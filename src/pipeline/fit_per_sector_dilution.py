"""PHOEBE fit with a separate dilution (flux offset and scale) for each TESS sector.

The daily pipeline fits one sector-consensus binned light curve whose errors are the
sector-to-sector scatter. In eclipse that scatter is 3-25x the white noise, because sectors
disagree on the eclipse depth (different apertures / pipelines dilute the eclipse by different
amounts). See docs/DD_EB_CROSSMATCH.md, cloud runs 2026-10-19 and 2026-10-20.

Here each sector is binned on its own (white-noise errors: std/sqrt(n) per bin) and compared with
one PHOEBE model through a sector-specific affine map, flux_s = a_s + b_s * model, solved by
weighted linear least squares at every model evaluation. b_s / mean(b) is the sector's relative
eclipse-depth scale (dilution). The 'common' chi2 uses one shared (a, b) for all sectors, i.e. the
same white-noise errors but no per-sector dilution, for comparison.

Needs per-point data with a sector column: time flux flux_err sector (e.g. fetch_tess_lc.py's
<des>_raw.txt, or J1122's outputs/j1122_tess_multisector_s10_s11_s37/..._LightCurve_multisector.txt).
The model grid is the times of the candidate's _binned.txt (same adaptive phase bins).
Global l3 is not fitted: it is degenerate with the per-sector affine terms.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
import run_phoebe_candidate as rpc  # noqa: E402

FREE = ["incl", "dt0", "rsum", "k", "teffratio"]
STEP = {"incl": 0.3, "dt0": 0.0005, "rsum": 0.004, "k": 0.03, "teffratio": 0.03, "ecosw": 0.001, "esinw": 0.01}
SMEAR = (-1 / 3, 0.0, 1 / 3)  # exposure average: 3 sub-exposures at the centres of its thirds


def robust_sigma(x):
    x = np.asarray(x, float)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def sector_bins(raw, edges, t0, period, min_points=3):
    """Bin each sector on the same phase edges. Returns {sector: (bin_index, flux, err, n)}."""
    t, f, s = raw[:, 0], raw[:, 1], raw[:, 3]
    ph_all = rpc.phase_values(t, t0, period)
    out = {}
    for sec in np.unique(s):
        m = s == sec
        idx = np.searchsorted(edges, ph_all[m]) - 1
        ok = (idx >= 0) & (idx < len(edges) - 1)
        fi, ii, pi = f[m][ok], idx[ok], ph_all[m][ok]
        rows = []
        for j in np.unique(ii):
            sel = ii == j
            v, pv = fi[sel], pi[sel]
            if len(v) >= min_points:
                med = np.median(v)
                keep = np.abs(v - med) < 5 * max(robust_sigma(v), 1e-6)
                v, pv = v[keep], pv[keep]
                rows.append((pv.mean(), v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v)))
        r = np.array(rows)
        # floor at the sector's median bin error, so bins with a few points cannot dominate
        r[:, 2] = np.maximum(r[:, 2], np.median(r[:, 2]))
        # each sector bin is compared with the model at its own mean phase (interpolated from
        # the model grid), not at the grid point, so bin-centre offsets do not bias the fit
        out[int(sec)] = (r[:, 0], r[:, 1], r[:, 2], r[:, 3].astype(int))
    return out


def affine_chi2(pred, bins, common=False, scale_prior=None):
    """chi2 of flux_s = a_s + b_s * model, per sector (or one shared a, b if common).
    pred is a callable sector -> model flux at that sector's bins (exposure-averaged).
    scale_prior = sigma adds a penalty ((mean_s b_s - 1) / sigma)^2, i.e. on average the
    model's own dilution (its l3) is right and sectors only scatter around it; pred must
    then be normalised to 1 out of eclipse."""
    if scale_prior:
        ss = list(bins)
        ncol = 2 * len(ss)
        Xs, ys = [], []
        for i, s in enumerate(ss):
            j, f, e, _ = bins[s]
            X = np.zeros((len(j), ncol))
            X[:, 2 * i] = 1 / e
            X[:, 2 * i + 1] = pred(s) / e
            Xs.append(X); ys.append(f / e)
        prow = np.zeros((1, ncol))
        prow[0, 1::2] = 1 / (len(ss) * scale_prior)
        X = np.vstack(Xs + [prow]); y = np.concatenate(ys + [[1 / scale_prior]])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(np.sum((y - X @ coef) ** 2)), {s: (coef[2 * i], coef[2 * i + 1]) for i, s in enumerate(ss)}
    if common:
        X = np.concatenate([np.c_[np.ones(len(b[0])), pred(s)] / b[2][:, None] for s, b in bins.items()])
        y = np.concatenate([f / e for _, f, e, _ in bins.values()])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(np.sum((y - X @ coef) ** 2)), {s: tuple(coef) for s in bins}
    total, coefs = 0.0, {}
    for s, (j, f, e, _) in bins.items():
        X = np.c_[np.ones(len(j)), pred(s)] / e[:, None]
        y = f / e
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        total += float(np.sum((y - X @ coef) ** 2))
        coefs[s] = tuple(coef)
    return total, coefs


def plot_sectors(path, name, bins, pred, coefs, res, phase2=0.5):
    """Per-sector data and model (after each sector's affine map) around both eclipses, with
    residuals in units of each bin's error."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p2 = phase2 if phase2 < 0.8 else phase2 - 1
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex="col", gridspec_kw={"height_ratios": [3, 1.3]})
    for col, (centre, lab) in enumerate(((0.0, "primary"), (p2, "secondary"))):
        ax, rx = axes[0, col], axes[1, col]
        for k, (s, (ph, f, e, _)) in enumerate(bins.items()):
            a, b = coefs[s]
            mod = a + b * pred(s)
            o = np.argsort(ph)
            m = np.abs(ph - centre) < 0.08
            off = -0.004 * k
            ax.errorbar(ph[m], f[m] + off, yerr=e[m], fmt=".", ms=3, lw=0.5, color=f"C{k}", label=f"S{s} (scale {b / res['mean_depth_scale_b']:.2f})")
            mo = o[np.abs(ph[o] - centre) < 0.08]
            ax.plot(ph[mo], mod[mo] + off, color=f"C{k}", lw=1.2)
            rx.plot(ph[m], (f[m] - mod[m]) / e[m], ".", ms=3, color=f"C{k}")
        ax.set_title(f"{lab} eclipse")
        rx.axhline(0, color="0.5", lw=0.8)
        rx.set_xlabel("phase"); rx.set_ylabel("resid / err")
    axes[0, 0].set_ylabel("flux (sectors offset)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"{name}: per-sector dilution fit, chi2_red {res['chi2_red_per_sector']:.2f} "
                 f"(shared dilution {res['chi2_red_common_same_model']:.2f})", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--designation", required=True)
    ap.add_argument("--raw", required=True, help="per-point file: time flux flux_err sector")
    ap.add_argument("--start-from", required=True, help="a _phoebe_result*.json or _localopt*.json to start from")
    ap.add_argument("--teff", type=float, required=True)
    ap.add_argument("--q", type=float, default=0.8)
    ap.add_argument("--ntri", type=int, default=400)
    ap.add_argument("--maxfev", type=int, default=400)
    ap.add_argument("--no-irrad", action="store_true")
    ap.add_argument("--tag", default="_persector")
    ap.add_argument("--eccentric", action="store_true", help="also fit ecosw, esinw (start values from --start-from)")
    ap.add_argument("--no-smear", action="store_true",
                    help="do not average the model over the exposure for sectors with cadence > 5 min")
    ap.add_argument("--l3", type=float, default=0.0, help="third-light fraction in the model (e.g. from TIC contamination)")
    ap.add_argument("--scale-prior", type=float, default=None,
                    help="sigma of a prior that the mean per-sector depth scale is 1 (needs a realistic --l3)")
    ap.add_argument("--walkers", type=int, default=0, help="emcee walkers after the optimisation (0 = no emcee)")
    ap.add_argument("--steps", type=int, default=0, help="emcee steps per walker")
    args = ap.parse_args()

    name = args.designation
    cdir, eph, times, flux, sig = rpc.load_inputs(name)
    period = eph["period"]
    t0 = eph["t0"]
    free = FREE + (["ecosw", "esinw"] if args.eccentric else [])
    rpc.PARAMS[:] = list(free)
    m1, _ = rpc.ms_mass_radius(args.teff)
    distortion = "roche" if period < 3 else "sphere"
    sf = json.loads(Path(args.start_from).read_text())
    if "v" in sf:
        d = dict(zip(sf["params"], sf["v"]))
    else:  # a _phoebe_result json
        d = dict(sf)
        d["dt0"] = sf["t0"] - t0
    v0 = np.array([d[n] for n in free], float)

    # phase bins: midpoints between the candidate's adaptive _binned.txt phases
    gph = np.sort(rpc.phase_values(times, t0, period))
    mids = 0.5 * (gph[1:] + gph[:-1])
    edges = np.r_[gph[0] - (mids[0] - gph[0]), mids, gph[-1] + (gph[-1] - mids[-1])]
    raw = np.genfromtxt(args.raw, comments="#")
    raw = raw[np.isfinite(raw).all(axis=1)]
    bins = sector_bins(raw, edges, t0, period)
    # PHOEBE is evaluated exactly at every sector bin's mean phase (no interpolation). For long
    # cadences (> 5 min) the model is averaged over the exposure (sub-exposure offsets SMEAR).
    secs = raw[:, 3]
    cadence = {s: float(np.median(np.diff(np.sort(raw[secs == s, 0])))) for s in bins}
    offs = {s: (np.array(SMEAR) * cadence[s] / period if (cadence[s] > 5 / 1440 and not args.no_smear)
                else np.array([0.0])) for s in bins}
    want = {s: (bins[s][0][:, None] + offs[s][None, :]) for s in bins}   # (nbin, nsub)
    allph = np.concatenate([w.ravel() for w in want.values()])
    uph, inv = np.unique(np.round(allph, 8), return_inverse=True)
    starts = np.cumsum([0] + [w.size for w in want.values()])
    where = {s: inv[starts[i]:starts[i + 1]].reshape(want[s].shape) for i, s in enumerate(bins)}
    gt = t0 + uph * period
    m = rpc.Model(gt, np.ones_like(gt), np.full_like(gt, 1e-3), period, t0, args.teff, m1, args.q, distortion,
                  args.ntri, l3=args.l3, irrad=not args.no_irrad)
    m.b.set_value("pblum_mode@lc01", "component-coupled")
    oot_idx = np.abs(rpc.phase_values(gt, t0, period)) > 0.15   # for normalising to 1 out of eclipse

    def interp_model(pred_grid):
        pg = pred_grid / np.median(pred_grid[oot_idx]) if oot_idx.any() else pred_grid
        return lambda s: pg[where[s]].mean(axis=1)
    nbin = sum(len(j) for j, *_ in bins.values())
    log_path = cdir / f"{name}_phoebe{args.tag}.log"
    logf = log_path.open("a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"=== {name} per-sector dilution fit: sectors {sorted(bins)}, {nbin} sector-bins, "
        f"{len(uph)} model phases, Teff1={args.teff:.0f}, l3={args.l3}, scale prior {args.scale_prior}, "
        f"free {free}, cadence (min) {({s: round(c * 1440, 1) for s, c in cadence.items()})}, "
        f"smeared {[s for s in bins if len(offs[s]) > 1]}, start {args.start_from}")
    best = {"chi2": np.inf, "v": v0}
    n = {"evals": 0}

    def model_pred(v):
        m.set(v)
        m.b.run_compute(compute="phoebe01", model="ps", overwrite=True, progressbar=False)
        return interp_model(np.asarray(m.b.get_value("fluxes@lc01@ps@model"), float))

    def objective(v):
        try:
            pred = model_pred(v)
        except Exception:
            return 1e12
        c, _ = affine_chi2(pred, bins, scale_prior=args.scale_prior)
        n["evals"] += 1
        if c < best["chi2"]:
            best.update(chi2=c, v=np.array(v), pred=pred)
        if n["evals"] % 25 == 0:
            log(f"  eval {n['evals']}: best chi2 {best['chi2']:.1f}")
        return c

    t_start = time.time()
    c0 = objective(v0)
    pred0 = best["pred"]
    c0_common, _ = affine_chi2(pred0, bins, common=True)
    log(f"  start: chi2 per-sector {c0:.1f}, common {c0_common:.1f} ({nbin} bins)")
    simplex = np.array([v0] + [v0 + np.eye(len(v0))[i] * STEP[p] for i, p in enumerate(free)])
    if args.maxfev > 1:
        minimize(objective, v0, method="Nelder-Mead",
                 options={"initial_simplex": simplex, "maxfev": args.maxfev, "xatol": 1e-5, "fatol": 0.05})
    v, pred = best["v"], best["pred"]
    c_ps, coefs = affine_chi2(pred, bins, scale_prior=args.scale_prior)
    c_common, coef_c = affine_chi2(pred, bins, common=True)
    nfree_ps = len(free) + 2 * len(bins)
    nfree_c = len(free) + 2
    bmean = np.mean([b for _, b in coefs.values()])
    res = {
        "designation": name, "teff1": args.teff, "params": free, "v": list(map(float, v)),
        "l3": args.l3, "scale_prior": args.scale_prior, "smeared_sectors": [s for s in bins if len(offs[s]) > 1],
        "mean_depth_scale_b": float(np.mean([b for _, b in coefs.values()])),
        "start_from": args.start_from, "n_sector_bins": nbin, "sectors": sorted(bins),
        "chi2_start_per_sector": c0, "chi2_start_common": c0_common,
        "chi2_per_sector": c_ps, "dof_per_sector": nbin - nfree_ps, "chi2_red_per_sector": c_ps / (nbin - nfree_ps),
        "chi2_common_same_model": c_common, "chi2_red_common_same_model": c_common / (nbin - nfree_c),
        "sector_affine": {str(s): {"a": a, "b": b, "depth_scale": b / bmean} for s, (a, b) in coefs.items()},
        "evals": n["evals"], "runtime_s": time.time() - t_start,
    }
    (cdir / f"{name}_persector{args.tag}.json").write_text(json.dumps(res, indent=1))
    plot_sectors(cdir / f"{name}_persector{args.tag}.png", name, bins, pred, coefs, res, phase2=eph.get("eclipses", {}).get("phase2", 0.5))
    log(f"  done: chi2 per-sector {c_ps:.1f} (red {res['chi2_red_per_sector']:.2f}), "
        f"common {c_common:.1f} (red {res['chi2_red_common_same_model']:.2f}); "
        "depth scales " + ", ".join(f"S{s} {x['depth_scale']:.3f}" for s, x in res['sector_affine'].items()))
    log(f"  params {dict(zip(free, np.round(v, 5)))}; mean b {res['mean_depth_scale_b']:.3f}")
    if args.walkers and args.steps:
        run_mcmc(args, name, cdir, free, v, c_ps, res, bins, model_pred, log)


def run_mcmc(args, name, cdir, free, v_best, chi2_best, res, bins, model_pred, log):
    """emcee around the optimum. The chi2 is divided by max(1, chi2_red) so the posterior width
    reflects the actual misfit, as in run_phoebe_candidate.py. Depth scales are stored as blobs."""
    import emcee
    dof = res["dof_per_sector"]
    scale = max(1.0, chi2_best / dof)
    secs = list(bins)
    lo = {"incl": 60.0, "dt0": -0.05, "rsum": 0.01, "k": 0.02, "teffratio": 0.1, "ecosw": -0.2, "esinw": -0.7}
    hi = {"incl": 120.0, "dt0": 0.05, "rsum": 0.8, "k": 3.0, "teffratio": 2.0, "ecosw": 0.2, "esinw": 0.7}
    nblob = len(secs)

    def lnp(v):
        if any(not (lo[p] < x < hi[p]) for p, x in zip(free, v)):
            return -np.inf, np.full(nblob, np.nan)
        try:
            pred = model_pred(v)
        except Exception:
            return -np.inf, np.full(nblob, np.nan)
        c, coefs = affine_chi2(pred, bins, scale_prior=args.scale_prior)
        bs = np.array([coefs[s][1] for s in secs])
        return -0.5 * c / scale, bs / bs.mean()

    rng = np.random.default_rng(3)
    step = np.array([0.2 * STEP[p] for p in free])
    p0 = v_best + step * rng.normal(size=(args.walkers, len(free)))
    # checkpoint the chain every 5 steps; a container restart (10-23) lost three chains at step 20/60.
    # On a rerun with the same tag, the walkers restart from the last saved positions and the
    # saved steps are kept (prepended) in the final chain.
    ck = cdir / f"{name}_persector_chainpart{args.tag}.npz"
    prev_chain = prev_blobs = None
    if ck.exists():
        z = np.load(ck)
        prev_chain, prev_blobs = z["chain"], z["blobs"]
        p0 = prev_chain[-1]
        log(f"  emcee: resuming from {ck.name} ({len(prev_chain)} saved steps)")
    nsteps = args.steps - (0 if prev_chain is None else len(prev_chain))
    t = time.time()
    sampler = emcee.EnsembleSampler(args.walkers, len(free), lnp)
    for i, _ in enumerate(sampler.sample(p0, iterations=max(nsteps, 0))):
        if (i + 1) % 5 == 0:
            log(f"  emcee step {i + 1}/{nsteps}, acc {np.mean(sampler.acceptance_fraction):.2f}")
            ch, bl = sampler.get_chain(), sampler.get_blobs()
            if prev_chain is not None:
                ch, bl = np.concatenate([prev_chain, ch]), np.concatenate([prev_blobs, bl])
            np.savez(ck, chain=ch, blobs=bl)
    if nsteps > 0:
        full_chain, full_blobs = sampler.get_chain(), sampler.get_blobs()
        if prev_chain is not None:
            full_chain = np.concatenate([prev_chain, full_chain])
            full_blobs = np.concatenate([prev_blobs, full_blobs])
    else:
        full_chain, full_blobs = prev_chain, prev_blobs
    burn = len(full_chain) // 2
    chain = full_chain[burn:].reshape(-1, len(free))
    blobs = full_blobs[burn:].reshape(-1, nblob)
    np.save(cdir / f"{name}_persector_chain{args.tag}.npy", full_chain)
    pct = {p: list(map(float, np.percentile(chain[:, i], [16, 50, 84]))) for i, p in enumerate(free)}
    ok = np.isfinite(blobs).all(axis=1)
    pct.update({f"depth_scale_S{s}": list(map(float, np.percentile(blobs[ok, i], [16, 50, 84])))
                for i, s in enumerate(secs)})
    res["mcmc"] = {"walkers": args.walkers, "steps": args.steps, "burn": burn, "chi2_scale": scale,
                   "acceptance": float(np.mean(sampler.acceptance_fraction)) if nsteps > 0 else float("nan"),
                   "percentiles_16_50_84": pct,
                   "runtime_s": time.time() - t}
    (cdir / f"{name}_persector{args.tag}.json").write_text(json.dumps(res, indent=1))
    log(f"  emcee {args.walkers}x{args.steps} (chi2 scale {scale:.2f}, acc {res['mcmc']['acceptance']:.2f}): "
        + "; ".join(f"{p} {v[1]:.4g} +{v[2] - v[1]:.2g} -{v[1] - v[0]:.2g}" for p, v in pct.items()))


if __name__ == "__main__":
    main()
