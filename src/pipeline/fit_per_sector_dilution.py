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
STEP = {"incl": 0.3, "dt0": 0.0005, "rsum": 0.004, "k": 0.03, "teffratio": 0.03}


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


def affine_chi2(pred, bins, common=False):
    """chi2 of flux_s = a_s + b_s * model, per sector (or one shared a, b if common).
    pred is a callable phase -> model flux."""
    if common:
        X = np.concatenate([np.c_[np.ones(len(j)), pred(j)] / e[:, None] for j, _, e, _ in bins.values()])
        y = np.concatenate([f / e for _, f, e, _ in bins.values()])
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(np.sum((y - X @ coef) ** 2)), {s: tuple(coef) for s in bins}
    total, coefs = 0.0, {}
    for s, (j, f, e, _) in bins.items():
        X = np.c_[np.ones(len(j)), pred(j)] / e[:, None]
        y = f / e
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        total += float(np.sum((y - X @ coef) ** 2))
        coefs[s] = tuple(coef)
    return total, coefs


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
    args = ap.parse_args()

    name = args.designation
    cdir, eph, times, flux, sig = rpc.load_inputs(name)
    period = eph["period"]
    t0 = eph["t0"]
    rpc.PARAMS[:] = list(FREE)
    m1, _ = rpc.ms_mass_radius(args.teff)
    distortion = "roche" if period < 3 else "sphere"
    sf = json.loads(Path(args.start_from).read_text())
    if "v" in sf:
        d = dict(zip(sf["params"], sf["v"]))
    else:  # a _phoebe_result json
        d = dict(sf)
        d["dt0"] = sf["t0"] - t0
    v0 = np.array([d[n] for n in FREE], float)

    # phase bins: midpoints between the candidate's adaptive _binned.txt phases
    gph = np.sort(rpc.phase_values(times, t0, period))
    mids = 0.5 * (gph[1:] + gph[:-1])
    edges = np.r_[gph[0] - (mids[0] - gph[0]), mids, gph[-1] + (gph[-1] - mids[-1])]
    raw = np.genfromtxt(args.raw, comments="#")
    raw = raw[np.isfinite(raw).all(axis=1)]
    bins = sector_bins(raw, edges, t0, period)
    # PHOEBE is evaluated exactly at every sector bin's mean phase (no interpolation)
    allph = np.concatenate([j for j, *_ in bins.values()])
    uph, inv = np.unique(np.round(allph, 7), return_inverse=True)
    starts = np.cumsum([0] + [len(j) for j, *_ in bins.values()])
    where = {s: inv[starts[i]:starts[i + 1]] for i, s in enumerate(bins)}
    gt = t0 + uph * period
    m = rpc.Model(gt, np.ones_like(gt), np.full_like(gt, 1e-3), period, t0, args.teff, m1, args.q, distortion,
                  args.ntri, l3=0.0, irrad=not args.no_irrad)
    m.b.set_value("pblum_mode@lc01", "component-coupled")
    phase_index = {}
    for s, (j, *_r) in bins.items():
        phase_index.update(zip(np.round(j, 7), where[s]))

    def interp_model(pred_grid):
        return lambda phase: pred_grid[[phase_index[x] for x in np.round(phase, 7)]]
    nbin = sum(len(j) for j, *_ in bins.values())
    log_path = cdir / f"{name}_phoebe{args.tag}.log"
    logf = log_path.open("a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"=== {name} per-sector dilution fit: sectors {sorted(bins)}, {nbin} sector-bins, "
        f"Teff1={args.teff:.0f}, start {args.start_from}")
    best = {"chi2": np.inf, "v": v0}
    n = {"evals": 0}

    def objective(v):
        try:
            m.set(v)
            m.b.run_compute(compute="phoebe01", model="ps", overwrite=True, progressbar=False)
            pred = interp_model(np.asarray(m.b.get_value("fluxes@lc01@ps@model"), float))
        except Exception:
            return 1e12
        c, _ = affine_chi2(pred, bins)
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
    simplex = np.array([v0] + [v0 + np.eye(len(v0))[i] * STEP[p] for i, p in enumerate(FREE)])
    minimize(objective, v0, method="Nelder-Mead",
             options={"initial_simplex": simplex, "maxfev": args.maxfev, "xatol": 1e-5, "fatol": 0.05})
    v, pred = best["v"], best["pred"]
    c_ps, coefs = affine_chi2(pred, bins)
    c_common, coef_c = affine_chi2(pred, bins, common=True)
    nfree_ps = len(FREE) + 2 * len(bins)
    nfree_c = len(FREE) + 2
    bmean = np.mean([b for _, b in coefs.values()])
    res = {
        "designation": name, "teff1": args.teff, "params": FREE, "v": list(map(float, v)),
        "start_from": args.start_from, "n_sector_bins": nbin, "sectors": sorted(bins),
        "chi2_start_per_sector": c0, "chi2_start_common": c0_common,
        "chi2_per_sector": c_ps, "dof_per_sector": nbin - nfree_ps, "chi2_red_per_sector": c_ps / (nbin - nfree_ps),
        "chi2_common_same_model": c_common, "chi2_red_common_same_model": c_common / (nbin - nfree_c),
        "sector_affine": {str(s): {"a": a, "b": b, "depth_scale": b / bmean} for s, (a, b) in coefs.items()},
        "evals": n["evals"], "runtime_s": time.time() - t_start,
    }
    (cdir / f"{name}_persector{args.tag}.json").write_text(json.dumps(res, indent=1))
    log(f"  done: chi2 per-sector {c_ps:.1f} (red {res['chi2_red_per_sector']:.2f}), "
        f"common {c_common:.1f} (red {res['chi2_red_common_same_model']:.2f}); "
        "depth scales " + ", ".join(f"S{s} {x['depth_scale']:.3f}" for s, x in res['sector_affine'].items()))
    log(f"  params {dict(zip(FREE, np.round(v, 5)))}")


if __name__ == "__main__":
    main()
