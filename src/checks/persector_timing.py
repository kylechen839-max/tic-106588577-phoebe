"""Per-sector eclipse timing check for a per-sector dilution fit.

Evaluates the PHOEBE model of a finished fit_per_sector_dilution.py result once, on a fine phase
grid, then for each sector shifts the model in time (dt_s) and refits that sector's affine map.
Reports chi2 per sector at dt_s = 0 and at the best dt_s, plus the in-eclipse chi2 split, so a
residual pattern (ingress/egress antisymmetry) can be attributed to timing (period / ETV) or not.

Usage: persector_timing.py DESIGNATION RAW RESULT_JSON TEFF [--distortion roche] [--no-irrad]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
import run_phoebe_candidate as rpc  # noqa: E402
from fit_per_sector_dilution import sector_bins  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("designation")
    ap.add_argument("raw")
    ap.add_argument("result")
    ap.add_argument("teff", type=float)
    ap.add_argument("--distortion", default="roche")
    ap.add_argument("--no-irrad", action="store_true")
    ap.add_argument("--ntri", type=int, default=400, help="use the fit's own --ntri (J2312 fits: 1500)")
    ap.add_argument("--q", type=float, default=None, help="the fit's --q (sets the Kepler sma; default: q in the result)")
    ap.add_argument("--grid-step", type=float, default=0.0001, help="phase step of the model grid near eclipses")
    ap.add_argument("--period", type=float, default=None, help="the fit's --period, if it overrode the ephemeris")
    ap.add_argument("--maxshift-min", type=float, default=12.0)
    args = ap.parse_args()

    name = args.designation
    cdir, eph, times, flux, sig = rpc.load_inputs(name)
    P, t0 = (args.period or eph["period"]), eph["t0"]
    res = json.loads(Path(args.result).read_text())
    free = res["params"]
    v = np.array(res["v"], float)
    d = dict(zip(free, v))
    rpc.PARAMS[:] = list(free)
    m1, _ = rpc.ms_mass_radius(args.teff)

    gph = np.sort(rpc.phase_values(times, t0, P))
    mids = 0.5 * (gph[1:] + gph[:-1])
    edges = np.r_[gph[0] - (mids[0] - gph[0]), mids, gph[-1] + (gph[-1] - mids[-1])]
    raw = np.genfromtxt(args.raw, comments="#")
    raw = raw[np.isfinite(raw).all(axis=1)]
    bins = sector_bins(raw, edges, t0, P)

    # fine grid near eclipses, coarse elsewhere (phase in [-0.5, 0.5))
    st = args.grid_step
    fine = np.r_[np.arange(-0.05, 0.05, st), np.arange(0.45, 0.55, st)]
    coarse = np.arange(-0.5, 0.5, 0.005)
    ph = np.unique(np.round(np.r_[fine, coarse], 6))
    ph = ((ph + 0.5) % 1.0) - 0.5
    ph = np.unique(ph)
    gt = t0 + ph * P
    cache = cdir / f"{name}_timing_model_{Path(args.result).stem}_{args.grid_step:g}_n{args.ntri}.npz"
    if cache.exists():
        z = np.load(cache)
        mph, mf = z["mph"], z["mf"]
    else:
        mph, mf = None, None
    if mph is None:
      m = rpc.Model(gt, np.ones_like(gt), np.full_like(gt, 1e-3), P, t0, args.teff, m1, args.q if args.q is not None else d.get("q", 0.8),
                    args.distortion, args.ntri, l3=res.get("l3", 0.0), irrad=not args.no_irrad)
      m.b.set_value("pblum_mode@lc01", "component-coupled")
      m.set(v)
      m.b.run_compute(compute="phoebe01", model="ps", overwrite=True, progressbar=False)
      mf = np.asarray(m.b.get_value("fluxes@lc01@ps@model"), float)
      mph = rpc.phase_values(m.b.get_value("times@lc01@ps@model"), t0, P)
      np.savez(cache, mph=mph, mf=mf)
    o = np.argsort(mph)
    mph, mf = mph[o], mf[o]
    mf = mf / np.median(mf[np.abs(mph) > 0.15])
    xp = np.r_[mph - 1, mph, mph + 1]
    fp = np.r_[mf, mf, mf]

    def model_at(phs):
        return np.interp(phs, xp, fp)

    phase2 = eph.get("eclipses", {}).get("phase2", 0.5)
    p2 = phase2 if phase2 < 0.5 else phase2 - 1
    shifts = np.arange(-args.maxshift_min, args.maxshift_min + 1e-9, 0.1) / 1440.0
    out = {"designation": name, "result": args.result, "period": P, "sectors": {}}
    print(f"{name}: per-sector timing scan, model from {args.result}")
    for s, (sph, f, e, _) in bins.items():
        best = None
        rows = []
        for dt in shifts:
            mod = model_at(sph - dt / P)
            X = np.c_[np.ones(len(f)), mod] / e[:, None]
            coef, *_ = np.linalg.lstsq(X, f / e, rcond=None)
            r = (f - (coef[0] + coef[1] * mod)) / e
            c = float(np.sum(r ** 2))
            rows.append(c)
            if best is None or c < best[1]:
                best = (dt, c, r)
        rows = np.array(rows)
        c0 = rows[np.argmin(np.abs(shifts))]
        dtb, cb, rb = best
        # 1-sigma from delta chi2 = 1 after rescaling chi2 to red = 1
        dof = len(f) - 3
        scale = cb / dof
        ok = rows <= cb + scale
        lo, hi = shifts[ok].min(), shifts[ok].max()
        r0 = None
        mod0 = model_at(sph)
        X = np.c_[np.ones(len(f)), mod0] / e[:, None]
        coef, *_ = np.linalg.lstsq(X, f / e, rcond=None)
        r0 = (f - (coef[0] + coef[1] * mod0)) / e
        prim = np.abs(sph) < 0.04
        sec = np.abs(((sph - p2 + 0.5) % 1.0) - 0.5) < 0.04
        ing = (sph > -0.04) & (sph < -0.005)
        egr = (sph > 0.005) & (sph < 0.04)
        info = {
            "n_bins": int(len(f)), "chi2_dt0": float(c0), "chi2_best": float(cb),
            "dt_best_min": float(dtb * 1440), "dt_1sig_min": [float(lo * 1440), float(hi * 1440)],
            "chi2_prim_dt0": float(np.sum(r0[prim] ** 2)), "chi2_prim_best": float(np.sum(rb[prim] ** 2)),
            "chi2_sec_dt0": float(np.sum(r0[sec] ** 2)), "chi2_sec_best": float(np.sum(rb[sec] ** 2)),
            "mean_resid_ingress_dt0": float(np.mean(r0[ing])) if ing.any() else None,
            "mean_resid_egress_dt0": float(np.mean(r0[egr])) if egr.any() else None,
            "mean_resid_ingress_best": float(np.mean(rb[ing])) if ing.any() else None,
            "mean_resid_egress_best": float(np.mean(rb[egr])) if egr.any() else None,
            "t_mid": float(np.median(raw[raw[:, 3] == s, 0])),
        }
        out["sectors"][str(s)] = info
        fm = lambda x: "n/a" if x is None else f"{x:+.2f}"
        print(f"  S{s}: chi2 {c0:.1f} -> {cb:.1f} at dt {dtb * 1440:+.2f} min "
              f"(1sig {lo * 1440:+.1f}..{hi * 1440:+.1f}); primary {info['chi2_prim_dt0']:.1f} -> "
              f"{info['chi2_prim_best']:.1f}; secondary {info['chi2_sec_dt0']:.1f} -> {info['chi2_sec_best']:.1f}; "
              f"ingress/egress mean resid {fm(info['mean_resid_ingress_dt0'])}/{fm(info['mean_resid_egress_dt0'])} -> "
              f"{fm(info['mean_resid_ingress_best'])}/{fm(info['mean_resid_egress_best'])}")
    # linear ephemeris correction from the per-sector shifts (dt = a + b * E)
    ss = sorted(out["sectors"], key=int)
    tm = np.array([out["sectors"][s]["t_mid"] for s in ss])
    dt = np.array([out["sectors"][s]["dt_best_min"] for s in ss])
    err = np.array([max(1e-3, 0.5 * (out["sectors"][s]["dt_1sig_min"][1] - out["sectors"][s]["dt_1sig_min"][0])) for s in ss])
    E = np.round((tm - t0) / P)
    if len(set(E)) >= 2:
        W = 1 / err ** 2
        A = np.c_[np.ones_like(E), E]
        cov = np.linalg.inv(A.T @ (A * W[:, None]))
        coef = cov @ (A.T @ (W * dt))
        chi2_lin = float(np.sum(W * (dt - A @ coef) ** 2))
        out["linear_ephemeris"] = {"dT0_min": float(coef[0]), "dP_min_per_cycle": float(coef[1]),
                                   "dP_days": float(coef[1] / 1440), "dP_err_days": float(np.sqrt(cov[1, 1]) / 1440),
                                   "chi2": chi2_lin, "n": len(ss)}
        print(f"  linear fit to shifts: dT0 {coef[0]:+.2f} min, dP {coef[1] / 1440:+.3e} "
              f"+/- {np.sqrt(cov[1, 1]) / 1440:.1e} d (chi2 {chi2_lin:.1f} for {len(ss)} sectors)")
    tag = Path(args.result).stem.replace(f"{name}_persector", "")
    outp = cdir / f"{name}_timing{tag}.json"
    outp.write_text(json.dumps(out, indent=1))
    print(f"  wrote {outp}")


if __name__ == "__main__":
    main()
