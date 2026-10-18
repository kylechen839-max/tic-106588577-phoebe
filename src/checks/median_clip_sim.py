"""Simulate how the running-median clip in fetch_tess_lc.normalise() treats real eclipses.

The raw TESS points are not kept in the repo, so this injects each candidate's own binned
eclipse profile (interpolated in phase) into synthetic light curves at the TESS cadences,
adds white noise, applies the same clip (median_filter size 9, 6 sigma on the deviation,
plus the 5 sigma upward cut), and reports the fraction of in-eclipse points removed and
the change in the binned eclipse minimum. Offline; numpy + scipy only.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter

ROOT = Path(__file__).resolve().parents[2]


def robust_sigma(x):
    x = np.asarray(x, float)
    s = 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))
    return float(s) if np.isfinite(s) and s > 0 else float(np.nanstd(x))


def clip_old(ff):
    """The clip used before 2026-10-18."""
    s = robust_sigma(ff)
    keep = ff < 1 + 5 * s
    dev = ff - median_filter(ff, size=9, mode="nearest")
    keep &= np.abs(dev) < 6 * robust_sigma(dev)
    return keep & (ff > 0.05)


def clip_new(ff):
    """Current fetch_tess_lc.clip_mask (copied so this runs without lightkurve)."""
    s = robust_sigma(ff)
    keep = ff < 1 + 5 * s
    dev = ff - median_filter(ff, size=9, mode="nearest")
    sd = robust_sigma(dev)
    sgn = np.sign(dev) * (np.abs(dev) > 3 * sd)
    prev = np.r_[0, sgn[:-1]]
    nxt = np.r_[sgn[1:], 0]
    in_run = (sgn != 0) & ((prev == sgn) | (nxt == sgn))
    keep &= (np.abs(dev) < 6 * sd) | in_run
    return keep & (ff > 0.05)


clip = clip_old


def simulate(des, cadence_s, noise, ndays=27.0, seed=1):
    d = ROOT / "outputs" / "candidates" / des
    rows = np.loadtxt(d / f"{des}_binned.txt")
    eph = json.loads((d / f"{des}_ephemeris.json").read_text())
    ph, fl = rows[:, 3] % 1, rows[:, 1]
    o = np.argsort(ph)
    ph, fl = ph[o], fl[o]
    ph = np.r_[ph - 1, ph, ph + 1]; fl = np.r_[fl, fl, fl]
    P = eph["period"]
    rng = np.random.default_rng(seed)
    t = np.arange(0, ndays, cadence_s / 86400.0)
    phase = (t / P) % 1
    truth = np.interp(phase, ph, fl)
    f = truth + rng.normal(0, noise, t.size)
    # 13.7-d orbits, as normalise() treats each continuous segment separately
    keep = np.zeros(t.size, bool)
    for seg in np.split(np.arange(t.size), np.where(np.diff(np.floor(t / 13.7)) > 0)[0] + 1):
        med = np.median(f[seg])
        keep[seg] = clip(f[seg] / med)
    ecl = eph["eclipses"]
    out = {}
    for name, p0, hw in (("primary", ecl["phase1"], ecl["halfwidth1"]), ("secondary", ecl["phase2"], ecl["halfwidth2"])):
        dph = ((phase - p0 + 0.5) % 1) - 0.5
        inn = np.abs(dph) < hw
        core = np.abs(dph) < 0.25 * hw
        # depth bias in the eclipse core: mean residual (kept - truth) of the core points
        bias_all = np.mean(f[core] - truth[core])
        bias_kept = np.mean(f[core & keep] - truth[core & keep]) if np.any(core & keep) else np.nan
        out[name] = {"frac_removed_in_eclipse": float(1 - keep[inn].mean()),
                     "frac_removed_core": float(1 - keep[core].mean()),
                     "core_bias_ppt": float(1e3 * (bias_kept - bias_all))}
    out["frac_removed_oot"] = float(1 - keep[~np.any([np.abs(((phase - p + 0.5) % 1) - 0.5) < h for p, h in
                                                      ((ecl["phase1"], ecl["halfwidth1"]), (ecl["phase2"], ecl["halfwidth2"]))], axis=0)].mean())
    return out


if __name__ == "__main__":
    des_list = sys.argv[1:]
    if des_list and des_list[0] == "--new":
        clip = clip_new
        des_list = des_list[1:]
    res = []
    for des in des_list:
        for cad in (1800, 600, 200, 120):
            for noise in (0.002, 0.005, 0.01):
                r = simulate(des, cad, noise)
                res.append({"designation": des, "cadence_s": cad, "noise": noise, **{f"{k}_{kk}": vv for k, v in r.items() if isinstance(v, dict) for kk, vv in v.items()}, "frac_removed_oot": r["frac_removed_oot"]})
                p, s = r["primary"], r["secondary"]
                print(f"{des} cad={cad:5d}s noise={noise:.3f}: removed in-ecl {100*p['frac_removed_in_eclipse']:5.1f}%/{100*s['frac_removed_in_eclipse']:5.1f}%  "
                      f"core {100*p['frac_removed_core']:5.1f}%/{100*s['frac_removed_core']:5.1f}%  core bias {p['core_bias_ppt']:+6.2f}/{s['core_bias_ppt']:+6.2f} ppt  oot {100*r['frac_removed_oot']:.2f}%")
    outp = ROOT / "outputs" / "checks"; outp.mkdir(parents=True, exist_ok=True)
    (outp / ("median_clip_sim_new.json" if clip is clip_new else "median_clip_sim.json")).write_text(json.dumps(res, indent=1))
