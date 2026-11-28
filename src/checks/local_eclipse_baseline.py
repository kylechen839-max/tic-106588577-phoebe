"""Per-eclipse local baseline for long-period EB FFI light curves.

For each observed eclipse event (primary and secondary) in a per-point file
(time flux flux_err sector), fit a low-order polynomial in time to the
out-of-eclipse points inside a window around the event and divide it out.
Only points inside the windows are written, so sector-to-sector offsets and
slopes left by per-orbit detrending cannot enter the eclipse fit.

usage: local_eclipse_baseline.py RAW EPHEMERIS_JSON OUT [--win 0.08] [--pad 1.3] [--order 1]
"""
import argparse, json
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("raw"); ap.add_argument("eph"); ap.add_argument("out")
ap.add_argument("--win", type=float, default=0.08, help="half window in phase around each eclipse")
ap.add_argument("--pad", type=float, default=1.3, help="mask |dphase| < pad*halfwidth as in-eclipse")
ap.add_argument("--order", type=int, default=1)
ap.add_argument("--min-side", type=int, default=20, help="min OOT points on each side of an event")
a = ap.parse_args()

raw = np.genfromtxt(a.raw, comments="#"); raw = raw[np.isfinite(raw).all(1)]
e = json.load(open(a.eph)); P, t0 = e["period"], e["t0"]; ec = e["eclipses"]
t, f, s, sec = raw[:, 0], raw[:, 1], raw[:, 2], raw[:, 3]
keep = np.zeros(len(t), bool); fn = f.copy(); log = []
for ph_c, hw, lab in [(ec["phase1"], ec["halfwidth1"], "P"), (ec["phase2"], ec["halfwidth2"], "S")]:
    n = np.round((t - t0) / P - ph_c)
    for k in np.unique(n):
        tc = t0 + (k + ph_c) * P
        dph = (t - tc) / P
        w = (np.abs(dph) < a.win) & (n == k)
        if not w.any():
            continue
        oot = w & (np.abs(dph) > a.pad * hw)
        nl, nr = (oot & (dph < 0)).sum(), (oot & (dph > 0)).sum()
        if nl < a.min_side or nr < a.min_side:
            log.append((lab, int(k), int(np.median(sec[w])), int(nl), int(nr), "dropped (one-sided)"))
            continue
        x = t[oot] - tc
        c = np.polyfit(x, f[oot], a.order)
        for _ in range(3):  # 4-sigma clip of the baseline fit
            r = f[oot] - np.polyval(c, x); m = np.abs(r) < 4 * 1.4826 * np.median(np.abs(r - np.median(r)))
            c = np.polyfit(x[m], f[oot][m], a.order)
        fn[w] = f[w] / np.polyval(c, t[w] - tc)
        keep |= w
        log.append((lab, int(k), int(np.median(sec[w])), int(nl), int(nr), f"slope {c[0]*1e3:+.2f} ppt/d, level {c[-1]:.4f}"))
np.savetxt(a.out, np.c_[t, fn, s, sec][keep], fmt="%.8f", header="time flux flux_err sector (per-eclipse local baseline)")
for row in log:
    print(*row)
print(f"kept {keep.sum()} of {len(t)} points")
