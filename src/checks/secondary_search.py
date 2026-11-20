#!/usr/bin/env python3
"""Secondary-eclipse search with the out-of-eclipse variation fitted at the same time.

Model, fitted by linear least squares to all out-of-primary points of the chosen sectors:
  flux = per-orbit polynomial (slow trends) + per-sector Fourier series in orbital phase (the wave)
         - depth * trapezoid(phase - phi2)
The trapezoid has the primary's total duration and ingress. For each trial phi2 (0.2-0.8) the
depth and its error are fitted; errors are scaled by the red-noise factor beta (rms of binned
residuals vs white expectation). The scan over phases away from 0.5 gives the empirical false-alarm
distribution, and injected secondaries at phase 0.5 test what depth is recovered.

  .venv-tess/bin/python src/checks/secondary_search.py LC.txt P T0 OUTPREFIX [sectors e.g. 58,85]
LC.txt columns: time flux flux_err sector.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

lcf, P, T0, outp = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
secs = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else None
NHARM, DEG = int(__import__("os").environ.get("NHARM", 3)), 2
DUR, ING = 0.068, 0.025       # primary total duration and ingress, in phase (measured from the fold)

a = np.genfromtxt(lcf, comments="#")
t, f, s = a[:, 0], a[:, 1], a[:, 3].astype(int)
if secs:
    m = np.isin(s, secs)
    t, f, s = t[m], f[m], s[m]
ph = ((t - T0) / P) % 1.0
keep = np.abs(((ph + 0.5) % 1) - 0.5) > DUR / 2 + 0.01        # drop the primary eclipse
t, f, s, ph = t[keep], f[keep], s[keep], ph[keep]
# clip isolated outliers against a running median
rm = np.array([np.median(f[max(0, i - 5):i + 6]) for i in range(len(f))])
d = f - rm
sig0 = 1.4826 * np.median(np.abs(d))
ok = np.abs(d) < 5 * sig0
t, f, s, ph = t[ok], f[ok], s[ok], ph[ok]


def trap(x):
    x = np.abs(((x + 0.5) % 1) - 0.5)
    full = DUR / 2 - ING
    return np.clip((DUR / 2 - x) / ING, 0, 1) * (x < DUR / 2) + (x <= full) * (1 - np.clip((DUR / 2 - x) / ING, 0, 1))


def base_columns():
    cols = []
    seg = np.r_[0, np.where(np.diff(t) > 0.75)[0] + 1, len(t)]
    for i in range(len(seg) - 1):
        idx = slice(seg[i], seg[i + 1])
        tt = np.zeros(len(t)); m = np.zeros(len(t), bool); m[idx] = True
        x = (t - t[idx].mean()) / max(np.ptp(t[idx]), 1e-3)
        for k in range(DEG + 1):
            cols.append(np.where(m, x ** k, 0.0))
    for sec in np.unique(s):
        m = s == sec
        for h in range(1, NHARM + 1):
            cols += [np.where(m, np.cos(2 * np.pi * h * ph), 0), np.where(m, np.sin(2 * np.pi * h * ph), 0)]
    return np.vstack(cols).T


B = base_columns()


def fit(y, phi2):
    box = -trap(ph - phi2)
    A = np.c_[B, box]
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c
    sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    cov = np.linalg.pinv(A.T @ A) * sig ** 2
    return c[-1], np.sqrt(cov[-1, -1]), r, sig


d05, e05, r, sig = fit(f, 0.5)
# red-noise factor: binned residual rms over the white expectation, for bins of the eclipse duration
nb = max(2, int(DUR * P / np.median(np.diff(t))))
bins = [r[i:i + nb] for i in range(0, len(r) - nb, nb)]
beta = max(1.0, np.std([b.mean() for b in bins]) / (sig / np.sqrt(nb)))
grid = np.arange(0.2, 0.8001, 0.0025)
scan = np.array([fit(f, p)[:2] for p in grid])
away = np.abs(grid - 0.5) > DUR
fa = scan[away, 0]
print(f"{lcf} sectors {sorted(set(s))}: {len(t)} points")
print(f"secondary depth at phase 0.5: {d05*1e3:.3f} ± {e05*1e3:.3f} ppt (white), ± {e05*beta*1e3:.3f} ppt (×beta {beta:.2f})")
print(f"scan away from 0.5: depth rms {np.std(fa)*1e3:.3f} ppt, max {fa.max()*1e3:.3f}, min {fa.min()*1e3:.3f} ppt; "
      f"fraction of phases with depth >= the 0.5 value: {np.mean(fa >= d05):.3f}")
best = grid[np.argmax(scan[:, 0] / scan[:, 1])]
print(f"highest-S/N dip in the scan: phase {best:.4f}, depth {scan[np.argmax(scan[:, 0] / scan[:, 1]), 0]*1e3:.3f} ppt")
inj = []
for dep in (0.0005, 0.001, 0.002, 0.003):
    rec = fit(f - dep * trap(ph - 0.5), 0.5)[0]
    inj.append((dep, rec))
    print(f"injected {dep*1e3:.1f} ppt at 0.5 -> recovered {rec*1e3:.3f} ppt")
fig, ax = plt.subplots(2, 1, figsize=(10, 7))
ax[0].plot(grid, scan[:, 0] * 1e3, "k-")
ax[0].fill_between(grid, (scan[:, 0] - beta * scan[:, 1]) * 1e3, (scan[:, 0] + beta * scan[:, 1]) * 1e3, color="0.8")
ax[0].axvline(0.5, color="r", lw=0.8)
ax[0].axhline(0, color="0.5", lw=0.5)
ax[0].set_ylabel("fitted secondary depth (ppt)")
ax[0].set_xlabel("trial secondary phase")
ax[0].set_title(f"depth at 0.5 = {d05*1e3:.2f} ± {e05*beta*1e3:.2f} ppt; scan rms away from 0.5 = {np.std(fa)*1e3:.2f} ppt")
A = np.c_[B, -trap(ph - 0.5)]
c, *_ = np.linalg.lstsq(A, f, rcond=None)
det = f - B @ c[:-1]
edges = np.linspace(0, 1, 201)
k = np.digitize(ph, edges) - 1
xb = 0.5 * (edges[1:] + edges[:-1])
yb = np.array([np.mean(det[k == i]) if (k == i).sum() > 3 else np.nan for i in range(200)])
eb = np.array([np.std(det[k == i]) / np.sqrt((k == i).sum()) if (k == i).sum() > 3 else np.nan for i in range(200)])
ax[1].errorbar(xb, (yb - np.nanmedian(yb)) * 1e3, eb * 1e3, fmt="o", ms=3)
xx = np.linspace(0.2, 0.8, 600)
ax[1].plot(xx, -c[-1] * trap(xx - 0.5) * 1e3, "r-", label="fitted secondary")
ax[1].set_xlim(0.2, 0.8)
ax[1].set_xlabel("phase")
ax[1].set_ylabel("flux - trends - wave (ppt)")
ax[1].legend()
fig.tight_layout()
fig.savefig(outp + ".png", dpi=110)
with open(outp + ".txt", "w") as fh:
    fh.write(f"depth05_ppt {d05*1e3:.4f}\nerr_white_ppt {e05*1e3:.4f}\nbeta {beta:.3f}\nerr_ppt {e05*beta*1e3:.4f}\n"
             f"scan_rms_ppt {np.std(fa)*1e3:.4f}\nscan_max_ppt {fa.max()*1e3:.4f}\nscan_frac_ge {np.mean(fa >= d05):.4f}\n")
    for dep, rec in inj:
        fh.write(f"inject_{dep*1e3:.1f}ppt_recovered {rec*1e3:.4f}\n")
