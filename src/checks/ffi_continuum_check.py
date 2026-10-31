#!/usr/bin/env python3
"""Is the out-of-eclipse shape of an FFI light curve real, or left by the detrending?

Re-does the FFI photometry for each sector and detrends it three ways:
  current  fetch_tess_ffi_s3.detrend() (per-orbit poly(2) + Fourier(4), poly divided out)
  linear   per-orbit straight line fitted to out-of-eclipse points only
  median   per-orbit median of out-of-eclipse points only (no trend removed)
then folds the out-of-eclipse points and fits a0 + Σ_{h=1,2} [a_h cos(2πhφ) + b_h sin(2πhφ)].
Ellipsoidal variation shows as a_2 < 0 (minima at both conjunctions); reflection as a_1 < 0.
If the shape is the same for all three detrends and all sectors, it is in the star.
With <cand>/<des>_ooe_models.npz (src/checks/ooe_model_curves.py) the PHOEBE sphere and Roche
curves are fitted the same way and overplotted.

  .venv-tess/bin/python src/checks/ffi_continuum_check.py --designation J231201.40+532028.6 --ra 348.0058494 --dec 53.3413022
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
import fetch_tess_ffi_s3 as ffi  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]

ap = argparse.ArgumentParser()
ap.add_argument("--designation", required=True)
ap.add_argument("--ra", type=float, required=True)
ap.add_argument("--dec", type=float, required=True)
ap.add_argument("--mask", type=float, default=None, help="eclipse half-width in phase (default: from the binned LC)")
args = ap.parse_args()
name = args.designation
cdir = ROOT / "outputs" / "candidates" / name
eph = json.loads((cdir / "ffi" / f"{name}_ephemeris.json").read_text())
P, t0 = eph["period"], eph["t0"]

# eclipse windows from the binned FFI light curve: where it dips > 0.3% below 1
binned = np.genfromtxt(cdir / "ffi" / f"{name}_binned.txt", comments="#")
ph_b, f_b = binned[:, 3], binned[:, 1]   # columns: time flux err phase n nsector
low = ph_b[f_b < 0.997]
sec_ph = np.median(low[np.abs(low) > 0.25] % 1.0) if (np.abs(low) > 0.25).any() else 0.5
hw = args.mask or float(max(np.abs(low[np.abs(low) < 0.25]).max(), np.abs((low[np.abs(low) > 0.25] % 1.0) - sec_ph).max()) + 0.01)


def phase(t):
    return ((t - t0) / P + 0.5) % 1.0 - 0.5


def ooe(ph):
    d2 = np.abs(((ph - sec_ph) + 0.5) % 1.0 - 0.5)
    return (np.abs(ph) > hw) & (d2 > hw)


def harm(ph, y):
    A = np.vstack([np.ones_like(ph)] + [fn(2 * np.pi * h * ph) for h in (1, 2) for fn in (np.cos, np.sin)]).T
    m = np.ones(len(y), bool)
    for _ in range(4):
        c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
        r = y - A @ c
        s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m])))
        m = np.abs(r) < 4 * s
    cov = np.linalg.inv(A[m].T @ A[m]) * s ** 2
    return c / c[0], np.sqrt(np.diag(cov)) / c[0], s


files = sorted((cdir / "ffi" / "ffi_cutouts").glob("tess-s*_astrocut.fits"))
rows, curves = [], {}
for fn in files:
    sec = int(fn.name[6:10])
    t, f, bg, *_ = ffi.photometry(fn, args.ra, args.dec)
    keep, trend = ffi.detrend(t, f, bg, P)
    ph = phase(t)
    o = ooe(ph)
    breaks = np.where(np.diff(t) > 0.75)[0] + 1
    lin = np.full(len(t), np.nan)
    med = np.full(len(t), np.nan)
    for idx in np.split(np.arange(len(t)), breaks):
        g = idx[keep[idx] & o[idx]]
        if len(g) < 30:
            continue
        c = np.polyfit(t[g] - t[g].mean(), f[g], 1)
        lin[idx] = f[idx] / np.polyval(c, t[idx] - t[g].mean())
        med[idx] = f[idx] / np.median(f[g])
    for meth, y in (("current", f / trend), ("linear", lin), ("median", med)):
        g = keep & o & np.isfinite(y)
        c, e, s = harm(ph[g], y[g])
        rows.append(dict(sector=sec, method=meth, n=int(g.sum()), a1=c[1], a1e=e[1], b1=c[2], a2=c[3], a2e=e[3], b2=c[4], rms=s))
        bins = np.linspace(-0.5, 0.5, 41)
        k = np.digitize(ph[g], bins)
        yb = np.array([np.median(y[g][k == i]) if (k == i).sum() > 5 else np.nan for i in range(1, len(bins))])
        curves[(sec, meth)] = (0.5 * (bins[1:] + bins[:-1]), yb / np.nanmedian(y[g]))

mods = {}
mf = cdir / f"{name}_ooe_models.npz"
if mf.exists():
    z = np.load(mf)
    for lab in z.files:
        if lab == "phase":
            continue
        g = ooe(z["phase"])
        c, e, s = harm(z["phase"][g], z[lab][g])
        mods[lab] = (z["phase"], z[lab])
        rows.append(dict(sector=0, method=f"model_{lab}", n=int(g.sum()), a1=c[1], a1e=0, b1=c[2], a2=c[3], a2e=0, b2=c[4], rms=0))

out = cdir / "ffi" / f"{name}_continuum_check"
with open(f"{out}.csv", "w") as fh:
    keys = list(rows[0])
    fh.write(",".join(keys) + "\n")
    for r in rows:
        fh.write(",".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
print(f"eclipse mask: |phase| > {hw:.3f} around 0 and {sec_ph:.4f}")
print(f"{'sector':>6} {'method':>16} {'n':>6} {'a1 (ppt)':>14} {'a2 (ppt)':>14} {'b1':>8} {'b2':>8}")
for r in rows:
    print(f"{r['sector']:>6} {r['method']:>16} {r['n']:>6} {1e3*r['a1']:7.2f}±{1e3*r['a1e']:4.2f} {1e3*r['a2']:7.2f}±{1e3*r['a2e']:4.2f} {1e3*r['b1']:8.2f} {1e3*r['b2']:8.2f}")

secs = sorted({s for s, _ in curves})
fig, axes = plt.subplots(1, len(secs), figsize=(4.5 * len(secs), 3.8), sharey=True, squeeze=False)
for ax, sec in zip(axes[0], secs):
    for meth, col in (("current", "k"), ("linear", "C0"), ("median", "C1")):
        x, y = curves[(sec, meth)]
        ax.plot(x, y, "o-", ms=3, color=col, label=meth)
    for lab, (x, y), col in zip(mods, mods.values(), ("C3", "C2", "C4")):
        g = ooe(x)
        ax.plot(np.where(g, x, np.nan), np.where(g, y, np.nan), "-", color=col, lw=1.5, label=f"PHOEBE {lab}")
    for c in (0, sec_ph, sec_ph - 1):
        ax.axvspan(c - hw, c + hw, color="0.85")
    ax.set_title(f"S{sec}")
    ax.set_xlabel("phase")
axes[0][0].set_ylabel("normalised flux (out of eclipse)")
lo = min(np.nanmin(y) for _, y in curves.values())
hi = max(np.nanmax(y) for _, y in curves.values())
axes[0][0].set_ylim(lo - 0.002, hi + 0.002)
axes[0][-1].legend(fontsize=7)
fig.suptitle(f"{name}: out-of-eclipse shape vs detrending (eclipses shaded)")
fig.tight_layout()
fig.savefig(f"{out}.png", dpi=120)
print(f"wrote {out}.csv/.png")
