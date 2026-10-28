"""Light-curve products from lightkurve_pipeline_s3.py output, mirroring the
outputs of Kyle's Colab notebook (light curve, Lomb-Scargle periodogram, fold).

Usage: lk_pipeline_products.py LKDIR DESIGNATION P T0 [good sectors, comma-separated]
"""
import sys, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

lkdir, des, P, T0 = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
good = [int(s) for s in sys.argv[5].split(",")] if len(sys.argv) > 5 else None

secs = {}
for fn in sorted(glob.glob(os.path.join(lkdir, f"{des}_S*_lk.txt"))):
    s = int(re.search(r"_S(\d+)_lk", fn).group(1))
    t, f, e, fn_ = np.loadtxt(fn, unpack=True)
    ok = np.isfinite(t) & np.isfinite(fn_)
    med = np.median(f[ok])
    secs[s] = dict(t=t[ok], f=fn_[ok], e=e[ok] / med, raw=f[ok])
good = good or sorted(secs)
out = lambda name: os.path.join(lkdir, f"{des}_{name}")

def phase(t):
    return ((t - T0) / P + 0.25) % 1 - 0.25          # primary at 0, range -0.25..0.75

def binned(ph, f, nb=200):
    edges = np.linspace(-0.25, 0.75, nb + 1)
    idx = np.digitize(ph, edges) - 1
    c, m, s, n = [], [], [], []
    for i in range(nb):
        sel = idx == i
        if sel.sum() < 3:
            continue
        c.append(0.5 * (edges[i] + edges[i + 1])); m.append(np.median(f[sel]))
        s.append(1.2533 * np.std(f[sel]) / np.sqrt(sel.sum())); n.append(sel.sum())
    return map(np.array, (c, m, s, n))

# 1. per-sector light curves
fig, axs = plt.subplots(len(secs), 1, figsize=(11, 2.2 * len(secs)))
for ax, (s, d) in zip(np.atleast_1d(axs), secs.items()):
    ax.plot(d["t"], d["f"], ".", ms=1.5, color="k" if s in good else "0.6")
    p2p = np.median(np.abs(np.diff(d["f"]))) * 1.4826 / np.sqrt(2) * 1e3
    ax.set_title(f"S{s}  {len(d['t'])} cadences  p2p rms {p2p:.2f} ppt"
                 + ("" if s in good else "  (not used in combined fold)"), fontsize=9, loc="left")
    ax.set_ylabel("norm. flux")
    ax.set_ylim(np.percentile(d["f"], 0.2) - 0.005, np.percentile(d["f"], 99.8) + 0.005)
np.atleast_1d(axs)[-1].set_xlabel("BTJD")
fig.suptitle(f"{des}: lightkurve RegressionCorrector light curves (TESScut 10x10, PCA-5 background)")
fig.tight_layout(); fig.savefig(out("LightCurve_sectors.png"), dpi=130); plt.close(fig)

# 2. combined exports and periodogram (good sectors)
T = np.concatenate([secs[s]["t"] for s in good]); F = np.concatenate([secs[s]["f"] for s in good])
E = np.concatenate([secs[s]["e"] for s in good]); S = np.concatenate([[s] * len(secs[s]["t"]) for s in good])
np.savetxt(out("LightCurve.txt"), np.c_[T, F, E, S], fmt="%.6f %.6f %.6f %d",
           header="Time(BTJD) Flux(norm per sector) Flux_Err Sector  [lightkurve pipeline; sectors %s]" % good)

freq = np.linspace(1 / 20, 1 / 0.2, 200000)
ls_pow = LombScargle(T, F, E).power(freq)
np.savetxt(out("Periodogram.txt"), np.c_[1 / freq, freq, ls_pow], fmt="%.7f %.7f %.6e",
           header="Period(d) Frequency(1/d) LS_power")
top = np.argsort(ls_pow)[::-1]
peaks = []
for i in top:
    if all(abs(1 / freq[i] - p) / p > 0.02 for p in peaks):
        peaks.append(1 / freq[i])
    if len(peaks) == 5:
        break
fig, ax = plt.subplots(figsize=(10, 3.5))
ax.plot(1 / freq, ls_pow, "k-", lw=0.6); ax.set_xscale("log")
for k, lab in [(1, "P"), (0.5, "P/2"), (2, "2P")]:
    ax.axvline(P * k, color="C1", ls=":", lw=1); ax.text(P * k, ax.get_ylim()[1] * 0.92, lab, color="C1")
ax.set_xlabel("Period (d)"); ax.set_ylabel("LS power")
ax.set_title("Lomb-Scargle, sectors %s; top peaks: %s d" % (good, ", ".join(f"{p:.4f}" for p in peaks)), fontsize=9)
fig.tight_layout(); fig.savefig(out("Periodogram.png"), dpi=130); plt.close(fig)

# 3. folds
ph = phase(T)
c, m, s_, n = binned(ph, F)
np.savetxt(out("Folded_LC.txt"), np.c_[ph, F, E, S], fmt="%.6f %.6f %.6f %d",
           header=f"Phase Flux Flux_Err Sector  [P={P} d, T0={T0} BTJD, primary at phase 0]")
np.savetxt(out("Folded_LC_binned.txt"), np.c_[c, m, s_, n], fmt="%.5f %.6f %.6f %d",
           header="Phase Flux_median Flux_err Npts  [200 bins]")

fig, axs = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True)
ax = axs[0]
ax.plot(ph, F, ".", ms=1, color="0.75", rasterized=True)
ax.errorbar(c, m, s_, fmt="o", ms=2.5, color="k", lw=0.8)
ax.set_ylabel("norm. flux"); ax.set_title(f"{des}  fold P = {P} d, T0 = {T0}  (sectors {good})", fontsize=10)
lo = np.percentile(F, 0.3); ax.set_ylim(lo - 0.003, np.percentile(F, 99.8) + 0.003)
ax = axs[1]
for k, sct in enumerate(sorted(secs)):
    d = secs[sct]; cc, mm, ss, nn = binned(phase(d["t"]), d["f"], 100)
    ax.plot(cc, mm - 0.008 * k, "-", lw=1, color=f"C{k}", label=f"S{sct}" + ("" if sct in good else " (noisy)"))
ax.set_ylabel("norm. flux (offset)"); ax.set_xlabel("phase"); ax.legend(fontsize=8, ncol=5)
for a in axs:
    a.axvline(0.5, color="C3", ls=":", lw=0.8)
fig.tight_layout(); fig.savefig(out("Folded_LC.png"), dpi=130); plt.close(fig)

# summary numbers
oot = (np.abs(c) > 0.06) & (np.abs(c - 0.5) > 0.06)
core = np.abs(c) < 0.01; sec = np.abs(c - 0.5) < 0.01
print("peaks:", [round(p, 5) for p in peaks])
print("primary depth (binned core) %.2f%%" % ((np.median(m[oot]) - m[core].mean()) * 100))
print("phase-0.5 level minus OOT median %.2f ppt" % ((m[sec].mean() - np.median(m[oot])) * 1e3))
print("OOT peak-to-peak (binned) %.2f%%" % ((m[oot].max() - m[oot].min()) * 100))
for sct, d in secs.items():
    cc, mm, ss, nn = binned(phase(d["t"]), d["f"], 100)
    o = (np.abs(cc) > 0.06) & (np.abs(cc - 0.5) > 0.06); cr = np.abs(cc) < 0.01
    print(f"S{sct} primary depth {(np.median(mm[o]) - mm[cr].mean()) * 100:.2f}%")
