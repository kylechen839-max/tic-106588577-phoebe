"""J2312 W4 radial profile with the nebula removed, in the layout of the PC's j2312_w4_profile.png (section 5 of the
figure report), so the two can be compared directly.

Curves (all peak-normalized, azimuthal means in 1.375" annuli):
  - J2312 with only a tilted plane removed (fitted 30-60"): what j2312_w4_profile.png shows
  - J2312 minus the getsf-style nebula model (j2312_getsf_figure.py, W = 120")
  - J2312 minus the 12 um-map nebula model (method B, j2312_model_figure.py)
  - the empirical PSF (median of isolated field stars in the tile), and the PSF blurred by 6" and 10" Gaussians
FWHM = twice the radius where each curve crosses 0.5 (the same read-off as j2312_w4_profile.png).
Output: outputs/j2312_contamination/report/fig5_profile_nebsub.png
"""
import runpy
from pathlib import Path
import numpy as np
from scipy import ndimage
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).parent
b = runpy.run_path(str(here / "j2312_model_figure.py"))
g = runpy.run_path(str(here / "j2312_getsf_figure.py"))
m = b["m"]; PIX = m.PIX
R, X, Y = b["R"], b["X"], b["Y"]
s4, P4 = b["s4"], b["P4"]
one = np.ones_like(s4)

ring = (R >= 30) & (R < 60)
cp, _ = m.lstsq([one, X, Y], s4, ring)
plane_only = s4 - (cp[0] + cp[1] * X + cp[2] * Y)
nebB = s4 - b["neb"]
nebG = g["s4"] - g["neb"]

edges = np.arange(0, 41.25 + 1e-6, PIX); rc = 0.5 * (edges[1:] + edges[:-1])
def prof(a, RR):
    p = np.array([np.nanmean(a[(RR >= e0) & (RR < e1)]) for e0, e1 in zip(edges[:-1], edges[1:])])
    return p / p[0]
def fwhm(p):
    i = np.argmax(p < 0.5)
    return 2 * (rc[i - 1] + (p[i - 1] - 0.5) / (p[i - 1] - p[i]) * (rc[i] - rc[i - 1]))

curves = {"psf": prof(P4, R)}
for fe in (6, 10):
    curves[f"psf{fe}"] = prof(ndimage.gaussian_filter(P4, fe / 2.3548 / PIX), R)
curves["plane"] = prof(plane_only, R)
curves["getsf"] = prof(nebG, g["R"])
curves["B"] = prof(nebB, R)
fw = {k: fwhm(v) for k, v in curves.items()}
for k, v in fw.items():
    print(f"{k:6s} FWHM {v:.1f}\"")

fig, ax = plt.subplots(figsize=(8.8, 5.4))
ax.plot(rc, curves["psf"], "-", color="0.25", lw=2.2, label=f"point source (median field star), FWHM {fw['psf']:.1f}\"")
ax.plot(rc, curves["psf6"], ":", color="C1", lw=1.6, label=f"point source blurred by 6\" (FWHM {fw['psf6']:.1f}\")")
ax.plot(rc, curves["psf10"], "--", color="C1", lw=1.6, label=f"point source blurred by 10\" (FWHM {fw['psf10']:.1f}\")")
ax.plot(rc, curves["plane"], "-", color="0.6", lw=1.6, label=f"J2312, flat background only (as in section 5), FWHM {fw['plane']:.1f}\"")
ax.plot(rc, curves["getsf"], "s-", color="C3", ms=3.5, lw=1.4, label=f"J2312 − getsf-style nebula, FWHM {fw['getsf']:.1f}\"")
ax.plot(rc, curves["B"], "o-", color="C0", ms=4, lw=2.2, label=f"J2312 − 12 µm-map nebula, FWHM {fw['B']:.1f}\"")
ax.axhline(0.5, color="0.85", lw=1); ax.axhline(0, color="0.7", lw=0.6)
ax.set_xlim(0, 40); ax.set_xlabel("radius from the star (arcsec)"); ax.set_ylabel("brightness / peak")
ax.set_title("J2312 at 22 µm with the nebula removed is still broader than a point source", fontsize=11)
ax.legend(fontsize=8)
fig.tight_layout()
out = m.D / "report" / "fig5_profile_nebsub.png"
fig.savefig(out, dpi=140)
print(out)
