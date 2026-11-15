"""Show where the 22 um ridge is around J2312: annotated W4/W3 images plus an E-W brightness cut.
Native tile 3476p530_ac51 (downloaded by j2312_w4_profile_plot.py). No astropy.wcs (blocked on this PC).
python src/checks/j2312_ridge_plot.py"""
import sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import ndimage
sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[2]; D = ROOT / "outputs" / "j2312_contamination"
TILE, RA, DEC, PIX = "3476p530_ac51", 348.0058494, 53.3413022, 1.375
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"

def load(b):
    with fits.open(D / "native" / f"{TILE}-{b}-int-3.fits") as h:
        return h[0].data.astype(float), h[0].header

def sin_pix(h, ra, dec):
    a0, d0, a, d = map(np.radians, (h["CRVAL1"], h["CRVAL2"], ra, dec))
    x = np.degrees(np.cos(d) * np.sin(a - a0)); y = np.degrees(np.sin(d) * np.cos(d0) - np.cos(d) * np.sin(d0) * np.cos(a - a0))
    return h["CRPIX1"] - 1 + x / h["CDELT1"], h["CRPIX2"] - 1 + y / h["CDELT2"]

half = 240  # arcsec
fig = plt.figure(figsize=(15, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15])
for k, (b, lab) in enumerate((("w4", "WISE 22 µm (W4)"), ("w3", "WISE 12 µm (W3)"))):
    img, h = load(b); x0, y0 = sin_pix(h, RA, DEC); n = int(half / PIX)
    sub = img[int(y0) - n:int(y0) + n + 1, int(x0) - n:int(x0) + n + 1]
    sm = ndimage.gaussian_filter(sub, 6 / PIX)
    ext = [half, -half, -half, half]  # East to the left
    ax = fig.add_subplot(gs[k])
    lo, hi = np.nanpercentile(sm, [10, 99.0 if b == "w4" else 98.5])
    ax.imshow(sm, origin="lower", cmap="magma", vmin=lo, vmax=hi, extent=ext)
    neb = ndimage.gaussian_filter(sub, 20 / PIX)
    lv = np.nanpercentile(neb, [70, 85, 93])
    ax.contour(neb, levels=lv, colors="w", linewidths=0.6, alpha=0.7, extent=ext, origin="lower")
    ax.plot(0, 0, "+", color="cyan", ms=14, mew=1.5)
    if b == "w4":
        ax.annotate("", xy=(10, -200), xytext=(10, 100), arrowprops=dict(arrowstyle="<->", color="w", lw=1.6))
        ax.text(25, -150, "N–S ridge\n(~5′ long)", color="w", fontsize=10, ha="left")
    ax.set_xlim(half, -half); ax.set_ylim(-half, half)
    ax.set_xlabel("← East   offset (arcsec)   West →"); ax.set_ylabel("← South   offset (arcsec)   North →" if k == 0 else "")
    ax.set_title(f"{lab}: 8′ × 8′, smoothed 6″; contours = large-scale nebula", fontsize=9.5, loc="left")
    if b == "w4":
        w4, cx, cy = sub, n, n
# E-W brightness cuts across the ridge in W4 at several declination offsets
ax = fig.add_subplot(gs[2])
sm4 = ndimage.gaussian_filter(w4, 8 / PIX)
xs = (np.arange(sm4.shape[1]) - cx) * PIX
bg = np.nanpercentile(sm4, 20)
for dy, col in ((60, "#86b6ef"), (-60, "#2a78d6"), (-120, "#1c5cab"), (-180, "#0d366b")):
    row = int(cy + dy / PIX)
    strip = np.nanmean(sm4[row - 5:row + 6, :], axis=0) - bg
    ax.plot(xs, strip, color=col, lw=1.8, label=f"{abs(dy)}″ {'north' if dy > 0 else 'south'} of the star")
ax.axvline(0, color=INK2, lw=0.8, ls=":"); ax.text(-4, -0.06, "star's RA", fontsize=8, color=INK2)
ax.set_xlim(half, -half); ax.set_xlabel("← East   offset (arcsec)   West →"); ax.set_ylabel("W4 brightness above background (DN)")
ax.set_title("East–west cuts across the ridge (avoiding the star itself)", fontsize=9.5, loc="left")
ax.legend(fontsize=8, frameon=False, loc="upper right"); [ax.spines[s].set_visible(False) for s in ("top", "right")]
fig.suptitle("Where the 22 µm ridge is: a north–south band of dust running through J2312 (cyan +)", x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
fig.tight_layout(); fig.savefig(D / "j2312_ridge.png", dpi=150, facecolor=SURF)
for dy in (60, -60, -120, -180):
    row = int(cy + dy / PIX); strip = np.nanmean(sm4[row - 5:row + 6, :], axis=0) - bg
    print(dy, "peak at x =", round(xs[np.nanargmax(strip[(np.abs(xs) < 120)]) + np.argmax(np.abs(xs) < 120)], 1), "arcsec")
