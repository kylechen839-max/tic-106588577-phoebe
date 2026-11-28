"""MIPS 24 aperture photometry on an SHA pbcd mosaic with a rotated CD matrix (astropy.wcs is blocked here).
python outputs/spitzer_mips24/sha_point.py FILE RA DEC W4_mJy OUTPNG"""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits
f, ra, dec, w4, png = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
h = fits.getheader(f); img = fits.getdata(f).astype(float)
a0, d0 = np.radians(h["CRVAL1"]), np.radians(h["CRVAL2"]); a, d = np.radians(ra), np.radians(dec)
c = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
X = np.degrees(np.cos(d) * np.sin(a - a0) / c); Y = np.degrees((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / c)
CD = np.array([[h["CD1_1"], h["CD1_2"]], [h["CD2_1"], h["CD2_2"]]])
dx, dy = np.linalg.solve(CD, [X, Y]); x0, y0 = h["CRPIX1"] - 1 + dx, h["CRPIX2"] - 1 + dy
pix = np.sqrt(abs(np.linalg.det(CD))) * 3600; sr = (pix / 206264.806) ** 2
print(f"pixel ({x0:.1f},{y0:.1f}) of {img.shape[::-1]}, {pix:.2f}\"/px, value {img[int(round(y0)), int(round(x0))]}")
yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]; r = np.hypot(xx - x0, yy - y0) * pix
sky = img[(r >= 20) & (r < 32)]; bg = np.nanmedian(sky); sig = 1.4826 * np.nanmedian(np.abs(sky - bg))
ap = r < 7; flux = np.nansum(img[ap] - bg) * sr * 1e9 * 1.61; err = sig * np.sqrt(ap.sum()) * sr * 1e9 * 1.61
pk = (np.nanmax(img[r < 5]) - bg) / sig
print(f"MIPS24 r=7\" flux {flux:.2f} ± {err:.2f} mJy, peak S/N {pk:.1f}, sky {bg:.2f}±{sig:.3f} MJy/sr; W4 {w4} mJy -> ratio {flux / w4:.2f}")
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
for ax, half in zip(axs, (60, 150)):
    m = r < half * 1.42
    ys, xs = np.where(m); sub = img[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    lo, hi = np.nanpercentile(sub, [5, 99.7])
    ax.imshow(sub, origin="lower", cmap="magma", vmin=lo, vmax=hi)
    for rr, st in ((6, "-"), (12, ":")):
        ax.add_patch(plt.Circle((x0 - xs.min(), y0 - ys.min()), rr / pix, fill=False, color="cyan" if st == "-" else "w", ls=st))
    ax.set_title(f"MIPS 24 µm (SHA, rotated frame), ~{2 * half}\" field", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(png, dpi=110)
