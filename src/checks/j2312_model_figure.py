"""J2312 W4: figure of the nebula-subtraction model (method B of j2312_psf_extraction.py) and what it leaves.

Top row: W4 data | nebula model (scaled, smoothed W3 template + plane) | star model (PSF x A) | residual.
Bottom row: data minus nebula model (the star + any excess) | radial profile of that against the empirical point
source | W4 flux at the star for each background model vs the binary photosphere.
Output: outputs/j2312_contamination/j2312_model_figure.png
"""
import importlib.util
from pathlib import Path
import numpy as np
from scipy import ndimage
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

spec = importlib.util.spec_from_file_location("m", str(Path(__file__).with_name("j2312_psf_extraction.py")))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
PIX = m.PIX

img4, w4 = m.load("w4"); img3, w3 = m.load("w3")
psf4, _, fw4 = m.build_psf(img4, w4); psf3, _, fw3 = m.build_psf(img3, w3)
hp = int(160 / PIX)
s4, dx, dy = m.cutout(img4, w4, hp); s3, dx3, dy3 = m.cutout(img3, w3, hp)
yy, xx = np.mgrid[:s4.shape[0], :s4.shape[1]]
X, Y = (xx - hp - dx) * PIX, (yy - hp - dy) * PIX
R = np.hypot(X, Y); one = np.ones_like(s4)
P4 = m.placed_psf(psf4, s4.shape, dx, dy); P3 = m.placed_psf(psf3, s3.shape, dx3, dy3)

# a0 reference (catalogue-like) for flux calibration
sky = np.nanmedian(s4[(R >= 50) & (R < 70)])
a0 = np.sum((P4 * (s4 - sky))[R < 55]) / np.sum((P4 * P4)[R < 55])
mjy = m.F_CAT / a0                      # mJy per unit PSF amplitude

# method B, as in j2312_psf_extraction.py
coef3, _ = m.lstsq([P3, one, X, Y], s3, R < 40)
neb3 = s3 - coef3[0] * P3
sm = ndimage.median_filter(neb3, size=int(30 / PIX)); resid = neb3 - sm
n3 = 1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid)))
pts = ndimage.binary_dilation(resid > 5 * n3, iterations=int(8 / PIX)) & (R > 20)
neb3[pts] = sm[pts]
tmpl = ndimage.gaussian_filter(neb3, np.sqrt(fw4 ** 2 - fw3 ** 2) / 2.3548 / PIX)
tmpl -= np.nanmedian(tmpl[(R > 140) & (R < 160)])
c, model = m.lstsq([P4, tmpl, one, X, Y], s4, R < 150)
neb = c[1] * tmpl + c[2] + c[3] * X + c[4] * Y
star = c[0] * P4
res = s4 - neb - star
starplus = s4 - neb
F_B = c[0] * mjy
print(f"method B: {F_B:.2f} mJy ({F_B / m.F_PHOT:.2f}x), s={c[1]:.4f}")

# to mJy per pixel for display, zeroed at the 140-160" sky level
k = mjy
z0 = np.nanmedian(s4[(R > 140) & (R < 160)])
v = np.nanpercentile((s4 - z0)[R < 150] * k, [2, 99.7])
noise = 1.4826 * np.nanmedian(np.abs(res[(R > 30) & (R < 150)])) * k

fig = plt.figure(figsize=(16, 8.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05])
ext = [-160, 160, -160, 160]
top = [((s4 - z0) * k, "W4 data (22 µm)", "magma", v),
       ((neb - z0) * k, f"Nebula model\n(12 µm map, blurred, × {c[1] * 5.277e-5 / 1.998e-6:.2f} in flux, + tilt)", "magma", v),
       (star * k, f"Point source at the star\n{F_B:.1f} mJy", "magma", v),
       (np.where(R < 150, res, np.nan) * k, "Residual: data − nebula − point source", "coolwarm", (-4 * noise, 4 * noise))]
for i, (a, t, cm, lim) in enumerate(top):
    ax = fig.add_subplot(gs[0, i])
    im = ax.imshow(a, origin="lower", cmap=cm, vmin=lim[0], vmax=lim[1], extent=ext)
    ax.set_title(t, fontsize=10); ax.plot(0, 0, "+", color="c", ms=9)
    ax.add_patch(plt.Circle((0, 0), 150, fill=False, ls=":", color="w" if cm == "magma" else "k", lw=0.7))
    ax.set_xticks([-150, 0, 150]); ax.set_yticks([-150, 0, 150]); ax.tick_params(labelsize=8)
    if i == 0: ax.set_ylabel("arcsec")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7)
    cb.set_label("mJy / pixel" if i < 3 else "mJy / pixel (±4σ)", fontsize=7)

ax = fig.add_subplot(gs[1, 0])
z = 60
lim2 = np.nanpercentile(starplus[R < z] * k, [1, 99.8])
ax.imshow(starplus * k, origin="lower", cmap="magma", vmin=lim2[0], vmax=lim2[1], extent=ext)
ax.set_xlim(-z, z); ax.set_ylim(-z, z); ax.plot(0, 0, "+", color="c", ms=9)
ax.set_title("Data − nebula model (zoom)\n= stars + whatever is left", fontsize=10); ax.tick_params(labelsize=8)
ax.add_patch(plt.Circle((0, 0), fw4 / 2, fill=False, color="w", lw=0.8, ls="--"))
ax.text(-z + 4, -z + 6, f"dashed: point-source FWHM {fw4:.0f}″", color="w", fontsize=8)

# radial profile
ax = fig.add_subplot(gs[1, 1:3])
edges = np.arange(0, 62, 2.75); rc = 0.5 * (edges[1:] + edges[:-1])
prof = lambda a: np.array([np.nanmean(a[(R >= e0) & (R < e1)]) for e0, e1 in zip(edges[:-1], edges[1:])])
p_obs = prof(starplus * k)
p_psf = prof(P4 * k)
p_phot = p_psf * (m.F_PHOT / m.F_CAT) * a0          # stars alone (PHOEBE binary photosphere)
p_fit = p_psf * c[0]
err = noise / np.sqrt(np.array([max(((R >= e0) & (R < e1)).sum() / 30, 1) for e0, e1 in zip(edges[:-1], edges[1:])]))
ax.errorbar(rc, p_obs, yerr=err, fmt="o", color="k", ms=4, label="J2312: data − nebula model")
ax.plot(rc, p_fit, "-", color="C1", lw=2, label=f"point source fitted to it ({F_B:.1f} mJy)")
ax.plot(rc, p_phot, "--", color="C0", lw=2, label=f"stars alone, no dust ({m.F_PHOT:.1f} mJy, PHOEBE)")
ax.fill_between(rc, p_phot, p_fit, color="C1", alpha=0.15, label="excess (dust)")
ax.axhline(0, color="0.6", lw=0.6)
ax.set_xlabel("distance from the star (arcsec)"); ax.set_ylabel("mean brightness (mJy / pixel)")
import importlib.util as iu
sp = iu.spec_from_file_location("t", str(Path(__file__).with_name("tierA_w4_shape.py"))); t = iu.module_from_spec(sp); sp.loader.exec_module(t)
fw_raw = t.gfit(s4, hp + dx, hp + dy, 20)[3]; fw_sub = t.gfit(starplus, hp + dx, hp + dy, 20)[3]; fw_psf = t.gfit(P4 * 1e3, hp + dx, hp + dy, 20)[3]
print(f"Gaussian FWHM within 20\": data+plane {fw_raw:.1f}, data-nebula {fw_sub:.1f}, PSF {fw_psf:.1f}")
ax.set_title(f"Radial profile at 22 µm. Gaussian FWHM within 20″: J2312 after nebula removal {fw_sub:.1f}″, "
             f"point source {fw_psf:.1f}″ (before removal {fw_raw:.1f}″)", fontsize=9)
ax.legend(fontsize=8); ax.set_xlim(0, 60)
ax2 = ax.twiny(); ax2.set_xlim(0, 60 * 742 / 1000); ax2.set_xlabel("projected distance at 742 pc (1000 AU)", fontsize=8)
ax2.tick_params(labelsize=8)

# flux comparison
ax = fig.add_subplot(gs[1, 3])
names = ["catalogue\nflat sky", "plane", "quadratic", "12 µm\ntemplate", "size\nseparation"]
vals = [15.36, 14.55, 12.3, F_B, 12.4]
cols = ["0.6", "0.6", "C1", "C1", "C1"]
y = np.arange(len(names))
ax.barh(y, vals, xerr=m.E_CAT, color=cols, alpha=0.85, capsize=3)
ax.axvline(m.F_PHOT, color="C0", ls="--", lw=2); ax.text(m.F_PHOT + 0.2, len(names) - 0.45, "stars\nalone", color="C0", fontsize=8, va="top")
for yi, vv in zip(y, vals):
    ax.text(vv + 1.1, yi, f"{vv / m.F_PHOT:.1f}×, {(vv - m.F_PHOT) / m.E_CAT:.1f}σ", va="center", fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8); ax.invert_yaxis()
ax.set_xlim(0, 24); ax.set_xlabel("W4 flux at the star (mJy)")
ax.set_title("W4 at the star, by background model\n(orange: nebula modelled)", fontsize=10)

fig.suptitle("J2312 at 22 µm: nebula model subtracted, point source fitted, and what is left", fontsize=12)
fig.tight_layout()
out = m.D / "j2312_model_figure.png"
fig.savefig(out, dpi=150)
print(out)
