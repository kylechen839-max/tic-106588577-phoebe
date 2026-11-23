"""J2312 W4: the getsf-style (W4-only multi-scale) model as a figure, same layout as j2312_model_figure.py.

Model (j2312_multiscale.py, background window W = 120"): image = smooth background + extended/filament structures
+ point sources. The nebula model is background + extended structures; J2312's flux is a PSF fit on the
point-source image, corrected by the field-star flux recovery of the same pipeline.
Top row: W4 data | nebula model | J2312 point source | residual (data - nebula - J2312 point source).
Bottom row: data - nebula model (zoom) | radial profile | W4 flux at the star vs the stars alone.
Output: outputs/j2312_contamination/j2312_getsf_figure.png
"""
import importlib.util
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt


def mod(name, f):
    s = importlib.util.spec_from_file_location(name, str(Path(__file__).with_name(f)))
    x = importlib.util.module_from_spec(s); s.loader.exec_module(x); return x


m = mod("m", "j2312_psf_extraction.py")
ms = mod("ms", "j2312_multiscale.py")
t = mod("t", "tierA_w4_shape.py")
PIX = m.PIX
W = 120

img, w = m.load("w4")
psf, _, fw = m.build_psf(img, w)
x0, y0 = (float(v) for v in w.world_to_pixel(m.C))
xi, yi = int(round(x0)), int(round(y0))
B0 = ms.BOX
sub = img[yi - B0: yi + B0 + 1, xi - B0: xi + B0 + 1].copy()
tx, ty = B0 + x0 - xi, B0 + y0 - yi

bg = ms.bg_iter(sub, W)
src, fil = ms.decompose(sub - bg)
a_src = ms.psf_amp(src, psf, tx, ty, const=False)
stars = ms.field_stars(sub, tx, ty)
rec = float(np.median([ms.psf_amp(src, psf, x, y, const=False) / ms.psf_amp(sub, psf, x, y) for x, y in stars]))

# 320" display cut centred on J2312
hp = int(160 / PIX)
cut = lambda a: a[B0 - hp: B0 + hp + 1, B0 - hp: B0 + hp + 1]
dx, dy = tx - B0, ty - B0
s4, nb, fl = cut(sub), cut(bg), cut(fil)
yy, xx = np.mgrid[:s4.shape[0], :s4.shape[1]]
X, Y = (xx - hp - dx) * PIX, (yy - hp - dy) * PIX
R = np.hypot(X, Y)
P4 = m.placed_psf(psf, s4.shape, dx, dy)

# flux scale: catalogue-like a0 fit -> AllWISE W4 flux
sky = np.nanmedian(s4[(R >= 50) & (R < 70)])
a0 = np.sum((P4 * (s4 - sky))[R < 55]) / np.sum((P4 * P4)[R < 55])
k = m.F_CAT / a0
A = a_src / rec
F = A * k
sig = (F - m.F_PHOT) / m.E_CAT
print(f"getsf-style W={W}\": source-image amplitude {a_src / a0:.3f} a0, field recovery {rec:.3f}, F = {F:.2f} mJy, "
      f"{F / m.F_PHOT:.2f}x photosphere, {sig:.1f} sigma")

neb = nb + fl
star = A * P4
res = s4 - neb - star
starplus = s4 - neb
z0 = np.nanmedian(nb[(R > 140) & (R < 160)])
v = np.nanpercentile((s4 - z0)[R < 150] * k, [2, 99.7])
noise = 1.4826 * np.nanmedian(np.abs(res[(R > 30) & (R < 150)] - np.nanmedian(res[(R > 30) & (R < 150)]))) * k

fig = plt.figure(figsize=(16, 8.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05])
ext = [-160, 160, -160, 160]
top = [((s4 - z0) * k, "W4 data (22 µm)", "magma", v),
       ((neb - z0) * k, f"Nebula model (getsf-style, 22 µm only)\nsmooth background + extended structure", "magma", v),
       (star * k, f"Point source at the star\n{F:.1f} mJy", "magma", v),
       (np.where(R < 150, res, np.nan) * k, "Residual: data − nebula − point source\n(other spots are field stars)", "coolwarm",
        (-4 * noise, 4 * noise))]
for i, (a, ttl, cm, lim) in enumerate(top):
    ax = fig.add_subplot(gs[0, i])
    im = ax.imshow(a, origin="lower", cmap=cm, vmin=lim[0], vmax=lim[1], extent=ext)
    ax.set_title(ttl, fontsize=10); ax.plot(0, 0, "+", color="c", ms=9)
    ax.add_patch(plt.Circle((0, 0), 150, fill=False, ls=":", color="w" if cm == "magma" else "k", lw=0.7))
    ax.set_xticks([-150, 0, 150]); ax.set_yticks([-150, 0, 150]); ax.tick_params(labelsize=8)
    if i == 0: ax.set_ylabel("arcsec")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=7)
    cb.set_label("mJy / pixel" if i < 3 else "mJy / pixel (±4σ)", fontsize=7)

ax = fig.add_subplot(gs[1, 0])
zz = 60
lim2 = np.nanpercentile(starplus[R < zz] * k, [1, 99.8])
ax.imshow(starplus * k, origin="lower", cmap="magma", vmin=lim2[0], vmax=lim2[1], extent=ext)
ax.set_xlim(-zz, zz); ax.set_ylim(-zz, zz); ax.plot(0, 0, "+", color="c", ms=9)
ax.set_title("Data − nebula model (zoom)\n= stars + whatever is left", fontsize=10); ax.tick_params(labelsize=8)
ax.add_patch(plt.Circle((0, 0), fw / 2, fill=False, color="w", lw=0.8, ls="--"))
ax.text(-zz + 4, -zz + 6, f"dashed: point-source FWHM {fw:.0f}″", color="w", fontsize=8)

ax = fig.add_subplot(gs[1, 1:3])
edges = np.arange(0, 62, 2.75); rc = 0.5 * (edges[1:] + edges[:-1])
prof = lambda a: np.array([np.nanmean(a[(R >= e0) & (R < e1)]) for e0, e1 in zip(edges[:-1], edges[1:])])
p_obs = prof(starplus * k); p_fit = prof(star * k); p_phot = p_fit * m.F_PHOT / F
err = noise / np.sqrt(np.array([max(((R >= e0) & (R < e1)).sum() / 30, 1) for e0, e1 in zip(edges[:-1], edges[1:])]))
fw_sub = t.gfit(starplus, hp + dx, hp + dy, 20)[3]; fw_psf = t.gfit(P4 * 1e3, hp + dx, hp + dy, 20)[3]
ax.errorbar(rc, p_obs, yerr=err, fmt="o", color="k", ms=4, label="J2312: data − nebula model")
ax.plot(rc, p_fit, "-", color="C1", lw=2, label=f"point source fitted ({F:.1f} mJy)")
ax.plot(rc, p_phot, "--", color="C0", lw=2, label=f"stars alone, no dust ({m.F_PHOT:.1f} mJy, PHOEBE)")
ax.fill_between(rc, p_phot, p_fit, color="C1", alpha=0.15, label="excess (dust)")
ax.axhline(0, color="0.6", lw=0.6)
ax.set_xlabel("distance from the star (arcsec)"); ax.set_ylabel("mean brightness (mJy / pixel)")
ax.set_title(f"Radial profile at 22 µm. Gaussian FWHM within 20″: J2312 after nebula removal {fw_sub:.1f}″, "
             f"point source {fw_psf:.1f}″", fontsize=9)
ax.legend(fontsize=8); ax.set_xlim(0, 60)
ax2 = ax.twiny(); ax2.set_xlim(0, 60 * 742 / 1000); ax2.set_xlabel("projected distance at 742 pc (1000 AU)", fontsize=8)
ax2.tick_params(labelsize=8)
print(f"FWHM after getsf nebula removal {fw_sub:.1f}, PSF {fw_psf:.1f}")

ax = fig.add_subplot(gs[1, 3])
names = ["stars alone\n(PHOEBE)", "getsf-style\npoint source"]
vals = [m.F_PHOT, F]
ax.barh([0, 1], vals, xerr=[0, m.E_CAT], color=["C0", "C1"], alpha=0.85, capsize=4)
ax.text(F + 1.2, 1, f"excess {sig:.1f}σ", va="center", fontsize=9)
ax.set_yticks([0, 1]); ax.set_yticklabels(names, fontsize=9); ax.invert_yaxis()
ax.set_xlim(0, 24); ax.set_xlabel("W4 flux at the star (mJy)")
ax.set_title("W4 at the star", fontsize=10)

fig.suptitle("J2312 at 22 µm: getsf-style nebula model subtracted, point source fitted, and what is left", fontsize=12)
fig.tight_layout()
out = m.D / "j2312_getsf_figure.png"
fig.savefig(out, dpi=150)
print(out)
