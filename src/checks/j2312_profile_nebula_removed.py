"""J2312 W4 profile and compact/extended split AFTER removing the 12 um-map nebula model (method B).

Re-implements method B of j2312_psf_extraction.py without astropy.wcs (blocked on Kyle's PC):
native Atlas tiles 3476p530_ac51 (RA---SIN, no rotation, 1.375"/px), empirical PSFs from isolated
field stars, W3 star (PSF + plane < 40") removed, other compact W3 sources masked and filled,
template blurred to W4 resolution, then W4 = A*PSF4 + s*template + plane fitted within 150".
The nebula term (s*template + plane) is subtracted and J2312's profile compared with field stars;
the compact + extended fit is repeated on the nebula-subtracted image.
python src/checks/j2312_profile_nebula_removed.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import ndimage

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "outputs" / "j2312_contamination"
TILE, RA, DEC, PIX = "3476p530_ac51", 348.0058494, 53.3413022, 1.375
H = 45            # PSF stamp half-size (px)
BIG = int(160 / PIX)
PHOT, ECAT = 4.31, 0.88


def load(b):
    with fits.open(D / "native" / f"{TILE}-{b}-int-3.fits") as h:
        return h[0].data.astype(float), h[0].header


def sin_pix(h, ra, dec):
    a0, d0, a, d = map(np.radians, (h["CRVAL1"], h["CRVAL2"], ra, dec))
    x = np.degrees(np.cos(d) * np.sin(a - a0))
    y = np.degrees(np.sin(d) * np.cos(d0) - np.cos(d) * np.sin(d0) * np.cos(a - a0))
    return h["CRPIX1"] - 1 + x / h["CDELT1"], h["CRPIX2"] - 1 + y / h["CDELT2"]


def plane_sub(st, cx, cy, r_in=30, r_out=60):
    yy, xx = np.mgrid[:st.shape[0], :st.shape[1]]
    r = np.hypot(xx - cx, yy - cy) * PIX
    mm = (r > r_in) & (r < r_out)
    c, *_ = np.linalg.lstsq(np.c_[np.ones(mm.sum()), xx[mm], yy[mm]], st[mm], rcond=None)
    return st - (c[0] + c[1] * xx + c[2] * yy)


def build_psf(img, x0, y0, amp_min):
    bg = ndimage.zoom(ndimage.median_filter(img[::4, ::4], size=15), 4, order=1)[:img.shape[0], :img.shape[1]]
    res = img - bg
    noise = 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res)))
    mx = ndimage.maximum_filter(res, size=31)
    sat = 0.98 * np.nanmax(img)
    pk = np.argwhere((res == mx) & (res > max(20 * noise, amp_min)))
    st_all = []
    for y, x in pk:
        if not (H + 5 < y < img.shape[0] - H - 5 and H + 5 < x < img.shape[1] - H - 5):
            continue
        if np.hypot(x - x0, y - y0) * PIX < 300 or img[y, x] >= sat:
            continue
        st = plane_sub(img[y - H:y + H + 1, x - H:x + H + 1].copy(), H, H)
        yy, xx = np.mgrid[:st.shape[0], :st.shape[1]]
        rr = np.hypot(xx - H, yy - H) * PIX
        if np.nanmax(ndimage.gaussian_filter(st, 2)[rr > 25]) > 0.12 * st[H, H]:
            continue
        w = np.clip(st, 0, None) * (rr < 5)
        cy, cx = (w * yy).sum() / w.sum(), (w * xx).sum() / w.sum()
        st = ndimage.shift(st, (H - cy, H - cx), order=3, mode="nearest")
        st_all.append(st / np.nansum(st[rr < 35]))
    psf = np.nanmedian(st_all[:60], axis=0)
    return psf / np.nansum(psf[rr < 35]), np.array(st_all[:60])


def radial(a, cx, cy, edges):
    yy, xx = np.mgrid[:a.shape[0], :a.shape[1]]
    r = np.hypot(xx - cx, yy - cy) * PIX
    return np.array([np.nanmean(a[(r >= p) & (r < q)]) for p, q in zip(edges[:-1], edges[1:])])


def fwhm(p, rc):
    p = p / p[0]
    return 2 * np.interp(0.5, p[::-1], rc[::-1])


def place(psf, shape, cx, cy):
    out = np.zeros(shape)
    ix, iy = int(round(cx)), int(round(cy))
    sh = ndimage.shift(psf, (cy - iy, cx - ix), order=3, mode="constant")
    out[iy - H:iy + H + 1, ix - H:ix + H + 1] = sh
    return out


edges = np.arange(0, 41, PIX)
rc = 0.5 * (edges[1:] + edges[:-1])
img4, h4 = load("w4"); img3, h3 = load("w3")
x4, y4 = sin_pix(h4, RA, DEC); x3, y3 = sin_pix(h3, RA, DEC)
assert abs(x4 - x3) < 0.01 and abs(y4 - y3) < 0.01  # same pixel grid
ix, iy = int(round(x4)), int(round(y4))
s4 = img4[iy - BIG:iy + BIG + 1, ix - BIG:ix + BIG + 1].copy()
s3 = img3[iy - BIG:iy + BIG + 1, ix - BIG:ix + BIG + 1].copy()
cx, cy = x4 - ix + BIG, y4 - iy + BIG
yy, xx = np.mgrid[:s4.shape[0], :s4.shape[1]]
R = np.hypot(xx - cx, yy - cy) * PIX
X, Y, one = (xx - cx) * PIX / 100, (yy - cy) * PIX / 100, np.ones_like(s4)

amp4 = plane_sub(img4[iy - H:iy + H + 1, ix - H:ix + H + 1].copy(), H, H)[H, H]
amp3 = plane_sub(img3[iy - H:iy + H + 1, ix - H:ix + H + 1].copy(), H, H)[H, H]
psf4, stars4 = build_psf(img4, x4, y4, amp4 / 3)
psf3, stars3 = build_psf(img3, x3, y3, amp3 / 3)
fw4 = fwhm(radial(psf4, H, H, edges), rc); fw3 = fwhm(radial(psf3, H, H, edges), rc)
print(f"PSFs: W4 {len(stars4)} stars, FWHM {fw4:.1f}\"; W3 {len(stars3)} stars, FWHM {fw3:.1f}\"")

# W3: remove J2312 (PSF + plane within 40"), mask other compact sources, blur to W4
P3 = place(psf3, s3.shape, cx, cy)
mk = R < 40
c3, *_ = np.linalg.lstsq(np.c_[P3[mk], one[mk], X[mk], Y[mk]], s3[mk], rcond=None)
neb3 = s3 - c3[0] * P3
sm = ndimage.median_filter(neb3, size=int(30 / PIX))
resid = neb3 - sm
n3 = 1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid)))
pts = ndimage.binary_dilation(resid > 5 * n3, iterations=int(8 / PIX)) & (R > 20)
neb3[pts] = sm[pts]
tmpl = ndimage.gaussian_filter(neb3, np.sqrt(max(fw4 ** 2 - fw3 ** 2, 1)) / 2.3548 / PIX)
tmpl -= np.nanmedian(tmpl[(R > 140) & (R < 160)])

# W4: A*PSF4 + s*template + plane within 150"
P4 = place(psf4, s4.shape, cx, cy)
mf = R < 150
cB, *_ = np.linalg.lstsq(np.c_[P4[mf], tmpl[mf], one[mf], X[mf], Y[mf]], s4[mf], rcond=None)
neb = cB[1] * tmpl + cB[2] + cB[3] * X + cB[4] * Y
c0, *_ = np.linalg.lstsq(np.c_[P4[R < 55]], (s4 - np.nanmedian(s4[(R >= 50) & (R < 70)]))[R < 55], rcond=None)
FB = 15.36 * cB[0] / c0[0]     # scale to the AllWISE catalogue flux via a catalogue-like fit (flat sky 50-70")
print(f"method B (re-implemented): W4 at star {FB:.2f} mJy = {FB / PHOT:.2f}x photosphere, "
      f"{(FB - PHOT) / ECAT:.1f} sigma; template scale s = {cB[1]:.4f}; {pts.sum()} W3 px masked")

star_img = s4 - neb                          # nebula removed (star + anything point-like or unmodelled)
plane_img = plane_sub(s4.copy(), cx, cy)     # the earlier, plane-only version

p_neb = radial(star_img, cx, cy, edges)
p_pl = radial(plane_img, cx, cy, edges)
p_psf = radial(psf4, H, H, edges)
fstars = [fwhm(radial(s, H, H, edges), rc) for s in stars4]
lo, hi = np.nanpercentile([radial(s, H, H, edges) / radial(s, H, H, edges)[0] for s in stars4], [16, 84], axis=0)
print(f"FWHM: field stars median {np.median(fstars):.1f}\" (95% {np.percentile(fstars, 95):.1f}\"); "
      f"J2312 plane-only {fwhm(p_pl, rc):.1f}\"; J2312 nebula removed {fwhm(p_neb, rc):.1f}\"")

# compact + extended split on the nebula-removed image (photosphere forced compact)
st = star_img[int(round(cy)) - H:int(round(cy)) + H + 1, int(round(cx)) - H:int(round(cx)) + H + 1]
sx, sy = cx - int(round(cx)) + H, cy - int(round(cy)) + H
st = ndimage.shift(st, (H - sy, H - sx), order=3, mode="nearest")
yyp, xxp = np.mgrid[:psf4.shape[0], :psf4.shape[1]]
rr = np.hypot(xxp - H, yyp - H) * PIX
m35 = rr < 35
data = st / np.nansum(st[m35])
fc_min = PHOT / FB
grid = []
for g in (6, 8, 10, 13, 16, 20, 25, 30, 40):
    ext = ndimage.gaussian_filter(psf4, g / 2.3548 / PIX); ext /= np.nansum(ext[m35])
    for fc in np.arange(fc_min, 1.0001, 0.01):
        mod = fc * psf4 + (1 - fc) * ext
        grid.append((np.sqrt(np.nanmean((data[m35] - mod[m35]) ** 2)), g, fc))
grid.sort()
r_pt = np.sqrt(np.nanmean((data[m35] - psf4[m35]) ** 2))
okk = [x for x in grid if x[0] <= 1.10 * grid[0][0]]
print(f"split: point source alone rms {r_pt:.2e}; best rms {grid[0][0]:.2e} at halo {grid[0][1]}\", compact fraction {grid[0][2]:.2f}")
for lab, fc in (("best", grid[0][2]), ("least compact within 10%", min(x[2] for x in okk)), ("most compact within 10%", max(x[2] for x in okk))):
    Fc = fc * FB
    print(f"  {lab}: compact {Fc:.1f} mJy = photosphere + {Fc - PHOT:.1f} mJy ({Fc / PHOT:.2f}x), extended {FB - Fc:.1f} mJy")

# figure
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw=dict(width_ratios=[1, 1, 1.3]))
ext_ = [-BIG * PIX, BIG * PIX, -BIG * PIX, BIG * PIX]
v = np.nanpercentile(s4, [3, 99.5])
for ax, a, t in ((axs[0], s4 - np.nanmedian(s4[(R > 140) & (R < 160)]), "W4 data"),
                 (axs[1], star_img - np.nanmedian(star_img[(R > 140) & (R < 160)]), "W4 − nebula model (12 µm map)")):
    ax.imshow(a, origin="lower", cmap="magma", vmin=0, vmax=v[1] - v[0], extent=ext_)
    ax.set_xlim(80, -80); ax.set_ylim(-80, 80); ax.set_title(t, fontsize=10, loc="left")
    ax.set_xlabel("← E   arcsec   W →"); ax.plot(0, 0, "+", color="cyan", ms=12)
ax = axs[2]
ax.fill_between(rc, lo, hi, color="#a9a8a2", alpha=0.35, lw=0, label=f"field stars 16–84% (n={len(stars4)})")
ax.plot(rc, p_psf / p_psf[0], color="#52514e", lw=2, label=f"point source (PSF), FWHM {np.median(fstars):.1f}″")
ax.plot(rc, p_pl / p_pl[0], color="#86b6ef", lw=2, ls="--", label=f"J2312, plane background only, FWHM {fwhm(p_pl, rc):.1f}″")
ax.plot(rc, p_neb / p_neb[0], color="#2a78d6", lw=2.5, marker="o", ms=4, label=f"J2312, nebula model removed, FWHM {fwhm(p_neb, rc):.1f}″")
ax.axhline(0.5, color="#e4e3df", lw=1, zorder=0); ax.set_xlim(0, 40); ax.set_ylim(-0.05, 1.05)
ax.set_xlabel("radius (arcsec)"); ax.set_ylabel("brightness / centre")
ax.set_title("22 µm radial profile", fontsize=10, loc="left"); ax.legend(fontsize=8, frameon=False)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.suptitle("J2312 at 22 µm after removing the 12 µm-map nebula model", x=0.01, ha="left", fontsize=13, fontweight="bold")
fig.tight_layout(); fig.savefig(D / "j2312_profile_nebula_removed.png", dpi=150)
