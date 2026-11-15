"""J2312: W3/W4 radial profiles vs field stars, to show the 22 um broadening in one figure.

Native AllWISE Atlas tile 3476p530_ac51 (IRSA AWS mirror, 1.375"/px, RA---SIN). The local background
(plane fitted at 30-60") is removed around each source; profiles are normalised to their own peak.
Field stars: isolated, unsaturated, >= 5' from J2312, peak >= 1/3 of J2312's (as in tierA_w4_shape.py).
Model: empirical PSF convolved with a circular Gaussian of FWHM g, i.e. a point source plus an
extended component of that size. astropy.wcs is not used (blocked on Kyle's PC).

  python src/checks/j2312_w4_profile_plot.py
"""
import sys
import urllib.request
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
NAT = D / "native"
NAT.mkdir(exist_ok=True)
TILE, RA, DEC, PIX, HALF = "3476p530_ac51", 348.0058494, 53.3413022, 1.375, 45
URL = "https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/allwise/images/p3am_cdd/34/3476/{0}/{0}-{1}-int-3.fits"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#a9a8a2"


def load(band):
    f = NAT / f"{TILE}-{band}-int-3.fits"
    if not f.exists():
        urllib.request.urlretrieve(URL.format(TILE, band), f)
    with fits.open(f) as h:
        return h[0].data.astype(float), h[0].header


def sin_pix(h, ra, dec):
    a0, d0, a, d = map(np.radians, (h["CRVAL1"], h["CRVAL2"], ra, dec))
    x = np.degrees(np.cos(d) * np.sin(a - a0))
    y = np.degrees(np.sin(d) * np.cos(d0) - np.cos(d) * np.sin(d0) * np.cos(a - a0))
    return h["CRPIX1"] - 1 + x / h["CDELT1"], h["CRPIX2"] - 1 + y / h["CDELT2"]


def stamp(img, x, y):
    ix, iy = int(round(x)), int(round(y))
    st = img[iy - HALF:iy + HALF + 1, ix - HALF:ix + HALF + 1].copy()
    yy, xx = np.mgrid[:st.shape[0], :st.shape[1]]
    cx, cy = x - ix + HALF, y - iy + HALF
    r = np.hypot(xx - cx, yy - cy) * PIX
    m = (r > 30) & (r < 60)
    A = np.c_[np.ones(m.sum()), xx[m], yy[m]]
    c, *_ = np.linalg.lstsq(A, st[m], rcond=None)
    st -= c[0] + c[1] * xx + c[2] * yy
    # recentre on the light-weighted centroid within 5"
    w = np.clip(st, 0, None) * (r < 5)
    cx, cy = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    st = ndimage.shift(st, (HALF - cy, HALF - cx), order=3, mode="nearest")
    return st


def profile(st, edges):
    yy, xx = np.mgrid[:st.shape[0], :st.shape[1]]
    r = np.hypot(xx - HALF, yy - HALF) * PIX
    p = np.array([np.nanmean(st[(r >= a) & (r < b)]) for a, b in zip(edges[:-1], edges[1:])])
    return p / p[0]


def fwhm(p, rc):
    return 2 * np.interp(0.5, p[::-1], rc[::-1])


edges = np.arange(0, 46, PIX)
rc = 0.5 * (edges[1:] + edges[:-1])
out = {}
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, band, label in zip(axs, ("w3", "w4"), ("W3 (12 µm)", "W4 (22 µm)")):
    img, h = load(band)
    x0, y0 = sin_pix(h, RA, DEC)
    st0 = stamp(img, x0, y0)
    pj = profile(st0, edges)
    amp = st0[HALF, HALF]
    # field stars
    bg = ndimage.median_filter(img[::4, ::4], size=15)
    bg = ndimage.zoom(bg, 4, order=1)[:img.shape[0], :img.shape[1]]
    res = img - bg
    noise = 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res)))
    mx = ndimage.maximum_filter(res, size=31)
    sat = 0.98 * np.nanmax(img)
    pk = np.argwhere((res == mx) & (res > max(20 * noise, amp / 3)))
    profs = []
    for y, x in pk:
        if not (HALF + 5 < y < img.shape[0] - HALF - 5 and HALF + 5 < x < img.shape[1] - HALF - 5):
            continue
        if np.hypot(x - x0, y - y0) * PIX < 300 or img[y, x] >= sat:
            continue
        st = stamp(img, float(x), float(y))
        yy, xx = np.mgrid[:st.shape[0], :st.shape[1]]
        rr = np.hypot(xx - HALF, yy - HALF) * PIX
        if np.nanmax(ndimage.gaussian_filter(st, 2)[rr > 25]) > 0.12 * st[HALF, HALF]:
            continue
        profs.append(profile(st, edges))
    profs = np.array(profs[:60])
    psf = np.nanmedian(profs, axis=0)
    lo, hi = np.nanpercentile(profs, [16, 84], axis=0)
    fw_star = [fwhm(p, rc) for p in profs]
    out[band] = dict(n=len(profs), fwhm_psf=fwhm(psf, rc), fwhm_j2312=fwhm(pj, rc),
                     fwhm_field_med=float(np.median(fw_star)), fwhm_field_95=float(np.percentile(fw_star, 95)))
    ax.fill_between(rc, lo, hi, color=GREY, alpha=0.35, lw=0, label=f"field stars, 16-84% (n = {len(profs)})")
    ax.plot(rc, psf, color=INK2, lw=2, label=f"median field star (PSF), FWHM {out[band]['fwhm_psf']:.1f}\"")
    if band == "w4":
        for g, ls in ((10, ":"), (13, "--")):
            # broaden the radial PSF profile: model a 2-D image from it, convolve, re-measure
            yy, xx = np.mgrid[-HALF:HALF + 1, -HALF:HALF + 1]
            rr = np.hypot(xx, yy) * PIX
            img2 = np.interp(rr, rc, psf, right=0)
            conv = ndimage.gaussian_filter(img2, g / 2.3548 / PIX)
            yy2, xx2 = np.mgrid[:conv.shape[0], :conv.shape[1]]
            r2 = np.hypot(xx2 - HALF, yy2 - HALF) * PIX
            pm = np.array([np.nanmean(conv[(r2 >= a) & (r2 < b)]) for a, b in zip(edges[:-1], edges[1:])])
            pm /= pm[0]
            ax.plot(rc, pm, color=ORANGE, lw=1.5, ls=ls,
                    label=f"PSF blurred by a {g}\" Gaussian (FWHM {fwhm(pm, rc):.1f}\")")
    ax.plot(rc, pj, color=BLUE, lw=2.5, marker="o", ms=4, label=f"J2312, FWHM {out[band]['fwhm_j2312']:.1f}\"")
    ax.axhline(0.5, color=GRID, lw=1, zorder=0)
    ax.set_xlim(0, 40); ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("radius from source centre (arcsec)")
    ax.set_title(label, loc="left", fontsize=12, fontweight="bold", color=INK)
    ax.legend(fontsize=8.5, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(color=GRID, lw=0.6)
axs[0].set_ylabel("brightness / peak (local background removed)")
fig.suptitle("J2312 is slightly broader than field stars at 12 µm and much broader at 22 µm",
             x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
fig.tight_layout()
fig.savefig(D / "j2312_w4_profile.png", dpi=170, facecolor=SURF)
for b, v in out.items():
    print(b, {k: (round(x, 2) if isinstance(x, float) else x) for k, x in v.items()})
