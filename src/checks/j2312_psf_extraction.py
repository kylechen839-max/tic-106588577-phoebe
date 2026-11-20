"""J2312: W4 point-source flux with the surrounding nebula modelled, from native AllWISE Atlas tiles.

Tiles: 3476p530_ac51 W3/W4 int (1.375"/px), downloaded from the IRSA AWS mirror into
outputs/j2312_contamination/native/ (67 MB each; never commit them).

1. Empirical PSFs: median stack of isolated, bright, unsaturated point sources in the tile
   (>= 5' from J2312), each centred to sub-pixel accuracy and normalised to unit peak-fit amplitude.
2. Method A (PSF + background model): fit A * PSF at the star plus one of
     a0  flat sky fixed to the 50-70" median (mimics the AllWISE catalogue),
     a1  free constant, a2 plane, a3 quadratic (all within 55"),
     a4  quadratic + a circular Gaussian halo centred on the star (width on a grid),
   so the nebula is fitted, not just subtracted from far away.
3. Method B (W3 nebula template): remove J2312 and other point sources from the W3 image, smooth
   the residual to W4 resolution, and fit W4 = A * PSF + s * template + plane within 150". The
   nebula's own W3 morphology then predicts the W4 nebula under the star.
Every flux is expressed relative to a0 and scaled to the AllWISE W4 flux (6.840 mag = 15.36 mJy),
so the zero point and PSF normalisation cancel. The photosphere is the binary photosphere of the
adopted PHOEBE fit (W4 ratio 3.561 -> 4.31 mJy).
"""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "outputs" / "j2312_contamination"
NAT = D / "native"
TILE = "3476p530_ac51"
C = SkyCoord(348.0058494, 53.3413022, unit="deg")
PIX = 1.375
F_CAT = 8.363e3 * 10 ** (-0.4 * 6.840)   # mJy, AllWISE w4mpro
E_CAT = F_CAT * 0.4 * np.log(10) * 0.062
F_PHOT = F_CAT / 3.561                    # binary photosphere (tierA_geometry.csv)
HALF = 45                                 # PSF stamp half-size (px) = 62"


def load(band):
    with fits.open(NAT / f"{TILE}-{band}-int-3.fits") as h:
        img = h[0].data.astype(float)
        w = WCS(h[0].header)
    return img, w


def background(img):
    small = img[: img.shape[0] // 8 * 8, : img.shape[1] // 8 * 8].reshape(img.shape[0] // 8, 8, img.shape[1] // 8, 8)
    med = np.nanmedian(small, axis=(1, 3))
    med = ndimage.median_filter(med, size=9)
    bg = ndimage.zoom(med, 8, order=1)
    out = np.full_like(img, np.nanmedian(img))
    out[: bg.shape[0], : bg.shape[1]] = bg
    return out


def centroid_shift(stamp):
    yy, xx = np.mgrid[: stamp.shape[0], : stamp.shape[1]]
    c = HALF
    m = np.hypot(xx - c, yy - c) < 4
    s = np.clip(stamp, 0, None) * m
    return (s * yy).sum() / s.sum() - c, (s * xx).sum() / s.sum() - c


def build_psf(img, w, nmax=40):
    bg = background(img)
    res = img - bg
    noise = 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res)))
    mx = ndimage.maximum_filter(res, size=31)
    pk = np.argwhere((res == mx) & (res > 20 * noise))
    x0, y0 = w.world_to_pixel(C)
    sat = 0.98 * np.nanmax(img)
    cand = []
    for y, x in pk:
        if not (HALF + 5 < y < img.shape[0] - HALF - 5 and HALF + 5 < x < img.shape[1] - HALF - 5):
            continue
        if np.hypot(x - x0, y - y0) * PIX < 300 or img[y, x] >= sat:
            continue
        st = res[y - HALF: y + HALF + 1, x - HALF: x + HALF + 1]
        yy, xx = np.mgrid[: st.shape[0], : st.shape[1]]
        rr = np.hypot(xx - HALF, yy - HALF)
        outer = ndimage.gaussian_filter(st, 2)[rr > 25]
        if np.nanmax(outer) > 0.12 * res[y, x]:   # a neighbour 34-62" away
            continue
        cand.append((res[y, x], y, x))
    cand.sort(reverse=True)
    stamps = []
    for amp, y, x in cand[:nmax]:
        st = res[y - HALF: y + HALF + 1, x - HALF: x + HALF + 1].copy()
        dy, dx = centroid_shift(st)
        st = ndimage.shift(st, (-dy, -dx), order=3, mode="nearest")
        yy, xx = np.mgrid[: st.shape[0], : st.shape[1]]
        rr = np.hypot(xx - HALF, yy - HALF) * PIX
        st -= np.nanmedian(st[(rr > 50) & (rr < 60)])
        stamps.append(st / np.nansum(st[rr < 40]))
    psf = np.nanmedian(stamps, axis=0)
    yy, xx = np.mgrid[: psf.shape[0], : psf.shape[1]]
    rr = np.hypot(xx - HALF, yy - HALF) * PIX
    psf -= np.nanmedian(psf[(rr > 50) & (rr < 60)])
    psf /= np.nansum(psf[rr < 40])
    # FWHM from the radial profile
    prof = [np.nanmean(psf[(rr >= a) & (rr < a + PIX)]) for a in np.arange(0, 40, PIX)]
    half = np.interp(0.5 * prof[0], prof[::-1], np.arange(0, 40, PIX)[::-1])
    return psf, len(stamps), 2 * half


def cutout(img, w, half_px):
    x0, y0 = w.world_to_pixel(C)
    xi, yi = int(round(float(x0))), int(round(float(y0)))
    sub = img[yi - half_px: yi + half_px + 1, xi - half_px: xi + half_px + 1].copy()
    return sub, float(x0) - xi, float(y0) - yi


def placed_psf(psf, shape, dx, dy):
    big = np.zeros(shape)
    cy, cx = shape[0] // 2, shape[1] // 2
    big[cy - HALF: cy + HALF + 1, cx - HALF: cx + HALF + 1] = psf
    return ndimage.shift(big, (dy, dx), order=3)


def lstsq(cols, data, mask):
    A = np.column_stack([c[mask] for c in cols])
    coef, *_ = np.linalg.lstsq(A, data[mask], rcond=None)
    model = sum(c * k for c, k in zip(cols, coef))
    return coef, model


def main():
    lines = []
    img4, w4 = load("w4")
    img3, w3 = load("w3")
    psf4, n4, fw4 = build_psf(img4, w4)
    psf3, n3, fw3 = build_psf(img3, w3)
    lines.append(f"empirical PSFs: W4 {n4} stars, FWHM {fw4:.1f}\"; W3 {n3} stars, FWHM {fw3:.1f}\"")

    hp = int(160 / PIX)
    s4, dx, dy = cutout(img4, w4, hp)
    yy, xx = np.mgrid[: s4.shape[0], : s4.shape[1]]
    cy = cx = hp
    X, Y = (xx - cx - dx) * PIX, (yy - cy - dy) * PIX
    R = np.hypot(X, Y)
    P4 = placed_psf(psf4, s4.shape, dx, dy)
    one = np.ones_like(s4)
    noise4 = 1.4826 * np.nanmedian(np.abs(s4[(R > 100) & (R < 160)] - np.nanmedian(s4[(R > 100) & (R < 160)])))

    res = {}
    # a0: flat sky fixed at the 50-70" median (catalogue-like)
    sky = np.nanmedian(s4[(R >= 50) & (R < 70)])
    m = R < 55
    a0 = np.sum((P4 * (s4 - sky))[m]) / np.sum((P4 * P4)[m])
    res["a0 flat sky 50-70\" (catalogue-like)"] = (a0, s4 - sky - a0 * P4)
    terms = {"a1 free constant": [one], "a2 plane": [one, X, Y], "a3 quadratic": [one, X, Y, X * X, Y * Y, X * Y]}
    for k, t in terms.items():
        coef, model = lstsq([P4] + t, s4, m)
        res[k] = (coef[0], s4 - model)
    best = None
    for sig in (8, 12, 16, 20, 25, 30, 40, 50):
        G = np.exp(-R ** 2 / (2 * sig ** 2))
        coef, model = lstsq([P4, G, one, X, Y, X * X, Y * Y, X * Y], s4, R < 120)
        chi = np.nansum(((s4 - model)[R < 120]) ** 2)
        if coef[1] >= 0 and (best is None or chi < best[0]):
            best = (chi, sig, coef[0], s4 - model, coef[1])
    res[f"a4 quad + Gaussian halo (sigma {best[1]}\"), DEGENERATE*"] = (best[2], best[3])

    # Method B: W3 nebula template
    s3, dx3, dy3 = cutout(img3, w3, hp)
    P3 = placed_psf(psf3, s3.shape, dx3, dy3)
    coef3, model3 = lstsq([P3, one, X, Y], s3, R < 40)
    neb3 = s3 - coef3[0] * P3
    # remove other compact W3 sources (they are mostly stars, faint at 22 um)
    sm = ndimage.median_filter(neb3, size=int(30 / PIX))
    resid = neb3 - sm
    n3 = 1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid)))
    pts = ndimage.binary_dilation(resid > 5 * n3, iterations=int(8 / PIX)) & (R > 20)
    neb3[pts] = sm[pts]
    ksig = np.sqrt(max(fw4 ** 2 - fw3 ** 2, 1.0)) / 2.3548 / PIX
    tmpl = ndimage.gaussian_filter(neb3, ksig)
    tmpl -= np.nanmedian(tmpl[(R > 140) & (R < 160)])
    coefB, modelB = lstsq([P4, tmpl, one, X, Y], s4, R < 150)
    res["B W3 nebula template + plane within 150\""] = (coefB[0], s4 - modelB)
    lines.append(f"W3 template scale s = {coefB[1]:.4f} (W4 DN per W3 DN); W3 star fit removed, {pts.sum()} px of other W3 sources masked")

    lines.append("")
    lines.append(f"{'method':58s} {'F/F(a0)':>8s} {'W4 mJy':>8s} {'F/F_phot':>9s} {'excess sigma':>12s}")
    rows = []
    for k, (amp, r) in res.items():
        f = amp / a0
        F = F_CAT * f
        sigma = (F - F_PHOT) / E_CAT
        rows.append((k, f, F, F / F_PHOT, sigma, r))
        lines.append(f"{k:58s} {f:8.3f} {F:8.2f} {F / F_PHOT:9.2f} {sigma:12.1f}")
    lines.append("")
    lines.append(f"photosphere {F_PHOT:.2f} mJy; catalogue error {E_CAT:.2f} mJy (0.062 mag); fit-region pixel noise {noise4:.3f} DN")
    lines.append("* a4 is degenerate: a 12\" Gaussian trades off against the 17\" PSF. The same fit on the W3 image, where J2312 is mostly photosphere, also moves 45% of the flux into the halo, so a4 is not used.")
    lines.append("A shared limit: nebula structure narrower than the W4 PSF and centred on the star is indistinguishable from the star in all methods.")
    (D / "j2312_psf_extraction_log.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    show = [rows[0], rows[3], rows[4], rows[5]]
    fig, axs = plt.subplots(1, 5, figsize=(17, 3.6))
    v = np.nanpercentile(s4[R < 150], [2, 99.5])
    axs[0].imshow(s4, origin="lower", cmap="magma", vmin=v[0], vmax=v[1]); axs[0].set_title("W4 data (5.3′)", fontsize=9)
    for ax, (k, f, F, ratio, sg, r) in zip(axs[1:], show):
        lim = 4 * noise4
        rd = np.where(R < (150 if k.startswith("B") else 120 if k.startswith("a4") else 55), r, np.nan)
        ax.imshow(rd, origin="lower", cmap="coolwarm", vmin=-lim, vmax=lim)
        ax.set_title(f"residual: {k.split(' (')[0].split(' within')[0].replace(', DEGENERATE*', ' (degenerate)')[:34]}\nW4 = {F:.1f} mJy, {ratio:.2f}× phot", fontsize=8)
    for ax in axs:
        ax.set_xticks([]); ax.set_yticks([])
        ax.add_patch(plt.Circle((hp + dx, hp + dy), 16.5 / PIX, fill=False, color="w", lw=0.6))
    fig.suptitle("J2312 W4: star (PSF) fitted together with a nebula model; residuals after removing both", fontsize=10)
    fig.tight_layout(); fig.savefig(D / "j2312_psf_extraction.png", dpi=130)


if __name__ == "__main__":
    main()
