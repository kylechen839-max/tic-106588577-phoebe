"""J2312: W4-only multi-scale separation of point sources, filaments and background, after getsf (Men'shchikov 2021,
A&A 649, A89). getsf itself is not reachable from the cloud container, so this is a simplified Python version of
its three steps, using W4 alone (no W3 colour assumption, unlike method B in j2312_psf_extraction.py):

1. Background: iterative median filtering with a window W (getsf: a few times the largest source size); pixels more
   than 2 sigma above the current background are replaced by it before the next pass, 15 passes.
2. Single-scale decomposition of image - background into difference-of-Gaussian scales (FWHM 4" to 120", factor
   1.25); in each scale, pixels above 3 sigma are labelled into connected structures; a structure is a "source" if it
   is round (second-moment axis ratio < 1.5) and not much larger than the scale, otherwise a "filament".
3. Source image = sum of source structures over scales. The J2312 flux is a PSF fit on the source image.
Flux is calibrated two ways: (a) relative to the catalogue-like a0 fit (F_CAT), as before; (b) per-image recovery
factor from isolated field stars run through the same pipeline (median of source-image/raw-image PSF amplitude).
Run for W in 60/90/120/180" to show the sensitivity to the one free choice.
"""
import importlib.util
from pathlib import Path
import numpy as np
from scipy import ndimage

spec = importlib.util.spec_from_file_location("m", str(Path(__file__).with_name("j2312_psf_extraction.py")))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
PIX, HALF = m.PIX, m.HALF
OUT = m.D
SCALES = 4.0 * 1.25 ** np.arange(0, 16)    # FWHM in arcsec, 4" .. 114"
BOX = int(1200 / PIX)                      # work on a 2400" box around J2312


def bg_iter(img, win_arcsec, npass=15):
    w = int(win_arcsec / PIX) | 1
    small = 4                                  # median filter on a 4x-binned image for speed
    def medf(a):
        ny, nx = a.shape[0] // small * small, a.shape[1] // small * small
        b = np.nanmedian(a[:ny, :nx].reshape(ny // small, small, nx // small, small), axis=(1, 3))
        b = ndimage.median_filter(b, size=max(3, w // small | 1))
        z = ndimage.zoom(b, small, order=1)
        o = np.full_like(a, np.nanmedian(a)); o[:z.shape[0], :z.shape[1]] = z
        return o
    B = medf(img)
    for _ in range(npass):
        R = img - B
        s = 1.4826 * np.nanmedian(np.abs(R - np.nanmedian(R)))
        B = medf(np.where(R > 2 * s, B, img))
    return B


def decompose(R):
    src = np.zeros_like(R); fil = np.zeros_like(R)
    prev = R
    sig = SCALES / 2.3548 / PIX
    for k, s in enumerate(sig):
        cur = ndimage.gaussian_filter(R, s)
        D = prev - cur
        prev = cur
        n = 1.4826 * np.nanmedian(np.abs(D - np.nanmedian(D)))
        lab, nl = ndimage.label(D > 3 * n)
        if nl == 0:
            continue
        idx = np.arange(1, nl + 1)
        yy, xx = np.indices(D.shape)
        w = np.where(lab > 0, D, 0)
        sw = ndimage.sum(w, lab, idx)
        mx = ndimage.sum(w * xx, lab, idx) / sw; my = ndimage.sum(w * yy, lab, idx) / sw
        cxx = ndimage.sum(w * xx * xx, lab, idx) / sw - mx ** 2
        cyy = ndimage.sum(w * yy * yy, lab, idx) / sw - my ** 2
        cxy = ndimage.sum(w * xx * yy, lab, idx) / sw - mx * my
        tr, det = cxx + cyy, cxx * cyy - cxy ** 2
        l1 = tr / 2 + np.sqrt(np.maximum(tr ** 2 / 4 - det, 0)); l2 = np.maximum(tr / 2 - np.sqrt(np.maximum(tr ** 2 / 4 - det, 0)), 1e-6)
        axis = np.sqrt(l1 / l2)
        size = 2.3548 * np.sqrt(np.sqrt(l1 * l2)) * PIX    # equivalent Gaussian FWHM of the structure, arcsec
        is_src = (axis < 1.5) & (size < 2.5 * max(SCALES[k], 17.0))
        lut = np.zeros(nl + 1, bool); lut[1:] = is_src
        smask = lut[lab] & (lab > 0)
        fmask = (~lut[lab]) & (lab > 0)
        src += np.where(smask, D, 0)
        fil += np.where(fmask, D, 0)
    return src, fil


def psf_amp(img, psf, x0, y0, rfit=40, const=True):
    xi, yi = int(round(x0)), int(round(y0))
    hp = int(80 / PIX)
    st = img[yi - hp: yi + hp + 1, xi - hp: xi + hp + 1]
    P = m.placed_psf(psf, st.shape, x0 - xi, y0 - yi)
    yy, xx = np.indices(st.shape)
    msk = np.hypot(xx - hp - (x0 - xi), yy - hp - (y0 - yi)) * PIX < rfit
    cols = [P, np.ones_like(st)] if const else [P]
    coef, _ = m.lstsq(cols, st, msk & np.isfinite(st))
    return coef[0]


def field_stars(img, x0, y0, n=25):
    bg = m.background(img); res = img - bg
    noise = 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res)))
    mx = ndimage.maximum_filter(res, size=31)
    pk = np.argwhere((res == mx) & (res > 10 * noise))
    sat = 0.98 * np.nanmax(img)
    out = []
    for y, x in pk:
        if not (130 < y < img.shape[0] - 130 and 130 < x < img.shape[1] - 130): continue
        if np.hypot(x - x0, y - y0) * PIX < 120 or img[y, x] >= sat: continue
        st = ndimage.gaussian_filter(res[y - 30:y + 31, x - 30:x + 31], 2)
        rr = np.hypot(*np.indices(st.shape) - 30)
        if np.nanmax(st[rr > 18]) > 0.2 * res[y, x]: continue
        out.append((x, y))
    return out[:n]


def main():
    lines = []
    for band in ("w4", "w3"):
        img, w = m.load(band)
        psf, npsf, fw = m.build_psf(img, w)
        x0, y0 = (float(v) for v in w.world_to_pixel(m.C))
        xi, yi = int(round(x0)), int(round(y0))
        sub = img[yi - BOX: yi + BOX + 1, xi - BOX: xi + BOX + 1].copy()
        tx, ty = BOX + x0 - xi, BOX + y0 - yi
        # catalogue-like reference amplitude (a0: flat sky 50-70")
        hp = int(160 / PIX)
        s4 = sub[BOX - hp: BOX + hp + 1, BOX - hp: BOX + hp + 1]
        yy, xx = np.indices(s4.shape); R = np.hypot(xx - hp - (tx - BOX), yy - hp - (ty - BOX)) * PIX
        P = m.placed_psf(psf, s4.shape, tx - BOX, ty - BOX)
        sky = np.nanmedian(s4[(R >= 50) & (R < 70)])
        a0 = np.sum((P * (s4 - sky))[R < 55]) / np.sum((P * P)[R < 55])
        stars = field_stars(sub, tx, ty)
        lines.append(f"{band.upper()}: PSF from {npsf} stars, FWHM {fw:.1f}\"; {len(stars)} isolated field stars for flux recovery")
        for W in (60, 90, 120, 180):
            B = bg_iter(sub, W)
            src, fil = decompose(sub - B)
            a_src = psf_amp(src, psf, tx, ty, const=False)
            rec = [psf_amp(src, psf, x, y, const=False) / psf_amp(sub, psf, x, y) for x, y in stars]
            rec = np.median(rec) if rec else np.nan
            f_a0 = a_src / a0
            f_cal = a_src / rec / a0
            line = f"  W={W:3d}\"  source/a0 {f_a0:.3f}  field-star recovery {rec:.3f}  recovery-corrected/a0 {f_cal:.3f}"
            if band == "w4":
                F = m.F_CAT * f_cal
                line += f"  -> {F:.2f} mJy, {F / m.F_PHOT:.2f}x photosphere, {(F - m.F_PHOT) / m.E_CAT:.1f} sigma"
            lines.append(line)
            if band == "w4" and W == 90:
                keep = (src.copy(), fil.copy(), B.copy(), sub.copy(), tx, ty)
    lines.append("")
    lines.append(f"W4 photosphere {m.F_PHOT:.2f} mJy, catalogue {m.F_CAT:.2f} mJy, catalogue error {m.E_CAT:.2f} mJy.")
    lines.append("W3 rows are a control: J2312 W3 is mostly photosphere, so recovery-corrected/a0 should be close to 1 if the")
    lines.append("method removes nebula without eating point-source flux.")
    lines.append("Shared limit: nebula structure as compact as the PSF and centred on the star is classed as 'source' here too.")
    txt = "\n".join(lines)
    print(txt)
    (OUT / "j2312_multiscale_log.txt").write_text(txt + "\n")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    src, fil, B, sub, tx, ty = keep
    h = int(200 / PIX)
    cut = lambda a: a[int(ty) - h:int(ty) + h + 1, int(tx) - h:int(tx) + h + 1]
    fig, ax = plt.subplots(1, 4, figsize=(15, 4))
    v = np.nanpercentile(cut(sub - B), [1, 99.5])
    for a, im, t in zip(ax, [sub - np.nanmedian(cut(B)), B - np.nanmedian(cut(B)), fil, src],
                        ["W4 image", "background (W=90\")", "filaments / extended", "point sources"]):
        a.imshow(cut(im), origin="lower", cmap="magma", vmin=v[0], vmax=v[1], extent=[-200, 200, -200, 200])
        a.set_title(t); a.plot(0, 0, "c+", ms=10); a.set_xlabel("arcsec")
    fig.suptitle("J2312 W4 multi-scale separation (getsf-style, W4 only)")
    fig.tight_layout(); fig.savefig(OUT / "j2312_multiscale.png", dpi=130)


if __name__ == "__main__":
    main()
