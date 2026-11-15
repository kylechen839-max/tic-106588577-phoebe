"""Tier A: is the W4 source a point source centred on the star? Two literature checks plus a field-calibrated shape test.

1. Koenig & Leisawitz (2014, ApJ 791, 131) W4 cut: keep W4 only if w4rchi2 < 0.2*w4snr - 2 (catalogue values,
   AllWISE p3as_psd via IRSA SCS, copied into ALLWISE below).
2. W3/W4 photocentre test (after Patel et al. 2017, AJ 153, 54): Gaussian+plane fits in each band within 20";
   offset of each band's centre from the AllWISE position, and W3-to-W4 offset.
3. Gaussian FWHM with a plane background, 20" and 30" fit radii (W3 at 20" as a control: the W3 light is mostly
   photosphere, so it should be point-like). The 40-star empirical PSF of tile 3476p530 gives Gaussian FWHM
   16.1" (W4, 20"), 20.2" (W4, 30") and 13.3" (W3, 20"), matching the field-source medians below.
Tests 2 and 3 are calibrated against isolated field point sources in the same tile (same PSF, same processing), so
no noise model is assumed. Bright W4 sources are rare, so every isolated field source brighter than A_target/3 is used.
Position noise scales as 1/SNR, so offsets are compared as offset x SNR (SNR = Gaussian amplitude / pixel noise);
the result is the target's percentile among the field sources.

Images: AllWISE Atlas tiles (1.375"/px) read with HTTP range requests from the IRSA AWS mirror, a 1201-row strip
per band, so no full tile is stored. Output: outputs/crossmatch/tierA_w4_shape/.
Usage: .venv-phoebe/bin/python src/checks/tierA_w4_shape.py
"""
import re, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy import ndimage
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch" / "tierA_w4_shape"
BUCKET = "https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/"
PREFIX = "wise/allwise/images/p3am_cdd/"
PIX = 1.375
HALF_STRIP = 600       # rows/cols either side of the target (825")
STAMP = 30             # stamp half-size (px) for fits

# AllWISE p3as_psd, IRSA SCS cone 5", queried 2026-11-15 (J2312 from outputs/j2312_contamination/*_allwise_3arcmin.csv)
ALLWISE = """designation,w3snr,w3rchi2,w4snr,w4rchi2,cc_flags
J073937.23-363011.8,52.5,1.432,21.1,1.207,0000
J221843.61+544715.0,32.8,2.026,14.3,1.155,0000
J231201.40+532028.6,51.3,,17.5,1.139,0000
J024542.12+480837.8,57.4,1.381,16.6,1.165,hh00
J223553.38+500414.8,72.2,0.8108,23.1,1.258,0h00
J160415.99-562627.3,33.3,1.094,12.0,1.172,0000
J141909.42-565518.1,43.6,0.7659,13.0,0.9723,0h00
J045912.77+165543.8,49.2,1.023,12.8,1.056,0000
J220912.52+582025.0,54.7,1.119,11.6,0.8963,hh00
J033226.58+520357.1,51.0,0.9544,12.1,1.130,0000
J173224.94-364153.9,29.3,0.7261,10.7,1.168,0000
J223949.47+583254.4,59.3,1.490,11.2,0.8882,hh00
"""


def get(url, rng=None):
    req = urllib.request.Request(url, headers={"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {})
    for _ in range(4):
        try:
            return urllib.request.urlopen(req, timeout=120).read()
        except Exception as e:
            err = e
    raise err


def list_prefixes(prefix):
    x = get(f"{BUCKET}?list-type=2&delimiter=/&prefix={prefix}").decode()
    return re.findall(r"<Prefix>([^<]*)</Prefix>", x)[1:]


def nearest_tile(ra, dec):
    """Atlas tile that contains the target with the largest margin (tiles are 1.564 deg squares on a tangent plane;
    the nearest centre is not always the containing tile). Folders are named by RA*10."""
    c = SkyCoord(ra, dec, unit="deg")
    half = 1.7 / max(np.cos(np.radians(dec)), 0.05)
    lo, hi = (ra - half) * 10, (ra + half) * 10
    best = (9e9, None)
    for top in sorted({f"{int(v) % 3600:04d}"[:2] for v in np.arange(lo, hi + 10, 10)}):
        for sub in list_prefixes(f"{PREFIX}{top}/"):
            v = int(sub.rstrip("/").split("/")[-1])
            if not any(lo <= v + k * 3600 <= hi for k in (-1, 0, 1)):
                continue
            for t in list_prefixes(sub):
                n = t.rstrip("/").split("/")[-1]
                tra, tdec = int(n[:4]) / 10, int(n[5:8]) / 10 * (1 if n[4] == "p" else -1)
                off = SkyCoord(tra, tdec, unit="deg").spherical_offsets_to(c)
                margin = 0.782 - max(abs(off[0].deg), abs(off[1].deg))
                best = min(best, (-margin, t))
    return best[1]


def strip(tile_prefix, band, ra, dec):
    url = f"{BUCKET}{tile_prefix}{tile_prefix.rstrip('/').split('/')[-1]}-w{band}-int-3.fits"
    raw = get(url, (0, 28799))
    hdr = fits.Header.fromstring(raw.decode("ascii", "ignore")[: raw.decode("ascii", "ignore").index("END" + " " * 77) + 80])
    nhdr = (raw.decode("ascii", "ignore").index("END" + " " * 77) // 2880 + 1) * 2880
    w = WCS(hdr)
    nx, ny = hdr["NAXIS1"], hdr["NAXIS2"]
    x0, y0 = (float(v) for v in w.world_to_pixel(SkyCoord(ra, dec, unit="deg")))
    y_lo, y_hi = max(0, int(y0) - HALF_STRIP), min(ny - 1, int(y0) + HALF_STRIP)
    assert 60 < x0 < nx - 60 and 60 < y0 < ny - 60, f"target at ({x0:.0f},{y0:.0f}) is off tile {url}"
    data = get(url, (nhdr + y_lo * nx * 4, nhdr + (y_hi + 1) * nx * 4 - 1))
    img = np.frombuffer(data, dtype=">f4").reshape(-1, nx).astype(float)
    assert img.shape[0] == y_hi - y_lo + 1
    x_lo, x_hi = max(0, int(x0) - HALF_STRIP), min(nx - 1, int(x0) + HALF_STRIP)
    return img[:, x_lo:x_hi + 1], x0 - x_lo, y0 - y_lo


def gfit(img, x0, y0, rmax, fix_centre=False):
    """Circular Gaussian + plane within rmax (arcsec) of (x0, y0). Returns (amp, xc, yc, fwhm_arcsec)."""
    xi, yi = int(round(x0)), int(round(y0))
    st = img[yi - STAMP:yi + STAMP + 1, xi - STAMP:xi + STAMP + 1]
    st = np.nan_to_num(st, nan=np.nanmedian(st))
    yy, xx = np.mgrid[:st.shape[0], :st.shape[1]]
    cx0, cy0 = STAMP + x0 - xi, STAMP + y0 - yi
    X, Y = xx - cx0, yy - cy0
    m = np.hypot(X, Y) * PIX < rmax

    def f(p):
        A, xc, yc, s, c0, gx, gy = p
        return (A * np.exp(-((xx - xc) ** 2 + (yy - yc) ** 2) / (2 * s * s)) + c0 + gx * X + gy * Y - st)[m]

    a0 = max(st[int(round(cy0)), int(round(cx0))] - np.median(st[m]), 1e-3 * np.nanstd(st[m]) + 1e-9)
    p0 = [a0, cx0, cy0, 5.0, np.median(st[m]), 0, 0]
    p = least_squares(f, p0, bounds=([0, cx0 - 8, cy0 - 8, 2, -np.inf, -np.inf, -np.inf],
                                     [np.inf, cx0 + 8, cy0 + 8, 12, np.inf, np.inf, np.inf])).x
    if p[3] < 2.05 or p[3] > 11.95 or abs(p[1] - cx0) > 7.9 or abs(p[2] - cy0) > 7.9 or p[0] <= 0:
        raise ValueError("fit hit a bound")
    return p[0], p[1] - STAMP + xi, p[2] - STAMP + yi, 2.3548 * abs(p[3]) * PIX


def pixel_noise(img):
    res = img - ndimage.median_filter(img, size=41)
    return 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res)))


def field_sources(img4, x0, y0):
    bg = ndimage.median_filter(img4, size=41)
    res = img4 - bg
    noise = pixel_noise(img4)
    sm = ndimage.gaussian_filter(res, 2)
    mx = ndimage.maximum_filter(sm, size=25)
    pk = np.argwhere((sm == mx) & (sm > 6 * noise))
    out = []
    for y, x in pk:
        if not (STAMP + 2 < y < img4.shape[0] - STAMP - 2 and STAMP + 2 < x < img4.shape[1] - STAMP - 2):
            continue
        if np.hypot(x - x0, y - y0) * PIX < 60:
            continue
        # isolation: no other peak > 30% of this one within 40"
        near = [(yy, xx) for yy, xx in pk if 0 < np.hypot(xx - x, yy - y) * PIX < 40 and sm[yy, xx] > 0.3 * sm[y, x]]
        if near:
            continue
        out.append((x, y))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    aw = pd.read_csv(pd.io.common.StringIO(ALLWISE))
    pos = pd.read_csv(ROOT / "outputs/crossmatch/dd1_candidates_scored.csv")[["designation", "ra", "dec"]]
    aw = aw.merge(pos, on="designation")
    aw["kl14_limit"] = 0.2 * aw.w4snr - 2
    aw["kl14_pass"] = aw.w4rchi2 < aw.kl14_limit
    rows, fld = [], []
    for r in aw.itertuples():
        tile = nearest_tile(r.ra, r.dec)
        im3, x3, y3 = strip(tile, 3, r.ra, r.dec)
        im4, x4, y4 = strip(tile, 4, r.ra, r.dec)
        try:
            A3, xc3, yc3, fw3 = gfit(im3, x3, y3, 20)
            A4, xc4, yc4, fw4c = gfit(im4, x4, y4, 20)
            _, _, _, fw4 = gfit(im4, x4, y4, 30)
        except ValueError as e:
            print(r.designation, "target fit failed:", e, flush=True)
            rows.append(dict(designation=r.designation, tile=tile.rstrip("/").split("/")[-1], w4snr=r.w4snr,
                             w4rchi2=r.w4rchi2, kl14_limit=round(r.kl14_limit, 2), kl14_pass=r.kl14_pass,
                             note="Gaussian fit within 20\" hit a bound (no compact W4 peak at the star)"))
            continue
        off3 = np.hypot(xc3 - x3, yc3 - y3) * PIX
        off4 = np.hypot(xc4 - x4, yc4 - y4) * PIX
        off34 = np.hypot(xc4 - xc3, yc4 - yc3) * PIX
        n4 = pixel_noise(im4)
        snr4 = A4 / n4
        sat = 0.98 * np.nanmax(im4)
        # field calibration: same fits on isolated W4 point sources brighter than A4/3
        f = []
        for x, y in field_sources(im4, x4, y4):
            try:
                a4, fx4, fy4, _ = gfit(im4, x, y, 20)
                if not (a4 > A4 / 3) or im4[int(round(fy4)), int(round(fx4))] >= sat:
                    continue
                a3, fx3, fy3, _ = gfit(im3, fx4, fy4, 20)
                _, _, _, ffw = gfit(im4, fx4, fy4, 30)
                _, _, _, ffw20 = gfit(im4, fx4, fy4, 20)
                f.append((a4, np.hypot(fx4 - fx3, fy4 - fy3) * PIX * a4 / n4, ffw, ffw20))
            except Exception:
                continue
        f = np.array(f) if f else np.zeros((0, 4))
        n = len(f)
        rows.append(dict(designation=r.designation, tile=tile.rstrip("/").split("/")[-1], w4snr=r.w4snr, w4rchi2=r.w4rchi2,
                         kl14_limit=round(r.kl14_limit, 2), kl14_pass=r.kl14_pass,
                         W3_offset_arcsec=round(off3, 2), W4_offset_arcsec=round(off4, 2), W3W4_offset_arcsec=round(off34, 2),
                         W4_fwhm30=round(fw4, 1), W4_snr_img=round(snr4, 1), n_field=n,
                         W3W4_offset_x_snr=round(off34 * snr4, 1),
                         field_offset_x_snr_p50=round(np.percentile(f[:, 1], 50), 1) if n else np.nan,
                         field_offset_x_snr_p95=round(np.percentile(f[:, 1], 95), 1) if n else np.nan,
                         W3W4_offset_pctile=round(100 * np.mean(f[:, 1] < off34 * snr4), 0) if n else np.nan,
                         field_fwhm_p50=round(np.percentile(f[:, 2], 50), 1) if n else np.nan,
                         field_fwhm_p95=round(np.percentile(f[:, 2], 95), 1) if n else np.nan,
                         W4_fwhm_pctile=round(100 * np.mean(f[:, 2] < fw4), 0) if n else np.nan,
                         W3_fwhm20=round(fw3, 1), W4_fwhm20=round(fw4c, 1),
                         field_fwhm20_p50=round(np.percentile(f[:, 3], 50), 1) if n else np.nan,
                         field_fwhm20_p95=round(np.percentile(f[:, 3], 95), 1) if n else np.nan,
                         W4_fwhm20_pctile=round(100 * np.mean(f[:, 3] < fw4c), 0) if n else np.nan))
        print({k: (float(v) if isinstance(v, np.floating) else v) for k, v in rows[-1].items()}, flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "tierA_w4_shape.csv", index=False)
    print(t.drop(columns=["tile"]).to_string(index=False))


if __name__ == "__main__":
    main()
