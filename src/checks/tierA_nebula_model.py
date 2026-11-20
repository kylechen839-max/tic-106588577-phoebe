"""Tier A: J2312's method B (W3 nebula template + PSF + plane, j2312_psf_extraction.py) applied to every tier A star,
plus the W4 width after the nebula model is removed.

Per star: AllWISE Atlas W3/W4 strips (range reads from the IRSA AWS mirror, +-1400 px), empirical PSFs from isolated
field stars in the strip (>= 300" from the target), then
  a0: PSF + flat sky (50-70") -> catalogue-like reference, scaled to AllWISE w4mpro (flux ratios cancel zero points)
  a3: PSF + quadratic (< 55")
  B : PSF + scaled, smoothed W3 template + plane (< 150")
The excess sigma uses the catalogue W4 error only. The photosphere is the binary photosphere of tierA_geometry.csv
(F_cat / W4_ratio_bin). Width: Gaussian + plane FWHM within 20" of (data - B nebula model) vs the same fit on the PSF.
Output: outputs/crossmatch/tierA_w4_shape/tierA_nebula_model.csv
Usage: .venv-phoebe/bin/python src/checks/tierA_nebula_model.py [designation ...]
"""
import importlib.util, sys
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy import ndimage


def mod(name, f):
    s = importlib.util.spec_from_file_location(name, str(Path(__file__).with_name(f)))
    x = importlib.util.module_from_spec(s); s.loader.exec_module(x); return x


m = mod("m", "j2312_psf_extraction.py")
t = mod("t", "tierA_w4_shape.py")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/crossmatch/tierA_w4_shape"
PIX = m.PIX
H = 1400                                  # +-1925": enough isolated bright W4 stars for the PSF


def strip_wcs(tile, band, c):
    url = f"{t.BUCKET}{tile}{tile.rstrip('/').split('/')[-1]}-w{band}-int-3.fits"
    raw = t.get(url, (0, 28799)).decode("ascii", "ignore")
    end = raw.index("END" + " " * 77)
    hdr = fits.Header.fromstring(raw[:end + 80]); nhdr = (end // 2880 + 1) * 2880
    w = WCS(hdr); nx, ny = hdr["NAXIS1"], hdr["NAXIS2"]
    x0, y0 = (float(v) for v in w.world_to_pixel(c))
    y_lo, y_hi = max(0, int(y0) - H), min(ny - 1, int(y0) + H)
    img = np.frombuffer(t.get(url, (nhdr + y_lo * nx * 4, nhdr + (y_hi + 1) * nx * 4 - 1)), dtype=">f4").reshape(-1, nx)
    x_lo, x_hi = max(0, int(x0) - H), min(nx - 1, int(x0) + H)
    return img[:, x_lo:x_hi + 1].astype(float), w[y_lo:y_hi + 1, x_lo:x_hi + 1]


def run(row, geo):
    c = SkyCoord(row.ra, row.dec, unit="deg")
    m.C = c                                   # build_psf / cutout use the module-level target position
    tile = t.nearest_tile(row.ra, row.dec)
    img4, w4 = strip_wcs(tile, 4, c); img3, w3 = strip_wcs(tile, 3, c)
    psf4, n4, fw4 = m.build_psf(img4, w4); psf3, n3, fw3 = m.build_psf(img3, w3)
    hp = int(160 / PIX)
    s4, dx, dy = m.cutout(img4, w4, hp); s3, dx3, dy3 = m.cutout(img3, w3, hp)
    yy, xx = np.mgrid[:s4.shape[0], :s4.shape[1]]
    X, Y = (xx - hp - dx) * PIX, (yy - hp - dy) * PIX
    R = np.hypot(X, Y); one = np.ones_like(s4)
    P4 = m.placed_psf(psf4, s4.shape, dx, dy); P3 = m.placed_psf(psf3, s3.shape, dx3, dy3)
    sky = np.nanmedian(s4[(R >= 50) & (R < 70)])
    a0 = np.sum((P4 * (s4 - sky))[R < 55]) / np.sum((P4 * P4)[R < 55])
    c3q, _ = m.lstsq([P4, one, X, Y, X * X, Y * Y, X * Y], s4, R < 55)
    coef3, _ = m.lstsq([P3, one, X, Y], s3, R < 40)
    neb3 = s3 - coef3[0] * P3
    sm = ndimage.median_filter(neb3, size=int(30 / PIX)); resid = neb3 - sm
    nn = 1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid)))
    pts = ndimage.binary_dilation(resid > 5 * nn, iterations=int(8 / PIX)) & (R > 20)
    neb3[pts] = sm[pts]
    tmpl = ndimage.gaussian_filter(neb3, np.sqrt(max(fw4 ** 2 - fw3 ** 2, 1.0)) / 2.3548 / PIX)
    tmpl -= np.nanmedian(tmpl[(R > 140) & (R < 160)])
    cB, _ = m.lstsq([P4, tmpl, one, X, Y], s4, R < 150)
    neb = cB[1] * tmpl + cB[2] + cB[3] * X + cB[4] * Y
    def width(a):
        try:
            return t.gfit(a, hp + dx, hp + dy, 20)[3]
        except ValueError:                    # no compact peak within 20": the fit hits a bound
            return np.nan
    fw_sub = width(s4 - neb)
    fw_raw = width(s4)
    big = np.zeros((200, 200)); big[100 - m.HALF:100 + m.HALF + 1, 100 - m.HALF:100 + m.HALF + 1] = psf4 * 1e3
    fw_psf = t.gfit(big, 100, 100, 20)[3]
    Fcat = 8.363e3 * 10 ** (-0.4 * row.W4mag); Ecat = Fcat * 0.4 * np.log(10) * row.e_W4mag
    Fph = Fcat / geo.W4_ratio_bin
    out = dict(designation=row.designation, tile=tile.rstrip("/").split("/")[-1], n_psf4=n4, psf_fwhm_w4=round(fw4, 1),
               F_cat=round(Fcat, 2), E_cat=round(Ecat, 2), F_phot=round(Fph, 2),
               s_flux=round(cB[1] * 5.277e-5 / 1.998e-6, 2))
    for kk, a in (("a0", a0), ("a3", c3q[0]), ("B", cB[0])):
        F = Fcat * a / a0
        out[f"F_{kk}"] = round(F, 2); out[f"ratio_{kk}"] = round(F / Fph, 2); out[f"sigma_{kk}"] = round((F - Fph) / Ecat, 1)
    out.update(fwhm20_raw=round(fw_raw, 1), fwhm20_nebsub=round(fw_sub, 1), fwhm20_psf=round(fw_psf, 1))
    return out


def main():
    geo = pd.read_csv(ROOT / "outputs/crossmatch/tierA_geometry/tierA_geometry.csv").set_index("designation")
    cand = pd.read_csv(ROOT / "outputs/crossmatch/candidates_all_in_one.csv").set_index("designation")
    pos = pd.read_csv(ROOT / "outputs/crossmatch/dd1_candidates_scored.csv").set_index("designation")[["ra", "dec"]]
    names = sys.argv[1:] or list(geo.index)
    rows = []
    for d in names:
        r = cand.loc[d].copy(); r["ra"], r["dec"] = pos.loc[d, "ra"], pos.loc[d, "dec"]; r.name = d
        r["designation"] = d
        try:
            rows.append(run(r, geo.loc[d]))
        except Exception as e:
            rows.append(dict(designation=d, note=f"failed: {type(e).__name__}: {e}"[:120]))
        print(rows[-1], flush=True)
    tb = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    tb.to_csv(OUT / "tierA_nebula_model.csv", index=False)
    print(tb.drop(columns=[c for c in ("tile",) if c in tb]).to_string(index=False))


if __name__ == "__main__":
    main()
