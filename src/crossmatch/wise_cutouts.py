"""WISE W1-W4 image cutouts (SkyView AllWISE atlas) for candidates, with a centring check.

For each band the centroid of the brightest blob within 15" of the target and the
peak-to-background contrast are measured, so W3/W4 emission that is offset from the
star (a neighbour or nebulosity) can be flagged.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astroquery.skyview import SkyView

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch" / "wise_cutouts"
BANDS = ["WISE 3.4", "WISE 4.6", "WISE 12", "WISE 22"]


def plane_background(img, xx, yy, ann):
    """Plane a + b x + c y fitted to the 18-35" annulus with 3-sigma clipping; returns (bg image, noise)."""
    ok = ann & np.isfinite(img)
    for _ in range(5):
        A = np.c_[np.ones(ok.sum()), xx[ok], yy[ok]]
        coef, *_ = np.linalg.lstsq(A, img[ok], rcond=None)
        resid = img[ok] - A @ coef
        noise = 1.4826 * np.median(np.abs(resid - np.median(resid)))
        keep = np.abs(resid) < 3 * noise
        if keep.all():
            break
        idx = np.flatnonzero(ok)
        ok = ok.copy(); ok.flat[idx[~keep]] = False
    return coef[0] + coef[1] * xx + coef[2] * yy, noise


def measure(img, wcs, coord, pix_arcsec):
    """Centroid offset and peak S/N of the W-band emission at the star.

    v2 (2026-10-04): the background is a plane fitted to an 18-35" annulus and the noise is the clipped
    scatter there. The old version used a median beyond 40", which in structured fields (gradients,
    nebulosity) both biased the peak and inflated the noise, so real compact sources (e.g. J0719) failed.
    """
    ny, nx = img.shape
    x0, y0 = wcs.world_to_pixel(coord)
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    r = np.hypot(xx - x0, yy - y0) * pix_arcsec
    bg, noise = plane_background(img, xx, yy, (r > 18) & (r < 35))
    inner = r < 15
    sub = np.where(inner, img - bg, 0)
    sub[sub < 3 * noise] = 0
    if sub.sum() <= 0 or not noise > 0:
        return {"offset_arcsec": np.nan, "peak_snr": 0.0}
    cx, cy = (sub * xx).sum() / sub.sum(), (sub * yy).sum() / sub.sum()
    peak = np.nanmax((img - bg)[r < 6])
    return {"offset_arcsec": float(np.hypot(cx - x0, cy - y0) * pix_arcsec), "peak_snr": float(peak / noise)}


def run(names, out_csv="wise_centring_v2.csv"):
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
    rows = []
    for name in names:
        if name in set(s.designation):
            r = s[s.designation == name].iloc[0]
            c = SkyCoord(r.ra, r.dec, unit="deg")
        else:  # e.g. J0719 (TIC 106588577), not in the scored list: decode the J-name
            c = SkyCoord(f"{name[1:3]}h{name[3:5]}m{name[5:10]}s {name[10:13]}d{name[13:15]}m{name[15:]}s")
        try:
            imgs = SkyView.get_images(position=c, survey=BANDS, pixels=60, width=2 * u.arcmin, height=2 * u.arcmin)
        except Exception as e:
            print(name, "skyview failed", e)
            continue
        fig, axs = plt.subplots(1, 4, figsize=(12, 3.3))
        res = {"designation": name}
        for ax, band, hdul in zip(axs, BANDS, imgs):
            img = hdul[0].data.astype(float)
            w = WCS(hdul[0].header)
            m = measure(img, w, c, 120 / 60)
            res[f"{band.split()[1]}_offset"] = m["offset_arcsec"]
            res[f"{band.split()[1]}_snr"] = m["peak_snr"]
            lo, hi = np.nanpercentile(img, [5, 99.5])
            ax.imshow(img, origin="lower", cmap="magma", vmin=lo, vmax=hi)
            x0, y0 = w.world_to_pixel(c)
            ax.add_patch(plt.Circle((x0, y0), 6 / 2, fill=False, color="cyan", lw=1))
            ax.set_title(f"{band} µm  off={m['offset_arcsec']:.1f}\"  S/N={m['peak_snr']:.0f}", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{name} (2'x2', cyan = 12\" W4 beam)", fontsize=9)
        fig.tight_layout(); fig.savefig(OUT / f"{name}_wise.png", dpi=110); plt.close(fig)
        res["w3w4_centred"] = bool(np.nanmax([res["12_offset"], res["22_offset"]]) < 4)
        rows.append(res)
        print(res)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / out_csv, index=False)
    return df


if __name__ == "__main__":
    names = sys.argv[1:]
    if not names:
        s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
        names = s[s.tier.isin(["A", "B"])].designation.tolist()
    run(names)
