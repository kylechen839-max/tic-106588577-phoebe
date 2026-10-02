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


def measure(img, wcs, coord, pix_arcsec):
    ny, nx = img.shape
    x0, y0 = wcs.world_to_pixel(coord)
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x0, yy - y0) * pix_arcsec
    bg = np.nanmedian(img[r > 40])
    noise = 1.4826 * np.nanmedian(np.abs(img[r > 40] - bg))
    inner = r < 15
    sub = np.where(inner, img - bg, 0)
    sub[sub < 3 * noise] = 0
    if sub.sum() <= 0:
        return {"offset_arcsec": np.nan, "peak_snr": 0.0}
    cx, cy = (sub * xx).sum() / sub.sum(), (sub * yy).sum() / sub.sum()
    peak = np.nanmax(img[r < 6]) - bg
    return {"offset_arcsec": float(np.hypot(cx - x0, cy - y0) * pix_arcsec), "peak_snr": float(peak / noise)}


def run(names):
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
    rows = []
    for name in names:
        r = s[s.designation == name].iloc[0]
        c = SkyCoord(r.ra, r.dec, unit="deg")
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
    df.to_csv(OUT / "wise_centring.csv", index=False)
    return df


if __name__ == "__main__":
    names = sys.argv[1:]
    if not names:
        s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
        names = s[s.tier == "A"].designation.tolist()
    run(names)
