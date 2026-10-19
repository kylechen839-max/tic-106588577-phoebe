"""Wide-field W3/W4 check for nebulosity or neighbours just outside the 2' cutout.

For one candidate: fetches 10'x10' AllWISE atlas images (SkyView), lists AllWISE sources
within 3' (IRSA), and measures W3/W4 aperture photometry at the star with a near sky
annulus (50-70", as AllWISE) and a far one (150-200"). If the near-annulus W4 flux is
well above the far-annulus flux, structured background is leaking into the photometry.
Needs IRSA and SkyView access (Kyle's Mac; the cloud container is blocked).

Usage: .venv-tess/bin/python src/crossmatch/w4_neighbours.py J231201.40+532028.6
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
from astroquery.ipac.irsa import Irsa
from astroquery.skyview import SkyView

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch" / "wise_cutouts"
BANDS = ["WISE 12", "WISE 22"]
PIX = 600 / 300  # arcsec per pixel for a 10' image at 300 px


def aperture(img, wcs, coord, r_ap, sky_in, sky_out):
    ny, nx = img.shape
    x0, y0 = wcs.world_to_pixel(coord)
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x0, yy - y0) * PIX
    sky = img[(r >= sky_in) & (r < sky_out)]
    bg = np.nanmedian(sky)
    flux = np.nansum(img[r < r_ap] - bg)
    return flux, bg, 1.4826 * np.nanmedian(np.abs(sky - bg))


def run(name):
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
    r = s[s.designation == name].iloc[0]
    c = SkyCoord(r.ra, r.dec, unit="deg")

    full = Irsa.query_region(c, catalog="allwise_p3as_psd", spatial="Cone", radius=3 * u.arcmin).to_pandas()
    full["sep_arcsec"] = c.separation(SkyCoord(full.ra, full.dec, unit="deg")).arcsec
    cols = ["designation", "sep_arcsec", "w1mpro", "w3mpro", "w3snr", "w4mpro", "w4snr", "w4rchi2", "ext_flg", "nb", "cc_flags"]
    cat = full.sort_values("sep_arcsec")[cols]
    cat.to_csv(OUT / f"{name}_allwise_3arcmin.csv", index=False)
    print(cat.to_string(index=False))

    imgs = SkyView.get_images(position=c, survey=BANDS, pixels=300, width=10 * u.arcmin, height=10 * u.arcmin)
    fig, axs = plt.subplots(1, 2, figsize=(9, 4.6))
    for ax, band, hdul in zip(axs, BANDS, imgs):
        img = hdul[0].data.astype(float)
        w = WCS(hdul[0].header)
        near = aperture(img, w, c, 16.5, 50, 70)
        far = aperture(img, w, c, 16.5, 150, 200)
        print(f"{band}: flux near-sky {near[0]:.1f} (sky {near[1]:.3f}), far-sky {far[0]:.1f} (sky {far[1]:.3f}), "
              f"near/far {near[0] / far[0]:.2f}")
        lo, hi = np.nanpercentile(img, [5, 99.5])
        ax.imshow(img, origin="lower", cmap="magma", vmin=lo, vmax=hi)
        x0, y0 = w.world_to_pixel(c)
        ax.add_patch(plt.Rectangle((x0 - 60 / PIX, y0 - 60 / PIX), 120 / PIX, 120 / PIX, fill=False, color="cyan", lw=0.8))
        ax.add_patch(plt.Circle((x0, y0), 16.5 / PIX, fill=False, color="cyan", lw=0.8))
        nb = cat[(cat.w4snr >= 3) & (cat.sep_arcsec > 5)]
        if len(nb):
            pos = SkyCoord(*[full.loc[nb.index, k] for k in ("ra", "dec")], unit="deg")
            sx, sy = w.world_to_pixel(pos)
            ax.scatter(sx, sy, s=40, facecolors="none", edgecolors="lime", lw=0.8)
        ax.set_title(f"{band} µm, near/far-sky flux {near[0] / far[0]:.2f}", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"{name} (10'x10', cyan = 2' cutout and 16.5\" aperture)", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / f"{name}_wise_10arcmin.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "J231201.40+532028.6")
