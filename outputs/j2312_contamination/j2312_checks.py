"""J2312 contamination checks: same method as src/crossmatch/w4_neighbours.py (081ab9c)
plus TIC v8.2 contamination ratio and Gaia DR3 neighbours within 60"."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa
from astroquery.skyview import SkyView
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
NAME = "J231201.40+532028.6"
c = SkyCoord(348.0058494, 53.3413022, unit="deg")  # from dd1_candidates_scored.csv
BANDS = ["WISE 12", "WISE 22"]
PIX = 600 / 300


def aperture(img, wcs, coord, r_ap, sky_in, sky_out):
    ny, nx = img.shape
    x0, y0 = wcs.world_to_pixel(coord)
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x0, yy - y0) * PIX
    sky = img[(r >= sky_in) & (r < sky_out)]
    bg = np.nanmedian(sky)
    flux = np.nansum(img[r < r_ap] - bg)
    return flux, bg, 1.4826 * np.nanmedian(np.abs(sky - bg))


# (a) AllWISE 3' cone
full = Irsa.query_region(c, catalog="allwise_p3as_psd", spatial="Cone", radius=3 * u.arcmin).to_pandas()
full["sep_arcsec"] = c.separation(SkyCoord(full.ra, full.dec, unit="deg")).arcsec
cols = ["designation", "ra", "dec", "sep_arcsec", "w1mpro", "w2mpro", "w3mpro", "w3snr", "w4mpro", "w4snr",
        "w4rchi2", "ext_flg", "nb", "cc_flags"]
cat = full.sort_values("sep_arcsec")[cols]
cat.to_csv(OUT / f"{NAME}_allwise_3arcmin.csv", index=False)
print("== AllWISE within 3'")
print(cat.drop(columns=["ra", "dec"]).to_string(index=False))

# (a) 10' W3/W4 images, near vs far sky
imgs = SkyView.get_images(position=c, survey=BANDS, pixels=300, width=10 * u.arcmin, height=10 * u.arcmin)
fig, axs = plt.subplots(1, 2, figsize=(10, 5.2))
for ax, band, hdul in zip(axs, BANDS, imgs):
    hdul.writeto(OUT / f"{NAME}_{band.replace(' ', '')}_10arcmin.fits", overwrite=True)
    img = hdul[0].data.astype(float)
    w = WCS(hdul[0].header)
    near = aperture(img, w, c, 16.5, 50, 70)
    far = aperture(img, w, c, 16.5, 150, 200)
    print(f"== {band}: flux near-sky {near[0]:.1f} (sky {near[1]:.3f}, rms {near[2]:.3f}), "
          f"far-sky {far[0]:.1f} (sky {far[1]:.3f}, rms {far[2]:.3f}), near/far {near[0] / far[0]:.3f}")
    lo, hi = np.nanpercentile(img, [5, 99.5])
    ax.imshow(img, origin="lower", cmap="magma", vmin=lo, vmax=hi)
    x0, y0 = w.world_to_pixel(c)
    ax.add_patch(plt.Rectangle((x0 - 60 / PIX, y0 - 60 / PIX), 120 / PIX, 120 / PIX, fill=False, color="cyan", lw=0.8))
    ax.add_patch(plt.Circle((x0, y0), 16.5 / PIX, fill=False, color="cyan", lw=0.8))
    for rr in (50, 70, 150, 200):
        ax.add_patch(plt.Circle((x0, y0), rr / PIX, fill=False, color="white", lw=0.4, ls=":"))
    nb = full[(full.w4snr >= 3) & (full.sep_arcsec > 5)]
    if len(nb):
        sx, sy = w.world_to_pixel(SkyCoord(nb.ra, nb.dec, unit="deg"))
        ax.scatter(sx, sy, s=40, facecolors="none", edgecolors="lime", lw=0.8)
    ax.set_title(f"{band} µm, near/far-sky flux {near[0] / far[0]:.2f}", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f"{NAME} 10'x10' (N up, E left; cyan = 2' box + 16.5\" ap; dotted = sky annuli; lime = W4 SNR>=3)", fontsize=8)
fig.tight_layout(); fig.savefig(OUT / f"{NAME}_wise_10arcmin.png", dpi=120); plt.close(fig)

# (b) TIC v8.2
tic = Catalogs.query_region(c, radius=60 * u.arcsec, catalog="TIC").to_pandas()
tic = tic.sort_values("dstArcSec")
tcols = ["ID", "dstArcSec", "Tmag", "GAIAmag", "Teff", "rad", "contratio", "disposition", "objType"]
tic[tcols].to_csv(OUT / f"{NAME}_tic_60arcsec.csv", index=False)
print("== TIC within 60\"")
print(tic[tcols].to_string(index=False))

# (b) Gaia DR3 within 60"
job = Gaia.launch_job_async(
    "SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, parallax_over_error, ruwe, "
    "DISTANCE(POINT(ra, dec), POINT(%f, %f))*3600 AS sep FROM gaiadr3.gaia_source "
    "WHERE 1=CONTAINS(POINT(ra, dec), CIRCLE(%f, %f, %f)) ORDER BY sep" % (c.ra.deg, c.dec.deg, c.ra.deg, c.dec.deg, 60 / 3600))
g = job.get_results().to_pandas()
g.to_csv(OUT / f"{NAME}_gaiadr3_60arcsec.csv", index=False)
print("== Gaia DR3 within 60\"")
print(g.to_string(index=False))
G0 = g.phot_g_mean_mag.iloc[0]
for lo_, hi_ in ((0, 12), (12, 60)):
    s = g[(g.sep > max(lo_, 1)) & (g.sep <= hi_)]
    frac = (10 ** (-0.4 * (s.phot_g_mean_mag - G0))).sum()
    print(f"Gaia {lo_}-{hi_}\": n={len(s)}, G<{G0 + 5:.1f}: {(s.phot_g_mean_mag < G0 + 5).sum()}, "
          f"summed G flux / target = {frac:.4f}")
