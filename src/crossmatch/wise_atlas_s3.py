"""10' AllWISE W3/W4 images and radial profiles from the IRSA AWS mirror (no IRSA/SkyView access needed).

  python src/crossmatch/wise_atlas_s3.py NAME RA DEC TILE OUTDIR
TILE is the AllWISE Atlas coadd_id covering the target, e.g. 0231p636_ac51.
Same near/far-sky aperture test as outputs/j2312_contamination/j2312_checks.py.
"""
import sys, urllib.request
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u

name, ra, dec, tile, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], Path(sys.argv[5])
out.mkdir(parents=True, exist_ok=True)
c = SkyCoord(ra, dec, unit="deg")
URL = "https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/allwise/images/p3am_cdd/{0}/{1}/{2}/{2}-{3}-int-3.fits"
PIX = 1.375

def aperture(img, x0, y0, r_ap, sky_in, sky_out):
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - x0, yy - y0) * PIX
    sky = img[(r >= sky_in) & (r < sky_out)]
    bg = np.nanmedian(sky)
    return np.nansum(img[r < r_ap] - bg), bg

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for ax, band, label in zip(axs, ("w3", "w4"), ("WISE 12 µm", "WISE 22 µm")):
    local = out / f"_{tile}-{band}.fits"
    if not local.exists():
        urllib.request.urlretrieve(URL.format(tile[:2], tile[:4], tile, band), local)
    with fits.open(local) as h:
        cut = Cutout2D(h[0].data.astype(float), c, 10 * u.arcmin, wcs=WCS(h[0].header), mode="partial")
    fits.PrimaryHDU(cut.data, header=cut.wcs.to_header()).writeto(out / f"{name}_{band.upper()}_10arcmin.fits", overwrite=True)
    img = cut.data
    x0, y0 = cut.wcs.world_to_pixel(c)
    near = aperture(img, x0, y0, 16.5, 50, 70)
    far = aperture(img, x0, y0, 16.5, 150, 200)
    ratio = near[0] / far[0]
    print(f"{band}: flux near-sky {near[0]:.1f} (sky {near[1]:.2f}), far-sky {far[0]:.1f} (sky {far[1]:.2f}), near/far {ratio:.3f}")
    lo, hi = np.nanpercentile(img, [5, 99.5])
    ax.imshow(img, origin="lower", cmap="magma", vmin=lo, vmax=hi)
    for rr, st in ((16.5, "-"), (50, ":"), (70, ":"), (150, ":"), (200, ":")):
        ax.add_patch(plt.Circle((x0, y0), rr / PIX, fill=False, color="cyan" if st == "-" else "w", ls=st, lw=1 if st == "-" else 0.6))
    ax.set_title(f"{label}, near/far-sky flux {ratio:.2f}", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    # radial profile
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - x0, yy - y0) * PIX
    edges = np.arange(0, 205, 5)
    sky = np.nanmedian(img[(r > 250)]) if np.any(r > 250) else np.nanpercentile(img, 20)
    prof = np.array([np.nanmedian(img[(r >= a) & (r < b)]) for a, b in zip(edges[:-1], edges[1:])]) - sky
    axs[2].plot(edges[:-1] + 2.5, prof / prof[0], label=label)
    local.unlink()
axs[2].axhline(0, color="k", lw=0.5); axs[2].set_xlabel("radius (arcsec)"); axs[2].set_ylabel("median flux / peak (sky-subtracted)")
axs[2].legend(); axs[2].set_title("Radial profiles")
fig.suptitle(f"{name}: AllWISE 10'x10' (N up, E left; cyan = 16.5\" aperture; dotted = sky annuli)", fontsize=11)
fig.tight_layout(); fig.savefig(out / f"{name}_wise_10arcmin.png", dpi=110)
