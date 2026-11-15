"""Spitzer MIPS 24 um check of WISE W4 excesses, from SEIP super-mosaics (IRSA, program 30002711).

  python src/crossmatch/mips24_check.py OUTDIR

MIPS 24 has a 6" PSF (WISE W4: 12"), so a compact excess should appear as a point source with the
W4-like flux, while nebular W4 emission spreads out. Aperture photometry: r = 7", sky 20-32",
aperture correction 1.61 (MIPS Instrument Handbook, 24 um point source). Mosaics are RA---TAN with no
rotation; sky-to-pixel is done directly because astropy.wcs is blocked on this PC.
"""
import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1])
s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv").set_index("designation")
TARGETS = ["J220912.52+582025.0"]  # J2239 lies in unobserved (NaN) parts of both SEIP tiles and the SMOG strip
W4 = {  # AllWISE W4 flux (mJy) from the SEIP source list's matched AllWISE row (II/368)
    "J220912.52+582025.0": (10.854, 0.930), "J223949.47+583254.4": (11.556, 1.032)}


def tan_world_to_pixel(h, ra, dec):
    a0, d0 = np.radians(h["CRVAL1"]), np.radians(h["CRVAL2"])
    a, d = np.radians(ra), np.radians(dec)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    x = np.degrees(np.cos(d) * np.sin(a - a0) / cosc)
    y = np.degrees((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc)
    return h["CRPIX1"] - 1 + x / h["CD1_1"], h["CRPIX2"] - 1 + y / h["CD2_2"]


rows = []
fig, axs = plt.subplots(len(TARGETS), 2, figsize=(10, 5 * len(TARGETS)), squeeze=False)
for k, name in enumerate(TARGETS):
    ra, dec = s.loc[name, "ra"], s.loc[name, "dec"]
    best = None
    for f in sorted(glob.glob(str(OUT / "raw" / "*.MIPS.1.mosaic.fits"))):
        h = fits.getheader(f)
        x, y = tan_world_to_pixel(h, ra, dec)
        if 20 < x < h["NAXIS1"] - 20 and 20 < y < h["NAXIS2"] - 20:
            img = fits.getdata(f).astype(float)
            if np.isfinite(img[int(round(y)), int(round(x))]):
                cov = np.isfinite(img[int(y) - 20:int(y) + 20, int(x) - 20:int(x) + 20]).mean()
                if best is None or cov > best[0]:
                    best = (cov, f, h, img, x, y)
    if best is None:
        print(name, "not covered by downloaded mosaics"); continue
    _, f, h, img, x0, y0 = best
    pix = abs(h["CD1_1"]) * 3600
    sr = (pix / 206264.806) ** 2
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - x0, yy - y0) * pix
    sky = img[(r >= 20) & (r < 32)]
    bg = np.nanmedian(sky)
    sig = 1.4826 * np.nanmedian(np.abs(sky - bg))
    ap = r < 7
    flux = np.nansum(img[ap] - bg) * sr * 1e9 * 1.61  # MJy/sr -> mJy, aperture-corrected
    err = sig * np.sqrt(ap.sum()) * sr * 1e9 * 1.61
    peak = (np.nanmax(img[r < 5]) - bg) / sig
    big = r < 60
    neb = np.nanmedian(img[(r > 40) & (r < 60)]), np.nanpercentile(img[big], 95)
    w4, ew4 = W4[name]
    rows.append(dict(designation=name, mosaic=Path(f).name, flux24_mJy=flux, err_mJy=err, peak_snr=peak,
                     sky_MJysr=bg, sky_rms=sig, wise_w4_mJy=w4, mips_over_w4=flux / w4))
    print(f"{name} [{Path(f).name}]: MIPS24 r=7\" flux {flux:.2f} ± {err:.2f} mJy (peak S/N {peak:.1f}); "
          f"sky {bg:.2f} ± {sig:.3f} MJy/sr; WISE W4 {w4:.2f} ± {ew4:.2f} mJy -> MIPS/W4 = {flux / w4:.2f}")
    for j, half in enumerate((60, 150)):
        ax = axs[k, j]
        n = int(half / pix)
        sub = img[int(y0) - n:int(y0) + n + 1, int(x0) - n:int(x0) + n + 1]
        lo, hi = np.nanpercentile(sub, [5, 99.7])
        ax.imshow(sub, origin="lower", cmap="magma", vmin=lo, vmax=hi,
                  extent=[-n * pix, n * pix, -n * pix, n * pix])
        ax.add_patch(plt.Circle((x0 - int(x0) - 0.5 * 0, y0 - int(y0)), 7, fill=False, color="cyan", lw=1))
        ax.add_patch(plt.Circle((0, 0), 6, fill=False, color="cyan", lw=1))
        ax.add_patch(plt.Circle((0, 0), 12, fill=False, color="w", lw=0.8, ls=":"))
        ax.set_title(f"{name[:5]}{name[10:13]}  MIPS 24 µm, {2 * half}\" box\ncyan r=6\" (MIPS beam), dotted r=12\" (W4 beam)", fontsize=9)
        ax.set_xlabel("arcsec E-W"); ax.set_ylabel("arcsec N-S")
fig.tight_layout(); fig.savefig(OUT / "mips24_cutouts.png", dpi=110)
pd.DataFrame(rows).to_csv(OUT / "mips24_photometry.csv", index=False)
