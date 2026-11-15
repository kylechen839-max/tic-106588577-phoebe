"""J2312-style contamination checks for one tier-A candidate, run where IRSA, MAST and Gaia are reachable.

  python src/crossmatch/contamination_checks.py DESIGNATION OUTDIR

Same method as outputs/j2312_contamination/ (j2312_checks.py + w4_bgsub.py), with the 10' W3/W4 images
cut from the AllWISE Atlas tiles on the IRSA AWS mirror (SkyView returns no WISE images as of 2026-11-11):
  (a) AllWISE sources within 3'; 16.5" aperture flux with near (50-70") vs far (150-200") sky
  (b) TIC v8.2 within 60" (contamination ratio), Gaia DR3 neighbours within 60"
  (c) W3/W4 excess re-scaled for four background estimates (far, 50-70", 25-40", halo extrapolated in)
Images are written as PNG/FITS so they can be checked by eye before any verdict.
astropy.wcs is not used (Windows Application Control blocks its compiled library on Kyle's PC); the Atlas
tiles are RA---SIN/DEC--SIN with no rotation, so sky-to-pixel is done directly below.
"""
import json
import re
import sys
import urllib.request
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.gaia import Gaia
from astroquery.ipac.irsa import Irsa
from astroquery.mast import Catalogs

ROOT = Path(__file__).resolve().parents[2]
sys.stdout.reconfigure(encoding="utf-8")
NAME, OUT = sys.argv[1], Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
s = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv")
row = s[s.designation == NAME].iloc[0]
c = SkyCoord(row.ra, row.dec, unit="deg")
BUCKET = "https://nasa-irsa-wise.s3.us-west-2.amazonaws.com"
PREFIX = "wise/allwise/images/p3am_cdd/"
PIX = 1.375  # arcsec per Atlas pixel


def list_prefixes(prefix):
    url = f"{BUCKET}/?list-type=2&delimiter=/&prefix={prefix}"
    xml = urllib.request.urlopen(url, timeout=60).read().decode()
    return re.findall(r"<Prefix>([^<]+)</Prefix>", xml)[1:]


def nearest_tile():
    best = (9e9, None)
    for ra in (row.ra - 1.6, row.ra, row.ra + 1.6):
        p2 = f"{int((ra % 360) * 10):04d}"[:2]
        for sub in list_prefixes(f"{PREFIX}{p2}/"):
            for q in list_prefixes(sub):
                n = q.rstrip("/").split("/")[-1]
                tra, tdec = int(n[:4]) / 10, int(n[5:8]) / 10 * (1 if n[4] == "p" else -1)
                sep = SkyCoord(tra, tdec, unit="deg").separation(c).deg
                best = min(best, (sep, n))
    return best[1]


def sin_world_to_pixel(hdr, ra, dec):
    """0-based pixel (x, y) for a SIN-projected, unrotated FITS image."""
    a0, d0 = np.radians(hdr["CRVAL1"]), np.radians(hdr["CRVAL2"])
    a, d = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    x = np.degrees(np.cos(d) * np.sin(a - a0))
    y = np.degrees(np.sin(d) * np.cos(d0) - np.cos(d) * np.sin(d0) * np.cos(a - a0))
    return hdr["CRPIX1"] - 1 + x / hdr["CDELT1"], hdr["CRPIX2"] - 1 + y / hdr["CDELT2"]


def cutout(data, hdr, size_arcsec):
    """Square cutout centred on the target; returns (image, x0, y0, offset, header)."""
    xc, yc = sin_world_to_pixel(hdr, c.ra.deg, c.dec.deg)
    half = int(round(size_arcsec / 2 / PIX))
    ix, iy = int(round(float(xc))), int(round(float(yc)))
    x1, y1 = max(ix - half, 0), max(iy - half, 0)
    img = data[y1:iy + half + 1, x1:ix + half + 1].astype(float)
    h2 = fits.Header({k: hdr[k] for k in ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "BUNIT")})
    h2["CRPIX1"], h2["CRPIX2"] = hdr["CRPIX1"] - x1, hdr["CRPIX2"] - y1
    return img, float(xc) - x1, float(yc) - y1, (x1, y1), h2


def aperture(img, x0, y0, r_ap, sky_in, sky_out):
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - x0, yy - y0) * PIX
    sky = img[(r >= sky_in) & (r < sky_out)]
    bg = np.nanmedian(sky)
    return np.nansum(img[r < r_ap] - bg), bg


log = {"designation": NAME}
# (a) AllWISE 3' cone
full = Irsa.query_region(c, catalog="allwise_p3as_psd", spatial="Cone", radius=3 * u.arcmin).to_pandas()
full["sep_arcsec"] = c.separation(SkyCoord(full.ra, full.dec, unit="deg")).arcsec
cols = ["designation", "ra", "dec", "sep_arcsec", "w1mpro", "w2mpro", "w3mpro", "w3sigmpro", "w3snr", "w4mpro",
        "w4sigmpro", "w4snr", "w4rchi2", "ext_flg", "nb", "cc_flags", "ph_qual"]
cat = full.sort_values("sep_arcsec")[cols]
cat.to_csv(OUT / f"{NAME}_allwise_3arcmin.csv", index=False)
tgt = cat.iloc[0]
nb = cat[(cat.sep_arcsec > 5)]
print("== AllWISE target:", tgt[["sep_arcsec", "w3mpro", "w3snr", "w4mpro", "w4snr", "w4rchi2", "ext_flg", "cc_flags", "ph_qual"]].to_dict())
print("== AllWISE neighbours within 60\" with W4 S/N >= 3:")
print(nb[(nb.sep_arcsec < 60) & (nb.w4snr >= 3)][["designation", "sep_arcsec", "w1mpro", "w3mpro", "w3snr", "w4mpro", "w4snr"]].to_string(index=False))
log["allwise_target"] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else str(v)) for k, v in tgt.items()}
log["n_w4_neighbours_lt60"] = int(((nb.sep_arcsec < 60) & (nb.w4snr >= 3)).sum())

# 10' W3/W4 images
tile = nearest_tile()
print("== Atlas tile", tile)
log["tile"] = tile
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
prof_out, bg_rows = {}, []
for ax, band, label in zip(axs, ("w3", "w4"), ("WISE 12 µm", "WISE 22 µm")):
    tmp = OUT / f"_{tile}-{band}.fits"
    if not tmp.exists():
        urllib.request.urlretrieve(f"{BUCKET}/{PREFIX}{tile[:2]}/{tile[:4]}/{tile}/{tile}-{band}-int-3.fits", tmp)
    with fits.open(tmp) as h:
        hdr = h[0].header
        img, x0, y0, (ox, oy), h2 = cutout(h[0].data, hdr, 600)
    tmp.unlink()
    fits.PrimaryHDU(img, header=h2).writeto(OUT / f"{NAME}_{band.upper()}_10arcmin.fits", overwrite=True)
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - x0, yy - y0) * PIX
    sky_far = np.nanmedian(img[(r >= 150) & (r < 200)])
    sky_near = np.nanmedian(img[(r >= 50) & (r < 70)])
    sky_in = np.nanmedian(img[(r >= 25) & (r < 40)])
    edges = np.arange(0, 202, 4)
    rc = edges + 2
    prof = np.array([np.nanmedian(img[(r >= a) & (r < a + 4)]) for a in edges])
    m = (rc >= 25) & (rc <= 70)
    p = np.polyfit(rc[m], prof[m], 1)
    sky_lin = np.polyval(p, 16.5 / 2)
    Fnear = np.nansum(img[r < 16.5] - sky_near)
    for lab, bg in (("far 150-200\"", sky_far), ("near 50-70\" (AllWISE)", sky_near), ("inner 25-40\"", sky_in),
                    ("halo extrapolated in", sky_lin)):
        bg_rows.append(dict(band=band.upper(), background=lab, sky=bg, flux_over_allwise=np.nansum(img[r < 16.5] - bg) / Fnear))
    q = [np.nanmedian(img[(r >= 50) & (r < 70) & (np.sign(xx - x0) == sx) & (np.sign(yy - y0) == sy)]) -
         np.nanmedian(img[(r >= 150) & (r < 200) & (np.sign(xx - x0) == sx) & (np.sign(yy - y0) == sy)])
         for sx in (-1, 1) for sy in (-1, 1)]
    near_far = aperture(img, x0, y0, 16.5, 50, 70)[0] / aperture(img, x0, y0, 16.5, 150, 200)[0]
    print(f"== {band.upper()}: near/far flux {near_far:.3f}; near-minus-far sky by quadrant {np.round(q, 3)}")
    log[f"{band}_near_over_far"] = float(near_far)
    log[f"{band}_quadrant_near_minus_far"] = [float(v) for v in q]
    lo, hi = np.nanpercentile(img, [5, 99.5])
    ax.imshow(img, origin="lower", cmap="magma", vmin=lo, vmax=hi)
    for rr, st in ((16.5, "-"), (50, ":"), (70, ":"), (150, ":"), (200, ":")):
        ax.add_patch(plt.Circle((x0, y0), rr / PIX, fill=False, color="cyan" if st == "-" else "w", ls=st,
                                lw=1 if st == "-" else 0.6))
    ax.add_patch(plt.Rectangle((x0 - 60 / PIX, y0 - 60 / PIX), 120 / PIX, 120 / PIX, fill=False, color="cyan", lw=0.6))
    w4n = full[(full.w4snr >= 3) & (full.sep_arcsec > 5)]
    if len(w4n):
        sx, sy = sin_world_to_pixel(hdr, w4n.ra.values, w4n.dec.values)
        ax.scatter(sx - ox, sy - oy, s=40, facecolors="none", edgecolors="lime", lw=0.8)
    ax.set_title(f"{label}: near/far-sky flux {near_far:.2f}", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    sky0 = np.nanpercentile(img, 20)
    axs[2].plot(rc, (prof - sky0) / (prof[0] - sky0), label=label)
axs[2].set_yscale("symlog", linthresh=0.01); axs[2].axhline(0, color="k", lw=0.5)
axs[2].set_xlabel("radius (arcsec)"); axs[2].set_ylabel("median flux / centre (sky = 20th pct)"); axs[2].legend()
fig.suptitle(f"{NAME}: 10'×10' AllWISE Atlas (N up, E left; cyan = 16.5\" aperture, 2' box; dotted = sky annuli; lime = W4 S/N≥3 sources)", fontsize=10)
fig.tight_layout(); fig.savefig(OUT / f"{NAME}_wise_10arcmin.png", dpi=110); plt.close(fig)

# (c) excess against each background, using the pipeline's photosphere ratios
bg = pd.DataFrame(bg_rows)
e4 = max(float(tgt.w4sigmpro) if pd.notna(tgt.w4sigmpro) else 1.0857 / tgt.w4snr, 0.03)
e3 = max(float(tgt.w3sigmpro) if pd.notna(tgt.w3sigmpro) else 1.0857 / tgt.w3snr, 0.03)
out = []
for _, b in bg.iterrows():
    e = e4 if b.band == "W4" else e3
    for lab, R0 in (("free", row[f"{b.band}_ratio"]), ("fixedTeff", row.get(f"fixed_{b.band}_ratio", np.nan))):
        if pd.notna(R0):
            R = R0 * b.flux_over_allwise
            out.append(dict(band=b.band, background=b.background, flux_over_allwise=b.flux_over_allwise, photosphere=lab,
                            excess_ratio=R, excess_sigma=(1 - 1 / R) / (0.4 * np.log(10) * e) if R > 0 else np.nan))
ex = pd.DataFrame(out)
ex.to_csv(OUT / f"{NAME}_bgsub_excess.csv", index=False)
print("== Excess vs background estimate (photometric sigma only):")
print(ex.round(3).to_string(index=False))
w4 = ex[(ex.band == "W4") & (ex.photosphere == "free")].set_index("background")
spread = 2.5 * np.log10(w4.flux_over_allwise.max() / w4.flux_over_allwise.min()) / 2
sig_sys = np.hypot(0.4 * np.log(10) * e4, 0.4 * np.log(10) * spread)
for k in w4.index:
    log[f"W4_excess_{k}"] = [float(w4.loc[k, "excess_ratio"]), float(w4.loc[k, "excess_sigma"]),
                             float((1 - 1 / w4.loc[k, "excess_ratio"]) / sig_sys)]
print(f"== W4 background spread ±{spread:.2f} mag; sigma incl. that systematic: " +
      ", ".join(f"{k}: {(1 - 1 / v) / sig_sys:.1f}" for k, v in w4.excess_ratio.items()))

# (b) TIC and Gaia neighbours
tic = Catalogs.query_region(c, radius=60 * u.arcsec, catalog="TIC").to_pandas().sort_values("dstArcSec")
tcols = ["ID", "dstArcSec", "Tmag", "GAIAmag", "Teff", "rad", "contratio", "disposition", "objType"]
tic[tcols].to_csv(OUT / f"{NAME}_tic_60arcsec.csv", index=False)
print("== TIC within 60\" (first 6):"); print(tic[tcols].head(6).to_string(index=False))
log["tic_contratio"] = None if pd.isna(tic.contratio.iloc[0]) else float(tic.contratio.iloc[0])
job = Gaia.launch_job_async(
    "SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, parallax_over_error, ruwe, "
    f"DISTANCE(POINT(ra, dec), POINT({c.ra.deg}, {c.dec.deg}))*3600 AS sep FROM gaiadr3.gaia_source "
    f"WHERE 1=CONTAINS(POINT(ra, dec), CIRCLE({c.ra.deg}, {c.dec.deg}, {60 / 3600})) ORDER BY sep")
g = job.get_results().to_pandas()
g.to_csv(OUT / f"{NAME}_gaiadr3_60arcsec.csv", index=False)
G0 = g.phot_g_mean_mag.iloc[0]
for lo_, hi_ in ((0, 12), (12, 60)):
    sub = g[(g.sep > max(lo_, 1)) & (g.sep <= hi_)]
    frac = float((10 ** (-0.4 * (sub.phot_g_mean_mag - G0))).sum())
    nbright = int((sub.phot_g_mean_mag < G0 + 5).sum())
    print(f"== Gaia {lo_}-{hi_}\": n={len(sub)}, within 5 mag of target: {nbright}, summed G flux / target = {frac:.4f}")
    log[f"gaia_{lo_}_{hi_}"] = dict(n=len(sub), n_within_5mag=nbright, flux_frac=frac,
                                    brightest_G=float(sub.phot_g_mean_mag.min()) if len(sub) else None)
(OUT / f"{NAME}_summary.json").write_text(json.dumps(log, indent=1, default=str))
