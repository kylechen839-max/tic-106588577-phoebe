# Spitzer 24 µm check of the tier-A W4 excesses (Kyle's PC, 2026-11-14)

Searched with IRSA SIA (`spitzer_seip`, `spitzer_sha`, 18″ radius) at all 12 tier-A positions. The
query was tested on Vega (50 SEIP images).

| Object | Spitzer coverage |
|---|---|
| **J2312+53** | **none** (0 images in SEIP or SHA), confirming the 11-14 cloud search |
| J0739, J2218, J0245, J1732, J2235, J1604, J1419, J0459, J0332 | none |
| **J2209+58** | SEIP super-mosaic (program 30002711) IRAC 4.5/5.8/8.0 + **MIPS 24 µm** |
| J2239+58 | inside the SEIP tile outlines and an SMOG strip (program 50398), but the star falls in **unobserved (blank) pixels** of every MIPS 24 mosaic. No usable 24 µm data |

SEIP outlines are tile polygons, larger than the area actually observed. That's why J2239 is "covered"
in the archive search but has no data.

## J2209+58 at 24 µm

- **SEIP source list** (VizieR II/368, via the Harvard mirror; Strasbourg VizieR was returning 404 and IRSA TAP
  502 on 11-14): `SSTSL2 J220912.54+582025.1`, 0.14″ from the star. MIPS 24 PSF flux 9.995 ± 0.131 mJy
  (S/N 76), flagged not extended (`24fEx` = 0).
- **Our aperture photometry** (`src/crossmatch/mips24_check.py`; r = 7″, sky 20–32″, aperture correction 1.61):
  9.83 ± 0.13 mJy, peak S/N 50. The two agree.
- **AllWISE W4** (22 µm, 12″ beam): 10.85 ± 0.93 mJy. **MIPS/W4 = 0.91**, the same within the W4 error.
- **By eye** (`mips24_cutouts.png`): a compact source centred on the star, on faint, patchy background.

**What this shows (verified):** at 24 µm with a 6″ beam, half the WISE W4 beam, essentially all of
the W4 flux is in a point source at the star. The excess is not extended nebulosity at ≳ 6″ scales.
This is the kind of check J2312 needs, but no such data exist for J2312.

**Size of the 24 µm excess (photosphere-dependent):**
- 2.9× when the photosphere is scaled from the pipeline's W4 ratio (3.39 mJy at 23.7 µm). That's ~16σ with a 4%
  MIPS calibration error added.
- 3.8× when it's the Rayleigh–Jeans extrapolation of IRAC 5.8 µm (2.56 mJy).
- IRAC 8.0 µm is already 1.3× the Rayleigh–Jeans extrapolation of 5.8 µm (30.0 vs 23.2 mJy), consistent with the W3 excess (1.28×).

**Caveat:** J2209's stellar temperature and distance are unreliable. Gaia parallax ≈ 0, RUWE 43, at b = +1.9°. A
point source at 6″ could still be unresolved dust ≲ 0.06 pc across at ≳ 2 kpc. So this confirms a
compact 24 µm excess, not a circumbinary disk.

Files: `mips24_photometry.csv`, `mips24_cutouts.png`, `sha_point.py` (rotated-frame photometry used to test
J2239), and `src/crossmatch/mips24_check.py`. The raw mosaics in `raw/` (8–17 MB) are not committed;
the URLs are in this README's history and in the IRSA SIA results.
