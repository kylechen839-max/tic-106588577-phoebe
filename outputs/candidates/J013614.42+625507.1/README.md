# J013614.42+625507.1 (J0136+62), TIC 390117626

DD 1.0 object; tier B on `mac-run-2026-10-04` (569ba9d), where the TESS 10k EB catalogue (Kostov 2025) lists it with P = 2.78951 d.
Built 2026-10-26 in the cloud from the public AWS mirrors, because MAST and IRSA are blocked there.

## TESS light curve (`src/pipeline/fetch_tess_ffi_s3.py`)
- 15×15 px FFI cutouts from the TESScut cubes on s3://stpubdata: sectors 18, 24, 25 (30-min) and 58, 85 (200-s).
  No SPOC 2-min or QLP/TESS-SPOC products are on the mirror, so this is our own aperture photometry (5–7 px aperture, per-cadence background).
- Cadences with scattered light (background > 8% of the star's flux above its quiet level) are dropped; per-orbit quadratic trends are removed while fitting a 4-harmonic phase model, so the periodic signal is kept.
  Then the usual `fetch_tess_lc.py` clipping, period refinement and sector-consensus binning.
- Result: P = 2.789497 d, T0 = BTJD 2900.40752; 20,493 points.
  - Primary eclipse ~2.5–3% deep, half-width ~0.06 in phase; odd/even depths agree (0.3σ), so P is not doubled.
  - No secondary eclipse detected.
  - Out-of-eclipse variation of ~1.3% peak-to-peak, once per orbit, minimum near phase 0.5. Its cause is not established.
- Per-sector primary depth: S18 0.023, S24 0.025, S25 0.034, S58 0.029, S85 0.028. The spread suggests different dilution per sector, since the aperture differs and the field is crowded.
- Difference images (out-of-eclipse minus in-eclipse, S58 and S85) centre within 0.2 px of the target, so the eclipse is on this star or within ~1 TESS pixel (21″) of it.

## WISE (`src/crossmatch/wise_atlas_s3.py`)
- 10′ AllWISE W3/W4 from the IRSA AWS mirror (tile 0231p636_ac51).
  - The star sits on the northern edge of extended 12/22 µm nebulosity that runs south of it.
  - Near/far sky flux ratio: 0.84 (W3), 0.88 (W4).
- The 2′ cutout (from mac-run) gives image S/N of 4 (W3) and 5 (W4) at the star's position, so the W3/W4 excess (7.9σ / 14.5σ against the photosphere) may be largely nebular. Same caveat as J2312.

## Files
- `_raw.txt`, `_binned.txt`, `_ephemeris.json`, `_lc_diagnostic.png`
- `_ffi_stamps.png`, `_diff_image.png`
- `_wise_10arcmin.png`, `_W3/W4_10arcmin.fits`, `_wise_2arcmin.png`
- The FFI cutouts (`ffi_cutouts/`, 120 MB) are not committed; rerun the script to rebuild them.

## Excess after subtracting the local background (2026-10-27)
Same method as `outputs/j2312_contamination/README.md`. Aperture flux relative to the 50–70″ sky: W4 0.90 / 0.72 / 0.56 and W3 0.85 / 0.74 / 0.71 for the 35–50″, 25–35″ and 20–30″ annuli.

Against the blackbody photosphere (W4 62.7×, 14.4σ; W3 2.28×, 7.2σ):
- 25–35″ sky: W4 45× (10.3σ), W3 1.68× (3.8σ).
- 20–30″ sky: W4 35× (8.0σ), W3 1.61× (3.4σ).
