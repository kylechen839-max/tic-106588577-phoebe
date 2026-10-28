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

## Clean light curve (2026-10-28)
Per-sector quality from `_raw.txt`. Noise is point-to-point, and "30-min equiv." scales it to 30-min cadence.

| Sector | Cadence | Eclipses | Noise (30-min equiv.) | Primary depth (±0.02 phase) | Eclipse S/N |
|---|---|---:|---:|---:|---:|
| 18 | 30 min | 7 | 3.1 ppt | 1.54% | 30 |
| 24 | 30 min | 6 | 2.9 ppt | 2.06% | 39 |
| 25 | 30 min | 6 | 3.1 ppt | 2.51% | 46 |
| 58 | 200 s | 10 | 2.1 ppt | 2.25% | 78 |
| 85 | 200 s | 6 | 2.1 ppt | 2.22% | 61 |

The 30-min sectors disagree in depth (1.5–2.5%, likely dilution from different apertures in a crowded field), so a five-sector stack would mix depths.
- **Best single sector:** S58, which has the most eclipses and the lowest noise.
- S58 and S85 agree in depth (2.25% vs 2.22%), so they are also stacked.

Files:
- `_clean_S58_timeseries.txt`;
- `_clean_S58_phasebinned.txt` and `_clean_S58_S85_phasebinned.txt` (200 phase bins, median error 1.1 and 0.87 ppt);
- `_clean_lc.png`.

The ~1.5% out-of-eclipse wave (maximum near phase 0.17, minimum near 0.5) is the same in both sectors. No secondary eclipse is visible at the 1 ppt level.

## Eclipse parameters (2026-10-28, from the S58 + S85 clean light curve)

| Parameter | Value |
|---|---|
| Period | 2.789497 d |
| T0 (primary) | BTJD 2900.40752 |
| Primary depth | 2.7% at the core; 2.32 ± 0.03% mean over \|φ\| < 0.0175 (relative to a local quadratic baseline) |
| Primary duration | ~0.07 in phase (~4.7 h) |
| Secondary depth at φ 0.5 | 0.10 ± 0.03% formal (3.6σ). **Not significant**: the same local-baseline test gives ±0.2% "dips" and "bumps" at random phases (0.18, 0.35, 0.42, 0.75) because of the out-of-eclipse wave. Upper limit ≈ 0.2%. |
| Primary − secondary | ≈ 2.2% (≈ 2.6% at the core). Secondary/primary < 0.08. |
| Out-of-eclipse wave | 1.4% peak to peak; maximum at φ 0.16, minimum at φ 0.49. Harmonic amplitudes (ppt): a₁ 5.88, b₁ −0.06, a₂ −3.19, b₂ 0.63. |
| Odd / even primary depth | 2.30% / 2.34% (consistent) |

Interpretation: either the companion is much cooler or fainter than the primary, or P = 5.579 d with two near-equal eclipses. Equal odd/even depths allow both. No Teff, distance or radii are available in the cloud copy of the tables; they are in the mac-run record.
