# J1122 (J112238.89-592027.5, TIC 451126962): cleaned 7-sector QLP light curve (Kyle's PC, 2026-11-06)

Downloaded from MAST (QLP HLSP) on Kyle's PC, since the cloud cannot reach MAST. Sectors 10, 11
(30 min), 37, 99 (10 min), 64, 90, 100 (200 s). Script `j1122_qlp.py`; console output `run_log.txt`.

## Cleaning
- QLP `sap_flux`, cadences with QLP `quality == 0` only. This drops 6-28% per sector: QLP's
  scattered-light bits (4096, 2^29, 2^30) are not in lightkurve's default mask, and the S100 bumps at
  BTJD 4085-4087 and 4099-4100 are among the dropped cadences.
- Not QLP's detrended flux. KSPSAP (S10-S37) and DET_FLUX (S64-S100) were tried first. Their splines
  partly fit the ~8-h primary and leave +0.5% shoulders beside it. They also cut S37's depth to 1.8%.
- Detrending: each segment (split at gaps > 0.5 d) is divided by a 1.0-d running median of the
  out-of-eclipse flux. Eclipses are masked at |phase| < 0.04 and |phase - 0.5| < 0.04, with the
  trend interpolated across them. No shoulders remain.
- `flux_err`: QLP has no SAP error column, so it is each segment's out-of-eclipse point-to-point
  scatter.
- Clipping: 5σ (MAD) on residuals from a 500-bin phase median, per sector. Removed 25 points.
- Ephemeris P 4.7917045 d, T0 BTJD 3762.8427 (docs/J1122_EPHEMERIS_CORRECTION.md).

## Files
- `_qlp_clean.csv`: 41,325 points. Columns: time_btjd, flux, flux_err, sector, orbit (segment),
  exptime_s, phase.
- `_qlp_all_with_flags.csv`: same, plus the 25 clipped points (`clipped`).
- `_qlp_folded_binned400.csv`: 400-bin phase median with error.
- `_qlp_sector_log.csv`, `_qlp_sectors.png`, `_qlp_folded.png`.

## Numbers
| sector | cadence | points | rms (ppt) | primary depth (\|φ\| < 0.005) | dip at φ 0.5 (\|Δφ\| < 0.005) |
|---|---|---|---|---|---|
| 10 | 1800 s | 915 | 7.7 | 4.24% | 0.40% |
| 11 | 1800 s | 955 | 7.2 | 3.97% | 0.54% |
| 37 | 600 s | 2,835 | 12.8 | 3.85% | 0.70% |
| 64 | 200 s | 10,673 | 20.3 | 4.59% | 0.03% |
| 90 | 200 s | 10,855 | 19.1 | 5.42% | 0.31% |
| 99 | 600 s | 7,088 | 4.9 | 4.61% | 0.47% |
| 100 | 200 s | 8,004 | 4.6 | 4.27% | 0.51% |

The table's depths are relative to a local 0.04-0.08 phase baseline. The stack below uses a
0.04-0.12 baseline, which is why the secondary values differ slightly.

- **Primary (7-sector stack):** 4.68 ± 0.06% at the core. The per-sector spread of 3.9-5.4% is
  larger than the errors, which points to different dilution per sector. That fits the per-sector
  fitter's free depth scales.
- **Secondary at φ 0.5:** 0.32 ± 0.04% over |Δφ| < 0.01, and 0.19 ± 0.03% over |Δφ| < 0.02. The
  formal significance is 8σ. But the same test at 24 other out-of-eclipse phases gives an rms of 0.11%
  (largest 0.29%), so the empirical significance is about 3σ. The two quietest sectors both show it:
  S99 0.47 ± 0.04% and S100 0.33 ± 0.04% (|Δφ| < 0.02). This is a tentative secondary of ~0.3%.
- S64 and S90 are the noisiest (~20 ppt per 200-s point) and add little.
