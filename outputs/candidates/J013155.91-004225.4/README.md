# J013155.91-004225.4: cleaned multi-sector TESS light curve (run 2026-10-26 on Kyle's PC)

Target: TIC 248986951, DD2 object, G 15.3, Tmag 13.8, parallax 31.5 mas (nearby M dwarf).
VSX lists it as SDSS J013155.94-004224.7, type EA, P = 0.351176 d, max r = 16.34. That
VSX entry has no epoch and no amplitude. In the DD crossmatch it is tier C (W4 upper
limit, VSX only).

## What was made

`J013155.91-004225.4_tess_clean.csv`: 65,055 points from 5 sectors. Columns: time_btjd,
flux (per-orbit median normalised), flux_err, sector, orbit, author, exptime_s.

| sector | product | cadence | flux | points kept | rms per point | CROWDSAP |
|---|---|---|---|---|---|---|
| 3 | SPOC | 120 s | PDCSAP | 12,948 / 12,951 | 14.2 ppt | 0.81 |
| 30 | SPOC | 120 s | PDCSAP | 16,082 / 16,083 | 13.3 ppt | 0.86 |
| 42 | SPOC | 120 s | PDCSAP | 13,033 / 13,035 | 14.8 ppt | 0.84 |
| 43 | SPOC | 120 s | PDCSAP | 14,037 / 14,039 | 12.9 ppt | 0.87 |
| 70 | TESS-SPOC | 200 s | PDCSAP | 8,955 / 8,961 | 11.4 ppt | 0.85 |
| 97 | QLP | 200 s | SAP | **excluded** | 14.1 ppt | n/a |

Sector 97 (QLP, no KSPSAP column) has large spikes and its periodogram is dominated by
multi-day systematics (7.0 / 3.6 / 2.8 d). It is kept, flagged, in
`_tess_all_with_flags.csv` (`use` = False) but left out of the clean file. Sector 108
exists only as TESScut FFIs and was not extracted.

Cleaning follows docs/PROJECT_CONTEXT.md: SPOC, then TESS-SPOC, then QLP per sector; the
default quality bitmask; NaNs removed; per-orbit median normalisation. Outlier rejection is
a 5σ (MAD) cut on residuals from a phase-binned median, iterated 3 times, done per sector.
It removed 11 points in sectors 3-70. Scripts are in `scripts/` (written fresh for this
run; the repo pipeline was not run because the PC's permission policy blocks running code
from the clone). The full console output is in `run_log.txt`.

## What the light curve shows

- **No eclipses at the VSX period.** Folded at 0.351176 d, every sector is flat. The BLS
  depth at that period is 0.8-1.7 ppt at SNR 2-3.5, i.e. not detected. The 2-min data
  would resolve an M-dwarf eclipse (~1 h).
- **A 0.4878 d quasi-sinusoidal signal in every sector.** Lomb-Scargle periods are
  0.48803, 0.48768, 0.48773, 0.48814 and 0.48751 d (mean 0.48782 ± 0.00026 d). Semi-amplitude
  is 3.5-4.8 ppt. FAP is below 1e-90 in each sector. The waveform changes from sector to
  sector: one hump per cycle, asymmetric, with amplitude and shape varying
  (`_per_sector_fold_LS.png`).
- No box-shaped eclipse at any period. BLS just picks up the sinusoid (0.487 d and its
  half and double). At 2P = 0.974 d the two halves have the same shape, with no odd/even
  difference that would suggest an ellipsoidal binary.

## Reading (interpretation, not proven)

The changing single-humped waveform at 0.488 d looks like starspot rotation of a
fast-rotating, active M dwarf, not an eclipsing or ellipsoidal binary. The VSX EA
classification and 0.351 d period are not confirmed by TESS. CROWDSAP 0.81-0.87 means
13-19% of the aperture flux is from other stars. PDCSAP corrects for that, but it
cannot be ruled out that the modulation belongs to a neighbour. A pixel-level check
(TESScut difference image) would settle it. Either way, this star is not a usable EB for the
PHOEBE pipeline from TESS data.
