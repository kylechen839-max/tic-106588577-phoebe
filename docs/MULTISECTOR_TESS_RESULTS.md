# Multi-Sector TESS Results

## Extraction Result

The original Downloads notebook found six TESScut sectors for `J071951.40-240400.6`, but the PHOEBE input used only Sector 07. I extracted all six sectors with the same Lightkurve aperture/regression-correction approach:

- Sector 07: 1086 points
- Sector 33: 3485 points
- Sector 34: 3474 points
- Sector 61: 10654 points
- Sector 87: 11168 points
- Sector 88: 11625 points

Combined raw extraction:

- `41492` points
- time baseline: `1491.66061` to `3717.92634` BTJD

## Ephemeris Finding

Folding all sectors with the old one-sector binary period did not align the eclipses:

- old period: `1.011031055532147 d`

The multi-sector period search preferred:

- Lomb-Scargle best period: `1.011848612783501 d`
- BLS best period: `1.0118536926383312 d`
- BLS t0: `1492.5441099999998`

I used the BLS period for the combined folded product:

- adopted period: `1.0118536926383312 d`
- adopted t0: `1492.5441099999998`

This aligns all six sectors well in phase.

## Binned PHOEBE Input

To avoid making PHOEBE compute all `41492` points, I folded and median-binned the multi-sector light curve into `1200` phase bins:

- `data/J071951.40-240400.6_LightCurve_multisector_binned.txt`
- median raw points per bin: `35`
- minimum raw points per bin: `25`
- maximum raw points per bin: `42`

The first three columns are PHOEBE-compatible:

```text
Time    Flux    Flux_Err
```

Extra columns record phase and bin count for diagnostics.

## PHOEBE Test

A full lc_geometry -> Powell rerun on the multi-sector input was not stable. The lc_geometry solution did not identify the secondary eclipse and pushed:

- `requivsumfrac = 0.7522599371673773`

This drove the semidetached model into repeated Roche-boundary and projection failures during Powell. I stopped that run.

As a safer diagnostic, I loaded the existing one-sector Powell-ready bundle, updated:

- dataset to the multi-sector binned light curve,
- period to `1.0118536926383312 d`,
- t0 to `1492.5441099999998`,
- secondary fill factor to `0.995` of `requiv_max` for numerical stability.

Result:

- chi-squared: `37097.54537827828`
- RMS residual: `0.0011397137096130843`
- mean absolute residual: `0.0007543885310056284`
- max absolute residual: `0.005439725591059208`

This RMS is lower than the previous one-sector `after_powell` RMS of `0.0014644575928165597`, but the chi-squared is not directly comparable because the multi-sector binned errors are constructed from bin scatter.

## Controlled Local Search

I ran a controlled local search around the previous Powell geometry instead of adopting the unstable multi-sector lc_geometry radius estimate. The search used a 240-bin version of the multi-sector light curve for speed, then evaluated the accepted candidate on the full 1200-bin input.

Best accepted candidate:

- period: `1.0118536926383312 d`
- `t0_supconj@binary`: `1492.54561 d`
- `incl@binary`: `159.59668610123063 deg`
- `requiv@primary`: `1.5360075942651437 solRad`
- `teffratio`: `0.7208571428571429`
- `gravb_bol@primary`: `0.9`
- `gravb_bol@secondary`: `1.0`
- secondary fill factor: `0.995`

Full 1200-bin result:

- chi-squared: `35867.291156438136`
- RMS residual: `0.001118510879687342`
- mean absolute residual: `0.0007281583647954803`
- max absolute residual: `0.0054177624871794006`

Compared with the old Powell geometry evaluated on the same multi-sector binned light curve:

- chi-squared improved from `37097.54537827828` to `35867.291156438136`
- RMS improved from `0.0011397137096130843` to `0.001118510879687342`
- mean absolute residual improved from `0.0007543885310056284` to `0.0007281583647954803`

The improvement is modest. The diagnostic figure still shows coherent residual structure, with the primary eclipse core underfit and a smaller mismatch near the secondary dip. This means the multi-sector data and refined period help, but the current detached/near-contact PHOEBE geometry still does not fully describe the eclipse shapes.

## Interpretation

The multi-sector extraction is clearly useful. It shows that the old period was the main reason sectors did not align over the multi-year baseline.

The remaining residuals are coherent in phase, not random TESS scatter. The strongest residual periodogram peak is near `3.09` cycles per orbital phase, and the residual plot shows the current model underfits the primary eclipse depth and has shape mismatch near the secondary dip.

## Figures

- `results/test/tess_multisector/period_refinement_search.png`
- `results/test/tess_multisector/old_vs_refined_period_fold.png`
- `results/test/tess_multisector/multisector_binned_fold.png`
- `results/test/tess_multisector/multisector_bin_counts.png`
- `results/test/tess_multisector/old_geometry_multisector_detached_diagnostic.png`
- `results/test/tess_multisector/multisector_residual_periodogram.png`
- `results/test/multisector_controlled_search/controlled_search_best_diagnostic.png`
- `results/test/multisector_controlled_search/controlled_search_result.json`
