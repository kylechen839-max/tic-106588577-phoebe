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

## Interpretation

The multi-sector extraction is clearly useful. It shows that the old period was the main reason sectors did not align over the multi-year baseline.

The remaining residuals are coherent in phase, not random TESS scatter. The strongest residual periodogram peak is near `3.09` cycles per orbital phase, and the residual plot shows the current model underfits the primary eclipse depth and has shape mismatch near the secondary dip.

Recommended modeling path:

1. Use the refined multi-sector binned light curve as the new input.
2. Keep the refined period fixed initially.
3. Do not fully adopt lc_geometry's `requivsumfrac` from the multi-sector curve.
4. Fit a controlled parameter set around the old Powell geometry:
   - `incl@binary`
   - `t0_supconj@binary`
   - `requiv@primary` or `requivsumfrac@binary`
   - `teffratio@binary`
   - possibly `gravb_bol@primary` and `gravb_bol@secondary`
5. Keep the secondary slightly detached during optimization, for example fill factor `0.995` to `0.999`, to avoid PHOEBE projection failures at exact Roche contact.

## Figures

- `results/test/tess_multisector/period_refinement_search.png`
- `results/test/tess_multisector/old_vs_refined_period_fold.png`
- `results/test/tess_multisector/multisector_binned_fold.png`
- `results/test/tess_multisector/multisector_bin_counts.png`
- `results/test/tess_multisector/old_geometry_multisector_detached_diagnostic.png`
- `results/test/tess_multisector/multisector_residual_periodogram.png`
