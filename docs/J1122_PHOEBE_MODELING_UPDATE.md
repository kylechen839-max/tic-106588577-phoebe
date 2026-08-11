# J1122 PHOEBE Modeling Update

Date: 2026-08-11
Target: J112238.89-592027.5
Input light curve: `outputs/j1122_tess_multisector_s10_s11_s37/J112238.89-592027.5_LightCurve_multisector_binned_220.txt`
Period fixed: 4.790978516308605 d
Initial t0 fixed: 1573.0869 BTJD

## What changed

The earlier J1122 PHOEBE run relied on `lc_geometry` and a very small emcee probe. That produced a poor preliminary model with RMS about 0.062 on the 220-bin light curve. I replaced that with a controlled detached-search workflow in `src/run_j1122_detached_search.py`.

The new workflow:

- builds a fresh `phoebe.default_binary()` for J1122 instead of reusing the TIC 106588577 bundle;
- uses the multisector 220-bin TESS light curve;
- fixes the period and circular orbit for this first stable pass;
- uses blackbody atmospheres and manual limb darkening for local robustness;
- applies a conservative sigma floor of 0.0004 so sparse bins do not dominate chi-square;
- searches detached spherical models over inclination and radii;
- saves both lowest-RMS and trend-weighted diagnostic models.

## Best current result

Best high-resolution RMS recompute:

- model file: `outputs/j1122_detached_search/j1122_detached_best_hires.phoebe`
- diagnostic plot: `outputs/j1122_detached_search/best_rms_hires_diagnostic.png`
- inclination: 87.0 deg
- sma: 20.0
- requiv primary: 0.7
- requiv secondary: 0.455
- teff ratio: 0.75
- chi2: 679.37
- RMS: 0.000902
- mean absolute residual: 0.000524

This is much better than the previous preliminary J1122 result and is now in the right shallow-eclipse regime.

## Caveats

The model is improved, but not final. The primary-eclipse residual still has a systematic negative median around the eclipse shoulders/core, so PHOEBE is not yet fully accounting for the primary eclipse morphology. A wider-radius focused test was run, but it did not beat the compact high-resolution RMS model.

The current result should be treated as a stable detached baseline for J1122, not a final physical solution. The next real improvement should fit `t0`, radii, inclination, and temperature ratio with a bounded optimizer/emcee initialized around this compact solution, ideally with a continuous model grid and possibly sector-by-sector normalization/systematics terms.

## Generated products

- `outputs/j1122_detached_search/j1122_detached_search_result.json`
- `outputs/j1122_detached_search/j1122_detached_search_all_results.json`
- `outputs/j1122_detached_search/j1122_detached_best.phoebe`
- `outputs/j1122_detached_search/j1122_detached_best_hires.phoebe`
- `outputs/j1122_detached_search/best_rms_hires_diagnostic.png`
- `outputs/j1122_detached_search/best_objective_hires_diagnostic.png`
- `outputs/j1122_detached_search/j1122_detached_wide_best_result.json`
- `outputs/j1122_detached_search/j1122_detached_wide_best_diagnostic.png`
- `outputs/j1122_detached_search/j1122_detached_wide_best.phoebe`
