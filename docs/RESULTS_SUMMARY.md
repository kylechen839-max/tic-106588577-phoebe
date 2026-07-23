# Results Summary

## Final Results

- `results/final/tic_106588577_final.phoebe`
- `results/final/tic_106588577_final_vm.phoebe`
- `results/final/tic_106588577_emcee_best.phoebe`
- `results/figures/powell_diagnostic_lightcurve.png`
- `results/figures/emcee_diagnostic_lightcurve.png`
- `results/figures/emcee_best_diagnostic_lightcurve.png`

## Test Results

- `results/test/mode_b_lightcurve.png`
- `results/test/mode_b_residuals.png`
- `results/test/emcee_mode_comparison.png`
- `results/test/emcee_mode_residuals.png`
- `results/test/emcee_posterior_models_autofig.json`
- `results/test/native_emcee_test1/emcee_experiment_failure.json`
- `results/test/native_emcee_probe1/emcee_experiment_result.json`

## Current Interpretation

The current best emcee-derived result is:

- Model bundle: `results/final/tic_106588577_emcee_best.phoebe`
- Diagnostic figure: `results/figures/emcee_best_diagnostic_lightcurve.png`
- Chi-squared: `127319.63687960047`
- Target `after_powell` chi-squared: `127974.04617631363`
- Improvement: `654.4092967131585`

The Powell diagnostic values were:

- Chi-squared: `127974.04617631363`
- RMS residual: `0.0014644575928165597`

The earlier emcee/ellc posterior diagnostic produced a valid figure, but those posterior curves fit the observed light curve poorly compared with the native PHOEBE model. Keep those outputs as test results until the ellc/emcee model setup is reconciled with the PHOEBE native model.
