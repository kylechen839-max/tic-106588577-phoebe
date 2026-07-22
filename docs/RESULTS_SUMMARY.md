# Results Summary

## Final Results

- `results/final/tic_106588577_final.phoebe`
- `results/final/tic_106588577_final_vm.phoebe`
- `results/figures/powell_diagnostic_lightcurve.png`
- `results/figures/emcee_diagnostic_lightcurve.png`

## Test Results

- `results/test/mode_b_lightcurve.png`
- `results/test/mode_b_residuals.png`
- `results/test/emcee_mode_comparison.png`
- `results/test/emcee_mode_residuals.png`
- `results/test/emcee_posterior_models_autofig.json`

## Current Interpretation

The Powell diagnostic remains the best current visual result. Its VM-generated diagnostic values were:

- Chi-squared: `127974.04617631363`
- RMS residual: `0.0014644575928165597`

The emcee/ellc posterior diagnostic produced a valid figure, but those posterior curves fit the observed light curve poorly compared with the Powell model. Keep those outputs as test results until the ellc/emcee model setup is reconciled with the PHOEBE native model.
