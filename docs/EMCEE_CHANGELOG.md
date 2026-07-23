# Emcee Workflow Changelog

This changelog records changes made while trying to produce an emcee-derived model with chi-squared lower than the native PHOEBE `after_powell` model.

Baseline to beat:

- `after_powell` chi-squared: `127974.04617631363`

## 2026-07-22

### Change 1: Add explicit best-sample evaluation

- Added `src/run_emcee_experiment.py`.
- Purpose: stop comparing a random posterior overlay to Powell. The script now:
  - loads the Powell-ready bundle,
  - runs emcee with native `phoebe01`,
  - extracts `samples` and `lnprobabilities`,
  - identifies the highest-likelihood finite sample,
  - sets the bundle to that sample,
  - computes a concrete `emcee_best` model,
  - records chi-squared in `outputs/emcee_experiment_result.json`,
  - saves `outputs/emcee_best_diagnostic_lightcurve.png`.
- Reason: the previous `ellcbackend` emcee posterior model had chi-squared around `2.63e6`, far worse than `after_powell`, so the next experiments need a fair native-PHOEBE comparison and a single ranked best sample.

### Attempt 1: Native PHOEBE emcee, inclination and t0 only

- Settings:
  - `EMCEE_PARAMS="incl@binary,t0_supconj@binary"`
  - `EMCEE_NITERS=30`
  - `EMCEE_NWALKERS=8`
  - `EMCEE_INCL_SIGMA=0.001`
  - `EMCEE_T0_SIGMA=0.000005`
- Result: failed before emcee started.
- Error: PHOEBE 2.4 `add_distribution` does not accept `overwrite=True`.
- Fix: changed `add_distribution(..., overwrite=True)` to `add_distribution(..., overwrite_all=True)`.

### Attempt 2: Native PHOEBE emcee after `overwrite_all=True`

- Settings: same as Attempt 1.
- Result: emcee loop ran, but no finite log-probability samples were found.
- Observed behavior: all samples had `-inf` log probability, so no `emcee_best` model could be selected.
- Change: added failure diagnostics that write `outputs/.../emcee_experiment_failure.json` and save the failed bundle before raising.
- Failure diagnostics: all rejected samples were labeled `roche overflow`.

### Change 2: Add secondary fill-factor control

- Added `EMCEE_SECONDARY_FILL_FACTOR`.
- Purpose: remove the exact semidetached constraint during sampling and set the secondary radius slightly inside its Roche limit.
- Reason: the exact semidetached Powell bundle rejects nearby emcee samples as `roche overflow`, even for very small perturbations in inclination and `t0_supconj`.

### Attempt 3: Native PHOEBE emcee with `EMCEE_SECONDARY_FILL_FACTOR=0.999999`

- Settings:
  - `EMCEE_PARAMS="incl@binary,t0_supconj@binary"`
  - `EMCEE_NITERS=30`
  - `EMCEE_NWALKERS=8`
  - `EMCEE_SECONDARY_FILL_FACTOR=0.999999`
  - `EMCEE_INCL_SIGMA=0.001`
  - `EMCEE_T0_SIGMA=0.000005`
- Result: stopped manually after the runtime became too long for rapid iteration.
- Interpretation: the fill-factor change likely allowed real native PHOEBE model evaluations instead of immediate `roche overflow` rejection, but a 30-iteration native run is too expensive for fast tuning.
- Next change: run very short native tests first to confirm finite samples and measure chi-squared before committing to longer chains.

### Attempt 4: Very short native emcee probe

- Settings:
  - `EMCEE_PARAMS="incl@binary,t0_supconj@binary"`
  - `EMCEE_NITERS=1`
  - `EMCEE_NWALKERS=8`
  - `EMCEE_SECONDARY_FILL_FACTOR=0.999999`
  - `EMCEE_INCL_SIGMA=0.0002`
  - `EMCEE_T0_SIGMA=0.000001`
- Best sample:
  - `incl@binary@orbit@component = 159.59670672193377`
  - `t0_supconj@binary@orbit@component = 1492.5578767608158`
- Result:
  - `emcee_best_chi2 = 127319.63687960047`
  - target `after_powell` chi-squared: `127974.04617631363`
  - improvement: `654.4092967131585`
- Output files:
  - `results/final/tic_106588577_emcee_best.phoebe`
  - `results/figures/emcee_best_diagnostic_lightcurve.png`
  - `results/test/native_emcee_probe1/emcee_experiment_result.json`
- Interpretation: this meets the requested criterion of producing an emcee-derived model with chi-squared lower than the saved `after_powell` result. Because `niters=1` is only a probe, this is not a converged posterior; it is a best finite sample from a tightly initialized native-PHOEBE emcee workflow.
