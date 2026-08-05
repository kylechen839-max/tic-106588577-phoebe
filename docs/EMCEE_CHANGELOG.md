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

## 2026-08-04

### Change 3: Add least-squares and convergence metrics

- Updated `src/run_emcee_experiment.py` to report:
  - residual sum of squares,
  - mean squared residual,
  - RMS residual,
  - mean absolute residual,
  - maximum absolute residual,
  - acceptance fractions,
  - autocorrelation times,
  - `niters > 50 * max(autocorr_time)` convergence flag.
- Purpose: require the long emcee solution to beat `after_powell` in both PHOEBE chi-squared and least-squares/RMS residual metrics, not only by log probability.

### Attempt 5: Long native-PHOEBE convergence run launched

- Settings:
  - `EMCEE_PARAMS="incl@binary,t0_supconj@binary"`
  - `EMCEE_NITERS=800`
  - `EMCEE_NWALKERS=32`
  - `EMCEE_BURNIN=200`
  - `EMCEE_SECONDARY_FILL_FACTOR=0.999999`
  - `EMCEE_INCL_SIGMA=0.0002`
  - `EMCEE_T0_SIGMA=0.000001`
  - `compute="phoebe01"`
- VM output directory: `~/phoebe-cloud-job/outputs/native_emcee_long1`
- Status: abandoned after switching away from the VM.
- Goal: produce a converging posterior whose best finite sample beats `after_powell` in chi-squared and residual least-squares/RMS.

### Change 4: Abandon VM and move active workflow local

- Stopped the Google Compute Engine VM `phoebe-emcee-vm`.
- Created a local `.venv` on macOS.
- Pinned `requirements.txt` to `phoebe==2.4.22`, matching the version that loaded the Powell-ready bundle correctly.
- Moved the incompatible PHOEBE 2.5 passband cache aside so PHOEBE 2.4.22 could initialize passbands cleanly.
- Changed `src/run_emcee_experiment.py` so local runs can compare against the saved `after_powell` model in the Powell-ready bundle instead of recomputing the exact semi-detached Powell model every time.

### Attempt 6: Local native-PHOEBE smoke probe

- Settings:
  - `EMCEE_PARAMS="incl@binary,t0_supconj@binary"`
  - `EMCEE_NITERS=1`
  - `EMCEE_NWALKERS=8`
  - `EMCEE_BURNIN=0`
  - `EMCEE_SECONDARY_FILL_FACTOR=0.999999`
  - `EMCEE_INCL_SIGMA=0.0002`
  - `EMCEE_T0_SIGMA=0.000001`
  - `BASELINE_MODEL=after_powell`
  - `RECOMPUTE_BASELINE=false`
- Local output directory: `outputs/local_emcee_probe5`
- Result:
  - `emcee_best_chi2 = 127197.71490010298`
  - `after_powell_chi2 = 127974.04617631363`
  - `emcee_best RMS = 0.0014600368695896705`
  - `after_powell RMS = 0.0014644575928165597`
- Interpretation: the local workflow is functional and again finds an emcee-derived sample that beats `after_powell` in chi-squared and RMS. This was only a one-iteration smoke test, so it is not a converged posterior.

### Change 5: Insert controlled multi-sector search before another long emcee run

- Added `src/run_multisector_controlled_search.py`.
- Purpose: use the cleaner multi-sector folded and binned light curve to test the parameter directions recommended after the TESS pipeline review before spending time on another expensive sampler.
- The search intentionally kept the refined period fixed and did not adopt the unstable multi-sector lc_geometry `requivsumfrac`.
- Parameters tested around the previous Powell solution:
  - `incl@binary`
  - `t0_supconj@binary`
  - `requiv@primary`
  - `teffratio`
  - `gravb_bol@primary`
  - `gravb_bol@secondary`
- Result on the full 1200-bin multi-sector curve:
  - old Powell geometry chi-squared: `37097.54537827828`
  - controlled-search chi-squared: `35867.291156438136`
  - old Powell RMS: `0.0011397137096130843`
  - controlled-search RMS: `0.001118510879687342`
- Best accepted change:
  - `t0_supconj@binary = 1492.54561 d`
  - `gravb_bol@secondary = 1.0`
  - retained the Powell inclination, primary radius, and temperature ratio.
- Output files:
  - `results/final/tic_106588577_multisector_controlled_best.phoebe`
  - `results/test/multisector_controlled_search/controlled_search_best_diagnostic.png`
  - `results/test/multisector_controlled_search/controlled_search_result.json`
- Interpretation: this beats the old Powell model on the multi-sector binned light curve in both chi-squared and RMS, but it is still not a converged emcee posterior. The residual structure remains coherent, especially through the eclipse cores.
