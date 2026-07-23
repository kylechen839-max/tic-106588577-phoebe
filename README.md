# TIC 106588577 PHOEBE Modeling

This repository contains the PHOEBE modeling workflow, cloud/VM job scripts, input light curve, saved bundles, and diagnostic results for TIC 106588577.

## Repository Layout

- `data/` - input light curve used for the PHOEBE model.
- `notebooks/` - exploratory notebooks and PHOEBE reference notebooks.
- `src/` - Python scripts for running the PHOEBE/Powell/emcee workflow and generating figures.
- `cloud/` - Docker and requirements files used during the Cloud Run experiment.
- `bundles/` - saved PHOEBE bundle after the Powell-ready stage.
- `results/final/` - final saved `.phoebe` bundles.
- `results/figures/` - current diagnostic figures.
- `results/test/` - earlier test figures and intermediate outputs.
- `docs/` - project notes and result summaries.

## Current Best Diagnostic Figure

The current best emcee-derived comparison figure is:

`results/figures/emcee_best_diagnostic_lightcurve.png`

This figure compares the observed light curve with the best finite native-PHOEBE emcee sample found so far and plots residuals.

## Current Notes

- The best emcee-derived model currently beats the saved `after_powell` chi-squared target.
- The VM successfully runs PHOEBE with the conda-forge `ellc` installation.
- The earlier emcee/ellc posterior diagnostic does not match the observed light curve as well as the native PHOEBE result, so it should be treated as a test result.

## Recreate Diagnostic Figures

From an environment with PHOEBE installed:

```bash
cd tic-106588577-phoebe
python src/make_powell_diagnostic_figure.py
python src/make_emcee_diagnostic_figure.py
```

The scripts expect the final bundle to be available at:

`outputs/tic_106588577_final.phoebe`

If running directly from this repository, copy or symlink the final bundle first:

```bash
mkdir -p outputs
cp results/final/tic_106588577_final_vm.phoebe outputs/tic_106588577_final.phoebe
```
