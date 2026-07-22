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

The current usable comparison figure is:

`results/figures/powell_diagnostic_lightcurve.png`

This figure compares the observed light curve with the `after_powell` PHOEBE model and plots residuals.

## Current Notes

- The Powell model remains the best-fitting diagnostic model in the saved results.
- The VM successfully runs PHOEBE with the conda-forge `ellc` installation.
- The emcee/ellc posterior diagnostic currently does not match the observed light curve as well as the Powell model, so it should be treated as a test result rather than the final scientific model.

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
