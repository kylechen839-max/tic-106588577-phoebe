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
- VM work has been abandoned for the active workflow. Current runs are local on macOS using `.venv` and `phoebe==2.4.22`.
- The earlier emcee/ellc posterior diagnostic does not match the observed light curve as well as the native PHOEBE result, so it should be treated as a test result.
- Local tests use the saved `after_powell` model in `bundles/tic_106588577_powell_ready.phoebe` as the baseline comparison.
- The original TESS extraction notebook used only Sector 07 even though six TESScut sectors were available. See `docs/TESS_PIPELINE_REVIEW.md` for the multi-sector extraction path.
- The six-sector TESS test refined the period to `1.0118536926383312 d` and produced a binned PHOEBE input. See `docs/MULTISECTOR_TESS_RESULTS.md`.

## DD × Eclipsing-Binary Candidate Search (2026-10-02)

- `src/crossmatch/` cross-matches Disk Detective DD1.0 (412,729 WZ subjects) and DD2.1 (1,419) against 16 eclipsing-binary surveys: Gaia DR3, TESS EBs, TESS OBA EBs, ZTF, ATLAS, ASAS-SN, ASAS-3, VSX, WISE, CRTS N/S, NSVS, Kepler, OGLE-IV, OGLE LMC and DEBCat. It finds 800 DD1 and 9 DD2 matches, then applies blackbody/SED, WISE-quality, confusion, TESS-usability and object-type sanity checks. Result: **15 tier-A** and 116 tier-B DD1 candidates; DD2 has none (W4 upper limits only). See `docs/DD_EB_CROSSMATCH.md`.
- `src/pipeline/fetch_tess_lc.py` builds multi-sector, eclipse-adaptive, sector-consensus binned TESS inputs. It runs an odd/even (P vs 2P) test and detects eccentric orbits.
- `src/pipeline/run_phoebe_candidate.py` is the generic PHOEBE pipeline: grid → Nelder–Mead → emcee → high-res recompute, with optional third-light and eccentricity fitting. It then runs post-fit photometry checks and a binary-blackbody SED check (photometric vs Gaia distance, and whether the W3/W4 excess survives the binary photosphere).
- `src/pipeline/eclipse_timing.py` produces per-sector O−C eclipse timings.
- **J1122 period correction**: the period used in all earlier J1122 work (4.7909785 d) is wrong. Seven-sector eclipse timing gives **P = 4.7917045 ± 0.0000058 d**. The old ephemeris shifted S10/S11/S37 eclipses by −1.2/−0.9/+1.3 h, which explains the sector-dependent eclipse bias. TIC also lists 37% TESS contamination for J1122. See `docs/J1122_EPHEMERIS_CORRECTION.md`.
- J0719 (TIC 106588577) ephemeris check: P = 1.0118518 ± 0.0000003 d, consistent with the README period to 2×10⁻⁶ d (minor).

## Recreate Diagnostic Figures

From an environment with PHOEBE installed:

```bash
cd tic-106588577-phoebe
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip setuptools wheel
.venv/bin/python -m pip install -r requirements.txt
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
