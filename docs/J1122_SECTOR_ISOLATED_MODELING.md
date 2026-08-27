# J1122 Sector-Isolated Eclipse Modeling

Date: 2026-08-26
Target: J112238.89-592027.5

## Goal

Reduce the eclipse residual bias in the J1122 PHOEBE model without using an empirical residual correction. This pass tested whether the mismatch was caused by PHOEBE parameters, sector-to-sector TESS differences, or sparse folded-bin construction.

## Tests Run

### Per-sector PHOEBE datasets

Script: `src/run_j1122_sector_dataset_model.py`

This rebuilt the PHOEBE bundle with separate datasets for sectors 10, 11, and 37, each with independent dataset scaling.

Result: worse than the merged/consensus input.

| Model | chi2 | RMS | primary mean abs | trend |
|---|---:|---:|---:|---:|
| per-sector start | 2741.21 | 0.001270 | 0.002059 | 0.011822 |

Interpretation: independent PHOEBE dataset scaling does not remove the eclipse bias. The per-sector plot shows the primary eclipse points are sparse and sector-dependent.

### Robust sector-consensus folded inputs

Script: `src/run_j1122_robust_binning_phoebe.py`

This rebuilt the phase-binned light curve using equal-weight sector medians, so one sector cannot dominate a sparse phase bin.

Best result: require each phase bin to have support from all three sectors.

| Input variant | chi2 | RMS | primary mean abs | trend |
|---|---:|---:|---:|---:|
| previous targeted merged input | 639.71 | 0.000880 | 0.001243 | 0.008486 |
| sector consensus, at least 1 sector | 487.83 | 0.000894 | 0.001283 | 0.008551 |
| sector consensus, at least 2 sectors | 375.75 | 0.000858 | 0.001286 | 0.008744 |
| sector consensus, all 3 sectors | 317.45 | 0.000671 | 0.000769 | 0.006631 |
| all 3 sectors with 3-sigma continuum clipping | 292.95 | 0.000663 | 0.000778 | 0.006730 |

Interpretation: the largest improvement came from the light-curve construction, not from changing binary geometry. Requiring all three sectors in a bin reduced the primary-eclipse bias by about 38 percent relative to the previous best targeted PHOEBE model.

### Local refinement on all-three-sector consensus

Script: `src/run_j1122_min3_refinement.py`

This explored small changes in `t0_supconj`, inclination, secondary radius, and temperature ratio using the all-three-sector consensus light curve.

| Selected result | chi2 | RMS | primary mean abs | trend |
|---|---:|---:|---:|---:|
| all-three-sector start | 317.45 | 0.000671 | 0.000769 | 0.006631 |
| best final chi2 | 311.99 | 0.000670 | 0.000775 | 0.006709 |
| best final primary mean abs | 319.02 | 0.000673 | 0.000769 | 0.006814 |

Interpretation: parameter refinement gives only marginal chi-square improvement and does not materially reduce the primary-eclipse bias beyond the robust all-three-sector input.

## Current Best Practical Result

The best eclipse-bias reduction is:

- input: `outputs/j1122_robust_binning_phoebe/sector_consensus_min3_lightcurve.txt`
- model candidate: original bounded-Powell physical solution
- plot: `outputs/j1122_robust_binning_phoebe/sector_consensus_min3_start_diagnostic.png`
- bundle: `outputs/j1122_robust_binning_phoebe/sector_consensus_min3_best.phoebe`

This is not a new exotic binary solution. It is the same compact detached model evaluated against a more reliable folded light curve that requires agreement among sectors.

## Main Conclusion

The eclipse mismatch is now reduced as much as this pass could achieve without forcing PHOEBE into an artificial parameter solution. The strongest evidence is that:

- dropping any one sector entirely makes the fit worse;
- fitting separate sector datasets makes the fit worse;
- requiring all sectors to contribute to each folded bin sharply reduces the primary-eclipse bias;
- additional local PHOEBE parameter changes barely improve the remaining bias.

The remaining structure is likely limited by sparse/sector-dependent TESS sampling near eclipse, not by an obvious missing detached-binary parameter in the current PHOEBE model family.
