# J1122 Powell and emcee Continuation

Date: 2026-08-11
Target: J112238.89-592027.5
Input light curve: `outputs/j1122_tess_multisector_s10_s11_s37/J112238.89-592027.5_LightCurve_multisector_binned_220.txt`

## Purpose

This run continued from the compact detached J1122 baseline and tested a more local optimization path:

- keep the detached PHOEBE model;
- fix period, semi-major axis, eccentricity, and gravity brightening for this pass;
- optimize inclination, `t0_supconj`, primary radius, secondary radius, and `teffratio`;
- use bounded external Powell calls where every likelihood evaluation is a PHOEBE model;
- run a short local emcee smoke test from the Powell endpoint.

The driver is `src/run_j1122_powell_emcee_continuation.py`.

## High-Resolution Comparison

All values below are recomputed with `ntriangles=1500`.

| Model | incl | t0 | r1 | r2 | teffratio | chi2 | RMS | mean abs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Detached start | 87.0000 | 1573.086900 | 0.7000 | 0.4550 | 0.7500 | 679.37 | 0.000902 | 0.000524 |
| Bounded Powell | 87.0064 | 1573.087569 | 0.6989 | 0.4551 | 0.6980 | 660.60 | 0.000894 | 0.000514 |
| Short emcee best | 86.9571 | 1573.088059 | 0.6970 | 0.4696 | 0.7022 | 670.65 | 0.000897 | 0.000520 |

The bounded Powell solution is currently the best physical PHOEBE solution from this batch.

## Interpretation

Powell made a modest but real improvement over the detached baseline. Most of the improvement came from shifting `t0_supconj` by about 0.00067 d and lowering `teffratio` from 0.75 to about 0.70. The radii and inclination stayed close to the compact detached baseline.

The short emcee run moved locally and produced a reasonable nearby solution, but it did not beat the high-resolution Powell endpoint. The acceptance fraction was about 0.21, which is usable for a smoke test, but 40 iterations is far too short to establish posterior convergence. The corner plot also warned that there were too few points for reliable contours.

## Remaining Mismatch

The residual plot is better than the earlier broad/preliminary models, especially out of eclipse, but the primary eclipse still has a sharp residual feature near the core. This means PHOEBE is now accounting for more of the light curve through physical parameters, but the eclipse morphology mismatch is not fully solved.

Likely next physical tests:

- allow eccentricity/argument of periastron if timing/asymmetry persists;
- test sector-by-sector normalization or per-sector flux offsets before fitting deeper physics;
- fit limb darkening or atmosphere choices if the residual is concentrated in ingress/egress/core shape;
- run a longer emcee only after the deterministic model family is improved.

## Generated Products

- `outputs/j1122_powell_emcee_continuation/j1122_powell_emcee_result.json`
- `outputs/j1122_powell_emcee_continuation/j1122_powell_history.json`
- `outputs/j1122_powell_emcee_continuation/j1122_start_hires_result.json`
- `outputs/j1122_powell_emcee_continuation/j1122_after_bounded_powell_hires_result.json`
- `outputs/j1122_powell_emcee_continuation/j1122_after_bounded_powell_hires_diagnostic.png`
- `outputs/j1122_powell_emcee_continuation/j1122_emcee_best_diagnostic.png`
- `outputs/j1122_powell_emcee_continuation/j1122_emcee_trace.png`
- `outputs/j1122_powell_emcee_continuation/j1122_emcee_corner.png`
- `outputs/j1122_powell_emcee_continuation/j1122_after_bounded_powell_hires.phoebe`
- `outputs/j1122_powell_emcee_continuation/j1122_emcee_best.phoebe`
