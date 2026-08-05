# J1122 Preliminary Multi-Sector Results

Target: `J112238.89-592027.5`

## TESS Extraction

The local Downloads notebook gave:

- target name: `J112238.89-592027.5`
- search string: `11 22 38.89 -59 20 27.5`
- custom aperture: `y=4:6`, `x=4:6`
- background aperture: `y=8:10`, `x=5:10`

MAST returned seven TESScut sectors: 10, 11, 37, 64, 90, 99, and 100. Sectors 10, 11, and 37 were extracted locally and combined. Sector 64 download stalled during this run, so the current PHOEBE input is a three-sector product.

Extracted points:

- Sector 10: 1108
- Sector 11: 1139
- Sector 37: 3378
- Total: 5625

Period search on the three-sector product:

- Lomb-Scargle period: `2.3875184586409772 d`
- BLS period: `4.790978516308605 d`
- BLS t0: `1573.0869`

The BLS period was adopted for the preliminary eclipsing-binary fold.

## PHOEBE/emcee Probe

A J1122-specific runner was added at `src/run_j1122_phoebe_emcee.py`.

The full Powell stage was attempted but stopped because native PHOEBE evaluation was too slow locally: one Powell iteration took about 2.7 minutes. The successful run therefore used:

- 220 phase bins
- detached binary
- lc_geometry-derived baseline
- no Powell optimization
- emcee probe with 8 walkers and 10 iterations
- sampled parameters: `incl@binary`, `t0_supconj@binary`

This produced a finite emcee model, but it is not converged and not a good physical fit yet.

Results:

- baseline model: `j1122_after_lcgeo`
- baseline chi-squared: `6613385.147542531`
- baseline RMS: `0.06250797676066748`
- emcee-best chi-squared: `6083007.995540482`
- emcee-best RMS: `0.061844010138341086`

Interpretation: the local PHOEBE/emcee pipeline is operational for J1122, but the binary model is preliminary and poor. A publishable binary solution needs either a longer Powell run, a better initial geometry, more aggressive binning for optimization, or a simpler/non-PHOEBE exploratory model before returning to PHOEBE.

## WISE Excess Diagnostic

The exact local AllWISE row was found in:

`/Users/kylechen/Downloads/WZ_subjects_extflag+1.5to+2.txt`

WISE/2MASS row:

`J112238.89-592027.5 170.6620748 -59.3409997 ...`

The free-free and dust blackbody diagnostic was run with:

```bash
python3 src/wise_excess_diagnostic.py \
  --designation J112238.89-592027.5 \
  --out-dir results/test/wise_excess_j1122
```

Key result:

- W3 flux/photosphere: `9.78`
- W4 flux/photosphere: `40.18`
- W3 excess significance: `31.44 sigma`
- W4 excess significance: `12.60 sigma`
- optically thin free-free chi-squared: `225.50`
- partially thick wind free-free chi-squared: `505.00`
- best plausible free-free chi-squared: `225.50`
- single-temperature dust chi-squared: `16.44`
- dust temperature: roughly `290-356 K`, depending on whether fitting only W3/W4 residuals or all four WISE bands

Interpretation: J1122 has a real mid-infrared excess in W3/W4. Ordinary free-free emission does not fit the WISE SED well. A dust blackbody is more plausible than free-free, but the case is less clean than J0719 because W2 is already above the simple W1/W2 photospheric baseline and W1 is below it. Before claiming a debris disk, inspect WISE quality flags/images and fit a fuller Gaia+2MASS+WISE photosphere.
