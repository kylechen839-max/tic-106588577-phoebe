# TESS Pipeline Review

Reviewed local notebook:

- `/Users/kylechen/Downloads/J071951_40_240400_6_LightkurveBasicPipelineForTESS(Kyle).ipynb`

## What It Currently Does

- Searches TESScut with `lk.search_tesscut("07 19 51.40 -24 04 00.6")`.
- Builds a threshold aperture for the target.
- Builds a separate background aperture.
- Uses `RegressionCorrector` with PCA background regressors.
- Exports `J071951.40-240400.6_LightCurve.txt` with `Time`, `Flux`, and `Flux_Err`.

## Important Finding

The saved notebook output shows six available TESScut sectors:

- Sector 07, 2019
- Sector 33, 2020
- Sector 34, 2021
- Sector 61, 2023
- Sector 87, 2024
- Sector 88, 2025

The notebook then runs:

```python
chosen_pf = pfs[0]
pf = chosen_pf.download(cutout_size=20)
```

So the current PHOEBE input uses only the first search result, Sector 07. The three downloaded copies of `J071951.40-240400.6_LightCurve.txt` are identical one-sector exports:

- 1086 data rows
- time range `1491.661` to `1516.077` BTJD
- span about `24.416` days
- median cadence about `0.021` days

## Why This Can Help PHOEBE

Combining sectors can help distinguish real binary-shape residuals from TESS sector-level systematics. If the same residual pattern appears in all sectors after separate sector normalization, it is more likely astrophysical or model-related. If the pattern changes by sector, it is more likely extraction/systematics-related.

## Added Tool

Added:

- `src/extract_multisector_tess_lightcurve.py`
- `requirements-tess.txt`

The script keeps the notebook's Lightkurve approach but processes all sectors separately, normalizes each sector by its own median, sigma-clips outliers, and writes a PHOEBE-compatible combined file whose first three columns remain:

```text
Time    Flux    Flux_Err
```

It also includes a fourth `Sector` column for diagnostics. The current PHOEBE loader reads the first three columns, so this remains compatible.

## Suggested Test

Create a separate environment for Lightkurve so the active PHOEBE environment is not disturbed:

```bash
cd "/Users/kylechen/Documents/Emory/Research/EB x DD/Scripts/tic-106588577-phoebe"
python3 -m venv .venv-tess
.venv-tess/bin/python -m pip install --upgrade pip setuptools wheel
.venv-tess/bin/python -m pip install -r requirements-tess.txt
.venv-tess/bin/python src/extract_multisector_tess_lightcurve.py
```

Then test the combined extraction in PHOEBE:

```bash
LIGHTCURVE_PATH=outputs/tess_multisector/J071951.40-240400.6_LightCurve_multisector.txt \
.venv/bin/python src/tic_106588577_phoebe.py
```

Start by comparing the Powell-stage residuals from the one-sector and multi-sector inputs before rerunning long emcee. If the multi-sector residuals are cleaner or less coherent, it is worth rebuilding the Powell-ready bundle from that input.
