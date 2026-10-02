# HANDOFF: DD × Eclipsing-Binary Search and PHOEBE Modelling

Handoff written 2026-10-02 (end of session). Upload this file to the Project. It is self-contained: goals, environment, what was done, what is still running, and exact next steps. For background detail, see the other docs in the repo (`docs/PROJECT_CONTEXT.md`, `docs/DD_EB_CROSSMATCH.md`, `docs/J1122_EPHEMERIS_CORRECTION.md`).

## Prompt to start the next session

```
Read HANDOFF.md. Continue from "Next steps", starting with step 1 (collect the round-2 PHOEBE
results). The repo is at ~/Documents/repos/tic-106588577-phoebe (GitHub:
kylechen839-max/tic-106588577-phoebe, branch main). Commit and push to main when a step is done.
```

## Goals (from the user)

1. Find eclipsing-binary surveys and cross-match them with Disk Detective DD1.0 and DD2.x objects to find the next targets after J0719 (TIC 106588577) and J1122.
2. Run sanity checks on each cross-match to judge how strong each candidacy is.
3. Build PHOEBE models, with the pipeline developed so far, for the candidates that pass the blackbody and photometry checks.

## Environment

| Item | Value |
|---|---|
| Repo | `~/Documents/repos/tic-106588577-phoebe` (on branch `main`; all session work pushed) |
| PHOEBE Python | `~/phoebe-env/bin/python` (3.11, phoebe 2.4.22, emcee, corner, scipy). No pandas/astroquery. **Never pip-install into it.** |
| Analysis Python | `.venv-tess/bin/python` (astroquery, lightkurve, pandas, tabulate). Gitignored; recreate with `~/phoebe-env/bin/python -m venv .venv-tess && .venv-tess/bin/pip install astroquery lightkurve pandas matplotlib scipy tabulate` |
| DD1 inputs | `~/Downloads/WZSubs.zip` → unzip into `data/catalogs/WZSubs/` (gitignored, 86 MB) |
| DD2 inputs | `data/catalogs/DDv2.1_Wise_Ids.txt` (committed) |
| Machine | macOS, 12 cores; one PHOEBE evaluation takes ~1–2.5 s; up to ~11 fits in parallel |
| `outputs/` | gitignored; commit selected files with `git add -f` |

## Results so far

### 1. Cross-match (`docs/DD_EB_CROSSMATCH.md`)

- 16 EB surveys:
  - Gaia DR3 ECL, VSX, ZTF, ATLAS
  - ASAS-SN, ASAS-3, OGLE LMC, TESS OBA EBs
  - TESS EBs (Prša), WISE variables, CRTS N/S, NSVS
  - OGLE-IV (Wang+2024), Kepler, DEBCat
- **800 DD1** and **9 DD2** objects match at least one survey. DD2 has no usable candidates: all 9 are M dwarfs with W4 upper limits only.
- Sanity checks:
  - match quality, number of independent surveys, period agreement and plausibility;
  - AllWISE quality flags;
  - **blackbody check**: W3/W4 excess ≥ 5σ against both a free blackbody and one fixed at the TIC Teff;
  - **photometry check**: catalogue-Teff photosphere fits J, H, Ks, W1;
  - dust vs free-free, Gaia crowding, TESS contamination, SIMBAD type, distance < 3 kpc;
  - WISE image centring.
- **DD1: 12 tier A, 119 tier B, 669 tier C.** The tables are in `docs/DD_EB_CROSSMATCH.md` and `outputs/crossmatch/tierA_table.md`.
- **Tier A:** J073937.23−363011.8, **J221843.61+544715.0** (best debris-like SED), J024542.12+480837.8, J173224.94−364153.9, J223949.47+583254.4, J231201.40+532028.6, J223553.38+500414.8, J160415.99−562627.3, J141909.42−565518.1, J045912.77+165543.8, J220912.52+582025.0 (RUWE 43, caution), J033226.58+520357.1.
- **Rejected by eye as nebular W3/W4:** J1022−62, J2041+46, J1246−65. Also rejected: CW Cep (free-free), and the LMC systems J0454 and J0537.
- **Gaia check:** all 408 Gaia ECL matches have `vari_eclipsing_binary` solutions. Where another survey also has a period, Gaia agrees for 96% (`src/crossmatch/gaia_veb_join.py`).

### 2. J1122 period was wrong (`docs/J1122_EPHEMERIS_CORRECTION.md`)

- Old P = 4.7909785 d → **new P = 4.7917045 ± 0.0000058 d, T0 = BTJD 3762.8427** (7-sector O−C).
- The old ephemeris shifted the S10/S11/S37 eclipses by −1.15/−0.92/+1.34 h. That is the likely cause of the earlier "sector-dependent eclipse bias".
- TIC lists 37% TESS contamination and Teff 7705 K for J1122; the earlier models assumed l3 = 0 and 6000 K.
- **7-sector PHOEBE fit:** χ²_red = 1.77, rms 1.86 ppt, i = 89.8°, r₁/r₂ = 0.099/0.033, T2/T1 = 0.57, fitted l3 = 0.16. It fails the binary-blackbody (JHKW1 χ² 116) and distance checks.
- J0719's ephemeris is fine: P = 1.0118518 ± 0.0000003 d.

### 3. PHOEBE results (7 finished; J1419 passes every check)

**J141909.42−565518.1 passes all checks:**
- eccentric EB: P = 16.634 d, e = 0.145, i = 89.2°, T2/T1 = 0.97;
- light-curve fit: χ²_red 0.92;
- blackbody/photometry: binary photosphere fits J, H, Ks, W1 (χ² 0.1); d_phot/d_Gaia = 0.91;
- IR excess: W4 at 8.5σ (W4-only, so the dust temperature is unconstrained); compact, centred W4 source.

Details are in `docs/DD_EB_CROSSMATCH.md`.

Fits that do not pass yet:

| Target | χ²_red | d_phot/d_Gaia | Verdict |
|---|---:|---:|---|
| J2312+53 | 4.4 | 0.84 | Closest to passing; small S-shaped eclipse residuals |
| J2218+54 | 47 | 0.65 | Shape good; structured eclipse residuals (O'Connell asymmetry, spot?); excess survives (18σ/14σ) |
| J1122 | 1.8 | 0.51 | Good LC fit; fails photosphere and distance checks |
| J0459+16 | 40 | 0.44 | Ellipsoidal shape needs fitted q |
| J0332+52 | 135 | 0.57 | Wrong local minimum (data show a total eclipse) |
| J2239+58 | 250 | 0.53 | Wrong local minimum |

**Systematic finding:** every system is 1.4–5× more luminous at its Gaia distance (well-measured parallaxes) than a main-sequence binary with the TIC Teff predicts. Possible causes:
- TIC Teffs biased low by reddening;
- evolved components;
- a luminous third star.

The third-star case would mean the IR excess may not belong to the EB. Results now include a `distance_scaling` block (implied luminosity, sma and mass).

## Running at handoff (local processes on the Mac)

Started ~19:50 UTC on 2026-10-02; expected to take ~2 h. The processes continue in the background on the Mac. They stop if the Mac sleeps or restarts; the cached grid and optimiser results allow a cheap relaunch.

| Run | Command tag | Output when done |
|---|---|---|
| J033226.58+520357.1, J045912.77+165543.8, J221843.61+544715.0, J231201.40+532028.6 | `--tag _v2` | `outputs/candidates/<name>/<name>_phoebe_result_v2.json` |
| J160415.99−562627.3, J232537.66+613847.9 (eccentric, q fitted; relaunched ~21:30 UTC from cached grids) | none | `<name>_phoebe_result.json` |

Check status and relaunch if needed:

```bash
cd ~/Documents/repos/tic-106588577-phoebe
pgrep -fl run_phoebe_candidate.py                      # still running?
find outputs/candidates -name "*_phoebe_result*.json"  # finished fits
tail -3 outputs/candidates/<name>/<name>_phoebe_v2.log # progress
# relaunch a dead run (cached grid/_localopt json files are reused automatically):
nohup ~/phoebe-env/bin/python src/pipeline/run_phoebe_candidate.py --designation <name> [--tag _v2] \
  > outputs/candidates/logs/<name>_stdout.log 2>&1 &
```

The round-2 code differs from round 1:
- q is fitted for P < 3 d;
- k = 1.3 is added to the grid;
- multi-start Nelder–Mead from the top 3 distinct grid points;
- dt0 starts from the measured eclipse phase;
- at least 2 × ndim emcee walkers;
- `_localopt*.json` cache.

## Next steps (priority order)

1. **Collect round 2.**
   - Run `.venv-tess/bin/python src/pipeline/summarize_phoebe.py`. It writes `outputs/candidates/phoebe_summary.md` and `phoebe_fits_montage.png` (the montage includes both `_v2` and round-1 results).
   - Update the "PHOEBE modelling status" section of `docs/DD_EB_CROSSMATCH.md` and say which candidates pass all checks (`passes_all` in the JSON).
   - Force-add the result JSONs and plots, plus the `.phoebe` bundles for passing models, then commit and push.
2. **Resolve the luminosity excess.** For each fitted system:
   - Re-derive Teff with an extinction-corrected SED (Gaia BP/RP + 2MASS with A_V free, e.g. from Bayestar/Green dust maps or a fit). Then rerun the PHOEBE fit with `--teff`.
   - If it is still overluminous, try `--fit-l3`.
   - A large fitted l3 together with an overluminous SED points to a third star, which would also affect whether the IR excess belongs to the EB.
   - **J2239+58 `_v2` is done and did not improve**: χ²_red 243; the primary is 0.047 too shallow (−55σ). Multi-start reaches the same minimum.
     - Likely cause: the deep (0.39) V-shaped primary needs the hotter star eclipsed by a larger, cooler companion near i = 90°.
     - Try `--eccentric yes`, since the secondary is at φ = 0.504 with unequal widths.
     - Seed T2/T1 lower, and let k go above 1.4 (the current bound).
3. **J1419−56 (passes every check):** look for higher-resolution 24 µm data (Spitzer MIPS, AKARI), check for apsidal motion with `eclipse_timing.py`, and test the result's sensitivity to the assumed Teff and q.
4. **J2218+54 (strongest SED):** add a cool spot for the O'Connell asymmetry and consider ld/gravity-darkening refinements. It is worth requesting RVs.
5. **J1122:** redo the custom-aperture extraction for all 7 sectors on the new ephemeris (`src/extract_multisector_tess_lightcurve.py --period 4.7917045 --t0 3762.8427`). Compare its eclipse depth with QLP and refit. Do not reuse the old `outputs/j1122_*` binned files; they use the wrong period.
6. **Improve the WISE image S/N** in `src/crossmatch/wise_cutouts.py`: use a plane-fit background in an 18–35″ annulus. J0719 falsely fails the current version. Then visually review the tier-B objects flagged `wise_image`.
7. **Fit the remaining tier-A objects:** J073937 (hot, no TIC Teff; use `--teff`), J024542, J173224, J223553, J220912.
8. **More surveys:** recent TESS-FFI EB catalogues (IJspeert+2024, Howard+2025), OGLE disk EBs, Pan-STARRS1 3π variables, K2 EBs.

## Pipeline commands

```bash
# cross-match and scoring (.venv-tess)
.venv-tess/bin/python src/crossmatch/build_dd_catalogs.py
.venv-tess/bin/python src/crossmatch/crossmatch_eb_surveys.py
.venv-tess/bin/python src/crossmatch/sanity_checks.py           # needs wise_cutouts csvs for the image check
.venv-tess/bin/python src/crossmatch/wise_cutouts.py [names]    # default: tier A
.venv-tess/bin/python src/crossmatch/gaia_veb_join.py
.venv-tess/bin/python src/crossmatch/make_tables.py
# per target
.venv-tess/bin/python src/pipeline/fetch_tess_lc.py --designation J... --ra .. --dec .. --tic .. --period ..
.venv-tess/bin/python src/pipeline/eclipse_timing.py --designation J...
~/phoebe-env/bin/python src/pipeline/run_phoebe_candidate.py --designation J... \
    [--teff K] [--fit-l3 --l3 0.3] [--eccentric auto|yes|no] [--fit-q auto|yes|no] [--tag _x]
.venv-tess/bin/python src/pipeline/summarize_phoebe.py
```

Target rows (ra, dec, TIC, period) for tier A are in `outputs/crossmatch/tierA_jobs.txt`. For targets outside the scored list, first run `src/crossmatch/support_rows_for_targets.py <designation>`.

## Pass criteria (`passes_all` in the result JSON)

| Check | Criterion |
|---|---|
| Fit quality | χ²_red < 5 |
| Eclipse bias | mean residual in each eclipse < 3σ |
| Binary-blackbody photosphere | fit to J, H, Ks, W1 with χ² < 60 |
| Distance | photometric/Gaia distance 0.6–1.8 |
| Excess | W3 or W4 excess ≥ 5σ against the binary photosphere |

## Gotchas

- `pkill -f` with "+" in a name silently matches nothing ("+" is a regex metacharacter); kill by PID.
- zsh errors on unmatched globs; use `find`.
- PHOEBE rejects short model labels ("m", "x").
- `~/Documents` may be iCloud-synced; `git add` occasionally hits "mmap failed", so retry per file.
- CDS XMatch can time out on 400k-row uploads; the script retries.
- `fetch_tess_lc.py`:
  - prefers SPOC, then TESS-SPOC, then QLP;
  - some QLP files lack `kspsap_flux` (falls back to `sap_flux` / `flux`);
  - doubles the period only if odd/even depths differ by > 4σ **and** there is no phase-0.5 secondary.
