# Project Context and Handoff (DD × EB search, PHOEBE modelling)

Last updated: 2026-10-02. Written so that a new Claude session (Claude Code, or a claude.ai Project) can pick up this work without the original chat.

## How to recreate this session in a project

**claude.ai Project:**
1. Create a Project and paste the "Project instructions" block below into its custom instructions.
2. Upload these files as project knowledge:
   - `docs/PROJECT_CONTEXT.md` (this file)
   - `docs/DD_EB_CROSSMATCH.md`
   - `docs/J1122_EPHEMERIS_CORRECTION.md`
   - `docs/ECB_TESTDATA_LITERATURE_COMPARISON.md`
   - `docs/J1122_SECTOR_ISOLATED_MODELING.md`
   - `README.md`
   - `outputs/crossmatch/dd1_candidates_scored.csv` (optional, 600 KB)
3. Start the first chat with: "Read PROJECT_CONTEXT.md and continue from 'Open tasks'."

**Claude Code:** open the repo (`~/Documents/repos/tic-106588577-phoebe`). `CLAUDE.md` at the repo root points here and loads automatically.

### Project instructions (paste into the Project)

```
This project models eclipsing binaries (EBs) that are also Disk Detective (DD) infrared-excess
targets, using TESS photometry and PHOEBE 2.4.22. Goals: (1) find DD1.0 / DD2.x objects that are
catalogued EBs, (2) run sanity checks on each match (WISE quality, blackbody SED excess, dust vs
free-free, Gaia crowding, TESS contamination, object type, WISE images), (3) build PHOEBE models for
the candidates that pass the blackbody and photometry checks, with the established pipeline
(grid -> Nelder-Mead -> emcee -> high-res recompute -> post-fit photometry + binary-blackbody SED check).
Repo: github.com/kylechen839-max/tic-106588577-phoebe. Read docs/PROJECT_CONTEXT.md first.
Be explicit about what is verified vs. assumed; report chi2/residual numbers, not adjectives.
```

## Background

- Repo: `tic-106588577-phoebe`. The original target is TIC 106588577 = WISE J071951.40−240400.6 (J0719), a hot EB with a debris-disk-like W3/W4 excess. J112238.89−592027.5 (J1122) is the second target.
- The PHOEBE workflow (fold/bin → grid → bounded optimiser → emcee → high-res) was validated on V501 Mon, which recovered published parameters, and on W UMa (`docs/ECB_TESTDATA_LITERATURE_COMPARISON.md`).
- **DD1.0** = Disk Detective 1.0 subject lists `WZ_subjects_*.txt` (AllWISE + 2MASS photometry, IPAC tables), from `~/Downloads/WZSubs.zip`.
- **DD2.x** = `DDv2.1_Wise_Ids.txt` (AllWISE designations only), from `~/Downloads/DDv2.1_Wise_Ids.txt.lf`.
- `EB_WiseID.txt.lf` (7,827 IDs) is an older EB list of unknown origin. Only 6 entries are DD1 objects: J0719, J1122, J0136+62, J2039+46, J2052+44 and J2219+61.

## Environment

| Item | Value |
|---|---|
| Repo checkout | `~/Documents/repos/tic-106588577-phoebe`, branch `dd-eb-crossmatch-2026-10-02` |
| PHOEBE env | `~/phoebe-env/bin/python` (Python 3.11, phoebe 2.4.22, emcee, corner, scipy). **No pandas, no astroquery. Do not pip-install into it**, because it could upgrade numpy and break PHOEBE. |
| Analysis env | `.venv-tess/` in the repo, created from phoebe-env's Python: astroquery, lightkurve, pandas, matplotlib, tabulate (gitignored) |
| Machine | macOS, 12 cores, about 1–2.5 s per PHOEBE evaluation at 400 triangles with ~200 bins. Run up to ~11 fits in parallel. |
| Network services | CDS XMatch (VizieR tables), MAST (lightkurve, TESScut), SkyView (WISE cutouts) |

## Pipeline map

| Step | Script | Output |
|---|---|---|
| Build DD lists | `src/crossmatch/build_dd_catalogs.py` | `outputs/crossmatch/dd{1,2}_targets.csv` |
| Cross-match 16 EB surveys | `src/crossmatch/crossmatch_eb_surveys.py` | `dd*_x_<survey>.csv`, `dd*_eb_matches_long.csv` |
| Sanity checks and tiers | `src/crossmatch/sanity_checks.py` | `dd*_candidates_scored.csv`, `dd*_sed_plots/`, support caches |
| WISE image check | `src/crossmatch/wise_cutouts.py [names]` | `outputs/crossmatch/wise_cutouts/` |
| Tables | `src/crossmatch/make_tables.py` | `tierA_table.md`, `tierB_top40_table.md` |
| Gaia EB-pipeline check | `src/crossmatch/gaia_veb_join.py` | `dd1_gaia_veb.csv` (Gaia veb period, model, depths) |
| Extra targets (J0719, J1122) | `src/crossmatch/support_rows_for_targets.py` | `extra_scored.csv`, `extra_support_allwise.csv` |
| TESS light curves | `src/pipeline/fetch_tess_lc.py --designation --ra --dec --tic --period [--fixed-period]` | `outputs/candidates/<name>/` raw, binned, ephemeris, diagnostic |
| Eclipse timing (O−C) | `src/pipeline/eclipse_timing.py --designation [--old-period --old-t0]` | `<name>_timing.json`, `_oc.png` |
| PHOEBE fit | `~/phoebe-env/bin/python src/pipeline/run_phoebe_candidate.py --designation <name> [--teff --fit-l3 --l3 --eccentric auto/yes/no --tag]` | `<name>_phoebe_result*.json`, `_phoebe_best*.phoebe`, fit and SED plots, grid cache `_grid*.json` |
| Summary | `.venv-tess/bin/python src/pipeline/summarize_phoebe.py` | `outputs/candidates/phoebe_summary.md`, montage |

Everything under `outputs/` is gitignored. Selected results are force-added (`git add -f`), as earlier commits did.

## Key findings so far

1. **Cross-match:** 800 DD1 and 9 DD2 objects match at least one EB survey. After the sanity checks: **DD1 12 tier A / 119 tier B / 669 tier C; DD2 none** (all 9 DD2 matches are M dwarfs with W4 upper limits only).
2. **Best new candidate:** J221843.61+544715.0.
   - 1.669 d detached EB found by 4 independent surveys, 8 TESS sectors.
   - Photospheric through W2, then a ~160 K excess at 18σ (W3) and 14σ (W4), f_IR ≈ 4×10⁻³.
   - Compact W3/W4 source on the star.
3. **Rejected by the WISE image check (nebular contamination):** J1022−62, J2041+46, J1246−65. **Rejected as non-debris:** CW Cep (free-free/Be), the LMC systems J0454−67 and J0537−66 (14 kpc, beam confusion).
4. **J1122 period was wrong.** The old 4.7909785 d gives way to **4.7917045 ± 0.0000058 d** from 7-sector O−C. The old ephemeris shifted S10/S11/S37 eclipses by −1.15/−0.92/+1.34 h, which probably explains the "sector-dependent eclipse bias". TIC also gives J1122 37% TESS contamination and Teff 7705 K (earlier models used l3 = 0 and 6000 K).
5. **J0719 ephemeris is fine:** P = 1.0118518 ± 0.0000003 d (2×10⁻⁶ d from the README value).
6. **Gaia DR3:** the cross-match used the variability classifier (`ECL`). Every match has a `vari_eclipsing_binary` solution, and the Gaia period agrees with other surveys for 96% of them (aliases in the rest).

## Decisions and conventions (keep these unless there is a reason to change)

- Match radius 3″, with separations kept. Independent-survey count maps VSX entries back to their parent survey (Gaia DR3, ASASSN-V, ZTF, ATO, TIC, …).
- **Blackbody check:** the W3 or W4 excess must be ≥ 5σ against both a free blackbody (fit to J,H,Ks,W1) and a blackbody fixed at the TIC Teff, and the band must have AllWISE qph A/B. **Photometry check:** a catalogue-Teff photosphere fits J,H,Ks,W1 with χ² < 60. A 0.03 mag error floor is used throughout.
- **Tier A** requires every critical check plus photosphere, dust over free-free, at least 2 independent surveys, d < 3 kpc, not a YSO, and passing the WISE image check. A failed image check only demotes to B, because the automated S/N is fooled by nebular gradients (J0719 fails it but its image looks fine).
- **Light curves:** SPOC, then TESS-SPOC, then QLP pipeline products (TESScut fallback); per-orbit median normalisation; running-median outlier filter. Bins are sector-consensus medians, ~30 bins across each eclipse, uniform elsewhere. Period doubling only when odd/even depths differ by > 4σ **and** there is no secondary at phase 0.5.
- **PHOEBE:**
  - blackbody atmospheres, linear limb darkening (0.5, or 0.35 above 9000 K), q fixed at 0.8, M1 from a main-sequence Teff relation, sma from Kepler's law;
  - Roche distortion for P < 3 d, otherwise spheres;
  - parameters incl, dt0, r₁+r₂ (R/a), k = r₂/r₁, T2/T1, plus optional l3 and ecosω/esinω (auto-enabled when |φ₂ − 0.5| > 0.01);
  - local optimiser is **Nelder–Mead in normalised coordinates** (scipy bounded Powell was unreliable here);
  - emcee runs 12 walkers × 30 steps with errors inflated so the best χ²_red ≈ 1.
- **Post-fit "passes_all"** requires: χ²_red < 5, |eclipse residual bias| < 3σ, binary-blackbody fit to JHKW1 χ² < 60, photometric/Gaia distance 0.6–1.8, and W3/W4 excess still ≥ 5σ against the binary photosphere.
- Gotchas:
  - `pkill -f` with names containing "+" fails silently because "+" is a regex metacharacter. Escape it or kill by PID.
  - zsh errors on unmatched globs; use `find`.
  - PHOEBE rejects short model labels such as "m" and "x".

## Daily next steps (rewritten by each daily cloud run; last: 2026-10-06)

Cloud runs: if `.venv-phoebe` is missing, run `python3 -m venv .venv-phoebe && .venv-phoebe/bin/pip install phoebe==2.4.22 emcee corner scipy matplotlib numpy pandas tabulate`, then `.venv-phoebe/bin/python src/pipeline/build_tess_passband.py`. **One usage window per day** (Kyle, 2026-10-05): plan fits to finish within ~4.5 h. Do not relaunch after a container restart or once the limit resets; checkpoint and stop. Fits die when the session goes idle. The container can also restart mid-run; relaunch from the cached `_grid`/`_localopt` files. Archive hosts are blocked. See `docs/DAILY_LOG.md`.

**Candidates passing all checks (2026-10-06): J1419 (round 1), J1604 (`_gteff`), J2312 (`_gteffnoirr`).** Light-curve notes: J1419 (all checks in round 1; the converged `_v2` fit has a marginal +3.1σ secondary bias; fitting q does not fix it), J1604 (`_gteff`), J2312 (`_gteffnoirr`, robust to Teff ±10%).

1. ~~SED check for the three~~ **Done 2026-10-06:** `outputs/crossmatch/dd1_support_cache/allwise.csv` is now built offline from the DD1.0 WZ_subjects tables (`src/crossmatch/build_allwise_cache_from_wzsubs.py`; reproduces the J1419 round-1 SED check exactly). **J1604 `_gteff` and J2312 `_gteffnoirr` (and its ±10% Teff runs) now pass all checks.** WZSubs.zip is not committed (29 MB); unzip `~/Downloads/WZSubs.zip` to `data/catalogs/` to rebuild.
2. **DD1.0 labels, partial (Kyle, 2026-10-06):** `data/catalogs/DD_1.0ObjectsFromMAST.csv` (DD1.0 MAST table from Alissa Bans's Teams chat; 30,659 subjects, with volunteer fractions, SciTeamFollowUp, W1–W4/2MASS photometry) is now in the repo. `src/crossmatch/join_dd1_labels.py` writes `outputs/crossmatch/dd1_labels_join.csv`. Only 115 of 800 cross-matches are in it (22 tier B, 93 tier C, **no tier A**); all 22 tier B pass the DD1.0 excess cut and have majority "good", and 10 are sci-team YES. The WZSubs tables have photometry but no volunteer labels, so the 12 tier-A objects still have no labels; ask Kyle/Alissa for a fuller label export. With the WZSubs photometry, all 12 tier A, 116/119 tier B and 663/669 tier C pass the DD1.0 [W1]−[W4] excess cut. `w4rchi2`/`cc_flags`/`ext_flg` are still missing.
3. **J1604 robustness:** `--teff 12726 --tag _gteffm10` (stopped in Nelder–Mead on 10-06; grid cached) and `--teff 15554 --tag _gteffp10`. Each takes ~2–3 h.
4. **J1419 secondary bias:** q is ruled out (`_v3`). Inspect the secondary-eclipse residuals in `_phoebe_fit_v2.png`; try `--fit-l3 --powell-maxiter 12 --tag _v4`, or accept as marginal.
5. **J2218+54 (strongest SED):** best is χ²_red 11.8 at 9,500 K with q fitted (`_t9500`), bias −6.6/+7.6σ. Try `--teff 9500 --eccentric yes --tag _t9500ecc`. A spot (O'Connell asymmetry) needs a code option.
6. **J2239:** primary core still 8.5σ too shallow (`_noirr2`, χ²_red 11.5). Try `--no-irrad --kmax 2.5 --eccentric yes --fit-l3 --tag _noirr3`.
7. **Low priority:** J0739 `_v2` gave no improvement (stays a fail); J2235 `--powell-maxiter 12 --tag _v2` is still unrun.
8. **If the network opens:** new EB surveys (IJspeert+2024, Howard+2025, OGLE disk, ZTF), J1122 custom aperture, WISE plane-fit background.

## Status at handoff

Latest daily cloud results: `docs/DAILY_LOG.md` and the "Cloud run" section of `docs/DD_EB_CROSSMATCH.md`. See **`docs/HANDOFF.md`** for the Mac-session status (PHOEBE round-1 results, round-2 runs in flight, exact next steps). It supersedes the status and open-task lists below.

## Open tasks (priority order)

1. Finish the PHOEBE runs. Run `summarize_phoebe.py`, add the results to `docs/DD_EB_CROSSMATCH.md` and `docs/J1122_EPHEMERIS_CORRECTION.md`, and list which candidates pass the blackbody and photometry checks.
2. ~~Join the Gaia `ECL` matches to I/358/veb.~~ Done: all 408 have a veb solution; periods agree for 96% of those with another survey period (`src/crossmatch/gaia_veb_join.py`).
3. Improve the WISE image S/N, for example with a plane-fit background in a 18–35″ annulus. Then re-run cutouts and visually review every tier-B object flagged `wise_image`.
4. Look for more EB surveys:
   - TESS-FFI EB catalogues newer than Prša+2022 (for example IJspeert+2024, Howard+2025);
   - OGLE Galactic-disk EBs;
   - Pan-STARRS1 3π variables;
   - Kepler/K2 EBs (Kruse);
   - recent ZTF EB catalogues.
5. Fit the remaining tier-A objects: J073937 (no TIC Teff, hot), J024542, J173224, J223553 (long period, shallow), J220912 (RUWE 43, caution).
6. Treat the top tier-B objects: YSOs as a separate primordial-disk class, single-survey EBs as needing TESS confirmation of eclipses.
7. For J1122: redo the custom-aperture extraction on all 7 sectors with the new ephemeris, and compare QLP vs custom eclipse depth (dilution).
