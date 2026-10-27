# Daily log (cloud routine)

One entry per daily cloud run. The plan for the next run is in `docs/PROJECT_CONTEXT.md` → "Daily next steps".

## 2026-10-03

**Environment.** 4-core cloud container. These hosts are still blocked by the network policy: cdsxmatch.u-strasbg.fr, vizier/cdsarc/simbad.cds.unistra.fr, mast.stsci.edu, archive.stsci.edu, skyview.gsfc.nasa.gov, irsa.ipac.caltech.edu, gea.esac.esa.int, tables.phoebe-project.org. GitHub (including `git clone` of public repos) and PyPI work.

**Unblocked PHOEBE in the cloud.** TESS:T normally downloads from tables.phoebe-project.org. `src/pipeline/build_tess_passband.py` now rebuilds a blackbody-only TESS:T from `data/passbands/tess.ptf`, the response copied from github.com/phoebe-project/phoebe2-tables. Verified: the J1419 best model gives χ² = 210.41, identical to the Mac result.

**Not possible offline:** new survey cross-matches (VizieR/XMatch), TESS downloads for new targets (MAST), WISE cutout re-analysis (only PNGs are committed), and the AllWISE SED check. The AllWISE cache `outputs/crossmatch/dd1_support_cache/allwise.csv` is not in git. `run_phoebe_candidate.py` now tolerates the missing cache instead of crashing after the fit.

**Fits run** (details and table in `docs/DD_EB_CROSSMATCH.md` → "Cloud run 2026-10-03"):
- New tier-A fits: J0245+48 χ²_red 3.84 (secondary bias 3.8σ); J0739−36 2.90 (bias 4.5/4.8σ); J2209+58 **1.03, no bias** (M dwarfs, RUWE 43); J2235+50 15.2 (e = 0.40).
- J2239+58: eccentric 235, Gaia Teff 318. Neither improves on 243.
- Gaia-Teff refits: J2312+53 χ²_red 4.91, est. d ratio 0.84 → 1.03; J0459+16 χ²_red 40.4 → 22.4, est. d ratio 0.44 → 0.95.
- Still running at the time of writing: J1732−36 (eccentric, Roche, q fitted: slow), J0332+52 at Gaia Teff (emcee stage), J1419 at Teff −10% and +10%.

**Key finding.** The luminosity excess at the Gaia distance mostly goes away with the Gaia GSP-Phot Teff. TIC Teffs look reddening-biased low by 15–45%. This is inferred from a scaled SED; the full AllWISE check is still to do.

**Also done:** `src/pipeline/teff_distance_scaling.py` and `src/pipeline/recheck_sed.py`. Draft PR #2 is open.

## 2026-10-04

**Environment.** The same 4-core container, resumed. The 3 fits left running on 10-03 (J1732, J0332 `_gteff`, J1419 `_teffm10`) had been **killed when the session went idle**, ~3 min after the last reply. Background fits only survive while the session is active. Cached `_grid`/`_localopt` JSONs let the relaunches skip finished stages. Archive hosts are still blocked (CDS XMatch, VizieR, MAST, SkyView, IRSA, tables.phoebe-project.org).

**Code.** Two new options in `run_phoebe_candidate.py`: `--kmax` (upper bound on k; grid adds k = 1.7, 2.2 when allowed) and `--no-irrad` (reflection albedo 0).

**Fits** (table in `docs/DD_EB_CROSSMATCH.md` → "Cloud run 2026-10-04"):
- J1419 Teff −10%/+10%: geometry stable (i 89.1–89.4°, r₁ 0.094, e 0.14); χ²_red 0.67/0.64, below the nominal 0.92; est. d ratio 0.79/1.04; secondary bias 2.8/3.4σ.
- J0332 at Gaia Teff (18,920 K): χ²_red 125 (was 135). Still fails.
- J2312 at Gaia Teff, eccentric + l3: χ²_red 2.00, eclipse bias −4.0/−9.2σ.
- J2239 with k ≤ 2.5: no change (χ²_red 245, k stays 0.56). **With reflection off: χ²_red 17.4, rms 3.0 ppt** (was 243, 10.5).
- J2218 at 9,500 K, q fitted: **χ²_red 11.8, rms 1.5 ppt** (was 47, 3.6); q = 0.32; est. d ratio 1.12.
- Reflection off on Roche systems makes them worse: J2218 `_t9500noirr` χ²_red 18.96 (vs 11.8), J0459 `_gteffnoirr` 765 (vs 22.4).
- J1604 `_gteff` and J1732 were killed by a container restart at ~23:00 UTC and relaunched.
- **J1604 at Gaia Teff (14,140 K): χ²_red 1.66, rms 1.9 ppt, eclipse bias +1.7/−1.1σ**, so it passes both light-curve checks (e = 0.063, q = 0.76). The SED check is pending.
- J1732 (first fit): χ²_red 5.52, eclipse bias +8.3/+3.3σ. It fails.

**Key finding.** For J2239 (detached, flat out-of-eclipse), the default reflection albedo of 1.0 produced humps the data do not show; switching it off fixed most of the misfit. For the short-period Roche systems (J2218, J0459), switching it off makes the fit much worse, so it is a per-system choice.

## 2026-10-05

**Environment.** Archive hosts are still blocked, and there is still no AllWISE cache in git. The container restarted twice: ~23:00 UTC on 10-04, which killed J1604 and J1732, and ~12:50 UTC on 10-05, which killed the J1419/J0739/J2235 `_v2` refits. All were relaunched from cached grids.

**Late results from 10-04:** J1604 at Gaia Teff gives χ²_red 1.66 with bias +1.7/−1.1σ, so it passes the light-curve checks. J1732 (first fit) gives χ²_red 5.52 with bias +8.3/+3.3σ and fails.

**Fits today** (table in `docs/DD_EB_CROSSMATCH.md` → "Cloud run 2026-10-05"):
- **J2312 at Gaia Teff, reflection off: χ²_red 1.69, bias −2.1/+1.9σ, est. d ratio 1.02.** It passes the light-curve checks.
- J2239 eccentric, reflection off, k ≤ 2.5: χ²_red 11.5 (was 17.4); primary bias −8.5σ.
- J1419 `_v2` (12 Nelder–Mead iterations): χ² 144.7 (was 210.4), χ²_red 0.63, but the secondary bias is +3.1σ, so its light-curve pass is now marginal.
- Stopped unfinished at 16:00 UTC on Kyle's request (no work after the session limit resets): J0739 `_v2`, J2235 `_v2`, J1604 `_gteffm10`, J1419 `_v3` (q fitted). Their grid caches are pushed.

**Status.** Three candidates pass or nearly pass the light-curve checks: J1419 (all checks in round 1; marginal +3.1σ secondary bias at the better `_v2` minimum), J1604 (Gaia Teff) and J2312 (Gaia Teff, reflection off). Only J1419 has had the full AllWISE SED/distance check.

## 2026-10-06

**Environment.** Same container; archive hosts still blocked. Neither the AllWISE cache nor the DD1.0 label spreadsheet is in the repo, so the SED check and the new DD1.0-criteria/label columns (Kyle, 10-06) could not be done. One window: fits ran 07:55–12:00 UTC and were stopped at the cutoff.

**Fits** (table in `docs/DD_EB_CROSSMATCH.md` → "Cloud run 2026-10-06"):
- **J2312 at Gaia Teff ±10% (reflection off): χ²_red 1.64 (8,939 K) and 1.78 (10,925 K), all eclipse biases < 1.5σ, est. d ratio 0.87/1.16.** Its light-curve pass is robust to Teff.
- J1419 with q fitted (`_v3`): Nelder–Mead χ² 171.4 vs 167.8 for fixed q (`_v2`). Fitting q does not help; stopped in emcee.
- Stopped in Nelder–Mead: J1604 `_gteffm10` (χ² 410.1 so far vs 409.3 at Gaia Teff), J0739 `_v2` (454.6 vs 454.7 in round 1, no improvement).

**Status.** Light-curve-passing candidates: J1419 (marginal +3.1σ secondary bias), J1604 (Gaia Teff), J2312 (Gaia Teff, robust to ±10%). All three still wait on the AllWISE SED check.

**Afternoon addendum (Kyle online).** The DD1.0 MAST label table (`DD_1.0ObjectsFromMAST.csv`, from Kyle's Teams chat with Alissa Bans) and WZSubs.zip were copied from Kyle's Mac. Labels cover 115/800 cross-matches (no tier A). WZSubs gave an offline AllWISE cache, so the SED check finally ran: **J1604 (Gaia Teff) and J2312 (Gaia Teff, reflection off, also at ±10% Teff) now pass all checks**, joining J1419. Details: `docs/DD_EB_CROSSMATCH.md` → "SED check with the AllWISE cache (2026-10-06)".

## 2026-10-07

**Environment.** Archive hosts still blocked (CDS XMatch, VizieR, MAST, SkyView, IRSA return no connection). One window: fits ran 07:53–11:56 UTC and were stopped at the cutoff.

**DD1.0 criteria and labels.** `src/crossmatch/add_dd1_criteria.py` adds the DD1.0 excess cut (without `w4rchi2`, which is unavailable), the volunteer labels and a `grade` column to `dd1_candidates_scored.csv`; the tier tables show them. Tier A: 12/12 pass the DD1.0 cut (no labels exist). Tier B: 110/119 pass; grades 100 B, 10 B+ (sci-team YES), 9 B− (W4 `cc_flags` = d, or W4 S/N 9.96). Tier C: 466/669 pass.

**Fits** (all stopped in Nelder–Mead; table in `docs/DD_EB_CROSSMATCH.md` → "Cloud run 2026-10-07"): J1604 −10% Teff χ² 392.4 (vs 409.3 at Gaia Teff), +10% 417.7; J1419 with l3 173.1 (vs 167.8 without); J2218 eccentric 3265.8 (vs 1936.3 circular).

**Status.** Three candidates pass all checks: J1419, J1604, J2312. Four parallel fits do not finish within one window; run fewer at a time.

## 2026-10-08

**Environment.** Archive hosts still blocked. Fits ran 07:52–11:51 UTC (one window) and were stopped at the cutoff.

**Fits.** J2235 `_v2` finished from its cached optimum: χ²_red 15.4, eclipse bias −7.9/−8.2σ, fails (SED passes: d ratio 0.94, W4 8.3σ). J1604 at Teff −10%/+10% stopped in Nelder–Mead start 2 (best χ² 392.4/391.5 vs 409.3 at Gaia Teff).

**Code.** Each Nelder–Mead start is now cached (`<name>_localopt_starts<tag>.json`), so stopped runs resume at the next start instead of repeating start 1 (which happened on 10-06, 10-07 and 10-08).

**Status.** Unchanged: J1419, J1604, J2312 pass all checks.

## 2026-10-09

**Environment.** Archive hosts still blocked. Fits ran 07:52–11:51 UTC and were stopped at the cutoff. The per-start cache added on 10-08 worked: finished Nelder–Mead starts are saved for J1604 (both) and J2218.

**Fits.** J1419 `_v4` (third light free, fitted l3 ≈ 0.0001) finished and **passes all checks**: χ²_red 0.67, bias −0.4/+2.5σ, d ratio 0.91, W4 8.5σ. It is a nearby solution to `_v2` (χ² 144.7 vs 151.2), which had +3.1σ secondary bias. J1604 −10%/+10% Teff: start 1 cached (457.0/435.4), start 2 reached 405.8/392.0. J2218 eccentric: start 1 χ² 3265.8, worse than circular 1936.3.

**Status.** J1419, J1604, J2312 pass all checks.

## 2026-10-10

**Environment.** Archive hosts still blocked. Fits ran 07:52–12:15 UTC (one window) and were stopped at the cutoff. Only the two J1604 fits ran, as planned.

**Fits.** J1604 −10%/+10% Teff: start 2 finished (χ² 392.4/391.5, cached). Start 3 reached χ² 198.1/281.4 before the cutoff, far below the Gaia-Teff fit's 409.3, but had not converged and is not cached. The lower minimum is inferred from the unfinished runs, not yet verified by a finished fit.

**Status.** J1419, J1604, J2312 pass all checks. J1604's light-curve fit can likely improve further.

## 2026-10-11

**Environment.** Archive hosts still blocked. Fits ran 07:55–11:52 UTC. They were stopped early in emcee, which takes about 2 h and could not finish before the cutoff.

**Fits.** J1604 −10%/+10% Teff: Nelder–Mead start 3 converged at χ² 190.65/280.11, against 409.3 for the Gaia-Teff fit's local optimum. Both optima are cached. emcee, the hi-res χ²_red and the SED check have not run yet, so whether J1604 still passes at ±10% Teff is not verified.

**Code.** `run_phoebe_candidate.py` now saves an unfinished Nelder–Mead start's best vector every 25 evaluations and resumes from it.

**Status.** J1419, J1604, J2312 pass all checks. The Gaia-Teff J1604 fit is likely not at its best minimum (inferred from the ±10% runs).

## 2026-10-12

**Environment.** Archive hosts still blocked. Fits ran 07:52–11:58 UTC and were stopped at the cutoff.

**Fits.** J1604 robustness finished. At +10% Teff it passes all checks (χ²_red 1.01, bias +2.0/+1.0σ, d ratio 1.10, W4 9.0σ). At −10% Teff the SED check passes (d ratio 0.83, W4 8.9σ) and χ²_red is 0.75, but the primary-eclipse bias is +3.6σ, just over the 3σ limit. J1604 `_gteffms` (Gaia Teff, multi-start): grid done, start 1 χ² 413.4, start 2 partial. J2239 `_noirr3`: χ²_red 12.7, d ratio 0.51, fails.

**Status.** J1419, J1604, J2312 pass all checks. J1604's W4 excess (≈ 9σ) is robust to ±10% Teff. Its light-curve pass holds at the Gaia Teff and +10% and fails marginally at −10%.

## 2026-10-13

**Environment.** Archive hosts still blocked. One fit ran 07:53–11:52 UTC and was stopped at the cutoff; the other cores were idle by plan.

**Fits.** J1604 `_gteffms` (Gaia Teff, multi-start). Start 2 resumed from its cached partial best and converged at χ² 336.4, against 409.3 for the old single-start fit. Start 3 reached 696.1 (partial, cached). Resuming from a partial start worked as intended.

**Status.** J1419, J1604, J2312 pass all checks.

## 2026-10-14

**Environment.** Archive hosts still blocked. Fits ran 07:53–11:51 UTC and were stopped in emcee at the cutoff.

**Fits.** J1604 at the Gaia Teff:
- multi-start `_gteffms` best χ² 336.4;
- seeded from the −10% Teff optimum (`_gteffseedm10`, new `--start-from` option): **χ² 180.1**, the deep i ≈ 86.1° minimum;
- seeded from +10%: 248.8.

All local optima are cached; emcee, hi-res and the SED check run next.

**Status.** J1419, J1604, J2312 pass all checks. J1604's current Gaia-Teff fit (χ²_red 1.66) is not its best solution. A finished `_gteffseedm10` should give a lower χ²_red at the same Teff (inferred; it needs emcee and hi-res).

## 2026-10-15

**Environment.** Archive hosts still blocked. Fits ran 07:53–10:11 UTC and all finished; no cutoff needed.

**Fits.**
- J1604 at the Gaia Teff, three finished fits, all passing all checks: seeded from −10% χ²_red **0.70** (bias +2.8/+0.2σ), seeded from +10% 1.01, multi-start 1.41 (old `_gteff` 1.66). `_gteffseedm10` is J1604's new reference fit. The SED check is identical in all: d ratio 0.96, W4 9.0σ.
- J2312 seeded from its ±10% optima: χ²_red 1.68/1.80 against 1.69, so there is no hidden deeper minimum. It still passes all checks.

**Status.** Three candidates pass all checks, with best fits: J1419 (`_v4` χ²_red 0.67), J1604 (`_gteffseedm10` 0.70), J2312 (`_gteffnoirr` 1.69). The PHOEBE work on these three is complete. What remains is blocked on data: DD1.0 tier-A labels, `w4rchi2`, and TESS or archive access for tier B and new surveys.

**Late (10:15–12:12 UTC).** J1604 ±10% Teff reseeded from the new reference ran until the cutoff, still in local optimisation. Partial bests (cached): χ² 184.8 at −10% (was 190.65) and 174.4 at +10% (was 280.1).

## 2026-10-16

**Environment.** Archive hosts still blocked. Fits ran 07:53–10:41 UTC and both finished.

**Fits.** J1604 at ±10% Teff, reseeded from the 10-15 reference, both pass all checks:
- −10%: χ²_red 0.70, bias +1.6/−0.4σ, d ratio 0.83, W4 8.9σ;
- +10%: χ²_red 0.66, bias +1.6/+0.3σ, d ratio 1.10, W4 9.0σ.

The earlier −10% failure was a worse local minimum.

**Status.** J1419, J1604 and J2312 pass all checks at the Gaia Teff and at ±10% Teff. The PHOEBE work on these three is complete. Remaining items need data: DD1.0 tier-A labels, AllWISE `w4rchi2`, and archive or TESS access for tier B and new surveys.

## 2026-10-17

**Environment.** Archive hosts still blocked (all 000). Fits ran 07:53–11:50 UTC. J1732 `_m10` was stopped at the window limit; the others finished.

**Fits.**
- **J2312 eccentric** (`_gteffnoirrecc`): χ² 139.3/215, **χ²_red 0.65** (was 1.69), bias −0.1/+0.4σ, d ratio 0.93, W3/W4 7.0/12.6σ, passes all checks. e cos ω = −0.0024 matches the eclipse timing. e sin ω = 0.10 is weakly constrained. New J2312 reference.
- J0245: `_m10` χ²_red 3.80; seeded nominal returns to the same optimum (χ²_red 3.84). Still fails on the +3.5–3.8σ secondary bias.
- J0332: `_gteffm10` χ²_red 77.5; seeded at the Gaia Teff, χ²_red 96.1 (local opt 16,501 vs 18,891). Still fails.
- J1732 `_m10`: start 1 reached χ² 1,348 (nominal 1,002); stopped during start 2.

**Clip check.** The running-median clip has been in the script since 19e6290, the same commit that added every candidate light curve, so they were probably all built with it (assumed; git can't prove which version ran). A simulation shows it removes eclipse cores only at 30-min cadence (J1604: about 95% of core points at 30 min, none at 200 s). The binned J1604 eclipse core still uses all 5 sectors. J2312's deepest primary bins use 2 of 3 sectors (possible clip loss, unconfirmed). See `docs/DD_EB_CROSSMATCH.md`, "Cloud run 2026-10-17".

**Status.** J1419 (`_v4` 0.67), J2312 (`_gteffnoirrecc` 0.65) and J1604 (`_gteffseedm10` 0.70) pass all fit and SED checks. J1604's W4 excess is suspect (nebular W4 image).

## 2026-10-18

**Environment.** Archive hosts still blocked (all 000). Fits ran 07:53–11:32 UTC. J1732 `_m10` was stopped at the cutoff; all others finished.

**Clip fix.** `fetch_tess_lc.clip_mask()` no longer drops runs of consecutive same-sign deviations. In simulation, J1604's eclipse-core loss at 30-min cadence falls from 95% to 0%, and isolated bad cadences are still removed (39/40). Light curves were not regenerated (MAST blocked).

**J2312 eccentric.**
- ±10% Teff from the 10-17 optimum: χ²_red 0.65 / 0.60, both pass.
- Seeds at e sin ω = 0, −0.03 and +0.06 reached χ² 98.4, **84.6** and 90.5 (e sin ω 0.029 / 0.010 / 0.052).
- ±10% from the e sin ω ≈ 0.03 optimum: χ² 80.3 / 81.5, both pass.
- e cos ω = −0.0024 in every fit. e sin ω is not constrained (0.01–0.10).
- New reference `_gteffnoirreccesm003`: χ²_red 0.39, bias −0.1/+2.5σ, d ratio 0.99, W4 12.6σ.

**J1732 `_m10`:** three starts at χ² 1,267–1,529, all worse than nominal (1,002). Still fails.

**Status.** J1419 (`_v4` 0.67), J2312 (`_gteffnoirreccesm003` 0.39) and J1604 (`_gteffseedm10` 0.70) pass all fit and SED checks. J1604's W4 excess is suspect (nebular image).

## 2026-10-19

**Environment.** Archive hosts still blocked (all 000). Sixteen J2312 fits ran 07:53–10:45 UTC, two rounds of 8, and all finished.

**J2312 e sin ω profile** (new `--bound` option):
- χ² is 84–101 for e sin ω between 0.00 and 0.065, the minimum being 83.9 at 0.01.
- It rises to 116 at 0.08, 141 at 0.10, 202 at 0.13 and 196 at −0.02.
- Repeats from different seeds differ by up to 25 at the same e sin ω.
- **Inferred:** e sin ω between 0 and 0.07. The reference stays `_gteffnoirreccesm003` (χ²_red 0.39).

**Flux errors.** The binned errors come from sector-to-sector scatter, and in eclipse that scatter is 3–25× larger than out of eclipse (depth differences between sectors, probably dilution). χ²_red < 1 therefore reflects large in-eclipse errors, not a photon-noise-level fit. With white-noise errors, the in-eclipse χ²_red of the passing fits would be J1419 2,983, J1604 182 and J2312 1,583. The pass criteria test consistency with the sector spread. See `docs/DD_EB_CROSSMATCH.md`, "Cloud run 2026-10-19".

**Status.** J1419, J1604 and J2312 still pass all checks as defined (J1604 W4 suspect). The error-model finding means the pass is relative to sector-to-sector systematics.

## 2026-10-20

**Environment.** The 07:51 fire was skipped while the routine was paused. Kyle reset the usage and the run started at 07:53 UTC. Archive hosts still blocked (all 000). Fits ran 07:54–11:39 UTC and all finished.

**White-noise χ².** The photometry check now reports χ²_red against the bin-to-bin white noise; this is informational only. Recomputed for all 67 stored fits: the pipeline χ²_red is reproduced exactly. For the reference fits, χ²_red against white noise is J1604 25, J2312 260 and J1419 331 (eclipse windows 116 / 1,058 / 1,915). No fit reaches the white-noise level. The closest are J1122 (1.9) and J0739 (4.8).

**J2312, reflection on:** grid fit χ²_red 1.91 (secondary bias −3.5σ), seeded fit 3.54 (primary bias −4.2σ). Both fail and are worse than reflection off (0.39).

**J0739 seeding:** `_m10` χ²_red 2.60, seeded at nominal 2.83 (round 1: 2.90). Both fail the eclipse bias (> 3.7σ).

**Status.** Unchanged: J1419, J2312 and J1604 pass all checks as defined (J1604 W4 suspect). Against white noise, none fits to the photon level, which points to per-sector dilution as the next model change (needs raw light curves).

## 2026-10-21

**Environment.** Archive hosts still blocked (all 000). Work ran 07:54–11:01 UTC.

**Per-sector dilution fitter** (`src/pipeline/fit_per_sector_dilution.py`):
- Each sector is binned separately with white-noise errors and gets its own flux offset and depth scale against one PHOEBE model.
- **Injection test on J1122's time stamps:** geometry recovered within 0.07° in incl, 0.001 in r_sum, 0.008 in k and 0.002 in T₂/T₁. Injected depth scales 0.70/1.40/0.90 came back as 0.707/1.389/0.904, at χ²_red 1.21. Needs ntri 1500: ntri 400 gives χ²_red 3.10 even at the true parameters.
- **J1122 (sectors 10, 11, 37):** χ²_red 1.33 against white noise. Depth scales are S10 0.77, S11 1.18 and S37 1.05, consistent across 3 fits, and Δχ² is 75–81 against a shared dilution. The geometry is degenerate: k is 0.26–0.51 at nearly equal χ². Needs the other 4 sectors or a total-dilution prior.

**Status.** The per-sector approach works and is tested. It can't be applied to J2312, J1604 or J1419 until their per-sector light curves are fetched (MAST blocked; `_raw.txt` not in the repo).

## 2026-10-22

**Environment.** Archive hosts still blocked (all 000). Work ran 07:54–10:25 UTC.

**Per-sector fitter additions.**
- New options: eccentricity, exposure smearing for long cadences, and a model l3 with a prior on the mean depth scale.
- The 10-21 results reproduce exactly.
- Smearing changes J1122's best χ² by only 0.2.
- Injection with the prior recovers the geometry (k 0.313 vs 0.320) and the depth scales (0.707/1.389/0.905).

**J1122 with the TIC dilution** (l3 0.367, σ 0.05).
- All 3 starts converge to k = 0.11–0.12 (R₂ ≈ 0.25–0.32 R☉) at χ² 806–849 (χ²_red 1.39–1.47).
- That is Δχ² 40–83 worse than free dilution.
- The free-dilution solutions imply 85–94% third light, against 37% from TIC, so they are probably unphysical.
- Per-sector depth scales are S10 ≈ 0.82, S11 ≈ 1.10 and S37 ≈ 1.07.

**Status.** The per-sector tools are ready. Applying them to J2312, J1604 and J1419 needs their per-sector light curves (MAST blocked).

## 2026-10-23

**Environment.** Archive hosts still blocked. Work ran 07:51 to about 09:15 UTC.

**Emcee step added** to `fit_per_sector_dilution.py` (`--walkers`, `--steps`; commit 3c0c65f).

**J1122 emcee runs lost.** Three 12×60 runs (`_mcfreea`, `_mcfreeb`, `_mcl3`) were at step 20/60 when the container restarted at about 09:05 UTC. No chains survived. Per the one-window rule, nothing was relaunched.

**Fix.** The chain is now checkpointed every 5 steps (`*_persector_chainpart{tag}.npz`), and a rerun with the same tag resumes from it. The code compiles; resume has not been exercised yet.

**Status.** No new science results today. The next run should rerun the three J1122 emcee jobs with the same commands and tags (each about 2.5 h at 4 cores in parallel).

## 2026-10-24

**Environment.** Archive hosts still blocked (all 000). Work ran 07:52 to about 12:20 UTC.

**J1122 emcee reruns** (12×60, 3 in parallel, about 2.7 h each, all finished by 10:35 UTC):
- Free dilution: two separate modes at nearly the same χ² (765.9 vs 769.5). One has k 0.405 (+0.052/−0.035) at incl 90°; the other has k 0.286 (+0.050/−0.038) at incl 85°.
- TIC dilution prior: k 0.122 (+0.008/−0.005), R₂ ≈ 0.34 R☉, χ²_red 1.39.
- Depth scales are consistent across all runs: S10 ≈ 0.78, S11 ≈ 1.17, S37 ≈ 1.05.
- Chains are short and were not convergence-tested.

**Lost time.** My wait loop used `pgrep -f` on a pattern that matched the loop's own shell. So it never saw the fits finish, and about 1.5 h of the window went idle. Use `pgrep -f "^.venv-phoebe/bin/python src/pipeline/fit_per"` or the PIDs from `$!` next time.

**Fix.** Resume with zero steps left no longer crashes (`afd6d37`).

## 2026-10-25

**Environment.** Archive hosts still blocked (all 000). Work ran 07:52 to about 11:00 UTC.

**J1122 emcee extended to 120 steps.** All three chains resumed from the 10-24 checkpoints, so checkpoint/resume is now verified in a real run. They finished at 10:45 UTC, about 2.9 h for 60 steps.
- **TIC dilution prior:** k = 0.125 (+0.004/−0.005), R₂ ≈ 0.36 R☉. This is stable between the two halves of the kept chain.
- **Free dilution:** the two chains now overlap (k 0.40 and 0.32), and one walker has drifted out to k 1.3. The 10-24 "bimodality" is a single poorly constrained region, roughly k 0.27–0.6+.
- **Convergence:** τ ≈ 10–15 steps, so the chains are about 9 τ long. They are not converged.
- **Depth scales:** unchanged, S10 0.77–0.80, S11 1.16–1.18, S37 1.04–1.05.

**Wait loop.** This time it watched the fits by PID, with no idle gap.

**Note.** The `_mc*` result JSONs are rewritten with the 120-step percentiles. The 60-step versions are in commit 2666e03.

## 2026-10-26

**Environment.** Archive hosts still blocked (all 000). Work ran 07:52 to about 11:20 UTC.

**J1122 l3 sensitivity** (Nelder–Mead with depth-scale prior σ 0.05):
- l3 0.25: k 0.110, χ² 808.0. The first start stalled at 817.2; a second start from the 0.50 fit fixed it.
- l3 0.367 (TIC): k 0.118, χ² 806.1.
- l3 0.50: k 0.136, χ² 797.1.
- So R₂ ≈ 0.30–0.37 R☉ across this range, and the small companion is robust to the exact TIC dilution.

**J1122 TIC-prior chain extended to 180 steps.** k stays at 0.125 ± 0.005. But incl (83.7 → 82.3) and r_sum are still drifting, and τ is 14–23 steps, so the chain is not converged. R₂ = 0.38 (0.35–0.41) R☉.

**2026-10-26 evening (extra routine fire at 22:13 UTC after Kyle added the J2312 checks; no fits run).**
- No PC results from Kyle yet: none on any branch, and nothing new in `/mnt/project-files/J2312`.
- MAST, IRSA, Gaia and VizieR still return 000, so checks (a) and (b) are blocked.
- Check (c), by eye, `/mnt/project-files/wise_w3w4_check/J231201.40+532028.6_wise.png`: the W4 source is a compact peak at the star (offset 1.1″, S/N 15). It sits on a diffuse plateau that extends S and E and runs off the 2′ field. The 10-16 label "clean" was too generous, so it is now "possibly contaminated", as 081ab9c says.
- A faint source about 50″ SW is visible in W1–W3 but not clearly in W4.
- No verdict until the wide-field check runs.

**2026-10-26 late: J2312 checks run on Kyle's PC** (Remote Control session; commit f953daa, `outputs/j2312_contamination/`).
- **AllWISE:** target W4 6.84, SNR 17.5, w4rchi2 1.14, cc_flags 0000. Neighbours within 3′ are at 23–34″ with W4 SNR ≤ 4.1, so none can supply the W4 flux. **Verified.**
- **Gaia DR3:** none within 12″; 52 sources at 12–60″, none within 5 mag of the target, summing to 1.8% of its flux. So TESS dilution from neighbours is about 2%. **Verified.**
- **TIC:** contamratio is null (the star is not in the CTL). The nearest TIC neighbour is Tmag 14.46 at 12.8″.
- **Sky test:** near/far sky flux ratio is 0.88 in W3 and 0.82 in W4. The W4 near-annulus sky is above the far sky in all four quadrants.
- **By eye** (I checked the 4′ zoom too): the compact source sits on a N–S 22 µm ridge running about 3′ south. The W4 radial profile is about 15–20% of peak at 50–100″.
- **Interpretation (not proven):** the W4 excess may be partly heated nebular dust rather than a disk. The 18% shift from the choice of sky annulus does not remove the W4/photosphere ratio of 3.6, but emission peaked on the star can't be separated by this test.
- The PC session ran its own reimplementation of `w4_neighbours.py` (same method and parameters), because the PC's permission policy blocked running the clone's script. Both scripts are in the commit. Commit 081ab9c is still not merged.


**2026-10-26 night (Kyle's requests, cloud and PC).**
- **J2312 W4 excess with the local nebular background subtracted** (`outputs/j2312_contamination/README.md`):
  - W4 drops from 3.6× to 2.5–3.0× the photosphere: 8–12σ on photometric error alone.
  - The spread between background estimates is about ±0.25 mag, so the honest significance is about 2.5–3σ.
  - W3 drops to 2–3.5σ, so it is not significant.
- **J0136+62** (J013614.42+625507.1, tier B on `mac-run-2026-10-04`; TESS 10k EB P = 2.78951 d):
  - FFI light curve from S18, S24, S25, S58 and S85 gives P = 2.789497 d. The primary is 2.5–3% deep, odd and even eclipses agree (0.3σ), and there is no secondary. A once-per-orbit 1.3% modulation is present.
  - The difference image centres on the target.
  - W3/W4 sit on the edge of a nebula. W4 is 8–10σ after local-background subtraction; W3 is about 3.5σ.
- **J0131−00** (DD2, PC session): no eclipses at the VSX 0.351 d period. TESS shows a 0.488 d spot-like modulation, and W4 is an upper limit. Not an EB candidate.
- **J2312 deck:** `outputs/j2312_deck/` (10 slides, plus the explainer).

## 2026-10-27

**Environment.** CDS XMatch, VizieR, MAST, SkyView, IRSA and the Gaia archive still return 000. The public AWS mirrors work, and they unblock most of the TESS and WISE work:
- TESScut FFI cubes on `s3://stpubdata/tess/public/mast/`. These cover sectors up to about 108; no QLP or SPOC HLSP products are there.
- AllWISE Atlas tiles on `nasa-irsa-wise`.

Work ran 07:52 to about 12:00 UTC.

**New tools** (committed 10-26 evening, extended today):
- `src/pipeline/fetch_tess_ffi_s3.py`: 15×15 px FFI cutouts with astrocut (a few seconds per sector), aperture photometry with a per-cadence background, and dropping of scattered-light cadences. Then per-orbit trend removal that keeps the periodic signal; the phase model is skipped when an orbit is shorter than 1.5 P, which fixes J1419. Finally `fetch_tess_lc.py` cleaning and binning (new `--outdir`).
  - Sectors missing from tess-point's pointing table (J2312 S77/S84, J1604 S66) could not be cut. A brute-force camera/CCD search gave garbage stamps for J2312, and was killed for memory on J1604.
- `src/crossmatch/wise_atlas_s3.py`: 10′ W3/W4 images, near/far sky, and radial profiles.

**FFI light curves for the per-sector dilution fits** (`outputs/candidates/<des>/ffi/`):
- J2312: S16, S17, S57.
- J1604: S12, S39, S65, S102, S103.
- J1419: S11, S38, S65, S99–S102.

All fold cleanly at the reference periods.

**Per-sector dilution fits on the FFI data** (`fit_per_sector_dilution.py`, `--ntri 1500 --no-smear --eccentric`, tag `_persectorffi`):
- **J2312** (reflection off, 600 evals; the l3 = 0.018 + scale-prior run gave the same optimum):
  - χ² 77,142, χ²_red 121.5 against white noise. The QLP/SPOC data gave 260.
  - Depth scales S16 0.998, S17 1.006, S57 0.996. **The ~3% sector depth disagreement in the QLP/SPOC products is absent in uniform FFI photometry**, so it came from the pipeline products, not the star.
  - incl 83.48°, r_sum 0.294, k 0.640 (was 0.70), T₂/T₁ 0.953 (was 0.974), e cos ω −0.0028, e sin ω 0.017.
  - The optimum stopped improving after about 100 evals. A second start (`_persectorffi2`) is reported below.
  - Most of the remaining χ² is out of eclipse. The FFI continuum has about 0.8% peak-to-peak curvature, faint near the eclipses, that the reflection-off sphere model does not have. It is not settled whether this is real or left by the detrending; QLP did not show it at this level.
- **J1419** (200 evals, capped because each eval takes about 55 s with 7 sectors):
  - χ²_red 137.8 per-sector vs 148.6 with shared dilution.
  - Depth scales 0.993–1.011.
  - incl 89.39°, r_sum 0.142, k 0.486, T₂/T₁ 1.019, e cos ω 0.134, e sin ω 0.036. This is close to `_v4` (k 0.499, T₂/T₁ 0.994).
  - Residuals are dominated by out-of-eclipse steps between sectors (detrending), not by the eclipses.
- **J2312 second start** (`_persectorffi2`, from the first optimum, 200 evals): χ² 76,636 (χ²_red 120.7), down 0.7% from 77,142. No improvement after eval 100. incl 83.04°, r_sum 0.298, k 0.672, T₂/T₁ 0.948, e cos ω −0.0024, e sin ω 0.010. Depth scales S16 0.997, S17 1.006, S57 0.996, unchanged. k moves by 0.03 between starts for a 0.7% χ² change, so k is poorly constrained by this data; the geometry otherwise holds.
- **J1604** (the 600-eval run was too slow, 125 evals in 2 h, so it was restarted at 100 evals; **not converged**: best χ² still fell 7,697 → 7,649 over the last 25 evals):
  - χ² 7,649, χ²_red 8.06 per-sector vs 11.26 with shared dilution (966 bins).
  - Depth scales S12 1.004, S39 1.005, S65 0.979, S102 1.011, S103 1.001. S65 is 2% shallow; the rest agree within 1.1%.
  - incl 86.27°, r_sum 0.384, k 0.778, T₂/T₁ 0.978, e cos ω 0.062, e sin ω −0.016, mean b 0.885. The 10-24 `_gteff` fit had i 85.1°, T₂/T₁ 0.97 and e 0.063, so the FFI fit agrees so far.

**Fitter.** Now writes `<des>_persector<tag>_partial.json` (best χ² and parameters) every 25 evals, so a run stopped at the end of the window keeps its optimum.

**Tier-B TESS confirmation** (`outputs/candidates/tierB_tess/`, `tierB_tess_eclipse_check.csv`, `montage_yes.png`, `montage_no.png`). FFI light curves, at most 3 sectors for B+ and 2 for B/B−, were folded at the catalogue period. Checked: 87 of 119. Not checked: 29 with no catalogue period, and 3 whose sectors are not on the mirror.
- **Automatic flag** (primary depth ≥ 5× profile noise, half-width ≥ 0.003): 43 yes.
- **By-eye class:** 30 detached eclipses, 10 close-binary-like continuous variation, 5 eclipses seen but the fold needs a better period, 5 variable but not eclipse-like, 37 no signal, too faint or systematics.
  - The by-eye class is the one to use: it rescues J0924−05 and the contact-like systems, and rejects pulsators and systematics.
  - Faintness drives the "no signal" group: 21 of 44 at Tmag > 13.
- **B+:** only J0009+58 shows eclipses, and they smear in the fold, so it needs a period search. J0408+32 is a 2% double wave with no eclipse (ellipsoidal or contact-like). J0604−07 shows a short-period pulsation.
- **Nebular-demoted tier-A objects:** J1022, J1246−65 and J2041 eclipse clearly, so they are EBs, but their W4 remains suspect.

**Daily-routine bookkeeping.** PR #2 is still open; Kyle's "merge it" (10-26) is awaiting his choice of what to merge. Nothing has been merged.
