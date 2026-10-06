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
