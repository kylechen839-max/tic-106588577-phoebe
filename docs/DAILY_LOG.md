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
