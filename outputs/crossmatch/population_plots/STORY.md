# How we narrowed and classified the 800 candidates

A narrative for the talk, following the order of `07_selection_funnel.png`. Every number comes from
`outputs/crossmatch/dd1_candidates_scored.csv` and the check definitions in
`src/crossmatch/sanity_checks.py`. The tier rules are: **tier A + B** passes all 8 critical checks;
**tier A** also passes 6 stricter ones.

## 1. The question
Debris disks trace the leftover material of planet formation. Most stars form in binaries, but very
few circumbinary debris disks are known. Eclipsing binaries are the best place to look: their
inclination, radii and temperatures can be measured, so the stars' own infrared light can be predicted
and subtracted.

## 2. The starting sample: 800 matches
Disk Detective DD1.0 lists stars with a WISE 22 µm excess that volunteers judged to be clean point
sources. We matched them within 3″ against 16 eclipsing-binary catalogues: Gaia DR3, VSX, ASAS-SN, the
TESS EB catalogues, ZTF, ATLAS, CRTS, WISE, OGLE and others. That gave **800 matches**.

The 800 are not a clean sample (`01_cmd_by_simbad_class.png`). SIMBAD types them as:
- 382 eclipsing or spectroscopic binaries;
- 199 giants and long-period variables, many of which Gaia mislabels as eclipsing;
- 59 young stellar objects;
- 46 hot emission-line stars (Be, Wolf-Rayet, symbiotic);
- 114 with no useful class.

75% lie within 10° of the Galactic plane (`03_galactic_map.png`), where the 22 µm sky is crowded and
nebulous. So the job is to separate real binaries with real excess from everything that imitates them.

## 3. Are these really eclipsing binaries? (800 → 483)
- **Position match** (≤ 1.5″, or listed by ≥ 2 surveys): 143 removed, 94 of them giants. Gaia's
  variability catalogue matches bright red giants loosely.
- **Plausible EB period** (0.15-30 d; contact binaries < 1.5 d): 174 removed. These are the
  long-period tail and spurious periods (`05_period_histogram.png`).

## 4. Is the infrared excess real? (483 → 311)
- **Clean W3/W4 photometry** (good quality flags, no upper limits): 129 removed.
- **Excess ≥ 5σ** against two photospheres, one with free temperature and one at the catalogue Teff,
  so a biased temperature cannot create the excess: 43 removed.

## 5. Is the excess from this star? (311 → 172)
- **No bright Gaia neighbour** within 6″, and < 25% of the flux from neighbours within 12″: 73 removed.
  WISE W4 has a 12″ beam, so crowded fields are the main contaminant.
- **Usable in TESS** (Tmag < 14.5, TIC contamination < 0.3), because we need a light curve to model:
  66 removed.

## 6. Is it the right kind of star? (172 → 131)
- **Not a galaxy or QSO**: removes nothing.
- **Not a giant, AGB or wind star** in SIMBAD: 41 removed, mostly Be, Wolf-Rayet and emission-line
  stars. Their excess comes from gas, not dust.
- **131 objects pass all 8 critical checks: tier A + B.**

## 7. The strict cut for tier A (131 → 12)
- **Compact W3/W4 source** in the WISE image itself, not just in the catalogue: 26 removed.
- **Photosphere fits J, H, Ks, W1** (no near-IR excess): 50 removed, including 13 YSOs whose disks are
  hot and young.
- **Dust fits better than free-free emission**: 2 removed.
- **≥ 2 independent EB surveys** (VSX copies other surveys, so it does not count twice): 40 removed.
  This is the biggest single tier-A cut.
- **Gaia distance < 3 kpc**: 1 removed. At larger distances the 12″ W4 beam spans > 0.17 pc of
  interstellar medium.
- **Not a known YSO**: removes nothing further.
- **12 objects: tier A.**

## 8. What the survivors have in common
- **All 12 tier A are catalogued binaries on the upper main sequence** (`02_cmd_by_tier.png`):
  photospheres of 7,000-18,000 K, so B-, A- and F-type stars.
- **Periods of 0.96-16.6 days**, median 2.6 d (`05_period_histogram.png`).
- **Modest W4 excesses** of 1.7-25× the photosphere (`06_w4_excess_vs_teff.png`). The 100-10,000×
  excesses all belong to giants and YSOs that the checks removed.
- **9 of 12 sit in the debris-disk regime** (`04_dust_temp_vs_fir.png`): fractional dust luminosity
  below 10⁻² at 70-240 K. Only 26% of the 800 do. The selection keeps the debris-disk-like systems,
  rather than picking them out of nowhere.

## 9. From candidates to systems
- For the strongest candidates we fit the TESS light curves with PHOEBE, after validating the pipeline
  on published binaries. The fitted temperatures and radii then predict the binary photosphere for the
  SED.
- **J2312+53**: two A-type stars in a 6.85-day orbit. W4 is 3.6× the photosphere from the catalogue,
  and 2.5-2.8× after subtracting the nearby nebular halo. That is still > 10σ on photometric error
  alone, but only about 2.5-3σ once the background uncertainty is included.

## 10. Caveats to have ready for questions
- **Nebulosity.** Tier A stars sit in the Galactic plane, and J2312 and J0136 both sit on nebulosity.
  Background subtraction lowers but does not remove their excess. Higher-resolution IR imaging is the
  real test.
- **Tier B in the LMC.** Four tier B objects are LMC OGLE binaries; three have Gaia distances of 14-16
  kpc. Tier B does not apply the distance cut, so they are almost certainly blended with LMC dust.
- **Fit-grid limits.** Dust temperatures at 40 and 1500 K, and photosphere temperatures at 2500 K, are
  grid limits, not measurements.
- **Follow-up needed**: radial velocities for masses; higher-resolution mid-IR imaging (Spitzer
  archive, JWST); far-IR or sub-mm to find the cold dust.
