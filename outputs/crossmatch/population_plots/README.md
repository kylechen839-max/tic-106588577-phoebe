# Population plots: all 800 DD1.0 × eclipsing-binary matches (2026-11-10)

Made on Kyle's PC at Dr. Bans's suggestion: plots that tie the 800 cross-matches together and show
whether they share characteristics. Source: `outputs/crossmatch/dd1_candidates_scored.csv`.
Script: `src/crossmatch/population_plots.py`. SIMBAD types are grouped into five classes, listed in
the script.

| file | what it shows | talking point |
|---|---|---|
| `01_cmd_by_simbad_class.png` | Gaia CMD (453 with parallax S/N > 5), coloured by SIMBAD class | The 800 are a mix: 382 binaries, 199 evolved / long-period variables (red giant branch, upper right), 59 YSOs (below the main sequence, lower right), 46 hot emission-line stars. All 12 tier A are catalogued binaries on the upper main sequence. |
| `02_cmd_by_tier.png` | Same CMD, coloured by tier | The checks remove the giants and keep hot, blue main-sequence stars. |
| `03_galactic_map.png` | Galactic l, b by tier | 75% of matches lie within 10° of the plane (tier A 83%). The LMC and SMC show up as clumps. |
| `04_dust_temp_vs_fir.png` | Dust temperature vs L_IR/L* by SIMBAD class | The key debris-disk diagnostic. 9 of 12 tier A sit at L_IR/L* < 10⁻² and 70-240 K, the debris-disk regime. Only 26% of all 800 do; YSOs and most binaries sit at 10⁻² to 1. |
| `05_period_histogram.png` | Catalogue period by tier | Tier A spans 0.96-16.6 d (median 2.6 d). Tier C carries the > 30 d tail, mostly giants. |
| `06_w4_excess_vs_teff.png` | W4 / photosphere vs photosphere temperature, by tier | Tier A are 7,000-18,000 K stars with modest W4 excess (1.7-25×). The huge excesses (100-10,000×) are cool stars at the 2,500 K fit limit: giants and YSOs. |
| `07_selection_funnel.png` | The real tier logic of `src/crossmatch/sanity_checks.py`, cumulative, stacked by SIMBAD class: 800 → 657 → 483 → 354 → 311 → 238 → 172 → 172 → **131** (all 8 critical checks = tier A + B), then → 105 → 55 → 53 → 13 → 12 → **12** (tier A) | The selection story in one chart; the colours show giants and wind stars disappearing first. The order within each group is a narrative choice; the end counts are not. |
| `09_first_failed_check.png` | The first check each object fails, by SIMBAD class (`first_failed_check.csv`) | What each check is for: the position match removes 94 giants with loose Gaia matches, the period check removes long periods, clean photometry removes W3/W4 upper limits, the Gaia-neighbour check removes crowded fields, the SIMBAD check removes Be/WR/emission-line stars, and the photosphere fit removes YSOs with near-IR excess. |
| `08_catalogue_overlap.png` | Number of EB catalogues per object | 58% of all matches (53% of tier A + B) are in ≥ 2 catalogues. VSX ingests other surveys, so these are not all independent. |

See `STORY.md` for the narrative of how the 800 were narrowed and classified.

Caveats:
- **Fit-grid limits.** Dust temperatures pile up at the grid limits of 40 K and 1500 K, and
  photosphere temperatures at 2500 K. Those points are limits, not measurements.
- **Magellanic Cloud stars in tier B.** Four tier B objects lie in the LMC: J045439.87-671306.0,
  J053749.52-664011.4, J054657.92-700933.2 and J050533.79-681109.5. Three have Gaia distances of 14-16
  kpc and W4 ratios of 300-700×. They are LMC OGLE binaries, almost certainly blended with LMC dust or
  star-forming regions. Tier B does not require the < 3 kpc check, so they were not removed. Worth
  dropping before quoting "119".
- **CMD coverage.** The CMD shows only the 453 stars with a good parallax.
- **Colours.** The validated default palette from the dataviz skill. The validator script could not
  be rerun here because Node is not installed.
