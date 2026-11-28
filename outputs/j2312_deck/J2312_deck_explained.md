# J2312 deck, slide by slide

This guide explains each slide of `J2312_research.pptx`: what the figure shows, what each number means and where it comes from.

The target is J231201.40+532028.6, called **J2312** for short. It is a bright star (Gaia G = 8.73) in Cassiopeia. It came out of our cross-match of the Disk Detective (DD 1.0) catalogue with eclipsing-binary catalogues: it is both an eclipsing binary and a star with too much infrared light.

---

## Slide 1: Title

The slide shows only the title and the two catalogues that list J2312 as an eclipsing binary. For reference, its TESS Input Catalog ID is TIC 252741805 and its Gaia brightness is G = 8.73.

- **Eclipsing-binary catalogues.** J2312 is listed as an eclipsing binary in two catalogues (our cross-match, `outputs/crossmatch/dd1_eb_matches_long.csv`):
  - **Gaia DR3 eclipsing binaries** (`vari_eclipsing_binary`, source_id 1995545679513855744). Gaia's model gives P = 6.8494 d, class ECL, matched at 0.03″.
  - **VSX** (AAVSO Variable Star Index) as **PMAK V105**, type **EA** (Algol-type, detached), P = 6.8495 d, matched at 0.05″.
  - Gaia DR3's variability classifier also labels it ECL. VSX gives its spectral type as A2, which fits the ~9,900 K temperature.
  - A VizieR search of 18 catalogues within 5″ (run 26 Oct from Kyle's PC) found it in no others. That search covered ASAS-SN, ZTF, ATLAS, the TESS EB catalogues (Prša+2022, IJspeert+2021), GCVS and SB9. At G = 8.7 it is probably too bright for most ground surveys; this is a guess, not checked.

## Slide 2: How we found it (cross-match)

The four cards are the steps of the Disk Detective × eclipsing-binary cross-match (details in `docs/DD_EB_CROSSMATCH.md`):
- **412,729**: all Disk Detective 1.0 objects (AllWISE + 2MASS photometry).
- **800**: those within 3″ of an eclipsing binary in 15 catalogues: 14 from VizieR (Gaia DR3 classifier, VSX, ZTF, ATLAS, ASAS-3, the TESS EB catalogues, OGLE, CRTS, WISE, NSVS, Kepler, DEBCat) plus ASAS-SN, matched locally.
- **12**: tier A, after scoring each match on 18 checks. The checks cover:
  - confirmation by ≥ 2 independent surveys;
  - a plausible, consistent period;
  - clean WISE photometry;
  - a W3/W4 excess above a stellar photosphere, better fitted by dust than by free-free emission;
  - no confusing neighbours;
  - usable TESS data;
  - not evolved or a wind source;
  - nearby;
  - a compact W3/W4 source in the images.
  119 matches are tier B and 669 tier C.
- **J2312**: one of the 12 tier-A objects, scoring 18 of 18.

## Slide 3: Summary

There are three headline numbers.

- **6.85 d** is the orbital period: the two stars circle each other every 6.85 days.
  - "Detached" means the stars are well separated and not exchanging mass.
  - "~9,900 K" is their surface temperature, from Gaia DR3 (GSP-Phot 9,932 K). That makes them hot A-type stars, hotter than the Sun's 5,800 K.
- **χ²ᵣ 0.39** is the reduced chi-squared of our PHOEBE model fit to the TESS light curve.
  - PHOEBE is the eclipsing-binary modelling code.
  - χ²ᵣ ≈ 1 means the model matches the data to within the error bars, and below 1 means it matches better than the error bars.
  - The low value partly reflects how we estimated those errors; slide 5 explains this.
- **12.6σ** is how far the measured 22 µm (WISE W4) flux sits above what the two stars alone should emit. It is 12.6 times the measurement uncertainty, which is far beyond chance.

**Bottom line:** the extra infrared light is real and does not come from a neighbouring star. However, the star sits on a patch of extended dust emission. The excess may therefore be from interstellar or nebular dust rather than a disk around the binary.

## Slide 4: Period

**Figure.** TESS light curves from sectors 16, 17 (QLP pipeline) and 57 (SPOC pipeline), folded on the period and stacked eclipse by eclipse.
- Blue curves are the primary eclipses, when the hotter or larger star is blocked.
- Red curves are the secondary eclipses, half an orbit later.

**Why 6.85 d and not 3.42 d.** If the two eclipses were identical, the true period could be half as long (3.42 d), with every dip being the same eclipse.
- The primary is deeper than the secondary by about 0.03 in normalised flux in sectors 16 and 17, and about 0.014 in sector 57.
- Combined over the sectors, this difference is 118σ.
- Two different eclipses mean the period is 6.85 d (6.849530 d).

**"Matches both catalogue periods."** Our 6.849530 d agrees with Gaia DR3 (6.8494 d) and VSX (6.8495 d). So neither catalogue lists the half period.

**"Secondary is 16 min early."** On a perfectly circular orbit the secondary falls exactly at phase 0.5. Here it falls at phase 0.4984, which is 16 ± 1 minutes early. That means the orbit is very slightly elongated (eccentric).

## Slide 5: Binary fit

**Figure.** The TESS light curve (black points, sectors combined and binned) with the PHOEBE model on top (red). The lower panel shows the residuals (data minus model). They are flat, so the model leaves no systematic misfit.

**Table:**
- **Inclination 82.8°** is the tilt of the orbit; 90° would be exactly edge-on. At 82.8° the stars partly eclipse each other.
- **R₁ / R₂ = 4.28 / 2.95 R☉** are the radii of the two stars in solar radii. Both are larger than main-sequence A stars of this temperature, so they are slightly evolved.
- **T₁ / T₂ = 9,932 / 9,700 K** are the surface temperatures.
  - T₁ is fixed to the Gaia value.
  - T₂ is fitted, as the ratio T₂/T₁ = 0.977.
  - "Near-twin" refers to these nearly equal temperatures.
- **Eccentricity ≈ 0.01** measures how elongated the orbit is; 0 is a circle.
  - It comes from the fitted e·cos ω = −0.0023 and e·sin ω = 0.010.
  - Allowing it improved the fit a lot.
- **χ²ᵣ 0.39** means χ² = 84.6 over 215 degrees of freedom.

**Assumptions in the model:**
- mass ratio q = 0.8;
- primary mass 2.5 M☉;
- no third light.

The light curve alone cannot pin down the masses. Radial velocities would be needed for that.

**Caveat line.** The error bars in the fit are the scatter *between* TESS sectors, not photon noise.
- The three sectors disagree on eclipse depth by up to ~3%, probably because of different pipelines and background light.
- Against photon noise alone χ²ᵣ would be about 260.
- The fix is a per-sector fit, with a separate dilution term for each sector. We built and tested that on J1122. For J2312 it needs the raw light curves from MAST, which the cloud setup can't download.

## Slide 6: Infrared excess (SED)

**Figure.** The spectral energy distribution (SED): brightness against wavelength.
- Black points are the measured photometry:
  - 2MASS J, H, Ks at 1.2–2.2 µm;
  - AllWISE W1–W4 at 3.4, 4.6, 12 and 22 µm.
- The blue curve is the light the two stars should emit. It is computed from the PHOEBE radii and temperatures and scaled to the Gaia distance.

**Points:**
- **Near-IR matches the photosphere.**
  - J, H, Ks, W1 and W2 sit on the blue curve, so the stellar model is right.
  - The photometric distance agrees with Gaia to within 1% (ratio 0.99).
- **W3 7.0σ, W4 12.6σ.** At 12 µm and 22 µm the star is brighter than the curve allows.
  - At 22 µm it is about 3.6 times brighter (slide title).
  - The excess holds if the temperatures are changed by ±10% (W4 stays at 12.4–12.8σ).
  - Warm dust is the usual cause: a debris or circumbinary disk, or nearby nebular dust.

## Slide 7: AllWISE sky images

**Figure.** 10′ × 10′ cutouts of the AllWISE all-sky images, centred on J2312, with north up and east left.
- Left: W3 (12 µm). Right: W4 (22 µm).
- The cyan circle is the 16.5″ photometry aperture on the star, and the cyan square is a 2′ box.
- The dotted rings are the near and far sky annuli used to compare background levels.
- The green circles are other AllWISE sources with W4 SNR ≥ 3. All are faint (W4 SNR ≤ 4.1).

**What it shows.**
- The star is a compact bright source in both bands.
- At 22 µm it sits on a north–south ridge of diffuse emission that extends ~3′ south.
- The sky is brighter near the star than farther out: near/far flux ratio 0.88 at W3 and 0.82 at W4.

This is the visual evidence for possible nebular contamination, and slide 8 quantifies it.

## Slide 8: Extended emission

**Figure.** 4′ × 4′ WISE images at 12 µm (W3) and 22 µm (W4), centred on the star, with radial profiles. A radial profile is the median brightness at each distance from the star, as a fraction of the peak.

**What it means.**
- WISE's beam at 22 µm is about 12″ wide.
- A point source (the star plus a compact disk) would fade to a few percent of peak within ~30″.
- At W4, the light stays at **15–20% of peak out to 100″**, so there is real extended dust emission around the star.
- The 10′ images on slide 7 show the star on a north–south ridge of 22 µm emission. The sky near the star is brighter than the sky farther out in every direction (near/far ratio 0.82 at W4).

This is the reason the excess may be contaminated. Part or all of the W4 flux in the photometry aperture could be this nebular dust rather than a disk.

## Slide 9: Neighbours

These checks rule out other stars as the source of the excess or of the eclipse dilution. They were run on 26 October from Kyle's PC against IRSA, SkyView, MAST and the Gaia archive.

- **≤ 4.1** is the highest W4 signal-to-noise of any other AllWISE source within 3′.
  - All are at the detection limit, and the nearest is 23–34″ away.
  - None is bright enough at 22 µm to supply the excess.
- **0** is the number of Gaia DR3 sources within 12″ (one WISE beam). Nothing hides inside the beam. The nearest TIC star is 12.8″ away at Tmag 14.5, about 6 magnitudes fainter.
- **1.8%** is the combined light of the 52 Gaia sources 12–60″ away, relative to J2312.
  - TESS pixels are large (21″), so neighbours dilute the eclipses.
  - About 2% dilution is small and does not change the fit.
  - TIC gives no official contamination ratio for this star, because it is outside the Candidate Target List.

## Slide 10: Verdict and next steps

**Excess real, origin uncertain.**
- Verified:
  - the binary model;
  - the 6.85 d period;
  - that no neighbouring star produces the 22 µm light.
- Not proven: that the excess comes from a disk. Hot stars embedded in nebular dust can heat it and produce the same signal.

**Next tests:**
1. **Excess vs aperture size.**
   - Measure the W4 excess in smaller and larger apertures.
   - Disk light is compact, so the excess should stay fixed.
   - Nebular light grows with aperture.
2. **Higher-resolution infrared images.** Spitzer MIPS 24 µm (about 6″ resolution) or Herschel, if they cover the field, would separate a point source from the extended emission.
3. **Per-sector light-curve fit.** This tightens the binary parameters once the raw TESS light curves can be pulled from MAST.
