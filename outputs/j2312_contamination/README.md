# J2312 contamination checks (run 2026-10-26 on Kyle's PC)

Target: J231201.40+532028.6 = TIC 252741805 = Gaia DR3 1995545679513855744, G = 8.73,
RA 348.00585, Dec +53.34130 (from `outputs/crossmatch/dd1_candidates_scored.csv`; the
RA 347.5058 quoted in the hand-off note is a typo, 23h12m01.40s = 348.0058 deg).

No earlier J2312 results existed anywhere under `C:\Users\kylec` before this run.

`j2312_checks.py` reimplements `src/crossmatch/w4_neighbours.py` (commit 081ab9c) with the
same method and parameters (SkyView 10' at 300 px = 2"/px, 16.5" aperture, near sky
50-70", far sky 150-200", AllWISE 3' cone), and adds the TIC and Gaia queries. It was
run instead of the repo script only because this PC's permission policy would not run
code from the fresh clone. `zoom.py` makes the 4' linear cutouts and radial profile.
Python 3.12, astroquery from PyPI, run live against IRSA, SkyView, MAST and the Gaia archive.

## Eclipsing-binary catalogues

Live VizieR cone search, radius 5", on 2026-10-26 (`eb_catalogues.py`, output in
`eb_catalogues_log.txt`). J2312 is listed as an eclipsing binary in:

| catalogue | VizieR | ID | entry |
|---|---|---|---|
| Gaia DR3 variability classifier | I/358/vclassre | 1995545679513855744 | Class ECL, score 0.48 (0.02") |
| Gaia DR3 eclipsing binaries | I/358/veb | 1995545679513855744 | P = 6.8494 ± 0.0003 d (freq 0.145999 c/d); two Gaussian eclipses, depths 0.29 / 0.25 mag at phases 0.239 / 0.739; ellipsoidal term on eclipse 1 |
| VSX (AAVSO) | B/vsx/vsx | PMAK V105 (OID 2223897) | EA, P = 6.84952 d, epoch HJD 2458758.51, V 8.75-0.31 mag (TESS-band depth), Sp A2 (0.07") |

Gaia DR3's variability summary (I/358/varisum) also flags it VEB = 1. The amplitude is 0.30 mag in G
from 52 epochs.

It is **not** in any of the other EB catalogues searched: ASAS-SN (II/366), TESS EBs of
Prša+2022 (J/ApJS/258/16), TESS OBA EBs of IJspeert+2021 (J/A+A/652/A120), ZTF
(J/ApJS/249/18), ATLAS (J/AJ/156/241), WISE variables (J/ApJS/237/28), CRTS
(J/ApJS/213/9), Kepler EBs (J/AJ/151/68), ASAS-3 (II/264), NSVS (J/AJ/138/466),
Avvakumova+2013 (J/AN/334/860), GCVS (B/gcvs), SB9 (B/sb9), or J/ApJS/247/56. This matches
the pipeline's `surveys = gaia_dr3;vsx`. The two catalogue periods agree to 0.0001 d (about
0.5σ). The star is at G = 8.7, which is bright enough to saturate ASAS-SN, ZTF, ATLAS
and CRTS. That probably explains why it is missing from those surveys (assumed, not checked).

## Results

**AllWISE within 3'** (`*_allwise_3arcmin.csv`): the target is W4 = 6.84, SNR 17.5,
w4rchi2 1.14, ext_flg 0, cc_flags 0000. Nearest neighbours are at 23-34" and are faint
(W1 13.4-15.4, W4 SNR <= 4.1, mostly at the W4 limit of ~8.4-9.0, with halo flags `hh`
from the target). The brightest W1 sources in the cone are J231205.47+531940.2 (60.6", W1 10.2,
W4 SNR 1.4) and J231210.14+532103.9 (85.8", W1 10.4, W4 SNR 1.4). No neighbour is a
plausible W4 contaminant at the 16.5" aperture.

**Near vs far sky** (DN, SkyView atlas images):

| band | flux (near sky) | flux (far sky) | near/far | sky near | sky far |
|---|---|---|---|---|---|
| W3 | 4539.7 | 5165.3 | 0.88 | 327.175 | 324.224 |
| W4 | 112.1 | 136.1 | 0.82 | 93.286 | 93.173 |

The W4 sky in the 50-70" annulus is above the 150-200" sky in all four quadrants
(S +0.151, E +0.131, N +0.124, W +0.052 DN), so the star sits at the centre of a
diffuse 22 um halo, not on a one-sided gradient.

**By eye** (`*_wise_10arcmin.png`, `*_wise_zoom_profile.png`): the W3 and W4 point
source is compact and clearly detected above the background. It sits on extended emission:
a N-S ridge through the star that runs ~3' to the south, inside a larger diffuse arc
seen at both 12 and 22 um. The median radial profile at W4 is still ~15-20% of the peak
level at 50-100" and ~5% at 150-200". At W3 it is ~5% and ~2%.

**TIC v8.2** (`*_tic_60arcsec.csv`): `contratio` is null for TIC 252741805 (Teff 8547 K,
R = 6.3 Rsun; the star is not in the CTL, which is where contratio is computed). The
scored CSV has it null too. Nearest TIC neighbour is 252741811, Tmag 14.46 at 12.8"
(dTmag 5.85).

**Gaia DR3 within 60"** (`*_gaiadr3_60arcsec.csv`): nothing within 12" except the target.
52 sources at 12-60". None is brighter than G = 13.7 (dG < 5). The brightest is G 14.86
at 12.9". Summed G flux at 12-60" is 1.8% of the target. This is my stand-in for the missing
TIC contratio: dilution of the TESS light curve is about 2% or less.

## Verdict

- **Neighbours / TESS dilution: not a problem (verified).** No catalogued IR or optical
  neighbour can supply the W4 flux, and optical contamination is ~2%.
- **W4 excess: possibly contaminated by extended dust emission (assumed, not proven).**
  The excess comes from a real compact source, but the star is centred on a 22 um halo plus
  a N-S ridge. Changing the sky annulus moves W4 by ~18%. That alone does not remove a
  W4/photosphere ratio of 3.6. But the annulus test cannot remove nebular emission that
  peaks at the star's position inside the 16.5" aperture. An 8500 K giant heating its
  surrounding ISM (an IR cirrus hotspot / bow-shock-like halo) would look like this, with
  no circumstellar disk needed. Higher-resolution IR imaging (Spitzer MIPS 24 if it exists,
  or Herschel) or an SED test of whether the excess flux scales with aperture size would
  settle it.

Verified here: catalogue numbers, aperture fluxes, profile and the visual inspection
above. Assumed or interpretive: the physical origin of the halo, and that Gaia flux ratio
stands in for the TIC contratio.

## Excess after subtracting the local W4/W3 background (2026-10-27, cloud)

The AllWISE catalogue fluxes subtract the sky measured 50–70″ from the star. The extended emission is brighter closer in, so I re-measured the 16.5″-aperture flux in the 10′ images with the sky taken from tighter annuli.
- Flux relative to the 50–70″ sky: W4 0.97 (35–50″), 0.84 (25–35″), 0.75 (20–30″); W3 0.96, 0.91, 0.87.
- The 20–30″ annulus still holds a little PSF wing, so it over-subtracts slightly.

Scaling the catalogue flux by these factors gives the excess above the binary photosphere (reference fit; ratio 3.56 and 12.6σ at W4, 1.24 and 6.9σ at W3). The catalogue error is kept:

| Sky annulus | W4 ratio | W4 excess | W3 ratio | W3 excess |
|---|---:|---:|---:|---:|
| 50–70″ (catalogue) | 3.56 | 12.6σ | 1.24 | 6.9σ |
| 25–35″ | 2.98 | 9.7σ | 1.12 | 3.5σ |
| 20–30″ | 2.66 | 8.2σ | 1.07 | 2.1σ |

Most of the W4 excess survives local background subtraction (still 2.7–3.0× the photosphere, 8–10σ). The W3 excess mostly does not.

## Excess after subtracting the W4 background halo (2026-10-27)

`w4_bgsub.py`, output in `w4_bgsub_log.txt` and `*_w4_background_subtraction.png`. The 16.5"
aperture flux is re-measured with four estimates of the background under the star, scaled
to the AllWISE magnitude (which uses the 50-70" sky), and divided by the same photospheres
as `dd1_candidates_scored.csv`. That is the free blackbody (10,325 K) and, for W4, the
catalogue Teff (8,547 K). Errors use the pipeline's rule: catalogue error with a 0.03 mag
floor. W4 is 0.062 mag (SNR 17.5).

| background under the star | W4 mag | W4 / photosphere (free, fixed Teff) | W4 excess σ | W3 / photosphere | W3 excess σ |
|---|---|---|---|---|---|
| far sky 150-200" (no halo removed) | 6.63 | 4.40, 4.05 | 13.5, 13.2 | 1.43 | 10.9 |
| 50-70" sky (AllWISE catalogue) | 6.84 | 3.63, 3.34 | 12.7, 12.3 | 1.26 | 7.4 |
| 25-40" sky | 6.98 | 3.18, 2.93 | 12.0, 11.5 | 1.17 | 5.3 |
| halo extrapolated to the star (linear fit 25-70") | 7.14 | 2.75, 2.53 | 11.1, 10.6 | 1.06 | 2.1 |

- **W4:** after removing the halo, W4 is fainter by up to 0.30 mag (flux × 0.76). The excess
  drops from 3.6× to about 2.5-2.8× the photosphere, but it stays above 10σ in every case.
  The W4 excess is not explained by the background halo alone.
- **W3:** the excess is small to begin with. With the most aggressive subtraction it falls to
  1.06× (2.1σ), i.e. the W3 excess may be entirely background.
- The two inner estimates probably over-subtract. At 25-40" the W4 profile still contains the
  star's PSF wings, so those rows are lower limits on the excess. What they cannot exclude is
  halo emission that rises steeply inside 16.5", which is invisible under the PSF (assumed
  not to happen).
- The excess σ shown is photometric only. Adding the spread of the background estimates
  (±0.25 mag in W4) as a systematic error gives about 2.5-2.6σ for the bottom row and 2.9-3.0σ for the AllWISE row.

## Is the compact excess the same material as the halo? (2026-11-14, cloud)

`src/checks/j2312_excess_colour.py`, output in `j2312_excess_colour_log.txt` and `j2312_excess_colour.png`. It uses the 10′ atlas images above, calibrated with the star's own AllWISE magnitudes in the 16.5″ aperture.

| Component | F_ν(12 µm)/F_ν(22 µm) | Blackbody T |
|---|---:|---:|
| Diffuse halo, rings 25–150″ above the 150–200″ sky | 0.75–0.99 | ~300 K-equivalent (not thermal; 12 µm has PAH bands) |
| Halo, pixel W3-vs-W4 regression 25–150″ | 0.55 | — |
| Compact excess, catalogue sky 50–70″ | 0.35 ± 0.05 | 200 K |
| Compact excess, halo extrapolated under the star | 0.12 ± 0.06 | 147 K |

Source widths are Gaussian FWHM in these images: W3 16.3″ (mostly photosphere, so this is the effective PSF) and W4 27.4″. The W4/W3 width ratio is 1.68, below the 22/12 µm diffraction ratio of ~1.9, so W4 is not resolved; it is offset 2.6″.

- **Verified:** the compact excess is 2–8× redder (weaker at 12 µm relative to 22 µm) than the nebula around it, for either background choice. The surrounding halo is also present on all sides and is more nearly flat in F_ν.
- **Interpretation (inferred):** if the excess were simply more of the surrounding nebula piled up under the star, it would have the nebula's colour (≈ 0.55–1.0). Instead it looks like a distinct component of ~150–200 K thermal grains. For L ≈ 230 L☉ that is ~30–60 AU (blackbody grains), unresolved at 742 pc (0.05″). Interstellar dust heated by the star at thousands of AU, the "IR cirrus hotspot" alternative, would be ≲ 30 K in equilibrium. Its 22 µm emission would come from stochastically heated small grains, which would carry the halo's PAH-bright 12 µm colour. This favours circumstellar or circumbinary dust over a heated-ISM hotspot, but does not prove a disk.
- **Caveats:**
  - the halo-subtracted W3 excess is only ~2σ (0.90 ± 0.44 mJy);
  - the halo inside 16.5″ is extrapolated, not measured;
  - atlas images are calibrated via the star (one point);
  - the width test only excludes extension larger than ~the W4 PSF (12″ ≈ 9,000 AU).
- **What would settle it:**
  - Spitzer MIPS 24 µm (6″) or IRAC, if J2312 was observed;
  - AKARI 18 µm, an independent epoch and beam;
  - Herschel 70/100 µm to see the cold side of the SED;
  - mid-IR spectroscopy showing silicate rather than PAH features.

  All of these need IRSA, from Kyle's PC.

## Higher-resolution / longer-wavelength archives (2026-11-14, IRSA cone searches)

These were IRSA SCS queries made through a web fetch, since the container can't reach IRSA directly.

| Archive (IRSA table) | Radius | Result |
|---|---|---|
| Spitzer SEIP source list (`slphotdr4`) | 3′ | 0 sources, so the field is almost certainly not covered by IRAC/MIPS |
| AKARI IRC 9/18 µm (`akari_irc`) | 36″ | 0 sources; the 18 µm PSC limit (~0.1 Jy) is ≫ the W4 flux (~15 mJy), so this is not constraining |
| AKARI FIS (`akari_fis`) | 72″ | 0 sources |
| Herschel PACS PSC 70 / 100 µm (`ppsc_70`, `ppsc_100`) | 3′ | 0 sources (coverage unknown) |
| IRAS PSC (`iraspsc`) | 72″ | 23097+5304 at 17.6″ (ellipse 47″×17″): 12/25/60 µm upper limits 0.25/0.50/0.58 Jy (fqual 1); 100 µm 1.62 Jy (fqual 3), cirr3 = 28 MJy/sr, cc_100 = D |

The IRAS 100 µm detection in a cirrus-flagged field matches the cold extended emission around the star, not a compact disk. No archival image sharper than WISE exists at 12–24 µm for J2312. Settling the disk question needs new mid-IR imaging (e.g. JWST/MIRI or a ground-based 10–20 µm camera).

## W4 point-source flux with the nebula modelled (2026-11-14, cloud)

`src/checks/j2312_psf_extraction.py` and `src/checks/j2312_w4_width.py` use the native AllWISE Atlas tiles (3476p530_ac51 W3/W4, 1.375″/px) from the IRSA AWS mirror. The tiles are 67 MB each in `native/` and are not committed. Empirical PSFs are median stacks of 40 isolated field stars per band (FWHM W4 17.0″, W3 11.3″). The star's PSF amplitude is fitted together with a background model, and every flux is scaled to the AllWISE catalogue value through the catalogue-like fit (a0), so zero points cancel. The photosphere is the binary photosphere (4.31 mJy).

| Background model fitted with the star | W4 (mJy) | × photosphere | Excess σ |
|---|---:|---:|---:|
| a0: flat sky at 50–70″ (catalogue-like) | 15.4 | 3.56 | 12.6 |
| a1/a2: free constant / plane (< 55″) | 14.5–14.6 | 3.37–3.38 | 11.6–11.7 |
| a3: quadratic (< 55″) | 12.3 | 2.84 | 9.1 |
| **B: W3 image as a nebula template, smoothed to W4 resolution and scaled, plus a plane (< 150″)** | **11.8** | **2.73** | **8.5** |
| a4: quadratic + free-width Gaussian halo at the star | — | — | degenerate, not used |

- **Method B is the most complete model, and its residual map is flat at the star** (`j2312_psf_extraction.png`). The nebula's own 12 µm shape predicts its 22 µm light under the star; after removing it, a point source of 2.7× the photosphere remains, an excess of ~7.5 mJy at 8.5σ. This agrees with the halo-extrapolation estimate above (2.75×).
- **a4 is degenerate:** a 12″-wide Gaussian can replace the 17″ PSF. On the W3 image, where J2312 is mostly photosphere, the same fit also moves ~45% of the flux into a "halo", so it cannot measure nebular light.
- ~~**Is the W4 source wider than a point source?** It depends on how much of the curved halo the fit includes. Within 20″ the Gaussian FWHM is 20.8″ (plane) or 17.1″ (quadratic); 76 field sources of similar brightness give 19.0 / 20.0 / 21.5″ (16/50/84%). So the core is point-like.~~ **Corrected 2026-11-15:** the field numbers quoted here were measured with a 30″ fit radius, not 20″, so the comparison was wrong. Like for like (20″, plane), J2312 is 20.8″ against 16.1″ for the 40-star PSF and 16.5″ (median) / 18.4″ (95%) for field sources; see the next section. Only the quadratic-background fit (17.1″) looks point-like, because the quadratic absorbs the wings.
- **Remaining limit (unchanged):** nebular emission concentrated inside ~10″ of the star and following the star's PSF cannot be separated at WISE resolution by any subtraction. Only sharper imaging (e.g. JWST/MIRI) removes it.

## Multi-scale separation and width test (2026-11-15, cloud)

Prompted by a literature check on building a 22 µm nebula map (Meisner & Finkbeiner 2014; Anderson et al. 2014; getsf, Men'shchikov 2021; Koenig & Leisawitz 2014; Patel et al. 2017).

**W4-only multi-scale separation** (`src/checks/j2312_multiscale.py`, `j2312_multiscale.png`, `j2312_multiscale_log.txt`). This is a simplified Python version of getsf's steps; getsf itself is not reachable from the cloud. Iterative median-filter background with window W, difference-of-Gaussian scales 4–114″, structures sorted into round "sources" and elongated "filaments", PSF fit on the source image. Flux recovery is calibrated on 12 isolated field stars (W4) run through the same pipeline. No W3 colour assumption is used.

| Background window W | Recovered W4 point flux | × photosphere | Excess σ (catalogue error) |
|---|---:|---:|---:|
| 60″ | 14.6 mJy | 3.39 | 11.7 |
| 90″ | 12.4 mJy | 2.87 | 9.2 |
| 120″ | 12.3 mJy | 2.86 | 9.1 |
| 180″ | 12.5 mJy | 2.89 | 9.3 |

A 60″ window is narrower than the halo, so it leaves halo in the "source" image. For W ≥ 90″ the answer is stable at 12.3–12.5 mJy, in line with the quadratic (12.3) and W3-template (11.8) fits. The W3 control (J2312 at 12 µm, mostly photosphere) gives 0.80–0.87 of the catalogue-like flux for W ≥ 90″, so ~15% of the catalogue W3 flux is nebula.

**The W4 source is broader than a point source** (`src/checks/tierA_w4_shape.py`, same tile). Gaussian + plane, 20″ fit radius:

| | W3 FWHM | W4 FWHM |
|---|---:|---:|
| 40-star empirical PSF | 13.3″ | 16.1″ |
| 13 field point sources ≥ ⅓ J2312's amplitude, median / 95% | — | 16.5″ / 18.4″ |
| **J2312** | **14.2″** | **20.8″** (above all 13) |

At 12 µm J2312 is point-like; at 22 µm it is not. Taken at face value, the extra width is a component of ≈ √(20.8² − 16.1²) ≈ 13″ FWHM, about 10,000 AU at 742 pc. That is far too large for a debris disk (blackbody estimate 30–60 AU, 0.05″), so part of the W4 light comes from dust on ~10⁴ AU scales around the star: nebula concentrated toward the star, or interstellar dust warmed by the star. All the flux numbers above put this component into the "point source", because none of the models can separate a 13″ component from a 16″ PSF (the a4 fit tried and was degenerate).

**What this changes (inferred):** the 2.7–2.9× is an upper limit on the compact excess rather than a measurement of it. The case for a circumbinary disk in J2312 is weaker than the 11-14 summary said. The same broadening is seen in 6 other tier A stars, but not in J2209 or J2218 (see `docs/DD_EB_CROSSMATCH.md`). J2209's MIPS 24 µm source is compact (`outputs/spitzer_mips24/`), which independently confirms the test.

**Width after the nebula model is removed (2026-11-15, `src/checks/j2312_model_figure.py`, `j2312_model_figure.png`):** after the method-B nebula model (W3 template + plane) is subtracted, J2312's W4 Gaussian FWHM within 20″ drops from 20.8″ to **18.8″**, compared with 16.1″ for the PSF. The field reference has a 95th percentile of 18.4″ and is measured without nebula subtraction. So about half of the extra width comes from the nebula ridge through the star, which the model removes. What remains is marginally broader than a point source (an intrinsic component of ≈ √(18.8² − 16.1²) ≈ 10″). The excess is still best read as an upper limit on a compact disk, but the evidence for an extended component is now weaker than the 20.8″ alone implied. The radial profile in the figure matches the fitted point source within the error bars out to ~30″.
