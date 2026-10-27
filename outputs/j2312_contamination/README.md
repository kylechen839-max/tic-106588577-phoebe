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
  (±0.25 mag in W4) as a systematic error gives about 2.6-2.7σ for the bottom row and 3.0-3.1σ for the AllWISE row.
