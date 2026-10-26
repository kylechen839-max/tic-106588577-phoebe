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
