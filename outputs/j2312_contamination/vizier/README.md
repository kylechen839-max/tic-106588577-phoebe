# J2312: what VizieR has (searched 2026-11-17, Harvard mirror)

All-catalogue search within 30″ returned 117 tables (`j2312_vizier_all_30arcsec.txt`). The IR and stellar rows are in
`j2312_vizier_detail.txt`. Scripts: `j2312_vizier_ir.py`, `j2312_vizier_detail.py`.

## Infrared beyond WISE
| Catalogue | Result |
|---|---|
| AKARI IRC (9/18 µm), AKARI FIS (65–160 µm) | no source within 60″ / 90″ |
| MSX, Herschel PACS point sources | none within 60″ |
| Spitzer | none (checked 11-14) |
| **IRAS PSC: IRAS 23097+5304**, 17″ away, error ellipse 47″×17″ | 12 µm < 0.25 Jy and 25 µm < 0.50 Jy (upper limits only); **60 µm 0.58 Jy, 100 µm 1.62 Jy** (good quality) |

The IRAS 60/100 µm source is not the star's dust:
- Its 60 and 100 µm shapes correlate poorly with a point source (CC_60 = C, CC_100 = D), and the 100 µm cirrus flag is high (Cirr3 = 28).
- A ~150 K dust component scaled to the 7.5 mJy W4 excess predicts only ~7 mJy at 60 µm, about 80× less than IRAS measured.
- So IRAS is seeing cold, extended nebula (the same ridge) in its arcminute beam.
- The 12/25 µm upper limits are about 50× above the WISE fluxes, so IRAS gives no useful constraint on J2312 itself.

**Conclusion: no archival IR data sharper or longer-wavelength than WISE can test disk vs envelope.**

## Stellar data (new)
- HD 236067 (BD+52 3388): HD type **A2**, Tycho/ASCC/TIC entries at the same position.
- **Gaia DR3 1995545679513855744**: G = 8.73, parallax 1.327 ± 0.015 mas (~**754 pc**), RUWE 0.89, flagged VARIABLE. No RVS radial velocity.
- **Gaia Apsis/FLAME** parameters. These assume a single star; light from the binary companion biases them.
  - GSP-Phot Teff 9930 K, log g 3.43, R 5.8 R☉, A_G 0.58.
  - ESP-HS Teff 9910 ± 76 K, log g 3.60, v sin i 57 km/s.
  - L 308 L☉, mass 3.6 M☉, **age 0.22 Gyr (0.20–0.24)**.
  - Hα equivalent width +0.38 Å, i.e. absorption: no sign of emission.
- No LAMOST spectrum or published MK spectral type turned up. Gaia BP/RP (XP) low-resolution spectra exist.

## What follows (derived, not measured)
- At 754 pc, JWST/MIRI's ~0.7″ beam at 21 µm is ~530 AU. A disk at ≲ 300 AU would be unresolved; a bow shock or dust wave at ≳ 1,000 AU (≥ 1.3″) would be resolved.
- For blackbody grains around L = 308 L☉ (T ≈ 278 K · L^¼ · d^-½), 150–200 K dust sits at ~35–60 AU. Small, superheated ISM grains could be that warm much farther out.
- The tangential velocity from the Gaia proper motion is ~16 km/s, heliocentric and not corrected to the local ISM.
