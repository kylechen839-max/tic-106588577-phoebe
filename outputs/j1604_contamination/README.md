# J1604-56 (J160415.99-562627.3) contamination checks (Kyle's PC, 2026-11-14)

The same method as `outputs/j2312_contamination/`, run with `src/crossmatch/contamination_checks.py`. The 10′ W3/W4
cutouts come from the AllWISE Atlas tile 2399m561_ac51 on the IRSA AWS mirror (SkyView returned no WISE images
on 11-11/11-14). Console output is in `run_log.txt`.

**By eye** (`*_wise_10arcmin.png`): the star sits inside a bright, extended 12/22 µm nebula that is centred on
it and fills the 10′ field. This is stronger than J2312's halo. The W4 radial profile is still ~50% of the
centre at 25″.

| background under the star | W4 / photosphere (free; fixed Teff) | W4 σ (photometric) | W3 / photosphere |
|---|---|---|---|
| far sky 150–200″ | 6.18; 5.88 | 10.0 | 1.88 |
| 50–70″ (AllWISE catalogue) | 3.93; 3.73 | 8.9 | 1.08 |
| 25–40″ | 2.54; 2.41 | 7.2 | 0.94 |
| halo extrapolated in (25–70″ linear) | 1.19; 1.13 | 1.9 | 0.82 |

- Near/far-sky flux: 0.58 (W3) and 0.64 (W4). The near annulus is brighter than the far one in all four quadrants.
- **The W3 excess vanishes** with any local background.
- **The W4 excess falls from 3.9× to 1.2–2.5×.** The background models span ±0.9 mag, so with that systematic the
  significance is ~1σ.
- Caveat: the W4 profile at 25–70″ still contains the star's PSF wings, so the "halo extrapolated" row over-subtracts.
  The 25–40″ row (2.5×) is the more defensible lower value.
- AllWISE: the target is W4 7.26 (S/N 12), ph_qual AAAA, cc 0000. 8 sources within 60″ have W4 S/N ≥ 3, mostly
  knots in the nebula.
- TIC: contratio not given. Gaia DR3: no star within 5 mag inside 12″ (flux 0.1%). At 12–60″, 517 stars and
  4 within 5 mag (21% of the target's G flux), a crowded plane field.

**Verdict (inferred):** J1604's WISE excess is probably dominated by the surrounding nebula. It is weaker evidence
for a disk than J2312, despite J1604's clean light-curve fit. Settling it needs data sharper than WISE (no Spitzer
coverage exists; see `outputs/spitzer_mips24/`).
