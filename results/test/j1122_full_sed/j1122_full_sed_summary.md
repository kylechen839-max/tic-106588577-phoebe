# J1122 Full SED and WISE Quality Check

Target: `J112238.89-592027.5` at RA=170.6620417, Dec=-59.3409722.

## WISE Image/Catalog Quality

Nearest AllWISE source: `J112238.89-592027.5` at separation 0.093 arcsec.

| Field | Value |
| --- | --- |
| `qph` / photometric quality | `AAAA` |
| `ccf` / contamination-confusion flags | `0000` |
| `ex` / extended-source flag | `0` |
| `snr1`, `snr2`, `snr3`, `snr4` | 48.5, 51.1, 34.7, 13.0 |
| `chi2W1`, `chi2W2`, `chi2W3`, `chi2W4` | 2.84, 1.59, 0.86, 1.06 |
| `nb`, `na` | 1, 0 |
| `sat1`-`sat4` | 0.0, 0.0, 0.0, 0.0 |
| `fmoon` | 0 |

Interpretation: the AllWISE catalog flags are favorable for using the W3/W4 points (`qph=AAAA`, `ccf=0000`, `ex=0`, no saturation). W4 is still lower-resolution and the local cutout has structured background, so the W4 excess should be treated as credible but not immune from follow-up checks.

## Nearby Sources

The matched 2MASS source is essentially coincident, but there are fainter 2MASS/Gaia neighbors inside the WISE W3/W4 beam. They are much fainter in near-IR/optical, so they are unlikely to explain the full W3/W4 excess alone, but they remain a blending caveat.

Closest 2MASS rows:

| sep arcsec | 2MASS | J | H | K | Qflg |
| ---: | --- | ---: | ---: | ---: | --- |
| 0.037 | 11223889-5920275 | 12.109000205993652 | 11.901000022888184 | 11.788999557495117 | AAA |
| 5.591 | 11223943-5920313 | 15.508000373840332 | 15.307000160217285 | 13.652999877929688 | CDU |
| 8.645 | 11223917-5920358 | 15.470999717712402 | 15.177000045776367 | 14.300000190734863 | ABU |
| 11.767 | 11223780-5920192 | 16.42099952697754 | 16.024999618530273 | 16.527000427246094 | BCU |
| 14.340 | 11223889-5920418 | 12.918999671936035 | 12.170000076293945 | 11.95300006866455 | AAA |

Closest Gaia DR3 rows:

| sep arcsec | source_id | G | BP | RP | parallax mas | RUWE |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 0.029 | 5339178269248016256 | 13.047956 | 13.277884 | 12.660706 | 0.496 | 0.991 |
| 1.878 | 5339178269205724032 | 20.234625 | None | None | None | None |
| 3.916 | 5339178269248012672 | 20.060041 | 19.712313 | 18.763662 | -0.1149 | 1.051 |
| 5.083 | 5339178269248013824 | 19.217934 | None | None | -0.121 | 0.968 |
| 5.910 | 5339178269205731200 | 17.54265 | 17.889229 | 16.610361 | 0.0861 | 1.055 |
| 6.952 | 5339178269205740160 | 20.250408 | 20.70131 | 19.50531 | 0.1868 | 1.066 |
| 7.247 | 5339178269207555584 | 20.77823 | 20.815016 | 19.788513 | 3.1916 | 1.029 |
| 7.938 | 5339178269205721856 | 18.913815 | 19.756845 | 17.944519 | 0.5445 | 0.957 |

## SED Fit

The stellar baseline is a single blackbody fit to `J`, `H`, `Ks`, and `W1`; W2/W3/W4 are then tested as excess bands. This is intentionally conservative because W2 can already contain warm excess.

Best stellar blackbody temperature: **4685 K**.

| Band | observed Jy | stellar Jy | obs/star | excess sigma |
| --- | ---: | ---: | ---: | ---: |
| J | 0.0228501 | 0.0217107 | 1.05 | 2.1 |
| H | 0.0177787 | 0.0183652 | 0.97 | -1.4 |
| Ks | 0.0128331 | 0.0142305 | 0.90 | -4.1 |
| W1 | 0.00853329 | 0.00798656 | 1.07 | 3.2 |
| W2 | 0.00535292 | 0.00487634 | 1.10 | 4.6 |
| W3 | 0.00757012 | 0.000958818 | 7.90 | 30.6 |
| W4 | 0.00851848 | 0.000280341 | 30.39 | 12.5 |

## Excess Models

| Model | Key parameter | chi2 W2/W3/W4 | chi2 all bands |
| --- | ---: | ---: | ---: |
| Stellar only | 4685 K | 1113.05 | 1145.90 |
| Stellar + dust blackbody | Tdust=316 K | 7.29 | 39.50 |
| Stellar + free-free power law | alpha=-0.10 | 710.30 | 868.26 |

## Finding

The W3/W4 excess remains highly significant after fitting the near-IR photosphere. A single-temperature dust blackbody near **316 K** explains the W2-W4 excess much better than an ordinary free-free spectrum constrained to alpha=[-0.1, 0.6]. On these data, free-free emission does **not** plausibly account for all of the mid-IR excess by itself.

The most plausible interpretation from this SED-only check is a warm circumstellar/circumbinary dust component, not a pure free-free continuum. The main caveat is WISE beam blending/background structure, especially in W4; the next high-value check would be higher-resolution mid-IR imaging or confirming the excess with independent W3/W4-quality photometry.

Generated files:

- `j1122_sed_fit.png`
- `j1122_wise_cutouts.png`
- `j1122_photometry.csv`
- `j1122_catalog_quality.json`
