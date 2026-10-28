# J1122 (J112238.89-592027.5) light curve from Kyle's lightkurve pipeline (2026-10-28)

Sectors 10, 11, 37, 64, 90, 99, 100 from AWS TESScut cubes (10x10 px), `src/pipeline/lightkurve_pipeline_s3.py`.

**Kyle's default settings fail here (verified).** `create_threshold_mask(threshold=2, reference_pixel=(5,5))` picks a brighter neighbour 3-5 px away in every sector, because J1122 (Tmag 12.7, TIC contamination 0.58) is not above the threshold. Those curves show no eclipse (primary depth -0.11% to +0.05%). They are `*_S*_lk.txt` here and should not be used.

**Used instead:** `--aperture eclipse --period 4.7917045 --t0 3762.8427 --sky-subtract --suffix _eclap`. The aperture is the 6-8 pixels within 1.6 px of the target whose in-minus-out-of-eclipse difference has S/N > 3. The background regression is Kyle's (PCA-5 of pixels [8:10,5:10], RegressionCorrector), after subtracting a per-cadence sky level. The difference-image centroid sits 1-5" from the target in every sector, so the eclipse is on J1122.

Files:
- `J112238.89-592027.5_S<nn>_lk_eclap.txt`: per sector (time, flux e/s, err, flux/median).
- `J112238.89-592027.5_eclap_LightCurve.txt`, `_Periodogram.txt`, `_Folded_LC.txt`, `_Folded_LC_binned.txt`, plus plots: all 7 sectors, each divided by its median.
- `J112238.89-592027.5_lk_eclap_7sector_raw.txt`: input to the per-sector PHOEBE fit. Each sector is divided by its out-of-eclipse median, then the deviations are multiplied by 0.3705, so the mean primary depth equals the QLP 7-sector binned depth (4.28%).

**Caveat on depths.** The sky subtraction over-subtracts by a different amount in each sector, so the raw per-sector primary depths (7.0-18.1%) are inflated and not comparable with each other. Use the shape and timing; the per-sector fitter fits each sector's dilution on its own. Point-to-point rms: S10 7.1, S11 5.5, S37 10.1, S64 44.8, S90 10.1, S99 16.7, S100 23.2 ppt. S64 is the noisiest.
