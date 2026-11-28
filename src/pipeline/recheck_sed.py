"""Re-run the post-fit binary-blackbody SED check on finished PHOEBE results.

Use when a fit ran without AllWISE photometry (blackbody_check status "no AllWISE row", e.g. in a cloud
container where VizieR is blocked and outputs/crossmatch/dd1_support_cache/ is absent). Once the
AllWISE cache is available, this recomputes blackbody_check, distance_scaling and passes_all in place.

Usage: .venv-phoebe/bin/python src/pipeline/recheck_sed.py [result.json ...]   (default: all results)
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_phoebe_candidate as r  # noqa: E402

paths = [Path(p) for p in sys.argv[1:]] or sorted((r.ROOT / "outputs" / "candidates").glob("*/*_phoebe_result*.json"))
for p in paths:
    res = json.loads(p.read_text())
    name = res["designation"]
    rows = r.read_rows(r.SCORED, name) or (r.read_rows(r.EXTRA_SCORED, name) if r.EXTRA_SCORED.exists() else [])
    row = rows[0] if rows else {}
    plx = r.fnum(row.get("plx"))
    gaia_dist = r.fnum(row.get("gaia_dist")) or (1000.0 / plx if plx and plx > 0 else np.nan)
    try:
        bb = r.blackbody_check(name, res["R1_rsun"], res["R2_rsun"], res["teff1"], res["teff2"], gaia_dist)
    except FileNotFoundError as e:
        print(f"{name}: AllWISE cache missing ({e.filename})")
        continue
    if "status" in bb:
        print(f"{name}: {bb['status']}")
        continue
    res["blackbody_check"] = {k: v for k, v in bb.items() if k not in ("lam", "flux", "sigma", "model")}
    ratio = bb.get("distance_ratio_phot_over_gaia")
    if ratio and np.isfinite(ratio):
        mtot = res["mass1_assumed"] * (1 + res.get("q", res["q_assumed"]))
        res["distance_scaling"] = {
            "luminosity_factor_vs_assumed": 1 / ratio**2,
            "implied_sma_rsun": res["sma_rsun"] / ratio,
            "implied_total_mass_msun": mtot / ratio**3,
            "assumed_total_mass_msun": mtot,
        }
    phot = res["photometry_check"]
    res["passes_all"] = bool(phot["pass_fit"] and phot["pass_no_eclipse_bias"] and bb.get("pass_photosphere")
                             and bb.get("pass_distance") is not False and bb.get("pass_excess_survives"))
    p.write_text(json.dumps(res, indent=2, default=float) + "\n")
    print(f"{name}: chi2_JHKW1={bb['binary_photosphere_chi2_JHKW1']:.1f} d_ratio={ratio:.2f} passes_all={res['passes_all']}")
