#!/usr/bin/env python3
"""Small targeted PHOEBE grid for J1122 eclipse morphology."""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np

import run_j1122_detached_search as base
import run_j1122_morphology_continuation as morph


OUTPUT_DIR = base.REPO_DIR / "outputs" / "j1122_targeted_eclipse_grid"


def candidate(label, **updates):
    c = morph.load_start()
    c.update(
        {
            "ld_primary": 0.5,
            "ld_secondary": 0.5,
            "ecc": 0.0,
            "per0": 0.0,
        }
    )
    c.update(updates)
    return label, c


def build_candidates():
    start = morph.load_start()
    t0a = 1573.08737
    candidates = [
        candidate("start_physical", **start),
        candidate("partial_powell_ld", t0=t0a, ld_primary=0.5114, ld_secondary=0.8865, gravb_primary=0.984, gravb_secondary=0.984),
        candidate("high_secondary_ld", t0=t0a, ld_primary=0.50, ld_secondary=0.92, gravb_primary=0.8, gravb_secondary=0.8),
        candidate("high_secondary_ld_g1", t0=t0a, ld_primary=0.50, ld_secondary=0.92, gravb_primary=1.0, gravb_secondary=1.0),
        candidate("moderate_ld", t0=t0a, ld_primary=0.60, ld_secondary=0.75, gravb_primary=0.8, gravb_secondary=0.8),
        candidate("low_primary_high_secondary_ld", t0=t0a, ld_primary=0.35, ld_secondary=0.92),
        candidate("high_primary_high_secondary_ld", t0=t0a, ld_primary=0.75, ld_secondary=0.92),
        candidate("earlier_t0_high_ld", t0=t0a - 0.00025, ld_primary=0.50, ld_secondary=0.92),
        candidate("later_t0_high_ld", t0=t0a + 0.00025, ld_primary=0.50, ld_secondary=0.92),
        candidate("lower_incl_high_ld", incl=start["incl"] - 0.020, t0=t0a, ld_primary=0.50, ld_secondary=0.92),
        candidate("higher_incl_high_ld", incl=start["incl"] + 0.020, t0=t0a, ld_primary=0.50, ld_secondary=0.92),
        candidate("slightly_smaller_r2_high_ld", t0=t0a, r2=start["r2"] - 0.006, ld_primary=0.50, ld_secondary=0.92),
        candidate("slightly_larger_r2_high_ld", t0=t0a, r2=start["r2"] + 0.006, ld_primary=0.50, ld_secondary=0.92),
        candidate("cooler_secondary_high_ld", t0=t0a, teffratio=0.675, ld_primary=0.50, ld_secondary=0.92),
        candidate("warmer_secondary_high_ld", t0=t0a, teffratio=0.725, ld_primary=0.50, ld_secondary=0.92),
        candidate("ecc_001_per90", t0=t0a, ld_primary=0.50, ld_secondary=0.92, ecc=0.01, per0=90.0),
        candidate("ecc_001_per270", t0=t0a, ld_primary=0.50, ld_secondary=0.92, ecc=0.01, per0=270.0),
        candidate("ecc_002_per90", t0=t0a, ld_primary=0.50, ld_secondary=0.92, ecc=0.02, per0=90.0),
        candidate("ecc_002_per270", t0=t0a, ld_primary=0.50, ld_secondary=0.92, ecc=0.02, per0=270.0),
    ]
    return candidates


def evaluate_candidates(ntriangles):
    b = morph.build_work_bundle(ntriangles)
    results = []
    for index, (label, cand) in enumerate(build_candidates()):
        model = f"targeted_{index:03d}_{label}"
        try:
            result = morph.evaluate(b, cand, model, keep_model=True)
            result["label"] = label
            results.append(result)
            print(
                f"{index:02d} {label}: obj={result['objective']:.6f} "
                f"chi2={result['chi2']:.2f} rms={result['rms']:.6f} "
                f"primary_abs={result['primary_mean_abs']:.6f}",
                flush=True,
            )
        except Exception as exc:
            results.append({"label": label, "candidate": cand, "error": repr(exc)})
            print(f"{index:02d} {label}: failed {exc}", flush=True)
    return b, results


def valid_results(results):
    return [r for r in results if "objective" in r and np.isfinite(r["objective"])]


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    work_ntriangles = 1000
    final_ntriangles = 2200
    b, results = evaluate_candidates(work_ntriangles)
    valid = valid_results(results)
    if not valid:
        raise RuntimeError("No valid targeted candidates")

    by_objective = sorted(valid, key=lambda r: r["objective"])
    by_chi2 = sorted(valid, key=lambda r: r["chi2"])
    top_labels = []
    for result in by_objective[:4] + by_chi2[:4]:
        if result["label"] not in top_labels:
            top_labels.append(result["label"])

    final_results = []
    for label in top_labels:
        source = next(r for r in valid if r["label"] == label)
        final_label = f"j1122_targeted_final_{label}"
        b_hi = morph.build_work_bundle(final_ntriangles)
        final = morph.evaluate(b_hi, source["candidate"], final_label, keep_model=True)
        final["label"] = label
        morph.base.save_plot(b_hi, final, OUTPUT_DIR, final_label)
        b_hi.save(str(OUTPUT_DIR / f"{final_label}.phoebe"))
        final_results.append(final)
        print(
            f"final {label}: obj={final['objective']:.6f} "
            f"chi2={final['chi2']:.2f} rms={final['rms']:.6f} "
            f"primary_abs={final['primary_mean_abs']:.6f}",
            flush=True,
        )

    best_final_objective = min(final_results, key=lambda r: r["objective"])
    best_final_chi2 = min(final_results, key=lambda r: r["chi2"])
    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": base.TARGET,
        "input_lc": str(base.INPUT_LC.relative_to(base.REPO_DIR)),
        "work_ntriangles": work_ntriangles,
        "final_ntriangles": final_ntriangles,
        "all_work_results": results,
        "final_results": final_results,
        "best_final_objective": best_final_objective,
        "best_final_chi2": best_final_chi2,
        "status": "success",
    }
    with (OUTPUT_DIR / "j1122_targeted_eclipse_grid_result.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    with (OUTPUT_DIR / "j1122_targeted_eclipse_grid_all_work_results.json").open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
        handle.write("\n")
    morph.base.save_plot(b, by_objective[0], OUTPUT_DIR, "j1122_targeted_best_work")
    b.save(str(OUTPUT_DIR / "j1122_targeted_best_work.phoebe"))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
