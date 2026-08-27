#!/usr/bin/env python3
"""Refine J1122 PHOEBE parameters on the all-three-sector consensus light curve."""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np

import run_j1122_detached_search as base
import run_j1122_morphology_continuation as morph


OUTPUT_DIR = base.REPO_DIR / "outputs" / "j1122_min3_refinement"
MIN3_LC = base.REPO_DIR / "outputs" / "j1122_robust_binning_phoebe" / "sector_consensus_min3_lightcurve.txt"


def load_rows():
    return np.genfromtxt(MIN3_LC, comments="#")


def build_bundle(rows, ntriangles):
    b = base.build_bundle(rows[:, 0], rows[:, 1], rows[:, 2])
    b.set_value_all("ntriangles", ntriangles)
    return b


def start_candidate():
    c = morph.load_start()
    c.update({"ld_primary": 0.5, "ld_secondary": 0.5, "ecc": 0.0, "per0": 0.0})
    return c


def candidates():
    start = start_candidate()
    trial_list = []
    index = 0
    for dt0 in [-0.00045, -0.00025, 0.0, 0.00025, 0.00045]:
        for dincl in [-0.035, 0.0, 0.035]:
            for dr2 in [-0.012, -0.006, 0.0, 0.006, 0.012]:
                for teffratio in [0.675, start["teffratio"], 0.725]:
                    c = dict(start)
                    c.update(
                        {
                            "t0": start["t0"] + dt0,
                            "incl": start["incl"] + dincl,
                            "r2": start["r2"] + dr2,
                            "teffratio": teffratio,
                        }
                    )
                    trial_list.append((f"min3_refine_{index:04d}", c))
                    index += 1
    return trial_list


def evaluate_grid(rows):
    b = build_bundle(rows, 1200)
    results = []
    for index, (label, candidate) in enumerate(candidates()):
        try:
            result = morph.evaluate(b, candidate, label, keep_model=False)
            result["label"] = label
            results.append(result)
            if index % 15 == 0:
                print(
                    f"{label}: obj={result['objective']:.6f} chi2={result['chi2']:.2f} "
                    f"rms={result['rms']:.6f} primary_abs={result['primary_mean_abs']:.6f}",
                    flush=True,
                )
        except Exception as exc:
            results.append({"label": label, "candidate": candidate, "error": repr(exc)})
    return results


def valid(results):
    return [r for r in results if "objective" in r and np.isfinite(r["objective"])]


def recompute_best(rows, results):
    selected = []
    for key in ["objective", "chi2", "primary_mean_abs", "trend_score"]:
        for result in sorted(valid(results), key=lambda r: r[key])[:3]:
            if result["label"] not in {r["label"] for r in selected}:
                selected.append(result)
    final = []
    for result in selected:
        label = f"final_{result['label']}"
        b = build_bundle(rows, 2400)
        recomputed = morph.evaluate(b, result["candidate"], label, keep_model=True)
        recomputed["label"] = result["label"]
        base.save_plot(b, recomputed, OUTPUT_DIR, label)
        b.save(str(OUTPUT_DIR / f"{label}.phoebe"))
        final.append(recomputed)
        print(
            f"{label}: obj={recomputed['objective']:.6f} chi2={recomputed['chi2']:.2f} "
            f"rms={recomputed['rms']:.6f} primary_abs={recomputed['primary_mean_abs']:.6f}",
            flush=True,
        )
    return final


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    results = evaluate_grid(rows)
    final = recompute_best(rows, results)
    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": base.TARGET,
        "input_lc": str(MIN3_LC.relative_to(base.REPO_DIR)),
        "n_trials": len(results),
        "all_results": results,
        "final_results": final,
        "best_final_objective": min(final, key=lambda r: r["objective"]),
        "best_final_chi2": min(final, key=lambda r: r["chi2"]),
        "best_final_primary_abs": min(final, key=lambda r: r["primary_mean_abs"]),
        "status": "success",
    }
    with (OUTPUT_DIR / "j1122_min3_refinement_result.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    with (OUTPUT_DIR / "j1122_min3_refinement_all_results.json").open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
