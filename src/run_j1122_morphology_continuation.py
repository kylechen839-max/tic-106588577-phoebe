#!/usr/bin/env python3
"""Target eclipse-shape residuals in the J1122 PHOEBE detached model."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import numpy as np
from scipy.optimize import minimize

import run_j1122_detached_search as base


OUTPUT_DIR = base.REPO_DIR / "outputs" / "j1122_morphology_continuation"

START = {
    "incl": 87.00638549741235,
    "t0": 1573.0875692510388,
    "sma": 20.0,
    "r1": 0.6989251602154395,
    "r2": 0.4550842850064301,
    "teffratio": 0.6979549275328235,
    "gravb_primary": 0.8,
    "gravb_secondary": 0.8,
    "ld_primary": 0.5,
    "ld_secondary": 0.5,
    "ecc": 0.0,
    "per0": 0.0,
}

CIRCULAR_PARAMS = [
    "incl",
    "t0",
    "sma",
    "r1",
    "r2",
    "teffratio",
    "ld_primary",
    "ld_secondary",
    "gravb_primary",
    "gravb_secondary",
]

ECCENTRIC_PARAMS = CIRCULAR_PARAMS + ["ecc", "per0"]

BOUNDS = {
    "incl": (85.8, 88.8),
    "t0": (base.T0 - 0.0025, base.T0 + 0.0025),
    "sma": (14.0, 30.0),
    "r1": (0.45, 1.05),
    "r2": (0.28, 0.8),
    "teffratio": (0.5, 0.95),
    "ld_primary": (0.05, 0.95),
    "ld_secondary": (0.05, 0.95),
    "gravb_primary": (0.4, 1.0),
    "gravb_secondary": (0.4, 1.0),
    "ecc": (0.0, 0.08),
    "per0": (0.0, 360.0),
}


def load_start():
    path = OUTPUT_DIR.parent / "j1122_powell_emcee_continuation" / "j1122_after_bounded_powell_hires_result.json"
    if not path.exists():
        return dict(START)
    candidate = json.loads(path.read_text())["candidate"]
    out = dict(START)
    out.update(candidate)
    return out


def build_work_bundle(ntriangles):
    times, fluxes, sigmas = base.load_light_curve(base.INPUT_LC)
    b = base.build_bundle(times, fluxes, sigmas)
    b.set_value_all("ntriangles", ntriangles)
    return b


def vector_to_candidate(vector, params, template):
    candidate = dict(template)
    for name, value in zip(params, vector):
        candidate[name] = float(value)
    if "ecc" not in params:
        candidate["ecc"] = 0.0
        candidate["per0"] = 0.0
    return candidate


def candidate_to_vector(candidate, params):
    return np.array([candidate[name] for name in params], dtype=float)


def in_bounds(vector, params):
    for name, value in zip(params, vector):
        lo, hi = BOUNDS[name]
        if value < lo or value > hi:
            return False
    candidate = vector_to_candidate(vector, params, START)
    if candidate["r1"] + candidate["r2"] > 0.16 * candidate["sma"]:
        return False
    return True


def set_candidate_extended(b, candidate):
    base.set_candidate(b, candidate)
    b.set_value("ld_coeffs@primary@lc01", [candidate["ld_primary"]])
    b.set_value("ld_coeffs@secondary@lc01", [candidate["ld_secondary"]])
    b.set_value("ecc@binary", candidate["ecc"])
    b.set_value("per0@binary", candidate["per0"])


def eclipse_shape_penalty(result):
    return (
        3.5 * result["primary_mean_abs"]
        + 2.0 * abs(result["primary_mean"])
        + 2.0 * result["secondary_mean_abs"]
        + 1.5 * abs(result["secondary_mean"])
        + 1.2 * result["shoulder_mean_abs"]
        + 0.7 * result["out_of_eclipse_mean_abs"]
    )


def evaluate(b, candidate, model, keep_model=False):
    set_candidate_extended(b, candidate)
    b.run_compute(compute="phoebe01", model=model, overwrite=True, progressbar=False)
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"), dtype=float)
    metrics = base.residual_metrics(times, observed, predicted, candidate["t0"])
    chi2 = float(b.calculate_chi2(model=model, dataset="lc01"))
    result = {
        "model": model,
        "candidate": dict(candidate),
        "chi2": chi2,
        "rms": metrics["rms"],
        "mean_abs": metrics["mean_abs"],
        "max_abs": metrics["max_abs"],
        **metrics,
    }
    result["eclipse_shape_penalty"] = float(eclipse_shape_penalty(result))
    result["objective"] = float(result["rms"] + 0.85 * result["eclipse_shape_penalty"])
    if not keep_model:
        try:
            b.remove_model(model)
        except Exception:
            pass
    return result


def optimize_stage(stage_name, start_candidate, params, ntriangles, maxiter):
    b = build_work_bundle(ntriangles)
    history = []
    counter = {"n": 0}
    maxfev = int(os.environ.get("J1122_MORPH_MAXFEV", "170"))

    def objective(vector):
        index = counter["n"]
        counter["n"] += 1
        if not in_bounds(vector, params):
            return 1e6
        candidate = vector_to_candidate(vector, params, start_candidate)
        try:
            result = evaluate(b, candidate, f"{stage_name}_eval_{index:04d}")
        except Exception as exc:
            print(f"{stage_name} {index} failed: {exc}", flush=True)
            return 1e6
        history.append(result)
        if index % 10 == 0:
            print(
                f"{stage_name} {index}: obj={result['objective']:.6f} "
                f"chi2={result['chi2']:.2f} rms={result['rms']:.6f} "
                f"primary_abs={result['primary_mean_abs']:.6f} "
                f"candidate={result['candidate']}",
                flush=True,
            )
        return result["objective"]

    x0 = candidate_to_vector(start_candidate, params)
    result = minimize(
        objective,
        x0,
        method="Powell",
        bounds=[BOUNDS[name] for name in params],
        options={
            "maxiter": maxiter,
            "maxfev": maxfev,
            "xtol": float(os.environ.get("J1122_MORPH_XTOL", "1e-4")),
            "ftol": float(os.environ.get("J1122_MORPH_FTOL", "1e-4")),
            "disp": True,
        },
    )
    best_candidate = vector_to_candidate(result.x, params, start_candidate)
    best = evaluate(b, best_candidate, stage_name, keep_model=True)
    base.save_plot(b, best, OUTPUT_DIR, stage_name)
    b.save(str(OUTPUT_DIR / f"{stage_name}.phoebe"))
    return {
        "optimizer_result": {
            "success": bool(result.success),
            "message": str(result.message),
            "fun": float(result.fun),
            "nit": int(result.nit),
            "nfev": int(result.nfev),
            "x": result.x.tolist(),
            "params": params,
        },
        "best": best,
        "history": history,
    }


def high_resolution_recompute(candidate, label):
    b = build_work_bundle(int(os.environ.get("J1122_MORPH_FINAL_NTRIANGLES", "1800")))
    result = evaluate(b, candidate, label, keep_model=True)
    base.save_plot(b, result, OUTPUT_DIR, label)
    b.save(str(OUTPUT_DIR / f"{label}.phoebe"))
    with (OUTPUT_DIR / f"{label}_result.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    return result


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ntriangles = int(os.environ.get("J1122_MORPH_NTRIANGLES", "650"))
    maxiter = int(os.environ.get("J1122_MORPH_MAXITER", "8"))
    start = load_start()

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": base.TARGET,
        "input_lc": str(base.INPUT_LC.relative_to(base.REPO_DIR)),
        "start": start,
        "bounds": BOUNDS,
        "status": "started",
    }

    try:
        result["start_hires"] = high_resolution_recompute(start, "j1122_morph_start_hires")
        circular = optimize_stage(
            "j1122_morph_circular_powell",
            start,
            CIRCULAR_PARAMS,
            ntriangles,
            maxiter,
        )
        result["circular_powell"] = {k: v for k, v in circular.items() if k != "history"}
        with (OUTPUT_DIR / "j1122_morph_circular_history.json").open("w", encoding="utf-8") as handle:
            json.dump(circular["history"], handle, indent=2)
            handle.write("\n")

        eccentric = optimize_stage(
            "j1122_morph_eccentric_powell",
            circular["best"]["candidate"],
            ECCENTRIC_PARAMS,
            ntriangles,
            maxiter,
        )
        result["eccentric_powell"] = {k: v for k, v in eccentric.items() if k != "history"}
        with (OUTPUT_DIR / "j1122_morph_eccentric_history.json").open("w", encoding="utf-8") as handle:
            json.dump(eccentric["history"], handle, indent=2)
            handle.write("\n")

        result["circular_hires"] = high_resolution_recompute(
            circular["best"]["candidate"],
            "j1122_morph_circular_hires",
        )
        result["eccentric_hires"] = high_resolution_recompute(
            eccentric["best"]["candidate"],
            "j1122_morph_eccentric_hires",
        )
        result["status"] = "success"
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = repr(exc)
        raise
    finally:
        with (OUTPUT_DIR / "j1122_morphology_result.json").open("w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
