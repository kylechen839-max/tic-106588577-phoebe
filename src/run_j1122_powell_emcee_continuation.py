#!/usr/bin/env python3
"""Continue J1122 from the detached baseline with bounded Powell and emcee."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import corner
import emcee
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import run_j1122_detached_search as base


OUTPUT_DIR = base.REPO_DIR / "outputs" / "j1122_powell_emcee_continuation"
PARAMS = ["incl", "t0", "r1", "r2", "teffratio"]
BOUNDS = {
    "incl": (85.8, 88.5),
    "t0": (base.T0 - 0.0018, base.T0 + 0.0018),
    "r1": (0.55, 0.9),
    "r2": (0.32, 0.62),
    "teffratio": (0.58, 0.95),
}
START = {
    "incl": 87.0,
    "t0": base.T0,
    "sma": 20.0,
    "r1": 0.7,
    "r2": 0.455,
    "teffratio": 0.75,
    "gravb_primary": 0.8,
    "gravb_secondary": 0.8,
}


def vector_to_candidate(vector):
    candidate = dict(START)
    candidate.update({name: float(value) for name, value in zip(PARAMS, vector)})
    return candidate


def candidate_to_vector(candidate):
    return np.array([candidate[name] for name in PARAMS], dtype=float)


def in_bounds(vector):
    return all(BOUNDS[name][0] <= value <= BOUNDS[name][1] for name, value in zip(PARAMS, vector))


def build_work_bundle(ntriangles):
    times, fluxes, sigmas = base.load_light_curve(base.INPUT_LC)
    b = base.build_bundle(times, fluxes, sigmas)
    b.set_value_all("ntriangles", ntriangles)
    return b


def evaluate_vector(b, vector, model, keep_model=False):
    if not in_bounds(vector):
        return None
    candidate = vector_to_candidate(vector)
    try:
        result = base.evaluate(b, candidate, model, keep_model=keep_model)
    except Exception as exc:
        print(f"{model} failed: {exc}", flush=True)
        return None
    return result


def run_powell():
    ntriangles = int(os.environ.get("J1122_POWELL_NTRIANGLES", "500"))
    maxiter = int(os.environ.get("J1122_POWELL_MAXITER", "12"))
    b = build_work_bundle(ntriangles)
    history = []
    counter = {"n": 0}

    def objective(vector):
        index = counter["n"]
        counter["n"] += 1
        result = evaluate_vector(b, vector, f"powell_eval_{index:04d}")
        if result is None:
            return 1e12
        history.append(result)
        # Optimize chi-square, but keep the trend metrics in the saved history.
        value = result["chi2"]
        if index % 8 == 0:
            print(
                f"powell {index}: chi2={result['chi2']:.3f} rms={result['rms']:.6f} "
                f"candidate={result['candidate']}",
                flush=True,
            )
        return value

    start = candidate_to_vector(START)
    result = minimize(
        objective,
        start,
        method="Powell",
        bounds=[BOUNDS[name] for name in PARAMS],
        options={
            "maxiter": maxiter,
            "xtol": float(os.environ.get("J1122_POWELL_XTOL", "2e-4")),
            "ftol": float(os.environ.get("J1122_POWELL_FTOL", "2e-4")),
            "disp": True,
        },
    )
    best_candidate = vector_to_candidate(result.x)
    best = base.evaluate(b, best_candidate, "j1122_after_bounded_powell", keep_model=True)
    base.save_plot(b, best, OUTPUT_DIR, "j1122_after_bounded_powell")
    b.save(str(OUTPUT_DIR / "j1122_after_bounded_powell.phoebe"))
    return {
        "optimizer_result": {
            "success": bool(result.success),
            "message": str(result.message),
            "fun": float(result.fun),
            "nit": int(result.nit),
            "nfev": int(result.nfev),
            "x": result.x.tolist(),
        },
        "best": best,
        "history": history,
    }


def log_prior(vector):
    return 0.0 if in_bounds(vector) else -np.inf


def run_emcee(powell_candidate):
    ntriangles = int(os.environ.get("J1122_EMCEE_NTRIANGLES", "400"))
    nwalkers = int(os.environ.get("J1122_EMCEE_NWALKERS", "16"))
    niters = int(os.environ.get("J1122_EMCEE_NITERS", "40"))
    discard = int(os.environ.get("J1122_EMCEE_DISCARD", str(max(20, niters // 3))))
    thin = int(os.environ.get("J1122_EMCEE_THIN", "2"))
    b = build_work_bundle(ntriangles)
    start = candidate_to_vector(powell_candidate)
    scales = np.array([0.08, 0.00025, 0.025, 0.025, 0.025], dtype=float)
    rng = np.random.default_rng(base.RNG_SEED + 22)
    pos = start + scales * rng.normal(size=(nwalkers, len(PARAMS)))
    for i in range(nwalkers):
        tries = 0
        while not in_bounds(pos[i]) and tries < 200:
            pos[i] = start + scales * rng.normal(size=len(PARAMS))
            tries += 1
        if not in_bounds(pos[i]):
            pos[i] = start

    counter = {"n": 0}
    cache = {}

    def log_probability(vector):
        lp = log_prior(vector)
        if not np.isfinite(lp):
            return -np.inf
        index = counter["n"]
        counter["n"] += 1
        result = evaluate_vector(b, vector, f"emcee_eval_{index:05d}")
        if result is None:
            return -np.inf
        cache[index] = {"vector": vector.tolist(), "result": result}
        return lp - 0.5 * result["chi2"]

    sampler = emcee.EnsembleSampler(nwalkers, len(PARAMS), log_probability)
    sampler.run_mcmc(pos, niters, progress=True)

    chain = sampler.get_chain()
    log_prob = sampler.get_log_prob()
    finite = np.isfinite(log_prob)
    if not np.any(finite):
        raise RuntimeError("No finite emcee samples")
    best_flat = int(np.nanargmax(np.where(finite, log_prob, -np.inf)))
    best_iter, best_walker = np.unravel_index(best_flat, log_prob.shape)
    best_vector = chain[best_iter, best_walker]
    best_candidate = vector_to_candidate(best_vector)

    b_hi = build_work_bundle(int(os.environ.get("J1122_FINAL_NTRIANGLES", "1500")))
    best = base.evaluate(b_hi, best_candidate, "j1122_emcee_best", keep_model=True)
    base.save_plot(b_hi, best, OUTPUT_DIR, "j1122_emcee_best")
    b_hi.save(str(OUTPUT_DIR / "j1122_emcee_best.phoebe"))

    flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
    flat_log_prob = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
    np.save(OUTPUT_DIR / "j1122_emcee_chain.npy", chain)
    np.save(OUTPUT_DIR / "j1122_emcee_log_prob.npy", log_prob)

    fig, axes = plt.subplots(len(PARAMS), 1, figsize=(10, 8), sharex=True)
    for i, name in enumerate(PARAMS):
        axes[i].plot(chain[:, :, i], alpha=0.55, lw=0.8)
        axes[i].set_ylabel(name)
    axes[-1].set_xlabel("iteration")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "j1122_emcee_trace.png", dpi=180)
    plt.close(fig)

    if flat_samples.shape[0] > 0:
        fig = corner.corner(flat_samples, labels=PARAMS, show_titles=True)
        fig.savefig(OUTPUT_DIR / "j1122_emcee_corner.png", dpi=180)
        plt.close(fig)

    summary = {
        "nwalkers": nwalkers,
        "niters": niters,
        "discard": discard,
        "thin": thin,
        "best_iter": int(best_iter),
        "best_walker": int(best_walker),
        "best_log_prob": float(log_prob[best_iter, best_walker]),
        "best": best,
        "acceptance_fraction_mean": float(np.mean(sampler.acceptance_fraction)),
        "acceptance_fraction_min": float(np.min(sampler.acceptance_fraction)),
        "acceptance_fraction_max": float(np.max(sampler.acceptance_fraction)),
        "posterior_median": dict(zip(PARAMS, np.nanmedian(flat_samples, axis=0).tolist()))
        if flat_samples.shape[0]
        else None,
        "posterior_std": dict(zip(PARAMS, np.nanstd(flat_samples, axis=0).tolist()))
        if flat_samples.shape[0]
        else None,
        "max_flat_log_prob": float(np.nanmax(flat_log_prob)) if flat_log_prob.size else None,
    }
    return summary


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": base.TARGET,
        "input_lc": str(base.INPUT_LC.relative_to(base.REPO_DIR)),
        "params": PARAMS,
        "bounds": BOUNDS,
        "start": START,
        "status": "started",
    }
    try:
        powell = run_powell()
        result["powell"] = {k: v for k, v in powell.items() if k != "history"}
        with (OUTPUT_DIR / "j1122_powell_history.json").open("w", encoding="utf-8") as handle:
            json.dump(powell["history"], handle, indent=2)
            handle.write("\n")
        result["emcee"] = run_emcee(powell["best"]["candidate"])
        result["status"] = "success"
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = repr(exc)
        raise
    finally:
        with (OUTPUT_DIR / "j1122_powell_emcee_result.json").open("w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
