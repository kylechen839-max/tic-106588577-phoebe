#!/usr/bin/env python3
"""Fast local PHOEBE search for a better detached J1122 light-curve model."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import phoebe


REPO_DIR = Path(__file__).resolve().parents[1]
TARGET = "J112238.89-592027.5"
PERIOD = 4.790978516308605
T0 = 1573.0869
INPUT_LC = (
    REPO_DIR
    / "outputs"
    / "j1122_tess_multisector_s10_s11_s37"
    / f"{TARGET}_LightCurve_multisector_binned_220.txt"
)
OUTPUT_DIR = REPO_DIR / "outputs" / "j1122_detached_search"
RNG_SEED = 112238


def load_light_curve(path: Path):
    data = np.genfromtxt(path, comments="#").T
    times, fluxes, sigmas = data[0], data[1], data[2]
    sigmas = np.where(np.isfinite(sigmas) & (sigmas > 0), sigmas, np.nanmedian(sigmas[sigmas > 0]))
    sigma_floor = float(os.environ.get("J1122_SIGMA_FLOOR", "0.0004"))
    sigmas = np.maximum(sigmas, sigma_floor)
    return times, fluxes, sigmas


def phase_values(times, t0=T0):
    phase = ((times - t0) / PERIOD) % 1.0
    return np.where(phase > 0.8, phase - 1.0, phase)


def build_bundle(times, fluxes, sigmas):
    b = phoebe.default_binary()
    b.add_dataset(
        "lc",
        times=times,
        fluxes=fluxes,
        sigmas=sigmas,
        dataset="lc01",
        passband="TESS:T",
    )
    b.set_value_all("atm", "blackbody")
    b.set_value_all("ld_mode_bol", "manual")
    b.set_value_all("ld_func_bol", "linear")
    b.set_value_all("ld_coeffs_bol", [0.5])
    b.set_value_all("ld_mode", "manual")
    b.set_value_all("ld_func", "linear")
    b.set_value_all("ld_coeffs", [0.5])
    b.set_value("period@binary", PERIOD)
    b.set_value("t0_supconj@binary", T0)
    b.set_value("ecc@binary", 0.0)
    b.set_value("q@binary", 1.0)
    b.set_value("sma@binary", 20.0)
    b.set_value("teff@primary", 6000.0)
    b.set_value("irrad_frac_refl_bol@primary", 0.5)
    b.set_value("irrad_frac_refl_bol@secondary", 0.5)
    b.set_value("pblum_mode@lc01", "dataset-scaled")
    b.set_value("distortion_method@primary@phoebe01", "sphere")
    b.set_value("distortion_method@secondary@phoebe01", "sphere")
    try:
        b.set_value_all("ntriangles", 300)
    except Exception:
        pass
    return b


def set_candidate(b, candidate):
    b.set_value("incl@binary", candidate["incl"])
    b.set_value("t0_supconj@binary", candidate["t0"])
    b.set_value("sma@binary", candidate["sma"])
    b.set_value("requiv@primary", candidate["r1"])
    b.set_value("requiv@secondary", candidate["r2"])
    b.set_value("teff@secondary", 6000.0 * candidate["teffratio"])
    b.set_value("gravb_bol@primary", candidate["gravb_primary"])
    b.set_value("gravb_bol@secondary", candidate["gravb_secondary"])


def residual_metrics(times, observed, predicted, t0):
    phase = phase_values(times, t0)
    residuals = observed - predicted
    primary = np.abs(phase) < 0.055
    secondary = np.abs(phase - 0.5) < 0.055
    out = ((phase > 0.12) & (phase < 0.38)) | ((phase > 0.62) & (phase < 0.78)) | (
        (phase < -0.09) & (phase > -0.19)
    )
    shoulders = ((np.abs(phase) >= 0.055) & (np.abs(phase) < 0.12)) | (
        (np.abs(phase - 0.5) >= 0.055) & (np.abs(phase - 0.5) < 0.12)
    )

    def safe_mean(mask):
        return float(np.mean(residuals[mask])) if np.any(mask) else float("nan")

    def safe_mae(mask):
        return float(np.mean(np.abs(residuals[mask]))) if np.any(mask) else float("nan")

    trend_score = (
        3.0 * abs(safe_mean(primary))
        + 2.5 * abs(safe_mean(secondary))
        + 1.5 * safe_mae(primary)
        + 1.3 * safe_mae(secondary)
        + 1.0 * safe_mae(shoulders)
        + 0.8 * safe_mae(out)
    )
    return {
        "rms": float(np.sqrt(np.mean(residuals**2))),
        "mean_abs": float(np.mean(np.abs(residuals))),
        "max_abs": float(np.max(np.abs(residuals))),
        "primary_mean": safe_mean(primary),
        "secondary_mean": safe_mean(secondary),
        "primary_mean_abs": safe_mae(primary),
        "secondary_mean_abs": safe_mae(secondary),
        "out_of_eclipse_mean_abs": safe_mae(out),
        "shoulder_mean_abs": safe_mae(shoulders),
        "trend_score": float(trend_score),
    }


def evaluate(b, candidate, model, keep_model=False):
    set_candidate(b, candidate)
    b.run_compute(compute="phoebe01", model=model, overwrite=True, progressbar=False)
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"), dtype=float)
    metrics = residual_metrics(times, observed, predicted, candidate["t0"])
    objective = metrics["rms"] + 0.7 * metrics["trend_score"]
    result = {
        "model": model,
        "candidate": dict(candidate),
        "chi2": float(b.calculate_chi2(model=model, dataset="lc01")),
        "objective": float(objective),
        **metrics,
    }
    if not keep_model:
        try:
            b.remove_model(model)
        except Exception:
            pass
    return result


def candidate_defaults():
    return {
        "incl": 88.0,
        "t0": T0,
        "sma": 20.0,
        "r1": 0.55,
        "r2": 0.55,
        "teffratio": 0.75,
        "gravb_primary": 0.8,
        "gravb_secondary": 0.8,
    }


def candidate_grid():
    base = candidate_defaults()
    candidates = []
    for incl in [84.0, 85.5, 87.0, 88.2, 89.0]:
        for r1 in [0.35, 0.5, 0.7, 0.9, 1.1]:
            for ratio in [0.65, 0.78, 0.9, 1.05]:
                c = dict(base)
                c.update({"incl": incl, "r1": r1, "r2": r1 * ratio})
                candidates.append(c)

    rng = np.random.default_rng(RNG_SEED)
    for _ in range(int(os.environ.get("J1122_RANDOM_CANDIDATES", "0"))):
        r1 = float(rng.uniform(0.28, 1.35))
        c = dict(base)
        c.update(
            {
                "incl": float(rng.uniform(82.0, 89.7)),
                "t0": float(T0 + rng.normal(0.0, 0.0035)),
                "sma": float(rng.uniform(18.0, 28.0)),
                "r1": r1,
                "r2": float(np.clip(r1 * rng.uniform(0.55, 1.35), 0.22, 1.55)),
                "teffratio": float(rng.uniform(0.48, 1.05)),
                "gravb_primary": float(rng.uniform(0.5, 1.0)),
                "gravb_secondary": float(rng.uniform(0.5, 1.0)),
            }
        )
        candidates.append(c)
    return candidates


def local_candidates(center, n=25):
    rng = np.random.default_rng(RNG_SEED + 1)
    candidates = []
    for _ in range(n):
        c = dict(center)
        c.update(
            {
                "incl": float(np.clip(center["incl"] + rng.normal(0.0, 0.9), 80.0, 89.9)),
                "t0": float(center["t0"] + rng.normal(0.0, 0.0012)),
                "sma": float(np.clip(center["sma"] + rng.normal(0.0, 1.3), 16.0, 32.0)),
                "r1": float(np.clip(center["r1"] + rng.normal(0.0, 0.11), 0.18, 1.8)),
                "r2": float(np.clip(center["r2"] + rng.normal(0.0, 0.11), 0.18, 1.8)),
                "teffratio": float(np.clip(center["teffratio"] + rng.normal(0.0, 0.06), 0.35, 1.2)),
                "gravb_primary": float(np.clip(center["gravb_primary"] + rng.normal(0.0, 0.08), 0.3, 1.0)),
                "gravb_secondary": float(np.clip(center["gravb_secondary"] + rng.normal(0.0, 0.08), 0.3, 1.0)),
            }
        )
        candidates.append(c)
    return candidates


def run_candidates(b, candidates, label):
    results = []
    best = None
    for index, candidate in enumerate(candidates):
        model = f"{label}_{index:04d}"
        try:
            result = evaluate(b, candidate, model)
            results.append(result)
            if best is None or result["objective"] < best["objective"]:
                best = result
                print(
                    f"{label} best {index}: obj={result['objective']:.6f} "
                    f"rms={result['rms']:.6f} chi2={result['chi2']:.1f} "
                    f"candidate={result['candidate']}",
                    flush=True,
                )
        except Exception as exc:
            results.append({"model": model, "candidate": candidate, "error": repr(exc)})
    if best is None:
        raise RuntimeError(f"No valid candidates in {label}")
    return best, results


def best_valid(results, key):
    valid = [result for result in results if key in result and np.isfinite(result[key])]
    if not valid:
        raise RuntimeError(f"No valid results for key={key}")
    return min(valid, key=lambda result: result[key])


def save_plot(b, result, output_dir, stem):
    model = result["model"]
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"), dtype=float)
    phase = phase_values(times, result["candidate"]["t0"])
    residuals = observed - predicted
    order = np.argsort(phase)
    bins = np.linspace(-0.2, 0.8, 80)
    mids, meds = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (phase >= lo) & (phase < hi)
        if np.any(mask):
            mids.append(0.5 * (lo + hi))
            meds.append(np.median(residuals[mask]))

    fig, (ax, rx) = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    ax.scatter(phase, observed, s=12, color="0.25", alpha=0.55, linewidths=0, label="binned TESS")
    ax.plot(phase[order], predicted[order], color="#dc2626", lw=2.2, label="PHOEBE model")
    ax.legend(loc="best", fontsize=9)
    ax.set_ylabel("Flux")
    ax.set_title(
        f"{TARGET} detached search: chi2={result['chi2']:.1f}, "
        f"rms={result['rms']:.6f}, trend={result['trend_score']:.6f}"
    )
    rx.axhline(0.0, color="0.45", lw=1)
    rx.scatter(phase, residuals, s=12, color="0.25", alpha=0.55, linewidths=0)
    rx.plot(mids, meds, color="#2563eb", lw=1.8, label="residual median")
    rx.legend(loc="best", fontsize=9)
    rx.set_xlabel("Phase")
    rx.set_ylabel("Obs-model")
    ax.set_xlim(-0.2, 0.8)
    fig.savefig(output_dir / f"{stem}_diagnostic.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    times, fluxes, sigmas = load_light_curve(INPUT_LC)
    b = build_bundle(times, fluxes, sigmas)

    baseline = evaluate(b, candidate_defaults(), "baseline_detached", keep_model=True)
    coarse_best, coarse_results = run_candidates(b, candidate_grid(), "coarse")
    local_results = []
    if os.environ.get("J1122_SKIP_LOCAL", "false").lower() not in ("1", "true", "yes"):
        _, local_results = run_candidates(b, local_candidates(coarse_best["candidate"]), "local")
    all_results = coarse_results + local_results
    best_objective = evaluate(
        b,
        best_valid(all_results, "objective")["candidate"],
        "j1122_detached_best_objective",
        keep_model=True,
    )
    best_rms = evaluate(
        b,
        best_valid(all_results, "rms")["candidate"],
        "j1122_detached_best_rms",
        keep_model=True,
    )

    save_plot(b, baseline, OUTPUT_DIR, "baseline_detached")
    save_plot(b, best_objective, OUTPUT_DIR, "j1122_detached_best_objective")
    save_plot(b, best_rms, OUTPUT_DIR, "j1122_detached_best_rms")
    b.save(str(OUTPUT_DIR / "j1122_detached_best.phoebe"))

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target": TARGET,
        "input_lc": str(INPUT_LC.relative_to(REPO_DIR)),
        "period": PERIOD,
        "t0_initial": T0,
        "sigma_floor": float(os.environ.get("J1122_SIGMA_FLOOR", "0.0004")),
        "baseline": baseline,
        "coarse_best": coarse_best,
        "best_objective": best_objective,
        "best_rms": best_rms,
        "n_coarse_results": len(coarse_results),
        "n_local_results": len(local_results),
    }
    with (OUTPUT_DIR / "j1122_detached_search_result.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    with (OUTPUT_DIR / "j1122_detached_search_all_results.json").open("w", encoding="utf-8") as handle:
        json.dump({"coarse": coarse_results, "local": local_results}, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
