import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import phoebe


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_BUNDLE = REPO_DIR / "results" / "final" / "tic_106588577_multisector_controlled_best.phoebe"
INPUT_LC = REPO_DIR / "data" / "J071951.40-240400.6_LightCurve_multisector_binned.txt"
OUTPUT_DIR = REPO_DIR / "outputs" / "two_spot_physical_search"
PERIOD = 1.0118536926383312
RNG_SEED = 106588579

PRIMARY_SPOT = {
    "colat": 76.346839700304,
    "long": 7.3981441145742295,
    "radius": 37.54240561077953,
    "relteff": 0.7385109974350407,
}


def load_input(path):
    data = np.loadtxt(path, comments="#").T
    times, fluxes, sigmas = data[0], data[1], data[2]
    sigmas = np.maximum(sigmas, np.nanmedian(sigmas))
    return times, fluxes, sigmas


def phase_values(times, period, t0):
    phase = ((times - t0) / period) % 1.0
    return np.where(phase > 0.8, phase - 1.0, phase)


def rebin(times, fluxes, sigmas, t0, n_bins):
    phase = phase_values(times, PERIOD, t0)
    bins = np.linspace(-0.2, 0.8, n_bins + 1)
    rows = []
    for lower, upper in zip(bins[:-1], bins[1:]):
        mask = (phase >= lower) & (phase < upper)
        if not np.any(mask):
            continue
        bin_phase = 0.5 * (lower + upper)
        bin_flux = np.nanmedian(fluxes[mask])
        scatter = 1.4826 * np.nanmedian(np.abs(fluxes[mask] - bin_flux))
        if not np.isfinite(scatter) or scatter <= 0:
            scatter = np.nanmedian(sigmas[mask])
        rows.append((t0 + bin_phase * PERIOD, bin_flux, scatter / np.sqrt(mask.sum())))
    out = np.asarray(rows, dtype=float)
    return out[:, 0], out[:, 1], np.maximum(out[:, 2], np.nanmedian(out[:, 2]))


def prepare_bundle(times, fluxes, sigmas):
    b = phoebe.load(str(BASE_BUNDLE))
    b.set_value("times@lc01@lc@dataset", times)
    b.set_value("fluxes@lc01@lc@dataset", fluxes)
    b.set_value("sigmas@lc01@lc@dataset", sigmas)
    b.add_feature("spot", component="primary", feature="spot_primary_main", overwrite=True)
    b.add_feature("spot", component="secondary", feature="spot_secondary_trial", overwrite=True)
    return b


def set_named_spot(b, feature, component, candidate):
    prefix = f"@{feature}@{component}@spot@feature"
    b.set_value(f"colat{prefix}", candidate["colat"])
    b.set_value(f"long{prefix}", candidate["long"])
    b.set_value(f"radius{prefix}", candidate["radius"])
    b.set_value(f"relteff{prefix}", candidate["relteff"])


def set_spots(b, candidate):
    set_named_spot(b, "spot_primary_main", "primary", PRIMARY_SPOT)
    set_named_spot(b, "spot_secondary_trial", "secondary", candidate)


def residual_metrics(phase, residuals):
    primary = np.abs(phase) < 0.045
    secondary = np.abs(phase - 0.5) < 0.035
    shoulders = ((np.abs(phase) >= 0.045) & (np.abs(phase) < 0.11)) | (
        (np.abs(phase - 0.5) >= 0.035) & (np.abs(phase - 0.5) < 0.09)
    )
    broad = (phase > 0.15) & (phase < 0.75)

    def mean(mask):
        return float(np.mean(residuals[mask]))

    def rms(mask):
        return float(np.sqrt(np.mean(residuals[mask] ** 2)))

    def mean_abs(mask):
        return float(np.mean(np.abs(residuals[mask])))

    trend_score = (
        3.0 * abs(mean(primary))
        + 2.0 * abs(mean(secondary))
        + 1.5 * mean_abs(primary)
        + 1.2 * mean_abs(secondary)
        + 0.8 * mean_abs(shoulders)
        + 0.8 * mean_abs(broad)
    )
    return {
        "rms": float(np.sqrt(np.mean(residuals**2))),
        "mean_abs": float(np.mean(np.abs(residuals))),
        "max_abs": float(np.max(np.abs(residuals))),
        "primary_core_mean": mean(primary),
        "secondary_core_mean": mean(secondary),
        "primary_core_rms": rms(primary),
        "secondary_core_rms": rms(secondary),
        "shoulder_mean_abs": mean_abs(shoulders),
        "broad_mean_abs": mean_abs(broad),
        "trend_score": float(trend_score),
    }


def evaluate(b, candidate, model):
    set_spots(b, candidate)
    b.run_compute(compute="phoebe01", model=model, overwrite=True)
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    t0 = float(b.get_value("t0_supconj@binary@orbit@component"))
    phase = phase_values(times, PERIOD, t0)
    residuals = observed - predicted
    metrics = residual_metrics(phase, residuals)
    objective = metrics["rms"] + 0.35 * metrics["trend_score"]
    return {
        "model": model,
        "secondary_spot": dict(candidate),
        "primary_spot_fixed": dict(PRIMARY_SPOT),
        "chi2": float(b.calculate_chi2(model=model, dataset="lc01")),
        "objective": float(objective),
        **metrics,
    }


def candidate_grid():
    rng = np.random.default_rng(RNG_SEED)
    candidates = [
        {"colat": 90.0, "long": 180.0, "radius": 18.0, "relteff": 0.85},
        {"colat": 90.0, "long": 180.0, "radius": 28.0, "relteff": 0.8},
        {"colat": 70.0, "long": 170.0, "radius": 24.0, "relteff": 0.85},
        {"colat": 110.0, "long": 190.0, "radius": 24.0, "relteff": 0.85},
    ]
    for _ in range(84):
        candidates.append(
            {
                "colat": float(rng.uniform(35, 145)),
                "long": float(rng.uniform(120, 240) % 360),
                "radius": float(rng.uniform(6, 40)),
                "relteff": float(rng.uniform(0.55, 1.15)),
            }
        )
    return candidates


def local_candidates(center):
    rng = np.random.default_rng(RNG_SEED + 1)
    candidates = []
    for _ in range(56):
        candidates.append(
            {
                "colat": float(np.clip(center["colat"] + rng.normal(0, 18), 20, 160)),
                "long": float((center["long"] + rng.normal(0, 18)) % 360),
                "radius": float(np.clip(center["radius"] + rng.normal(0, 8), 4, 55)),
                "relteff": float(np.clip(center["relteff"] + rng.normal(0, 0.12), 0.45, 1.2)),
            }
        )
    return candidates


def save_plot(b, result, output_dir, prefix):
    model = result["model"]
    t0 = float(b.get_value("t0_supconj@binary@orbit@component"))
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    phase = phase_values(times, PERIOD, t0)
    residuals = observed - predicted
    order = np.argsort(phase)
    bins = np.linspace(-0.2, 0.8, 101)
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
    ax.scatter(phase, observed, s=8, color="0.25", alpha=0.45, linewidths=0)
    ax.plot(phase[order], predicted[order], color="#dc2626", lw=2)
    ax.set_ylabel("Flux")
    ax.set_title(
        f"{prefix}: chi2={result['chi2']:.1f}, rms={result['rms']:.6f}, "
        f"trend={result['trend_score']:.6f}"
    )
    rx.axhline(0, color="0.45", lw=1)
    rx.scatter(phase, residuals, s=8, color="0.25", alpha=0.45, linewidths=0)
    rx.plot(mids, meds, color="#2563eb", lw=1.8, label="binned residual median")
    rx.legend(loc="upper right", fontsize=8)
    rx.set_xlabel("Phase")
    rx.set_ylabel("Obs-model")
    ax.set_xlim(-0.2, 0.8)
    fig.savefig(output_dir / f"{prefix}_diagnostic.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    full_times, full_fluxes, full_sigmas = load_input(INPUT_LC)
    base = phoebe.load(str(BASE_BUNDLE))
    t0 = float(base.get_value("t0_supconj@binary@orbit@component"))
    search_times, search_fluxes, search_sigmas = rebin(full_times, full_fluxes, full_sigmas, t0, 180)

    search_b = prepare_bundle(search_times, search_fluxes, search_sigmas)
    results, failures = [], []
    best = None
    candidates = candidate_grid()
    for idx, candidate in enumerate(candidates):
        try:
            result = evaluate(search_b, candidate, f"two_spot_search_{idx:04d}")
            results.append(result)
            if best is None or result["objective"] < best["objective"]:
                best = result
                print(
                    "accepted",
                    idx,
                    f"objective={result['objective']:.7f}",
                    f"rms={result['rms']:.7f}",
                    f"trend={result['trend_score']:.7f}",
                    f"chi2={result['chi2']:.2f}",
                    result["secondary_spot"],
                )
            elif idx % 20 == 0:
                print("checked", idx, f"best_objective={best['objective']:.7f}")
        except Exception as exc:
            failures.append({"idx": idx, "candidate": candidate, "error": str(exc)})

    for offset, candidate in enumerate(local_candidates(best["secondary_spot"]), start=len(candidates)):
        try:
            result = evaluate(search_b, candidate, f"two_spot_search_{offset:04d}")
            results.append(result)
            if result["objective"] < best["objective"]:
                best = result
                print(
                    "accepted",
                    offset,
                    f"objective={result['objective']:.7f}",
                    f"rms={result['rms']:.7f}",
                    f"trend={result['trend_score']:.7f}",
                    f"chi2={result['chi2']:.2f}",
                    result["secondary_spot"],
                )
            elif offset % 20 == 0:
                print("checked", offset, f"best_objective={best['objective']:.7f}")
        except Exception as exc:
            failures.append({"idx": offset, "candidate": candidate, "error": str(exc)})

    ranked = sorted(results, key=lambda item: item["objective"])
    full_b = prepare_bundle(full_times, full_fluxes, full_sigmas)
    final = evaluate(full_b, ranked[0]["secondary_spot"], "two_spot_physical_best")
    save_plot(full_b, final, OUTPUT_DIR, "two_spot_physical_best")
    full_b.save(str(OUTPUT_DIR / "tic_106588577_two_spot_physical_best.phoebe"))

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_bundle": str(BASE_BUNDLE.relative_to(REPO_DIR)),
        "period": PERIOD,
        "random_seed": RNG_SEED,
        "search_bins": int(search_times.size),
        "full_bins": int(full_times.size),
        "best_search_result": ranked[0],
        "top_10_search_results": ranked[:10],
        "final_full_result": final,
        "n_success": len(results),
        "n_failures": len(failures),
        "failures": failures[:50],
    }
    with open(OUTPUT_DIR / "two_spot_physical_search_result.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result["final_full_result"], indent=2))


if __name__ == "__main__":
    main()
