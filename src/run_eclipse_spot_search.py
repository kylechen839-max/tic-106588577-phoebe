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
OUTPUT_DIR = REPO_DIR / "outputs" / "eclipse_spot_search"
PERIOD = 1.0118536926383312
RNG_SEED = 106588581


def load_input(path):
    data = np.loadtxt(path, comments="#").T
    times, fluxes, sigmas = data[0], data[1], data[2]
    sigmas = np.maximum(sigmas, np.nanmedian(sigmas))
    return times, fluxes, sigmas


def phase_values(times, t0):
    phase = ((times - t0) / PERIOD) % 1.0
    return np.where(phase > 0.8, phase - 1.0, phase)


def rebin(times, fluxes, sigmas, t0, n_bins):
    phase = phase_values(times, t0)
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
    b.add_feature("spot", component="primary", feature="primary_eclipse_spot", overwrite=True)
    return b


def set_candidate(b, candidate):
    b.set_value("incl@binary@orbit@component", candidate["incl"])
    b.set_value("t0_supconj@binary@orbit@component", candidate["t0"])
    b.set_value("colat@primary_eclipse_spot@primary@spot@feature", candidate["colat"])
    b.set_value("long@primary_eclipse_spot@primary@spot@feature", candidate["long"])
    b.set_value("radius@primary_eclipse_spot@primary@spot@feature", candidate["radius"])
    b.set_value("relteff@primary_eclipse_spot@primary@spot@feature", candidate["relteff"])


def residual_metrics(phase, residuals):
    primary_core = np.abs(phase) < 0.025
    primary_wings = (np.abs(phase) >= 0.025) & (np.abs(phase) < 0.07)
    secondary_core = np.abs(phase - 0.5) < 0.035
    out_of_eclipse = (
        ((phase > 0.10) & (phase < 0.42))
        | ((phase > 0.57) & (phase < 0.78))
        | ((phase > -0.2) & (phase < -0.08))
    )

    def mean(mask):
        return float(np.mean(residuals[mask]))

    def mean_abs(mask):
        return float(np.mean(np.abs(residuals[mask])))

    def rms(mask):
        return float(np.sqrt(np.mean(residuals[mask] ** 2)))

    primary_shape = abs(mean(primary_core)) + abs(mean(primary_wings))
    trend_score = (
        5.0 * primary_shape
        + 1.6 * mean_abs(primary_core)
        + 1.0 * mean_abs(primary_wings)
        + 1.4 * mean_abs(secondary_core)
        + 1.0 * mean_abs(out_of_eclipse)
    )
    return {
        "rms": float(np.sqrt(np.mean(residuals**2))),
        "mean_abs": float(np.mean(np.abs(residuals))),
        "max_abs": float(np.max(np.abs(residuals))),
        "primary_core_mean": mean(primary_core),
        "primary_wing_mean": mean(primary_wings),
        "primary_core_rms": rms(primary_core),
        "primary_wing_rms": rms(primary_wings),
        "secondary_core_mean": mean(secondary_core),
        "secondary_core_rms": rms(secondary_core),
        "out_of_eclipse_mean_abs": mean_abs(out_of_eclipse),
        "trend_score": float(trend_score),
    }


def evaluate(b, candidate, model):
    set_candidate(b, candidate)
    b.run_compute(compute="phoebe01", model=model, overwrite=True)
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    t0 = float(b.get_value("t0_supconj@binary@orbit@component"))
    phase = phase_values(times, t0)
    residuals = observed - predicted
    metrics = residual_metrics(phase, residuals)
    objective = metrics["rms"] + 0.30 * metrics["trend_score"]
    return {
        "model": model,
        "candidate": dict(candidate),
        "chi2": float(b.calculate_chi2(model=model, dataset="lc01")),
        "objective": float(objective),
        **metrics,
    }


def candidates():
    base = phoebe.load(str(BASE_BUNDLE))
    base_t0 = float(base.get_value("t0_supconj@binary@orbit@component"))
    base_incl = float(base.get_value("incl@binary@orbit@component"))
    out = []
    for colat in [55, 70, 85, 100, 115, 130]:
        for long in [340, 350, 0, 10, 20]:
            for radius in [6, 10, 14, 18, 24]:
                for relteff in [0.45, 0.6, 0.75, 0.9]:
                    out.append(
                        {
                            "incl": base_incl,
                            "t0": base_t0,
                            "colat": float(colat),
                            "long": float(long % 360),
                            "radius": float(radius),
                            "relteff": float(relteff),
                        }
                    )
    rng = np.random.default_rng(RNG_SEED)
    for _ in range(90):
        out.append(
            {
                "incl": float(np.clip(base_incl + rng.normal(0, 1.2), 156.0, 166.0)),
                "t0": float(base_t0 + rng.normal(0, 0.00035)),
                "colat": float(rng.uniform(45, 135)),
                "long": float(rng.choice([rng.uniform(330, 360), rng.uniform(0, 30)])),
                "radius": float(rng.uniform(4, 26)),
                "relteff": float(rng.uniform(0.35, 0.95)),
            }
        )
    return out


def save_plot(b, result, output_dir, prefix):
    model = result["model"]
    t0 = float(b.get_value("t0_supconj@binary@orbit@component"))
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    phase = phase_values(times, t0)
    residuals = observed - predicted
    order = np.argsort(phase)
    bins = np.linspace(-0.2, 0.8, 121)
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
    search_times, search_fluxes, search_sigmas = rebin(full_times, full_fluxes, full_sigmas, t0, 260)

    search_b = prepare_bundle(search_times, search_fluxes, search_sigmas)
    results, failures, best = [], [], None
    for idx, candidate in enumerate(candidates()):
        try:
            result = evaluate(search_b, candidate, f"eclipse_spot_search_{idx:04d}")
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
                    result["candidate"],
                )
            elif idx % 50 == 0:
                print("checked", idx, f"best_objective={best['objective']:.7f}")
        except Exception as exc:
            failures.append({"idx": idx, "candidate": candidate, "error": str(exc)})

    ranked = sorted(results, key=lambda item: item["objective"])
    full_b = prepare_bundle(full_times, full_fluxes, full_sigmas)
    final = evaluate(full_b, ranked[0]["candidate"], "eclipse_spot_best")
    save_plot(full_b, final, OUTPUT_DIR, "eclipse_spot_best")
    full_b.save(str(OUTPUT_DIR / "tic_106588577_eclipse_spot_best.phoebe"))

    output = {
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
        "failures": failures[:25],
    }
    with open(OUTPUT_DIR / "eclipse_spot_search_result.json", "w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    main()
