import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import phoebe


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_BUNDLE = REPO_DIR / "bundles" / "tic_106588577_powell_ready.phoebe"
INPUT_LC = REPO_DIR / "data" / "J071951.40-240400.6_LightCurve_multisector_binned.txt"
OUTPUT_DIR = REPO_DIR / "outputs" / "multisector_controlled_search"
PERIOD = 1.0118536926383312
T0 = 1492.5441099999998


def load_input(path):
    data = np.loadtxt(path, comments="#").T
    times, fluxes, sigmas = data[0], data[1], data[2]
    sigmas = np.maximum(sigmas, np.nanmedian(sigmas))
    return times, fluxes, sigmas


def phase_values(times):
    phase = ((times - T0) / PERIOD) % 1.0
    return np.where(phase > 0.8, phase - 1.0, phase)


def rebin_for_search(times, fluxes, sigmas, n_bins=240):
    phase = phase_values(times)
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
        rows.append((T0 + bin_phase * PERIOD, bin_flux, scatter / np.sqrt(mask.sum())))
    data = np.asarray(rows, dtype=float)
    return data[:, 0], data[:, 1], np.maximum(data[:, 2], np.nanmedian(data[:, 2]))


def prepare_bundle(times, fluxes, sigmas, fill_factor):
    b = phoebe.load(str(BASE_BUNDLE))
    for twig in b.filter(context="constraint", constraint_func="semidetached").twigs:
        try:
            b.remove_constraint(twig)
        except Exception:
            pass

    b.set_value("period@binary@orbit@component", PERIOD)
    b.set_value("t0_supconj@binary@orbit@component", T0)
    b.set_value("times@lc01@lc@dataset", times)
    b.set_value("fluxes@lc01@lc@dataset", fluxes)
    b.set_value("sigmas@lc01@lc@dataset", sigmas)
    b.set_value("pblum_mode@lc01", "dataset-scaled")

    requiv_max = b.get_value("requiv_max@secondary@star@component")
    b.set_value("requiv@secondary@star@component", fill_factor * requiv_max)
    return b


def set_candidate(b, candidate):
    if "incl" in candidate:
        b.set_value("incl@binary@orbit@component", candidate["incl"])
    if "t0" in candidate:
        b.set_value("t0_supconj@binary@orbit@component", candidate["t0"])
    if "requiv_primary" in candidate:
        b.set_value("requiv@primary@star@component", candidate["requiv_primary"])
    if "teffratio" in candidate:
        teff_primary = b.get_value("teff@primary@star@component")
        b.set_value("teff@secondary@star@component", candidate["teffratio"] * teff_primary)
    if "gravb_primary" in candidate:
        b.set_value("gravb_bol@primary@star@component", candidate["gravb_primary"])
    if "gravb_secondary" in candidate:
        b.set_value("gravb_bol@secondary@star@component", candidate["gravb_secondary"])


def evaluate(b, candidate, model):
    set_candidate(b, candidate)
    b.run_compute(compute="phoebe01", model=model, overwrite=True)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    residuals = observed - predicted
    return {
        "model": model,
        "candidate": dict(candidate),
        "chi2": float(b.calculate_chi2(model=model, dataset="lc01")),
        "rms": float(np.sqrt(np.mean(residuals**2))),
        "mean_abs": float(np.mean(np.abs(residuals))),
        "max_abs": float(np.max(np.abs(residuals))),
    }


def merge(base, **updates):
    candidate = dict(base)
    candidate.update(updates)
    return candidate


def run_stage(b, base_candidate, name, values_by_key):
    results = []
    for key, values in values_by_key.items():
        for value in values:
            candidate = merge(base_candidate, **{key: value})
            model = f"{name}_{key}_{len(results):03d}"
            try:
                result = evaluate(b, candidate, model)
                print(name, key, value, result["rms"], result["chi2"])
                results.append(result)
            except Exception as exc:
                print(name, key, value, "FAILED", exc)
                results.append(
                    {
                        "model": model,
                        "candidate": candidate,
                        "error": str(exc),
                    }
                )
    finite = [r for r in results if "rms" in r and np.isfinite(r["rms"])]
    best = min(finite, key=lambda r: r["rms"]) if finite else None
    return best, results


def save_diagnostic(b, model, output_dir):
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    phase = phase_values(times)
    residuals = observed - predicted
    order = np.argsort(phase)

    fig, (ax, rx) = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    ax.scatter(phase, observed, s=8, color="0.25", alpha=0.45, linewidths=0)
    ax.plot(phase[order], predicted[order], color="#dc2626", lw=2)
    rx.axhline(0, color="0.45", lw=1)
    rx.scatter(phase, residuals, s=8, color="0.25", alpha=0.45, linewidths=0)
    rx.plot(phase[order], residuals[order], color="#dc2626", lw=1)
    ax.set_ylabel("Flux")
    rx.set_ylabel("Obs-model")
    rx.set_xlabel("Phase")
    ax.set_xlim(-0.2, 0.8)
    ax.set_title("Controlled search best model on multi-sector binned light curve")
    fig.savefig(output_dir / "controlled_search_best_diagnostic.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    times, fluxes, sigmas = load_input(INPUT_LC)
    search_times, search_fluxes, search_sigmas = rebin_for_search(
        times,
        fluxes,
        sigmas,
        n_bins=240,
    )

    base_candidate = {
        "fill_factor": 0.995,
        "incl": 159.59668610123063,
        "t0": T0,
        "requiv_primary": 1.5360075942651437,
        "teffratio": 0.7208571428571429,
        "gravb_primary": 0.9,
        "gravb_secondary": 0.9,
    }
    all_results = []

    b = prepare_bundle(search_times, search_fluxes, search_sigmas, base_candidate["fill_factor"])
    baseline = evaluate(b, base_candidate, "search_baseline")
    all_results.append(baseline)
    current = dict(baseline["candidate"])

    stages = [
        ("stage1", {"incl": [156, 158, 160, 162, 164, 166, 168, 170, 172, 174]}),
        ("stage2", {"requiv_primary": [1.35, 1.45, 1.5360075942651437, 1.65, 1.8, 2.0]}),
        ("stage3", {"teffratio": [0.66, 0.69, 0.7208571428571429, 0.75, 0.78]}),
        ("stage4", {"gravb_primary": [0.6, 0.75, 0.9, 1.0], "gravb_secondary": [0.6, 0.75, 0.9, 1.0]}),
        ("stage5", {"t0": [T0 - 0.003, T0 - 0.0015, T0, T0 + 0.0015, T0 + 0.003]}),
    ]

    for name, values_by_key in stages:
        best, results = run_stage(b, current, name, values_by_key)
        all_results.extend(results)
        if best and best["rms"] < evaluate(b, current, f"{name}_current_check")["rms"]:
            current = dict(best["candidate"])
            print("accepted", name, current)

    full_b = prepare_bundle(times, fluxes, sigmas, current["fill_factor"])
    final = evaluate(full_b, current, "controlled_search_best")
    save_diagnostic(full_b, "controlled_search_best", OUTPUT_DIR)
    full_b.save(str(OUTPUT_DIR / "tic_106588577_multisector_controlled_best.phoebe"))

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "period": PERIOD,
        "t0": T0,
        "search_bins": int(search_times.size),
        "full_bins": int(times.size),
        "baseline_search": baseline,
        "all_search_results": all_results,
        "accepted_candidate": current,
        "final_full_result": final,
    }
    with open(OUTPUT_DIR / "controlled_search_result.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result["final_full_result"], indent=2))


if __name__ == "__main__":
    main()
