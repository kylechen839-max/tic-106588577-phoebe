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
OUTPUT_DIR = REPO_DIR / "outputs" / "residual_trend_search"
PERIOD = 1.0118536926383312
T0 = 1492.5441099999998
RNG_SEED = 106588577


BASE_CANDIDATE = {
    "fill_factor": 0.995,
    "incl": 159.59668610123063,
    "t0": 1492.54561,
    "requiv_primary": 1.5360075942651437,
    "teffratio": 0.7208571428571429,
    "gravb_primary": 0.9,
    "gravb_secondary": 1.0,
    "ld_primary": 0.5,
    "ld_secondary": 0.5,
    "sync_primary": 1.0,
    "sync_secondary": 1.0,
    "q": 1.0,
    "sma": 10.12605462848939,
}


def load_input(path):
    data = np.loadtxt(path, comments="#").T
    times, fluxes, sigmas = data[0], data[1], data[2]
    sigmas = np.maximum(sigmas, np.nanmedian(sigmas))
    return times, fluxes, sigmas


def phase_values(times, t0):
    phase = ((times - t0) / PERIOD) % 1.0
    return np.where(phase > 0.8, phase - 1.0, phase)


def rebin(times, fluxes, sigmas, n_bins):
    phase = phase_values(times, T0)
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
    out = np.asarray(rows, dtype=float)
    return out[:, 0], out[:, 1], np.maximum(out[:, 2], np.nanmedian(out[:, 2]))


def prepare_bundle(times, fluxes, sigmas):
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
    return b


def set_if_possible(b, twig, value):
    try:
        b.set_value(twig, value)
        return None
    except Exception as exc:
        return str(exc)


def set_candidate(b, candidate):
    errors = []
    for twig, key in [
        ("q@binary@orbit@component", "q"),
        ("sma@binary@orbit@component", "sma"),
        ("incl@binary@orbit@component", "incl"),
        ("t0_supconj@binary@orbit@component", "t0"),
        ("requiv@primary@star@component", "requiv_primary"),
        ("gravb_bol@primary@star@component", "gravb_primary"),
        ("gravb_bol@secondary@star@component", "gravb_secondary"),
        ("syncpar@primary@star@component", "sync_primary"),
        ("syncpar@secondary@star@component", "sync_secondary"),
    ]:
        err = set_if_possible(b, twig, candidate[key])
        if err:
            errors.append(f"{twig}: {err}")

    err = set_if_possible(b, "ld_coeffs@primary@lc01@lc@dataset", [candidate["ld_primary"]])
    if err:
        errors.append(f"ld primary: {err}")
    err = set_if_possible(b, "ld_coeffs@secondary@lc01@lc@dataset", [candidate["ld_secondary"]])
    if err:
        errors.append(f"ld secondary: {err}")

    teff_primary = b.get_value("teff@primary@star@component")
    err = set_if_possible(
        b,
        "teff@secondary@star@component",
        candidate["teffratio"] * teff_primary,
    )
    if err:
        errors.append(f"teff secondary: {err}")

    try:
        requiv_max = b.get_value("requiv_max@secondary@star@component")
        err = set_if_possible(
            b,
            "requiv@secondary@star@component",
            candidate["fill_factor"] * requiv_max,
        )
        if err:
            errors.append(f"requiv secondary: {err}")
    except Exception as exc:
        errors.append(f"requiv secondary max: {exc}")

    if errors:
        raise ValueError("; ".join(errors))


def residual_trend_metrics(phase, residuals):
    primary = np.abs(phase) < 0.045
    secondary = np.abs(phase - 0.5) < 0.035
    shoulders = ((np.abs(phase) >= 0.045) & (np.abs(phase) < 0.11)) | (
        (np.abs(phase - 0.5) >= 0.035) & (np.abs(phase - 0.5) < 0.09)
    )

    def rms(mask):
        return float(np.sqrt(np.mean(residuals[mask] ** 2)))

    def mean_abs(mask):
        return float(np.mean(np.abs(residuals[mask])))

    def mean(mask):
        return float(np.mean(residuals[mask]))

    # This is the specific thing we are trying to drive down: coherent offsets
    # in the eclipse cores and shoulders, rather than global scatter alone.
    trend_score = (
        3.0 * abs(mean(primary))
        + 2.0 * abs(mean(secondary))
        + 1.5 * mean_abs(primary)
        + 1.2 * mean_abs(secondary)
        + 0.8 * mean_abs(shoulders)
    )

    return {
        "primary_core_mean": mean(primary),
        "secondary_core_mean": mean(secondary),
        "primary_core_rms": rms(primary),
        "secondary_core_rms": rms(secondary),
        "shoulder_mean_abs": mean_abs(shoulders),
        "trend_score": float(trend_score),
    }


def evaluate(b, candidate, model):
    set_candidate(b, candidate)
    b.run_compute(compute="phoebe01", model=model, overwrite=True)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    phase = phase_values(times, candidate["t0"])
    residuals = observed - predicted
    trend = residual_trend_metrics(phase, residuals)
    rms = float(np.sqrt(np.mean(residuals**2)))
    chi2 = float(b.calculate_chi2(model=model, dataset="lc01"))
    score = rms + trend["trend_score"]
    return {
        "model": model,
        "candidate": dict(candidate),
        "chi2": chi2,
        "rms": rms,
        "mean_abs": float(np.mean(np.abs(residuals))),
        "max_abs": float(np.max(np.abs(residuals))),
        "objective": float(score),
        **trend,
    }


def sample_candidate(rng, center, scale):
    cand = dict(center)
    ranges = {
        "incl": (156.0, 164.0),
        "t0": (T0 - 0.012, T0 + 0.012),
        "requiv_primary": (1.0, 2.5),
        "fill_factor": (0.93, 0.999),
        "teffratio": (0.48, 0.95),
        "gravb_primary": (0.3, 1.0),
        "gravb_secondary": (0.3, 1.0),
        "ld_primary": (0.05, 0.95),
        "ld_secondary": (0.05, 0.95),
        "sync_primary": (0.75, 1.3),
        "sync_secondary": (0.75, 1.3),
        "q": (0.55, 1.8),
        "sma": (8.8, 13.5),
    }
    widths = {
        "incl": 2.5,
        "t0": 0.004,
        "requiv_primary": 0.35,
        "fill_factor": 0.03,
        "teffratio": 0.12,
        "gravb_primary": 0.2,
        "gravb_secondary": 0.2,
        "ld_primary": 0.25,
        "ld_secondary": 0.25,
        "sync_primary": 0.16,
        "sync_secondary": 0.16,
        "q": 0.35,
        "sma": 1.0,
    }
    for key, (lo, hi) in ranges.items():
        value = center[key] + rng.normal(0, widths[key] * scale)
        cand[key] = float(np.clip(value, lo, hi))
    return cand


def save_plot(b, result, output_dir, prefix):
    candidate = result["candidate"]
    model = result["model"]
    times = np.asarray(b.get_value("times@lc01@lc@dataset"), dtype=float)
    observed = np.asarray(b.get_value("fluxes@lc01@lc@dataset"), dtype=float)
    predicted = np.asarray(
        b.get_value(f"fluxes@lc01@phoebe01@{model}@lc@model"),
        dtype=float,
    )
    phase = phase_values(times, candidate["t0"])
    residuals = observed - predicted
    order = np.argsort(phase)
    bins = np.linspace(-0.2, 0.8, 101)
    mids, meds = [], []
    for lower, upper in zip(bins[:-1], bins[1:]):
        mask = (phase >= lower) & (phase < upper)
        if np.any(mask):
            mids.append(0.5 * (lower + upper))
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
    rx.axhline(0, color="0.45", lw=1)
    rx.scatter(phase, residuals, s=8, color="0.25", alpha=0.45, linewidths=0)
    rx.plot(mids, meds, color="#2563eb", lw=1.8, label="binned residual median")
    rx.legend(loc="upper right", fontsize=8)
    ax.set_ylabel("Flux")
    rx.set_ylabel("Obs-model")
    rx.set_xlabel("Phase")
    ax.set_xlim(-0.2, 0.8)
    title = (
        f"{prefix}: chi2={result['chi2']:.1f}, rms={result['rms']:.6f}, "
        f"trend={result['trend_score']:.6f}"
    )
    ax.set_title(title)
    fig.savefig(output_dir / f"{prefix}_diagnostic.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    times, fluxes, sigmas = load_input(INPUT_LC)
    search_times, search_fluxes, search_sigmas = rebin(times, fluxes, sigmas, 180)

    rng = np.random.default_rng(RNG_SEED)
    b = prepare_bundle(search_times, search_fluxes, search_sigmas)
    results = []
    failures = []

    candidates = [dict(BASE_CANDIDATE)]
    center = dict(BASE_CANDIDATE)
    for scale, count in [(1.0, 80), (0.55, 80), (0.28, 80)]:
        for _ in range(count):
            candidates.append(sample_candidate(rng, center, scale))

    best = None
    for idx, candidate in enumerate(candidates):
        model = f"trend_search_{idx:04d}"
        try:
            result = evaluate(b, candidate, model)
            results.append(result)
            if best is None or result["objective"] < best["objective"]:
                best = result
                center = dict(result["candidate"])
                print(
                    "accepted",
                    idx,
                    f"objective={result['objective']:.7f}",
                    f"rms={result['rms']:.7f}",
                    f"trend={result['trend_score']:.7f}",
                    f"chi2={result['chi2']:.2f}",
                )
            elif idx % 20 == 0:
                print(
                    "checked",
                    idx,
                    f"best_objective={best['objective']:.7f}",
                    f"best_rms={best['rms']:.7f}",
                )
        except Exception as exc:
            failures.append({"idx": idx, "candidate": candidate, "error": str(exc)})
            if idx % 20 == 0:
                print("failed", idx, exc)

    ranked = sorted(results, key=lambda item: item["objective"])
    full_b = prepare_bundle(times, fluxes, sigmas)
    final = evaluate(full_b, ranked[0]["candidate"], "residual_trend_best")
    save_plot(full_b, final, OUTPUT_DIR, "residual_trend_best")
    full_b.save(str(OUTPUT_DIR / "tic_106588577_residual_trend_best.phoebe"))

    result = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "period": PERIOD,
        "reference_t0": T0,
        "search_bins": int(search_times.size),
        "full_bins": int(times.size),
        "random_seed": RNG_SEED,
        "base_candidate": BASE_CANDIDATE,
        "best_search_result": ranked[0],
        "top_10_search_results": ranked[:10],
        "final_full_result": final,
        "n_success": len(results),
        "n_failures": len(failures),
        "failures": failures[:50],
    }
    with open(OUTPUT_DIR / "residual_trend_search_result.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    print(json.dumps(result["final_full_result"], indent=2))


if __name__ == "__main__":
    main()
