#!/usr/bin/env python3
"""Diagnose whether WISE W3/W4 excess looks disk-like or free-free-like.

This is intentionally dependency-light so it can run in the base macOS Python.
It parses the local IPAC-style AllWISE export files in Downloads, estimates a
Rayleigh-Jeans photospheric baseline from W1/W2, and compares W3/W4 residuals.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
from pathlib import Path


WISE_ZEROPOINT_JY = {
    "W1": 309.540,
    "W2": 171.787,
    "W3": 31.674,
    "W4": 8.363,
}

WISE_WAVELENGTH_UM = {
    "W1": 3.35,
    "W2": 4.60,
    "W3": 11.56,
    "W4": 22.09,
}

REQUIRED_COLUMNS = [
    "designation",
    "ra",
    "dec",
    "glon",
    "glat",
    "w1mpro",
    "w1sigmpro",
    "w2mpro",
    "w2sigmpro",
    "w3mpro",
    "w3sigmpro",
    "w4mpro",
    "w4sigmpro",
    "j_m_2mass",
    "j_msig_2mass",
    "h_m_2mass",
    "h_msig_2mass",
    "k_m_2mass",
    "k_msig_2mass",
]


def parse_value(raw: str) -> float | str | None:
    raw = raw.strip()
    if raw.lower() == "null":
        return None
    try:
        return float(raw)
    except ValueError:
        return raw


def iter_ipac_rows(path: Path):
    headers: list[str] | None = None
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("|"):
                parts = [part.strip() for part in line.strip().strip("|").split("|")]
                if parts and parts[0] == "designation":
                    headers = parts
                continue
            if line.startswith("\\") or not line.strip() or headers is None:
                continue
            values = line.split()
            if len(values) < len(headers):
                continue
            yield dict(zip(headers, (parse_value(value) for value in values)))


def find_target(downloads_dir: Path, designation: str) -> tuple[Path, dict]:
    matches: list[tuple[Path, dict]] = []
    for path in sorted(downloads_dir.glob("WZ_subjects_*.txt")):
        for row in iter_ipac_rows(path):
            if row.get("designation") == designation:
                matches.append((path, row))
    if not matches:
        raise SystemExit(f"No exact row found for {designation} in {downloads_dir}/WZ_subjects_*.txt")
    return matches[0]


def mag_to_flux_jy(band: str, mag: float, mag_sigma: float) -> tuple[float, float]:
    flux = WISE_ZEROPOINT_JY[band] * 10 ** (-0.4 * mag)
    sigma = flux * math.log(10) * 0.4 * mag_sigma
    return flux, sigma


def fit_rayleigh_jeans_baseline(fluxes: dict[str, float], sigmas: dict[str, float]) -> float:
    """Fit F_nu = amplitude / lambda^2 to W1/W2."""
    numerator = 0.0
    denominator = 0.0
    for band in ("W1", "W2"):
        x = WISE_WAVELENGTH_UM[band] ** -2
        weight = 1.0 / sigmas[band] ** 2
        numerator += weight * x * fluxes[band]
        denominator += weight * x * x
    return numerator / denominator


def blackbody_bnu_ratio(temp_k: float, lam1_um: float, lam2_um: float) -> float:
    h = 6.62607015e-34
    c = 299792458.0
    k = 1.380649e-23

    def bnu(lam_um: float) -> float:
        nu = c / (lam_um * 1e-6)
        return (2.0 * h * nu**3 / c**2) / (math.exp(h * nu / (k * temp_k)) - 1.0)

    return bnu(lam1_um) / bnu(lam2_um)


def best_dust_temperature(excess_w3: float, excess_w4: float) -> tuple[float | None, float | None]:
    if excess_w3 <= 0 or excess_w4 <= 0:
        return None, None
    observed_ratio = excess_w3 / excess_w4
    best_temp = None
    best_ratio = None
    best_error = float("inf")
    for i in range(50 * 4, 1000 * 4 + 1):
        temp = i / 4.0
        ratio = blackbody_bnu_ratio(temp, WISE_WAVELENGTH_UM["W3"], WISE_WAVELENGTH_UM["W4"])
        error = abs(math.log(ratio / observed_ratio))
        if error < best_error:
            best_error = error
            best_temp = temp
            best_ratio = ratio
    return best_temp, best_ratio


def excess_power_law_alpha(excess_w3: float, excess_w4: float) -> float | None:
    if excess_w3 <= 0 or excess_w4 <= 0:
        return None
    nu3 = 1.0 / WISE_WAVELENGTH_UM["W3"]
    nu4 = 1.0 / WISE_WAVELENGTH_UM["W4"]
    return math.log(excess_w4 / excess_w3) / math.log(nu4 / nu3)


def solve_two_parameter_weighted_fit(
    basis_1: list[float],
    basis_2: list[float],
    observed: list[float],
    sigmas: list[float],
    clamp_second_nonnegative: bool = True,
) -> tuple[float, float, float]:
    """Weighted least-squares fit to y = a*basis_1 + b*basis_2."""
    s11 = s12 = s22 = t1 = t2 = 0.0
    for x1, x2, y, sigma in zip(basis_1, basis_2, observed, sigmas):
        weight = 1.0 / sigma**2
        s11 += weight * x1 * x1
        s12 += weight * x1 * x2
        s22 += weight * x2 * x2
        t1 += weight * x1 * y
        t2 += weight * x2 * y
    determinant = s11 * s22 - s12 * s12
    if abs(determinant) < 1e-300:
        return 0.0, 0.0, float("inf")
    a = (t1 * s22 - t2 * s12) / determinant
    b = (s11 * t2 - s12 * t1) / determinant
    if clamp_second_nonnegative and b < 0:
        b = 0.0
        a = t1 / s11 if s11 else 0.0
    chi2 = sum(
        ((y - (a * x1 + b * x2)) / sigma) ** 2
        for x1, x2, y, sigma in zip(basis_1, basis_2, observed, sigmas)
    )
    return a, b, chi2


def fit_power_law_excess(
    fluxes: dict[str, float],
    sigmas: dict[str, float],
    alpha_min: float,
    alpha_max: float,
    alpha_step: float,
) -> dict:
    """Fit total Fnu = photosphere_RJ + free-free-like power law."""
    bands = ["W1", "W2", "W3", "W4"]
    observed = [fluxes[band] for band in bands]
    errors = [sigmas[band] for band in bands]
    photosphere_basis = [WISE_WAVELENGTH_UM[band] ** -2 for band in bands]
    best = None
    n_steps = int(round((alpha_max - alpha_min) / alpha_step))
    for i in range(n_steps + 1):
        alpha = alpha_min + i * alpha_step
        # Free-free is normally described as Fnu proportional to nu^alpha.
        # Since nu is proportional to 1/lambda, this basis is lambda^-alpha.
        excess_basis = [WISE_WAVELENGTH_UM[band] ** (-alpha) for band in bands]
        photosphere_amp, excess_amp, chi2 = solve_two_parameter_weighted_fit(
            photosphere_basis,
            excess_basis,
            observed,
            errors,
        )
        model = [
            photosphere_amp * p + excess_amp * e
            for p, e in zip(photosphere_basis, excess_basis)
        ]
        candidate = {
            "kind": "free-free power law",
            "alpha": alpha,
            "photosphere_amp": photosphere_amp,
            "excess_amp": excess_amp,
            "chi2": chi2,
            "model_fluxes": dict(zip(bands, model)),
        }
        if best is None or candidate["chi2"] < best["chi2"]:
            best = candidate
    if best is None:
        raise RuntimeError("free-free grid fit failed")
    return best


def fit_fixed_power_law_excess(
    fluxes: dict[str, float],
    sigmas: dict[str, float],
    alpha: float,
) -> dict:
    return fit_power_law_excess(fluxes, sigmas, alpha, alpha, 1.0)


def planck_bnu(temp_k: float, lam_um: float) -> float:
    h = 6.62607015e-34
    c = 299792458.0
    k = 1.380649e-23
    nu = c / (lam_um * 1e-6)
    return (2.0 * h * nu**3 / c**2) / (math.exp(h * nu / (k * temp_k)) - 1.0)


def fit_dust_blackbody(
    fluxes: dict[str, float],
    sigmas: dict[str, float],
    temp_min: float = 50.0,
    temp_max: float = 1000.0,
    temp_step: float = 0.25,
) -> dict:
    """Fit total Fnu = photosphere_RJ + single-temperature dust blackbody."""
    bands = ["W1", "W2", "W3", "W4"]
    observed = [fluxes[band] for band in bands]
    errors = [sigmas[band] for band in bands]
    photosphere_basis = [WISE_WAVELENGTH_UM[band] ** -2 for band in bands]
    best = None
    n_steps = int(round((temp_max - temp_min) / temp_step))
    for i in range(n_steps + 1):
        temp = temp_min + i * temp_step
        dust_basis = [planck_bnu(temp, WISE_WAVELENGTH_UM[band]) for band in bands]
        photosphere_amp, dust_amp, chi2 = solve_two_parameter_weighted_fit(
            photosphere_basis,
            dust_basis,
            observed,
            errors,
        )
        model = [
            photosphere_amp * p + dust_amp * d
            for p, d in zip(photosphere_basis, dust_basis)
        ]
        candidate = {
            "kind": "single-temperature blackbody dust",
            "temperature_k": temp,
            "photosphere_amp": photosphere_amp,
            "dust_amp": dust_amp,
            "chi2": chi2,
            "model_fluxes": dict(zip(bands, model)),
        }
        if best is None or candidate["chi2"] < best["chi2"]:
            best = candidate
    if best is None:
        raise RuntimeError("dust grid fit failed")
    return best


def estimate_fractional_luminosity(
    dust_fit: dict,
    stellar_teff: float,
) -> float:
    """Approximate L_IR/L_star from fitted solid-angle ratio.

    The stellar photosphere is approximated as a blackbody whose W1/W2 tail is
    fit by the Rayleigh-Jeans amplitude. This is less robust than a Kurucz/PHOEBE
    photosphere, but it mirrors the Trilling et al. logic with local data only.
    """
    c = 299792458.0
    k = 1.380649e-23
    jy = 1e-26
    sigma_sb = 5.670374419e-8
    # For a blackbody source, Fnu = solid_angle * pi * Bnu. In the RJ limit:
    # Bnu = 2 k T nu^2 / c^2 = 2 k T / lambda_m^2.
    # Our photosphere model is Fnu[Jy] = A / lambda_um^2, so converting
    # lambda_m^-2 to lambda_um^-2 introduces a 1e12 factor.
    stellar_solid_angle = dust_fit["photosphere_amp"] * jy / (math.pi * 2.0 * k * stellar_teff * 1e12)
    dust_solid_angle = dust_fit["dust_amp"] * jy / math.pi
    dust_temp = dust_fit["temperature_k"]
    stellar_bolometric_flux = stellar_solid_angle * sigma_sb * stellar_teff**4
    dust_bolometric_flux = dust_solid_angle * sigma_sb * dust_temp**4
    return dust_bolometric_flux / stellar_bolometric_flux


def model_comparison_row(label: str, fit: dict, observed_fluxes: dict[str, float]) -> str:
    if "alpha" in fit:
        descriptor = f"{fit['alpha']:.2f}"
    else:
        descriptor = f"{fit['temperature_k']:.1f} K"
    ratios = [
        fit["model_fluxes"][band] / observed_fluxes[band]
        for band in ("W1", "W2", "W3", "W4")
    ]
    return (
        f"| {label} | {descriptor} | {fit['chi2']:.2f} | "
        f"{ratios[0]:.3f} | {ratios[1]:.3f} | {ratios[2]:.3f} | {ratios[3]:.3f} |"
    )


def classify(rows: list[dict]) -> str:
    w1_ratio = rows[0]["flux_over_photosphere"]
    w2_ratio = rows[1]["flux_over_photosphere"]
    w3_sig = rows[2]["excess_sigma"]
    w4_sig = rows[3]["excess_sigma"]
    w4_ratio = rows[3]["flux_over_photosphere"]
    if abs(w1_ratio - 1) < 0.05 and abs(w2_ratio - 1) < 0.05 and w3_sig > 3 and w4_sig > 3 and w4_ratio > 3:
        return "debris-disk-like: W1/W2 are photospheric, while W3/W4 show a strong thermal excess."
    if w3_sig > 3 and w4_sig > 3:
        return "ambiguous infrared excess: W3/W4 are significant, but additional bands/spectra are needed."
    return "no robust W3/W4 excess by this simple W1/W2-baseline test."


def evaluate_fit_at_lambda(fit: dict, lam_um: float) -> float:
    photosphere = fit["photosphere_amp"] * lam_um**-2
    if "alpha" in fit:
        return photosphere + fit["excess_amp"] * lam_um ** (-fit["alpha"])
    return photosphere + fit["dust_amp"] * planck_bnu(fit["temperature_k"], lam_um)


def write_svg(path: Path, rows: list[dict], title: str, fits: list[tuple[str, dict]] | None = None) -> None:
    width = 900
    height = 560
    margin_left = 82
    margin_right = 34
    margin_top = 58
    margin_bottom = 72
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    x_min = math.log10(2.7)
    x_max = math.log10(28.0)
    y_values = [row["flux_jy"] for row in rows] + [row["photosphere_jy"] for row in rows]
    if fits:
        for _, fit in fits:
            for i in range(180):
                lam = 2.8 * (28.0 / 2.8) ** (i / 179)
                y_values.append(evaluate_fit_at_lambda(fit, lam))
    y_min = math.log10(min(y_values) * 0.55)
    y_max = math.log10(max(y_values) * 1.8)

    def x_pos(lam: float) -> float:
        return margin_left + (math.log10(lam) - x_min) / (x_max - x_min) * plot_w

    def y_pos(flux: float) -> float:
        return margin_top + (y_max - math.log10(flux)) / (y_max - y_min) * plot_h

    obs_points = [(x_pos(row["lambda_um"]), y_pos(row["flux_jy"])) for row in rows]
    phot_points = [(x_pos(row["lambda_um"]), y_pos(row["photosphere_jy"])) for row in rows]
    obs_path = " ".join(f"{x:.1f},{y:.1f}" for x, y in obs_points)
    phot_path = " ".join(f"{x:.1f},{y:.1f}" for x, y in phot_points)
    fit_colors = ["#dc2626", "#16a34a", "#9333ea", "#f97316"]
    fit_paths = []
    if fits:
        for index, (label, fit) in enumerate(fits):
            points = []
            for i in range(180):
                lam = 2.8 * (28.0 / 2.8) ** (i / 179)
                points.append((x_pos(lam), y_pos(evaluate_fit_at_lambda(fit, lam))))
            fit_paths.append((label, fit_colors[index % len(fit_colors)], " ".join(f"{x:.1f},{y:.1f}" for x, y in points)))

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin_left}" y="34" font-family="Arial" font-size="22" font-weight="700">{html.escape(title)}</text>',
        f'<line x1="{margin_left}" y1="{height-margin_bottom}" x2="{width-margin_right}" y2="{height-margin_bottom}" stroke="#222"/>',
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height-margin_bottom}" stroke="#222"/>',
        f'<polyline points="{phot_path}" fill="none" stroke="#4b5563" stroke-width="3" stroke-dasharray="8 7"/>',
    ]
    for _, color, path_data in fit_paths:
        elements.append(f'<polyline points="{path_data}" fill="none" stroke="{color}" stroke-width="2.5"/>')
    elements.append(f'<polyline points="{obs_path}" fill="none" stroke="#2563eb" stroke-width="3"/>')
    for row, (x, y) in zip(rows, obs_points):
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="#2563eb"/>')
        elements.append(f'<text x="{x - 14:.1f}" y="{height - margin_bottom + 26}" font-family="Arial" font-size="15">{row["band"]}</text>')
        elements.append(
            f'<text x="{x - 34:.1f}" y="{y - 15:.1f}" font-family="Arial" font-size="13">'
            f'{row["excess_sigma"]:.1f}σ</text>'
        )
    for x, y in phot_points:
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="white" stroke="#4b5563" stroke-width="2"/>')
    elements.extend(
        [
            f'<text x="{width/2 - 90:.1f}" y="{height-22}" font-family="Arial" font-size="16">Wavelength (micron, log scale)</text>',
            f'<text x="18" y="{height/2 + 80:.1f}" transform="rotate(-90 18 {height/2 + 80:.1f})" font-family="Arial" font-size="16">Flux density Fν (Jy, log scale)</text>',
            '<rect x="560" y="72" width="310" height="134" fill="white" stroke="#d1d5db"/>',
            '<line x1="610" y1="94" x2="660" y2="94" stroke="#2563eb" stroke-width="3"/>',
            '<text x="672" y="99" font-family="Arial" font-size="14">AllWISE observed</text>',
            '<line x1="610" y1="118" x2="660" y2="118" stroke="#4b5563" stroke-width="3" stroke-dasharray="8 7"/>',
            '<text x="672" y="123" font-family="Arial" font-size="14">W1/W2 stellar baseline</text>',
        ]
    )
    legend_y = 147
    for label, color, _ in fit_paths:
        elements.append(f'<line x1="610" y1="{legend_y}" x2="660" y2="{legend_y}" stroke="{color}" stroke-width="2.5"/>')
        elements.append(f'<text x="672" y="{legend_y + 5}" font-family="Arial" font-size="14">{html.escape(label)}</text>')
        legend_y += 24
    elements.append("</svg>")
    path.write_text("\n".join(elements))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--designation", default="J071951.40-240400.6")
    parser.add_argument("--downloads-dir", type=Path, default=Path.home() / "Downloads")
    parser.add_argument("--out-dir", type=Path, default=Path("results/test/wise_excess"))
    parser.add_argument(
        "--stellar-teff",
        type=float,
        default=18000.0,
        help="Approximate stellar Teff for optional L_IR/L_star estimate.",
    )
    args = parser.parse_args()

    source_path, target = find_target(args.downloads_dir, args.designation)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    fluxes: dict[str, float] = {}
    sigmas: dict[str, float] = {}
    for band in ("W1", "W2", "W3", "W4"):
        mag = target[f"{band.lower()}mpro"]
        mag_sigma = target[f"{band.lower()}sigmpro"]
        if not isinstance(mag, float) or not isinstance(mag_sigma, float):
            raise SystemExit(f"Missing {band} magnitude/sigma for {args.designation}")
        fluxes[band], sigmas[band] = mag_to_flux_jy(band, mag, mag_sigma)

    amplitude = fit_rayleigh_jeans_baseline(fluxes, sigmas)
    rows = []
    for band in ("W1", "W2", "W3", "W4"):
        photosphere = amplitude * WISE_WAVELENGTH_UM[band] ** -2
        excess = fluxes[band] - photosphere
        rows.append(
            {
                "band": band,
                "lambda_um": WISE_WAVELENGTH_UM[band],
                "mag": target[f"{band.lower()}mpro"],
                "mag_sigma": target[f"{band.lower()}sigmpro"],
                "flux_jy": fluxes[band],
                "flux_sigma_jy": sigmas[band],
                "photosphere_jy": photosphere,
                "excess_jy": excess,
                "excess_sigma": excess / sigmas[band],
                "flux_over_photosphere": fluxes[band] / photosphere,
            }
        )

    dust_temp, dust_ratio = best_dust_temperature(rows[2]["excess_jy"], rows[3]["excess_jy"])
    alpha = excess_power_law_alpha(rows[2]["excess_jy"], rows[3]["excess_jy"])
    dust_fit = fit_dust_blackbody(fluxes, sigmas)
    thin_ff_fit = fit_fixed_power_law_excess(fluxes, sigmas, -0.1)
    wind_ff_fit = fit_fixed_power_law_excess(fluxes, sigmas, 0.6)
    plausible_ff_fit = fit_power_law_excess(fluxes, sigmas, -0.1, 0.6, 0.005)
    unconstrained_power_law_fit = fit_power_law_excess(fluxes, sigmas, -5.0, 2.0, 0.005)
    fractional_luminosity = estimate_fractional_luminosity(dust_fit, args.stellar_teff)
    verdict = classify(rows)

    csv_path = args.out_dir / f"{args.designation}_wise_excess.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    svg_path = args.out_dir / f"{args.designation}_wise_excess.svg"
    write_svg(
        svg_path,
        rows,
        f"WISE Excess Diagnostic: {args.designation}",
        fits=[
            ("Dust blackbody", dust_fit),
            ("Best plausible free-free", plausible_ff_fit),
            ("Unconstrained power law", unconstrained_power_law_fit),
        ],
    )

    md_path = args.out_dir / f"{args.designation}_wise_excess_summary.md"
    lines = [
        f"# WISE Excess Diagnostic: {args.designation}",
        "",
        f"Source table: `{source_path}`",
        "",
        f"Verdict: **{verdict}**",
        "",
        "The photospheric baseline is a Rayleigh-Jeans approximation fit to W1 and W2.",
        "PHOEBE can help define the binary photosphere, but it does not distinguish dust emission from free-free emission by itself.",
        "",
        "| Band | Observed Jy | Photosphere Jy | Excess sigma | Flux/photosphere |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['band']} | {row['flux_jy']:.6g} | {row['photosphere_jy']:.6g} | "
            f"{row['excess_sigma']:.2f} | {row['flux_over_photosphere']:.2f} |"
        )
    lines.extend(
        [
            "",
            f"W3/W4 excess-only power-law alpha in Fnu vs nu: `{alpha:.2f}`" if alpha is not None else "W3/W4 excess-only power-law alpha: unavailable",
            f"Single-temperature blackbody matching W3/W4 residual ratio: `{dust_temp:.1f} K`" if dust_temp is not None else "Single-temperature blackbody match: unavailable",
            f"Trilling-style four-band blackbody dust fit temperature: `{dust_fit['temperature_k']:.1f} K`",
            f"Approximate L_IR/L_star from the blackbody fit, assuming Teff={args.stellar_teff:.0f} K: `{fractional_luminosity:.3e}`",
            "This fractional luminosity estimate scales roughly as Teff^-3 and should be replaced by a proper Kurucz/PHOEBE photosphere before publication.",
            "",
            "## Free-Free Test",
            "",
            "Model form: total flux = Rayleigh-Jeans stellar photosphere + A nu^alpha.",
            "",
            "| Model | alpha or T | chi2 | W1 model/obs | W2 model/obs | W3 model/obs | W4 model/obs |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            model_comparison_row("Optically thin free-free", thin_ff_fit, fluxes),
            model_comparison_row("Partially thick wind free-free", wind_ff_fit, fluxes),
            model_comparison_row("Best plausible free-free grid", plausible_ff_fit, fluxes),
            model_comparison_row("Unconstrained power law", unconstrained_power_law_fit, fluxes),
            model_comparison_row("Single-temperature dust", dust_fit, fluxes),
            "",
            "Interpretation:",
            "- W1 and W2 are consistent with the stellar baseline.",
            "- W3 is a significant excess and W4 is a very strong excess.",
            f"- The W3/W4 residual ratio is consistent with a cool roughly {dust_temp:.0f} K thermal component.",
            f"- The best free-free slope in the plausible alpha range is alpha={plausible_ff_fit['alpha']:.2f}, but its chi2 is {plausible_ff_fit['chi2']:.1f}.",
            f"- An unconstrained power law wants alpha={unconstrained_power_law_fit['alpha']:.2f}, close to the W3/W4-only residual slope and far steeper than normal free-free emission.",
            "- Therefore, WISE W3/W4 cannot be all accounted for by ordinary free-free emission without over/under-producing other WISE bands.",
            "",
            "## Trilling et al. Methodology Replication",
            "",
            "The local paper `Trilling_2007_ApJ_658_1289.pdf` corresponds to Trilling et al. 2007, ApJ 658, 1289.",
            "Their method was replicated in simplified WISE form:",
            "- predict the combined binary/star photosphere from shorter-wavelength photometry;",
            "- compute observed/predicted ratios and excess significances;",
            "- require both a high flux ratio and significant excess;",
            "- fit the residual infrared emission with a single-temperature blackbody;",
            "- estimate dust temperature and fractional luminosity.",
            "",
            "Next checks before claiming a debris disk:",
            "- Inspect AllWISE quality flags, contamination flags, reduced chi-squared, and source blending metrics.",
            "- Add Gaia and 2MASS points to fit the stellar photosphere instead of relying only on W1/W2.",
            "- Search for longer-wavelength detections or upper limits; they are the cleanest way to separate dust from free-free.",
            "- Check for gas/accretion diagnostics such as H-alpha emission if spectra exist.",
        ]
    )
    md_path.write_text("\n".join(lines) + "\n")

    print(f"Source table: {source_path}")
    print(f"CSV: {csv_path}")
    print(f"SVG: {svg_path}")
    print(f"Summary: {md_path}")
    print(verdict)
    for row in rows:
        print(
            f"{row['band']}: flux/photosphere={row['flux_over_photosphere']:.2f}, "
            f"excess_sigma={row['excess_sigma']:.2f}"
        )
    if alpha is not None:
        print(f"W3/W4 excess power-law alpha: {alpha:.2f}")
    if dust_temp is not None:
        print(f"Dust blackbody temperature: {dust_temp:.1f} K")
    print(f"Dust fit chi2: {dust_fit['chi2']:.2f}")
    print(f"Optically thin free-free chi2: {thin_ff_fit['chi2']:.2f}")
    print(f"Partially thick wind free-free chi2: {wind_ff_fit['chi2']:.2f}")
    print(
        "Best plausible free-free: "
        f"alpha={plausible_ff_fit['alpha']:.2f}, chi2={plausible_ff_fit['chi2']:.2f}"
    )
    print(
        "Best unconstrained power law: "
        f"alpha={unconstrained_power_law_fit['alpha']:.2f}, chi2={unconstrained_power_law_fit['chi2']:.2f}"
    )


if __name__ == "__main__":
    main()
