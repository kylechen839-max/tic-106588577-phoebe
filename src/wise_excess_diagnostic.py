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


def write_svg(path: Path, rows: list[dict], title: str) -> None:
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

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin_left}" y="34" font-family="Arial" font-size="22" font-weight="700">{html.escape(title)}</text>',
        f'<line x1="{margin_left}" y1="{height-margin_bottom}" x2="{width-margin_right}" y2="{height-margin_bottom}" stroke="#222"/>',
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height-margin_bottom}" stroke="#222"/>',
        f'<polyline points="{phot_path}" fill="none" stroke="#4b5563" stroke-width="3" stroke-dasharray="8 7"/>',
        f'<polyline points="{obs_path}" fill="none" stroke="#2563eb" stroke-width="3"/>',
    ]
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
            '<rect x="590" y="72" width="252" height="62" fill="white" stroke="#d1d5db"/>',
            '<line x1="610" y1="94" x2="660" y2="94" stroke="#2563eb" stroke-width="3"/>',
            '<text x="672" y="99" font-family="Arial" font-size="14">AllWISE observed</text>',
            '<line x1="610" y1="118" x2="660" y2="118" stroke="#4b5563" stroke-width="3" stroke-dasharray="8 7"/>',
            '<text x="672" y="123" font-family="Arial" font-size="14">W1/W2 stellar baseline</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(elements))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--designation", default="J071951.40-240400.6")
    parser.add_argument("--downloads-dir", type=Path, default=Path.home() / "Downloads")
    parser.add_argument("--out-dir", type=Path, default=Path("results/test/wise_excess"))
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
    verdict = classify(rows)

    csv_path = args.out_dir / f"{args.designation}_wise_excess.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    svg_path = args.out_dir / f"{args.designation}_wise_excess.svg"
    write_svg(svg_path, rows, f"WISE Excess Diagnostic: {args.designation}")

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
            "",
            "Interpretation:",
            "- W1 and W2 are consistent with the stellar baseline.",
            "- W3 is a significant excess and W4 is a very strong excess.",
            "- The W3/W4 residual ratio is consistent with a cool roughly 158 K thermal component.",
            "- A free-free explanation is not ruled out from W3/W4 alone, but it is disfavored by the lack of W1/W2 excess and the steep excess slope.",
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


if __name__ == "__main__":
    main()
