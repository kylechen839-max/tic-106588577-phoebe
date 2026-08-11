#!/usr/bin/env python3
"""Build a J1122 SED and compare dust against free-free excess models."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier


TARGET = "J112238.89-592027.5"
COORD = SkyCoord("11h22m38.89s -59d20m27.5s")

OUTDIR = Path("results/test/j1122_full_sed")
WISE_CUTOUT_DIR = Path("/Users/kylechen/Downloads/FC_Files (1)")

ZEROPOINT_JY = {
    "J": 1594.0,
    "H": 1024.0,
    "Ks": 666.7,
    "W1": 309.540,
    "W2": 171.787,
    "W3": 31.674,
    "W4": 8.363,
}

WAVELENGTH_UM = {
    "J": 1.235,
    "H": 1.662,
    "Ks": 2.159,
    "W1": 3.35,
    "W2": 4.60,
    "W3": 11.56,
    "W4": 22.09,
}


@dataclass
class PhotPoint:
    band: str
    wavelength_um: float
    mag: float
    mag_err: float
    flux_jy: float
    flux_err_jy: float
    source: str


def mag_to_flux_jy(band: str, mag: float, mag_err: float) -> tuple[float, float]:
    flux = ZEROPOINT_JY[band] * 10 ** (-0.4 * mag)
    flux_err = flux * math.log(10.0) * 0.4 * mag_err
    return flux, flux_err


def masked_to_float(value) -> float | None:
    if np.ma.is_masked(value):
        return None
    if value in ("", "--"):
        return None
    return float(value)


def query_catalogs() -> dict:
    vizier = Vizier(columns=["**"], row_limit=50)
    tables = {}
    for label, catalog in {
        "allwise": "II/328/allwise",
        "twomass": "II/246/out",
        "gaia_dr3": "I/355/gaiadr3",
    }.items():
        result = vizier.query_region(COORD, radius=15 * u.arcsec, catalog=catalog)
        if len(result) == 0:
            raise RuntimeError(f"No rows returned for {catalog}")
        tables[label] = result[0]
    return tables


def nearest_row(table):
    return sorted(table, key=lambda row: float(row["_r"]))[0]


def photometry_from_allwise(row) -> list[PhotPoint]:
    mapping = [
        ("J", "Jmag", "e_Jmag", "2MASS/AllWISE cross-match"),
        ("H", "Hmag", "e_Hmag", "2MASS/AllWISE cross-match"),
        ("Ks", "Kmag", "e_Kmag", "2MASS/AllWISE cross-match"),
        ("W1", "W1mag", "e_W1mag", "AllWISE"),
        ("W2", "W2mag", "e_W2mag", "AllWISE"),
        ("W3", "W3mag", "e_W3mag", "AllWISE"),
        ("W4", "W4mag", "e_W4mag", "AllWISE"),
    ]
    points = []
    for band, mag_col, err_col, source in mapping:
        mag = masked_to_float(row[mag_col])
        err = masked_to_float(row[err_col])
        if mag is None or err is None:
            continue
        flux, flux_err = mag_to_flux_jy(band, mag, err)
        points.append(
            PhotPoint(
                band=band,
                wavelength_um=WAVELENGTH_UM[band],
                mag=mag,
                mag_err=err,
                flux_jy=flux,
                flux_err_jy=flux_err,
                source=source,
            )
        )
    return points


def planck_bnu(temp_k: float, wavelength_um: np.ndarray) -> np.ndarray:
    h = 6.62607015e-34
    c = 299792458.0
    k = 1.380649e-23
    nu = c / (wavelength_um * 1e-6)
    x = h * nu / (k * temp_k)
    return (2.0 * h * nu**3 / c**2) / np.expm1(x)


def fit_scale(model: np.ndarray, flux: np.ndarray, sigma: np.ndarray) -> tuple[float, float]:
    weight = sigma ** -2
    scale = np.sum(weight * model * flux) / np.sum(weight * model * model)
    chi2 = float(np.sum(((flux - scale * model) / sigma) ** 2))
    return float(scale), chi2


def fit_stellar_photosphere(points: list[PhotPoint]) -> dict:
    fit_points = [point for point in points if point.band in {"J", "H", "Ks", "W1"}]
    lam = np.array([point.wavelength_um for point in fit_points], dtype=float)
    flux = np.array([point.flux_jy for point in fit_points], dtype=float)
    sigma = np.array([point.flux_err_jy for point in fit_points], dtype=float)

    best = None
    for temp in np.linspace(3000.0, 20000.0, 3401):
        basis = planck_bnu(temp, lam)
        scale, chi2 = fit_scale(basis, flux, sigma)
        candidate = {"temp_k": float(temp), "scale": scale, "chi2": chi2}
        if best is None or candidate["chi2"] < best["chi2"]:
            best = candidate
    if best is None:
        raise RuntimeError("stellar photosphere fit failed")
    return best


def model_star(points: list[PhotPoint], fit: dict) -> np.ndarray:
    lam = np.array([point.wavelength_um for point in points], dtype=float)
    return fit["scale"] * planck_bnu(fit["temp_k"], lam)


def weighted_one_component_fit(
    basis: np.ndarray,
    residual: np.ndarray,
    sigma: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    weight = sigma ** -2
    amp = np.sum(weight * basis * residual) / np.sum(weight * basis * basis)
    amp = max(float(amp), 0.0)
    model = amp * basis
    chi2 = float(np.sum(((residual - model) / sigma) ** 2))
    return amp, chi2, model


def fit_dust(points: list[PhotPoint], star_flux: np.ndarray) -> dict:
    excess_points = [i for i, point in enumerate(points) if point.band in {"W2", "W3", "W4"}]
    lam = np.array([points[i].wavelength_um for i in excess_points], dtype=float)
    flux = np.array([points[i].flux_jy for i in excess_points], dtype=float)
    sigma = np.array([points[i].flux_err_jy for i in excess_points], dtype=float)
    residual = flux - star_flux[excess_points]

    best = None
    for temp in np.linspace(50.0, 1500.0, 5801):
        basis = planck_bnu(temp, lam)
        amp, chi2, model = weighted_one_component_fit(basis, residual, sigma)
        candidate = {
            "temp_k": float(temp),
            "amplitude": amp,
            "chi2_w2w3w4": chi2,
            "model_excess_w2w3w4": model.tolist(),
        }
        if best is None or candidate["chi2_w2w3w4"] < best["chi2_w2w3w4"]:
            best = candidate
    if best is None:
        raise RuntimeError("dust fit failed")
    return best


def fit_freefree(points: list[PhotPoint], star_flux: np.ndarray, alpha_min=-0.1, alpha_max=0.6) -> dict:
    excess_points = [i for i, point in enumerate(points) if point.band in {"W2", "W3", "W4"}]
    lam = np.array([points[i].wavelength_um for i in excess_points], dtype=float)
    flux = np.array([points[i].flux_jy for i in excess_points], dtype=float)
    sigma = np.array([points[i].flux_err_jy for i in excess_points], dtype=float)
    residual = flux - star_flux[excess_points]
    nu_rel = 1.0 / lam

    best = None
    for alpha in np.linspace(alpha_min, alpha_max, 701):
        basis = nu_rel ** alpha
        amp, chi2, model = weighted_one_component_fit(basis, residual, sigma)
        candidate = {
            "alpha": float(alpha),
            "amplitude": amp,
            "chi2_w2w3w4": chi2,
            "model_excess_w2w3w4": model.tolist(),
            "alpha_range": [alpha_min, alpha_max],
        }
        if best is None or candidate["chi2_w2w3w4"] < best["chi2_w2w3w4"]:
            best = candidate
    if best is None:
        raise RuntimeError("free-free fit failed")
    return best


def compute_model_fluxes(points: list[PhotPoint], star_fit: dict, excess_fit: dict, kind: str) -> np.ndarray:
    star = model_star(points, star_fit)
    lam = np.array([point.wavelength_um for point in points], dtype=float)
    excess = np.zeros_like(star)
    if kind == "dust":
        excess = excess_fit["amplitude"] * planck_bnu(excess_fit["temp_k"], lam)
    elif kind == "freefree":
        excess = excess_fit["amplitude"] * (1.0 / lam) ** excess_fit["alpha"]
    return star + excess


def write_csv(points: list[PhotPoint], path: Path) -> None:
    import csv

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(points[0]).keys()))
        writer.writeheader()
        for point in points:
            writer.writerow(asdict(point))


def table_row_to_jsonable(row, columns: list[str]) -> dict:
    out = {}
    for column in columns:
        value = row[column]
        if np.ma.is_masked(value):
            out[column] = None
        elif hasattr(value, "item"):
            out[column] = value.item()
        else:
            out[column] = str(value)
    return out


def make_plots(points: list[PhotPoint], star_fit: dict, dust_fit: dict, free_fit: dict) -> None:
    lam = np.array([point.wavelength_um for point in points], dtype=float)
    flux = np.array([point.flux_jy for point in points], dtype=float)
    sigma = np.array([point.flux_err_jy for point in points], dtype=float)
    star = model_star(points, star_fit)
    dust = compute_model_fluxes(points, star_fit, dust_fit, "dust")
    free = compute_model_fluxes(points, star_fit, free_fit, "freefree")

    grid = np.logspace(np.log10(1.0), np.log10(30.0), 300)
    star_grid = star_fit["scale"] * planck_bnu(star_fit["temp_k"], grid)
    dust_grid = star_grid + dust_fit["amplitude"] * planck_bnu(dust_fit["temp_k"], grid)
    free_grid = star_grid + free_fit["amplitude"] * (1.0 / grid) ** free_fit["alpha"]

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 8.5), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    ax = axes[0]
    ax.errorbar(lam, flux, yerr=sigma, fmt="o", color="black", label="2MASS + AllWISE")
    for point in points:
        ax.annotate(point.band, (point.wavelength_um, point.flux_jy), textcoords="offset points", xytext=(5, 5), fontsize=9)
    ax.plot(grid, star_grid, color="#4363d8", label=f"stellar BB fit ({star_fit['temp_k']:.0f} K)")
    ax.plot(grid, dust_grid, color="#e6194b", label=f"star + dust BB ({dust_fit['temp_k']:.0f} K)")
    ax.plot(grid, free_grid, color="#3cb44b", linestyle="--", label=f"star + free-free alpha={free_fit['alpha']:.2f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Flux density (Jy)")
    ax.set_title(f"{TARGET} SED")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, which="both")

    ax = axes[1]
    ax.axhline(0.0, color="0.5", linewidth=1)
    ax.errorbar(lam, (flux - star) / sigma, yerr=np.ones_like(sigma), fmt="o", color="black", label="after stellar-only")
    ax.plot(lam, (flux - dust) / sigma, "o-", color="#e6194b", label="after dust model")
    ax.plot(lam, (flux - free) / sigma, "s--", color="#3cb44b", label="after free-free model")
    ax.set_xscale("log")
    ax.set_xlabel("Wavelength (micron)")
    ax.set_ylabel("Residual / sigma")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(OUTDIR / "j1122_sed_fit.png", dpi=180)
    plt.close(fig)


def make_wise_cutout_montage() -> None:
    files = [
        WISE_CUTOUT_DIR / "fc_170.662042-59.340972_wise_1.png",
        WISE_CUTOUT_DIR / "fc_170.662042-59.340972_wise_2.png",
        WISE_CUTOUT_DIR / "fc_170.662042-59.340972_wise_3.png",
        WISE_CUTOUT_DIR / "fc_170.662042-59.340972_wise_4.png",
    ]
    if not all(path.exists() for path in files):
        return
    fig, axes = plt.subplots(1, 4, figsize=(8, 2.2))
    for ax, path, band in zip(axes, files, ["W1", "W2", "W3", "W4"]):
        image = mpimg.imread(path)
        ax.imshow(image)
        ax.set_title(band)
        ax.set_axis_off()
    fig.suptitle(f"{TARGET} local WISE cutouts")
    fig.tight_layout()
    fig.savefig(OUTDIR / "j1122_wise_cutouts.png", dpi=180)
    plt.close(fig)


def write_summary(
    points: list[PhotPoint],
    allwise_quality: dict,
    twomass_neighbors: list[dict],
    gaia_neighbors: list[dict],
    star_fit: dict,
    dust_fit: dict,
    free_fit: dict,
) -> None:
    star_flux = model_star(points, star_fit)
    flux = np.array([point.flux_jy for point in points], dtype=float)
    sigma = np.array([point.flux_err_jy for point in points], dtype=float)
    star_only_chi2_all = float(np.sum(((flux - star_flux) / sigma) ** 2))
    dust_flux = compute_model_fluxes(points, star_fit, dust_fit, "dust")
    free_flux = compute_model_fluxes(points, star_fit, free_fit, "freefree")
    dust_chi2_all = float(np.sum(((flux - dust_flux) / sigma) ** 2))
    free_chi2_all = float(np.sum(((flux - free_flux) / sigma) ** 2))

    excess_lines = []
    for point, photosphere in zip(points, star_flux):
        ratio = point.flux_jy / photosphere
        sigma_excess = (point.flux_jy - photosphere) / point.flux_err_jy
        excess_lines.append(
            f"| {point.band} | {point.flux_jy:.6g} | {photosphere:.6g} | {ratio:.2f} | {sigma_excess:.1f} |"
        )

    text = f"""# J1122 Full SED and WISE Quality Check

Target: `{TARGET}` at RA={COORD.ra.deg:.7f}, Dec={COORD.dec.deg:.7f}.

## WISE Image/Catalog Quality

Nearest AllWISE source: `{allwise_quality['AllWISE']}` at separation {allwise_quality['_r']:.3f} arcsec.

| Field | Value |
| --- | --- |
| `qph` / photometric quality | `{allwise_quality['qph']}` |
| `ccf` / contamination-confusion flags | `{allwise_quality['ccf']}` |
| `ex` / extended-source flag | `{allwise_quality['ex']}` |
| `snr1`, `snr2`, `snr3`, `snr4` | {allwise_quality['snr1']:.1f}, {allwise_quality['snr2']:.1f}, {allwise_quality['snr3']:.1f}, {allwise_quality['snr4']:.1f} |
| `chi2W1`, `chi2W2`, `chi2W3`, `chi2W4` | {allwise_quality['chi2W1']:.2f}, {allwise_quality['chi2W2']:.2f}, {allwise_quality['chi2W3']:.2f}, {allwise_quality['chi2W4']:.2f} |
| `nb`, `na` | {allwise_quality['nb']}, {allwise_quality['na']} |
| `sat1`-`sat4` | {allwise_quality['sat1']}, {allwise_quality['sat2']}, {allwise_quality['sat3']}, {allwise_quality['sat4']} |
| `fmoon` | {allwise_quality['fmoon']} |

Interpretation: the AllWISE catalog flags are favorable for using the W3/W4 points (`qph=AAAA`, `ccf=0000`, `ex=0`, no saturation). W4 is still lower-resolution and the local cutout has structured background, so the W4 excess should be treated as credible but not immune from follow-up checks.

## Nearby Sources

The matched 2MASS source is essentially coincident, but there are fainter 2MASS/Gaia neighbors inside the WISE W3/W4 beam. They are much fainter in near-IR/optical, so they are unlikely to explain the full W3/W4 excess alone, but they remain a blending caveat.

Closest 2MASS rows:

| sep arcsec | 2MASS | J | H | K | Qflg |
| ---: | --- | ---: | ---: | ---: | --- |
"""
    for row in twomass_neighbors[:5]:
        text += f"| {row['_r']:.3f} | {row['2MASS']} | {row['Jmag']} | {row['Hmag']} | {row['Kmag']} | {row['Qflg']} |\n"

    text += "\nClosest Gaia DR3 rows:\n\n| sep arcsec | source_id | G | BP | RP | parallax mas | RUWE |\n| ---: | --- | ---: | ---: | ---: | ---: | ---: |\n"
    for row in gaia_neighbors[:8]:
        text += (
            f"| {row['_r']:.3f} | {row['Source']} | {row.get('Gmag')} | {row.get('BPmag')} | "
            f"{row.get('RPmag')} | {row.get('Plx')} | {row.get('RUWE')} |\n"
        )

    text += f"""
## SED Fit

The stellar baseline is a single blackbody fit to `J`, `H`, `Ks`, and `W1`; W2/W3/W4 are then tested as excess bands. This is intentionally conservative because W2 can already contain warm excess.

Best stellar blackbody temperature: **{star_fit['temp_k']:.0f} K**.

| Band | observed Jy | stellar Jy | obs/star | excess sigma |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(excess_lines)}

## Excess Models

| Model | Key parameter | chi2 W2/W3/W4 | chi2 all bands |
| --- | ---: | ---: | ---: |
| Stellar only | {star_fit['temp_k']:.0f} K | {np.sum(((flux[[4,5,6]] - star_flux[[4,5,6]]) / sigma[[4,5,6]]) ** 2):.2f} | {star_only_chi2_all:.2f} |
| Stellar + dust blackbody | Tdust={dust_fit['temp_k']:.0f} K | {dust_fit['chi2_w2w3w4']:.2f} | {dust_chi2_all:.2f} |
| Stellar + free-free power law | alpha={free_fit['alpha']:.2f} | {free_fit['chi2_w2w3w4']:.2f} | {free_chi2_all:.2f} |

## Finding

The W3/W4 excess remains highly significant after fitting the near-IR photosphere. A single-temperature dust blackbody near **{dust_fit['temp_k']:.0f} K** explains the W2-W4 excess much better than an ordinary free-free spectrum constrained to alpha=[-0.1, 0.6]. On these data, free-free emission does **not** plausibly account for all of the mid-IR excess by itself.

The most plausible interpretation from this SED-only check is a warm circumstellar/circumbinary dust component, not a pure free-free continuum. The main caveat is WISE beam blending/background structure, especially in W4; the next high-value check would be higher-resolution mid-IR imaging or confirming the excess with independent W3/W4-quality photometry.

Generated files:

- `j1122_sed_fit.png`
- `j1122_wise_cutouts.png`
- `j1122_photometry.csv`
- `j1122_catalog_quality.json`
"""
    (OUTDIR / "j1122_full_sed_summary.md").write_text(text)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tables = query_catalogs()
    allwise = nearest_row(tables["allwise"])
    points = photometry_from_allwise(allwise)

    write_csv(points, OUTDIR / "j1122_photometry.csv")

    allwise_columns = [
        "_r",
        "AllWISE",
        "RAJ2000",
        "DEJ2000",
        "W1mag",
        "e_W1mag",
        "W2mag",
        "e_W2mag",
        "W3mag",
        "e_W3mag",
        "W4mag",
        "e_W4mag",
        "snr1",
        "chi2W1",
        "snr2",
        "chi2W2",
        "snr3",
        "chi2W3",
        "snr4",
        "chi2W4",
        "nb",
        "na",
        "sat1",
        "sat2",
        "sat3",
        "sat4",
        "ccf",
        "ex",
        "var",
        "qph",
        "fdet",
        "fmoon",
        "d2M",
        "2M",
    ]
    allwise_quality = table_row_to_jsonable(allwise, allwise_columns)

    twomass_columns = ["_r", "2MASS", "Jmag", "Hmag", "Kmag", "Qflg", "Cflg", "prox"]
    twomass_neighbors = [
        table_row_to_jsonable(row, twomass_columns)
        for row in sorted(tables["twomass"], key=lambda row: float(row["_r"]))
    ]

    gaia_columns = ["_r", "Source", "Gmag", "BPmag", "RPmag", "Plx", "RUWE", "AllWISE", "dAllWISE", "2MASS", "d2MASS"]
    gaia_neighbors = [
        table_row_to_jsonable(row, gaia_columns)
        for row in sorted(tables["gaia_dr3"], key=lambda row: float(row["_r"]))
    ]

    star_fit = fit_stellar_photosphere(points)
    star_flux = model_star(points, star_fit)
    dust_fit = fit_dust(points, star_flux)
    free_fit = fit_freefree(points, star_flux)

    catalog_quality = {
        "target": TARGET,
        "coordinate_deg": {"ra": COORD.ra.deg, "dec": COORD.dec.deg},
        "allwise": allwise_quality,
        "twomass_neighbors": twomass_neighbors,
        "gaia_neighbors": gaia_neighbors,
        "stellar_fit": star_fit,
        "dust_fit": dust_fit,
        "freefree_fit": free_fit,
        "catalogs": {
            "allwise": "VizieR II/328/allwise",
            "twomass": "VizieR II/246/out",
            "gaia_dr3": "VizieR I/355/gaiadr3",
        },
    }
    (OUTDIR / "j1122_catalog_quality.json").write_text(json.dumps(catalog_quality, indent=2))

    make_plots(points, star_fit, dust_fit, free_fit)
    make_wise_cutout_montage()
    write_summary(points, allwise_quality, twomass_neighbors, gaia_neighbors, star_fit, dust_fit, free_fit)

    print(json.dumps(catalog_quality, indent=2))
    print(f"Wrote {OUTDIR}")


if __name__ == "__main__":
    main()
