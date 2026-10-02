"""Collect outputs/candidates/*/*_phoebe_result*.json into a markdown table and a montage."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
CAND = ROOT / "outputs" / "candidates"


def yn(x):
    return "–" if x is None else ("✓" if x else "✗")


def main():
    rows, pngs = [], []
    for f in sorted(CAND.glob("*/*_phoebe_result*.json")):
        r = json.loads(f.read_text())
        ph, bb = r["photometry_check"], r["blackbody_check"]
        ex = bb.get("excess_vs_binary_photosphere", {})
        rows.append("| " + " | ".join([
            r["designation"] + (" (" + f.stem.split("_phoebe_result")[1].strip("_") + ")" if f.stem.endswith("sector") else ""),
            f"{r['period']:.5f}", f"{r['incl']:.2f}", f"{r['r1_frac']:.3f} / {r['r2_frac']:.3f}",
            f"{r['teff1']:.0f} / {r['teff2']:.0f}",
            f"{r.get('ecc', 0):.3f}" if r.get("eccentric") else "0 (fixed)",
            f"{r.get('l3_frac', 0):.2f}" + (" fit" if r.get("l3_fitted") else ""),
            f"{r['R1_rsun']:.2f} / {r['R2_rsun']:.2f}",
            f"{ph['chi2_red']:.1f}", f"{ph['rms']*1e3:.2f}",
            f"{ph.get('primary_bias_sigma', float('nan')):.1f} / {ph.get('secondary_bias_sigma', float('nan')):.1f}",
            f"{bb.get('distance_ratio_phot_over_gaia', float('nan')):.2f}",
            f"{ex.get('W3', {}).get('sigma', float('nan')):.1f} / {ex.get('W4', {}).get('sigma', float('nan')):.1f}",
            yn(ph.get("pass_fit")), yn(ph.get("pass_no_eclipse_bias")), yn(bb.get("pass_photosphere")),
            yn(bb.get("pass_distance")), yn(bb.get("pass_excess_survives")), "**" + yn(r.get("passes_all")) + "**",
        ]) + " |")
        png = f.parent / f.name.replace("_phoebe_result", "_phoebe_fit").replace(".json", ".png")
        if png.exists():
            pngs.append(png)
    head = ("| target | P (d) | i (°) | r1/r2 (R/a) | T1/T2 (K) | e | l3 | R1/R2 (R☉, assumed M) | χ²_red | rms (ppt) "
            "| eclipse bias σ (pri/sec) | d_phot/d_Gaia | W3σ/W4σ vs binary | fit | no bias | phot. | dist. | excess | all |\n"
            "|" + "---|" * 20)
    table = head + "\n" + "\n".join(rows) + "\n"
    (CAND / "phoebe_summary.md").write_text(table)
    print(table)
    if pngs:
        n = len(pngs)
        cols = 3
        fig, axs = plt.subplots((n + cols - 1) // cols, cols, figsize=(6 * cols, 4.4 * ((n + cols - 1) // cols)))
        for ax in axs.ravel():
            ax.axis("off")
        for ax, p in zip(axs.ravel(), pngs):
            ax.imshow(mpimg.imread(p))
        fig.tight_layout()
        fig.savefig(CAND / "phoebe_fits_montage.png", dpi=90)


if __name__ == "__main__":
    main()
