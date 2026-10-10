"""Write markdown tables of tier-A and top tier-B DD1 x EB candidates."""
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "outputs" / "crossmatch"


def fmt(df):
    star = lambda a, b, f: f.format(a) if pd.notna(a) else (f.format(b) + "*" if pd.notna(b) else "–")
    out = pd.DataFrame({
        "DD1 designation": df.designation,
        "EB surveys (indep.)": [f"{a} ({n})" for a, n in zip(df.surveys, df.n_independent)],
        "P (d)": df.period_d.map(lambda x: f"{x:.4f}"),
        "TIC Teff": df.tic_teff.map(lambda x: "–" if pd.isna(x) else f"{x:.0f}"),
        "W3σ / W4σ": [star(a, c, "{:.1f}") + " / " + star(b, d, "{:.1f}") for a, b, c, d in
                      zip(df.fixed_W3_sigma, df.fixed_W4_sigma, df.W3_sigma, df.W4_sigma)],
        "W4 obs/phot": [star(a, b, "{:.1f}") for a, b in zip(df.fixed_W4_ratio, df.W4_ratio)],
        "T_dust (K)": [star(a, b, "{:.0f}") for a, b in zip(df.fixed_dust_temp, df.dust_temp)],
        "Tmag": df.Tmag.map(lambda x: f"{x:.1f}"),
        "d (pc)": df.gaia_dist.map(lambda x: "–" if pd.isna(x) else f"{x:.0f}"),
        "RUWE": df.ruwe.map(lambda x: f"{x:.2f}"),
        "TESS sect.": df.n_tess_sectors.map(lambda x: "–" if pd.isna(x) else f"{x:.0f}"),
        "SIMBAD": df.simbad_otype.fillna(""),
        "failed checks": df.fail_reasons.fillna(""),
        "DD1 cut": df.dd_excess_pass.map({True: "pass", False: "fail"}) if "dd_excess_pass" in df else "",
        "DD1 label (good frac, sci team)": [("–" if pd.isna(g) else f"{g:.2f}, {t}") for g, t in
                                            zip(df.get("dd1_good_fraction", pd.Series(index=df.index, dtype=float)),
                                                df.get("dd1_sciteam", pd.Series(index=df.index, dtype=object)))],
        "grade": df.get("grade", df.tier),
    })
    return out.to_markdown(index=False) + "\n"


if __name__ == "__main__":
    s = pd.read_csv(OUT / "dd1_candidates_scored.csv")
    (OUT / "tierA_table.md").write_text(fmt(s[s.tier == "A"]))
    (OUT / "tierB_top40_table.md").write_text(fmt(s[s.tier == "B"].head(40)))
    print((OUT / "tierA_table.md").read_text())
