"""Add Disk Detective 1.0 excess criteria and DD1.0 volunteer labels as columns of the scored DD1 table.

Columns added to outputs/crossmatch/dd1_candidates_scored.csv (re-running replaces them):
  dd_w1w4, dd_w1w4_sigma, dd_w4snr      [W1]-[W4] colour, its significance, W4 S/N (1.0857 / sigma_W4)
  dd_ccf_ok, dd_ext_ok                  AllWISE cc_flags W4 character == '0', ext_flg == 0
  dd_mstar                              Teff < 3900 K (Gaia GSP-Phot, else TIC)
  dd_w4rchi2                            not available offline (NaN): w4rchi2 is in neither WZSubs nor the cache
  dd_excess_pass                        DD1.0 cut without w4rchi2 (Kuchner+2016, Silverberg+2018)
  dd1_zooniverse_id, dd1_good_fraction, dd1_majority_label, dd1_classifiers, dd1_sciteam
                                        volunteer labels from data/catalogs/DD_1.0ObjectsFromMAST.csv
  dd_vs_ours                            'agree', 'DD1 only' or 'ours only' (DD1 cut vs our c_ir_excess)
  grade                                 tier plus '+' (sci team YES) or '-' (DD1 cut fails, or volunteer
                                        majority is not "good"); unlabelled objects keep the bare tier

Usage: .venv-tess/bin/python src/crossmatch/add_dd1_criteria.py
"""
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv"
s = pd.read_csv(P)
s = s[[c for c in s.columns if not (c.startswith("dd_") or c.startswith("dd1_") or c == "grade")]]
aw = pd.read_csv(ROOT / "outputs" / "crossmatch" / "dd1_support_cache" / "allwise.csv")
aw = aw.set_index("designation").reindex(s.designation)

w1w4 = (aw.W1mag - aw.W4mag).values
sig = w1w4 / np.hypot(aw.e_W1mag, aw.e_W4mag).values
snr = 1.0857 / aw.e_W4mag.values
ccf = s.ccf.astype(str).where(s.ccf.notna(), "").map(lambda c: c.zfill(4) if c.isdigit() else c)
teff = s.gaia_teff.fillna(s.tic_teff)
s["dd_w1w4"], s["dd_w1w4_sigma"], s["dd_w4snr"] = w1w4.round(3), sig.round(1), snr.round(1)
s["dd_ccf_ok"] = ccf.str.len().eq(4) & ccf.str[3].eq("0")
s["dd_ext_ok"] = s.ex.eq(0)
s["dd_mstar"] = teff < 3900
s["dd_w4rchi2"] = np.nan
s["dd_excess_pass"] = ((w1w4 > 0.25) & (sig > 5) & (snr >= 10) & s.dd_ccf_ok & s.dd_ext_ok
                       & (~s.dd_mstar | (w1w4 > 0.9)))

lab = pd.read_csv(ROOT / "data" / "catalogs" / "DD_1.0ObjectsFromMAST.csv", comment="#", low_memory=False)
fr = ["goodFraction", "multiFraction", "ovalFraction", "emptyFraction", "extendedFraction", "shiftFraction"]
lab["majority"] = lab[fr].idxmax(axis=1).str.replace("Fraction", "")
lab = lab.set_index("designation").reindex(s.designation)
s["dd1_zooniverse_id"] = lab.ZooniverseID.values
s["dd1_good_fraction"] = lab.goodFraction.round(3).values
s["dd1_majority_label"] = lab.majority.values
s["dd1_classifiers"] = lab.classifiers.values
s["dd1_sciteam"] = lab.SciTeamFollowUp.values

ours = s.c_ir_excess.astype(bool)
s["dd_vs_ours"] = np.select([s.dd_excess_pass == ours, s.dd_excess_pass], ["agree", "DD1 only"], "ours only")
yes = s.dd1_sciteam.astype(str).str.upper().str.contains("YES")
bad = (~s.dd_excess_pass) | (s.dd1_good_fraction.notna() & (s.dd1_good_fraction <= 0.5))
s["grade"] = s.tier + np.where(bad, "-", np.where(yes, "+", ""))
s.to_csv(P, index=False)

print(f"DD1 cut passes: {s.dd_excess_pass.sum()} of {len(s)}; labelled: {s.dd1_zooniverse_id.notna().sum()}")
print(pd.crosstab(s.tier, s.grade))
print(pd.crosstab(s.tier, s.dd_vs_ours))
fail = s[~s.dd_excess_pass & s.tier.isin(["A", "B"])]
why = lambda r: ",".join(k for k, v in [("colour", r.dd_w1w4 <= 0.25 or r.dd_w1w4_sigma <= 5), ("w4snr", r.dd_w4snr < 10),
                                       ("ccf", not r.dd_ccf_ok), ("ext", not r.dd_ext_ok)] if v)
print(fail.assign(why=fail.apply(why, axis=1))[["designation", "tier", "dd_w1w4", "dd_w1w4_sigma", "dd_w4snr", "ccf", "ex", "why"]].to_string(index=False))
