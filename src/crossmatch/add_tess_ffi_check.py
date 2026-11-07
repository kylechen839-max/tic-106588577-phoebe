"""Join the tier-B TESS FFI eclipse check (outputs/candidates/tierB_tess/) into the scored DD1 table.

Adds tess_bls_period (filling period_d where the catalogue has none), tess_ffi_class (one letter: E detached eclipses, C close-binary-like, P eclipses but period needs
work, V variable not eclipse-like, N no signal, D dipper-like) and tess_ffi_note, and appends
" (TESS)" to the grade of objects whose FFI light curve shows eclipses (E or C).
Run after add_dd1_criteria.py (which rewrites the grade), then make_tables.py.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / "outputs" / "crossmatch" / "dd1_candidates_scored.csv"
chk = pd.read_csv(ROOT / "outputs" / "candidates" / "tierB_tess" / "tierB_tess_eclipse_check.csv").set_index("designation")
s = pd.read_csv(P, low_memory=False)
eye = chk.eye_class.reindex(s.designation).values
s["tess_ffi_class"] = pd.Series(eye).str[0].values
s["tess_ffi_note"] = pd.Series(eye).str[2:].values
# BLS periods for objects with no catalogue period (src/checks/bls_ffi.py); kept only where the fold showed eclipses
bls = pd.read_csv(ROOT / "outputs" / "candidates" / "tierB_tess" / "noperiod" / "bls_results.csv").set_index("designation")
s["tess_bls_period"] = bls.P.where(bls.eye_class.str[0].isin(["E", "C"])).reindex(s.designation).values
s["period_d"] = s.period_d.fillna(pd.Series(s.tess_bls_period.values, index=s.index))
base = s.grade.astype(str).str.replace(" (TESS)", "", regex=False)
s["grade"] = base.where(~s.tess_ffi_class.isin(["E", "C"]), base + " (TESS)")
s.to_csv(P, index=False)
print(pd.crosstab(s.tier, s.tess_ffi_class.fillna("-")))
print(s.grade.value_counts().to_string())
