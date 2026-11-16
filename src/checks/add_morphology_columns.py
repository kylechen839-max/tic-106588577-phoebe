"""Add the tier A W4 morphology results (tierA_w4_shape.csv, tierA_nebula_model.csv) to the candidate tables.

New columns (blank for stars outside tier A):
  kl14_pass          Koenig & Leisawitz 2014 W4 cut (w4rchi2 < 0.2 w4snr - 2)
  W4_fwhm20          Gaussian + plane FWHM within 20" of the star, arcsec
  W4_fwhm20_field95  95th percentile of the same fit on field point sources in the tile
  W4_fwhm20_nebsub   the same fit after subtracting the method-B nebula model (W3 template + plane)
  W4_ratio_nebB, W4_sigma_nebB   W4 / binary photosphere and excess sigma with the method-B nebula model
  W4_ratio_quad, W4_sigma_quad   the same with a quadratic background within 55"
  w4_morphology      point-like / extended (nebula-subtracted width vs field 95th percentile),
                     no compact peak after nebula removal, or no field comparison
Writes outputs/crossmatch/candidates_all_in_one.csv and the Candidates sheet of DD1_EB_candidates.xlsx in place.
Usage: .venv-tess/bin/python src/checks/add_morphology_columns.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "outputs/crossmatch"
COLS = ["kl14_pass", "W4_fwhm20", "W4_fwhm20_field95", "W4_fwhm20_nebsub", "W4_ratio_nebB", "W4_sigma_nebB",
        "W4_ratio_quad", "W4_sigma_quad", "w4_morphology"]


def table():
    sh = pd.read_csv(X / "tierA_w4_shape/tierA_w4_shape.csv").set_index("designation")
    nb = pd.read_csv(X / "tierA_w4_shape/tierA_nebula_model.csv").set_index("designation")
    t = pd.DataFrame(index=sh.index)
    t["kl14_pass"] = sh["kl14_pass"]
    t["W4_fwhm20"] = sh["W4_fwhm20"]
    t["W4_fwhm20_field95"] = sh["field_fwhm20_p95"]
    t = t.join(nb[["fwhm20_nebsub", "ratio_B", "sigma_B", "ratio_a3", "sigma_a3"]].rename(columns={
        "fwhm20_nebsub": "W4_fwhm20_nebsub", "ratio_B": "W4_ratio_nebB", "sigma_B": "W4_sigma_nebB",
        "ratio_a3": "W4_ratio_quad", "sigma_a3": "W4_sigma_quad"}))

    def verdict(r):
        w = r.W4_fwhm20_nebsub
        if not np.isfinite(w):
            return "no compact peak after nebula removal"
        if not np.isfinite(r.W4_fwhm20_field95):
            return "no field comparison"
        return "point-like" if w <= r.W4_fwhm20_field95 else "extended"
    t["w4_morphology"] = t.apply(verdict, axis=1)
    return t[COLS[:-1] + ["w4_morphology"]]


def main():
    t = table()
    c = pd.read_csv(X / "candidates_all_in_one.csv", dtype=str, keep_default_na=False)   # keep existing values verbatim
    c = c.drop(columns=[k for k in COLS if k in c]).merge(t, left_on="designation", right_index=True, how="left")
    c.to_csv(X / "candidates_all_in_one.csv", index=False)

    p = X / "DD1_EB_candidates.xlsx"
    wb = openpyxl.load_workbook(p); ws = wb["Candidates"]
    head = [x.value for x in ws[1]]
    for k in COLS:                      # append missing headers, copying the last header's style
        if k not in head:
            src = ws.cell(1, len(head)); new = ws.cell(1, len(head) + 1, k)
            new._style = src._style; head.append(k)
    di = head.index("designation") + 1
    for row in range(2, ws.max_row + 1):
        d = ws.cell(row, di).value
        for k in COLS:
            v = t.at[d, k] if d in t.index else None
            if isinstance(v, (float, np.floating)) and not np.isfinite(v):
                v = None
            elif isinstance(v, (np.bool_, bool)):
                v = bool(v)
            elif isinstance(v, np.floating):
                v = float(v)
            ws.cell(row, head.index(k) + 1, v)
    wb.save(p)
    print(t.to_string())


if __name__ == "__main__":
    main()
