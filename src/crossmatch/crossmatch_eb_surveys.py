"""Cross-match the DD1.0 / DD2.x target lists against eclipsing-binary surveys.

Uses the CDS XMatch service for VizieR tables, and a local astropy match
for ASAS-SN (not served by XMatch). Writes one CSV per survey per DD list and
a merged table of every DD object with at least one EB-class match.
"""
from pathlib import Path
import argparse
import re
import time

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch"
RADIUS = 3.0  # arcsec; separations are kept so tighter cuts can be applied later

VSX_ECLIPSING = {"E", "EA", "EB", "EW", "ED", "EC", "ESD", "E-DO", "ESD/EC"}


def vsx_is_eclipsing(vtype):
    """True if any VSX type token is an eclipsing class (EP = exoplanet transit is excluded)."""
    tokens = re.split(r"[|+/]", str(vtype).strip())
    return any(t.strip().rstrip(":") in VSX_ECLIPSING for t in tokens)


ID_COLS = {"gaia_dr3": "Source", "tess_ebs": "TIC", "tess_oba_ebs": "TIC", "ztf_chen2020": "ID",
           "vsx": "Name", "atlas_heinze2018": "ATOID", "wise_chen2018": "WISE", "crts_drake2014": "CRTS",
           "crts_south2017": "CRTS", "kepler_ebs": "KIC", "asassn": "ASASSN-V", "asas3": "ASAS",
           "nsvs_hoffman2009": "ID", "debcat": "Name", "ogle4_wang2024": "ID", "ogle_lmc_ecl": "Star",
           "tess_10k_kostov2025": "TIC", "tess_obaf_ijspeert2024": "TIC", "tess_ea_shi2022": "TIC",
           "tess_gao2025": "TIC", "css_algol_papageorgiou2018": "CRTS", "lamost_ea_qian2018": "Name",
           "k2_armstrong2015": "EPIC", "ztf_2plus2_vaessen2024": "GaiaDR3", "css_contact_wang2024": "CRTS"}

# name: (vizier table, function(df) -> boolean mask for eclipsing class, period column)
SURVEYS = {
    "gaia_dr3": ("I/358/vclassre", lambda d: d["Class"].astype(str).str.strip() == "ECL", None),
    "tess_ebs": ("J/ApJS/258/16/tess-ebs", lambda d: np.ones(len(d), bool), "Per"),
    "tess_oba_ebs": ("J/A+A/652/A120/eb-cat", lambda d: np.ones(len(d), bool), None),
    "ztf_chen2020": ("J/ApJS/249/18/table2", lambda d: d["Type"].astype(str).str.strip().isin(["EA", "EW"]), "Per"),
    "vsx": ("B/vsx/vsx", lambda d: d["Type"].astype(str).map(vsx_is_eclipsing), "Period"),
    "atlas_heinze2018": ("J/AJ/156/241/table4", lambda d: d["Class"].astype(str).str.strip().isin(["CBF", "CBH", "DBF", "DBH"]), "fp-period"),
    "wise_chen2018": ("J/ApJS/237/28/table2", lambda d: d["Type"].astype(str).str.strip().isin(["EW", "EA"]), "Per"),
    "crts_drake2014": ("J/ApJS/213/9/table3", lambda d: d["Cl"].isin([1, 2, 3]), "Per"),
    "crts_south2017": ("J/MNRAS/469/3688/table4", lambda d: d["Cl"].isin([1, 2, 3]), "Per"),
    "kepler_ebs": ("J/AJ/151/68/table1", lambda d: np.ones(len(d), bool), "Per"),
    "asas3": ("II/264/asas3", lambda d: d["Class"].astype(str).map(vsx_is_eclipsing), "Per"),
    "nsvs_hoffman2009": ("J/AJ/138/466/variables",
                         lambda d: d["Vcl"].astype(str).str.strip().isin(["A", "W", "B"]), "Per"),
    "debcat": ("V/152/debcat", lambda d: np.ones(len(d), bool), "Per"),
    "ogle4_wang2024": ("J/ApJS/275/12/table1", lambda d: np.ones(len(d), bool), "Per"),
    "ogle_lmc_ecl": ("J/AcA/66/421/ecl", lambda d: d["Type"].astype(str).str.strip().isin(["C", "NC"]), "Per"),
    # added 2026-10-04 (recent TESS-FFI, CSS, LAMOST, K2 and ZTF catalogues)
    "tess_10k_kostov2025": ("J/ApJS/279/50/table3", lambda d: np.ones(len(d), bool), "Per"),
    "tess_obaf_ijspeert2024": ("J/A+A/691/A242/obaf-eb1", lambda d: np.ones(len(d), bool), None),
    "tess_ea_shi2022": ("J/ApJS/259/50/table1", lambda d: np.ones(len(d), bool), "Per"),
    "tess_gao2025": ("J/ApJS/276/57/table5", lambda d: d["Type"].astype(str).str.strip().isin(["EA", "EB", "EW"]), "Per"),
    "css_algol_papageorgiou2018": ("J/ApJS/238/4/binaries", lambda d: np.ones(len(d), bool), "Per"),
    "lamost_ea_qian2018": ("J/ApJS/235/5/table2", lambda d: np.ones(len(d), bool), "Per"),
    "k2_armstrong2015": ("J/A+A/579/A19/table3", lambda d: d["Type"].astype(str).str.strip() == "EB", "Per"),
    "ztf_2plus2_vaessen2024": ("J/A+A/682/A164/table1", lambda d: np.ones(len(d), bool), "PerA"),
    "css_contact_wang2024": ("J/ApJS/273/31/table1", lambda d: np.ones(len(d), bool), "Period"),
}


def xmatch(df, table, label):
    t = Table.from_pandas(df[["designation", "ra", "dec"]])
    for attempt in range(4):
        try:
            res = XMatch.query(cat1=t, cat2=f"vizier:{table}", max_distance=RADIUS * u.arcsec,
                               colRA1="ra", colDec1="dec")
            return res.to_pandas()
        except Exception as e:  # network hiccups on large uploads
            print(f"  {label}: attempt {attempt + 1} failed: {e}")
            time.sleep(20)
    raise RuntimeError(f"XMatch failed for {label}")


def asassn_local(df):
    cache = OUT / "asassn_eb_catalog.csv"
    if cache.exists():
        cat = pd.read_csv(cache)
    else:
        v = Vizier(columns=["ASASSN-V", "RAJ2000", "DEJ2000", "Vmag", "Amp", "Per", "Type", "Prob"],
                   row_limit=-1, timeout=600)
        frames = []
        for typ in ["EA", "EB", "EW", "EA:", "EB:", "EW:"]:
            r = v.query_constraints(catalog="II/366/catv2021", Type=f"={typ}")
            if len(r):
                frames.append(r[0].to_pandas())
                print(f"  ASAS-SN type {typ}: {len(r[0])}")
        cat = pd.concat(frames, ignore_index=True)
        cat.to_csv(cache, index=False)
    c1 = SkyCoord(df["ra"].values * u.deg, df["dec"].values * u.deg)
    c2 = SkyCoord(cat["RAJ2000"].values * u.deg, cat["DEJ2000"].values * u.deg)
    idx, sep, _ = c1.match_to_catalog_sky(c2)
    ok = sep.arcsec <= RADIUS
    m = cat.iloc[idx[ok]].reset_index(drop=True)
    m.insert(0, "angDist", sep.arcsec[ok])
    m.insert(0, "designation", df["designation"].values[ok])
    return m


def ijspeert_consensus_period(row, tol=0.01):
    """IJspeert+2024 gives one period per pipeline (SPOC/QLP) and year; XMatch does not return them.
    Keep the value that at least two pipeline-years agree on within 1%, or the only value if there is one."""
    pers = sorted(float(row[c]) for c in row.index if c.startswith("Per") and pd.notna(row[c]) and float(row[c]) > 0)
    if len(pers) == 1:
        return pers[0]
    best = [p for p in pers if sum(abs(q / p - 1) < tol for q in pers) >= 2]
    return float(np.median(best)) if best else np.nan


def add_ijspeert_periods(m):
    """Attach periods fetched by cone search (outputs/crossmatch/ijspeert2024_periods.csv)."""
    path = OUT / "ijspeert2024_periods.csv"
    if not path.exists() or not len(m):
        return m
    per = pd.read_csv(path)
    per["Per"] = per.apply(ijspeert_consensus_period, axis=1)
    return m.merge(per[["designation", "Per"]], on="designation", how="left")


def run(dd_name, df, surveys):
    rows = []
    for name in surveys:
        out = OUT / f"{dd_name}_x_{name}.csv"
        if out.exists():
            m = pd.read_csv(out)
        else:
            print(f"{dd_name} x {name} ...", flush=True)
            if name == "asassn":
                m = asassn_local(df)
            else:
                table, _, _ = SURVEYS[name]
                m = xmatch(df, table, f"{dd_name} x {name}")
            m.to_csv(out, index=False)
        if name == "asassn":
            ecl = np.ones(len(m), bool)
            per = "Per"
            typecol = "Type"
        else:
            _, sel, per = SURVEYS[name]
            ecl = np.asarray(sel(m)) if len(m) else np.zeros(0, bool)
            typecol = next((c for c in ["Class", "Type", "Cl", "Morph"] if c in m.columns), None)
        m = m[ecl]
        if name == "tess_obaf_ijspeert2024":
            m = add_ijspeert_periods(m)
            per = "Per" if "Per" in m.columns else None
        print(f"  {dd_name} x {name}: {len(m)} eclipsing-class matches within {RADIUS}\"")
        for _, r in m.iterrows():
            rows.append({
                "designation": r["designation"], "survey": name,
                "sep_arcsec": r["angDist"],
                "survey_id": str(r.get(ID_COLS[name], "")),
                "period_d": r.get(per, np.nan) if per else np.nan,
                "survey_class": r.get(typecol, "") if typecol else "",
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--surveys", nargs="*", default=list(SURVEYS) + ["asassn"])
    ap.add_argument("--lists", nargs="*", default=["dd2", "dd1"])
    args = ap.parse_args()
    for dd_name in args.lists:
        df = pd.read_csv(OUT / f"{dd_name}_targets.csv")
        res = run(dd_name, df, args.surveys)
        res.to_csv(OUT / f"{dd_name}_eb_matches_long.csv", index=False)
        print(f"{dd_name}: {res.designation.nunique()} unique DD objects with an EB-class match")


if __name__ == "__main__":
    main()
