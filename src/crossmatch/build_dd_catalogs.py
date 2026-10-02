"""Build position tables for the Disk Detective 1.0 and 2.x target lists.

DD1.0: the WZ_subjects_*.txt IPAC tables (AllWISE + 2MASS photometry).
DD2.x: DDv2.1_Wise_Ids.txt (AllWISE designations only); positions are
parsed from the designation and photometry is fetched later via XMatch.
"""
from pathlib import Path
import re

import numpy as np
import pandas as pd
from astropy.io import ascii

ROOT = Path(__file__).resolve().parents[2]
CAT = ROOT / "data" / "catalogs"
OUT = ROOT / "outputs" / "crossmatch"

PHOT = ["w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro", "w3mpro", "w3sigmpro",
        "w4mpro", "w4sigmpro", "j_m_2mass", "j_msig_2mass", "h_m_2mass",
        "h_msig_2mass", "k_m_2mass", "k_msig_2mass"]


def designation_to_radec(name):
    m = re.match(r"J(\d\d)(\d\d)(\d\d\.\d+)([+-])(\d\d)(\d\d)(\d\d\.\d+)", name.strip())
    if not m:
        return np.nan, np.nan
    h, mi, s, sg, d, dm, ds = m.groups()
    ra = 15.0 * (int(h) + int(mi) / 60 + float(s) / 3600)
    dec = int(d) + int(dm) / 60 + float(ds) / 3600
    return ra, -dec if sg == "-" else dec


def build_dd1():
    frames = []
    for f in sorted((CAT / "WZSubs").glob("WZ_subjects_*.txt")):
        if "copy" in f.name:
            continue
        t = ascii.read(f, format="ipac").to_pandas()
        t["source_file"] = f.name
        frames.append(t)
    df = pd.concat(frames, ignore_index=True)
    df["designation"] = df["designation"].str.strip()
    df = df.drop_duplicates("designation").reset_index(drop=True)
    return df


def build_dd2():
    names = [l.strip() for l in (CAT / "DDv2.1_Wise_Ids.txt").read_text().splitlines() if l.strip()]
    names = list(dict.fromkeys(names))
    radec = [designation_to_radec(n) for n in names]
    return pd.DataFrame({"designation": names,
                         "ra": [r for r, _ in radec], "dec": [d for _, d in radec]})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dd1 = build_dd1()
    dd1.to_csv(OUT / "dd1_targets.csv", index=False)
    dd2 = build_dd2()
    dd2.to_csv(OUT / "dd2_targets.csv", index=False)
    prev = set(l.strip() for l in (CAT / "EB_WiseID_previous.txt").read_text().splitlines() if l.strip())
    print(f"DD1 rows: {len(dd1)}  DD2 rows: {len(dd2)}  previous EB list: {len(prev)}")
    print(f"previous EB list in DD1: {len(prev & set(dd1.designation))}, in DD2: {len(prev & set(dd2.designation))}")
    print(f"DD2 also in DD1: {len(set(dd2.designation) & set(dd1.designation))}")


if __name__ == "__main__":
    main()
