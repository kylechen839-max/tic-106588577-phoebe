"""Fetch AllWISE / TIC / Gaia support rows for named targets outside the scored cross-match
(e.g. the earlier J0719 and J1122 targets) so src/pipeline/run_phoebe_candidate.py can run on them.
Writes outputs/crossmatch/extra_support_allwise.csv and outputs/crossmatch/extra_scored.csv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sanity_checks import xmatch, nearest, OUT  # noqa: E402


def main(names):
    dd1 = pd.read_csv(OUT / "dd1_targets.csv")
    t = dd1[dd1.designation.isin(names)][["designation", "ra", "dec"]].reset_index(drop=True)
    aw = xmatch(t, "vizier:II/328/allwise", 2.0)
    aw.to_csv(OUT / "extra_support_allwise.csv", index=False)
    tic = nearest(xmatch(t, "vizier:IV/39/tic82", 2.5))
    gaia = nearest(xmatch(t, "vizier:I/355/gaiadr3", 2.0))
    rows = []
    for _, r in t.iterrows():
        d = r.designation
        row = {"designation": d, "ra": r.ra, "dec": r.dec}
        if d in tic.index:
            row.update(tic=str(tic.loc[d, "TIC"]), tic_teff=tic.loc[d, "Teff"], Tmag=tic.loc[d, "Tmag"],
                       tic_contratio=tic.loc[d, "Rcont"])
        if d in gaia.index:
            row.update(gaia_teff=gaia.loc[d, "Teff"], gaia_dist=gaia.loc[d, "Dist"], plx=gaia.loc[d, "Plx"],
                       ruwe=gaia.loc[d, "RUWE"])
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "extra_scored.csv", index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main(sys.argv[1:] or ["J112238.89-592027.5", "J071951.40-240400.6"])
