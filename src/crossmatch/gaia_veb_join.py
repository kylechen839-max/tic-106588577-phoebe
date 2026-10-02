"""Join DD x Gaia-DR3 ECL matches to gaiadr3.vari_eclipsing_binary (Mowlavi+2023) via the Gaia archive.

Adds the Gaia EB-pipeline period, model type, eclipse depths and secondary phase, and checks the
Gaia period against the other surveys. Output: outputs/crossmatch/dd1_gaia_veb.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from astroquery.gaia import Gaia

OUT = Path(__file__).resolve().parents[2] / "outputs" / "crossmatch"


def main():
    long = pd.read_csv(OUT / "dd1_eb_matches_long.csv")
    g = long[long.survey == "gaia_dr3"][["designation", "survey_id"]].drop_duplicates()
    g["source_id"] = g.survey_id.astype("int64")
    up = Table.from_pandas(g[["source_id"]])
    q = ("SELECT u.source_id, v.frequency, v.model_type, v.num_model_parameters, v.global_ranking, "
         "v.derived_primary_ecl_phase, v.derived_primary_ecl_depth, v.derived_secondary_ecl_phase, "
         "v.derived_secondary_ecl_depth, v.derived_primary_ecl_duration, v.derived_secondary_ecl_duration "
         "FROM tap_upload.ids AS u LEFT JOIN gaiadr3.vari_eclipsing_binary AS v ON u.source_id = v.source_id")
    job = Gaia.launch_job_async(q, upload_resource=up, upload_table_name="ids")
    r = job.get_results().to_pandas()
    r.columns = [c.lower() for c in r.columns]
    df = g.merge(r, on="source_id", how="left")
    df["gaia_veb_period"] = 1 / df["frequency"]
    df["in_gaia_veb"] = df["frequency"].notna()
    # compare with the other surveys' periods (allowing x2 / x1/2)
    other = long[(long.survey != "gaia_dr3") & long.period_d.notna()].groupby("designation").period_d.median()
    df["other_period"] = df.designation.map(other)
    ratio = df.gaia_veb_period / df.other_period
    df["period_agrees"] = np.where(df.other_period.notna() & df.gaia_veb_period.notna(),
                                   np.any([np.abs(ratio / f - 1) < 0.01 for f in (0.5, 1, 2)], axis=0), np.nan)
    df.to_csv(OUT / "dd1_gaia_veb.csv", index=False)
    print(f"{len(df)} Gaia ECL matches; in vari_eclipsing_binary: {df.in_gaia_veb.sum()}")
    both = df.period_agrees.notna()
    print(f"with another survey period: {both.sum()}; Gaia period agrees: {int(np.nansum(df.period_agrees))}")
    sc = pd.read_csv(OUT / "dd1_candidates_scored.csv")
    a = sc[sc.tier.isin(["A", "B"])].designation
    print(df[df.designation.isin(a)][["designation", "in_gaia_veb", "gaia_veb_period", "other_period",
                                       "period_agrees", "model_type"]].to_string())


if __name__ == "__main__":
    main()
