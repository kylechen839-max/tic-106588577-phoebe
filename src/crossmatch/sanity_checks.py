"""Sanity checks and candidacy scoring for DD x eclipsing-binary cross-matches.

For every DD object with at least one EB-class survey match this computes:

1. match quality      - nearest separation, number of independent surveys
2. period consistency - survey periods agree (allowing factor-2 aliases)
3. WISE quality       - AllWISE qph / ccf / ex flags for W3 and W4
4. IR excess          - blackbody photosphere fit to J,H,Ks,W1 (same method as
                        src/j1122_full_sed.py); W3/W4 excess significance
5. dust vs free-free  - single-temperature dust blackbody vs free-free
                        power law (alpha in [-0.1, 0.6]) on the W2-W4 excess
6. confusion          - Gaia DR3 neighbours inside the W1 (6") and W4 (12") beams
7. astrometry         - Gaia RUWE, parallax S/N (flag only; binaries inflate RUWE)
8. TESS usability     - TIC contamination ratio, Tmag, number of TESScut sectors
9. SIMBAD object type - flags YSOs / galaxies / known non-EB variables

Outputs: outputs/crossmatch/<list>_candidates_scored.csv and a markdown
summary; per-target SED plots for the top candidates.
"""
from pathlib import Path
import argparse
import json
import math
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.table import Table
from astroquery.xmatch import XMatch

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch"

ZP = {"J": 1594.0, "H": 1024.0, "K": 666.7, "W1": 309.540, "W2": 171.787, "W3": 31.674, "W4": 8.363}
LAM = {"J": 1.235, "H": 1.662, "K": 2.159, "W1": 3.35, "W2": 4.60, "W3": 11.56, "W4": 22.09}
MIN_ERR = 0.03  # mag floor so catalogue errors of 0.01 do not dominate the fit
# SIMBAD main types whose IR excess is not expected to be debris-like
EVOLVED_OR_WIND = {"LongPeriodV*", "LongPeriodV*_Candidate", "Mira", "AGB*", "AGB*_Candidate", "post-AGB*",
                   "post-AGB*_Candidate", "RGB*", "RedSG*", "RedSG*_Candidate", "OH/IR*", "C*", "S*", "WolfRayet*",
                   "WolfRayet*_Candidate", "BlueSG", "BlueSG*", "BlueSG_Candidate", "Be*", "Be*_Candidate",
                   "EmLine*", "PlanetaryNeb", "Symbiotic*", "Symbiotic*_Candidate", "SemiRegV*", "Cepheid",
                   "classicalCep", "RVTauV*", "HighMassXBin", "XrayBin", "CV*", "Nova", "LBV*", "Supergiant"}
YSO_TYPES = {"TTauri*", "TTauri*_Candidate", "OrionV*", "YSO", "YSO_Candidate", "Ae*", "Ae*_Candidate",
             "HerbigHaroObj", "FUOr", "EXOr", "Outflow"}
CONTACT_CLASSES = {"EW", "EC", "CBF", "CBH", "EW:"}
# VSX re-ingests survey classifications; map its IDs back to the parent survey
VSX_PREFIX = [("Gaia DR3", "gaia_dr3"), ("ASASSN-V", "asassn"), ("ASAS J", "asas"), ("ZTF", "ztf_chen2020"),
              ("ATO ", "atlas_heinze2018"), ("TIC ", "tess_ebs"), ("WISE J", "wise_chen2018"),
              ("CSS_", "crts_drake2014"), ("SSS_", "crts_south2017"), ("MLS_", "crts_south2017")]
GOOD_OTYPE_BAD = ("G", "QSO", "AGN", "GiG", "GiC", "ClG", "PN", "HII", "Cl*", "OpC", "GlC", "Rad", "IR")


# ----------------------------------------------------------------------------- helpers
def xmatch(df, cat, radius):
    t = Table.from_pandas(df[["designation", "ra", "dec"]])
    for attempt in range(4):
        try:
            return XMatch.query(cat1=t, cat2=cat, max_distance=radius * u.arcsec,
                                colRA1="ra", colDec1="dec").to_pandas()
        except Exception as e:
            print(f"  xmatch {cat} attempt {attempt + 1} failed: {e}")
            time.sleep(20)
    raise RuntimeError(cat)


def nearest(df):
    return df.sort_values("angDist").drop_duplicates("designation").set_index("designation")


def planck(T, lam_um):
    h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
    nu = c / (np.asarray(lam_um) * 1e-6)
    return 2 * h * nu**3 / c**2 / np.expm1(np.clip(h * nu / (k * T), 1e-8, 700))


def periods_consistent(periods):
    p = np.array([x for x in periods if np.isfinite(x) and x > 0])
    if len(p) < 2:
        return np.nan
    ref = p[0]
    for q in p[1:]:
        r = q / ref
        if not any(abs(r - f) / f < 0.01 for f in (0.5, 1.0, 2.0)):
            return False
    return True


def fit_sed(mags, errs, fixed_teff=None):
    """mags/errs: dict band->value. Returns dict with fit results or None.

    With fixed_teff the photosphere temperature is held at the catalogue value
    (TIC/Gaia) and only the scale is fitted, so the excess does not depend on a
    reddening-biased free blackbody temperature.
    """
    bands = [b for b in ["J", "H", "K", "W1", "W2", "W3", "W4"] if np.isfinite(mags.get(b, np.nan))]
    phot = [b for b in ["J", "H", "K", "W1"] if b in bands]
    if len(phot) < 3:
        return None
    lam = np.array([LAM[b] for b in bands])
    e = np.array([max(errs.get(b, np.nan) if np.isfinite(errs.get(b, np.nan)) else 0.1, MIN_ERR) for b in bands])
    f = np.array([ZP[b] * 10 ** (-0.4 * mags[b]) for b in bands])
    s = f * 0.4 * math.log(10) * e
    ip = np.array([b in phot for b in bands])

    best = None
    grid = [fixed_teff] if fixed_teff else np.linspace(2500, 30000, 1101)
    for T in grid:
        m = planck(T, lam[ip])
        w = s[ip] ** -2
        A = np.sum(w * m * f[ip]) / np.sum(w * m * m)
        chi2 = np.sum(((f[ip] - A * m) / s[ip]) ** 2)
        if best is None or chi2 < best[2]:
            best = (T, A, chi2)
    T, A, chi2_star = best
    star = A * planck(T, lam)
    out = {"bb_teff": T, "bb_chi2_jhkw1": chi2_star, "bands": bands, "lam": lam, "flux": f, "sigma": s,
           "star": star, "bb_scale": A}
    for b in ("W2", "W3", "W4"):
        if b in bands:
            i = bands.index(b)
            out[f"{b}_ratio"] = f[i] / star[i]
            out[f"{b}_sigma"] = (f[i] - star[i]) / s[i]
        else:
            out[f"{b}_ratio"] = np.nan
            out[f"{b}_sigma"] = np.nan

    ex = np.array([b in ("W2", "W3", "W4") for b in bands])
    if ex.sum() >= 2:
        res, sx, lx = f[ex] - star[ex], s[ex], lam[ex]
        w = sx ** -2

        def one(basis):
            a = max(np.sum(w * basis * res) / np.sum(w * basis * basis), 0.0)
            return a, float(np.sum(((res - a * basis) / sx) ** 2))

        dust = min(((Td,) + one(planck(Td, lx)) for Td in np.linspace(40, 1500, 1461)), key=lambda r: r[2])
        ff = min(((al,) + one((1 / lx) ** al) for al in np.linspace(-0.1, 0.6, 71)), key=lambda r: r[2])
        out.update(dust_temp=dust[0], dust_amp=dust[1], dust_chi2=dust[2], ff_alpha=ff[0], ff_chi2=ff[2])
        # fractional IR luminosity: ratio of integrated blackbodies (L ~ A * T^4 / pi... use scale*sigma T^4)
        out["f_ir"] = (dust[1] * dust[0] ** 4) / (A * T ** 4) if A > 0 else np.nan
    else:
        out.update(dust_temp=np.nan, dust_amp=np.nan, dust_chi2=np.nan, ff_alpha=np.nan, ff_chi2=np.nan, f_ir=np.nan)
    return out


def plot_sed(name, fit, path, title_extra=""):
    lam, f, s = fit["lam"], fit["flux"], fit["sigma"]
    g = np.logspace(0, np.log10(30), 300)
    star = fit["bb_scale"] * planck(fit["bb_teff"], g)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.errorbar(lam, f, yerr=s, fmt="o", color="k", label="2MASS + AllWISE")
    for b, x, y in zip(fit["bands"], lam, f):
        ax.annotate(b, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax.plot(g, star, color="#4363d8", label=f"photosphere BB {fit['bb_teff']:.0f} K")
    if np.isfinite(fit.get("dust_temp", np.nan)) and fit["dust_amp"] > 0:
        ax.plot(g, star + fit["dust_amp"] * planck(fit["dust_temp"], g), color="#e6194b",
                label=f"+ dust BB {fit['dust_temp']:.0f} K")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wavelength (micron)"); ax.set_ylabel("Flux density (Jy)")
    ax.set_title(f"{name} {title_extra}", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.25, which="both")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def tess_sectors(ra, dec):
    from astroquery.mast import Tesscut
    from astropy.coordinates import SkyCoord
    try:
        t = Tesscut.get_sectors(coordinates=SkyCoord(ra, dec, unit="deg"))
        return len(t), ",".join(str(x) for x in t["sector"])
    except Exception:
        return np.nan, ""


# ----------------------------------------------------------------------------- main
def score(dd_name, top_plots=40, tess=True):
    long = pd.read_csv(OUT / f"{dd_name}_eb_matches_long.csv")
    targets = pd.read_csv(OUT / f"{dd_name}_targets.csv")
    agg = long.groupby("designation").agg(
        n_surveys=("survey", "nunique"),
        surveys=("survey", lambda s: ";".join(sorted(set(s)))),
        min_sep=("sep_arcsec", "min"),
        periods=("period_d", lambda p: ";".join(f"{x:.6g}" for x in p if np.isfinite(x))),
        classes=("survey_class", lambda c: ";".join(sorted(set(str(x) for x in c)))),
        ids=("survey_id", lambda c: ";".join(sorted(set(str(x) for x in c)))),
    ).reset_index()
    agg["period_consistent"] = [periods_consistent([float(x) for x in p.split(";") if x])
                                for p in agg["periods"].fillna("")]
    agg["period_d"] = [np.nanmedian([float(x) for x in p.split(";") if x]) if p else np.nan
                       for p in agg["periods"].fillna("")]
    def independent(group):
        src = set()
        for survey, sid in zip(group["survey"], group["survey_id"].astype(str)):
            if survey == "vsx":
                parent = next((p for pre, p in VSX_PREFIX if sid.startswith(pre)), "vsx_literature")
                src.add(parent)
            else:
                src.add(survey)
        return len(src)

    agg["n_independent"] = agg["designation"].map(long.groupby("designation").apply(independent))
    cand = agg.merge(targets[["designation", "ra", "dec"]], on="designation", how="left")
    print(f"{dd_name}: {len(cand)} candidates; querying support catalogues")

    cache = OUT / f"{dd_name}_support_cache"
    cache.mkdir(exist_ok=True)

    def cached(name, cat, radius):
        """XMatch support catalogue; only designations not queried before are sent."""
        p = cache / f"{name}.csv"
        done_p = cache / f"{name}_queried.txt"
        old = pd.read_csv(p) if p.exists() else pd.DataFrame(columns=["designation"])
        done = set(done_p.read_text().split()) if done_p.exists() else set(old["designation"])
        todo = cand[~cand["designation"].isin(done)]
        if len(todo):
            new = xmatch(todo, cat, radius)
            old = pd.concat([old, new], ignore_index=True) if len(old) else new
            old.to_csv(p, index=False)
            done |= set(todo["designation"])
            done_p.write_text("\n".join(sorted(done)))
        return old

    aw = nearest(cached("allwise", "vizier:II/328/allwise", 2.0))
    gaia_all = cached("gaia", "vizier:I/355/gaiadr3", 12.0)
    tic = nearest(cached("tic", "vizier:IV/39/tic82", 2.5))
    sim = nearest(cached("simbad", "simbad", 3.0))

    rows = []
    for _, c in cand.iterrows():
        d = c["designation"]
        r = dict(c)
        # --- AllWISE quality + photometry
        if d in aw.index:
            a = aw.loc[d]
            qph, ccf, ex = str(a["qph"]), str(a["ccf"]), a["ex"]
            r.update(allwise_sep=a["angDist"], qph=qph, ccf=ccf, ex=ex)
            r["w3w4_quality_ok"] = (len(qph) == 4 and qph[2] in "AB" and qph[3] in "AB"
                                    and len(ccf) == 4 and ccf[2] in "0d" and ccf[3] in "0d"
                                    and str(ex) in ("0", "0.0"))
            mags = {b: a.get(k, np.nan) for b, k in zip(["J", "H", "K", "W1", "W2", "W3", "W4"],
                                                        ["Jmag", "Hmag", "Kmag", "W1mag", "W2mag", "W3mag", "W4mag"])}
            errs = {b: a.get("e_" + k, np.nan) for b, k in zip(["J", "H", "K", "W1", "W2", "W3", "W4"],
                                                              ["Jmag", "Hmag", "Kmag", "W1mag", "W2mag", "W3mag", "W4mag"])}
            mags = {k: float(v) if pd.notna(v) else np.nan for k, v in mags.items()}
            errs = {k: float(v) if pd.notna(v) else np.nan for k, v in errs.items()}
            fit = fit_sed(mags, errs)
            teff_cat = np.nan
            if d in tic.index and pd.notna(tic.loc[d, "Teff"]):
                teff_cat = float(tic.loc[d, "Teff"])
            fit_fixed = fit_sed(mags, errs, fixed_teff=teff_cat) if np.isfinite(teff_cat) else None
            if fit_fixed:
                r.update(fixed_teff=teff_cat, fixed_chi2_jhkw1=fit_fixed["bb_chi2_jhkw1"],
                         fixed_W3_sigma=fit_fixed["W3_sigma"], fixed_W4_sigma=fit_fixed["W4_sigma"],
                         fixed_W4_ratio=fit_fixed["W4_ratio"], fixed_W2_sigma=fit_fixed["W2_sigma"],
                         fixed_dust_temp=fit_fixed["dust_temp"], fixed_f_ir=fit_fixed["f_ir"])
        else:
            r.update(allwise_sep=np.nan, qph="", ccf="", ex=np.nan, w3w4_quality_ok=False)
            fit = None
        if fit:
            for k in ("bb_teff", "bb_chi2_jhkw1", "W2_ratio", "W2_sigma", "W3_ratio", "W3_sigma", "W4_ratio",
                      "W4_sigma", "dust_temp", "dust_chi2", "ff_alpha", "ff_chi2", "f_ir"):
                r[k] = fit.get(k, np.nan)
        # --- Gaia
        g = gaia_all[gaia_all["designation"] == d].sort_values("angDist")
        if len(g):
            g0 = g.iloc[0]
            gm = g0["Gmag"]
            r.update(gaia_source=str(g0["Source"]), gaia_sep=g0["angDist"], Gmag=gm, bp_rp=g0.get("BP-RP"),
                     plx=g0["Plx"], plx_snr=g0["Plx"] / g0["e_Plx"] if g0["e_Plx"] else np.nan,
                     ruwe=g0["RUWE"], gaia_teff=g0.get("Teff"), gaia_dist=g0.get("Dist"))
            if pd.isna(r["gaia_dist"]) and pd.notna(g0["Plx"]) and g0["Plx"] > 0 and r["plx_snr"] > 5:
                r["gaia_dist"] = 1000.0 / g0["Plx"]  # GSP-Phot distance missing: inverse parallax
            others = g.iloc[1:]
            r["n_gaia_6as_bright"] = int(((others["angDist"] <= 6) & (others["Gmag"] < gm + 2.5)).sum())
            r["n_gaia_12as_bright"] = int(((others["angDist"] <= 12) & (others["Gmag"] < gm + 2.5)).sum())
            # flux fraction of neighbours within the W4 beam
            fl = 10 ** (-0.4 * others.loc[others["angDist"] <= 12, "Gmag"].dropna())
            r["gaia_12as_flux_frac"] = float(fl.sum() / 10 ** (-0.4 * gm)) if np.isfinite(gm) else np.nan
        # --- TIC
        if d in tic.index:
            t = tic.loc[d]
            r.update(tic=str(t["TIC"]), Tmag=t["Tmag"], tic_teff=t["Teff"], tic_rad=t["Rad"],
                     tic_contratio=t["Rcont"], tic_lumclass=t["LClass"])
        # --- SIMBAD
        if d in sim.index:
            s = sim.loc[d]
            r.update(simbad_id=s["main_id"], simbad_otype=s["main_type"] if "main_type" in s else s["otype"],
                     simbad_sptype=s["sp_type"])
        rows.append(r)
    df = pd.DataFrame(rows)

    # ------------------------------------------------------------------ checks
    def col(n, default=np.nan):
        return df[n] if n in df.columns else pd.Series(default, index=df.index)

    chk = pd.DataFrame(index=df.index)
    chk["c_match"] = (col("min_sep") <= 1.5) | (col("n_independent") >= 2)
    chk["c_multi_survey"] = col("n_independent") >= 2
    chk["c_period"] = col("period_consistent").astype(object).where(col("period_consistent").notna(), True).astype(bool)
    classes = col("classes", "").fillna("").astype(str)
    per = col("period_d")
    is_contact = classes.apply(lambda c: any(x.strip() in CONTACT_CLASSES for x in c.split(";")))
    chk["c_period_plausible"] = per.isna() | ((per > 0.15) & (per < 30) & (~is_contact | (per < 1.5)))
    chk["c_wise_quality"] = col("w3w4_quality_ok").fillna(False).astype(bool)
    # an excess band only counts when its AllWISE photometric quality is A or B (not an upper limit)
    q = col("qph", "").fillna("").astype(str)
    w3_ok = q.str.len().eq(4) & q.str[2].isin(["A", "B"])
    w4_ok = q.str.len().eq(4) & q.str[3].isin(["A", "B"])
    free_ex = ((col("W4_sigma") >= 5) & w4_ok) | ((col("W3_sigma") >= 5) & w3_ok)
    fixed_ex = ((col("fixed_W4_sigma") >= 5) & w4_ok) | ((col("fixed_W3_sigma") >= 5) & w3_ok)
    # blackbody check: excess must survive both a free-temperature and a catalogue-Teff photosphere
    chk["c_ir_excess"] = free_ex & (fixed_ex | col("fixed_W4_sigma").isna())
    # photometry check: catalogue-Teff photosphere describes J,H,Ks,W1 (no strong near-IR / hot-dust excess)
    chk["c_photosphere"] = (col("fixed_chi2_jhkw1") < 60) | (col("fixed_chi2_jhkw1").isna() & (col("bb_chi2_jhkw1") < 60))
    chk["c_dust_over_ff"] = col("dust_chi2") < col("ff_chi2")
    chk["c_unconfused"] = (col("n_gaia_6as_bright").fillna(0) == 0) & (col("gaia_12as_flux_frac").fillna(0) < 0.25)
    chk["c_tess_usable"] = (col("tic_contratio").fillna(0) < 0.3) & (col("Tmag").fillna(99) < 14.5)
    otype = col("simbad_otype", "").fillna("").astype(str)
    chk["c_not_extragalactic"] = ~otype.isin(GOOD_OTYPE_BAD)
    chk["c_not_evolved_or_wind"] = ~otype.isin(EVOLVED_OR_WIND)
    # W4 beam (12") spans >0.17 pc beyond 3 kpc: confusion with ISM / clusters dominates
    chk["c_nearby"] = col("gaia_dist").isna() | (col("gaia_dist") < 3000)
    df = pd.concat([df, chk], axis=1)
    df["is_yso"] = otype.isin(YSO_TYPES)
    weights = {"c_match": 1, "c_multi_survey": 1, "c_period": 1, "c_period_plausible": 1, "c_wise_quality": 2,
               "c_ir_excess": 2, "c_photosphere": 1, "c_dust_over_ff": 1, "c_unconfused": 2, "c_tess_usable": 2,
               "c_not_extragalactic": 1, "c_not_evolved_or_wind": 2, "c_nearby": 1}
    # image check (src/crossmatch/wise_cutouts.py): a compact W3 or W4 source (S/N >= 5) within 4" of the
    # star. Catalogue W3/W4 photometry with clean flags can still be dominated by extended nebulosity.
    cut_files = sorted((OUT / "wise_cutouts").glob("wise_centring*.csv"))
    if cut_files:
        cut = pd.concat([pd.read_csv(f) for f in cut_files]).drop_duplicates("designation", keep="last")
        cut = cut.set_index("designation")
        def img_ok(d):
            if d not in cut.index:
                return np.nan
            r = cut.loc[d]
            ok3 = (r["12_snr"] >= 5) and (r["12_offset"] < 4)
            ok4 = (r["22_snr"] >= 5) and (r["22_offset"] < 4)
            return bool(ok3 or ok4)
        df["c_wise_image"] = df["designation"].map(img_ok)
    else:
        df["c_wise_image"] = np.nan
    df["score"] = sum(df[k].astype(int) * w for k, w in weights.items())
    df["max_score"] = sum(weights.values())
    critical = ["c_wise_quality", "c_ir_excess", "c_unconfused", "c_tess_usable", "c_period_plausible",
                "c_not_extragalactic", "c_not_evolved_or_wind", "c_match"]
    crit_ok = df[critical].all(axis=1)
    # A: every critical check, plus the binary-photosphere/photometry checks, plus independent EB confirmation
    image_ok = df["c_wise_image"].astype(object).where(df["c_wise_image"].notna(), True).astype(bool)
    # The automated image S/N can be fooled by bright nebular gradients (J0719 has a compact, centred
    # W3/W4 source but measures S/N~4), so a failed image check only blocks tier A and flags the object
    # for visual review; it does not reject it.
    tier_a = (crit_ok & image_ok & df["c_photosphere"] & df["c_dust_over_ff"] & df["c_multi_survey"]
              & df["c_nearby"] & ~df["is_yso"])
    tier_b = crit_ok & ~tier_a
    df["tier"] = np.where(tier_a, "A", np.where(tier_b, "B", "C"))
    df["fail_reasons"] = [";".join([k[2:] for k in weights if not r[k]] +
                                   (["wise_image"] if r["c_wise_image"] == False else []))  # noqa: E712 (nan-safe)
                          for _, r in df.iterrows()]
    df = df.sort_values(["tier", "score", "n_independent", "W4_sigma"],
                        ascending=[True, False, False, False]).reset_index(drop=True)

    if tess:
        top = df["tier"].isin(["A", "B"])
        print(f"querying TESScut sectors for {top.sum()} tier A/B candidates")
        sec_cache_path = cache / "tess_sectors.json"
        sec_cache = json.loads(sec_cache_path.read_text()) if sec_cache_path.exists() else {}
        sec = []
        for r, t in zip(df.itertuples(), top):
            if not t:
                sec.append((np.nan, ""))
                continue
            if r.designation not in sec_cache:
                sec_cache[r.designation] = tess_sectors(r.ra, r.dec)
            sec.append(tuple(sec_cache[r.designation]))
        sec_cache_path.write_text(json.dumps(sec_cache))
        df["n_tess_sectors"] = [s[0] for s in sec]
        df["tess_sectors"] = [s[1] for s in sec]

    df.to_csv(OUT / f"{dd_name}_candidates_scored.csv", index=False)

    plot_dir = OUT / f"{dd_name}_sed_plots"
    plot_dir.mkdir(exist_ok=True)
    for _, r in df.head(top_plots).iterrows():
        if r["designation"] not in aw.index:
            continue
        a = aw.loc[r["designation"]]
        mags = {b: float(a[k]) if pd.notna(a[k]) else np.nan for b, k in zip(
            ["J", "H", "K", "W1", "W2", "W3", "W4"], ["Jmag", "Hmag", "Kmag", "W1mag", "W2mag", "W3mag", "W4mag"])}
        errs = {b: float(a["e_" + k]) if pd.notna(a["e_" + k]) else np.nan for b, k in zip(
            ["J", "H", "K", "W1", "W2", "W3", "W4"], ["Jmag", "Hmag", "Kmag", "W1mag", "W2mag", "W3mag", "W4mag"])}
        fit = fit_sed(mags, errs)
        if fit:
            plot_sed(r["designation"], fit, plot_dir / f"{r['designation']}_sed.png",
                     f"tier {r['tier']} score {r['score']}/{r['max_score']}")
    print(df["tier"].value_counts().to_dict())
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--lists", nargs="*", default=["dd2", "dd1"])
    ap.add_argument("--no-tess", action="store_true")
    args = ap.parse_args()
    for name in args.lists:
        score(name, tess=not args.no_tess)
