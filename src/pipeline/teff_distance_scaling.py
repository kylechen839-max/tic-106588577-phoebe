"""Estimate how the photometric/Gaia distance ratio of each PHOEBE result changes with the primary Teff.

Keeps the fitted geometry (R/a, T2/T1, q) and rescales sma with the main-sequence mass for a new Teff,
then scales the binary JHKW1 blackbody photosphere. Needs no AllWISE photometry: the stored
d_phot/d_Gaia ratio scales as sqrt(new_shape/old_shape). An approximation; a refit with --teff is the real test.

Usage: .venv-phoebe/bin/python src/pipeline/teff_distance_scaling.py
"""
import json, glob, math, numpy as np, pandas as pd, sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from run_phoebe_candidate import ms_mass_radius, kepler_sma, planck, LAM
d = pd.read_csv('outputs/crossmatch/dd1_candidates_scored.csv').set_index('designation')
lam = np.array([LAM[b] for b in ('J','H','K','W1')])
rows=[]
for p in sorted(glob.glob('outputs/candidates/*/*_phoebe_result*.json')):
    r = json.load(open(p)); n = r['designation']; bb = r['blackbody_check']
    ratio = bb.get('distance_ratio_phot_over_gaia')
    if not ratio: continue
    row = d.loc[n] if n in d.index else {}
    def shape(T1):
        m1,_ = ms_mass_radius(T1); q = r.get('q', r['q_assumed'])
        a = kepler_sma(r['period'], m1*(1+q))
        R1, R2 = r['r1_frac']*a, r['r2_frac']*a
        return R1**2*planck(T1,lam) + R2**2*planck(T1*r['teffratio'],lam)
    s0 = shape(r['teff1'])
    def newratio(T):
        return ratio*math.sqrt(np.mean(shape(T)/s0))
    # Teff needed for ratio=1
    Ts = np.arange(4000, 40001, 100); rs = np.array([newratio(T) for T in Ts])
    need = Ts[np.argmin(abs(rs-1))]
    g = row.get('gaia_teff', np.nan) if len(row) else np.nan
    rows.append(dict(target=n+p.split('result')[1].replace('.json',''), teff_used=r['teff1'], d_ratio=round(ratio,2),
        gaia_teff=g, ratio_at_gaia_teff=round(newratio(g),2) if g==g else None, teff_for_ratio1=int(need),
        bp_rp=row.get('bp_rp') if len(row) else None, l3=r['l3_frac']))
print(pd.DataFrame(rows).to_markdown(index=False))
