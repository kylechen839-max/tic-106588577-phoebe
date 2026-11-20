"""Build the AllWISE + 2MASS cache (outputs/crossmatch/dd1_support_cache/allwise.csv) from the DD1.0 WZ_subjects tables.

The cloud container cannot reach VizieR, but the WZ_subjects_*.txt IPAC tables (data/catalogs/WZSubs/, from
~/Downloads/WZSubs.zip) already carry AllWISE W1-W4 and 2MASS JHK with errors for every DD1.0 subject.
Writes only the cross-matched designations in outputs/crossmatch/dd1_candidates_scored.csv, angDist = 0.

Usage: .venv-tess/bin/python src/crossmatch/build_allwise_cache_from_wzsubs.py [WZSubs dir]
"""
import sys, glob, pandas as pd
from astropy.io import ascii
wz = sys.argv[1] if len(sys.argv) > 1 else 'data/catalogs/WZSubs'
want = set(pd.read_csv('outputs/crossmatch/dd1_candidates_scored.csv').designation)
rows = []
for f in sorted(glob.glob(f'{wz}/WZ_subjects_*.txt')):
    t = ascii.read(f, format='ipac').to_pandas()
    rows.append(t[t.designation.isin(want)])
t = pd.concat(rows).drop_duplicates('designation')
ren = {'j_m_2mass':'Jmag','j_msig_2mass':'e_Jmag','h_m_2mass':'Hmag','h_msig_2mass':'e_Hmag','k_m_2mass':'Kmag',
       'k_msig_2mass':'e_Kmag'}
for i in range(1,5): ren[f'w{i}mpro'] = f'W{i}mag'; ren[f'w{i}sigmpro'] = f'e_W{i}mag'
out = t.rename(columns=ren); out['angDist'] = 0.0
cols = ['designation','angDist','Jmag','e_Jmag','Hmag','e_Hmag','Kmag','e_Kmag','W1mag','e_W1mag','W2mag','e_W2mag',
        'W3mag','e_W3mag','W4mag','e_W4mag']
out[cols].to_csv('outputs/crossmatch/dd1_support_cache/allwise.csv', index=False)
print(f'{len(out)} of {len(want)} cross-matched designations written')
