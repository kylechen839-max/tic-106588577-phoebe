"""Join Disk Detective 1.0 volunteer labels and the DD1.0 excess criteria onto the EB cross-matches.

Input: data/catalogs/DD_1.0ObjectsFromMAST.csv (DD1.0 MAST table, shared by Alissa Bans 2025-08-14;
~30.7k subjects, not the full DD1.0 list) and outputs/crossmatch/dd1_candidates_scored.csv.
Output: outputs/crossmatch/dd1_labels_join.csv (one row per cross-match found in the MAST table).

DD1.0 criteria (Kuchner+2016, Silverberg+2018) computable from this table:
  [W1]-[W4] > 0.25 and > 5 sigma of the colour error; W4 S/N >= 10 (S/N ~ 1.0857/w4sigmpro);
  [W1]-[W4] > 0.9 if the star is an M dwarf (Tstar < 3900 K). w4rchi2/cc_flags/ext_flg are not in the table.
Volunteer label: majority "good" = goodFraction > 0.5.

Usage: .venv-tess/bin/python src/crossmatch/join_dd1_labels.py
"""
import numpy as np, pandas as pd
d = pd.read_csv('data/catalogs/DD_1.0ObjectsFromMAST.csv', comment='#', low_memory=False)
s = pd.read_csv('outputs/crossmatch/dd1_candidates_scored.csv')
keep = ['designation','ZooniverseID','SciTeamFollowUp','SciTeamComment','flags','classifiers','state',
        'goodFraction','multiFraction','ovalFraction','emptyFraction','extendedFraction','shiftFraction',
        'Tstar','Tdisk','Lir_Lstar','w1mpro','w1sigmpro','w4mpro','w4sigmpro']
m = s[['designation','tier']].merge(d[keep], on='designation', how='inner')
m['w1_w4'] = m.w1mpro - m.w4mpro
m['w1_w4_sig'] = m.w1_w4 / np.hypot(m.w1sigmpro, m.w4sigmpro)
m['w4_snr'] = 1.0857 / m.w4sigmpro
mdwarf = m.Tstar < 3900
m['dd1_excess_pass'] = (m.w1_w4 > 0.25) & (m.w1_w4_sig > 5) & (m.w4_snr >= 10) & (~mdwarf | (m.w1_w4 > 0.9))
fr = ['goodFraction','multiFraction','ovalFraction','emptyFraction','extendedFraction','shiftFraction']
m['majority_label'] = m[fr].idxmax(axis=1).str.replace('Fraction','')
m['volunteer_good'] = m.goodFraction > 0.5
m.sort_values(['tier','designation']).to_csv('outputs/crossmatch/dd1_labels_join.csv', index=False)
print(f'{len(m)} of {len(s)} cross-matches are in the DD1.0 MAST table')
print(m.groupby('tier')[['dd1_excess_pass','volunteer_good']].sum().assign(n=m.groupby('tier').size()))
print(m.SciTeamFollowUp.value_counts().to_string())
