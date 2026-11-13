import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd; from astropy.timeseries import LombScargle
d=pd.read_csv("lc_out/J013155.91-004225.4_tess_clean.csv"); d=d[(d.flux>0.9)&(d.flux<1.1)]
for s,a in d.groupby("sector"):
    ls=LombScargle(a.time_btjd,a.flux); f,p=ls.autopower(minimum_frequency=0.1,maximum_frequency=12,samples_per_peak=50)
    o=np.argsort(p)[::-1]; pk=[]
    for i in o:
        if all(abs(f[i]-x)>0.05 for x in pk): pk.append(f[i])
        if len(pk)==3: break
    amp=[np.ptp(ls.model(np.linspace(0,1/x,200)+a.time_btjd.iloc[0],x))*1e3 if False else None for x in pk]
    m=LombScargle(a.time_btjd,a.flux).model(np.linspace(a.time_btjd.min(),a.time_btjd.min()+1/pk[0],300),pk[0])
    print(f"S{s}: top P = "+", ".join(f"{1/x:.5f}" for x in pk)+f" d; FAP(top)={ls.false_alarm_probability(p.max()):.1e}; semi-amp {1e3*np.ptp(m)/2:.2f} ppt")
