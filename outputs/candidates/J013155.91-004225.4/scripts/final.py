import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
N="J013155.91-004225.4"; d=pd.read_csv(f"lc_out/{N}_tess_all_with_flags.csv")
d["use"]=(~d.clipped)&(d.sector!=97)
d.to_csv(f"lc_out/{N}_tess_all_with_flags.csv",index=False)
c=d[d.use].drop(columns=["clipped","use"]); c.to_csv(f"lc_out/{N}_tess_clean.csv",index=False)
Ps=[0.48803,0.48768,0.48773,0.48814,0.48751]; P=np.mean(Ps); print("mean P",P,"std",np.std(Ps,ddof=1))
secs=sorted(c.sector.unique()); fig,axs=plt.subplots(1,len(secs),figsize=(15,3.2),sharey=True)
for ax,s in zip(axs,secs):
    a=c[c.sector==s]; f,p=LombScargle(a.time_btjd,a.flux).autopower(minimum_frequency=0.5,maximum_frequency=6,samples_per_peak=50); Ps_=1/f[np.argmax(p)]
    ph=((a.time_btjd-a.time_btjd.min())/Ps_)%1; bins=np.linspace(0,1,51); m=pd.Series(a.flux.values).groupby(np.digitize(ph,bins)).median(); ok=m.index<=50
    ax.plot(ph,a.flux,".",ms=.4,c="0.75"); ax.plot(0.5*(bins[1:]+bins[:-1])[m.index.values[ok]-1],m[ok],"r-"); ax.set_ylim(0.985,1.015); ax.set_title(f"S{s} {a.author.iloc[0]}  P={Ps_:.4f} d",fontsize=8); ax.set_xlabel("phase",fontsize=7)
fig.suptitle(f"{N} (TIC 248986951): each sector folded at its own Lomb-Scargle period; red = 50-bin median",fontsize=9); fig.tight_layout(); fig.savefig(f"lc_out/{N}_per_sector_fold_LS.png",dpi=105)
print("clean points",len(c), "sectors",secs)
