import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
from astroquery.vizier import Vizier; from astropy.coordinates import SkyCoord; import astropy.units as u
v=Vizier(columns=["**"]); r=v.query_region(SkyCoord(22.982958,-0.707056,unit="deg"),radius=5*u.arcsec,catalog="B/vsx/vsx")
print(r[0].to_pandas().iloc[0].dropna().to_string())
d=pd.read_csv("lc_out/J013155.91-004225.4_tess_clean.csv"); P0=0.351176
fig,axs=plt.subplots(2,3,figsize=(13,6.5)); out=[]
for ax,(s,a) in zip(axs.flat,d.groupby("sector")):
    a=a[(a.flux>0.8)&(a.flux<1.2)]
    b=BoxLeastSquares(a.time_btjd.values,a.flux.values,a.flux_err.values)
    pg=b.power(np.exp(np.linspace(np.log(0.2),np.log(6),60000)),[0.02,0.04,0.08])
    k=np.argmax(pg.power); st=b.compute_stats(pg.period[k],pg.duration[k],pg.transit_time[k])
    # power at VSX period
    pv=b.power([P0],[0.02,0.04],objective="snr"); 
    i=np.argmax(pg.depth_snr); out.append((s,pg.period[k],pg.depth[k]*1e3,pg.depth_snr[k],pv.depth[0]*1e3,pv.depth_snr[0]))
    ph=((a.time_btjd-a.time_btjd.min())/P0)%1; bins=np.linspace(0,1,81); m=pd.Series(a.flux.values).groupby(np.digitize(ph,bins)).median()
    ax.plot(ph,a.flux,".",ms=.5,c="0.7"); ax.plot(0.5*(bins[1:]+bins[:-1])[m.index.values[m.index<=80]-1],m[m.index<=80],"r-")
    ax.set_ylim(0.97,1.03); ax.set_title(f"S{s} fold at VSX P; BLS best {pg.period[k]:.4f} d, {pg.depth[k]*1e3:.1f} ppt (SNR {pg.depth_snr[k]:.1f})",fontsize=7)
fig.tight_layout(); fig.savefig("lc_out/J013155.91-004225.4_per_sector_vsx_fold_bls.png",dpi=100)
print(pd.DataFrame(out,columns=["sector","bls_P","bls_depth_ppt","bls_snr","depth_at_vsxP_ppt","snr_at_vsxP"]).to_string(index=False))
