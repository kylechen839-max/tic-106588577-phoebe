import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
d=pd.read_csv("lc_out/J013155.91-004225.4_tess_clean.csv")
q=d[d.sector==97]; print("S97 flux pct 1/50/99:",np.percentile(q.flux,[1,50,99]).round(4),"n",len(q),"rms",(1.4826*np.median(np.abs(q.flux-1))).round(4))
a=d[(d.sector!=97)&(d.flux>0.9)&(d.flux<1.1)]
b=BoxLeastSquares(a.time_btjd.values,a.flux.values,a.flux_err.values)
for lo,hi,n in ((0.2,1.2,400000),):
    pg=b.power(np.linspace(lo,hi,n),[0.015,0.025,0.04,0.06]); k=np.argmax(pg.power)
    print("combined BLS:",pg.period[k],pg.duration[k],pg.depth[k]*1e3,pg.depth_snr[k],pg.transit_time[k])
    top=np.argsort(pg.power)[::-1]; seen=[]
    for i in top:
        if all(abs(pg.period[i]-s)>0.002 for s in seen): seen.append(pg.period[i]); print(f"  peak P={pg.period[i]:.6f} pow={pg.power[i]:.4g} depth={pg.depth[i]*1e3:.2f}ppt snr={pg.depth_snr[i]:.1f}")
        if len(seen)==6: break
P=pg.period[k]; fine=np.linspace(P-2e-5,P+2e-5,4001); pf=b.power(fine,pg.duration[k]); k2=np.argmax(pf.power); P=fine[k2]; T0=pf.transit_time[k2]; D=pg.duration[k]
st=b.compute_stats(P,D,T0); print(f"refined P={P:.7f} T0={T0:.5f} dur={D*24:.2f} h depth={pf.depth[k2]*1e3:.2f}±{pf.depth_err[k2]*1e3:.2f} ppt")
print("odd/even depth ppt:",np.array(st["depth_odd"])*1e3,np.array(st["depth_even"])*1e3,"half-phase depth ppt:",np.array(st["depth_half"])*1e3, "phase-0.5 'secondary':",np.array(st["depth_phased"])*1e3)
for PP,lab in ((P,"P"),(2*P,"2P")):
    ph=((a.time_btjd-T0)/PP+0.25)%1-0.25; bins=np.linspace(-0.25,0.75,201); m=pd.Series(a.flux.values).groupby(np.digitize(ph,bins)).median(); e=pd.Series(a.flux.values).groupby(np.digitize(ph,bins)).agg(lambda x:1.4826*np.median(np.abs(x-np.median(x)))/np.sqrt(len(x)))
    ok=(m.index>0)&(m.index<=200); x=0.5*(bins[1:]+bins[:-1])[m.index.values[ok]-1]
    plt.figure(figsize=(10,3.5)); plt.errorbar(x,m[ok],e[ok],fmt=".",ms=3); plt.title(f"S3+30+42+43+70 folded at {lab}={PP:.7f} d, T0={T0:.4f} BTJD (200 bins)",fontsize=8); plt.xlabel("phase"); plt.tight_layout(); plt.savefig(f"lc_out/J013155.91-004225.4_fold_{lab}.png",dpi=110)
    for s,g in a.groupby("sector"):
        ph2=((g.time_btjd-T0)/P+0.5)%1-0.5; inn=np.abs(ph2)<D/P/2; print(f"  {lab} S{s}: in-eclipse depth {1e3*(np.median(g.flux[~inn])-np.median(g.flux[inn])):.2f} ppt n_in={inn.sum()}") if lab=="P" else None
