"""J2218: DASCH DR7 (Harvard plates, 1888-1990) eclipse search and timing against the TESS template.
Data: api.starglass.cfa.harvard.edu/public/dasch/dr7 querycat + lightcurve (refcat apass, ref 422184361544715),
saved as outputs/candidates/J221843.61+544715.0/timing/dasch_dr7_apass_422184361544715.csv. Run from the repo root.
Cuts: magcal_local_rms<0.25, time_accuracy<0.01 d, limiting_mag_local>13, reject_flag==0, 4-sigma MAD clip."""
import numpy as np, pandas as pd, json
R='outputs/candidates/J221843.61+544715.0/ffi/'
d=np.loadtxt(R+'J221843.61+544715.0_raw.txt'); eph=json.load(open(R+'J221843.61+544715.0_ephemeris.json'))
P0,T0=eph['period'],eph['t0']
ref=np.isin(d[:,3],[56,57,76,77,83,84]); ph=((d[ref,0]-T0)/P0)%1
e=np.arange(0,1.0001,0.005); x=0.5*(e[1:]+e[:-1]); idx=np.digitize(ph,e)-1
y=np.array([np.median(d[ref,1][idx==i]) for i in range(len(x))]); y=y/np.median(y)
tm=-2.5*np.log10(y); TX=np.r_[x-1,x,x+1]; TY=np.r_[tm,tm,tm]
g=pd.read_csv('outputs/candidates/J221843.61+544715.0/timing/dasch_dr7_apass_422184361544715.csv')
g=g[np.isfinite(g.magcal_magdep)&(g.magcal_local_rms<0.25)&(g.time_accuracy_days<0.01)&(g.limiting_mag_local>13.0)&(g.reject_flag==0)]
m=g.magcal_magdep.values; s=np.maximum(g.magcal_local_rms.values,0.08)
med=np.median(m); k=np.abs(m-med)<4*1.4826*np.median(np.abs(m-med)); g=g[k]; m=m[k]; s=s[k]
t=g.date_jd.values-2457000
print('N',len(t),'years',2000+(t.min()+2457000-2451544.5)/365.25,2000+(t.max()+2457000-2451544.5)/365.25,'rms',np.std(m))
w=1/s**2
def score(phase):
    # best a + b*T over phase shifts; returns max delta chi2 and best shift, depth b
    best=(0,0,0)
    mm=m-np.sum(w*m)/w.sum(); chi0=np.sum(w*mm**2)
    for sh in np.arange(0,1,0.01):
        T=np.interp((phase-sh)%1,TX,TY); Tc=T-np.sum(w*T)/w.sum()
        b=np.sum(w*Tc*mm)/np.sum(w*Tc*Tc); dchi=b*np.sum(w*Tc*mm)
        if b>0 and dchi>best[0]: best=(dchi,sh,b)
    return best
out={}
# 1) constant period scan
Ps=P0+np.arange(-6e-4,6e-4,1.5e-6)
res=np.array([score(((t-T0)/P)%1) for P in Ps])
i=res[:,0].argmax(); print('const-P best: P=%.7f dchi2=%.1f shift=%.3f depth=%.3f'%(Ps[i],*res[i]))
srt=np.sort(res[:,0]); print('  dchi2 percentiles 50/99/max:',np.percentile(res[:,0],50),np.percentile(res[:,0],99),srt[-1])
out['const']=dict(P=Ps.tolist(),dchi2=res[:,0].tolist(),shift=res[:,1].tolist(),depth=res[:,2].tolist())
# modern-P only
print('at modern P0: dchi2=%.1f shift=%.3f depth=%.3f'%score(((t-T0)/P0)%1))

# ---- per-era timing and broad period search ----
def fine(phase,ww=None):
    mm=m-np.sum(w*m)/w.sum(); best=(0,0,0); prof=[]
    for sh in np.arange(-0.1,0.1,0.001):
        T=np.interp((phase-sh)%1,TX,TY); Tc=T-np.sum(w*T)/w.sum()
        b=np.sum(w*Tc*mm)/np.sum(w*Tc*Tc); dchi=b*np.sum(w*Tc*mm); prof.append(dchi)
        if dchi>best[0]: best=(dchi,sh,b)
    prof=np.array(prof); shs=np.arange(-0.1,0.1,0.001); ok=shs[prof>=prof.max()-1]
    return best, ok.min(), ok.max()
yr=2000+(t+2457000-2451544.5)/365.25
full=(m.copy(),s.copy(),t.copy(),w.copy())
for lab,sel in (('all',yr>0),('1888-1935',yr<1935),('1935-1950',(yr>=1935)&(yr<1950)),('1950-1990',yr>=1950)):
    m,s,t,w=[a[sel] for a in full]
    (dchi,sh,b),lo,hi=fine(((t-T0)/P0)%1)
    print(f'{lab:10s} n={sel.sum():4d} mean yr {yr[sel].mean():.0f}: dchi2 {dchi:5.1f}  shift {sh*P0*1440:+6.0f} min [{lo*P0*1440:+.0f},{hi*P0*1440:+.0f}]  depth scale {b:.2f}')
m,s,t,w=full
# broad periodogram
Ps=np.linspace(0.5,5,1)
fr=np.arange(1/2.5,1/1.2,2e-6)  # frequency grid
res=[]
mm=m-np.sum(w*m)/w.sum()
for f in fr[::5]:
    best=0
    for sh in np.arange(0,1,0.02):
        T=np.interp((((t-T0)*f)-sh)%1,TX,TY); Tc=T-np.sum(w*T)/w.sum(); b=np.sum(w*Tc*mm)/np.sum(w*Tc*Tc)
        if b>0: best=max(best,b*np.sum(w*Tc*mm))
    res.append(best)
res=np.array(res); fr=fr[::5]; o=np.argsort(res)[::-1][:8]
print('broad search 1.2-2.5 d, top peaks:'); [print(f'  P={1/fr[i]:.6f} dchi2={res[i]:.1f}') for i in o]
print('  median',np.median(res),'99.9pct',np.percentile(res,99.9),'n trials',len(res))
