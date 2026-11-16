"""J2218: time Gaia DR3 G epoch photometry and ZTF i against the TESS template, then fit linear and
quadratic O-C models to TESS (eclipse_split_timing.py output) + Gaia. Run from the repo root.
Data: Gaia DR3 epoch photometry (DataLink, DATAMODEL_GAIA) and IRSA ZTF light curves, saved in
outputs/candidates/J221843.61+544715.0/timing/. Gaia times are BJD(TCB)-2455197.5; BTJD = t - 1802.5."""
import numpy as np, json, csv, ast, sys
R='outputs/candidates/J221843.61+544715.0/ffi/'; T='outputs/candidates/J221843.61+544715.0/timing/'
d=np.loadtxt(R+'J221843.61+544715.0_raw.txt'); eph=json.load(open(R+'J221843.61+544715.0_ephemeris.json'))
P,t0=eph['period'],eph['t0']
ref=np.isin(d[:,3],[56,57,76,77,83,84])
ph=((d[ref,0]-t0)/P)%1
e=np.arange(0,1.0001,0.002); x=0.5*(e[1:]+e[:-1]); idx=np.digitize(ph,e)-1
y=np.array([np.median(d[ref,1][idx==i]) for i in range(len(x))])
X=np.r_[x-1,x,x+1]; Y=np.r_[y,y,y]
def fit(t,f,s,label,shifts=np.arange(-60,60.01,0.25)):
    out=[]
    for dt in shifts:
        m=np.interp(((t-t0)/P - dt/1440/P)%1, X, Y)
        A=np.vstack([np.ones_like(m),m]).T
        c,*_=np.linalg.lstsq(A/s[:,None],f/s,rcond=None); r=(f-A@c)/s; out.append(r@r)
    out=np.array(out); i=out.argmin(); red=out[i]/(len(t)-3); sc=max(red,1)
    w=shifts[(out-out[i])/sc<=1]
    print(f'{label}: n={len(t)} mean BTJD {t.mean():.0f} ({2457000+t.mean():.0f} JD)  dt = {shifts[i]:+.1f} min  1sig [{w.min():+.1f},{w.max():+.1f}]  chi2_red {red:.2f}')
    return dict(n=len(t),t_mean=float(t.mean()),dt=float(shifts[i]),lo=float(w.min()),hi=float(w.max()),chi2_red=float(red))
res={}
# Gaia
rows=list(csv.DictReader(open(T+'gaia_dr3_2005076418115691648_epoch_phot.csv')))
r=rows[0]
arr=lambda k: np.array([np.nan if v.strip() in ('null','') else float(v) for v in r[k].strip('()').split(',')])
gt=arr('g_transit_time')-1802.5; gf=arr('g_transit_flux'); ge=arr('g_transit_flux_error')
rej=np.array([v.strip()=='true' for v in r['variability_flag_g_reject'].strip('()').split(',')])
ok=np.isfinite(gt)&np.isfinite(gf)&~rej
print('Gaia span BTJD',gt[ok].min(),gt[ok].max())
gf=gf[ok]/np.median(gf[ok]); ge=np.maximum(ge[ok]/np.median(arr('g_transit_flux')[ok]),0.003)
res['gaia_G']=fit(gt[ok],gf,ge,'Gaia G')
# ZTF
z=list(csv.DictReader(open(T+'ztf_lc.csv')))
zt=np.array([float(a['hjd']) for a in z])-2457000; zm=np.array([float(a['mag']) for a in z]); ze=np.array([float(a['magerr']) for a in z]); zf=np.array([int(a['catflags']) for a in z])
ok=zf!=32768
fl=10**(-0.4*(zm[ok]-np.median(zm[ok]))); fe=np.maximum(fl*ze[ok]*0.921,0.005)
res['ztf_i']=fit(zt[ok],fl,fe,'ZTF i')
json.dump(res,open(T+'J221843.61+544715.0_external_timing.json','w'),indent=1)

# ---- O-C fit ----
R='outputs/candidates/J221843.61+544715.0/'
st=json.load(open(R+'J221843.61+544715.0_split_timing.json')); ex=res
P,t0=st['period'],st['t0']
rows=[]
for s,v in st['sectors'].items():
    # combine primary and secondary (weighted)
    w=np.array([1/v['primary']['err_min']**2,1/v['secondary']['err_min']**2]); y=np.array([v['primary']['dt_min'],v['secondary']['dt_min']])
    rows.append(('TESS S'+s,v['t_mid'],(w@y)/w.sum(),w.sum()**-0.5))
g=ex['gaia_G']; rows.append(('Gaia G',g['t_mean'],g['dt'],0.5*(g['hi']-g['lo'])))
names=[r[0] for r in rows]; t=np.array([r[1] for r in rows]); y=np.array([r[2] for r in rows]); e=np.array([r[3] for r in rows])
E=np.round((t-t0)/P)
def pfit(m,deg):
    c=np.polyfit(E[m],y[m],deg,w=1/e[m]); r=(y-np.polyval(c,E))/e; return c,r
for lab,m in (('TESS only',np.array(['TESS' in n for n in names])),('TESS+Gaia',np.ones(len(y),bool))):
    for deg in (1,2):
        c,r=pfit(m,deg); dof=m.sum()-deg-1
        print(f'{lab} deg{deg}: chi2 {np.sum(r[m]**2):.1f} dof {dof}  pred Gaia {np.polyval(c,E[-1]):+.1f} min (obs {y[-1]:+.1f}±{e[-1]:.1f})')
    if deg==2:
        q=c[0]  # min per cycle^2
        dPdE=2*q/1440  # days per cycle
        dPdt=dPdE/P*365.25
        print(f'   quadratic c={q:.3e} min/cyc^2 -> dP/dt = {dPdt*86400:.3f} s/yr = {dPdt:.2e} d/yr ; vertex at BTJD {t0+P*(-c[1]/(2*c[0])):.0f}')
for n,a,b,cc in zip(names,t,y,e): print(f'{n:9s} BTJD {a:7.1f}  O-C {b:+6.2f} ± {cc:.2f} min')
json.dump(dict(points=[dict(name=n,btjd=float(a),oc_min=float(b),err_min=float(cc)) for n,a,b,cc in zip(names,t,y,e)]),open(T+'J221843.61+544715.0_oc_points.json','w'),indent=1)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
c,_=pfit(np.ones(len(y),bool),2); cl,_=pfit(np.array(['TESS' in n for n in names]),1)
yr=lambda bt: 2000+(bt+2457000-2451544.5)/365.25
tt=np.linspace(200,3700,400); EE=(tt-t0)/P
fig,ax=plt.subplots(figsize=(7.5,4.2))
ax.plot(yr(tt),np.polyval(c,EE),'k-',lw=1,label=f'quadratic (period increasing, dP/dt ≈ +0.49 s/yr)')
ax.plot(yr(tt),np.polyval(cl,EE),'k:',lw=1,label='best linear ephemeris from TESS alone')
m=np.array(['TESS' in n for n in names])
ax.errorbar(yr(t[m]),y[m],e[m],fmt='o',color='C0',label='TESS FFI sectors (primary+secondary)')
ax.errorbar(yr(t[~m]),y[~m],e[~m],fmt='s',color='C3',label='Gaia DR3 G epoch photometry (2014–2017)')
ax.set_xlabel('Year'); ax.set_ylabel('O − C (minutes)'); ax.axhline(0,color='0.7',lw=0.5)
ax.set_title('J2218 eclipse timing: O − C against P = 1.669258 d'); ax.legend(fontsize=8); fig.tight_layout()
fig.savefig(T+'J221843.61+544715.0_oc_gaia_tess.png',dpi=150)
