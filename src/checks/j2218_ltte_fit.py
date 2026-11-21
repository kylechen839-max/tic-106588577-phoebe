"""J2218: DASCH per-window eclipse timing (bootstrap errors) and light-travel-time (LTTE) fits to the full O-C.
Run from the repo root after j2218_dasch_timing.py inputs exist. Windows with dchi2 < 9 are dropped; DASCH errors floored at 10 min."""
import os; os.chdir(os.path.dirname(os.path.abspath(__file__)) + '/../..')
import numpy as np, json
exec(open('src/checks/j2218_dasch_timing.py').read().split('out={}')[0].replace("print('N'","pass #"))
yr=2000+(t+2457000-2451544.5)/365.25
full=(m.copy(),s.copy(),t.copy())
shs=np.arange(-0.15,0.15,0.0005)
def prof(mm_,ww,ph):
    mm=mm_-np.sum(ww*mm_)/ww.sum(); out=[]
    for sh in shs:
        T=np.interp((ph-sh)%1,TX,TY); Tc=T-np.sum(ww*T)/ww.sum()
        b=np.sum(ww*Tc*mm)/np.sum(ww*Tc*Tc); out.append(b*np.sum(ww*Tc*mm) if b>0 else 0)
    return np.array(out)
rows=[]
edges=[1885,1905,1915,1925,1930,1935,1940,1945,1950,1960,1975,1991]
for lo,hi in zip(edges[:-1],edges[1:]):
    sel=(yr>=lo)&(yr<hi); n=sel.sum()
    if n<40: print(lo,hi,n,'skip'); continue
    mm,ss,tt=[a[sel] for a in full]; ww=1/ss**2
    ph=((tt-T0)/P0)%1; p=prof(mm,ww,ph); i=p.argmax()
    # bootstrap the timing
    rng=np.random.default_rng(1); bs=[]
    for k in range(150):
        j=rng.integers(0,n,n); pj=prof(mm[j],ww[j],ph[j]); bs.append(shs[pj.argmax()])
    bs=np.array(bs)*P0*1440
    lo16,hi84=np.percentile(bs,[16,84])
    rows.append(dict(yr=float(yr[sel].mean()),btjd=float(tt.mean()),n=int(n),dchi2=float(p[i]),oc=float(shs[i]*P0*1440),lo=float(lo16),hi=float(hi84)))
    print(f'{lo}-{hi}: n={n:4d} <yr>={yr[sel].mean():.1f} dchi2={p[i]:5.1f} O-C={shs[i]*P0*1440:+6.0f} min  boot 16-84% [{lo16:+.0f},{hi84:+.0f}]')
json.dump(rows,open('outputs/candidates/J221843.61+544715.0/timing/J221843.61+544715.0_dasch_decades.json','w'),indent=1)

# ---- LTTE fits ----

from scipy.optimize import least_squares
P0=1.6692583691009597
yr=lambda bt: 2000+(np.asarray(bt)+2457000-2451544.5)/365.25
D=json.load(open('outputs/candidates/J221843.61+544715.0/timing/J221843.61+544715.0_dasch_decades.json'))
pts=json.load(open('outputs/candidates/J221843.61+544715.0/timing/J221843.61+544715.0_oc_points.json'))['points']
T=[];Y=[];E=[];L=[]
for d in D:
    if d['dchi2']>=9:
        T.append(d['btjd']);Y.append(d['oc']);E.append(max(0.5*(d['hi']-d['lo']),10));L.append('DASCH %.0f'%d['yr'])
T.append(581);Y.append(50.5);E.append(5.8);L.append('ASAS-SN')
for p in pts: T.append(p['btjd']);Y.append(p['oc_min']);E.append(p['err_min']);L.append(p['name'])
T=np.array(T);Y=np.array(Y);E=np.array(E)
tyr=yr(T)
def lin(p): return p[0]+p[1]*(tyr-2000)
def ltte_circ(p): return lin(p)+p[2]*np.sin(2*np.pi*(tyr-p[4])/p[3])
def kepler(M,e):
    Ea=M.copy()
    for _ in range(50): Ea=Ea-(Ea-e*np.sin(Ea)-M)/(1-e*np.cos(Ea))
    return Ea
def ltte_ecc(p):
    c0,c1,A,P3,tp,e,w=p
    M=2*np.pi*(tyr-tp)/P3; Ea=kepler(np.mod(M,2*np.pi),e)
    nu=2*np.arctan2(np.sqrt(1+e)*np.sin(Ea/2),np.sqrt(1-e)*np.cos(Ea/2))
    r=(1-e**2)/(1+e*np.cos(nu))
    return lin(p)+A*(r*np.sin(nu+w)+e*np.sin(w))
res={}
r=least_squares(lambda p:(lin(p)-Y)/E,[0,0]); res['linear']=(np.sum(r.fun**2),2,r.x)
best=None
for P3 in np.arange(15,200,1.0):
    for ph in np.linspace(0,P3,12,endpoint=False):
        r=least_squares(lambda p:(ltte_circ(p)-Y)/E,[0,0,40,P3,1950+ph],bounds=([-500,-50,0,10,1800],[500,50,500,300,2200]))
        c=np.sum(r.fun**2)
        if best is None or c<best[0]: best=(c,5,r.x)
res['circular LTTE']=best
bestE=None
rng=np.random.default_rng(0)
for k in range(3000):
    x0=[0,0,rng.uniform(20,200),rng.uniform(20,200),rng.uniform(1850,2050),rng.uniform(0,0.8),rng.uniform(0,2*np.pi)]
    try: r=least_squares(lambda p:(ltte_ecc(p)-Y)/E,x0,bounds=([-500,-50,0,10,1700,0,0],[500,50,500,300,2300,0.95,2*np.pi]))
    except Exception: continue
    c=np.sum(r.fun**2)
    if bestE is None or c<bestE[0]: bestE=(c,7,r.x)
res['eccentric LTTE']=bestE
n=len(Y)
for k,(c,npar,x) in res.items():
    print(f'{k:15s} chi2={c:7.1f} dof={n-npar:2d}  BIC={c+npar*np.log(n):6.1f}  params={np.round(x,3).tolist()}')
json.dump(dict(T=T.tolist(),Y=Y.tolist(),E=E.tolist(),L=L,fits={k:dict(chi2=v[0],npar=v[1],x=list(map(float,v[2]))) for k,v in res.items()}),open('outputs/candidates/J221843.61+544715.0/timing/J221843.61+544715.0_ltte_fit.json','w'),indent=1)
# derived for eccentric
c,npar,x=bestE; A=x[2]/1440*173.14 # AU light-days: 1 day = 173.14 AU
a12=A/np.sqrt(1-(x[5]*np.cos(x[6]))**2)
print('ecc: a12 sin i ~ %.1f AU (A=%.0f min), P3=%.0f yr, e=%.2f'%(a12,x[2],x[3],x[5]))
for nm,(c,npar,x) in (('circ',res['circular LTTE']),('ecc',bestE)):
    A=x[2]/1440*173.14 if nm=='circ' else a12; P3=x[3]; f=A**3/P3**2
    M12=3.1; m3=0.1
    for _ in range(200): m3=(f*(M12+m3)**2)**(1/3)
    print(f'{nm}: f(m3)={f:.3f} Msun -> m3 sin i >= {m3:.2f} Msun (M12 3.1)')
