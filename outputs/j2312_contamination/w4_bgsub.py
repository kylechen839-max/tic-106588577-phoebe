"""J2312: W4 (and W3) excess after removing the extended background halo under the star."""
import warnings; warnings.filterwarnings("ignore")
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits; from astropy.wcs import WCS; from astropy.coordinates import SkyCoord
D=sys.argv[1]; N="J231201.40+532028.6"; c=SkyCoord(348.0058494,53.3413022,unit="deg"); PIX=2.0; RAP=16.5
CAT={"W3":dict(mag=8.062,snr=51.3,ratio_free=1.257643,ratio_fixed=None),
     "W4":dict(mag=6.840,snr=17.5,ratio_free=3.629307,ratio_fixed=3.342305)}
fig,axs=plt.subplots(1,2,figsize=(11,4))
for ax,(b,f) in zip(axs,(("W3","WISE12"),("W4","WISE22"))):
    h=fits.open(f"{D}/{N}_{f}_10arcmin.fits")[0]; img=h.data.astype(float); w=WCS(h.header); x0,y0=w.world_to_pixel(c)
    yy,xx=np.mgrid[:img.shape[0],:img.shape[1]]; r=np.hypot(xx-x0,yy-y0)*PIX
    edges=np.arange(0,202,4); prof=np.array([np.nanmedian(img[(r>=a)&(r<a+4)]) for a in edges]); rc=edges+2
    def ap(bg): return np.nansum(img[r<RAP]-bg)
    sky_far=np.nanmedian(img[(r>=150)&(r<200)]); sky_near=np.nanmedian(img[(r>=50)&(r<70)]); sky_in=np.nanmedian(img[(r>=25)&(r<40)])
    m=(rc>=25)&(rc<=70); p=np.polyfit(rc[m],prof[m],1); sky_lin=np.polyval(p,RAP/2)  # halo extrapolated in, evaluated at mean aperture radius ~8"
    Fnear=ap(sky_near); e=max(1.0857/CAT[b]["snr"],0.03)  # pipeline MIN_ERR floor
    print(f"\n== {b}: sky far {sky_far:.3f}, near(50-70\") {sky_near:.3f}, inner(25-40\") {sky_in:.3f}, halo-extrapolated(r~8\") {sky_lin:.3f}  [DN]")
    rows=[("far sky 150-200\" (no halo removal)",sky_far),("AllWISE near sky 50-70\" (catalogue)",sky_near),("inner sky 25-40\"",sky_in),("halo extrapolated to star (25-70\" linear fit)",sky_lin)]
    for lab,bg in rows:
        F=ap(bg); s=F/Fnear; mag=CAT[b]["mag"]-2.5*np.log10(s)
        out=f"  {lab:48s} flux/AllWISE {s:5.3f}  {b}={mag:6.3f}"
        for k in ("ratio_free","ratio_fixed"):
            R0=CAT[b][k]
            if R0: R=R0*s; sig=(1-1/R)/(0.4*np.log(10)*e); out+=f" | {k[6:]}: F/F* {R:5.2f}, {sig:5.1f} sigma"
        print(out)
    ax.plot(rc,prof-sky_far,"k.-",ms=3,label="median radial profile − far sky"); ax.axhline(sky_near-sky_far,c="C0",ls="--",label="50-70\" sky")
    ax.axhline(sky_in-sky_far,c="C1",ls="--",label="25-40\" sky"); rr=np.linspace(0,70,50); ax.plot(rr,np.polyval(p,rr)-sky_far,"C3",label="halo linear fit 25-70\"")
    ax.axvline(RAP,c="0.5",lw=.6); ax.set_yscale("symlog",linthresh=(sky_near-sky_far)/3); ax.set_xlim(0,200); ax.set_xlabel("radius (arcsec)"); ax.set_title(f"{b}",fontsize=9); ax.legend(fontsize=7)
fig.suptitle(f"{N}: background levels under the star",fontsize=9); fig.tight_layout(); fig.savefig(f"{D}/{N}_w4_background_subtraction.png",dpi=110)
