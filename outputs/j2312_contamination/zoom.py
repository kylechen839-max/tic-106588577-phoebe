import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits; from astropy.wcs import WCS; from astropy.coordinates import SkyCoord
c=SkyCoord(348.0058494,53.3413022,unit="deg"); fig,axs=plt.subplots(1,3,figsize=(13,4.3))
for ax,b in zip(axs[:2],["WISE12","WISE22"]):
    h=fits.open(f"out/J231201.40+532028.6_{b}_10arcmin.fits")[0]; img=h.data.astype(float); w=WCS(h.header); x0,y0=w.world_to_pixel(c)
    xi,yi=int(round(float(x0))),int(round(float(y0))); sub=img[yi-60:yi+61,xi-60:xi+61]; sky=np.nanpercentile(img,20)
    ax.imshow(sub,origin="lower",cmap="magma",vmin=sky,vmax=np.nanpercentile(sub,99.8)); ax.set_title(f"{b} 4'x4' linear",fontsize=8); ax.set_xticks([]);ax.set_yticks([])
    yy,xx=np.mgrid[:img.shape[0],:img.shape[1]]; r=np.hypot(xx-x0,yy-y0)*2; th=np.degrees(np.arctan2(yy-y0,-(xx-x0)))
    edges=np.arange(0,242,6); prof=[np.nanmedian(img[(r>=a)&(r<a+6)]) for a in edges]
    axs[2].plot(edges+3,(np.array(prof)-sky)/(prof[0]-sky),label=f"{b} median radial")
    if b=="WISE22":
        for lab,m in (("N",(th>45)&(th<135)),("S",(th<-45)&(th>-135)),("E",np.abs(th)<45),("W",np.abs(th)>135)):
            p=[np.nanmedian(img[(r>=a)&(r<a+6)&m]) for a in edges]; print(b,lab,"sky 50-70:",round(np.nanmedian(img[(r>=50)&(r<70)&m]),3),"150-200:",round(np.nanmedian(img[(r>=150)&(r<200)&m]),3))
axs[2].axhline(0,c="k",lw=.5); axs[2].set_yscale("symlog",linthresh=0.01); axs[2].set_xlabel("radius (arcsec)"); axs[2].legend(fontsize=7); axs[2].set_title("normalised radial profile (sky=20th pct of 10')",fontsize=8)
fig.tight_layout(); fig.savefig("out/J231201.40+532028.6_wise_zoom_profile.png",dpi=110)
