import warnings; warnings.filterwarnings("ignore")
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord; from astropy.wcs import WCS; import astropy.units as u
from astroquery.skyview import SkyView; from astroquery.ipac.irsa import Irsa
N="J013155.91-004225.4"; c=SkyCoord(22.982958333,-0.707055556,unit="deg"); out=sys.argv[1]
a=Irsa.query_region(c,catalog="allwise_p3as_psd",spatial="Cone",radius=1*u.arcmin).to_pandas()
a["sep"]=c.separation(SkyCoord(a.ra,a.dec,unit="deg")).arcsec; print(a.sort_values("sep")[["designation","sep","w1mpro","w2mpro","w3mpro","w3snr","w4mpro","w4snr","ph_qual","cc_flags","ext_flg"]].head(8).to_string(index=False))
for W,px in ((2,120),(10,300)):
    im=SkyView.get_images(position=c,survey=["WISE 12","WISE 22"],pixels=px,width=W*u.arcmin,height=W*u.arcmin)
    fig,axs=plt.subplots(1,2,figsize=(9,4.7))
    for ax,b,h in zip(axs,["W3 (12 µm)","W4 (22 µm)"],im):
        h.writeto(f"{out}/{N}_{b[:2]}_{W}arcmin.fits",overwrite=True)
        d=h[0].data.astype(float); w=WCS(h[0].header); lo,hi=np.nanpercentile(d,[5,99.7]); ax.imshow(d,origin="lower",cmap="magma",vmin=lo,vmax=hi)
        x,y=w.world_to_pixel(c); s=W*60/px; ax.add_patch(plt.Circle((x,y),16.5/s if b[1]=="4" else 8.25/s,fill=False,color="cyan",lw=.8))
        ax.set_title(f"{b}, {W}'x{W}'",fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"{N} AllWISE atlas (N up, E left; cyan = target, 1x PSF FWHM radius)",fontsize=9); fig.tight_layout(); fig.savefig(f"{out}/{N}_wise_W3W4_{W}arcmin.png",dpi=120); plt.close(fig)
