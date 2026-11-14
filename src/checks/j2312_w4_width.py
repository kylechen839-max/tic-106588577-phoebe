"""J2312: Gaussian FWHM of the W4 source vs field point sources of similar brightness (native AllWISE tile; see j2312_psf_extraction.py)."""
import importlib.util,numpy as np
from scipy import ndimage
from scipy.optimize import least_squares
from pathlib import Path
spec=importlib.util.spec_from_file_location('m',str(Path(__file__).with_name('j2312_psf_extraction.py')));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
H=m.HALF;P=m.PIX
def gfwhm(st,x0,y0,rmax=30,quad=False):
    st=np.nan_to_num(st,nan=np.nanmedian(st));yy,xx=np.mgrid[:st.shape[0],:st.shape[1]];X=xx-x0;Y=yy-y0;msk=np.hypot(X,Y)*P<rmax
    def f(p):
        A,xc,yc,s,c0,cx,cy=p[:7];b=c0+cx*X+cy*Y
        if quad: b=b+p[7]*X*X+p[8]*Y*Y+p[9]*X*Y
        return (A*np.exp(-((xx-xc)**2+(yy-yc)**2)/(2*s**2))+b-st)[msk]
    p0=[st[int(round(y0)),int(round(x0))]-np.median(st),x0,y0,5,np.median(st),0,0]+([0,0,0] if quad else [])
    p=least_squares(f,p0).x;return 2.3548*abs(p[3])*P,p[0]
img,w=m.load('w4');bg=m.background(img);res=img-bg
noise=1.4826*np.nanmedian(np.abs(res-np.nanmedian(res)))
mx=ndimage.maximum_filter(res,size=31);pk=np.argwhere((res==mx)&(res>8*noise))
x0,y0=w.world_to_pixel(m.C);xi,yi=int(round(float(x0))),int(round(float(y0)))
tgt=img[yi-H:yi+H+1,xi-H:xi+H+1];tx,ty=H+float(x0)-xi,H+float(y0)-yi
for rm in (20,30,40):
  for q in (False,True):
    print('J2312 fit radius',rm,'quad' if q else 'plane',round(gfwhm(tgt,tx,ty,rm,q)[0],1))
fa=[]
for y,x in pk:
    if not (H+5<y<img.shape[0]-H-5 and H+5<x<img.shape[1]-H-5): continue
    if np.hypot(x-x0,y-y0)*P<300: continue
    fw,A=gfwhm(img[y-H:y+H+1,x-H:x+H+1],H,H);fa.append((A,fw))
fa=np.array(fa);At=gfwhm(tgt,tx,ty)[0]
_,At=gfwhm(tgt,tx,ty)
sel=fa[(fa[:,0]>At/2)&(fa[:,0]<At*2)][:,1]
print('field sources with amp within x2 of J2312:',len(sel),'FWHM percentiles 16/50/84/95:',np.round(np.percentile(sel,[16,50,84,95]),1),'fraction >= J2312 width:',np.mean(sel>=25.0).round(3))
