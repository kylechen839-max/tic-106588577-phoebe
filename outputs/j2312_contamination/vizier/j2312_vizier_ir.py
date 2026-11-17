"""Search VizieR (Harvard mirror) for all catalogues with a source near J2312, and list the IR ones.
python outputs/j2312_contamination/vizier/j2312_vizier_ir.py"""
import sys, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(__file__).parent

c = SkyCoord(348.0058494, 53.3413022, unit="deg")
v = Vizier(columns=["**", "_r"], row_limit=50); v.VIZIER_SERVER = "vizier.cfa.harvard.edu"
res = v.query_region(c, radius=30 * u.arcsec)
lines = []
for name in res.keys():
    t = res[name]
    lines.append(f"{name}\t{len(t)}\tmin_r={min(t['_r']) if '_r' in t.colnames else 'NA'}\t{t.meta.get('description','')}")
(OUT / "j2312_vizier_all_30arcsec.txt").write_text("\n".join(lines), encoding="utf-8")
print(len(lines), "tables"); print("\n".join(lines))
