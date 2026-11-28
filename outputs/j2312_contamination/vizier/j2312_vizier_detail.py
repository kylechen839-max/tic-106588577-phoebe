"""Print the IR / spectral-type / stellar-parameter rows for J2312 from VizieR (Harvard mirror)."""
import sys, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(__file__).parent
c = SkyCoord(348.0058494, 53.3413022, unit="deg")
cats = {"II/125/main": 60, "II/167/catalog": 60, "II/275/fsr": 60, "II/275/assoc": 60, "II/297/irc": 60, "II/298/fis": 90,
        "V/114/msx6_main": 60, "VIII/106": 60, "B/mk/mktypes": 30, "III/135A/catalog": 30, "III/231/catalog": 10,
        "I/355/paramp": 3, "I/355/gaiadr3": 3, "I/352/gedr3dis": 3, "I/360/goldoba": 3, "I/358/veb": 3,
        "J/A+A/704/A155/auxcode": 3, "J/ApJ/984/58/table1": 3, "II/361/mdfc-v10": 5, "V/136/tycall": 5}
out = []
for cat, r in cats.items():
    v = Vizier(columns=["**", "_r"], row_limit=20); v.VIZIER_SERVER = "vizier.cfa.harvard.edu"
    try:
        res = v.query_region(c, radius=r * u.arcsec, catalog=cat)
    except Exception as e:
        out.append(f"=== {cat}: ERROR {e}"); continue
    if len(res) == 0:
        out.append(f"=== {cat} (r<{r}\"): no rows"); continue
    for k in res.keys():
        t = res[k]; t.sort("_r")
        out.append(f"=== {k} (r<{r}\"): {len(t)} rows, {t.meta.get('description','')}")
        for row in t[:3]:
            out.append("  " + "; ".join(f"{n}={row[n]}" for n in t.colnames if str(row[n]) not in ("--", "", "nan")))
txt = "\n".join(out); (OUT / "j2312_vizier_detail.txt").write_text(txt, encoding="utf-8"); print(txt)
