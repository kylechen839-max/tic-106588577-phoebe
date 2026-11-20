import warnings; warnings.filterwarnings("ignore")
from astroquery.vizier import Vizier; from astropy.coordinates import SkyCoord; import astropy.units as u
c=SkyCoord(348.0058494,53.3413022,unit="deg"); v=Vizier(columns=["**","+_r"],row_limit=5)
T={"Gaia DR3 vari classifier":"I/358/vclassre","Gaia DR3 eclipsing binaries":"I/358/veb","VSX":"B/vsx/vsx","ASAS-SN variables":"II/366/catv2021",
"TESS EBs (Prsa+2022)":"J/ApJS/258/16/tess-ebs","TESS OBA EBs (IJspeert+2021)":"J/A+A/652/A120/eb-cat","ZTF periodic (Chen+2020)":"J/ApJS/249/18/table2",
"ATLAS variables (Heinze+2018)":"J/AJ/156/241/table4","WISE variables (Chen+2018)":"J/ApJS/237/28/table2","CRTS (Drake+2014)":"J/ApJS/213/9/table3",
"Kepler EBs":"J/AJ/151/68/table1","ASAS-3":"II/264/asas3","NSVS (Hoffman+2009)":"J/AJ/138/466/variables","EB catalogue (Avvakumova+2013)":"J/AN/334/860/catalog",
"GCVS":"B/gcvs/gcvs_cat","SB9":"B/sb9/main","TESS EBs sectors 1-26 north (Justesen/Howard)":"J/ApJS/247/56/table1","Gaia DR3 vari summary":"I/358/varisum"}
for k,t in T.items():
    try: r=v.query_region(c,radius=5*u.arcsec,catalog=t)
    except Exception as e: print(f"{k} [{t}]: ERROR {e}"); continue
    if not len(r): print(f"{k} [{t}]: no match"); continue
    tb=r[0]; print(f"{k} [{t}]: MATCH"); print(tb.to_pandas().iloc[0].dropna().to_string()); print()
