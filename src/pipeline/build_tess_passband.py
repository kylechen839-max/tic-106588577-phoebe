"""Build and install a blackbody-only TESS:T passband for PHOEBE 2.4 without network access.

PHOEBE normally downloads TESS:T from tables.phoebe-project.org. Cloud containers block that host,
so this rebuilds the passband from the TESS response curve (data/passbands/tess.ptf, copied from
github.com/phoebe-project/phoebe2-tables). Only blackbody intensities are computed, which is all
run_phoebe_candidate.py uses (atm=blackbody, manual linear limb darkening). Verified 2026-10-03:
the J141909.42-565518.1 best model gives chi2 = 210.41, identical to the result made with the
official table.

Usage: .venv-phoebe/bin/python src/pipeline/build_tess_passband.py
"""
import os
from pathlib import Path

import astropy.units as u
import phoebe
from phoebe.atmospheres import passbands

ROOT = Path(__file__).resolve().parents[2]

if "TESS:T" in phoebe.list_installed_passbands():
    print("TESS:T already installed")
else:
    pb = passbands.Passband(ptf=str(ROOT / "data" / "passbands" / "tess.ptf"), pbset="TESS", pbname="T",
                            effwl=7972.3, wlunits=u.AA, calibrated=True,
                            reference="TESS response from phoebe2-tables (blackbody only, built locally)",
                            version=1.0, comments="offline blackbody-only rebuild")
    pb.compute_blackbody_response()
    out = os.path.expanduser("~/tess_T_bb.fits")
    pb.save(out)
    phoebe.install_passband(out, local=True)
    print("installed:", phoebe.list_installed_passbands())
