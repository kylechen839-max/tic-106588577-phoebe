"""How do the binned flux errors compare with the actual bin-to-bin noise?

For each fit given as designation:tag, loads the stored PHOEBE model ('final') and the
_binned.txt light curve, and reports:
  - chi2_red with the binned errors (as in the pipeline, up to the free scale factor);
  - the median error inside and outside the eclipses;
  - the white-noise level of the binned series (robust second-difference scatter out of eclipse);
  - chi2_red if every bin had that white-noise error, in and out of eclipse.
Offline (PHOEBE venv). Writes outputs/checks/flux_error_scale.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import phoebe

ROOT = Path(__file__).resolve().parents[2]


def rsig(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def check(des, tag):
    d = ROOT / "outputs" / "candidates" / des
    r = np.loadtxt(d / f"{des}_binned.txt")
    e = json.loads((d / f"{des}_ephemeris.json").read_text())["eclipses"]
    b = phoebe.load(str(d / f"{des}_phoebe_best{tag}.phoebe"))
    mod = np.asarray(b.get_value(qualifier="fluxes", context="model", model="final"), float)
    sig_used = np.asarray(b.get_value(qualifier="sigmas", context="dataset", dataset="lc01"), float)
    obs = r[:, 1]
    sig = r[:, 2]
    a = np.sum(obs * mod / sig_used**2) / np.sum(mod**2 / sig_used**2)
    res = obs - a * mod
    ph = r[:, 3] % 1
    ine = np.zeros(len(ph), bool)
    for p0, hw in ((e["phase1"], e["halfwidth1"]), (e["phase2"], e["halfwidth2"])):
        ine |= np.abs(((ph - p0 + 0.5) % 1) - 0.5) < hw
    o = np.argsort(ph)
    fo, oo = r[o, 1], ~ine[o]
    i = np.where(oo[1:-1] & oo[:-2] & oo[2:])[0] + 1
    white = rsig((fo[i] - 0.5 * (fo[i - 1] + fo[i + 1])) / np.sqrt(1.5))
    n = len(obs)
    out = {"designation": des, "tag": tag, "n_bins": n, "n_in_eclipse": int(ine.sum()),
           "chi2_red_used_sigmas": float(np.sum((res / sig_used) ** 2) / n),
           "median_err_in_eclipse_ppt": float(1e3 * np.median(sig_used[ine])),
           "median_err_out_ppt": float(1e3 * np.median(sig_used[~ine])),
           "white_noise_ppt": float(1e3 * white),
           "resid_rms_in_eclipse_ppt": float(1e3 * np.sqrt(np.mean(res[ine] ** 2))),
           "resid_rms_out_ppt": float(1e3 * np.sqrt(np.mean(res[~ine] ** 2))),
           "chi2_red_white_all": float(np.mean((res / white) ** 2)),
           "chi2_red_white_in_eclipse": float(np.mean((res[ine] / white) ** 2)),
           "chi2_red_white_out": float(np.mean((res[~ine] / white) ** 2))}
    return out


if __name__ == "__main__":
    rows = [check(*a.split(":")) for a in sys.argv[1:]]
    for o in rows:
        print(f"{o['designation']} {o['tag']}: chi2_red(used) {o['chi2_red_used_sigmas']:.2f}; err in/out "
              f"{o['median_err_in_eclipse_ppt']:.2f}/{o['median_err_out_ppt']:.2f} ppt; white {o['white_noise_ppt']:.2f} ppt; "
              f"resid rms in/out {o['resid_rms_in_eclipse_ppt']:.2f}/{o['resid_rms_out_ppt']:.2f} ppt; "
              f"chi2_red(white) all {o['chi2_red_white_all']:.1f}, in {o['chi2_red_white_in_eclipse']:.1f}, out {o['chi2_red_white_out']:.1f}")
    p = ROOT / "outputs" / "checks"; p.mkdir(parents=True, exist_ok=True)
    (p / "flux_error_scale.json").write_text(json.dumps(rows, indent=1))
