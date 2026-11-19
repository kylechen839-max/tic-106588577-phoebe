"""J2218 per-star limb-darkening scan: chi2 of the per-sector fits vs (u1, u2), linear law.

Reads the _persector_persectorffi_ld1*_ld2*.json results (and the shared-law ldlin060/070 runs) and
writes timing-independent summary plot + json next to them.
"""
import json, re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "J221843.61+544715.0"
O = Path("outputs/candidates") / D
rows = []
for f in sorted(O.glob(f"{D}_persector_persectorffi_ld1*_ld2*.json")):
    if f.name.endswith("_partial.json"):
        continue
    r = json.loads(f.read_text())
    p = dict(zip(r["params"], r["v"]))
    rows.append({"tag": f.stem.split("persectorffi_")[1], "u1": float(r["ld1"]), "u2": float(r["ld2"]),
                 "chi2": r["chi2_per_sector"] if "chi2_per_sector" in r else r["chi2"], **{k: p[k] for k in ("incl", "k", "teffratio", "q")}})
for tag, u in (("ldlin060", 0.6), ("ldlin070", 0.7)):
    r = json.loads((O / f"{D}_persector_persectorffi_{tag}.json").read_text())
    p = dict(zip(r["params"], r["v"]))
    rows.append({"tag": tag, "u1": u, "u2": u, "chi2": r["chi2_per_sector"] if "chi2_per_sector" in r else r["chi2"],
                 **{k: p[k] for k in ("incl", "k", "teffratio", "q")}})
rows.sort(key=lambda x: x["chi2"])
(O / f"{D}_ld_grid.json").write_text(json.dumps(rows, indent=1))
for x in rows:
    print(f"{x['tag']:14s} u1 {x['u1']:.2f} u2 {x['u2']:.2f} chi2 {x['chi2']:8.1f} incl {x['incl']:.2f} k {x['k']:.3f} T2/T1 {x['teffratio']:.3f} q {x['q']:.3f}")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
sc = ax[0].scatter([x["u1"] for x in rows], [x["u2"] for x in rows], c=[x["chi2"] for x in rows], s=160, cmap="viridis_r")
for x in rows:
    ax[0].annotate(f"{x['chi2']:.0f}", (x["u1"], x["u2"]), textcoords="offset points", xytext=(8, 6), fontsize=8)
ax[0].plot([0.3, 0.4], [0.5, 0.6], alpha=0)  # keep limits
ax[0].axvspan(0.3, 0.4, color="0.85", zorder=0, label="expected primary (9500 K)")
ax[0].axhspan(0.5, 0.6, color="0.92", zorder=0, label="expected secondary (6700 K)")
ax[0].set_xlabel("primary linear LD u1"); ax[0].set_ylabel("secondary linear LD u2")
ax[0].set_title("per-sector fit chi2 (1514 bins)"); ax[0].legend(fontsize=7, loc="center right")
fig.colorbar(sc, ax=ax[0], label="chi2")
ch = [x["chi2"] for x in rows]
for key, lab in (("q", "q"), ("k", "k = R2/R1"), ("teffratio", "T2/T1")):
    ax[1].scatter(ch, [x[key] for x in rows], label=lab)
ax[1].set_xlabel("chi2"); ax[1].set_ylabel("value"); ax[1].legend(fontsize=8)
ax[1].set_title("geometry vs fit quality")
fig.suptitle(f"{D}: per-star limb darkening (2026-11-19)")
fig.tight_layout()
fig.savefig(O / f"{D}_ld_grid.png", dpi=130)
