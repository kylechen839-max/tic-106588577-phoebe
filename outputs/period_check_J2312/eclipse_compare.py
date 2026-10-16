"""J2312: per-event primary vs secondary eclipse comparison at P = 6.849530274 d."""
import sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
raw, outdir = sys.argv[1], sys.argv[2]
P, T0 = 6.849530274, 2868.13335
W, OOT_IN, ECL, CORE = 0.08, 0.045, 0.037, 0.01
t, f, e, s = np.loadtxt(raw, unpack=True); s = s.astype(int)
out = []
def say(*a):
    line = " ".join(str(x) for x in a); print(line); out.append(line)

def quad_min(x, y, w):
    A = np.vander(x, 3)
    Aw = A * w[:, None]; yw = y * w
    coef, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    res = y - A @ coef; dof = max(len(x) - 3, 1)
    s2 = np.sum((res * w) ** 2) / dof            # scale errors by actual scatter
    cov = np.linalg.inv(Aw.T @ Aw) * s2
    a, b, c = coef
    if a > 0 and abs(-b / (2 * a)) < CORE:
        xv = -b / (2 * a); fmin = c - b * b / (4 * a)
        J = np.array([b * b / (4 * a * a), -b / (2 * a), 1.0])
    else:                                          # vertex outside window: use fitted minimum on the grid
        g = np.linspace(-CORE, CORE, 201); v = np.polyval(coef, g); i = np.argmin(v)
        xv, fmin = g[i], v[i]; J = np.array([xv * xv, xv, 1.0])
    return fmin, np.sqrt(J @ cov @ J), xv

def fwhm(x, y, depth, xmin):
    # bin to suppress noise, then interpolate half-depth crossings outward from the minimum
    bw = max(0.002, 1.5 * np.median(np.diff(np.sort(x))))
    edges = np.arange(-W, W + bw, bw); c = 0.5 * (edges[1:] + edges[:-1])
    idx = np.digitize(x, edges) - 1
    yb = np.array([np.median(y[idx == k]) if np.any(idx == k) else np.nan for k in range(len(c))])
    ok = np.isfinite(yb); c, yb = c[ok], yb[ok]
    lev = 1 - 0.5 * depth; i0 = np.argmin(np.abs(c - xmin))
    def cross(direction):
        i = i0
        while 0 <= i + direction < len(c):
            j = i + direction
            if yb[j] >= lev and yb[i] < lev:
                return c[i] + (lev - yb[i]) * (c[j] - c[i]) / (yb[j] - yb[i])
            i = j
        return np.nan
    return cross(+1) - cross(-1)

kmin = int(np.floor((t.min() - T0) / (P / 2))) - 1; kmax = int(np.ceil((t.max() - T0) / (P / 2))) + 1
events, skipped = [], []
for k in range(kmin, kmax + 1):
    tc = T0 + k * P / 2
    dph = (t - tc) / P; m = np.abs(dph) < W
    if m.sum() == 0: continue
    x, y, ey, sec = dph[m], f[m], e[m], np.bincount(s[m]).argmax()
    typ = "P" if k % 2 == 0 else "S"
    # coverage of the eclipse itself on a 0.001-phase grid
    grid = np.arange(-ECL, ECL, 0.001)
    cad = np.median(np.diff(np.sort(t[s == sec]))) / P
    cov = np.mean([np.any(np.abs(x - g - 0.0005) <= max(0.0005, cad / 2)) for g in grid])
    oot = np.abs(x) > OOT_IN
    left, right = np.sum(oot & (x < 0)), np.sum(oot & (x > 0))
    if cov < 0.6 or left < 3 or right < 3:
        skipped.append((k, typ, sec, round(cov, 2), left, right)); continue
    pl = np.polyfit(x[oot], y[oot], 1); base = np.polyval(pl, x)
    yn, en = y / base, ey / base
    core = np.abs(x) < CORE
    if core.sum() < 4:
        skipped.append((k, typ, sec, round(cov, 2), "core<4")); continue
    fmin, efmin, xv = quad_min(x[core], yn[core], 1 / en[core])
    depth = 1 - fmin
    events.append(dict(k=k, typ=typ, sec=sec, depth=depth, edepth=efmin, xv=xv,
                       fwhm=fwhm(x, yn, depth, xv), x=x, y=yn, e=en, cov=cov))

say(f"J231201.40+532028.6 eclipse comparison; P={P} d, T0=BTJD {T0}; raw file {raw}")
say(f"window +/-{W} phase; OOT {OOT_IN}-{W} linear normalisation; eclipse |dphase|<{ECL}; depth from quadratic in |dphase|<{CORE}")
say(f"events used: {len(events)}; skipped: {len(skipped)} -> " + "; ".join(str(x) for x in skipped))
say("")
say("Per-event: k type sector depth err  FWHM(phase)  vertex_offset(phase)  coverage")
for ev in events:
    say(f"  {ev['k']:4d} {ev['typ']} S{ev['sec']}  {ev['depth']:.4f} {ev['edepth']:.4f}  {ev['fwhm']:.4f}  {ev['xv']:+.4f}  {ev['cov']:.2f}")

def stats(v):
    v = np.asarray([x for x in v if np.isfinite(x)])
    if len(v) == 0: return np.nan, np.nan, 0
    return v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan), len(v)

say("")
say("Group | nP nS | depth_P +/- SE | depth_S +/- SE | dDepth sigma | FWHM_P +/- SE | FWHM_S +/- SE | dFWHM sigma | odd-even primary depth (sigma)")
summary = {}
for name, sel in [("S16", [16]), ("S17", [17]), ("S57", [57]), ("ALL", [16, 17, 57])]:
    E = [ev for ev in events if ev["sec"] in sel]
    Pd = stats([ev["depth"] for ev in E if ev["typ"] == "P"]); Sd = stats([ev["depth"] for ev in E if ev["typ"] == "S"])
    Pw = stats([ev["fwhm"] for ev in E if ev["typ"] == "P"]); Sw = stats([ev["fwhm"] for ev in E if ev["typ"] == "S"])
    dd = (Pd[0] - Sd[0]) / np.hypot(Pd[1], Sd[1]); dw = (Pw[0] - Sw[0]) / np.hypot(Pw[1], Sw[1])
    odd = stats([ev["depth"] for ev in E if ev["typ"] == "P" and (ev["k"] // 2) % 2])
    even = stats([ev["depth"] for ev in E if ev["typ"] == "P" and not (ev["k"] // 2) % 2])
    oe = (odd[0] - even[0]) / np.hypot(odd[1], even[1])
    if name == "ALL":   # sectors have different dilution (QLP vs SPOC): also do it with per-sector means removed
        pass
    say(f"{name} | {Pd[2]} {Sd[2]} | {Pd[0]:.4f} +/- {Pd[1]:.4f} | {Sd[0]:.4f} +/- {Sd[1]:.4f} | {dd:.1f} | "
        f"{Pw[0]:.4f} +/- {Pw[1]:.4f} | {Sw[0]:.4f} +/- {Sw[1]:.4f} | {dw:.1f} | "
        f"odd {odd[0]:.4f}+/-{odd[1]:.4f} (n={odd[2]}) even {even[0]:.4f}+/-{even[1]:.4f} (n={even[2]}) -> {oe:.1f}")
    summary[name] = (Pd, Sd, dd)

# per-sector differences combined (pooling sectors directly mixes QLP/SPOC dilution)
dif = [(summary[n][0][0] - summary[n][1][0], np.hypot(summary[n][0][1], summary[n][1][1])) for n in ("S16", "S17", "S57")]
w = np.array([1 / d[1] ** 2 for d in dif]); dm = np.sum(w * [d[0] for d in dif]) / w.sum(); de = 1 / np.sqrt(w.sum())
say(f"Inverse-variance combination of per-sector (P - S) depth differences: {dm:.4f} +/- {de:.4f} ({dm/de:.0f} sigma); "
    f"depth ratio S/P per sector: " + ", ".join(f"{n} {summary[n][1][0]/summary[n][0][0]:.3f}" for n in ("S16", "S17", "S57")))
say("Note: per-event depth errors are white-noise only; group SEs use the event-to-event scatter (S57 scatter ~0.01 is systematics).")
say("")
for n, sel in [("S16", 16), ("S17", 17), ("S57", 57)]:
    a = stats([ev["xv"] for ev in events if ev["typ"] == "P" and ev["sec"] == sel]); b = stats([ev["xv"] for ev in events if ev["typ"] == "S" and ev["sec"] == sel])
    say(f"{n}: secondary minus primary minimum time = {(b[0]-a[0])*P*1440:+.1f} +/- {np.hypot(a[1], b[1])*P*1440:.1f} min")
# timing: mean vertex offset per type (secondary displaced from 0.5 => eccentric => 2 distinct eclipses)
for typ in "PS":
    v = stats([ev["xv"] for ev in events if ev["typ"] == typ])
    say(f"ALL: mean quadratic-vertex offset {typ}: {v[0]:+.5f} +/- {v[1]:.5f} phase ({v[0]*P*24*60:+.1f} +/- {v[1]*P*24*60:.1f} min), n={v[2]}")

a = stats([ev["xv"] for ev in events if ev["typ"] == "P"]); b = stats([ev["xv"] for ev in events if ev["typ"] == "S"])
say(f"ALL: secondary at phase {0.5 + b[0] - a[0]:.5f} +/- {np.hypot(a[1], b[1]):.5f} -> offset {(b[0]-a[0])*P*1440:+.1f} +/- {np.hypot(a[1], b[1])*P*1440:.1f} min ({(b[0]-a[0])/np.hypot(a[1], b[1]):.0f} sigma) "
    "(a 3.42 d period with identical eclipses requires 0)")
# stacked profiles, per sector (dilution differs), chi2 summed
def binw(sec): return 0.002 if sec == 57 else 0.004
def stack(E, bw, align=False):
    edges = np.arange(-W, W + 1e-9, bw); c = 0.5 * (edges[1:] + edges[:-1])
    per = []
    for ev in E:
        x = ev["x"] - (ev["xv"] if align else 0); idx = np.digitize(x, edges) - 1
        per.append([ev["y"][idx == j].mean() if np.any(idx == j) else np.nan for j in range(len(c))])
    per = np.array(per); nev = np.sum(np.isfinite(per), 0)
    mu = np.nanmean(per, 0)
    se = np.where(nev >= 2, np.nanstd(per, 0, ddof=1) / np.sqrt(np.maximum(nev, 1)), np.nan) if len(E) > 1 else np.full(len(c), np.nan)
    return c, mu, se, nev
say("")
say("Stacked comparison: each event binned (0.002 phase S57, 0.004 QLP), bin value = mean over events, error = event-to-event std/sqrt(n_events) (needs >=2 events);")
say(f"chi2 of primary-minus-secondary over |dphase|<{ECL}; 'aligned' = each event shifted to its own fitted minimum (removes the timing offset, tests shape+depth only)")
stacks = {}
for align in (False, True):
    chi_tot = 0; dof_tot = 0
    for name, sel in [("S16", [16]), ("S17", [17]), ("S57", [57])]:
        Ep = [ev for ev in events if ev["sec"] in sel and ev["typ"] == "P"]
        Es = [ev for ev in events if ev["sec"] in sel and ev["typ"] == "S"]
        if len(Ep) < 2 or len(Es) < 2: say(f"{name}: <2 P or S events"); continue
        c, mp, sp, _ = stack(Ep, binw(sel[0]), align); _, ms, ss, _ = stack(Es, binw(sel[0]), align)
        d = mp - ms; ed = np.hypot(sp, ss); m = np.isfinite(d) & np.isfinite(ed) & (np.abs(c) < ECL)
        chi = np.sum((d[m] / ed[m]) ** 2); dof = m.sum()
        mo = np.isfinite(d) & np.isfinite(ed) & (np.abs(c) > OOT_IN)
        say(f"  {'aligned' if align else 'ephemeris'} {name} (nP={len(Ep)}, nS={len(Es)}): chi2 = {chi:.1f} / {dof} -> chi2_red = {chi/dof:.2f}"
            f"   [OOT control chi2_red = {np.sum((d[mo]/ed[mo])**2)/max(mo.sum(),1):.2f}, {mo.sum()} bins]")
        stacks[(name, align)] = (c, mp, sp, ms, ss, d, ed)
        chi_tot += chi; dof_tot += dof
    from scipy import stats as st
    p = st.chi2.sf(chi_tot, dof_tot)
    say(f"  {'aligned' if align else 'ephemeris'} ALL sectors: chi2 = {chi_tot:.1f} / {dof_tot} -> chi2_red = {chi_tot/dof_tot:.2f}, p(identical) = {p:.2e} (~{st.norm.isf(max(p/2,1e-300)):.1f} sigma)")

# plot: one column per sector; top = stacked P and S (ephemeris-centred), bottom = difference
secs = [n for n in ("S16", "S17", "S57") if (n, False) in stacks]
fig, ax = plt.subplots(2, len(secs), figsize=(4.2 * len(secs), 6.5), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]), squeeze=False)
for i, n in enumerate(secs):
    c, mp, sp, ms, ss, d, ed = stacks[(n, False)]
    ax[0, i].errorbar(c, mp, sp, fmt="o", ms=3, color="#1f4e9e", label="primary (phase 0)")
    ax[0, i].errorbar(c, ms, ss, fmt="s", ms=3, color="#d1495b", label="secondary (phase 0.5)")
    ax[0, i].set_title(f"{n}"); ax[1, i].errorbar(c, d, ed, fmt="o", ms=3, color="k"); ax[1, i].axhline(0, color="0.5")
    for a in ax[:, i]: a.axvspan(-ECL, ECL, color="0.92", zorder=0)
    ax[1, i].set_xlabel("phase from eclipse centre")
ax[0, 0].set_ylabel("normalised flux"); ax[1, 0].set_ylabel("primary - secondary"); ax[0, 0].legend(fontsize=8)
fig.suptitle(f"J2312 stacked eclipses at P = {P} d (errors: event-to-event scatter)")
fig.tight_layout(); fig.savefig(f"{outdir}/J2312_eclipse_comparison.png", dpi=150)
open(f"{outdir}/results.txt", "w").write("\n".join(out) + "\n")
