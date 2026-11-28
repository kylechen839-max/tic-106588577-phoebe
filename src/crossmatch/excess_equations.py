"""Render the IR-excess equations (sanity_checks.fit_sed) as slide-ready PNGs, with J0719 W4 worked numbers.

Uses matplotlib mathtext (no LaTeX install needed). Writes outputs/crossmatch/excess_equations/:
one PNG per step, a combined sheet, and excess_equations.tex with the same maths.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "crossmatch" / "excess_equations"
OUT.mkdir(parents=True, exist_ok=True)
sys.stdout.reconfigure(encoding="utf-8")
INK, INK2, ACC, SURF = "#0b0b0b", "#52514e", "#2a78d6", "#fcfcfb"
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "Arial", "figure.facecolor": SURF, "savefig.facecolor": SURF})

STEPS = [
    ("1  Magnitude to flux density",
     [r"$F_\nu = F_0\,10^{-0.4\,m}\qquad \sigma_F = 0.4\ln 10\;F_\nu\,\sigma_m$"],
     r"J0719, W4:  $m = 6.744 \pm 0.076,\; F_0 = 8.363$ Jy  $\Rightarrow\; F_\nu = 16.78 \pm 1.17$ mJy   ($\sigma_m \geq 0.03$ mag)"),
    (r"2  Stellar photosphere (fitted to J, H, K$_s$, W1)",
     [r"$F_{\rm phot}(\lambda) = A\,B_\nu(\lambda, T)\qquad A = \dfrac{\sum_i w_i\,B_\nu(\lambda_i,T)\,F_i}{\sum_i w_i\,B_\nu(\lambda_i,T)^2},\quad w_i = \sigma_i^{-2}$",
      r"$T = \mathrm{arg\,min}_T\;\chi^2(T),\qquad T \in [2500,\,30000]$ K"],
     r"J0719:  $T = 16\,550$ K.  Diagnostic version: $F_{\rm phot} = A/\lambda^2$ fitted to W1, W2  $\Rightarrow F_{\rm phot}({\rm W4}) = 1.57$ mJy"),
    ("3  Excess ratio",
     [r"$R = \dfrac{F_{\rm obs}}{F_{\rm phot}}$"],
     r"J0719, W4:  $R = 16.78 / 1.57 = 10.7$     ($R = 1$: pure starlight)"),
    ("4  Excess significance",
     [r"$S = \dfrac{F_{\rm obs} - F_{\rm phot}}{\sigma_F}\qquad S \leq \dfrac{1}{0.4\ln 10 \times 0.03} \approx 36$"],
     r"J0719, W4:  $S = (16.78 - 1.57)/1.17 = 12.95\sigma$.   Excess counts if $S_{\rm W3}$ or $S_{\rm W4} \geq 5$ for both photospheres"),
    (r"5  Dust vs. free-free (fit to the residual $r_i = F_{\rm obs} - F_{\rm phot}$ in W2, W3, W4)",
     [r"dust: $r_i \approx a\,B_\nu(\lambda_i, T_d)$     free-free: $r_i \approx a\,\nu_i^{\,\alpha},\;\alpha \in [-0.1,\,0.6]$",
      r"$\chi^2 = \sum_i \left(\dfrac{r_i - {\rm model}_i}{\sigma_i}\right)^{2}\qquad {\rm dust\ wins\ if}\;\;\chi^2_{\rm dust} < \chi^2_{\rm ff}$"],
     r"J0719:  $T_d = 151$ K, $\chi^2_{\rm dust} = 0.12$;  $\chi^2_{\rm ff} = 120$   (diagnostic: 158 K, $\chi^2$ 0.02 vs 103)"),
    ("6  Fractional dust luminosity",
     [r"$\dfrac{L_{\rm IR}}{L_\ast} = \dfrac{a\,T_d^{\,4}}{A\,T^{4}}$"],
     r"J0719:  $L_{\rm IR}/L_\ast = 1.1 \times 10^{-4}$     (debris disks: $\leq 10^{-2}$)"),
]


def draw(ax, head, eqs, ex, y0=0.92):
    ax.axis("off")
    ax.text(0.0, y0, head, fontsize=15, color=INK, fontweight="bold", va="top", transform=ax.transAxes)
    y = y0 - 0.30
    for e in eqs:
        ax.text(0.04, y, e, fontsize=21, color=INK, va="center", transform=ax.transAxes)
        y -= 0.30
    ax.text(0.04, y + 0.06, ex, fontsize=12.5, color=ACC, va="center", transform=ax.transAxes)


for i, (head, eqs, ex) in enumerate(STEPS, 1):
    h = 1.5 + 0.75 * len(eqs)
    fig = plt.figure(figsize=(12, h))
    ax = fig.add_axes([0.03, 0.02, 0.95, 0.96])
    draw(ax, head, eqs, ex)
    fig.savefig(OUT / f"eq{i}.png", dpi=220); plt.close(fig)

fig, axs = plt.subplots(len(STEPS), 1, figsize=(12, 3.05 * len(STEPS)),
                        gridspec_kw=dict(height_ratios=[1.5 + 0.75 * len(e) for _, e, _ in STEPS]))
for ax, (head, eqs, ex) in zip(axs, STEPS):
    draw(ax, head, eqs, ex)
fig.suptitle("How the infrared excess is calculated (src/crossmatch/sanity_checks.py)", fontsize=17, color=INK,
             fontweight="bold", x=0.03, ha="left", y=0.995)
fig.subplots_adjust(left=0.03, right=0.98, top=0.965, bottom=0.01, hspace=0.05)
fig.savefig(OUT / "excess_equations_sheet.png", dpi=170); plt.close(fig)

TEX = r"""\documentclass[11pt]{article}
\usepackage{amsmath,amssymb}
\usepackage[margin=1in]{geometry}
\begin{document}
\section*{How the infrared excess is calculated}
Source: \texttt{src/crossmatch/sanity_checks.py} (\texttt{fit\_sed}). Worked numbers: J0719 (WISE J071951.40$-$240400.6), W4.

\paragraph{1. Magnitude to flux density}
\begin{equation}
F_\nu = F_0\,10^{-0.4\,m}, \qquad \sigma_F = 0.4\ln 10\;F_\nu\,\sigma_m \quad (\sigma_m \ge 0.03\ \text{mag})
\end{equation}
$F_0 = 1594,\,1024,\,666.7,\,309.5,\,171.8,\,31.67,\,8.363$ Jy for J, H, K$_s$, W1--W4.
J0719 W4: $m = 6.744 \pm 0.076 \Rightarrow F_\nu = 16.78 \pm 1.17$ mJy.

\paragraph{2. Stellar photosphere} (blackbody fitted to J, H, K$_s$, W1)
\begin{equation}
F_{\rm phot}(\lambda) = A\,B_\nu(\lambda,T), \qquad
A = \frac{\sum_i w_i\,B_\nu(\lambda_i,T)\,F_i}{\sum_i w_i\,B_\nu(\lambda_i,T)^2}, \qquad w_i = \sigma_i^{-2},
\end{equation}
with $T$ minimising $\chi^2$ on a 2500--30000 K grid (25 K steps), and again with $T$ fixed at the catalogue
$T_{\rm eff}$. J0719: $T = 16\,550$ K. The diagnostic version uses $F_{\rm phot} = A/\lambda^2$ fitted to W1 and W2,
giving $F_{\rm phot}({\rm W4}) = 1.57$ mJy.

\paragraph{3. Excess ratio}
\begin{equation}
R = \frac{F_{\rm obs}}{F_{\rm phot}} \qquad \text{J0719 W4: } R = 16.78/1.57 = 10.7
\end{equation}

\paragraph{4. Excess significance}
\begin{equation}
S = \frac{F_{\rm obs} - F_{\rm phot}}{\sigma_F} \qquad \text{J0719 W4: } S = \frac{16.78 - 1.57}{1.17} = 12.95
\end{equation}
Only the photometric error enters; $S \lesssim 1/(0.4\ln 10 \times 0.03) \approx 36$. An excess counts if $S_{\rm W3}$ or
$S_{\rm W4} \ge 5$ (AllWISE quality A/B) for both photospheres.

\paragraph{5. Dust versus free-free} Residuals $r_i = F_{\rm obs} - F_{\rm phot}$ in W2, W3, W4:
\begin{gather}
\text{dust: } r_i \approx a\,B_\nu(\lambda_i, T_d),\ T_d \in [40,1500]\ \text{K}; \qquad
\text{free-free: } r_i \approx a\,\nu_i^{\alpha},\ \alpha \in [-0.1, 0.6] \\
\chi^2 = \sum_i \left(\frac{r_i - \text{model}_i}{\sigma_i}\right)^2, \qquad \text{dust preferred if } \chi^2_{\rm dust} < \chi^2_{\rm ff}
\end{gather}
J0719: $T_d = 151$ K, $\chi^2_{\rm dust} = 0.12$, $\chi^2_{\rm ff} = 120$ (diagnostic: 158 K, 0.02 vs.\ 103).

\paragraph{6. Fractional dust luminosity}
\begin{equation}
\frac{L_{\rm IR}}{L_\ast} = \frac{a\,T_d^{4}}{A\,T^{4}} \qquad \text{J0719: } 1.1\times10^{-4}\ \ (\text{debris disks} \lesssim 10^{-2})
\end{equation}
\end{document}
"""
(OUT / "excess_equations.tex").write_text(TEX, encoding="utf-8")
print("wrote", sorted(p.name for p in OUT.iterdir()))
