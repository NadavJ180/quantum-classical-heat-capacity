"""
make_figures.py -- the figures used in PT_Box_Summary.tex and PT_Box_Derivation.tex.

    python make_figures.py        (writes figures/*.pdf and figures/*.png)

Units: hbar = m = alpha = k_B = 1, so hbar^2 alpha^2 / 2m = 1/2.
Temperatures are plotted as k_B T / (hbar^2 alpha^2 / 2m), energies as E / (hbar^2 alpha^2 / 2m).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cv_plot import cv_quantum, cv_classical, cv_quantum_box_only, COLORS

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
if not os.path.isdir(FIG):
    os.makedirs(FIG)

EPS = 0.5      # hbar^2 alpha^2 / 2m with hbar = m = alpha = 1
C_Q, C_CL = COLORS[0], COLORS[1]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.5, "savefig.dpi": 200, "savefig.bbox": "tight",
})


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    fig.savefig(os.path.join(FIG, name + ".png"))
    plt.close(fig)


def exact_levels(L, n_box):
    """Bound state and first n_box box levels (alpha = 1), from
    kappa = tanh(L/2) tanh(kappa L/2)  and  k_j L = j pi + 2 arctan(k_j / tanh(L/2))."""
    c = np.tanh(L / 2)
    kappa = 1.0
    for _ in range(300):
        kappa = c * np.tanh(kappa * L / 2)
    j = np.arange(1, n_box + 1)
    k = j * np.pi / L
    for _ in range(300):
        k = (j * np.pi + 2 * np.arctan(k / c)) / L
    return -kappa**2 / 2, k**2 / 2


def plot_potential(ax, L=20.0):
    x = np.linspace(-L / 2, L / 2, 2000)
    V = -2 * EPS / np.cosh(x)**2
    Eb, Ebox = exact_levels(L, 12)
    Ebox = Ebox[Ebox < 1.6 * EPS]

    ax.plot(x, V / EPS, color=INK, lw=1.6)
    ax.plot([-L / 2, -L / 2], [-2.2, 1.75], color=INK, lw=3)
    ax.plot([L / 2, L / 2], [-2.2, 1.75], color=INK, lw=3)
    for E in Ebox:
        ax.plot([-L / 2, L / 2], [E / EPS, E / EPS], color=C_Q, lw=0.9)
    xb = np.arccosh(np.sqrt(2.0))                       # where V = E_b
    ax.plot([-xb, xb], [Eb / EPS, Eb / EPS], color=C_CL, lw=2.2)
    xs = x[np.abs(x) < 4]
    ax.plot(xs, Eb / EPS + 0.55 / np.cosh(xs), color=C_CL, lw=0.9, ls=":")  # bound state ~ 1/cosh(alpha x)

    ax.text(1.3, -1.0, r"bound state $E_b=-\hbar^2\alpha^2/2m$", fontsize=8, color=INK, va="center",
            bbox={"fc": "white", "ec": "none", "pad": 1})
    ax.text(0.5, -2.0, r"$\leftarrow$ classical floor $-\hbar^2\alpha^2/m$", fontsize=8, color=INK, va="center",
            bbox={"fc": "white", "ec": "none", "pad": 1})
    ax.text(-L / 2 + 0.4, 1.55, r"box levels $E_j>0$", fontsize=8, color=C_Q, va="center")
    ax.set_xlim(-L / 2 - 0.8, L / 2 + 0.8)
    ax.set_ylim(-2.3, 1.75)
    ax.set_xlabel(r"$\alpha x$")
    ax.set_ylabel(r"$V,\ E$   in units of $\hbar^2\alpha^2/2m$")
    ax.set_title(r"(a) potential and exact levels, $\alpha L=20$")


def plot_cv(ax, cases=(100.0, 1000.0)):
    th = np.linspace(0.02, 1.0, 1200)                   # k_B T / (hbar^2 alpha^2 / 2m)
    T = th * EPS
    for i, L in enumerate(cases):
        cq = cv_quantum(T, 1.0, L)
        ccl = cv_classical(T, 1.0, L)
        c = COLORS[i]
        ax.plot(th, cq, color=c, label=r"$C_v^{q}$, $\alpha L=%g$" % L)
        ax.plot(th, ccl, "--", color=c, label=r"$C_v^{cl}$, $\alpha L=%g$" % L)
        ax.fill_between(th, ccl, cq, where=cq > ccl, color=c, alpha=0.25, lw=0)
    ax.set_xlabel(r"$k_BT\ /\ (\hbar^2\alpha^2/2m)$")
    ax.set_ylabel(r"$C_v/k_B$")
    ax.set_title(r"(b) heat capacities; shaded: $C_v^{q}>C_v^{cl}$")
    ax.legend(loc="upper right")
    ax.set_xlim(0, 1.0)


# --- Figure 1: potential + levels | Cv  (summary and derivation) ---
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 2.6))
plot_potential(a1)
plot_cv(a2)
fig.tight_layout()
save(fig, "fig_summary")

# --- Figure 2: Cv for several alpha L with the difference (derivation) ---
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.0, 4.6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
th = np.linspace(0.02, 1.0, 1500)
T = th * EPS
for i, L in enumerate([30.0, 100.0, 400.0, 1000.0]):
    cq, ccl = cv_quantum(T, 1.0, L), cv_classical(T, 1.0, L)
    a1.plot(th, cq, color=COLORS[i], label=r"$\alpha L=%g$" % L)
    a1.plot(th, ccl, "--", color=COLORS[i])
    a2.plot(th, cq - ccl, color=COLORS[i])
    a2.fill_between(th, 0, cq - ccl, where=cq > ccl, color=COLORS[i], alpha=0.25, lw=0)
a1.plot([], [], color=INK2, label="quantum")
a1.plot([], [], "--", color=INK2, label="classical")
a1.legend(ncol=2)
a1.set_ylabel(r"$C_v/k_B$")
a2.axhline(0, color=INK2, lw=0.8)
a2.set_ylabel(r"$(C_v^{q}-C_v^{cl})/k_B$")
a2.set_xlabel(r"$k_BT\ /\ (\hbar^2\alpha^2/2m)$")
a2.set_xlim(0, 1.0)
fig.tight_layout()
save(fig, "fig_cv_cases")

# --- Figure 3: empty box only (derivation) ---
fig, ax = plt.subplots(figsize=(5.0, 2.8))
E1 = np.pi**2 / 2                                       # L = 1
T = E1 * np.logspace(-1, 3, 800)
cq = cv_quantum_box_only(T, 1.0)
r = np.sqrt(2 * np.pi / T)                              # lambda_th / L
ok = r < 1.2
ax.plot(T / E1, cq, color=C_Q, label="quantum, exact level sum")
ax.plot(T[ok] / E1, 0.5 + r[ok] * (1 - r[ok]) / (2 * (2 - r[ok])**2), ":", color=INK, lw=1.2,
        label=r"quantum, $Z=L/\lambda_{th}-1/2$")
ax.plot(T / E1, 0.5 + 0 * T, "--", color=C_CL, label="classical")
ax.fill_between(T / E1, 0.5, cq, where=cq > 0.5, color=C_Q, alpha=0.2, lw=0)
ax.set_xscale("log")
ax.set_xlabel(r"$k_BT\ /\ E_1$,   $E_1=\pi^2\hbar^2/2mL^2$")
ax.set_ylabel(r"$C_v/k_B$")
ax.legend(loc="lower right")
fig.tight_layout()
save(fig, "fig_box_only")

# --- Figure 4: scale-invariant stroke (derivation) ---
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 2.7))
T = np.linspace(0.005, 0.6, 1500)
for gamma, ls, lab in [(1.0, "-", r"$\gamma=1$"), (1.5, "-", r"$\gamma=1.5$")]:
    alpha, L = 1.0 / gamma, 400.0 * gamma
    cq, ccl = cv_quantum(T, alpha, L), cv_classical(T, alpha, L)
    lw = 1.5 if gamma == 1.0 else 2.6
    al = 1.0 if gamma == 1.0 else 0.45
    for ax, x in [(a1, T), (a2, gamma**2 * T)]:
        ax.plot(x, cq, color=C_Q, lw=lw, alpha=al, label=r"$C_v^{q}$, " + lab)
        ax.plot(x, ccl, "--", color=C_CL, lw=lw, alpha=al, label=r"$C_v^{cl}$, " + lab)
a1.set_xlabel(r"$k_BT$   ($\hbar=m=k_B=1$, before stroke $\alpha=1$)")
a1.set_title(r"(a) after $\alpha\to\alpha/\gamma$, $L\to\gamma L$ ($\alpha L=400$)")
a2.set_xlabel(r"$\gamma^2 k_BT$")
a2.set_title(r"(b) same curves against $\gamma^2T$: they coincide")
a1.set_ylabel(r"$C_v/k_B$")
a1.set_xlim(0, 0.3)
a2.set_xlim(0, 0.3)
a1.legend(fontsize=7)
fig.tight_layout()
save(fig, "fig_stroke")

print("figures written to", FIG)
