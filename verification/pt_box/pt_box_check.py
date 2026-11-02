"""
pt_box_check.py
=====================================================================
Every number and figure quoted in PT_Box_Cv.tex comes from this script.

System: a reflectionless Poschl-Teller well (lambda = 2) centred in an
infinite square well of length L,

    V(x) = -(hbar^2 alpha^2 / m) sech^2(alpha x)   for |x| < L/2,
    V(x) = +infinity                                 otherwise.

Units: hbar = m = k_B = alpha = 1, so the energy unit is
eps = hbar^2 alpha^2 / (2m) = 1/2. Everything reported is dimensionless:
theta = k_B T / eps, Lambda = alpha L, E / eps, and Cv / k_B.

Sections (each prints its results and writes them to RESULTS.md):
  1. Exact levels from the closed-form quantization conditions.
  2. Levels against two DVRs: the repo's sinc-DVR (src/DVR), whose grid
     edges act as approximate walls, and a box sine-DVR with exact walls.
  3. General integer lambda, needed for the xi ladder: asymptotic
     quantization against the sine-DVR.
  4. Quantum Cv: Var(E) over the exact levels, against the closed-form Z_q.
  5. Classical Cv: closed form, cross-checked against quadrature and
     finite differences.
  6. The xi formula along the ladder xi_n = sqrt(n(n-1)/2), converging
     to the classical Cv.
  7. Cv_q against Cv_cl for several Lambda: the anomaly windows.
  8. The empty box: its own Cv_q > 1/2 overshoot and the
     lambda_th/(8L) high-T asymptote.

Runs on both of the repo's environments (Python 3.7 / numpy 1.18 and
Python 3.12 / numpy 2.x): no walrus operator, no np.trapezoid.
=====================================================================
"""

import os
import sys
import warnings

import numpy as np
import scipy.linalg as la
from scipy.integrate import quad
from scipy.special import erf

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass
warnings.filterwarnings("ignore", category=RuntimeWarning)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src", "DVR"))
from DVR_Algorithm import colbert_miller_dvr_1d  # noqa: E402  (repo sinc-DVR)

FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

EPS = 0.5          # energy unit eps = hbar^2 alpha^2 / 2m with hbar = m = alpha = 1
U0 = 2.0 * EPS     # classical depth of the lambda = 2 well
LINES = []         # RESULTS.md body


def out(s=""):
    print(s)
    LINES.append(s)


# ---------------------------------------------------------------------
# Plot style (static print figures): thin lines, recessive grid,
# palette slots 1-3 of the validated categorical palette.
# ---------------------------------------------------------------------
C_Q, C_CL, C_3 = "#2a78d6", "#eb6834", "#1baf7a"
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


# =====================================================================
# 1. Exact levels (lambda = 2, finite L)
# =====================================================================
def tau_of(L):
    """tau = alpha tanh(alpha L / 2)."""
    return np.tanh(L / 2.0)


def bound_kappa(L, iters=200):
    """Even bound state: kappa = tau tanh(kappa L / 2) (fixed point)."""
    tau, kap = tau_of(L), 1.0
    for _ in range(iters):
        kap = tau * np.tanh(kap * L / 2.0)
    return kap


def box_k(L, jmax, iters=200):
    """Continuum: k_j L = j pi + 2 arctan(k_j / tau), j = 1..jmax (vectorized fixed point)."""
    tau = tau_of(L)
    j = np.arange(1, jmax + 1, dtype=float)
    k = j * np.pi / L
    for _ in range(iters):
        k = (j * np.pi + 2.0 * np.arctan(k / tau)) / L
    resid = np.max(np.abs(k * L - j * np.pi - 2.0 * np.arctan(k / tau)))
    return k, resid


def exact_levels(L, jmax):
    kap = bound_kappa(L)
    k, resid = box_k(L, jmax)
    return np.concatenate([[-kap ** 2 / 2.0], k ** 2 / 2.0]), resid


def V_pt(x, lam=2):
    return -EPS * lam * (lam - 1) / np.cosh(x) ** 2


out("# PT (lambda = 2) in a box -- verification results")
out()
out("Written by `pt_box_check.py`. Units: hbar = m = k_B = alpha = 1, energy unit eps = 1/2;")
out("theta = k_B T / eps, Lambda = alpha L. All quoted in PT_Box_Cv.tex.")
out()
out("## 1. Exact levels (Lambda = 60)")
L60 = 60.0
E60, res60 = exact_levels(L60, 3000)
Efree = 0.5 * (np.arange(1, 6) * np.pi / L60) ** 2
out(f"- quantization residual max |k L - j pi - 2 arctan(k/tau)| = {res60:.1e}")
out(f"- bound state E_b/eps = {E60[0] / EPS:.15f}   (L -> infinity: -1)")
out("- first box levels E_j/eps vs the empty box (j pi / L)^2 / 2 / eps:")
out()
out("| j | PT + box | empty box | shift |")
out("|---|---|---|---|")
for j in range(1, 6):
    out(f"| {j} | {E60[j] / EPS:.8e} | {Efree[j - 1] / EPS:.8e} | {(E60[j] - Efree[j - 1]) / EPS:+.3e} |")
out()

# =====================================================================
# 2. DVR comparison (lambda = 2, L = 60, 120 levels)
# =====================================================================
def sine_dvr(Vf, L, N, nlev):
    """Colbert-Miller DVR on the interval (-L/2, L/2): exact Dirichlet walls, N-1 interior points."""
    i = np.arange(1, N)
    x = -L / 2.0 + L * i / N
    pref = 0.5 * 0.5 * (np.pi / L) ** 2          # (hbar^2 / 2m) * pi^2 / (2 L^2)
    I, J = np.meshgrid(i, i, indexing="ij")
    with np.errstate(divide="ignore"):
        T = pref * (-1.0) ** (I - J) * (1.0 / np.sin(np.pi * (I - J) / (2 * N)) ** 2
                                        - 1.0 / np.sin(np.pi * (I + J) / (2 * N)) ** 2)
    T[np.diag_indices(N - 1)] = pref * ((2.0 * N * N + 1.0) / 3.0 - 1.0 / np.sin(np.pi * i / N) ** 2)
    return la.eigvalsh(T + np.diag(Vf(x)), subset_by_index=[0, nlev - 1])


out("## 2. Exact levels vs DVR (Lambda = 60, 120 levels)")
out()
out("Absolute errors |E_DVR - E_exact| / eps.")
out()
out("| alpha dx | repo sinc-DVR: bound | repo sinc-DVR: levels 1-20 | repo sinc-DVR: levels 1-119 | sine-DVR: all 120 |")
out("|---|---|---|---|---|")
NL = 120
dvr_err = {}
for dx in (0.2, 0.1, 0.05, 0.025):
    n = int(round(L60 / dx)) - 1
    xg = (np.arange(n) - (n - 1) / 2.0) * dx      # grid edges one dx inside the walls
    Es = colbert_miller_dvr_1d(V_pt, NL, xg[0], xg[-1], n)
    Eb = sine_dvr(V_pt, L60, n + 1, NL)
    es, eb = np.abs(Es - E60[:NL]), np.abs(Eb - E60[:NL])
    dvr_err[dx] = (es, eb)
    out(f"| {dx} | {es[0] / EPS:.1e} | {es[1:21].max() / EPS:.1e} | {es[1:].max() / EPS:.1e} | {eb.max() / EPS:.1e} |")
out()

fig, ax = plt.subplots(figsize=(5.6, 3.0))
idx = np.arange(NL)
shades = {0.1: "#86b6ef", 0.05: "#3987e5", 0.025: "#1c5cab"}
for dx in (0.1, 0.05, 0.025):
    ax.semilogy(idx[1:], dvr_err[dx][0][1:] / EPS, color=shades[dx], lw=1.3,
                label=rf"repo sinc-DVR, $\alpha\Delta x={dx}$")
ax.semilogy(idx[1:], np.maximum(dvr_err[0.05][1][1:], 1e-16) / EPS, color=C_CL, lw=1.3,
            label=r"box sine-DVR, $\alpha\Delta x=0.05$")
ax.set_xlabel(r"level index $j$")
ax.set_ylabel(r"$|E_j^{\rm DVR}-E_j^{\rm exact}|\,/\,(\hbar^2\alpha^2/2m)$")
ax.set_title(r"Box levels of PT ($\lambda=2$) in a box, $\alpha L=60$")
ax.legend(loc="center right", frameon=False)
save(fig, "fig_dvr_levels")

# =====================================================================
# 3. General integer lambda: asymptotic quantization vs sine-DVR
# =====================================================================
def levels_int(lam, L, nlev):
    """Integer lambda, large alpha L: bound -eps p^2 (p = 1..lam-1); box k L = j pi + 2 sum_p arctan(k/p)."""
    Eb = sorted([-EPS * p * p for p in range(1, lam)])
    nbox = nlev - (lam - 1)
    j = np.arange(1, nbox + 1, dtype=float)
    k = j * np.pi / L
    for _ in range(400):
        k = (j * np.pi + 2.0 * sum(np.arctan(k / p) for p in range(1, lam))) / L
    return np.concatenate([Eb, k ** 2 / 2.0])


out("## 3. Integer lambda > 2 (asymptotic quantization) vs sine-DVR (Lambda = 40, 80 levels)")
out()
out("| lambda | bound states E/eps | max |E_DVR - E_asym| / eps, 80 levels |")
out("|---|---|---|")
for lam in (2, 3, 4, 6):
    Ed = sine_dvr(lambda x, lam=lam: V_pt(x, lam), 40.0, 1600, 80)
    Ea = levels_int(lam, 40.0, 80)
    out(f"| {lam} | {', '.join(f'{e / EPS:g}' for e in Ea[:lam - 1])} | {np.max(np.abs(Ed - Ea)) / EPS:.1e} |")
out()

# =====================================================================
# 4. Quantum Cv: Var(E) over exact levels vs closed-form Z_q
# =====================================================================
def cv_from_levels(E, beta):
    """Cv/k_B = beta^2 Var(E), with the exponent shift for stability."""
    w = np.exp(-beta * (E - E[0]))
    Z = w.sum()
    m1 = (E * w).sum() / Z
    m2 = (E * E * w).sum() / Z
    return beta * beta * (m2 - m1 * m1)


def lnZq_closed(beta, L, lam=2, eps=EPS, hbar=1.0):
    """Z_q = L/lambda_th - 1/2 + sum_p exp(beta p^2 eps) erf(p sqrt(beta eps))."""
    lth = np.sqrt(2.0 * np.pi * hbar ** 2 * beta)
    p = np.arange(1, lam)
    return np.log(L / lth - 0.5 + np.sum(np.exp(beta * p * p * eps) * erf(p * np.sqrt(beta * eps))))


def cvq_closed(beta, L, lam=2, eps=EPS, hbar=1.0):
    """Cv/k_B = beta^2 [Z''/Z - (Z'/Z)^2] from the closed-form Z_q and its analytic derivatives."""
    A = L / np.sqrt(2.0 * np.pi * hbar ** 2)
    c = eps * np.arange(1, lam) ** 2
    G = np.exp(beta * c) * erf(np.sqrt(beta * c))
    G1 = c * G + np.sqrt(c / (np.pi * beta))
    G2 = c * G1 - 0.5 * np.sqrt(c / np.pi) * beta ** -1.5
    Z = A * beta ** -0.5 - 0.5 + G.sum()
    Z1 = -0.5 * A * beta ** -1.5 + G1.sum()
    Z2 = 0.75 * A * beta ** -2.5 + G2.sum()
    return beta * beta * (Z2 / Z - (Z1 / Z) ** 2)


THETA = np.geomspace(0.02, 5.0, 500)


def levels_for(L, theta_max):
    """Enough exact levels that beta E_jmax >= 45 at the hottest theta."""
    jmax = int(L * np.sqrt(2.0 * 45.0 * theta_max * EPS) / np.pi) + 20
    E, _ = exact_levels(L, jmax)
    return E


out("## 4. Quantum Cv: Var(E) over the exact levels vs the closed-form Z_q")
out()
out("| Lambda | max rel. diff ln Z | max abs diff Cv_q | levels used |")
out("|---|---|---|---|")
for Lam in (60.0, 100.0, 1000.0):
    E = levels_for(Lam, THETA.max())
    dz, dc = 0.0, 0.0
    for th in THETA:
        b = 1.0 / (th * EPS)
        lz = np.log(np.sum(np.exp(-b * (E - E[0])))) - b * E[0]
        dz = max(dz, abs(lz / lnZq_closed(b, Lam) - 1.0))
        dc = max(dc, abs(cv_from_levels(E, b) - cvq_closed(b, Lam)))
    out(f"| {Lam:g} | {dz:.1e} | {dc:.1e} | {len(E)} |")
out()

# =====================================================================
# 5. Classical Cv: closed form
# =====================================================================
def I_of(a):
    """I(a) * alpha = int (e^{a sech^2 u} - 1) du = int_{-1}^{1} (e^{a(1-t^2)} - 1)/(1-t^2) dt."""
    return quad(lambda t: np.expm1(a * (1.0 - t * t)) / (1.0 - t * t), -1.0, 1.0, limit=200)[0]


def I_series(a, kmax=400):
    """(2a/alpha) 2F2(1,1;3/2,2;a) as its power series: sqrt(pi) sum_k a^k / (k Gamma(k+1/2))."""
    from scipy.special import gammaln
    k = np.arange(1, kmax + 1)
    return np.sqrt(np.pi) * np.sum(np.exp(k * np.log(a) - np.log(k) - gammaln(k + 0.5)))


def J1(a):
    """int sech^2 e^{a sech^2} dx = sqrt(pi) e^a erf(sqrt a) / sqrt a."""
    return np.sqrt(np.pi) * np.exp(a) * erf(np.sqrt(a)) / np.sqrt(a)


def J2(a):
    """int sech^4 e^{a sech^2} dx = J1 (1 - 1/(2a)) + 1/a."""
    return J1(a) * (1.0 - 0.5 / a) + 1.0 / a


def cvcl_closed(beta, L):
    """Cv_cl/k_B = 1/2 + beta^2 [I''/(L+I) - (I'/(L+I))^2], I' = U0 J1, I'' = U0^2 J2."""
    a = beta * U0
    D = L + I_of(a)
    return 0.5 + beta * beta * (U0 ** 2 * J2(a) / D - (U0 * J1(a) / D) ** 2)


def lnZcl(beta, L):
    return -0.5 * np.log(2.0 * np.pi * beta) + np.log(L + I_of(beta * U0))


out("## 5. Classical closed forms vs quadrature / finite differences")
out()
worst = {"I": 0.0, "J1": 0.0, "J2": 0.0}
for a in (0.05, 0.5, 2.0, 10.0, 50.0):
    worst["I"] = max(worst["I"], abs(I_series(a) / I_of(a) - 1.0))
    q1 = quad(lambda x: np.exp(a / np.cosh(x) ** 2) / np.cosh(x) ** 2, -60, 60, limit=400, points=[0])[0]
    q2 = quad(lambda x: np.exp(a / np.cosh(x) ** 2) / np.cosh(x) ** 4, -60, 60, limit=400, points=[0])[0]
    worst["J1"] = max(worst["J1"], abs(J1(a) / q1 - 1.0))
    worst["J2"] = max(worst["J2"], abs(J2(a) / q2 - 1.0))
dfd = 0.0
for th in (0.05, 0.2, 1.0, 4.0):
    b = 1.0 / (th * EPS)
    h = 1e-3 * b
    fd = b * b * (lnZcl(b + h, 100.0) - 2 * lnZcl(b, 100.0) + lnZcl(b - h, 100.0)) / h ** 2
    dfd = max(dfd, abs(fd - cvcl_closed(b, 100.0)))
out(f"- I(a): power series (2F2) vs quadrature, max rel. diff over a in [0.05, 50]: {worst['I']:.1e}")
out(f"- J1(a) erf form vs quadrature: {worst['J1']:.1e};  J2(a): {worst['J2']:.1e}")
out(f"- Cv_cl closed form vs finite-difference beta^2 d^2 lnZ_cl (Lambda = 100): max abs diff {dfd:.1e} (FD truncation)")
out()

# =====================================================================
# 6. xi ladder
# =====================================================================
def cv_ladder(theta, Lam, n):
    """Cv^(xi_n) for xi_n^2 = n(n-1)/2: lambda = n, hbar_eff = 1/xi_n, eps_n = eps / xi_n^2."""
    xi2 = n * (n - 1) / 2.0
    return cvq_closed(1.0 / (theta * EPS), Lam, lam=n, eps=EPS / xi2, hbar=1.0 / np.sqrt(xi2))


out("## 6. xi formula along the integer ladder xi_n = sqrt(n(n-1)/2)  (Lambda = 100)")
out()
NS = (2, 3, 5, 10, 20, 50, 100, 300, 1000, 3000)
out("| theta | Cv_cl | " + " | ".join(f"n={n} (xi={np.sqrt(n * (n - 1) / 2.0):.3g})" for n in NS) + " |")
out("|---|---|" + "---|" * len(NS))
ladder_err = {}
for th in (0.2, 0.4, 1.0):
    ref = cvcl_closed(1.0 / (th * EPS), 100.0)
    errs = [cv_ladder(th, 100.0, n) - ref for n in NS]
    ladder_err[th] = errs
    out(f"| {th} | {ref:.6f} | " + " | ".join(f"{e:+.1e}" for e in errs) + " |")
out()
out("(entries: Cv^(xi_n) - Cv_cl; n = 2 is the physical system itself, xi = 1)")
out()

fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8))
th_plot = np.geomspace(0.05, 5.0, 300)
cl = np.array([cvcl_closed(1.0 / (t * EPS), 100.0) for t in th_plot])
ramp = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#104281"]
for n, col in zip((2, 3, 5, 10, 40), ramp):
    xi = np.sqrt(n * (n - 1) / 2.0)
    axs[0].semilogx(th_plot, [cv_ladder(t, 100.0, n) for t in th_plot], color=col, lw=1.2,
                    label=rf"$\xi_{{{n}}}={xi:.3g}$")
axs[0].semilogx(th_plot, cl, color=C_CL, lw=1.6, ls="--", label=r"$C_v^{\rm cl}$ (closed form)")
axs[0].set_xlabel(r"$k_BT\,/\,(\hbar^2\alpha^2/2m)$")
axs[0].set_ylabel(r"$C_v/k_B$")
axs[0].set_title(r"$C_v^{(\xi_n)}$ along the ladder, $\alpha L=100$")
axs[0].legend(frameon=False, fontsize=7)
xis = np.array([np.sqrt(n * (n - 1) / 2.0) for n in NS])
for th, col in zip((0.2, 0.4, 1.0), (C_Q, C_3, INK2)):
    axs[1].loglog(xis, np.maximum(np.abs(ladder_err[th]), 1e-9), "o-", ms=3, color=col, lw=1.0,
                  label=rf"$k_BT={th}\,\hbar^2\alpha^2/2m$")
axs[1].loglog(xis, 0.05 / xis, color=INK2, lw=0.8, ls=":", label=r"$\propto 1/\xi$")
axs[1].set_xlabel(r"$\xi_n$")
axs[1].set_ylabel(r"$|C_v^{(\xi_n)}-C_v^{\rm cl}|\,/\,k_B$")
axs[1].set_title("convergence to the classical limit")
axs[1].legend(frameon=False, fontsize=7)
fig.tight_layout()
save(fig, "fig_xi_ladder")

# =====================================================================
# 7. Cv_q vs Cv_cl: anomaly windows
# =====================================================================
def windows(th, d, tol=1e-3):
    """Maximal theta-intervals where d > tol."""
    pos = d > tol
    ws, start = [], None
    for i, p in enumerate(pos):
        if p and start is None:
            start = i
        if (not p or i == len(pos) - 1) and start is not None:
            stop = i if p else i - 1
            ws.append((th[start], th[stop]))
            start = None
    return ws


out("## 7. Cv_q vs Cv_cl (exact levels for Cv_q, closed form for Cv_cl)")
out()
out("| Lambda | max (Cv_q - Cv_cl) | at theta | Cv_q there | Cv_cl there | windows Cv_q > Cv_cl (theta) | Cv_q peak (theta) | Cv_cl peak (theta) | min (Cv_q - Cv_cl) (theta) | wall term lth/8L at max |")
out("|---|---|---|---|---|---|---|---|---|---|")
curves = {}
for Lam in (30.0, 100.0, 400.0, 1000.0):
    E = levels_for(Lam, THETA.max())
    cq = np.array([cv_from_levels(E, 1.0 / (t * EPS)) for t in THETA])
    cc = np.array([cvcl_closed(1.0 / (t * EPS), Lam) for t in THETA])
    d = cq - cc
    i = int(np.argmax(d))
    wall = np.sqrt(4.0 * np.pi / THETA[i]) / (8.0 * Lam)
    ws = "; ".join(f"[{a:.3g}, {b:.3g}]" for a, b in windows(THETA, d)) or "none"
    out(f"| {Lam:g} | {d[i]:+.3f} | {THETA[i]:.3f} | {cq[i]:.3f} | {cc[i]:.3f} | {ws} | "
        f"{cq.max():.2f} ({THETA[np.argmax(cq)]:.3f}) | {cc.max():.2f} ({THETA[np.argmax(cc)]:.3f}) | "
        f"{d.min():+.2f} ({THETA[np.argmin(d)]:.3f}) | {wall:.1e} |")
    curves[Lam] = (cq, cc)
out()
out("(the upper window edge 5 is the end of the theta grid; above the dimple window the")
out(" remaining small excess is the wall term lambda_th/(8L), see section 8)")
out()

fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6), sharey=False)
for ax, Lam in zip(axs, (100.0, 400.0, 1000.0)):
    cq, cc = curves[Lam]
    ax.semilogx(THETA, cq, color=C_Q, label=r"$C_v^{\rm q}$")
    ax.semilogx(THETA, cc, color=C_CL, ls="--", label=r"$C_v^{\rm cl}$")
    ax.fill_between(THETA, cc, cq, where=cq > cc + 1e-3, color=C_Q, alpha=0.18, lw=0)
    ax.set_title(rf"$\alpha L={Lam:g}$")
    ax.set_xlabel(r"$k_BT\,/\,(\hbar^2\alpha^2/2m)$")
    ax.set_xlim(0.03, 5)
axs[0].set_ylabel(r"$C_v/k_B$")
axs[0].legend(frameon=False, loc="upper right")
fig.tight_layout()
save(fig, "fig_cv_compare")

# =====================================================================
# 8. Empty box
# =====================================================================
out("## 8. The empty box (no dimple)")
nn = np.arange(1, 6000)
tg = np.geomspace(0.05, 400.0, 4000)        # T in units of E_g = pi^2 hbar^2 / (2 m L^2)
cvb = np.array([cv_from_levels(nn ** 2.0, 1.0 / t) for t in tg])
ib = int(np.argmax(cvb))
first = tg[np.argmax(cvb > 0.5)]
t100 = tg[np.searchsorted(tg, 100.0)]
asym = np.sqrt(4.0 / (np.pi * t100)) / 8.0
out(f"- Cv_q max = {cvb[ib]:.4f} at T/E_g = {tg[ib]:.2f}; Cv_q > 1/2 for all T/E_g > {first:.2f}")
out(f"- at T/E_g = {t100:.0f}: Cv_q - 1/2 = {cvb[np.searchsorted(tg, 100.0)] - 0.5:.4f}; "
    f"wall term lambda_th/(8L) = {asym:.4f}")
out()

# Potential figure
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.6), gridspec_kw={"width_ratios": [1.4, 1]})
Lshow = 100.0
x = np.linspace(-Lshow / 2, Lshow / 2, 4001)
axs[0].plot(x, V_pt(x) / EPS, color=INK)
axs[0].vlines([-Lshow / 2, Lshow / 2], -2.2, 3.0, color=INK, lw=2.4)
axs[0].set_xlabel(r"$\alpha x$")
axs[0].set_ylabel(r"$V(x)\,/\,(\hbar^2\alpha^2/2m)$")
axs[0].set_title(r"full box, $\alpha L=100$")
axs[0].set_ylim(-2.3, 3.0)
xz = np.linspace(-4, 4, 801)
axs[1].plot(xz, V_pt(xz) / EPS, color=INK)
axs[1].axhline(-1.0, color=C_Q, lw=1.3)
axs[1].axhline(-2.0, color=C_CL, lw=1.3, ls="--")
axs[1].text(1.3, -0.93, r"bound state $E_b$", color=INK, fontsize=7.5, va="bottom")
axs[1].text(1.3, -1.93, r"classical floor $V_{\min}$", color=INK, fontsize=7.5, va="bottom")
axs[1].text(-3.9, 0.06, r"$V(x)$", color=INK, fontsize=7.5, va="bottom")
axs[1].set_xlabel(r"$\alpha x$")
axs[1].set_title("zoom on the dimple")
axs[1].set_ylim(-2.3, 0.5)
fig.tight_layout()
save(fig, "fig_potential")

with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(LINES) + "\n")
print("\nwrote RESULTS.md and figures/")
