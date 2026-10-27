"""
cv_plot.py -- quantum vs classical heat capacity of a Poschl-Teller well in a box.

    V(x) = -(hbar^2 alpha^2 / 2m) * 2 / cosh^2(alpha x)   for |x| < L/2
    V(x) = +infinity                                    for |x| >= L/2

Edit the CONSTANTS block and run:   python cv_plot.py

Every formula below is written out (with its derivation) in
PT_Box_Summary.pdf and PT_Box_Derivation.pdf in this folder.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf


# ============================ CONSTANTS (edit these) ============================
HBAR = 1.0
M = 1.0
KB = 1.0

# One pair of curves (quantum solid, classical dashed) per entry.
# alpha = 1 / (well width), L = box length. Only alpha*L and k_B T / (hbar^2 alpha^2 / 2m) matter.
CASES = [
    {"alpha": 1.0, "L": 30.0},
    {"alpha": 1.0, "L": 100.0},
    {"alpha": 1.0, "L": 400.0},
    {"alpha": 1.0, "L": 1000.0},
]

T_MIN = 0.01        # k_B T range of the plot (same units as hbar^2 alpha^2 / m)
T_MAX = 0.5
N_T = 1000
LOG_T = False       # True: logarithmic temperature axis

SHOW_BOX_ONLY = True    # also show the empty box (no well): its own small Cv_q > Cv_cl
# ===============================================================================

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
GL_NODES, GL_WEIGHTS = np.polynomial.legendre.leggauss(1000)   # quadrature points on [-1, 1]


def cv_quantum(T, alpha, L):
    """Cv_q from Z_q = L/lambda_th - 1/2 + exp(a/2) erf(sqrt(a/2)),  a = hbar^2 alpha^2 / (m k_B T).

    Cv = k_B a^2 d^2(ln Z)/da^2 = k_B a^2 [Z''/Z - (Z'/Z)^2]   (' = d/da).
    Z, Z', Z'' are all divided by exp(a/2) (the bound-state weight) so nothing
    overflows at low T; the ratios Z''/Z and Z'/Z do not change.
    """
    a = HBAR**2 * alpha**2 / (M * KB * T)
    s = np.exp(-a / 2)
    box = alpha * L / np.sqrt(2 * np.pi * a)          # L / lambda_th
    g = erf(np.sqrt(a / 2))                           # bound-state term exp(a/2) erf(...), divided by exp(a/2)
    dg = g / 2 + s / np.sqrt(2 * np.pi * a)
    d2g = dg / 2 - s / (2 * np.sqrt(2 * np.pi) * a**1.5)

    Z = (box - 0.5) * s + g
    dZ = -box / (2 * a) * s + dg
    d2Z = 3 * box / (4 * a**2) * s + d2g

    return KB * a**2 * (d2Z / Z - (dZ / Z)**2)


def cv_classical(T, alpha, L):
    """Cv_cl from Z_cl = (L + I) / lambda_th,  I = integral of (exp(-V/k_B T) - 1) dx.

    With t = tanh(alpha x):  alpha*I = integral_{-1}^{1} (exp(a(1-t^2)) - 1) / (1-t^2) dt
    (= sum_l (2a)^l / (l (2l-1)!!), the closed form in the notes),
    alpha*J1 = d(alpha*I)/da,  alpha*J2 = d(alpha*J1)/da,
    Cv_cl = k_B/2 + k_B a^2 [J2/(L+I) - (J1/(L+I))^2].
    Everything is divided by exp(a) (classical well weight) so nothing overflows.
    """
    a = HBAR**2 * alpha**2 / (M * KB * T)

    # alpha*I*exp(-a) by Gauss-Legendre quadrature in t (smooth integrand, no end points)
    t, w = GL_NODES, GL_WEIGHTS
    I = np.zeros(len(a))
    for i in range(len(a)):
        f = np.exp(-a[i] * t**2) * (-np.expm1(-a[i] * (1 - t**2))) / (1 - t**2)
        I[i] = np.sum(w * f)

    J1 = np.sqrt(np.pi) * erf(np.sqrt(a)) / np.sqrt(a)
    J2 = J1 * (1 - 1 / (2 * a)) + np.exp(-a) / a
    W = alpha * L * np.exp(-a) + I

    return KB * (0.5 + a**2 * (J2 / W - (J1 / W)**2))


# ---------------- BOX-ONLY (delete from here ... ) ----------------
def cv_quantum_box_only(T, L):
    """Empty box, exact levels E_n = (hbar pi n / L)^2 / 2m.  Cv = k_B beta^2 Var(E).
    (The classical empty box has Cv_cl = k_B/2 at every T.)"""
    E1 = (HBAR * np.pi / L)**2 / (2 * M)
    cv = np.zeros(len(T))
    for i in range(len(T)):
        beta = 1 / (KB * T[i])
        n_max = int(np.sqrt(60 / (beta * E1))) + 3          # exp(-60) is negligible
        E = E1 * np.arange(1, n_max + 1)**2
        p = np.exp(-beta * (E - E[0]))
        p = p / np.sum(p)
        cv[i] = KB * beta**2 * (np.sum(p * E**2) - np.sum(p * E)**2)
    return cv


def plot_box_only():
    """Empty box: Cv_q vs Cv_cl = k_B/2 against k_B T / E_1 (same curve for every L)."""
    L = 1.0
    E1 = (HBAR * np.pi / L)**2 / (2 * M)
    T = E1 * np.logspace(-1, 3, 800)
    cq = cv_quantum_box_only(T, L)

    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(T / E1, cq / KB, color=COLORS[0], label="quantum")
    ax.plot(T / E1, 0.5 * np.ones(len(T)), "--", color=COLORS[1], label="classical (= 1/2)")
    ax.fill_between(T / E1, 0.5, cq / KB, where=cq / KB > 0.5, color=COLORS[0], alpha=0.15)
    i = np.argmax(cq)
    ax.annotate(r"max %.4f at $k_BT$ = %.2f $E_1$" % (cq[i] / KB, T[i] / E1), (T[i] / E1, cq[i] / KB),
                xytext=(30, -60), textcoords="offset points", fontsize=9,
                arrowprops={"arrowstyle": "->", "color": "0.3"})
    ax.set_xscale("log")
    ax.set_xlabel(r"$k_BT\,/\,E_1$,   $E_1=\pi^2\hbar^2/2mL^2$")
    ax.set_ylabel(r"$C_v/k_B$")
    ax.set_title("Empty box only (no well): the wall effect, at most +1/16")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
# ---------------- ... to here to remove the box-only option ----------------


def main():
    if LOG_T:
        T = np.logspace(np.log10(T_MIN), np.log10(T_MAX), N_T)
    else:
        T = np.linspace(T_MIN, T_MAX, N_T)

    fig, (ax, ax_d) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})

    print("case          alpha*L   max(Cv_q - Cv_cl)/k_B   at k_B T   at k_B T/(hbar^2 alpha^2/2m)")
    for i, case in enumerate(CASES):
        alpha, L = case["alpha"], case["L"]
        c = COLORS[i % len(COLORS)]
        cq = cv_quantum(T, alpha, L)
        ccl = cv_classical(T, alpha, L)
        name = r"$\alpha L$ = %g" % (alpha * L)
        if len(set(cc["alpha"] for cc in CASES)) > 1:
            name = r"$\alpha$ = %g, L = %g" % (alpha, L)

        ax.plot(T, cq / KB, color=c, label=name)
        ax.plot(T, ccl / KB, "--", color=c)
        ax.fill_between(T, ccl / KB, cq / KB, where=cq > ccl, color=c, alpha=0.2)

        ax_d.plot(T, (cq - ccl) / KB, color=c)
        ax_d.fill_between(T, 0, (cq - ccl) / KB, where=cq > ccl, color=c, alpha=0.2)

        # ---- BOX-ONLY (delete these 2 lines to remove) ----
        if SHOW_BOX_ONLY:
            ax_d.plot(T, cv_quantum_box_only(T, L) / KB - 0.5, ":", color=c, lw=1.2)

        j = np.argmax(cq - ccl)
        eps = HBAR**2 * alpha**2 / (2 * M)
        print("%-12d  %7g   %+10.3f              %8.4f   %8.4f" % (i + 1, alpha * L, (cq[j] - ccl[j]) / KB, T[j], KB * T[j] / eps))

    # legend for the line styles (colour = case)
    ax.plot([], [], "-", color="0.3", label="quantum")
    ax.plot([], [], "--", color="0.3", label="classical")
    # ---- BOX-ONLY (delete these 3 lines to remove) ----
    if SHOW_BOX_ONLY:
        ax_d.plot([], [], ":", color="0.3", label="empty box only")
        ax_d.legend(loc="upper right")
    ax.legend(ncol=2)
    ax.set_ylabel(r"$C_v/k_B$")
    ax.set_title(r"Pöschl–Teller well ($\lambda=2$) in a box: shaded where $C_v^{q}>C_v^{cl}$")
    ax_d.axhline(0, color="0.3", lw=0.8)
    ax_d.set_ylabel(r"$(C_v^{q}-C_v^{cl})/k_B$")
    ax_d.set_xlabel(r"$k_BT$   ($\hbar$ = %g, m = %g, $k_B$ = %g)" % (HBAR, M, KB))
    if LOG_T:
        ax_d.set_xscale("log")
    for a_ in (ax, ax_d):
        a_.grid(alpha=0.3)
    fig.tight_layout()

    # ---- BOX-ONLY (delete these 2 lines to remove) ----
    if SHOW_BOX_ONLY:
        plot_box_only()

    plt.show()


if __name__ == "__main__":
    main()
