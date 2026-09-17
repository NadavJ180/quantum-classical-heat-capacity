"""
homogeneous_scaling_check.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Advisor-requested sanity check, prompted by docs/refs/qvsclas book
chapter (1).pdf (Levy & Gelbwaser-Klimovsky, "Quantum Features and
Signatures of Quantum Thermal Machines", Eqs. 4.15-4.17): the simple
formula for an Otto cycle's heat/work as an INTEGRAL OVER Cv,

    Q_h = \\int_{T_c/q}^{T_h} Cv dT,   Q_c = -q \\int Cv dT,   W = (q-1) \\int Cv dT,

is only valid if the potential deformation during the adiabatic
strokes HOMOGENEOUSLY scales the energy levels, i.e. E_n(2) = q *
E_n(1) with q the SAME for every n (chapter Sec. 4.2.1). Before this
project can use that shortcut for the quartic/cubic/quadratic double
well, we need to confirm the potential family actually admits such a
scaling.

THE CLAIM BEING TESTED
---------------------------------------------------------------------
For a potential of the form V(x) = A x^4 + B x^3 + C x^2 (mass and
hbar fixed at 1, as everywhere else in this project), substituting
x = g*y shows

    A/g^6 * (g y)^4 + B/g^5 * (g y)^3 + C/g^4 * (g y)^2
        = (1/g^2) * [A y^4 + B y^3 + C y^2],

i.e. rescaling the coefficients by (A/g^6, B/g^5, C/g^4) is exactly
equivalent to stretching x by g and shrinking the whole potential by
1/g^2. Plugging this into the Schroedinger equation and substituting
x = g*y shows the eigenvalues must satisfy

    E_n(g) = E_n(g=1) / g^2   for EVERY n,

a pure, n-independent overall rescaling -- exactly the homogeneous-
scaling condition the Otto-cycle shortcut needs. This is a dimensional-
analysis identity, not an approximation: with mass/hbar held fixed, it
should hold to numerical precision for ANY (A, B, C), not just special
values. This script verifies that our own DVR pipeline actually
reproduces it (no bug silently breaking the scaling), and reproduces
the specific check the advisor asked for: for 2-3 of the potentials
already used in Cv_Coefficient_Sweep.py's own `b`-sweep, sweep g,
rescale (A, B, C) -> (A/g^6, B/g^5, C/g^4), and confirm that Cv(T)
plotted against T normalized by the level-0/1 gap (E1-E0) collapses
onto ONE g-independent curve.

WHAT'S IN THIS FOLDER
---------------------------------------------------------------------
homogeneous_scaling_check.py   This script (run it directly).
figures/                       Output plots (created on first run).
RESULTS.md                     Written by this script: the printed
                                eigenvalue-ratio table and collapse
                                error report, saved for the record.
=====================================================================
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from DVR.DVR_Algorithm import auto_configure_dvr, get_fully_converged_energy_levels
from Quantum_Classical_Combined import compute_quantum_heat_capacity_curve

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Fixed a, c (this project's base double well) -- b is the same
# SCAN_PARAM already used by Cv_Coefficient_Sweep.py. Three of the
# curves already plotted in cv_coefficient_sweep.png: the two swept
# extremes and the base value.
A_BASE, C_BASE = 0.25, -0.5
B_LINES = [-0.9, -0.5, -0.1]

# Compression/expansion factors to sweep. g=1 recovers the unscaled
# potential; g<1 and g>1 both tested (not just one direction).
G_VALUES = [0.75, 1.0, 1.25, 1.5]

NUM_STATES = 350
N_EIGEN_CHECK = 6          # how many low-lying levels to report in the ratio table
X_GRID = np.geomspace(1e-2, 15.0, 600)   # shared T / (E1-E0) axis for every g


def make_potential(A, B, C):
    def V(x):
        return A * x**4 + B * x**3 + C * x**2
    return V


def solve(A, B, C, num_states=NUM_STATES):
    V = make_potential(A, B, C)
    x_min, x_max, n_grid = auto_configure_dvr(V, num_states, mass=1.0, hbar=1.0)
    E = get_fully_converged_energy_levels(V, num_states, x_min, x_max, n_grid, mass=1.0, hbar=1.0)
    return E


def run_one_line(b, report_lines):
    A, C = A_BASE, C_BASE
    B = b
    report_lines.append(f"\n{'='*72}\nLine: (A, B, C) = ({A}, {B}, {C})   [b={b} in the project's own sweep]\n{'='*72}")

    per_g = {}
    for g in G_VALUES:
        Ag, Bg, Cg = A / g**6, B / g**5, C / g**4
        E = solve(Ag, Bg, Cg)
        delta01 = float(E[1] - E[0])
        per_g[g] = {"E": E, "delta01": delta01}
        print(f"  [line b={b}] g={g:g}: solved ({len(E)} states), "
              f"E0={E[0]:.6g} E1={E[1]:.6g} delta01={delta01:.6g}")

    # ---- Check A: eigenvalue-level scaling, g^2 * E_n(g) vs E_n(g=1) ----
    E_ref = per_g[1.0]["E"]
    report_lines.append("\nEigenvalue check: g^2 * E_n(g) / E_n(g=1)  (expect 1.0000 for every n, g)")
    header = "  n  | " + " | ".join(f"g={g:<6g}" for g in G_VALUES)
    report_lines.append(header)
    for n in range(N_EIGEN_CHECK):
        row = [f"{n:>3d}  |"]
        for g in G_VALUES:
            ratio = (g**2 * per_g[g]["E"][n]) / E_ref[n]
            row.append(f" {ratio:8.5f} |")
        report_lines.append(" ".join(row))

    # ---- Check B: Cv(T/(E1-E0)) collapse across g ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(G_VALUES)))

    cv_by_g = {}
    for g, color in zip(G_VALUES, colors):
        delta01 = per_g[g]["delta01"]
        T_phys = X_GRID * delta01
        beta_arr = 1.0 / T_phys
        cv = compute_quantum_heat_capacity_curve(per_g[g]["E"], beta_arr, xi=1.0)
        cv_by_g[g] = cv
        axes[0].plot(T_phys, cv, color=color, label=f"g={g:g}")
        axes[1].plot(X_GRID, cv, color=color, label=f"g={g:g}  (\u03b401={delta01:.3g})")

    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"physical $k_B T$")
    axes[0].set_ylabel(r"$C_v/k_B$")
    axes[0].set_title("NOT normalized -- curves should differ")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3, linestyle="--")

    axes[1].set_xscale("log")
    axes[1].set_xlabel(r"$k_B T \,/\, (E_1-E_0)$")
    axes[1].set_ylabel(r"$C_v/k_B$")
    axes[1].set_title("Normalized by (E1-E0) -- expect ONE curve")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3, linestyle="--")

    fig.suptitle(f"Homogeneous-scaling check: V(x)={A}x^4{B:+g}x^3{C:+g}x^2, coefficients scaled "
                 r"(A/g$^6$, B/g$^5$, C/g$^4$)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    b_slug = f"{b:g}".replace("-", "m").replace(".", "p")
    out_path = os.path.join(FIG_DIR, f"collapse_b_{b_slug}.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  saved {out_path}")

    # Collapse-quality metric: max |Cv_g(x) - Cv_{g=1}(x)| over the shared x grid
    ref_curve = cv_by_g[1.0]
    report_lines.append("\nCollapse quality: max_x |Cv_g(x) - Cv_{g=1}(x)|, x = T/(E1-E0)")
    for g in G_VALUES:
        max_dev = float(np.nanmax(np.abs(cv_by_g[g] - ref_curve)))
        report_lines.append(f"  g={g:<6g}: max deviation = {max_dev:.2e}")

    return out_path


def main():
    report_lines = [
        "Homogeneous energy scaling check",
        "=================================",
        "Potential family: V(x) = A x^4 + B x^3 + C x^2  (mass=hbar=1)",
        f"Lines tested (A, C fixed at project base = {A_BASE}, {C_BASE}; b swept as in Cv_Coefficient_Sweep.py): {B_LINES}",
        f"g values: {G_VALUES}",
        f"NUM_STATES per solve: {NUM_STATES}",
    ]
    figs = []
    for b in B_LINES:
        figs.append(run_one_line(b, report_lines))

    report_text = "\n".join(report_lines) + "\n"
    print("\n" + report_text)
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("# Homogeneous energy scaling check -- results\n\n```\n" + report_text + "```\n")
    print("Wrote RESULTS.md")
    print("Figures:", figs)


if __name__ == "__main__":
    main()
