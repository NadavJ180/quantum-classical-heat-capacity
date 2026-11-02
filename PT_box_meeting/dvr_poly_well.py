"""
dvr_poly_well.py
The Poschl-Teller well (lambda = 2) inside a high-order polynomial trap instead of
hard walls, computed with the repo's own DVR pipeline (src/Quick_Scan.py):

    V(x) = (hbar^2 alpha^2 / 2m) (2x / L)^POWER  -  (hbar^2 alpha^2 / m) / cosh^2(alpha x)

The polynomial is ~0 for |x| < L/2 and reaches the binding energy hbar^2 alpha^2 / 2m
at x = +-L/2, so L plays the role of the box length (alpha L = trap length in well widths).

Nothing in src/ is changed: this script puts the potential into the config module in
memory (the same names a potential block of config.py sets) and calls
Quick_Scan.run_quick_scan. The quantum Cv comes from the DVR spectrum of V; the classical
Cv is the pipeline's xi-scan (hbar -> hbar/xi, xi^2 V re-solved by the DVR at every xi).

Two figures: the quick-scan window (fig_dvr_poly_well.png) and, from the same saved results,
every temperature the saved DVR levels cover (fig_dvr_poly_well_wide.png).

Run from this folder (about 70 minutes and ~5 GB of memory for the settings below):
    python dvr_poly_well.py               # quick scan, saves results/ and the figures
    python dvr_poly_well.py --plot-only   # redraw the figures from the saved results
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import config
import Quick_Scan

# ------------------------------------------------------------------
# CONSTANTS (units: hbar = m = k_B = 1, as in config.py)
# ------------------------------------------------------------------
ALPHA = 1.0             # the well's width is 1/alpha
L = 100.0               # where the polynomial reaches hbar^2 alpha^2 / 2m (alpha L = 100)
POWER = 10              # order of the polynomial (even)
RESOLUTION = 1          # key of config.QUICK_SCAN_PRESETS (1 = the quick scan)
T_MIN, T_MAX = 0.08, 0.13   # k_B T window around the anomaly
SHOW_EXACT_CLASSICAL = True  # also draw the classical Cv from the exact 1-D integral (a check, not the pipeline)
WIDE_T_MIN = 0.01       # cold end of the wide figure (its hot end is where the saved levels run out)

PARAMS = {"L": L, "alpha": ALPHA, "power": POWER}
RESULTS_FILE = os.path.join(HERE, "results", "dvr_poly_well_quick.npz")
FIGURE_FILE = os.path.join(HERE, "figures", "fig_dvr_poly_well.png")
FIGURE_WIDE_FILE = os.path.join(HERE, "figures", "fig_dvr_poly_well_wide.png")


def poly_pt_potential(x, p=PARAMS):
    """(hbar^2 alpha^2 / 2m) (2x/L)^power - (hbar^2 alpha^2 / m) / cosh^2(alpha x)."""
    depth = config.HBAR**2 * p["alpha"]**2 / config.MASS
    u = np.minimum(np.abs(p["alpha"] * x), 300.0)          # cosh(300)^2 still fits in a float
    return 0.5 * depth * (2.0 * x / p["L"])**p["power"] - depth / np.cosh(u)**2


def use_this_potential():
    """Point config at this potential, exactly as swapping a potential block in config.py would."""
    config.POTENTIAL_PARAMS = PARAMS
    config.my_potential = poly_pt_potential
    config.SYSTEM_NAME = "1-D PT well in x%d trap" % POWER
    config.T_UNITS_LABEL = r"$k_B T$  ($\hbar = m = k_B = 1$)"
    config.POTENTIAL_FORMULA = (r"$V = \frac{1}{2}\alpha^2 (2x/L)^{%d} - \alpha^2/\cosh^2(\alpha x)$,  "
                                r"$L$ = %g, $\alpha$ = %g" % (POWER, L, ALPHA))
    # The DVR picks dx from the momentum at its energy ceiling; at xi = 1 that is too coarse
    # for the well itself (width 1/alpha), so cap it at a quarter width (config.py's DVR_DX_MAX).
    config.DVR_DX_MAX = 0.25 / ALPHA


def cv_classical_exact(T_arr, p=PARAMS):
    """Cv_cl = 1/2 + beta^2 Var(V) over exp(-beta V): the exact classical heat capacity
    (kinetic 1/2 from the Gaussian momentum integral, the rest from the x-integral)."""
    x = np.linspace(-0.8 * p["L"], 0.8 * p["L"], 400001)
    v = poly_pt_potential(x, p)
    cv = []
    for T in T_arr:
        w = np.exp(-(v - v.min()) / T)
        mean = np.sum(w * v) / np.sum(w)
        var = np.sum(w * (v - mean)**2) / np.sum(w)
        cv.append(0.5 + var / T**2)
    return np.array(cv)


def cv_from_levels(E, T_arr):
    """Cv_q = beta^2 Var(E) over the Boltzmann weights of the levels E (sorted, E[0] lowest)."""
    cv = []
    for T in T_arr:
        w = np.exp(-(E - E[0]) / T)
        mean = np.sum(w * E) / np.sum(w)
        cv.append(np.sum(w * (E - mean)**2) / np.sum(w) / T**2)
    return np.array(cv)


def run_scan():
    use_this_potential()
    out = Quick_Scan.run_quick_scan("single", RESOLUTION, (1.0 / T_MAX, 1.0 / T_MIN))
    if not out["results"]:
        raise SystemExit("The quick scan failed (see the console output above).")
    r = out["results"][0]
    os.makedirs(os.path.dirname(RESULTS_FILE), exist_ok=True)
    np.savez(RESULTS_FILE, T=out["T_arr"], cv_q=r["cv_quantum"], cv_cl=r["cv_classical"],
             eps=r["error_estimate"], energies=r["energies"], resolution=RESOLUTION,
             L=L, alpha=ALPHA, power=POWER, seconds=out["seconds"],
             pipeline_figure=str(out["figure_path"]))
    print("Saved", RESULTS_FILE)


def figure_with_potential(E0):
    """The layout of both figures: the trap and a zoom on the well on the left (with the DVR
    bound state at E0), the Cv panel and the Cv_q - Cv_cl panel on the right."""
    fig = plt.figure(figsize=(11, 6.2))
    ax_trap = fig.add_axes([0.06, 0.58, 0.30, 0.34])
    ax_well = fig.add_axes([0.06, 0.10, 0.30, 0.34])
    ax = fig.add_axes([0.45, 0.42, 0.52, 0.50])
    ax_d = fig.add_axes([0.45, 0.10, 0.52, 0.25], sharex=ax)

    depth = config.HBAR**2 * ALPHA**2 / config.MASS
    x = np.linspace(-0.7 * L, 0.7 * L, 4001)
    ax_trap.plot(x, poly_pt_potential(x), color="k", lw=1.2)
    ax_trap.set_ylim(-1.15 * depth, depth)
    ax_trap.set_title("the trap: polynomial (order %d) + PT well, $L$ = %g" % (POWER, L), fontsize=10)
    xw = np.linspace(-4 / ALPHA, 4 / ALPHA, 801)
    ax_well.plot(xw, poly_pt_potential(xw), color="k", lw=1.2)
    xb = np.arccosh(np.sqrt(depth / -E0)) / ALPHA          # where V = E_b
    ax_well.plot([-xb, xb], [E0, E0], color="#2a78d6", lw=3,
                 label="bound state\n(DVR): E = %.6f" % E0)
    ax_well.axhline(0, color="0.6", lw=0.8, ls=":")
    ax_well.set_ylim(-1.5 * depth, 0.1 * depth)                # room for the legend under the well
    ax_well.locator_params(axis="y", nbins=4)
    ax_well.set_title(r"zoom on the well: depth $\hbar^2\alpha^2/m$ = %g" % depth, fontsize=10)
    ax_well.legend(loc="lower left", fontsize=8, frameon=False)
    for a in (ax_trap, ax_well):
        a.set_xlabel("x")
        a.set_ylabel("V(x)")

    ax.set_ylabel(r"$C_v / k_B$")
    ax.grid(alpha=0.3)
    ax.tick_params(labelbottom=False)
    ax_d.axhline(0, color="k", lw=0.8)
    ax_d.set_xlabel(r"$k_B T$  ($\hbar = m = k_B = 1$; bound state at $-\hbar^2\alpha^2/2m$ = %g)" % (-0.5 * depth))
    ax_d.set_ylabel(r"$\Delta C_v / k_B$")
    ax_d.grid(alpha=0.3)
    return fig, ax, ax_d


def plot_results():
    data = np.load(RESULTS_FILE)
    order = np.argsort(data["T"])
    T, cq, ccl, eps = (data[k][order] for k in ("T", "cv_q", "cv_cl", "eps"))
    E = data["energies"]
    d = cq - ccl
    above = d > Quick_Scan.ERROR_MARGIN * eps

    fig, ax, ax_d = figure_with_potential(E[0])
    ax.plot(T, cq, color="#2a78d6", lw=2, label="quantum (DVR spectrum of V)")
    ax.plot(T, ccl, color="#eb6834", lw=2, ls="--", label=r"classical limit (DVR, $\hbar\to\hbar/\xi$)")
    ax.errorbar(T, ccl, yerr=eps, fmt="none", ecolor="#eb6834", capsize=0)   # separate, so the legend keeps this order
    if SHOW_EXACT_CLASSICAL:
        ax.plot(T, cv_classical_exact(T), color="0.3", lw=1, ls=":",
                label="check: exact classical integral")
    if above.any():
        ax.fill_between(T, 0, 1, where=above, color="#1baf7a", alpha=0.15,
                        transform=ax.get_xaxis_transform(), label=r"$C_v^q > C_v^{cl}$ (resolved)")
    ax.set_title(r"Quick scan (resolution %d): $\alpha L$ = %g, polynomial order %d"
                 % (int(data["resolution"]), ALPHA * L, POWER), fontsize=11)
    ax.legend(fontsize=9, frameon=False, loc="upper left")

    ax_d.plot(T, d, color="#1baf7a", lw=2, label=r"$C_v^q - C_v^{cl}$")
    ax_d.fill_between(T, d - 2 * eps, d + 2 * eps, color="#1baf7a", alpha=0.25, lw=0,
                      label=r"$\pm 2\times$ error estimate ($\leq$ %.3f)" % (2 * np.nanmax(eps)))
    ax_d.legend(fontsize=9, frameon=False, loc="lower center")

    os.makedirs(os.path.dirname(FIGURE_FILE), exist_ok=True)
    fig.savefig(FIGURE_FILE, dpi=160)
    print("Saved", FIGURE_FILE)

    i = int(np.nanargmax(d))
    print("max Cv_q - Cv_cl = %+.3f (error estimate %.1e) at k_B T = %.4f: Cv_q = %.3f, Cv_cl = %.3f"
          % (d[i], eps[i], T[i], cq[i], ccl[i]))
    if above.any():
        print("Cv_q > Cv_cl (resolved) for k_B T = %.4f .. %.4f" % (T[above].min(), T[above].max()))
    if SHOW_EXACT_CLASSICAL:
        diff = ccl - cv_classical_exact(T)
        print("check: xi-scan - exact classical = %+.1e .. %+.1e (xi-scan error estimate %.1e .. %.1e)"
              % (np.nanmin(diff), np.nanmax(diff), np.nanmin(eps), np.nanmax(eps)))
        low = np.abs(diff) > Quick_Scan.ERROR_MARGIN * eps
        if low.any():
            print("       |xi-scan - exact| > %g x error estimate at %d temperatures, k_B T = %.4f .. %.4f"
                  % (Quick_Scan.ERROR_MARGIN, low.sum(), T[low].min(), T[low].max()))


def plot_wide():
    """The same comparison without zooming in, from the saved results only (no new DVR run).
    Cv_q comes from the saved DVR levels; by the quick scan's own rule (levels up to
    E_0 + thermal_coverage k_B T) they cover k_B T up to (E_max - E_0) / thermal_coverage,
    which is the hot end. Cv_cl is the exact classical integral, because the xi-scan exists
    only in the quick-scan window; its points are drawn on top there."""
    data = np.load(RESULTS_FILE)
    E = data["energies"]
    coverage = config.QUICK_SCAN_PRESETS[int(data["resolution"])]["thermal_coverage"]
    T = np.linspace(WIDE_T_MIN, (E[-1] - E[0]) / coverage, 400)
    cq, ccl = cv_from_levels(E, T), cv_classical_exact(T)
    d = cq - ccl
    order = np.argsort(data["T"])[::3]                         # every third quick-scan point
    T_xi, cq_xi, ccl_xi = (data[k][order] for k in ("T", "cv_q", "cv_cl"))

    fig, ax, ax_d = figure_with_potential(E[0])
    ax.plot(T, cq, color="#2a78d6", lw=2, label="quantum (saved DVR levels)")
    ax.plot(T, ccl, color="#eb6834", lw=2, ls="--", label=r"classical (exact integral of $e^{-\beta V}$)")
    ax.plot(T_xi, ccl_xi, "o", color="#eb6834", ms=4,
            label=r"classical limit (DVR, $\hbar\to\hbar/\xi$), quick-scan window")
    ax.fill_between(T, 0, 1, where=d > 0, color="#1baf7a", alpha=0.15,
                    transform=ax.get_xaxis_transform(), label=r"$C_v^q > C_v^{cl}$")
    ax.set_title(r"$\alpha L$ = %g, polynomial order %d: all $k_B T$ the saved levels cover ($\leq$ %.3f)"
                 % (ALPHA * L, POWER, T[-1]), fontsize=10.5)
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    ax_d.plot(T, d, color="#1baf7a", lw=2, label=r"$C_v^q - C_v^{cl}$")
    ax_d.plot(T_xi, cq_xi - ccl_xi, "o", color="#1baf7a", ms=4, label=r"with the $\xi$-scan")
    ax_d.legend(fontsize=9, frameon=False, loc="lower left")

    os.makedirs(os.path.dirname(FIGURE_WIDE_FILE), exist_ok=True)
    fig.savefig(FIGURE_WIDE_FILE, dpi=160)
    print("Saved", FIGURE_WIDE_FILE)
    s = np.where(np.diff(np.sign(d)) != 0)[0]
    crossings = [T[j] - d[j] * (T[j + 1] - T[j]) / (d[j + 1] - d[j]) for j in s]   # linear interpolation
    print("wide figure: k_B T = %.3f .. %.3f; Cv_q - Cv_cl changes sign at k_B T = %s"
          % (T[0], T[-1], ", ".join("%.4f" % t for t in crossings)))


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        run_scan()
    plot_results()
    plot_wide()
    plt.show()
