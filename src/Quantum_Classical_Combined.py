"""
Quantum_Classical_Combined.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
General-purpose (system-agnostic) Cv pipeline: given an array of
energy levels for ANY system, this file computes the true quantum
Cv(T) curve directly from the spectrum, calls into
Classical_Limit_Numerical.py to find the numerical classical
limit Cv(T) across the same temperature range, and produces all the
diagnostic + summary plots (xi-convergence diagnostic, n-convergence
diagnostic, and the combined Cv(T) curve plot). Everything here is
driven by `run(energies, potential_func, ...)`: the quantum Cv comes
from the pre-computed spectrum (whatever produced it), while the
classical limit needs the potential itself, because the xi-scan
re-solves the scaled potential xi^2 V at every xi (see
Classical_Limit_Numerical.py).

This file does NOT contain any hard-coded physical systems (no Box,
no HO, no Double Well). Those live in their own driver files/sections
(Quantum_HO_Master.py, Cv_AutoTune.py, Cv_Coefficient_Sweep.py) so
this stays a reusable, system-agnostic pipeline -- callers may invoke
`run()` more than once per process (e.g. the auto-tune escalation
loop in Quantum_HO_Master.py re-runs it with a growing NUM_STATES/
XI_START each round), and each call's diagnostic plots simply
overwrite the previous ones on disk.
=====================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from Classical_Limit_Numerical import compute_cv, sweep_temperature_range
from figures.output_paths import save_figure


# =====================================================================
# True quantum Cv(T) curve, directly from a (finite) energy spectrum
# =====================================================================
def compute_quantum_heat_capacity_curve(energies, beta_arr, xi=1.0):
    """
    Evaluate the literal quantum heat capacity Cv(T) (no xi-rescaling
    by default) at every beta in `beta_arr`, using the given,
    truncated energy spectrum directly.

    Parameters
    ----------
    energies : array_like
        Energy eigenvalues for the system (e.g. from DVR).
    beta_arr : array_like
        Inverse temperatures to evaluate Cv at.
    xi : float, optional
        Spectrum scaling factor (default 1.0 -- the real, unscaled
        quantum Cv). Non-default values are mainly useful for
        diagnostics; physically you want xi=1.0 here.

    Returns
    -------
    cv_curve : ndarray, shape (len(beta_arr),)
        Cv/k_B at each beta.
    """
    return np.array([compute_cv(energies, b, xi) for b in beta_arr])


# =====================================================================
# Diagnostic plot: how Cv(xi) behaved at the hardest-to-converge T
# =====================================================================
def plot_xi_convergence_diagnostic(xi_result, beta_val, T_K_val, tol_xi, system_name):
    """
    Plot Cv as a function of the scaling factor xi for a single
    (hardest-converging) temperature. Each point is Cv evaluated with the
    spectrum of xi^2 V at temperature xi^2 T (see
    Classical_Limit_Numerical.py). Points are colored by whether they
    belong to the verified plateau, satisfy the per-step stability test
    without (yet) being part of it, or are still moving.

    Parameters
    ----------
    xi_result : dict
        Output of `converge_xi` for this temperature.
    beta_val : float
        The inverse temperature this diagnostic corresponds to.
    T_K_val : float
        The corresponding temperature (1/beta_val), used only for the title.
    tol_xi : float
        The tolerance used during xi-convergence (for coloring/legend).
    system_name : str
        Used in the plot title.

    Returns
    -------
    None (saves the figure to disk under figures/<system>/<params>/<category>/; see figures/output_paths.py).
    """
    BLUE, GREEN, ORANGE, YELLOW, GRAY = "#1f77b4", "#2ca02c", "#d62728", "#bcbd22", "#7f7f7f"
    xis, cvs, errs = xi_result["xi_values"], xi_result["cv_values"], xi_result["errors"]
    plateau = xi_result["plateau"]
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle(f"{system_name} \u2014 \u03be-Convergence Diagnostic (\u03be\u00b2V re-solved)\nHardest T: {T_K_val:.3g}  (\u03b2 = {beta_val:.4f})", fontsize=12, fontweight="bold")
    ax.plot(xis, cvs, color=BLUE, linewidth=1.5, marker="s", markersize=5, zorder=3, label="Cv(\u03be)")
    for i in range(len(xis)):
        if plateau is not None and plateau[0] <= i <= plateau[1]:
            c = GREEN
        elif errs[i] is not None and errs[i] < tol_xi and cvs[i] >= 0.5 - tol_xi:
            c = YELLOW
        else:
            c = GRAY
        ax.scatter([xis[i]], [cvs[i]], color=c, zorder=5, s=60)
    if xi_result["converged"]:
        xc, cc, ec = xi_result["xi_converged"], xi_result["cv_converged"], xi_result["error_estimate"]
        ax.axvline(xc, color=GREEN, linestyle=":", linewidth=1.3, label=f"\u03be_conv = {xc:.4g}")
        ax.scatter([xc], [cc], color=GREEN, zorder=6, s=100, label=f"Cv_conv = {cc:.4f}")
        ann, ann_colour = (f"Converged \u2713\n\u03be_conv = {xc:.4g}\nCv_conv/kB = {cc:.5f}\n"
                           f"est. distance to limit {ec:.1e}"), GREEN
    else:
        ann, ann_colour = f"NOT converged\n({xi_result['stop_reason']})", ORANGE
    ax.text(0.97, 0.03, ann, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=ann_colour,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=ann_colour, alpha=0.9))
    dot_legend = [mpatches.Patch(color=GREEN, label="verified plateau (reported value = last point)"),
                  mpatches.Patch(color=YELLOW, label="stable step, not (yet) a full plateau"),
                  mpatches.Patch(color=GRAY, label="still moving / first point")]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + dot_legend, fontsize=9, loc="center left")
    ax.set_xscale("log")
    ax.set_xlabel("Scaling factor  \u03be   (V \u2192 \u03be\u00b2V,  T \u2192 \u03be\u00b2T)", fontsize=11)
    ax.set_ylabel("Cv / kB", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    save_figure(fig, "convergence", "xi_convergence")


# =====================================================================
# Diagnostic plot: how Cv(n) behaved at the hardest-to-converge T
# =====================================================================
def plot_n_convergence_diagnostic(n_result, beta_val, T_K_val, tol_cv, system_name):
    """
    Plot Cv as a function of the number of included energy levels n
    for a single (hardest-converging) temperature, marking where (if
    anywhere) it became stable for the remainder of the spectrum.

    Parameters
    ----------
    n_result : dict
        Output of `converge_n` for this temperature.
    beta_val : float
        The inverse temperature this diagnostic corresponds to.
    T_K_val : float
        The corresponding temperature (1/beta_val), used only for the title.
    tol_cv : float
        The tolerance used during n-convergence (for annotation only).
    system_name : str
        Used in the plot title.

    Returns
    -------
    None (saves the figure to disk under figures/<system>/<params>/<category>/; see figures/output_paths.py).
    """
    BLUE, GREEN, ORANGE = "#1f77b4", "#2ca02c", "#d62728"
    ns, cvs = n_result["n_values"], n_result["cv_values"]
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle(f"{system_name} \u2014 n-Convergence Diagnostic\nHardest T: {T_K_val:.2f}  (\u03b2 = {beta_val:.4f})", fontsize=12, fontweight="bold")
    ax.plot(ns, cvs, color=BLUE, linewidth=1.5, marker="o", markersize=3, zorder=3, label="Cv(n levels)")
    if n_result["converged"]:
        nc = n_result["n_converged"]
        idx = ns.index(nc)
        ax.axvline(nc, color=GREEN, linestyle=":", linewidth=1.3, label=f"Converged at n = {nc}")
        ax.scatter([nc], [cvs[idx]], color=GREEN, zorder=6, s=100)
        ann, ann_colour = f"Converged \u2713  at n = {nc}\nCv/kB = {cvs[idx]:.5f}", GREEN
    else:
        ann, ann_colour = "NOT converged within N levels\nIncrease N_MAX", ORANGE
    ax.text(0.97, 0.05, ann, transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=ann_colour,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=ann_colour, alpha=0.9))
    ax.set_xlabel("Number of energy levels  n", fontsize=11)
    ax.set_ylabel("Cv / kB", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    save_figure(fig, "convergence", "n_convergence")


# =====================================================================
# Summary plot: quantum Cv(T), numerical classical limit, optional analytic overlay
# =====================================================================
def plot_cv_curves(T_arr, cv_quantum, cv_classical, xi_conv_arr, n_conv_arr, system_name,
                    cv_analytic_classical=None, T_units_label=r"$k_B T \,/\, E_0$"):
    """
    Plot the true quantum Cv(T) curve together with the numerically
    found classical-limit Cv(T) curve, plus (optionally) a reference
    analytic classical-limit curve/value. A secondary y-axis shows
    the converged xi and n values across temperature, which is a
    useful at-a-glance indicator of how hard convergence was at each T.

    Parameters
    ----------
    T_arr : ndarray
        Temperature values (= 1/beta_arr).
    cv_quantum : ndarray
        True quantum Cv(T), from `compute_quantum_heat_capacity_curve`.
    cv_classical : ndarray
        Numerical classical-limit Cv(T), from `sweep_temperature_range`.
    xi_conv_arr, n_conv_arr : ndarray
        Converged xi/n at each temperature (NaN where convergence failed).
    system_name : str
        Used in the plot title.
    cv_analytic_classical : float, ndarray, or None, optional
        If provided, overlaid as a reference "Analytic classical limit"
        curve (constant if a scalar is given).
    T_units_label : str, optional
        X-axis label (LaTeX-formatted), defaults to a generic
        dimensionless temperature label.

    Returns
    -------
    None (saves the figure to disk under figures/<system>/<params>/<category>/; see figures/output_paths.py).
    """
    BLUE, GREEN, ORANGE, PURPLE, RED = "#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#d62728"
    fig, ax1 = plt.subplots(figsize=(9, 6))
    fig.suptitle(f"{system_name} \u2014 Cv(T)", fontsize=13, fontweight="bold")
    ax1.plot(T_arr, cv_quantum, color=BLUE, linewidth=2, label="Quantum Cv(T)")
    ax1.plot(T_arr, cv_classical, color=GREEN, linewidth=2, linestyle="--", label="Numerical classical limit")
    if cv_analytic_classical is not None:
        cv_ref = np.full_like(T_arr, cv_analytic_classical) if np.isscalar(cv_analytic_classical) else cv_analytic_classical
        ax1.plot(T_arr, cv_ref, color=ORANGE, linewidth=1.5, linestyle=":", label="Analytic classical limit")
    ax1.set_xlabel(T_units_label, fontsize=12)
    ax1.set_ylabel(r"$C_v \,/\, k_B$", fontsize=12)
    ax1.set_xscale("log")
    ax1.legend(fontsize=10, loc="upper left")
    ax1.grid(True, linestyle="--", alpha=0.4)
    ax2 = ax1.twinx()
    valid_xi = ~np.isnan(xi_conv_arr)
    valid_n = ~np.isnan(n_conv_arr)
    ax2.plot(T_arr[valid_xi], xi_conv_arr[valid_xi], color=PURPLE, linewidth=1, linestyle="-.", alpha=0.6, label="\u03be_conv(T)")
    ax2.plot(T_arr[valid_n], n_conv_arr[valid_n], color=RED, linewidth=1, linestyle=":", alpha=0.6, label="n_conv(T)")
    ax2.set_ylabel("Converged \u03be  /  n  (secondary axis)", fontsize=10, color=PURPLE)
    ax2.tick_params(axis="y", colors=PURPLE)
    ax2.legend(fontsize=9, loc="upper right")
    plt.tight_layout()
    save_figure(fig, "cv", "cv_summary")


# =====================================================================
# Full pipeline: quantum Cv + numerical classical limit + all plots
# =====================================================================
def run(energies, potential_func, system_name,
        beta_min=0.02, beta_max=5.0, n_beta=200,
        xi_start=1.0, tol_xi=2e-3, min_stable_xi=3, xi_multiplier=1.25, max_xi_steps=30,
        tol_cv=1e-4, min_stable_n=3,
        mass=1.0, hbar=1.0, thermal_coverage=20.0, xi_max=np.inf, spectra=None,
        cv_analytic=None, T_units_label=r"$k_B T \,/\, E_0$"):
    """
    Run the full general-purpose Cv pipeline for ANY system: sweep the
    temperature range, find the numerical classical limit at every T
    (via Classical_Limit_Numerical, which re-solves xi^2 V at every xi),
    compute the true quantum Cv(T) curve from `energies`, and produce
    the xi-convergence diagnostic, n-convergence diagnostic (each shown
    at the single hardest-to-converge temperature), and the combined
    Cv(T) summary plot.

    Parameters
    ----------
    energies : array_like
        Energy eigenvalues for the system (e.g. from DVR), ascending.
        Used for the quantum Cv(T) only.
    potential_func : callable
        The system's potential V(x), unscaled. Used for the classical
        limit only (each xi step solves xi^2 V).
    system_name : str
        Human-readable system name, used in plot titles/console output.
    beta_min, beta_max, n_beta : float, float, int, optional
        Inverse-temperature sweep range and point count. `beta_arr` is
        LOG-spaced (`np.geomspace`), not linear: every plot in this
        project displays T on a log axis, and beta_min/beta_max
        typically span more than a decade, so linear spacing in beta
        leaves the high-T tail (small beta) drastically under-sampled
        relative to how it's actually viewed -- e.g. for a typical
        beta_min~0.075/beta_max=50/n_beta=1000, linear spacing puts
        just ~3 points across the T=5-13 decade (visibly polygonal,
        "connected straight lines" instead of a smooth curve, most
        noticeable wherever a curve's transition happens to land in
        that stretch), where log spacing puts ~150.
    xi_start, tol_xi, min_stable_xi, xi_multiplier, max_xi_steps :
        The xi ladder and plateau criterion, passed through to
        `sweep_temperature_range` / `converge_xi`.
    tol_cv, min_stable_n :
        Passed through to `converge_n` at every temperature.
    mass, hbar : float, optional
        Physical constants of the system (default 1.0).
    thermal_coverage : float, optional
        Each scaled spectrum keeps levels up to E_0 + thermal_coverage
        * k_B T (default 20.0; normally config.HOT_STATE_SAFETY).
    xi_max : float, optional
        Hard cap on xi (normally config.XI_MAX).
    spectra : ScaledSpectra or None, optional
        Cache of scaled spectra to reuse (e.g. from the previous
        auto-tune round); None builds a new one. Returned in the result.
    cv_analytic : float, ndarray, or None, optional
        If known, the analytic classical-limit Cv to overlay as a
        reference curve on the summary plot.
    T_units_label : str, optional
        X-axis label for the summary plot.

    Returns
    -------
    dict with keys:
        beta_arr, T_arr : ndarray
        cv_quantum, cv_classical : ndarray
        xi_conv, n_conv : ndarray
        sweep : dict
            Full output of `sweep_temperature_range` (includes
            per-temperature convergence traces for further inspection).
        spectra : ScaledSpectra
            The scaled-spectrum cache, for reuse by a later call.
    """
    beta_arr = np.geomspace(beta_min, beta_max, n_beta)
    T_arr = 1.0 / beta_arr

    rule = "\u2550" * 60
    print(f"\n{rule}\n  {system_name}\n{rule}")
    print(f"  {len(energies)} levels, E_min={energies[0]:.3g}, E_max={energies[-1]:.3g}")
    print(f"  \u03b2: {beta_min} \u2192 {beta_max}  ({n_beta} log-spaced points)")

    sweep = sweep_temperature_range(
        potential_func, beta_arr,
        xi_start, tol_xi, min_stable_xi, xi_multiplier, max_xi_steps,
        tol_cv, min_stable_n, mass=mass, hbar=hbar,
        thermal_coverage=thermal_coverage, xi_max=xi_max, spectra=spectra, verbose=True,
    )
    cv_classical = sweep["cv_classical"]
    xi_conv = sweep["xi_conv"]
    n_conv = sweep["n_conv"]

    # The quantum Cv uses the whole base spectrum. (It used to be cut at
    # the largest n the classical n-scan needed, but that n now refers to
    # the scaled spectra of xi^2 V, not to `energies`.) This is exact as
    # long as the top level is thermally inaccessible, which the auto-tune
    # loop enforces (Cv_AutoTune: E_max >= HOT_STATE_SAFETY * T_hot).
    cv_quantum = compute_quantum_heat_capacity_curve(energies, beta_arr, xi=1.0)

    valid_xi_mask = ~np.isnan(xi_conv)
    if valid_xi_mask.any():
        idx_hard_xi = int(np.nanargmax(xi_conv))
        print(f"  Hardest \u03be-convergence: T*={T_arr[idx_hard_xi]:.4g}, \u03be_conv={xi_conv[idx_hard_xi]:.3f}")
        plot_xi_convergence_diagnostic(sweep["xi_results"][idx_hard_xi], beta_arr[idx_hard_xi], T_arr[idx_hard_xi], tol_xi, system_name)

    valid_n_mask = ~np.isnan(n_conv)
    if valid_n_mask.any():
        idx_hard_n = int(np.nanargmax(n_conv))
        print(f"  Hardest n-convergence:  T*={T_arr[idx_hard_n]:.4g}, n_conv={int(n_conv[idx_hard_n])}")
        nr_hard = sweep["n_results"][idx_hard_n]
        if nr_hard is not None:
            plot_n_convergence_diagnostic(nr_hard, beta_arr[idx_hard_n], T_arr[idx_hard_n], tol_cv, system_name)

    plot_cv_curves(T_arr, cv_quantum, cv_classical, xi_conv, n_conv, system_name, cv_analytic_classical=cv_analytic, T_units_label=T_units_label)

    print(f"\n  \u03be-conv: {valid_xi_mask.sum()}/{n_beta}  (max \u03be={np.nanmax(xi_conv):.2f})" if valid_xi_mask.any() else "  \u03be-conv: not applicable")
    print(f"  n-conv: {valid_n_mask.sum()}/{n_beta}  (max n={int(np.nanmax(n_conv))})" if valid_n_mask.any() else "  n-conv: failed at all T")
    print(f"{rule}\n")

    return {"beta_arr": beta_arr, "T_arr": T_arr, "cv_quantum": cv_quantum, "cv_classical": cv_classical,
            "xi_conv": xi_conv, "n_conv": n_conv, "sweep": sweep, "spectra": sweep["spectra"]}
