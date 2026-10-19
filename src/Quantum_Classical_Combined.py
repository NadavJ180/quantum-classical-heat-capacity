"""
Quantum_Classical_Combined.py
Quantum and classical Cv(T) of ONE potential, plus the checks that say where
each curve can be trusted. Both pipelines (Quick_Scan.py, Quantum_HO_Master.py)
call `run` here; they only pass in different settings.
"""

import time

import numpy as np
import matplotlib.pyplot as plt

import config
from DVR.DVR_Algorithm import solve_energy_levels, number_of_grid_points
from DVR.DVR_Reference_Generator import reference_grid
from Classical_Limit_Numerical import (heat_capacity_curve, classical_limit_curves,
                                       classical_error_estimate, cv_vs_number_of_levels)
from figures.output_paths import save_figure


def potential_label(params):
    """Text such as 'double_well(a=0.25, b=-0.5, c=-0.5, d=0)' for titles."""
    values = ", ".join(f"{k}={v:g}" for k, v in params.items())
    return f"{config.POTENTIAL.__name__}({values})"


def run(potential_func, dx, num_states, xi_values, n_T, cv_tol):
    """Quantum Cv, classical Cv and trust checks of one potential on the config grid."""
    t_start = time.time()
    x_min, x_max, walls = config.X_MIN, config.X_MAX, config.WALLS
    mass, hbar = config.MASS, config.HBAR
    T_arr = np.geomspace(config.T_MIN, config.T_MAX, n_T)

    # --- Quantum Cv from the levels of V(x) ---
    # The extra levels (REFERENCE_LEVELS_FACTOR) are only used by the "more levels" check.
    more_states = int(round(num_states * config.REFERENCE_LEVELS_FACTOR))
    more_states = min(more_states, number_of_grid_points(x_min, x_max, dx, walls) - 1)
    print(f"  Quantum levels: {num_states} on [{x_min:g}, {x_max:g}], dx = {dx:g}, walls = {walls}",
          flush=True)
    energies_all = solve_energy_levels(potential_func, more_states, x_min, x_max, dx, mass, hbar, walls)
    energies = energies_all[:num_states]
    cv_quantum = heat_capacity_curve(energies, T_arr)

    # --- Trust checks: recompute the quantum Cv with ONE knob changed at a time ---
    checks = {}
    E = solve_energy_levels(potential_func, num_states, x_min, x_max,
                            dx / config.REFERENCE_DX_FACTOR, mass, hbar, walls)
    checks[f"finer grid (dx -> {dx / config.REFERENCE_DX_FACTOR:.3g})"] = \
        np.abs(heat_capacity_curve(E, T_arr) - cv_quantum)
    if not walls:
        x_lo, x_hi, _ = reference_grid(x_min, x_max, dx, config.REFERENCE_SPAN_FACTOR, 1.0)
        E = solve_energy_levels(potential_func, num_states, x_lo, x_hi, dx, mass, hbar, walls)
        checks[f"wider grid ([{x_lo:.3g}, {x_hi:.3g}])"] = \
            np.abs(heat_capacity_curve(E, T_arr) - cv_quantum)
    checks[f"more levels ({num_states} -> {more_states})"] = \
        np.abs(heat_capacity_curve(energies_all, T_arr) - cv_quantum)
    quantum_error = np.max(np.array(list(checks.values())), axis=0)

    # --- Classical limit: the same V(x) with hbar -> hbar/xi ---
    print("  Classical limit (hbar -> hbar/xi):", flush=True)
    curves = classical_limit_curves(potential_func, T_arr, xi_values, num_states, x_min, x_max, dx,
                                    mass, hbar, walls)
    cv_classical = curves[xi_values[-1]]
    classical_error = classical_error_estimate(curves, xi_values)

    # --- Where each curve can be trusted, and quantum minus classical ---
    # Each xi solve has the same points per wavelength and energy window as xi = 1,
    # so the classical curve also needs the quantum checks to pass.
    quantum_ok = quantum_error < cv_tol
    with np.errstate(invalid="ignore"):
        classical_ok = quantum_ok & (classical_error < cv_tol * cv_classical)
    difference = cv_quantum - cv_classical
    uncertainty = quantum_error + classical_error
    resolved = np.abs(difference) > config.ERROR_MARGIN * uncertainty
    above = classical_ok & resolved & (difference > 0)
    below = classical_ok & resolved & (difference < 0)
    undecided = ~above & ~below

    print(f"  Done in {time.time() - t_start:.1f} s", flush=True)
    return {
        "T_arr": T_arr, "energies": energies, "energies_all": energies_all,
        "cv_quantum": cv_quantum, "cv_classical": cv_classical, "curves": curves,
        "xi_values": list(xi_values), "checks": checks,
        "quantum_error": quantum_error, "classical_error": classical_error,
        "quantum_ok": quantum_ok, "classical_ok": classical_ok,
        "difference": difference, "uncertainty": uncertainty,
        "above": above, "below": below, "undecided": undecided,
        "dx": dx, "num_states": num_states, "cv_tol": cv_tol, "seconds": time.time() - t_start,
    }


def describe_T(T_arr, mask):
    """The temperatures where mask is True, as text like 'T = 0.05-3.2, 5-7.1'."""
    if not mask.any():
        return "none"
    if mask.all():
        return "all T"
    parts = []
    i = 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j + 1 < len(mask) and mask[j + 1]:
                j += 1
            parts.append(f"{T_arr[i]:.3g}-{T_arr[j]:.3g}" if j > i else f"{T_arr[i]:.3g}")
            i = j + 1
        else:
            i += 1
    return "T = " + ", ".join(parts)


def print_report(result):
    """Console summary: where each curve is trusted, and how quantum compares to classical."""
    T = result["T_arr"]
    print(f"\n  Trust checks (Cv recomputed with one knob changed; OK where it moves less than "
          f"cv_tol = {result['cv_tol']:g}):")
    for name, change in result["checks"].items():
        print(f"    {name:<40}: OK at {describe_T(T, change < result['cv_tol'])}"
              f"   (largest change {np.max(change):.1e})")
    print(f"  Quantum Cv trusted at   : {describe_T(T, result['quantum_ok'])}")
    xis = ", ".join(f"{xi:g}" for xi in result["xi_values"])
    print(f"  Classical Cv converged at: {describe_T(T, result['classical_ok'])}"
          f"   (xi = {xis}; estimated distance to hbar -> 0 below cv_tol * Cv)")
    print(f"  Quantum minus classical (resolved = more than {config.ERROR_MARGIN:g} x its uncertainty):")
    print(f"    ABOVE classical (trusted, resolved): {describe_T(T, result['above'])}")
    print(f"    below classical (trusted, resolved): {describe_T(T, result['below'])}")
    print(f"    undecided (a curve not trusted, or the difference within its uncertainty):"
          f" {describe_T(T, result['undecided'])}")
    print(f"      of which Cv_q is above the current classical estimate at:"
          f" {describe_T(T, result['undecided'] & (result['difference'] > 0))}")
    if result["classical_ok"].any():
        d = np.where(result["classical_ok"], result["difference"], np.nan)
        i = int(np.nanargmax(d))
        print(f"    largest Cv_q - Cv_cl where both are trusted: {d[i]:+.3e} at T = {T[i]:.3g}"
              f" (uncertainty {result['uncertainty'][i]:.1e})")


def print_grid_suggestion(suggestion, num_states, dx):
    """Print DVR_Algorithm.suggest_grid next to the values chosen in config.py."""
    s = suggestion
    print(f"  Suggestion for k_B T up to {config.T_MAX:g}: levels up to E = V_min + "
          f"{config.THERMAL_CUTOFF:g} k_B T_MAX = {s['E_top']:.4g} matter.")
    edge = "  <- V is still below E at a grid end" if s["touches_edge"] and not config.WALLS else ""
    print(f"    region where V < E      : x = {s['region'][0]:.3g} .. {s['region'][1]:.3g}"
          f"   (grid: {config.X_MIN:g} .. {config.X_MAX:g}){edge}")
    print(f"    levels below E (WKB)    : ~{s['num_levels']:.0f}   (num_states: {num_states})")
    print(f"    dx for 4 points per wavelength: {s['dx']:.3g}   (dx: {dx:g})")


def shade(ax, T_arr, mask, color, alpha, label=None):
    """Shade the temperature ranges where mask is True."""
    edges = np.concatenate(([T_arr[0]], np.sqrt(T_arr[1:] * T_arr[:-1]), [T_arr[-1]]))
    for i in range(len(T_arr)):
        if mask[i]:
            ax.axvspan(edges[i], edges[i + 1], color=color, alpha=alpha, lw=0, label=label)
            label = None


def plot_cv_summary(result, title, category, name):
    """Three panels: Cv(T) curves, Cv_q - Cv_cl, and the trust-check errors."""
    T = result["T_arr"]
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 11), sharex=True,
                                        gridspec_kw={"height_ratios": [2.2, 1.2, 1.4]})
    fig.suptitle(title, fontsize=12, fontweight="bold")

    greens = plt.get_cmap("Greens")
    xis = result["xi_values"]
    for k, xi in enumerate(xis[:-1]):
        ax1.plot(T, result["curves"][xi], color=greens(0.3 + 0.5 * k / max(len(xis) - 1, 1)),
                 lw=1, label=f"$\\hbar/\\xi$, $\\xi$ = {xi:g}")
    ax1.plot(T, result["cv_classical"], color="#2ca02c", lw=2.2, ls="--",
             label=f"classical estimate ($\\xi$ = {xis[-1]:g})")
    ax1.plot(T, result["cv_quantum"], color="#1f77b4", lw=2.2, label="quantum")
    shade(ax1, T, ~result["quantum_ok"], "gray", 0.25, "quantum not trusted")
    shade(ax1, T, result["quantum_ok"] & ~result["classical_ok"], "gold", 0.2, "classical not converged")
    shade(ax1, T, result["above"], "red", 0.3, "quantum ABOVE classical")
    ax1.set_ylabel(r"$C_v / k_B$")
    ax1.legend(fontsize=8, loc="upper left")

    d, band = result["difference"], config.ERROR_MARGIN * result["uncertainty"]
    ax2.plot(T, d, color="k", lw=1.5)
    ax2.fill_between(T, d - band, d + band, color="gray", alpha=0.3,
                     label=f"$\\pm${config.ERROR_MARGIN:g} x uncertainty")
    ax2.axhline(0, color="r", lw=1)
    ax2.set_ylabel(r"$C_{v,q} - C_{v,cl}$")
    ax2.legend(fontsize=8, loc="upper left")

    for (label, change), color in zip(result["checks"].items(), ["#1f77b4", "#ff7f0e", "#9467bd"]):
        ax3.plot(T, np.where(change > 0, change, np.nan), color=color, lw=1.3, label=f"quantum: {label}")
    with np.errstate(invalid="ignore", divide="ignore"):
        relative = result["classical_error"] / result["cv_classical"]
    ax3.plot(T, np.where(relative > 0, relative, np.nan), color="#2ca02c", lw=1.6, ls="--",
             label="classical: estimated distance to $\\hbar\\to 0$, / $C_v$")
    ax3.axhline(result["cv_tol"], color="r", ls=":", label=f"cv_tol = {result['cv_tol']:g}")
    ax3.set_yscale("log")
    ax3.set_ylabel("change in $C_v$")
    ax3.set_xlabel(config.T_UNITS_LABEL)
    ax3.legend(fontsize=7, loc="upper left")

    for ax in (ax1, ax2, ax3):
        ax.set_xscale("log")
        ax.grid(True, ls="--", alpha=0.4)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return save_figure(fig, category, name)


def plot_xi_convergence(result, n_temperatures=6):
    """Cv against 1/xi^2 at a few temperatures. The leading correction is of order
    hbar^2, so near 1/xi^2 = 0 the points lie on a straight line whose end is the
    classical Cv."""
    T, xis = result["T_arr"], result["xi_values"]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = plt.get_cmap("viridis")
    for k, i in enumerate(np.linspace(0, len(T) - 1, n_temperatures).astype(int)):
        x = [1.0 / xi**2 for xi in xis]
        y = [result["curves"][xi][i] for xi in xis]
        color = colors(k / max(n_temperatures - 1, 1))
        ax.plot(x, y, "o-", color=color, label=f"T = {T[i]:.3g}")
        if len(xis) >= 2:   # straight line through the last two points, down to 1/xi^2 = 0
            slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
            ax.plot([0, x[-1]], [y[-1] - slope * x[-1], y[-1]], ":", color=color)
    ax.set_xlabel(r"$1/\xi^2$   ($\hbar \to \hbar/\xi$;  0 = classical)")
    ax.set_ylabel(r"$C_v / k_B$")
    ax.set_title("Approach to the classical limit (dotted: $\\hbar^2$ law through the last two points)",
                 fontsize=10)
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend(fontsize=8)
    plt.tight_layout()
    return save_figure(fig, "convergence", "xi_convergence")


def plot_n_convergence(result, n_temperatures=6):
    """Cv against the number of levels kept, at a few temperatures. Flat = enough levels."""
    T = result["T_arr"]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = plt.get_cmap("viridis")
    for k, i in enumerate(np.linspace(0, len(T) - 1, n_temperatures).astype(int)):
        n_values, cv_values = cv_vs_number_of_levels(result["energies_all"], T[i])
        ax.plot(n_values, cv_values, color=colors(k / max(n_temperatures - 1, 1)), label=f"T = {T[i]:.3g}")
    ax.axvline(result["num_states"], color="r", ls=":", label=f"num_states = {result['num_states']}")
    ax.set_xlabel("number of levels kept, n")
    ax.set_ylabel(r"$C_v / k_B$")
    ax.set_xscale("log")
    ax.set_title("Quantum Cv against the number of levels kept", fontsize=10)
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend(fontsize=8)
    plt.tight_layout()
    return save_figure(fig, "convergence", "n_convergence")
