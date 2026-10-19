"""
Cv_Coefficient_Sweep.py
The same calculation for several values of one potential parameter
(config.SCAN_PARAM takes every value in config.SCAN_VALUES, everything else
fixed). Each variant's quantum Cv is compared with ITS OWN classical Cv: the
classical Cv depends on the potential too.
"""

import numpy as np
import matplotlib.pyplot as plt

import config
from Quantum_Classical_Combined import run, print_report, describe_T, potential_label
from figures.output_paths import save_figure


def variant_params(base_params, scan_param, values):
    """One copy of base_params per value, with scan_param set to that value."""
    variants = []
    for value in values:
        params = dict(base_params)
        params[scan_param] = value
        variants.append(params)
    return variants


def variant_colors(n):
    cmap = plt.get_cmap("tab10") if n <= 10 else plt.get_cmap("viridis")
    return [cmap(i if n <= 10 else i / max(n - 1, 1)) for i in range(n)]


def plot_coefficient_sweep(results, values, title, category, name):
    """Top: quantum (solid) and classical (dashed) Cv of every variant.
    Bottom: Cv_q - Cv_cl, solid where both curves are trusted, dotted elsewhere."""
    T = results[0]["T_arr"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1.2]})
    fig.suptitle(title, fontsize=12, fontweight="bold")
    for r, value, color in zip(results, values, variant_colors(len(values))):
        ax1.plot(T, r["cv_quantum"], color=color, lw=1.8, label=f"{config.SCAN_PARAM} = {value:g}")
        ax1.plot(T, r["cv_classical"], color=color, lw=1.5, ls="--")
        ax2.plot(T, np.where(r["classical_ok"], r["difference"], np.nan), color=color, lw=1.8)
        ax2.plot(T, r["difference"], color=color, lw=0.8, ls=":")
    ax1.plot([], [], color="0.3", label="quantum")
    ax1.plot([], [], color="0.3", ls="--", label="classical (same potential)")
    ax1.set_ylabel(r"$C_v / k_B$")
    ax1.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax2.axhline(0, color="r", lw=1)
    ax2.set_ylabel(r"$C_{v,q} - C_{v,cl}$" + "\n(dotted: not trusted)")
    ax2.set_xlabel(config.T_UNITS_LABEL)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.grid(True, ls="--", alpha=0.4)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return save_figure(fig, category, name)


def plot_variant_potentials(params_list, results, values, category, name, levels_to_draw=5):
    """V(x) of every variant, with its lowest levels drawn across the allowed region."""
    x = np.linspace(config.X_MIN, config.X_MAX, 4000)
    fig, ax = plt.subplots(figsize=(9, 6))
    y_low, y_high = np.inf, -np.inf
    for params, r, value, color in zip(params_list, results, values, variant_colors(len(values))):
        V = config.my_potential(x, params)
        ax.plot(x, V, color=color, lw=1.8, label=f"{config.SCAN_PARAM} = {value:g}")
        for E in r["energies"][:levels_to_draw]:
            inside = x[V <= E]
            if inside.size:
                ax.hlines(E, inside.min(), inside.max(), color=color, lw=0.8, alpha=0.7)
        y_low = min(y_low, V.min())
        y_high = max(y_high, r["energies"][min(levels_to_draw, len(r["energies"])) - 1])
    margin = 0.15 * (y_high - y_low)
    ax.set_ylim(y_low - margin, y_high + margin)
    shown = np.zeros(len(x), dtype=bool)               # zoom on x where some V is in view
    for params in params_list:
        shown |= config.my_potential(x, params) <= y_high + margin
    ax.set_xlim(x[shown].min(), x[shown].max())
    ax.set_xlabel("$x$")
    ax.set_ylabel("Energy")
    ax.set_title(f"{config.SYSTEM_NAME}: V(x) and lowest {levels_to_draw} levels", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, ls="--", alpha=0.4)
    plt.tight_layout()
    return save_figure(fig, category, name)


def run_coefficient_sweep(dx, num_states, xi_values, n_T, cv_tol, category="coefficient_sweep",
                          name="cv_coefficient_sweep"):
    """Run Quantum_Classical_Combined.run for every variant, then plot and summarise."""
    values = list(config.SCAN_VALUES)
    params_list = variant_params(config.POTENTIAL_PARAMS, config.SCAN_PARAM, values)
    results = []
    for k, params in enumerate(params_list):
        print(f"\n  [{k + 1}/{len(params_list)}] {potential_label(params)}")
        potential_func = lambda x, p=params: config.my_potential(x, p)
        result = run(potential_func, dx, num_states, xi_values, n_T, cv_tol)
        print_report(result)
        results.append(result)

    title = (f"{config.SYSTEM_NAME}: quantum vs classical $C_v$ across {config.SCAN_PARAM}\n"
             f"{potential_label(config.POTENTIAL_PARAMS)} with {config.SCAN_PARAM} varied,"
             f" dx = {dx:g}, {num_states} levels")
    cv_path = plot_coefficient_sweep(results, values, title, category, name)
    v_path = plot_variant_potentials(params_list, results, values, category, name + "_potentials")

    T = results[0]["T_arr"]
    print(f"\n  Sweep summary ({config.SCAN_PARAM}):")
    for value, r in zip(values, results):
        d = np.where(r["classical_ok"], r["difference"], np.nan)
        best = f"{np.nanmax(d):+.2e}" if np.isfinite(d).any() else "-"
        print(f"    {config.SCAN_PARAM} = {value:<8g} both trusted at {describe_T(T, r['classical_ok'])};"
              f" ABOVE at {describe_T(T, r['above'])}; largest Cv_q - Cv_cl (trusted) {best}")
    print(f"  Figures: {cv_path}\n           {v_path}")
    return {"values": values, "params": params_list, "results": results}
