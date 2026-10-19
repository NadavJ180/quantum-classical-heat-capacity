"""
Cv_Numerical_Benchmark.py
The quantum Cv(T) from the base spectrum against the same Cv(T) from the
reference spectrum (wider and finer grid, DVR_Reference_Generator.py). Optionally
the classical Cv too, recomputed with every xi solve on the reference grid.
The error shown is the absolute difference |Cv_base - Cv_ref| in units of k_B:
the quantum Cv goes to 0 at low T, where a relative error is meaningless.
"""

import numpy as np
import matplotlib.pyplot as plt

import config
from Classical_Limit_Numerical import heat_capacity_curve, classical_limit_curves
from figures.output_paths import save_figure


def plot_cv_comparison(T_arr, cv_base, cv_ref, curve_name, system_name, reference_text, name):
    """Top: both curves. Bottom: |Cv_base - Cv_ref| (log scale)."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                   gridspec_kw={"height_ratios": [2.2, 1]})
    fig.suptitle(f"{system_name}: {curve_name} Cv(T), base vs reference\n{reference_text}",
                 fontsize=12, fontweight="bold")
    ax1.plot(T_arr, cv_base, color="#1f77b4", lw=2, label="base grid")
    ax1.plot(T_arr, cv_ref, color="#d62728", lw=1.5, ls="--", label="reference grid")
    ax1.set_ylabel(r"$C_v / k_B$")
    ax1.legend()
    error = np.abs(cv_base - cv_ref)
    ax2.plot(T_arr, np.where(error > 0, error, np.nan), color="#1f77b4", lw=1.5)
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$|C_{v,\mathrm{base}} - C_{v,\mathrm{ref}}|$")
    ax2.set_xlabel(config.T_UNITS_LABEL)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.grid(True, ls="--", alpha=0.4)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    return save_figure(fig, "cv", name)


def run_cv_numerical_benchmark(result, reference, potential_func, system_name, run_classical=False):
    """result: Quantum_Classical_Combined.run output; reference: generate_reference_energies output."""
    T = result["T_arr"]
    reference_text = f"reference grid [{reference['x_min']:.3g}, {reference['x_max']:.3g}], dx = {reference['dx']:.3g}"
    cv_ref = heat_capacity_curve(reference["energies"][:result["num_states"]], T)
    error = np.abs(result["cv_quantum"] - cv_ref)
    print(f"  Quantum Cv:   max |Cv_base - Cv_ref| = {np.max(error):.3e}"
          f" (at T = {T[int(np.argmax(error))]:.3g})")
    plot_cv_comparison(T, result["cv_quantum"], cv_ref, "quantum", system_name, reference_text,
                       "quantum_cv_benchmark")
    if not run_classical:
        return {"cv_quantum_ref": cv_ref}

    print("  Classical Cv on the reference grid:")
    curves = classical_limit_curves(potential_func, T, result["xi_values"], result["num_states"],
                                    reference["x_min"], reference["x_max"], reference["dx"],
                                    config.MASS, config.HBAR, config.WALLS)
    cv_cl_ref = curves[result["xi_values"][-1]]
    error = np.abs(result["cv_classical"] - cv_cl_ref)
    print(f"  Classical Cv: max |Cv_base - Cv_ref| = {np.max(error):.3e}"
          f" (at T = {T[int(np.argmax(error))]:.3g})")
    plot_cv_comparison(T, result["cv_classical"], cv_cl_ref, "classical", system_name, reference_text,
                       "classical_cv_benchmark")
    return {"cv_quantum_ref": cv_ref, "cv_classical_ref": cv_cl_ref}
