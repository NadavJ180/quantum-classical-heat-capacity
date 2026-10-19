"""
error_energylevels.py
Compares two spectra level by level (base DVR vs reference DVR).
The error is an ENERGY difference |E_n - E_n,ref|: Boltzmann factors depend on
(E_n - E_0) / k_B T, so an energy error matters compared with k_B T, not with E_n
itself (E_n can be negative or zero, e.g. in a double well).
"""

import numpy as np
import matplotlib.pyplot as plt

from figures.output_paths import save_figure


def compute_energy_level_errors(E_numeric, E_reference):
    """|E_numeric[n] - E_reference[n]| for every level n, plus the worst level."""
    E_numeric = np.asarray(E_numeric, dtype=float)
    E_reference = np.asarray(E_reference, dtype=float)
    abs_error = np.abs(E_numeric - E_reference)
    worst = int(np.argmax(abs_error))
    return {"abs_error": abs_error, "max_abs_error": float(abs_error[worst]), "max_abs_idx": worst,
            "mean_abs_error": float(np.mean(abs_error))}


def plot_energy_level_comparison(E_numeric, E_reference, error_result, system_name, tolerance=None):
    """Top: E_n of both spectra. Bottom: |E_n - E_n,ref| against n (log scale)."""
    n = np.arange(len(E_numeric))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    fig.suptitle(f"{system_name}\nbase DVR vs reference DVR", fontsize=12, fontweight="bold")
    ax1.plot(n, E_reference, color="#d62728", lw=2, label="reference")
    ax1.plot(n, E_numeric, color="#1f77b4", lw=1, ls="--", label="base")
    ax1.set_ylabel("$E_n$")
    ax1.legend()
    error = error_result["abs_error"]
    ax2.plot(n, np.where(error > 0, error, np.nan), color="#1f77b4", lw=1.2)
    if tolerance is not None:
        ax2.axhline(tolerance, color="r", ls=":", label=f"tolerance = {tolerance:g}")
        ax2.legend()
    ax2.set_yscale("log")
    ax2.set_xlabel("level index n")
    ax2.set_ylabel(r"$|E_n - E_{n,\mathrm{ref}}|$")
    for ax in (ax1, ax2):
        ax.grid(True, ls="--", alpha=0.4)
    plt.tight_layout()
    return save_figure(fig, "energy_levels", "energy_level_comparison")


def print_accuracy_summary(error_result, system_name):
    print(f"  {system_name}: {len(error_result['abs_error'])} levels compared with the reference")
    print(f"    mean |dE| = {error_result['mean_abs_error']:.3e}")
    print(f"    max  |dE| = {error_result['max_abs_error']:.3e}  (at n = {error_result['max_abs_idx']})")
