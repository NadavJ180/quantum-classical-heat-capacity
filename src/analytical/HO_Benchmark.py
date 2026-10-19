"""
HO_Benchmark.py
The numerical Cv(T) of the harmonic oscillator against the exact results
(HO_Analytical.py). The full pipeline runs it when config.POTENTIAL is
potentials.harmonic.
"""

import numpy as np
import matplotlib.pyplot as plt

from analytical.HO_Analytical import analytic_cv_HO_quantum, analytic_cv_HO_classical
from figures.output_paths import save_figure


def run_ho_benchmark(numeric_results, hbar=1.0, omega=1.0, kB=1.0):
    """numeric_results: Quantum_Classical_Combined.run output for the harmonic oscillator."""
    T = numeric_results["T_arr"]
    cv_q_exact = analytic_cv_HO_quantum(T, hbar, omega, kB)
    cv_cl_exact = analytic_cv_HO_classical(kB)
    error_q = np.abs(numeric_results["cv_quantum"] - cv_q_exact)
    error_cl = np.abs(numeric_results["cv_classical"] - cv_cl_exact)
    trusted = numeric_results["classical_ok"]
    print(f"  Harmonic oscillator vs exact: max |Cv_q error| = {np.max(error_q):.2e}")
    print(f"    max |Cv_cl error| = {np.max(error_cl):.2e} at all T,"
          f" {np.max(error_cl[trusted]) if trusted.any() else np.nan:.2e} where the classical Cv is trusted")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                   gridspec_kw={"height_ratios": [2.2, 1]})
    fig.suptitle("Harmonic oscillator: numerical vs exact", fontsize=13, fontweight="bold")
    ax1.plot(T, numeric_results["cv_quantum"], color="#1f77b4", lw=2, label="quantum, DVR")
    ax1.plot(T, cv_q_exact, color="#d62728", lw=1.5, ls="--", label="quantum, Einstein formula")
    ax1.plot(T, numeric_results["cv_classical"], color="#2ca02c", lw=1.8, ls="-.", label="classical, hbar/xi")
    ax1.axhline(cv_cl_exact, color="#9467bd", ls=":", label="classical, exact ($k_B$)")
    ax1.set_ylabel(r"$C_v / k_B$")
    ax1.legend(fontsize=9)
    ax2.plot(T, np.where(error_q > 0, error_q, np.nan), color="#1f77b4", label="quantum")
    ax2.plot(T, np.where(error_cl > 0, error_cl, np.nan), color="#2ca02c", label="classical")
    ax2.set_yscale("log")
    ax2.set_ylabel("|numerical - exact|")
    ax2.set_xlabel(r"$k_B T / \hbar\omega$")
    ax2.legend(fontsize=9)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.grid(True, ls="--", alpha=0.4)
    plt.tight_layout()
    save_figure(fig, "cv", "analytic_cv_benchmark")
    return {"error_quantum": error_q, "error_classical": error_cl}
