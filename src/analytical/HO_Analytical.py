"""
HO_Analytical.py
Exact results for the 1-D harmonic oscillator, used to test the numerics.
"""

import numpy as np


def analytic_energy_levels_HO(num_levels, hbar=1.0, omega=1.0):
    """E_n = hbar omega (n + 1/2), n = 0 .. num_levels-1."""
    return hbar * omega * (np.arange(num_levels) + 0.5)


def analytic_cv_HO_quantum(T, hbar=1.0, omega=1.0, kB=1.0):
    """Einstein formula: Cv = k_B x^2 e^x / (e^x - 1)^2, with x = hbar omega / (k_B T)."""
    x = hbar * omega / (kB * np.asarray(T, dtype=float))
    return kB * x**2 * np.exp(x) / (np.exp(x) - 1.0)**2


def analytic_cv_HO_classical(kB=1.0):
    """Equipartition: k_B/2 from p^2/2m plus k_B/2 from m omega^2 x^2/2, at every T."""
    return kB
