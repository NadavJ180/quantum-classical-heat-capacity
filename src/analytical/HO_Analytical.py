"""
HO_Analytical.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Closed-form (analytic) reference solutions for the 1-D Quantum
Harmonic Oscillator (HO). These are exact formulas -- no DVR, no
diagonalization, no numerical convergence scanning -- used purely as
ground truth to benchmark the numerical (DVR-based) pipeline against.
This file is HO-only: it exists specifically to certify the pipeline's
numerics against a system with a known closed-form answer, one time,
rather than as part of the main (potential-agnostic) driver.

Provides:
    - Exact energy eigenvalues E_n = hbar*omega*(n + 1/2)
    - The exact quantum Cv(T) (the Einstein oscillator formula)
    - The exact classical (high-temperature) limit of Cv, which for a
      1-D harmonic oscillator is simply k_B (one quadratic kinetic +
      one quadratic potential degree of freedom, by equipartition)

This file does no plotting and makes no figures -- it just returns
numbers. Plotting lives in HO_Benchmark.py.
=====================================================================
"""

import numpy as np


# =====================================================================
# Exact HO energy eigenvalues
# =====================================================================
def analytic_energy_levels_HO(num_levels, hbar=1.0, omega=1.0):
    """
    Exact 1-D quantum harmonic oscillator energy eigenvalues,
    E_n = hbar * omega * (n + 1/2), for n = 0, 1, ..., num_levels-1.

    Parameters
    ----------
    num_levels : int
        Number of energy levels to return (n = 0 .. num_levels-1).
    hbar : float, optional
        Reduced Planck constant (default 1.0, dimensionless units).
    omega : float, optional
        Angular frequency of the oscillator (default 1.0).

    Returns
    -------
    E : ndarray, shape (num_levels,)
        Energy eigenvalues in ascending order, E[n] = hbar*omega*(n+1/2).
    """
    n = np.arange(num_levels, dtype=float)
    return hbar * omega * (n + 0.5)


# =====================================================================
# Exact HO quantum heat capacity (Einstein oscillator formula)
# =====================================================================
def analytic_cv_HO_quantum(T, hbar=1.0, omega=1.0, kB=1.0):
    """
    Exact quantum heat capacity of a single 1-D harmonic oscillator
    (the Einstein oscillator formula),

        Cv(T) = k_B * x^2 * exp(x) / (exp(x) - 1)^2,   x = hbar*omega / (k_B*T)

    Valid for scalar or array-like T (vectorized via NumPy).

    Parameters
    ----------
    T : float or array_like
        Temperature(s), in the same energy-equivalent units as hbar*omega/k_B.
    hbar : float, optional
        Reduced Planck constant (default 1.0).
    omega : float, optional
        Oscillator angular frequency (default 1.0).
    kB : float, optional
        Boltzmann constant (default 1.0).

    Returns
    -------
    Cv : float or ndarray
        Exact quantum heat capacity at each T, same shape as input T.
    """
    T = np.asarray(T, dtype=float)
    x = (hbar * omega) / (kB * T)
    return kB * (x**2 * np.exp(x)) / (np.exp(x) - 1.0)**2


# =====================================================================
# Exact HO classical (high-T) limit of Cv
# =====================================================================
def analytic_cv_HO_classical(kB=1.0):
    """
    Exact classical (high-temperature) limit of the 1-D harmonic
    oscillator heat capacity. By the equipartition theorem, each
    quadratic degree of freedom (kinetic + potential, here) 
    contributes (1/2) k_B, giving Cv_classical = k_B for a 1-D HO.

    This is a constant (temperature-independent) by definition of
    being the *classical* limit -- it is what the quantum Cv(T)
    curve approaches as T -> infinity.

    Parameters
    ----------
    kB : float, optional
        Boltzmann constant (default 1.0).

    Returns
    -------
    Cv_classical : float
        The classical limit value, equal to kB.
    """
    return kB
