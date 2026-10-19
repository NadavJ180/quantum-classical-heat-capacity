"""
Classical_Limit_Numerical.py
Heat capacity from a list of energy levels, and its classical limit.

The classical limit is hbar -> 0 at fixed V(x) and T. For every xi in a list
that you choose (config.py) the SAME potential is solved with a smaller Planck
constant, hbar/xi, and Cv(T) is computed from those levels. As xi grows, Cv(T)
approaches the classical heat capacity. This is the xi-scaling of
Gelbwaser-Klimovsky et al. (Eq. S7): solving xi^2 V at temperature xi^2 T gives
exactly the same Cv as solving V with hbar/xi at temperature T.
"""

import time

import numpy as np

from DVR.DVR_Algorithm import solve_energy_levels, number_of_grid_points


def compute_cv(energies, T):
    """Cv / k_B of a canonical ensemble with these levels at temperature T (k_B = 1):
    Cv = Var(E) / T^2, with Boltzmann probabilities p_n = exp(-E_n / T) / Z."""
    E = np.asarray(energies, dtype=float)
    w = np.exp(-(E - E.min()) / T)       # Boltzmann factors, measured from the lowest level
    p = w / w.sum()                       # occupation probability of each level
    E_mean = np.sum(p * E)
    variance = np.sum(p * (E - E_mean)**2)
    return variance / T**2


def heat_capacity_curve(energies, T_arr):
    """compute_cv at every temperature in T_arr."""
    return np.array([compute_cv(energies, T) for T in T_arr])


def classical_limit_curves(potential_func, T_arr, xi_values, num_levels, x_min, x_max, dx,
                           mass=1.0, hbar=1.0, walls=False):
    """Cv(T) for every xi in xi_values, from the levels of V(x) with hbar -> hbar/xi.
    Each xi uses dx/xi and num_levels*xi levels: the de Broglie wavelength and the
    level spacing both shrink by xi, so this keeps the same points per wavelength
    and the same energy window as xi = 1. Returns {xi: Cv(T) array}."""
    curves = {}
    for xi in xi_values:
        t0 = time.time()
        levels = int(round(num_levels * xi))
        points = number_of_grid_points(x_min, x_max, dx / xi, walls)
        print(f"    xi = {xi:<5g} hbar/xi = {hbar / xi:<7.4g} dx = {dx / xi:<8.4g} {levels} levels,"
              f" {points} points ({8e-9 * points**2:.2g} GB matrix) ...", end=" ", flush=True)
        energies = solve_energy_levels(potential_func, levels, x_min, x_max, dx / xi,
                                       mass, hbar / xi, walls)
        curves[xi] = heat_capacity_curve(energies, T_arr)
        print(f"{time.time() - t0:.1f} s", flush=True)
    return curves


def classical_error_estimate(curves, xi_values):
    """Estimated distance of the largest-xi curve from the xi -> infinity limit.
    The leading quantum correction to the classical Cv is of order hbar^2
    (Wigner-Kirkwood), so Cv(xi) = Cv_cl + a / xi^2 for large xi. From the last
    two xi values:  a / xi_2^2 = |Cv(xi_2) - Cv(xi_1)| / ((xi_2 / xi_1)^2 - 1)."""
    if len(xi_values) < 2:
        return np.full(len(next(iter(curves.values()))), np.nan)
    xi_1, xi_2 = xi_values[-2], xi_values[-1]
    return np.abs(curves[xi_2] - curves[xi_1]) / ((xi_2 / xi_1)**2 - 1.0)


def cv_vs_number_of_levels(energies, T):
    """Cv at temperature T when only the lowest n levels are kept, for n = 2 .. N.
    Where the curve is flat, the levels above n are not populated at this T."""
    E = np.asarray(energies, dtype=float)
    n_values = np.arange(2, len(E) + 1)
    cv_values = np.array([compute_cv(E[:n], T) for n in n_values])
    return n_values, cv_values
