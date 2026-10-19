"""
DVR_Algorithm.py
Energy levels of a particle of mass m in a 1-D potential V(x), from the
Discrete Variable Representation (DVR) of Colbert & Miller (1992).

    colbert_miller_dvr_1d : V(x) itself confines the particle (grid with open ends)
    box_dvr_1d            : infinite walls at x_min and x_max (particle in a box + V)
    solve_energy_levels   : picks one of the two and builds the grid from dx
    suggest_grid          : physics estimate of the grid a temperature needs (printed only)
"""

import numpy as np
import scipy.linalg as la


def colbert_miller_dvr_1d(potential_func, num_levels, x_min, x_max, num_points, mass=1.0, hbar=1.0):
    """Lowest `num_levels` eigenvalues of H = p^2/2m + V(x) on `num_points` evenly
    spaced points from x_min to x_max (sinc-DVR, Colbert & Miller eq. A7).
    The grid ends are not walls: V must already be large there."""
    if num_points <= num_levels:
        raise ValueError(f"{num_levels} levels need more than {num_points} grid points.")
    x = np.linspace(x_min, x_max, num_points)
    dx = x[1] - x[0]

    # Kinetic energy between grid points i and j depends only on n = |i - j|:
    #   T(0) = hbar^2/(2 m dx^2) * pi^2/3,   T(n) = hbar^2/(2 m dx^2) * 2 (-1)^n / n^2
    n = np.arange(1, num_points)
    kinetic_row = np.zeros(num_points)
    kinetic_row[0] = np.pi**2 / 3.0
    kinetic_row[1:] = 2.0 * (-1.0)**n / n**2
    kinetic_row *= hbar**2 / (2.0 * mass * dx**2)

    H = la.toeplitz(kinetic_row)                    # kinetic energy matrix
    V = np.asarray(potential_func(x), dtype=float)
    if not np.all(np.isfinite(V)):
        raise ValueError("V(x) is infinite or NaN on the grid; for hard walls use WALLS = True.")
    i = np.arange(num_points)
    H[i, i] += V                                    # potential energy on the diagonal
    return la.eigvalsh(H, overwrite_a=True, subset_by_index=[0, num_levels - 1])


def box_dvr_1d(potential_func, num_levels, x_min, x_max, num_points, mass=1.0, hbar=1.0):
    """Same as colbert_miller_dvr_1d, but with infinite walls at x_min and x_max
    (Colbert & Miller eq. A6). The `num_points` grid points are inside the box;
    the wavefunction is zero on the walls."""
    if num_points <= num_levels:
        raise ValueError(f"{num_levels} levels need more than {num_points} grid points.")
    N = num_points + 1                              # number of intervals between the walls
    L = x_max - x_min
    k = np.arange(1, N)                             # grid point labels 1 .. N-1
    x = x_min + L * k / N

    H = np.zeros((num_points, num_points))
    for row in range(num_points):
        i = k[row]
        with np.errstate(divide="ignore", invalid="ignore"):
            H[row] = (-1.0)**(i - k) * (1.0 / np.sin(np.pi * (i - k) / (2 * N))**2
                                        - 1.0 / np.sin(np.pi * (i + k) / (2 * N))**2)
        H[row, row] = (2.0 * N**2 + 1.0) / 3.0 - 1.0 / np.sin(np.pi * i / N)**2
    H *= hbar**2 / (2.0 * mass) * np.pi**2 / (2.0 * L**2)

    V = np.asarray(potential_func(x), dtype=float)
    if not np.all(np.isfinite(V)):
        raise ValueError("V(x) is infinite or NaN inside the box.")
    i = np.arange(num_points)
    H[i, i] += V
    return la.eigvalsh(H, overwrite_a=True, subset_by_index=[0, num_levels - 1])


def number_of_grid_points(x_min, x_max, dx, walls=False):
    """Grid points for spacing dx: both ends included, or (walls) only the inside of the box."""
    if walls:
        return int(round((x_max - x_min) / dx)) - 1
    return int(round((x_max - x_min) / dx)) + 1


def solve_energy_levels(potential_func, num_levels, x_min, x_max, dx, mass=1.0, hbar=1.0, walls=False):
    """The lowest `num_levels` energies on a grid of spacing dx between x_min and x_max.
    walls=True puts infinite walls at x_min and x_max."""
    num_points = number_of_grid_points(x_min, x_max, dx, walls)
    if walls:
        return box_dvr_1d(potential_func, num_levels, x_min, x_max, num_points, mass, hbar)
    return colbert_miller_dvr_1d(potential_func, num_levels, x_min, x_max, num_points, mass, hbar)


def suggest_grid(potential_func, T_max, x_min, x_max, mass=1.0, hbar=1.0, thermal_cutoff=20.0):
    """What the grid needs so that k_B T_max is described correctly. Levels more than
    thermal_cutoff * k_B T above the bottom of V carry a Boltzmann weight below
    exp(-thermal_cutoff), so the highest energy that matters is E_top = V_min + cutoff * T_max.
      - region: where V(x) < E_top (the classically allowed region at E_top)
      - num_levels: Bohr-Sommerfeld count of levels below E_top, (1/pi hbar) * integral of p dx
      - dx: 4 grid points per shortest de Broglie wavelength, 2 pi hbar / p_max
    Only looks at V inside [x_min, x_max]."""
    x = np.linspace(x_min, x_max, 20001)
    V = np.asarray(potential_func(x), dtype=float)
    v_min = V.min()
    E_top = v_min + thermal_cutoff * T_max
    allowed = x[V <= E_top]
    p = np.sqrt(2.0 * mass * np.clip(E_top - V, 0.0, None))   # classical momentum at E_top
    return {
        "E_top": E_top,
        "v_min": v_min,
        "region": (allowed.min(), allowed.max()),
        "touches_edge": bool(V[0] <= E_top or V[-1] <= E_top),
        "num_levels": np.sum(p) * (x[1] - x[0]) / (np.pi * hbar),
        "dx": np.pi * hbar / (2.0 * np.sqrt(2.0 * mass * (E_top - v_min))),
    }
