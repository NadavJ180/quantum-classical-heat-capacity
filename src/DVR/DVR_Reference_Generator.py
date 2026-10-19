"""
DVR_Reference_Generator.py
A reference spectrum: the same DVR on a wider and finer grid. It is the
"ground truth" for potentials with no exact solution. If the base grid and the
reference agree, the base grid is wide enough (the wavefunctions die out before
the grid ends) and fine enough (enough points per de Broglie wavelength).
"""

import time

from DVR.DVR_Algorithm import solve_energy_levels


def reference_grid(x_min, x_max, dx, span_factor, dx_factor, walls=False):
    """Grid that is span_factor times wider (same centre) and dx_factor times finer.
    With walls the box size is part of the physics, so only dx changes."""
    if walls:
        span_factor = 1.0
    centre = 0.5 * (x_min + x_max)
    half_span = 0.5 * (x_max - x_min) * span_factor
    return centre - half_span, centre + half_span, dx / dx_factor


def generate_reference_energies(potential_func, num_levels, x_min, x_max, dx,
                                span_factor=2.0, dx_factor=2.0, mass=1.0, hbar=1.0, walls=False):
    """The lowest num_levels energies on the reference grid."""
    x_min_ref, x_max_ref, dx_ref = reference_grid(x_min, x_max, dx, span_factor, dx_factor, walls)
    print(f"  Base grid:      [{x_min:.4g}, {x_max:.4g}], dx = {dx:.4g}")
    print(f"  Reference grid: [{x_min_ref:.4g}, {x_max_ref:.4g}], dx = {dx_ref:.4g}"
          f"  ({(x_max_ref - x_min_ref) / dx_ref:.0f} points)", flush=True)
    t0 = time.time()
    energies = solve_energy_levels(potential_func, num_levels, x_min_ref, x_max_ref, dx_ref,
                                   mass, hbar, walls)
    print(f"  Reference solved in {time.time() - t0:.1f} s")
    return {"energies": energies, "x_min": x_min_ref, "x_max": x_max_ref, "dx": dx_ref}
