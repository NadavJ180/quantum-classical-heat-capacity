"""
DVR_Limit_Finder.py
Where does the DVR stop being trustworthy? The grid is redone with a list of
grid spacings (and a list of grid widths), and every level is compared with the
reference spectrum. A level is "correct" if it is within `tolerance` of the
reference. This shows, for each dx (or width), how many levels can be trusted.

Physics: a level needs a few grid points per de Broglie wavelength, 2 pi hbar / p,
and higher levels have shorter wavelengths, so a coarser dx loses the top levels
first. A narrower grid cuts into the wavefunction tails, which are longest for
the highest levels, so it also loses the top levels first.
"""

import numpy as np
import matplotlib.pyplot as plt

from DVR.DVR_Algorithm import solve_energy_levels, number_of_grid_points
from DVR.DVR_Reference_Generator import reference_grid
from figures.output_paths import save_figure


def count_correct_levels(errors, tolerance):
    """How many of the lowest levels are all within tolerance of the reference."""
    bad = np.nonzero(errors > tolerance)[0]
    return len(errors) if len(bad) == 0 else int(bad[0])


def errors_vs_dx(potential_func, num_levels, x_min, x_max, dx_values, reference_energies,
                 mass=1.0, hbar=1.0, walls=False):
    """{dx: |E_n - E_n,ref| for n = 0 .. num_levels-1}, same grid width, each dx in dx_values."""
    errors = {}
    for dx in dx_values:
        if number_of_grid_points(x_min, x_max, dx, walls) <= num_levels:
            print(f"    dx = {dx:.4g}: too few grid points for {num_levels} levels, skipped")
            continue
        E = solve_energy_levels(potential_func, num_levels, x_min, x_max, dx, mass, hbar, walls)
        errors[dx] = np.abs(E - reference_energies[:num_levels])
    return errors


def errors_vs_span(potential_func, num_levels, x_min, x_max, dx, span_factors, reference_energies,
                   mass=1.0, hbar=1.0):
    """{span factor: |E_n - E_n,ref|}, same dx, grid width multiplied by each factor."""
    errors = {}
    for factor in span_factors:
        x_lo, x_hi, _ = reference_grid(x_min, x_max, dx, factor, 1.0)
        if number_of_grid_points(x_lo, x_hi, dx) <= num_levels:
            print(f"    span x{factor:g}: too few grid points for {num_levels} levels, skipped")
            continue
        E = solve_energy_levels(potential_func, num_levels, x_lo, x_hi, dx, mass, hbar)
        errors[factor] = np.abs(E - reference_energies[:num_levels])
    return errors


def plot_errors(errors, label, tolerance, title, name):
    """|E_n - E_n,ref| against n, one curve per tested value."""
    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = plt.get_cmap("viridis")
    for k, (value, err) in enumerate(errors.items()):
        n_good = count_correct_levels(err, tolerance)
        ax.plot(np.arange(len(err)), np.where(err > 0, err, np.nan), lw=1.3,
                color=colors(k / max(len(errors) - 1, 1)),
                label=f"{label} = {value:.4g}  ({n_good} correct levels)")
    ax.axhline(tolerance, color="r", ls=":", label=f"tolerance = {tolerance:g}")
    ax.set_yscale("log")
    ax.set_xlabel("level index n")
    ax.set_ylabel(r"$|E_n - E_{n,\mathrm{ref}}|$")
    ax.set_title(title, fontsize=10)
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend(fontsize=8)
    plt.tight_layout()
    return save_figure(fig, "dvr_limits", name)


def run_dvr_limit_analysis(potential_func, system_name, reference_energies, num_levels, x_min, x_max,
                           dx, dx_factors, span_factors, tolerance, mass=1.0, hbar=1.0, walls=False):
    """Both scans, their figures, and a console table."""
    dx_values = [dx * f for f in dx_factors]
    print(f"  Grid spacing scan (width fixed at [{x_min:g}, {x_max:g}]):")
    by_dx = errors_vs_dx(potential_func, num_levels, x_min, x_max, dx_values, reference_energies,
                         mass, hbar, walls)
    for value, err in by_dx.items():
        print(f"    dx = {value:<8.4g}: {count_correct_levels(err, tolerance):>5} correct levels"
              f"   (max error {np.max(err):.1e})")
    plot_errors(by_dx, "dx", tolerance, f"{system_name}: level errors vs grid spacing", "level_errors_vs_dx")

    if walls:
        print("  Grid width scan skipped: with WALLS = True the box size is part of the physics.")
        return {"by_dx": by_dx, "by_span": {}}
    print(f"  Grid width scan (dx fixed at {dx:g}):")
    by_span = errors_vs_span(potential_func, num_levels, x_min, x_max, dx, span_factors,
                             reference_energies, mass, hbar)
    for value, err in by_span.items():
        print(f"    span x{value:<6.3g}: {count_correct_levels(err, tolerance):>5} correct levels"
              f"   (max error {np.max(err):.1e})")
    plot_errors(by_span, "span factor", tolerance, f"{system_name}: level errors vs grid width",
                "level_errors_vs_span")
    return {"by_dx": by_dx, "by_span": by_span}
