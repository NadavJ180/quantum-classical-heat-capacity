"""
plot_potential.py  (src/figures/)
=====================================================================
Plots the potential V(x) for whatever system is currently configured
in config.py, with the ACTUAL computed energy levels overlaid -- not
a placeholder count. Every input (potential, mass, hbar, bounds,
number of levels) is read from config.py and DVR_Algorithm.py, the
same modules Quantum_HO_Master.py itself uses for Section 1. That
means this script always plots exactly the system currently being
researched: change the potential in config.py and re-run, nothing
here needs to be touched.

Because it calls the real grid auto-configuration and the real
3-pass converged DVR solve for all NUM_STATES levels, this script
does real (if modest) computation -- it is not instantaneous the way
a purely illustrative plot would be, but the eigenvalues it draws are
the same ones the rest of the pipeline is using.

Run from anywhere; paths are resolved relative to this file.

Usage:
    python plot_potential.py
Output:
    figures/<system>/<params>/energy_levels/potential.png
    (see src/figures/output_paths.py for the naming scheme)
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Make `src/` (for config.py) and `src/DVR/` (via the `DVR` package)
# importable regardless of the current working directory this script
# is launched from.
# ---------------------------------------------------------------------
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from config import MASS, HBAR, NUM_STATES, SYSTEM_NAME, POTENTIAL_PARAMS, my_potential  # noqa: E402
from DVR.DVR_Algorithm import (                                        # noqa: E402
    auto_configure_dvr,
    get_fully_converged_energy_levels,
)
from figures.output_paths import set_context, save_figure             # noqa: E402

# ============================== USER CONFIG ===============================
# Number of levels to draw. Defaults to ALL of NUM_STATES (i.e. exactly the
# scope currently being researched, per config.py) -- override to an int
# only if you want a deliberately reduced, less visually dense subset.
LEVELS_TO_DRAW = NUM_STATES
# ============================================================================


def  classical_turning_points(V, E, x):
    """Leftmost/rightmost x where V(x) <= E, for drawing a level segment
    that spans the classically allowed region at that energy."""
    allowed = x[V(x) <= E]
    if allowed.size == 0:
        return None
    return allowed.min(), allowed.max()


def main():
    set_context(SYSTEM_NAME, POTENTIAL_PARAMS)
    print(f"Configuring grid for '{SYSTEM_NAME}', {NUM_STATES} levels ...")
    x_min, x_max, n_grid = auto_configure_dvr(
        my_potential, NUM_STATES, mass=MASS, hbar=HBAR
    )

    print("Solving (3-pass converged) for the real energy levels "
          "-- this reuses the same call as Section 1, so it can take "
          "a little while for large NUM_STATES ...")
    energies = get_fully_converged_energy_levels(
        potential_func=my_potential,
        num_levels=NUM_STATES,
        x_min=x_min, x_max=x_max, num_points=n_grid,
        mass=MASS, hbar=HBAR,
    )

    levels = np.sort(np.asarray(energies))[:LEVELS_TO_DRAW]
    n_levels = len(levels)

    x = np.linspace(x_min, x_max, 4000)
    V = my_potential(x)

    # The auto-configured DVR grid is sized for numerical accuracy at
    # NUM_STATES levels, so V(x) at its far edges can be orders of
    # magnitude above anything physically interesting here. Left at
    # full scale, those steep wings dominate the y-axis and flatten
    # the actual well shape -- and any asymmetry between multiple
    # wells -- into an indistinguishable sliver at the bottom. Zoom to
    # a window that comfortably contains every drawn level instead,
    # since that's what this figure is actually meant to show.
    v_bottom = V.min()
    y_top = levels[-1] + 0.15 * (levels[-1] - v_bottom)
    y_bottom = v_bottom - 0.05 * (y_top - v_bottom)
    x_visible = x[V <= y_top]
    x_lo, x_hi = (x_visible.min(), x_visible.max()) if x_visible.size else (x_min, x_max)
    x_pad = 0.08 * (x_hi - x_lo)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, V, color="#1f4e8c", linewidth=2.2, label="$V(x)$", zorder=3)
    ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
    ax.set_ylim(y_bottom, y_top)

    # Many levels (NUM_STATES can be in the hundreds) -> thin, semi-
    # transparent, colormap-graded lines rather than individually labelled
    # ones. This renders as a density gradient that is still informative
    # (dense near the bottom, sparser as E grows) instead of a solid block.
    cmap = plt.cm.autumn_r
    label_every = max(1, n_levels // 10)  # label ~10 levels even if N is large
    for i, E in enumerate(levels):
        tp = classical_turning_points(my_potential, E, x)
        if tp is None:
            continue
        xl, xr = tp
        color = cmap(0.15 + 0.8 * i / max(n_levels - 1, 1))
        if i % 10 == 0:
            ax.hlines(E, xl, xr, color=color, linewidth=1.0,
                    alpha=0.55 if n_levels > 40 else 1.0, zorder=2)
        if n_levels <= 40 or i % label_every == 0 or i == n_levels - 1:
            ax.text(xr + 0.015 * (x_hi - x_lo), E, f"$n={i}$",
                     va="center", fontsize=7, color="#555555")

    ax.set_xlabel("$x$")
    ax.set_ylabel("Energy")
    ax.set_title(f"{SYSTEM_NAME} -- potential and computed spectrum\n"
                 f"({n_levels} levels shown)", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(loc="upper center")
    fig.tight_layout()

    out_path = save_figure(fig, "energy_levels", "potential")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()