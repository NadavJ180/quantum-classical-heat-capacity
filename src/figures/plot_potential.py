"""
plot_potential.py  (src/figures/)
=====================================================================
Plots the potential V(x) for whatever system is currently configured
in config.py, with the ACTUAL computed energy levels overlaid -- not
a placeholder count. Produces TWO complementary figures:

    potential_full_spectrum.png -- zoomed just enough to show every
        drawn level (up to the highest one), so you can see how many
        levels were actually needed and how the spacing changes going
        up the spectrum.
    potential_zoomed.png -- zoomed tightly on the well's own local
        structure near its minima (any interior barrier, for a multi-
        well potential), so the shape itself -- and any asymmetry
        between wells -- is actually visible, even though that usually
        means most of a large requested spectrum falls outside the
        frame.

Both use the same density-graded level-drawing heuristic (dense,
individually-labelled lines when there are few levels in view; a
thinned, colormap-graded subset when there are many), just scoped to
however many of the requested levels fall inside each figure's own
window.

The reusable part is `plot_potential_with_spectrum`: given a grid span,
potential function, and an already-computed spectrum, it builds and
saves both figures without doing any DVR work of its own. This is what
lets Quantum_HO_Master.py call it directly from its own Section 1,
reusing the base grid and energies it already computed there instead
of paying for a second 3-pass converged solve -- so running the master
pipeline once is enough; this script never needs to be run separately
just to get these figures.

Run as a standalone script (`python plot_potential.py`), `main()` does
its own grid auto-configuration and 3-pass converged solve first (for
when you want just these figures, without running the full pipeline).

Run from anywhere; paths are resolved relative to this file.

Usage:
    python plot_potential.py
Output:
    figures/<system>/<params>/energy_levels/potential_full_spectrum.png
    figures/<system>/<params>/energy_levels/potential_zoomed.png
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

from figures.output_paths import save_figure                          # noqa: E402


def  classical_turning_points(V, E, x):
    """Leftmost/rightmost x where V(x) <= E, for drawing a level segment
    that spans the classically allowed region at that energy."""
    allowed = x[V(x) <= E]
    if allowed.size == 0:
        return None
    return allowed.min(), allowed.max()


def _full_spectrum_window(x, V, v_bottom, levels):
    """
    Zoom just enough to show every one of the requested/drawn levels
    (up to the highest one) -- the "overview" companion to
    `_origin_zoom_window`. Still zoomed relative to the full
    auto-configured DVR grid (which is sized for numerical accuracy,
    not for being a good picture, and can be orders of magnitude
    taller than any requested level), just not zoomed as tightly as
    the origin view.
    """
    if len(levels):
        y_top = levels[-1] + 0.15 * (levels[-1] - v_bottom)
    else:
        y_top = v_bottom + 1.0
    y_bottom = v_bottom - 0.05 * (y_top - v_bottom)
    x_visible = x[V <= y_top]
    x_lo, x_hi = (x_visible.min(), x_visible.max()) if x_visible.size else (x[0], x[-1])
    return x_lo, x_hi, y_bottom, y_top


def _origin_zoom_window(x, V, v_bottom, levels):
    """
    Pick a (x_lo, x_hi, y_bottom, y_top) window that actually shows the
    well's shape near its own minima, rather than the full requested
    spectrum (at the far edges of the auto-configured DVR grid, V(x)
    can reach orders of magnitude above anything physically
    interesting, and plotted at that scale it swamps the y-axis and
    flattens the actual well shape, and any asymmetry between multiple
    wells, into an indistinguishable sliver at the bottom).

    Method: find every local extremum of the sampled V(x) by sign
    changes in its discrete derivative, and keep only the maxima (a
    local max between two wells is a barrier). If at least one interior
    barrier exists, zoom just above the highest one -- this reveals the
    multi-well structure near the origin even if that means most of a
    large requested spectrum falls outside the frame, which is exactly
    the point (their shape, not the full level count, is what this
    figure is showing). If there's no interior barrier at all (e.g. a
    single well like the HO), fall back to a window sized around a
    modest number of the lowest levels instead.
    """
    dV = np.diff(V)
    sign_changes = np.diff(np.sign(dV))
    maxima_idx = np.where(sign_changes < 0)[0] + 1  # + -> - : local max
    barriers = V[maxima_idx]

    if barriers.size:
        barrier_top = barriers.max()
        y_top = barrier_top + 0.3 * (barrier_top - v_bottom)
    else:
        cap = min(len(levels), 15)
        y_top = levels[cap - 1] if cap else v_bottom + 1.0

    y_bottom = v_bottom - 0.05 * (y_top - v_bottom)
    x_visible = x[V <= y_top]
    x_lo, x_hi = (x_visible.min(), x_visible.max()) if x_visible.size else (x[0], x[-1])
    return x_lo, x_hi, y_bottom, y_top


def _render_and_save(x, V, potential_func, levels, n_levels, system_name, window, title_suffix, filename):
    """
    Build one potential-plus-spectrum figure for a given (x_lo, x_hi,
    y_bottom, y_top) window and save it. Shared by both the full-
    spectrum and origin-zoomed figures -- they differ only in which
    window they're rendered into and their filename/title.

    Levels above the window's y_top are dropped entirely -- not just
    axes-clipped -- because an unclipped ax.text positioned far outside
    the visible range still inflates savefig(bbox_inches="tight") to
    cover it, which previously produced a runaway-tall PNG. Levels are
    also individually labelled while there are few of them in view, or
    thinned to a colormap-graded subset once there are many -- scoped
    to how many levels this window actually contains, not the full
    requested spectrum.
    """
    x_lo, x_hi, y_bottom, y_top = window
    x_pad = 0.08 * (x_hi - x_lo)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, V, color="#1f4e8c", linewidth=2.2, label="$V(x)$", zorder=3)
    ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
    ax.set_ylim(y_bottom, y_top)

    in_view = [(i, E) for i, E in enumerate(levels) if E <= y_top]
    n_candidates = len(in_view)
    thin = n_candidates > 40
    label_every = max(1, n_candidates // 10)
    cmap = plt.cm.autumn_r
    n_shown = 0
    for j, (i, E) in enumerate(in_view):
        tp = classical_turning_points(potential_func, E, x)
        if tp is None:
            continue
        xl, xr = tp
        n_shown += 1
        color = cmap(0.15 + 0.8 * i / max(n_levels - 1, 1))
        if not thin or j % 10 == 0:
            ax.hlines(E, xl, xr, color=color, linewidth=1.0,
                    alpha=0.55 if thin else 1.0, zorder=2)
        if not thin or j % label_every == 0 or j == n_candidates - 1:
            ax.text(xr + 0.015 * (x_hi - x_lo), E, f"$n={i}$",
                     va="center", fontsize=7, color="#555555")

    ax.set_xlabel("$x$")
    ax.set_ylabel("Energy")
    shown = f"{n_shown} of {n_levels}" if n_shown < n_levels else f"{n_levels}"
    ax.set_title(f"{system_name} -- {title_suffix}\n"
                 f"({shown} levels shown)", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(loc="upper center")
    fig.tight_layout()

    return save_figure(fig, "energy_levels", filename)


def plot_potential_with_spectrum(x_min, x_max, potential_func, energies, system_name,
                                  levels_to_draw=None):
    """
    Build and save both potential-plus-spectrum figures (full-spectrum
    overview and origin-zoomed) for an already solved system -- does no
    DVR work itself, so it's cheap to call right after a driver script
    has already computed `energies` on [x_min, x_max].

    Parameters
    ----------
    x_min, x_max : float
        The grid span the energies were computed on (only used here to
        sample V(x) for the blue curve; the actual displayed windows
        are zoomed in separately for each figure).
    potential_func : callable
        V(x) -> float or ndarray, e.g. config.my_potential.
    energies : array_like
        Computed energy levels (any order; sorted internally).
    system_name : str
        Used in the plot titles. The figures' save location
        (figures/<system>/<params>/...) is controlled separately by
        whatever already called figures.output_paths.set_context --
        this function does not call it, so the caller's context (set
        once by the driver script) is what's used.
    levels_to_draw : int or None, optional
        How many of the lowest computed levels to overlay (default:
        all of `energies`).

    Returns
    -------
    paths : dict
        {"full_spectrum": path, "zoomed": path} -- see
        `figures.output_paths.save_figure`.
    """
    levels = np.sort(np.asarray(energies))[:levels_to_draw]
    n_levels = len(levels)

    x = np.linspace(x_min, x_max, 4000)
    V = potential_func(x)
    v_bottom = V.min()

    full_path = _render_and_save(
        x, V, potential_func, levels, n_levels, system_name,
        window=_full_spectrum_window(x, V, v_bottom, levels),
        title_suffix="potential and computed spectrum (full)",
        filename="potential_full_spectrum",
    )
    zoomed_path = _render_and_save(
        x, V, potential_func, levels, n_levels, system_name,
        window=_origin_zoom_window(x, V, v_bottom, levels),
        title_suffix="potential and computed spectrum (zoomed on minima)",
        filename="potential_zoomed",
    )
    return {"full_spectrum": full_path, "zoomed": zoomed_path}


def main():
    from config import MASS, HBAR, NUM_STATES, SYSTEM_NAME, POTENTIAL_PARAMS, my_potential
    from DVR.DVR_Algorithm import auto_configure_dvr, get_fully_converged_energy_levels
    from figures.output_paths import set_context

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

    out_paths = plot_potential_with_spectrum(
        x_min, x_max, my_potential, energies, SYSTEM_NAME, levels_to_draw=NUM_STATES
    )
    print(f"Saved {out_paths['full_spectrum']}")
    print(f"Saved {out_paths['zoomed']}")


if __name__ == "__main__":
    main()
