"""
plot_potential.py  (src/figures/)
V(x) with the computed energy levels drawn on it, each level spanning its
classically allowed region (where V(x) < E). Two figures: all levels, and a
zoom on the bottom of the well(s).

    python figures/plot_potential.py     (uses config.py and config.QUICK_SCAN)
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

SRC_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from figures.output_paths import save_figure                          # noqa: E402


def classical_turning_points(V, E, x):
    """Leftmost and rightmost x where V(x) <= E (None if there is none)."""
    allowed = x[V(x) <= E]
    if allowed.size == 0:
        return None
    return allowed.min(), allowed.max()


def _draw(x, V, potential_func, levels, system_name, y_top, title, filename):
    v_bottom = V.min()
    y_bottom = v_bottom - 0.05 * (y_top - v_bottom)
    visible = x[V <= y_top]
    x_lo, x_hi = (visible.min(), visible.max()) if visible.size else (x[0], x[-1])
    pad = 0.08 * (x_hi - x_lo)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, V, color="#1f4e8c", lw=2.2, label="$V(x)$", zorder=3)
    in_view = [(n, E) for n, E in enumerate(levels) if E <= y_top]
    step = max(1, len(in_view) // 40)                 # draw at most ~40 levels
    cmap = plt.get_cmap("autumn_r")
    for k, (n, E) in enumerate(in_view):
        if k % step != 0:
            continue
        tp = classical_turning_points(potential_func, E, x)
        if tp is None:
            continue
        ax.hlines(E, tp[0], tp[1], color=cmap(0.15 + 0.8 * n / max(len(levels) - 1, 1)), lw=1.0, zorder=2)
        if len(in_view) <= 15 or k % (10 * step) == 0:
            ax.text(tp[1] + 0.015 * (x_hi - x_lo), E, f"$n={n}$", va="center", fontsize=7, color="#555555")
    ax.set_xlim(x_lo - pad, x_hi + pad)
    ax.set_ylim(y_bottom, y_top)
    ax.set_xlabel("$x$")
    ax.set_ylabel("Energy")
    ax.set_title(f"{system_name}: {title}\n({len(in_view)} of {len(levels)} levels in view)",
                 fontsize=11, fontweight="bold")
    ax.grid(alpha=0.3, ls="--")
    ax.legend(loc="upper center")
    fig.tight_layout()
    return save_figure(fig, "energy_levels", filename)


def plot_potential_with_spectrum(x_min, x_max, potential_func, energies, system_name):
    """Save both figures for levels already computed on [x_min, x_max]."""
    levels = np.sort(np.asarray(energies))
    x = np.linspace(x_min, x_max, 4000)
    V = potential_func(x)
    v_bottom = V.min()
    full_top = levels[-1] + 0.15 * (levels[-1] - v_bottom)
    # Zoom: just above the highest barrier between wells (a local maximum of V),
    # so the shape of the wells is visible; with no barrier, up to the 15th level.
    # At least the lowest 3 levels are always in view.
    barriers = [V[i] for i in range(1, len(V) - 1) if V[i - 1] < V[i] > V[i + 1]]
    if barriers:
        zoom_top = max(barriers) + 0.3 * (max(barriers) - v_bottom)
    else:
        zoom_top = levels[min(len(levels), 15) - 1]
    zoom_top = max(zoom_top, levels[min(len(levels), 3) - 1])
    zoom_top += 0.15 * (zoom_top - v_bottom)
    return {
        "full_spectrum": _draw(x, V, potential_func, levels, system_name, full_top,
                               "potential and computed levels", "potential_full_spectrum"),
        "zoomed": _draw(x, V, potential_func, levels, system_name, zoom_top,
                        "lowest levels", "potential_zoomed"),
    }


def main():
    import config
    from DVR.DVR_Algorithm import solve_energy_levels
    from figures.output_paths import set_context

    set_context(config.SYSTEM_NAME, config.POTENTIAL_PARAMS)
    s = config.QUICK_SCAN
    energies = solve_energy_levels(config.my_potential, s["num_states"], config.X_MIN, config.X_MAX,
                                   s["dx"], config.MASS, config.HBAR, config.WALLS)
    paths = plot_potential_with_spectrum(config.X_MIN, config.X_MAX, config.my_potential, energies,
                                         config.SYSTEM_NAME)
    print(f"Saved {paths['full_spectrum']}\nSaved {paths['zoomed']}")


if __name__ == "__main__":
    main()
