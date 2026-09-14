"""
Cv_Coefficient_Sweep.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Section 7 of the pipeline: compares several variants of the run's
potential that differ only in one named POTENTIAL_PARAMS coefficient
(e.g. the double well's cubic term `b`), to see how that coefficient
affects the shape/size of a (numerical) Schottky anomaly. Produces,
all under one dedicated `coefficient_sweep/` figure folder (see
figures/output_paths.py -- deliberately NOT the `cv`/`energy_levels`
folders the rest of the pipeline uses, since these figures are about
comparing variants, not diagnosing the base run):

    cv_coefficient_sweep.png    -- quantum Cv(T) per variant, overlaid
                                    with the base run's classical limit
    potential_comparison.png    -- V(x) + low-lying spectrum per variant,
                                    for correlating a Cv anomaly with the
                                    potential shape that causes it (side-
                                    by-side panels, falling back to one
                                    overlay plot, then to one figure per
                                    variant, as the variant count grows --
                                    see `plot_variant_potentials`)

WHY NO XI/N-CONVERGENCE SEARCH HERE
---------------------------------------------------------------------
The classical-limit curve is reused verbatim from the base run's
`base_cv_results["cv_classical"]` -- it is NOT recomputed per variant.
Each variant's quantum Cv(T) is instead a direct evaluation of the
literal, unscaled quantum heat capacity (xi=1.0, no rescaling): given
a spectrum, Cv(T) = k_B*beta^2*Var(E) is an exact formula, not an
approximation that needs a convergence search -- the xi/n-convergence
machinery elsewhere in this project exists only to locate the
classical-limit PLATEAU (a genuinely different, harder problem), which
this file deliberately does not attempt per variant.

IS IT VALID TO REUSE THE BASE RUN'S NUM_STATES FOR EVERY VARIANT?
---------------------------------------------------------------------
Mostly, but not unconditionally -- and this file checks for the one
case where it can fail rather than assuming it's always fine.

A direct quantum Cv(T) evaluation from a TRUNCATED spectrum (only the
lowest `num_states` levels) is exact for every temperature at which
the omitted, higher levels would have carried negligible Boltzmann
weight anyway. Concretely (the same criterion `Cv_AutoTune.py` already
uses for the base run's own escalation):

    E_max - E_0  >>  k_B * T_hot      (T_hot = 1 / beta_arr.min())

where E_max is the top level actually computed. `NUM_STATES` was
chosen (via the base run's auto-tune loop) to satisfy this for the
BASE potential's own spectrum. Reusing that same level COUNT for a
different coefficient value is safe only if that variant's spectrum
reaches at least as high an E_max for the same level count -- which
is NOT guaranteed in general, because the level-spacing scaling
itself depends on the coefficient being swept. For example, WKB
quantization of a quartic well V ~ a*x^4 gives E_n ~ a^(1/3) * n^(4/3):
sweeping the leading coefficient `a` upward only ever makes a fixed
NUM_STATES MORE conservative (E_max grows), but sweeping it DOWNWARD
(toward zero, or negative -- exactly the direction that produced the
non-confining potential this file's resilience handling was built
for) shrinks E_max for the same level count, and could in principle
silently reproduce the very truncation artifact -- a numerical
Schottky-like collapse at the hot end -- that Cv_AutoTune exists to
prevent for the base run, without anything flagging it here.

`solve_variant_with_hot_coverage` closes that gap: after solving a
variant's own spectrum, it checks the SAME E_max/k_B/T_hot criterion
(`hot_state_safety`, reused from config.HOT_STATE_SAFETY) and, if it
fails, escalates that ONE variant's own NUM_STATES (reusing config's
NUM_STATES_GROWTH/NUM_STATES_CAP/MAX_ESCALATION_ROUNDS) and re-solves,
exactly mirroring the base run's own escalation loop but scoped to a
single variant. A variant that still can't clear the bar after that
is kept (its low/mid-T behavior -- where the anomaly actually shows --
is normally unaffected by hot-end truncation) but flagged, both in the
console output and with a marker in the plot legend, rather than
silently trusted.

WHAT HAPPENS WHEN A VARIANT FAILS TO SOLVE AT ALL
---------------------------------------------------------------------
Not every coefficient value produces a genuinely confining potential
(e.g. a quartic leading coefficient that goes negative is unbounded
from below) -- the DVR grid auto-configurator or its 3-pass
convergence check will then fail outright for that one variant.
`run_coefficient_sweep` catches that per variant, prints a diagnostic
(including a cheap, potential-agnostic unboundedness check -- see
`_diagnose_variant_failure`), and keeps going: the final plots are
built from whichever variants DID converge, with a printed summary of
which ones didn't and why, rather than losing the whole sweep to one
bad coefficient value.
=====================================================================
"""

import functools
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from DVR.DVR_Algorithm import auto_configure_dvr, get_fully_converged_energy_levels
from Quantum_Classical_Combined import compute_quantum_heat_capacity_curve
from figures.output_paths import save_figure
# Reused rather than reimplemented so the potential-comparison panels
# below use exactly the same "zoom on the well's own structure" window
# logic and turning-point drawing as the main pipeline's own potential
# figures (figures/plot_potential.py, Section 1 & 4) -- these two are
# genuinely private helpers of that module, imported explicitly by
# name (not via `import *`), so keep this in sync if their signatures
# ever change there.
from figures.plot_potential import classical_turning_points, _origin_zoom_window

# Variant-count thresholds for `plot_variant_potentials`'s side-by-side
# -> overlay -> separate-figures cascade (see its docstring).
SIDE_BY_SIDE_MAX = 6
OVERLAY_MAX = 15


# =====================================================================
# Build the list of per-variant parameter dicts
# =====================================================================
def generate_variant_params(base_params, scan_param, scan_step, scan_count):
    """
    Build `2*scan_count + 1` copies of `base_params`, each with
    `scan_param` shifted by a multiple of `scan_step`, centered on
    (and including exactly) the value already in `base_params`.

    Parameters
    ----------
    base_params : dict
        e.g. config.POTENTIAL_PARAMS.
    scan_param : str
        Key of `base_params` to vary.
    scan_step : float
        Spacing between consecutive variants.
    scan_count : int
        Number of EXTRA variants added on each side of the base value
        (total variants returned = 2*scan_count + 1).

    Returns
    -------
    list of dict
        One shallow copy of `base_params` per variant, ordered by
        ascending offset (most negative first); the middle entry
        (offset 0) is exactly `base_params`'s own value.

    Raises
    ------
    KeyError
        If `scan_param` is not a key of `base_params`.
    """
    if scan_param not in base_params:
        raise KeyError(
            f"SCAN_PARAM '{scan_param}' is not a key of POTENTIAL_PARAMS "
            f"(available keys: {sorted(base_params.keys())})."
        )
    base_value = base_params[scan_param]
    variants = []
    for offset in range(-scan_count, scan_count + 1):
        params = dict(base_params)
        params[scan_param] = base_value + offset * scan_step
        variants.append(params)
    return variants


# =====================================================================
# Build the plot annotation: the formula, fixed coefficients numeric,
# the swept one symbolic
# =====================================================================
def format_potential_formula(formula_template, base_params, scan_param):
    """
    Substitute `base_params`' values into `formula_template` for
    display on the coefficient-sweep plot, EXCEPT `scan_param`, which
    is left as its own symbol (e.g. "b") since that is the one thing
    that varies across the plotted curves.

    Placeholders are `<<name>>` tokens (not Python's `{name}`) so they
    never collide with LaTeX's own braces (e.g. `\\frac{1}{2}`) that
    may appear elsewhere in the template -- see config.POTENTIAL_FORMULA.
    Each substituted numeric value carries its own explicit sign
    (`+0.25`, `-0.5`, ...) and the template is expected to place NO
    literal operator between consecutive placeholders, so terms can be
    reordered or added without needing to fix up stray `+`/`-` signs;
    a redundant leading "+" right after "=" (when the first term in the
    template happens to be positive, or is the symbolic one) is
    stripped for readability.

    Parameters
    ----------
    formula_template : str or None
        e.g. config.POTENTIAL_FORMULA. A template with no `<<...>>`
        tokens at all is returned unchanged (useful for potentials
        whose formula doesn't decompose into one additive term per
        coefficient, e.g. the harmonic oscillator). None returns None.
    base_params : dict
        e.g. config.POTENTIAL_PARAMS -- the FIXED values (every key
        except `scan_param`) to substitute numerically.
    scan_param : str
        The coefficient shown symbolically instead of numerically.

    Returns
    -------
    str or None
    """
    if formula_template is None:
        return None
    text = formula_template
    for key, value in base_params.items():
        token = f"<<{key}>>"
        if token not in text:
            continue
        replacement = f" + {key}" if key == scan_param else f"{value:+.3g}"
        text = text.replace(token, replacement)
    # Collapse a redundant leading "+" (e.g. "= +0.25..." -> "= 0.25...",
    # or "=  + b..." -> "= b...") -- cosmetic only, never changes meaning.
    text = re.sub(r"=\s*\+\s*", "= ", text, count=1)
    return text


# =====================================================================
# Cheap, potential-agnostic diagnostic for a variant that failed to solve
# =====================================================================
def _diagnose_variant_failure(potential_func, probe=1000.0):
    """
    Best-effort explanation for why a variant's DVR solve might have
    failed to converge, without any knowledge of the specific
    potential's functional form. Samples V(x) far from the origin and
    compares it to V(x) near the origin: a genuinely confining
    potential must rise well above its near-origin values out there,
    so if it doesn't (or is lower), the potential is likely unbounded
    from below (or otherwise non-confining) for this coefficient
    value -- the classic cause of `auto_configure_dvr`'s turning-point
    search or span-expansion loop never settling.

    Parameters
    ----------
    potential_func : callable
        V(x) -> float or ndarray, e.g. the failing variant's potential.
    probe : float, optional
        How far from the origin to sample (default 1000.0).

    Returns
    -------
    str or None
        A one-line diagnostic message, or None if the cheap check
        didn't find anything to report (the failure has some other
        cause -- e.g. a genuinely pathological grid-resolution issue).
    """
    try:
        xs = np.array([-probe, -1.0, 0.0, 1.0, probe])
        vs = np.asarray(potential_func(xs), dtype=float)
    except Exception:
        return None
    if not np.all(np.isfinite(vs)):
        return "V(x) is non-finite somewhere on the probe grid -- check the potential formula for this coefficient value."
    v_far = min(vs[0], vs[-1])
    v_near = max(vs[1:4])
    if v_far < v_near:
        return (
            f"V(x) at x=±{probe:g} is {v_far:.3g}, not above the near-origin "
            f"values (~{v_near:.3g}) -- this potential is likely unbounded from "
            f"below (or otherwise non-confining) for this coefficient value."
        )
    return None


# =====================================================================
# Per-variant DVR solve (one attempt, no escalation)
# =====================================================================
def solve_variant_spectrum(my_potential, variant_params, num_states, mass=1.0, hbar=1.0):
    """
    Solve for one variant potential's own energy spectrum, with the
    same rigor as the base run's Section 1 & 4: auto-configured grid,
    then the 3-pass (base/resolution/boundary-span) convergence check.

    Parameters
    ----------
    my_potential : callable
        config.my_potential -- must accept `p=` to select which
        parameter dict it reads its coefficients from.
    variant_params : dict
        This variant's POTENTIAL_PARAMS.
    num_states : int
        Number of energy levels to solve for.
    mass, hbar : float, optional

    Returns
    -------
    energies, x_min, x_max : ndarray, float, float
    """
    potential_func = functools.partial(my_potential, p=variant_params)
    x_min, x_max, n_grid = auto_configure_dvr(potential_func, num_states, mass=mass, hbar=hbar)
    energies = get_fully_converged_energy_levels(
        potential_func=potential_func, num_levels=num_states,
        x_min=x_min, x_max=x_max, num_points=n_grid,
        mass=mass, hbar=hbar,
    )
    return energies, x_min, x_max


# =====================================================================
# Per-variant DVR solve + direct quantum Cv(T) evaluation (no escalation)
# =====================================================================
def compute_variant_quantum_cv(my_potential, variant_params, num_states,
                                beta_arr, mass=1.0, hbar=1.0):
    """
    `solve_variant_spectrum` followed by a direct quantum Cv(T)
    evaluation over `beta_arr` -- no xi/n-convergence search, no
    hot-coverage escalation (see `solve_variant_with_hot_coverage` for
    the escalation-aware version `run_coefficient_sweep` actually uses).
    Kept as a simple, single-attempt building block.

    Parameters
    ----------
    my_potential : callable
    variant_params : dict
    num_states : int
    beta_arr : array_like
    mass, hbar : float, optional

    Returns
    -------
    cv_curve : ndarray, shape (len(beta_arr),)
    """
    energies, _, _ = solve_variant_spectrum(my_potential, variant_params, num_states, mass, hbar)
    return compute_quantum_heat_capacity_curve(energies, beta_arr, xi=1.0)


# =====================================================================
# Per-variant DVR solve with hot-end thermal-coverage escalation
# =====================================================================
def solve_variant_with_hot_coverage(my_potential, variant_params, num_states, beta_arr,
                                     mass, hbar, hot_state_safety,
                                     num_states_growth, num_states_cap, max_rounds):
    """
    `solve_variant_spectrum`, escalating this ONE variant's own
    `num_states` if its top computed level doesn't clear the same
    thermal-truncation safety margin `Cv_AutoTune.py` uses for the
    base run -- see this module's docstring, "IS IT VALID TO REUSE THE
    BASE RUN'S NUM_STATES FOR EVERY VARIANT?", for the physical
    rationale. Bounded by `max_rounds`/`num_states_cap`, same as the
    base run's own escalation loop.

    Parameters
    ----------
    my_potential : callable
    variant_params : dict
    num_states : int
        Starting guess (normally the base run's final, post-escalation
        NUM_STATES).
    beta_arr : array_like
        Shared inverse-temperature array; only its minimum (T_hot) is
        used here.
    mass, hbar : float
    hot_state_safety : float
        Target ratio of E_max to k_B*T_hot (config.HOT_STATE_SAFETY).
    num_states_growth : float
        Growth factor applied to `num_states` each escalation round
        (config.NUM_STATES_GROWTH).
    num_states_cap : int
        Upper bound on `num_states` (config.NUM_STATES_CAP).
    max_rounds : int
        Maximum escalation rounds (config.MAX_ESCALATION_ROUNDS).

    Returns
    -------
    dict with keys:
        energies, x_min, x_max : ndarray, float, float
        num_states_used : int
        coverage_ok : bool
            False if every round still failed the safety margin (the
            curve is still returned -- see module docstring).
    """
    beta_min = float(np.min(beta_arr))
    n = num_states
    energies = x_min = x_max = None
    coverage_ok = False
    for round_i in range(max_rounds):
        energies, x_min, x_max = solve_variant_spectrum(my_potential, variant_params, n, mass, hbar)
        coverage_ok = float(energies[-1]) >= hot_state_safety / beta_min
        if coverage_ok or n >= num_states_cap or round_i == max_rounds - 1:
            break
        n_next = min(int(n * num_states_growth), num_states_cap)
        print(f"    hot-end coverage marginal (E_max={energies[-1]:.3g}, need "
              f"≥{hot_state_safety / beta_min:.3g}) -- escalating this "
              f"variant's NUM_STATES {n} → {n_next}")
        n = n_next
    return {
        "energies": energies, "x_min": x_min, "x_max": x_max,
        "num_states_used": n, "coverage_ok": coverage_ok,
    }


# =====================================================================
# Shared color assignment (kept identical across the Cv and potential
# comparison plots so a variant's color means the same thing in both)
# =====================================================================
def _variant_colors(n):
    """tab10 (clearly distinguishable) for a typical small sweep;
    viridis (perceptually ordered) once there are more variants than
    tab10 has distinct colors for."""
    if n <= 10:
        cmap = cm.get_cmap("tab10")
        return [cmap(i) for i in range(n)]
    cmap = cm.get_cmap("viridis")
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


# =====================================================================
# Plot: reused classical limit + one quantum curve per variant
# =====================================================================
def plot_coefficient_sweep(T_arr, cv_classical_base, variant_values, variant_curves,
                            scan_param, system_name, formula_text=None,
                            marginal_values=None, T_units_label=r"$k_B T \,/\, E_0$"):
    """
    Plot the base run's classical-limit Cv(T) (reused, unchanged)
    together with one quantum Cv(T) curve per coefficient variant.
    Deliberately excludes everything else (xi/n convergence, secondary
    axes) -- this figure is only about comparing the variants' quantum
    curves against one shared classical-limit reference. `variant_values`/
    `variant_curves` may be shorter than the full requested sweep (see
    `run_coefficient_sweep`) -- variants that failed to converge are
    simply not plotted; even an empty list still produces a valid
    figure showing just the classical-limit reference.

    Parameters
    ----------
    T_arr : ndarray
        Temperature axis, shared with the base run.
    cv_classical_base : ndarray
        The base run's already-computed numerical classical limit.
    variant_values : list of float
        This variant's value of `scan_param`, one per curve, used for
        the legend and the color mapping (ascending order). May be
        shorter than the requested sweep if some variants failed.
    variant_curves : list of ndarray
        Quantum Cv(T) curves, one per entry in `variant_values`.
    scan_param : str
        Name of the varied coefficient (for the legend/title).
    system_name : str
    formula_text : str or None, optional
        The potential's formula (see `format_potential_formula`), shown
        as a second title line. None omits that line entirely.
    marginal_values : set of float or None, optional
        Values whose hot-end thermal coverage stayed marginal even
        after escalation (see `solve_variant_with_hot_coverage`) --
        flagged with a "*" in the legend and a footnote, rather than
        silently plotted as if fully trustworthy at every T.
    T_units_label : str, optional

    Returns
    -------
    None (saves the figure; see figures/output_paths.py).
    """
    # Reserved for the classical-limit reference line only -- deliberately
    # NOT reused for any variant curve (tab10's own green, #2ca02c, would
    # otherwise collide with whichever variant lands on that palette slot
    # and become invisible underneath the dashed reference line).
    CLASSICAL_LIMIT_COLOR = "#000000"
    marginal_values = marginal_values or set()
    fig, ax = plt.subplots(figsize=(9, 6))
    title_lines = [f"{system_name} — Quantum $C_v(T)$ vs {scan_param}"]
    if formula_text:
        title_lines.append(formula_text)
    title_lines.append("(shared classical limit reused from the base run)")
    fig.suptitle("\n".join(title_lines), fontsize=13, fontweight="bold")

    colors = _variant_colors(len(variant_values))
    for value, curve, color in zip(variant_values, variant_curves, colors):
        flag = " *" if value in marginal_values else ""
        ax.plot(T_arr, curve, color=color, linewidth=1.8,
                label=f"{scan_param} = {value:g}{flag}")

    ax.plot(T_arr, cv_classical_base, color=CLASSICAL_LIMIT_COLOR, linewidth=2.2, linestyle="--",
            label="Numerical classical limit (base run)")

    ax.set_xlabel(T_units_label, fontsize=12)
    ax.set_ylabel(r"$C_v \,/\, k_B$", fontsize=12)
    ax.set_xscale("log")
    legend_title = "* hot-end coverage marginal (see console)" if marginal_values else None
    ax.legend(fontsize=9, loc="upper left", title=legend_title, title_fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    # tight_layout() doesn't reserve room for a multi-line suptitle by
    # itself; push the axes down afterward so the (up to 3-line) title
    # never overlaps the plotted curves.
    fig.subplots_adjust(top=0.99 - 0.06 * len(title_lines))
    save_figure(fig, "coefficient_sweep", "cv_coefficient_sweep")


# =====================================================================
# Render one V(x) + low-lying-spectrum panel onto a given Axes
# =====================================================================
def _render_potential_panel(ax, record, levels_to_draw, color=None, min_levels_shown=3):
    """
    Draw one variant's potential and its lowest `levels_to_draw` energy
    levels onto `ax`, zoomed on the well's own structure (reusing
    `figures.plot_potential`'s origin-zoom window and turning-point
    logic) -- this is the view that actually shows the barrier/
    asymmetry structure responsible for a Schottky-anomaly-like bump,
    as opposed to a full-spectrum overview.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    record : dict
        One entry from `run_coefficient_sweep`'s per-variant records:
        must contain "potential_func", "energies", "x_min", "x_max".
    levels_to_draw : int
        How many of the lowest computed levels to overlay.
    color : str or None, optional
        Color for the V(x) curve itself (level markers stay a fixed
        contrasting color regardless, so they read the same across
        panels no matter which variant color they belong to).
    min_levels_shown : int, optional
        `_origin_zoom_window`'s barrier-based crop (calibrated for a
        single full-detail figure) can end up tighter than a shallow
        well's own ground state, leaving a panel with zero levels drawn
        -- which would defeat a side-by-side comparison. The window is
        widened (never shrunk) to include at least this many of the
        lowest levels, if they exist (default 3).

    Returns
    -------
    int
        Number of levels actually drawn.
    """
    potential_func = record["potential_func"]
    levels = np.sort(np.asarray(record["energies"]))[:levels_to_draw]
    x = np.linspace(record["x_min"], record["x_max"], 2000)
    V = np.asarray(potential_func(x), dtype=float)
    v_bottom = float(V.min())

    x_lo, x_hi, y_bottom, y_top = _origin_zoom_window(x, V, v_bottom, levels)
    if len(levels):
        floor_idx = min(min_levels_shown, len(levels)) - 1
        if levels[floor_idx] > y_top:
            y_top = levels[floor_idx] + 0.05 * (levels[floor_idx] - v_bottom)
            y_bottom = v_bottom - 0.05 * (y_top - v_bottom)
    x_pad = 0.08 * (x_hi - x_lo)

    ax.plot(x, V, color=color or "#1f4e8c", linewidth=2.0, zorder=3)
    ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
    ax.set_ylim(y_bottom, y_top)

    n_shown = 0
    for E in levels:
        if E > y_top:
            continue
        tp = classical_turning_points(potential_func, E, x)
        if tp is None:
            continue
        xl, xr = tp
        ax.hlines(E, xl, xr, color="#d62728", linewidth=0.9, alpha=0.75, zorder=2)
        n_shown += 1
    ax.set_xlabel("$x$", fontsize=9)
    ax.grid(alpha=0.3, linestyle="--")
    return n_shown


# =====================================================================
# Potential comparison: side-by-side panels (small variant counts)
# =====================================================================
def _plot_potentials_side_by_side(records, scan_param, system_name, levels_to_draw, formula_text=None):
    n = len(records)
    fig, axes = plt.subplots(1, n, figsize=(4.3 * n, 5), squeeze=False)
    axes = axes[0]
    colors = _variant_colors(n)
    for ax, record, color in zip(axes, records, colors):
        n_shown = _render_potential_panel(ax, record, levels_to_draw, color)
        ax.set_title(f"{scan_param} = {record['value']:g}\n({n_shown} levels shown)",
                     fontsize=10, fontweight="bold", color=color)
    axes[0].set_ylabel("Energy", fontsize=10)
    title_lines = [f"{system_name} — potential & spectrum vs {scan_param}"]
    if formula_text:
        title_lines.append(formula_text)
    fig.suptitle("\n".join(title_lines), fontsize=13, fontweight="bold")
    top = 0.90 if formula_text is None else 0.84
    plt.tight_layout(rect=[0, 0, 1, top])
    return save_figure(fig, "coefficient_sweep", "potential_comparison")


# =====================================================================
# Potential comparison: single overlay (medium variant counts)
# =====================================================================
def _plot_potentials_overlay(records, scan_param, system_name, formula_text=None):
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = _variant_colors(len(records))

    curves, windows = [], []
    for record in records:
        x = np.linspace(record["x_min"], record["x_max"], 2000)
        V = np.asarray(record["potential_func"](x), dtype=float)
        levels = np.sort(np.asarray(record["energies"]))[:5]
        windows.append(_origin_zoom_window(x, V, float(V.min()), levels))
        curves.append((x, V))

    # Shared window = the union of every variant's own zoom window, so
    # every well's own structure stays visible on one shared axes.
    x_lo = min(w[0] for w in windows)
    x_hi = max(w[1] for w in windows)
    y_bottom = min(w[2] for w in windows)
    y_top = max(w[3] for w in windows)

    for record, (x, V), color in zip(records, curves, colors):
        ax.plot(x, V, color=color, linewidth=1.8, label=f"{scan_param} = {record['value']:g}")

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_bottom, y_top)
    ax.set_xlabel("$x$", fontsize=11)
    ax.set_ylabel("Energy", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, linestyle="--")
    title_lines = [f"{system_name} — potential overlay vs {scan_param}"]
    if formula_text:
        title_lines.append(formula_text)
    title_lines.append("(individual levels omitted -- too many variants for side-by-side panels)")
    fig.suptitle("\n".join(title_lines), fontsize=13, fontweight="bold")
    top = 0.90 if formula_text is None else 0.84
    plt.tight_layout(rect=[0, 0, 1, top])
    return save_figure(fig, "coefficient_sweep", "potential_comparison")


# =====================================================================
# Potential comparison: one figure per variant (large variant counts)
# =====================================================================
def _plot_potentials_separate(records, scan_param, system_name, levels_to_draw, formula_template=None):
    paths = []
    for record in records:
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        n_shown = _render_potential_panel(ax, record, levels_to_draw)
        ax.set_ylabel("Energy", fontsize=11)
        # Fully numeric (not symbolic) here -- each figure is exactly one
        # concrete potential, so scan_param=None substitutes every
        # coefficient (including the swept one) with this record's own value.
        formula_text = format_potential_formula(formula_template, record["params"], None)
        title_lines = [f"{system_name} — {scan_param} = {record['value']:g}"]
        if formula_text:
            title_lines.append(formula_text)
        title_lines.append(f"({n_shown} levels shown)")
        ax.set_title("\n".join(title_lines), fontsize=12, fontweight="bold")
        plt.tight_layout()
        value_slug = f"{record['value']:g}".replace("-", "m").replace(".", "p")
        paths.append(save_figure(fig, "coefficient_sweep", f"potential_{scan_param}_{value_slug}"))
    return paths


# =====================================================================
# Potential comparison: dispatcher
# =====================================================================
def plot_variant_potentials(records, scan_param, system_name, formula_template=None,
                             base_params=None, levels_to_draw=25):
    """
    Compare the swept potentials and their low-lying spectra -- this
    is what lets a Cv anomaly's size be correlated with the actual
    change in well shape that caused it. Falls back through three
    layouts as the variant count grows, since all three become
    unreadable at some point in the opposite direction:

        variants <= SIDE_BY_SIDE_MAX (6)  -> one figure, one subplot
            per variant, each independently zoomed + its own levels
            drawn (most detail, needs the most horizontal space).
        SIDE_BY_SIDE_MAX < variants <= OVERLAY_MAX (15) -> one figure,
            all V(x) curves overlaid on shared axes, colors matching
            the Cv plot's legend (no individual levels -- would be
            unreadable overlaid).
        variants > OVERLAY_MAX -> one figure per variant (most figures,
            least clutter per figure).

    Every layout's title reflects WHICHEVER coefficient is actually
    being swept -- there is nothing hard-coded to any one coefficient
    name (`scan_param` drives every label, filename, and formula
    substitution) -- so switching `config.SCAN_PARAM` to a different
    key needs no changes here. The side-by-side/overlay titles show
    the formula with the swept coefficient symbolic (matching the Cv
    plot); the separate-figure layout shows each one fully numeric,
    since each of those figures is exactly one concrete potential.

    Parameters
    ----------
    records : list of dict
        One entry per CONVERGED variant (see `run_coefficient_sweep`),
        each with "value", "params", "potential_func", "energies",
        "x_min", "x_max".
    scan_param : str
    system_name : str
    formula_template : str or None, optional
        config.POTENTIAL_FORMULA -- see `format_potential_formula`.
        None omits the formula annotation from every layout.
    base_params : dict or None, optional
        config.POTENTIAL_PARAMS, needed (together with `formula_template`)
        for the side-by-side/overlay layouts' symbolic formula line.
    levels_to_draw : int, optional
        How many of each variant's lowest levels to overlay in the
        side-by-side / separate-figure layouts (default 25 -- enough
        to show the low-lying structure without cluttering a subplot;
        unused in overlay mode, which omits individual levels).

    Returns
    -------
    str or list of str or None
        Path(s) of the figure(s) saved; None if `records` is empty.
    """
    n = len(records)
    if n == 0:
        return None
    if n <= SIDE_BY_SIDE_MAX:
        formula_text = format_potential_formula(formula_template, base_params, scan_param) if base_params else None
        return _plot_potentials_side_by_side(records, scan_param, system_name, levels_to_draw, formula_text)
    if n <= OVERLAY_MAX:
        formula_text = format_potential_formula(formula_template, base_params, scan_param) if base_params else None
        return _plot_potentials_overlay(records, scan_param, system_name, formula_text)
    return _plot_potentials_separate(records, scan_param, system_name, levels_to_draw, formula_template)


# =====================================================================
# Orchestrator
# =====================================================================
def run_coefficient_sweep(base_cv_results, my_potential, base_params,
                           scan_param, scan_step, scan_count, num_states,
                           mass, hbar, system_name, formula_template=None,
                           hot_state_safety=20.0, num_states_growth=1.6,
                           num_states_cap=4000, max_escalation_rounds=4,
                           T_units_label=r"$k_B T \,/\, E_0$"):
    """
    Generate the coefficient variants, solve each one's own spectrum
    (escalating that variant's own NUM_STATES if its hot-end thermal
    coverage is marginal -- see `solve_variant_with_hot_coverage` and
    this module's docstring), compute its quantum Cv(T) directly (no
    xi/n-convergence search), and produce both the Cv comparison plot
    and the potential/spectrum comparison plot(s).

    A variant whose DVR solve fails outright (typically: this
    coefficient value makes the potential non-confining, e.g. a
    quartic leading term that went negative) does NOT abort the sweep
    -- its failure is caught, diagnosed (see `_diagnose_variant_failure`),
    and printed; the figures are still produced from whichever variants
    DID converge. If every variant fails, the Cv figure still saves
    (showing just the classical-limit reference) and a clear warning
    is printed.

    Parameters
    ----------
    base_cv_results : dict
        Output of `Quantum_Classical_Combined.run()` for the base
        potential. Must contain "T_arr" and "cv_classical".
    my_potential : callable
    base_params : dict
        config.POTENTIAL_PARAMS.
    scan_param : str
    scan_step : float
    scan_count : int
    num_states : int
        The base run's final NUM_STATES (post auto-tune escalation) --
        the STARTING guess for every variant's own solve.
    mass, hbar : float
    system_name : str
    formula_template : str or None, optional
        config.POTENTIAL_FORMULA -- see `format_potential_formula`.
        None omits the formula annotation from the Cv plot.
    hot_state_safety, num_states_growth, num_states_cap, max_escalation_rounds :
        Passed through to `solve_variant_with_hot_coverage` for each
        variant -- normally config.HOT_STATE_SAFETY/NUM_STATES_GROWTH/
        NUM_STATES_CAP/MAX_ESCALATION_ROUNDS, the same knobs the base
        run's own auto-tune loop uses.
    T_units_label : str, optional

    Returns
    -------
    dict with keys:
        variant_params, variant_curves : list of dict / ndarray
            Only the variants that converged, in swept order.
        failed : list of dict
            One {"params", "value", "error", "diagnosis"} entry per
            variant that failed to converge outright.
        marginal_values : list of float
            Values whose hot-end thermal coverage stayed marginal even
            after escalation (still plotted, but flagged).
        potential_figure_paths : str, list of str, or None
            Whatever `plot_variant_potentials` returned.
    """
    beta_arr = base_cv_results["beta_arr"]
    T_arr = base_cv_results["T_arr"]

    variant_params = generate_variant_params(base_params, scan_param, scan_step, scan_count)
    all_values = [p[scan_param] for p in variant_params]

    rule = "─" * 60
    print(f"\n{rule}\n  Coefficient sweep: {scan_param} over {all_values}\n{rule}")

    ok_params, ok_values, ok_curves, ok_records = [], [], [], []
    failed, marginal_values = [], []
    for i, params in enumerate(variant_params):
        value = params[scan_param]
        print(f"  [{i + 1}/{len(variant_params)}] {scan_param} = {value:g} ...", flush=True)
        try:
            result = solve_variant_with_hot_coverage(
                my_potential, params, num_states, beta_arr, mass, hbar,
                hot_state_safety, num_states_growth, num_states_cap, max_escalation_rounds,
            )
        except Exception as exc:
            potential_func = functools.partial(my_potential, p=params)
            diagnosis = _diagnose_variant_failure(potential_func)
            print(f"  ✗ {scan_param} = {value:g} did NOT converge: {exc}")
            if diagnosis:
                print(f"    diagnosis: {diagnosis}")
            failed.append({"params": params, "value": value, "error": str(exc), "diagnosis": diagnosis})
            continue

        energies = result["energies"]
        curve = compute_quantum_heat_capacity_curve(energies, beta_arr, xi=1.0)
        if not result["coverage_ok"]:
            marginal_values.append(value)
            print(f"    ⚠ {scan_param} = {value:g} kept, but hot-end coverage stayed marginal "
                  f"after escalation (NUM_STATES={result['num_states_used']}) -- treat T above "
                  f"~{1.0 / beta_arr.min():.3g} with caution for this curve.")

        ok_params.append(params)
        ok_values.append(value)
        ok_curves.append(curve)
        ok_records.append({
            "value": value, "params": params, "energies": energies,
            "x_min": result["x_min"], "x_max": result["x_max"],
            "potential_func": functools.partial(my_potential, p=params),
            "coverage_ok": result["coverage_ok"], "num_states_used": result["num_states_used"],
        })

    formula_text = format_potential_formula(formula_template, base_params, scan_param)
    plot_coefficient_sweep(
        T_arr, base_cv_results["cv_classical"], ok_values, ok_curves,
        scan_param, system_name, formula_text, set(marginal_values), T_units_label,
    )
    potential_figure_paths = plot_variant_potentials(
        ok_records, scan_param, system_name, formula_template, base_params,
    )

    if failed:
        print(f"  ⚠ {len(failed)}/{len(variant_params)} coefficient values did not converge "
              f"and are OMITTED from both plots: {[f['value'] for f in failed]}")
        if not ok_values:
            print("  ⚠ NO variant converged -- the Cv plot shows only the classical-limit reference, "
                  "and no potential-comparison figure was produced.")
    if marginal_values:
        print(f"  ⚠ {len(marginal_values)} coefficient value(s) kept with marginal hot-end coverage "
              f"(flagged with * in the Cv plot legend): {marginal_values}")
    print(f"  ✓ Coefficient sweep figures saved to the coefficient_sweep/ folder.\n{rule}\n")

    return {
        "variant_params": ok_params, "variant_curves": ok_curves,
        "failed": failed, "marginal_values": marginal_values,
        "potential_figure_paths": potential_figure_paths,
    }
