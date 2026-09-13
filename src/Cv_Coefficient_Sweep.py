"""
Cv_Coefficient_Sweep.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Produces exactly one figure: the quantum Cv(T) curve for several
variants of the run's potential that differ only in one named
POTENTIAL_PARAMS coefficient (e.g. the double well's cubic term `b`),
overlaid with the single classical-limit Cv(T) curve already computed
for the base potential in config.py, and annotated with the actual
potential formula (config.POTENTIAL_FORMULA, fixed coefficients shown
numerically, the swept one shown symbolically). This is what lets a
coefficient's effect on the size/shape of a (numerical) Schottky
anomaly be compared directly across several potentials, without
changing anything about the rest of Quantum_HO_Master.py's pipeline
for the base potential itself -- this file only ever adds Section 7,
its one new figure.

WHY NO XI/N-CONVERGENCE SEARCH HERE
---------------------------------------------------------------------
The classical-limit curve is reused verbatim from the base run's
`base_cv_results["cv_classical"]` -- it is NOT recomputed per variant.
Each variant's quantum Cv(T) is a direct evaluation of the literal,
unscaled quantum heat capacity (xi=1.0, no rescaling) over the base
run's own `beta_arr`, using that variant's own DVR spectrum computed
with the base run's own (possibly auto-tune-escalated) NUM_STATES.
No xi-scan or n-convergence search is needed for a direct quantum Cv
evaluation -- that machinery exists only to locate the classical-limit
plateau, which this file deliberately does not redo per variant.

WHAT HAPPENS WHEN A VARIANT FAILS TO SOLVE
---------------------------------------------------------------------
Not every coefficient value produces a genuinely confining potential
(e.g. a quartic leading coefficient that goes negative is unbounded
from below) -- the DVR grid auto-configurator or its 3-pass
convergence check will then fail for that one variant. `run_coefficient_sweep`
catches that per variant, prints a diagnostic (including a cheap,
potential-agnostic unboundedness check -- see `_diagnose_variant_failure`),
and keeps going: the final plot is built from whichever variants DID
converge, with a printed summary of which ones didn't and why, rather
than losing the whole sweep to one bad coefficient value.
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
# Per-variant DVR solve + direct quantum Cv(T) evaluation
# =====================================================================
def compute_variant_quantum_cv(my_potential, variant_params, num_states,
                                beta_arr, mass=1.0, hbar=1.0):
    """
    Solve for one variant potential's own energy spectrum (same rigor
    as the base run's Section 1: auto-configured grid + 3-pass
    convergence check), then directly evaluate the literal quantum
    Cv(T) over `beta_arr` -- no xi/n-convergence search.

    Parameters
    ----------
    my_potential : callable
        config.my_potential -- must accept `p=` to select which
        parameter dict it reads its coefficients from.
    variant_params : dict
        This variant's POTENTIAL_PARAMS.
    num_states : int
        Number of energy levels to compute (reused from the base run).
    beta_arr : array_like
        Inverse-temperature array (reused from the base run).
    mass, hbar : float, optional

    Returns
    -------
    cv_curve : ndarray, shape (len(beta_arr),)
    """
    potential_func = functools.partial(my_potential, p=variant_params)
    x_min, x_max, n_grid = auto_configure_dvr(potential_func, num_states, mass=mass, hbar=hbar)
    energies = get_fully_converged_energy_levels(
        potential_func=potential_func, num_levels=num_states,
        x_min=x_min, x_max=x_max, num_points=n_grid,
        mass=mass, hbar=hbar,
    )
    return compute_quantum_heat_capacity_curve(energies, beta_arr, xi=1.0)


# =====================================================================
# Plot: reused classical limit + one quantum curve per variant
# =====================================================================
def plot_coefficient_sweep(T_arr, cv_classical_base, variant_values, variant_curves,
                            scan_param, system_name, formula_text=None,
                            T_units_label=r"$k_B T \,/\, E_0$"):
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
    fig, ax = plt.subplots(figsize=(9, 6))
    title_lines = [f"{system_name} — Quantum $C_v(T)$ vs {scan_param}"]
    if formula_text:
        title_lines.append(formula_text)
    title_lines.append("(shared classical limit reused from the base run)")
    fig.suptitle("\n".join(title_lines), fontsize=13, fontweight="bold")

    # tab10 gives clearly distinguishable colors for a typical (small)
    # sweep; viridis is used instead once there are more variants than
    # tab10 has distinct colors for, where perceptual ordering matters
    # more than per-curve contrast.
    n = len(variant_values)
    if n <= 10:
        cmap = cm.get_cmap("tab10")
        colors = [cmap(i) for i in range(n)]
    else:
        cmap = cm.get_cmap("viridis")
        colors = [cmap(i / max(n - 1, 1)) for i in range(n)]
    for value, curve, color in zip(variant_values, variant_curves, colors):
        ax.plot(T_arr, curve, color=color, linewidth=1.8,
                label=f"{scan_param} = {value:g}")

    ax.plot(T_arr, cv_classical_base, color=CLASSICAL_LIMIT_COLOR, linewidth=2.2, linestyle="--",
            label="Numerical classical limit (base run)")

    ax.set_xlabel(T_units_label, fontsize=12)
    ax.set_ylabel(r"$C_v \,/\, k_B$", fontsize=12)
    ax.set_xscale("log")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    # tight_layout() doesn't reserve room for a multi-line suptitle by
    # itself; push the axes down afterward so the (up to 3-line) title
    # never overlaps the plotted curves.
    fig.subplots_adjust(top=0.99 - 0.06 * len(title_lines))
    save_figure(fig, "cv", "cv_coefficient_sweep")


# =====================================================================
# Orchestrator
# =====================================================================
def run_coefficient_sweep(base_cv_results, my_potential, base_params,
                           scan_param, scan_step, scan_count, num_states,
                           mass, hbar, system_name, formula_template=None,
                           T_units_label=r"$k_B T \,/\, E_0$"):
    """
    Generate the coefficient variants, compute each one's quantum
    Cv(T) directly (no new convergence search), and plot the ones that
    converged against the base run's already-computed classical limit.

    A variant whose DVR solve fails (typically: this coefficient value
    makes the potential non-confining, e.g. a quartic leading term that
    went negative) does NOT abort the sweep -- its failure is caught,
    diagnosed (see `_diagnose_variant_failure`), and printed; the
    figure is still produced from whichever variants DID converge. If
    every variant fails, the figure still saves (showing just the
    classical-limit reference) and a clear warning is printed.

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
        The base run's final NUM_STATES (post auto-tune escalation).
    mass, hbar : float
    system_name : str
    formula_template : str or None, optional
        config.POTENTIAL_FORMULA -- see `format_potential_formula`.
        None omits the formula annotation from the plot.
    T_units_label : str, optional

    Returns
    -------
    dict with keys:
        variant_params, variant_curves : list of dict / ndarray
            Only the variants that converged, in swept order.
        failed : list of dict
            One {"params", "value", "error", "diagnosis"} entry per
            variant that failed to converge.
    """
    beta_arr = base_cv_results["beta_arr"]
    T_arr = base_cv_results["T_arr"]

    variant_params = generate_variant_params(base_params, scan_param, scan_step, scan_count)
    all_values = [p[scan_param] for p in variant_params]

    rule = "─" * 60
    print(f"\n{rule}\n  Coefficient sweep: {scan_param} over {all_values}\n{rule}")

    ok_params, ok_values, ok_curves = [], [], []
    failed = []
    for i, params in enumerate(variant_params):
        value = params[scan_param]
        print(f"  [{i + 1}/{len(variant_params)}] {scan_param} = {value:g} ...", flush=True)
        try:
            curve = compute_variant_quantum_cv(my_potential, params, num_states, beta_arr, mass, hbar)
        except Exception as exc:
            potential_func = functools.partial(my_potential, p=params)
            diagnosis = _diagnose_variant_failure(potential_func)
            print(f"  ✗ {scan_param} = {value:g} did NOT converge: {exc}")
            if diagnosis:
                print(f"    diagnosis: {diagnosis}")
            failed.append({"params": params, "value": value, "error": str(exc), "diagnosis": diagnosis})
            continue
        ok_params.append(params)
        ok_values.append(value)
        ok_curves.append(curve)

    formula_text = format_potential_formula(formula_template, base_params, scan_param)
    plot_coefficient_sweep(
        T_arr, base_cv_results["cv_classical"], ok_values, ok_curves,
        scan_param, system_name, formula_text, T_units_label,
    )

    if failed:
        print(f"  ⚠ {len(failed)}/{len(variant_params)} coefficient values did not converge "
              f"and are OMITTED from the plot: {[f['value'] for f in failed]}")
        if not ok_values:
            print("  ⚠ NO variant converged -- the plot shows only the classical-limit reference.")
    print(f"  ✓ Coefficient sweep figure saved.\n{rule}\n")

    return {"variant_params": ok_params, "variant_curves": ok_curves, "failed": failed}

    return {"variant_params": variant_params, "variant_curves": variant_curves}
