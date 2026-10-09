"""
Quick_Scan.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
A fast first look at a potential before committing the full pipeline
(Quantum_HO_Master.py) to it. It produces ONE figure -- the quantum
Cv(T) (solid) and the classical-limit Cv(T) (dashed, same color) -- for

    "single"  the potential in config.py, or
    "sweep"   every variant of config.SCAN_PARAM's sweep (SCAN_STEP,
              SCAN_COUNT, SCAN_SYMMETRIC_VALUE, exactly as Section 7),
              each against its own classical limit,

plus a console verdict per potential on whether the quantum Cv is above
the classical one. No convergence figures, no reference grids, no
auto-tune loop, no potential/spectrum panels.

THE SAME METHOD AS THE PIPELINE, AT A LOWER RESOLUTION
---------------------------------------------------------------------
Nothing here is a new algorithm. Every number comes from the functions
the full pipeline uses:

    classical limit   Classical_Limit_Numerical.sweep_temperature_range:
                      the xi-scan with xi^2 V re-solved at every xi
                      (Gelbwaser Eq. S7), NaN where no plateau is reached
    quantum Cv        Quantum_Classical_Combined.compute_quantum_heat_
                      capacity_curve on the spectrum of V itself
    variants          Cv_Coefficient_Sweep.generate_variant_params
    figure            Cv_Coefficient_Sweep.plot_coefficient_sweep

Only the settings are coarser, chosen by config.QUICK_SCAN_RESOLUTION
(a key of config.QUICK_SCAN_PRESETS): fewer temperatures, a looser
tolerance on the classical value (tol_xi), a coarser xi ladder with a
shorter plateau, fewer levels per solve (thermal_coverage) and a looser
DVR check. The tolerance is what saves the time. The cold end needs the
largest xi, the xi needed grows as ~1/sqrt(tol_xi), the DVR grid of
xi^2 V grows ~xi and a solve costs ~grid^3. Resolution 3 is the full
pipeline's own settings.

THE QUANTUM SPECTRUM COSTS NOTHING EXTRA
---------------------------------------------------------------------
The xi ladder starts at xi = 1 (config.XI_START), so the scan's first
solve -- xi = 1 at the hottest temperature -- is the spectrum of V
itself: the same DVR with its 3-pass check, holding every level up to
E_0 + thermal_coverage * k_B T_hot. The quantum Cv is evaluated on that
cached spectrum. On the double well at resolution 1 it agrees with the
full pipeline's 500-level quantum Cv to 1.5e-7.

READING THE VERDICT: IS THE QUANTUM CV ABOVE THE CLASSICAL ONE?
---------------------------------------------------------------------
Each classical value carries the xi-scan's own estimate of its distance
to the xi -> inf limit, eps = |dCv| / (xi_mult^2 - 1) (see converge_xi).
With d = Cv_quantum - Cv_classical at each temperature:

    d >  2 eps     quantum ABOVE classical, resolved
    |d| <= 2 eps   not resolved at this resolution
    d < -2 eps     quantum below classical, resolved

Why 2 eps (ERROR_MARGIN): the estimate tracks the true error to within
~10% (audit/quick_scan/), but a coarse classical value sits BELOW the
limit by about eps, which shifts d upward by about eps. With a margin of
1, a temperature where the curves have merged could be called "above";
with 2, a wrong verdict would need the error to exceed twice its
estimate.

At high T the quantum Cv approaches the classical one, so d -> 0 and a
band of unresolved temperatures at the hot end is expected at any
resolution; the verdict says when that is the only unresolved region.
Unresolved points anywhere else, or gaps in a classical curve (its
ladder ran out), are the cue to raise QUICK_SCAN_RESOLUTION, or to zoom
QUICK_SCAN_BETA_RANGE onto that window and raise it there (cheaper,
especially when the window leaves out the cold end).

USAGE
---------------------------------------------------------------------
Set the potential and the QUICK_SCAN_* settings in config.py, then from
src/:

    python Quick_Scan.py                              # config's settings
    python Quick_Scan.py --mode single --resolution 2
    python Quick_Scan.py --resolution 2 --beta-range 0.5 5

The figure is saved to figures/<system>/<params>/quick_scan/ (see
figures/output_paths.py), named after the mode, resolution and any
zoom window, so a quick scan never overwrites the pipeline's figures.
=====================================================================
"""

import argparse
import functools
import math
import time

import numpy as np

import config
from Classical_Limit_Numerical import ScaledSpectra, sweep_temperature_range
from Quantum_Classical_Combined import compute_quantum_heat_capacity_curve
from Cv_AutoTune import resolve_beta_min
from Cv_Coefficient_Sweep import (generate_variant_params, format_potential_formula,
                                  plot_coefficient_sweep, _diagnose_variant_failure)
from figures.output_paths import set_context as set_figure_context

# A difference d = Cv_q - Cv_cl counts as resolved only where |d| exceeds
# this many times the classical value's error estimate (module docstring).
ERROR_MARGIN = 2.0


# =====================================================================
# Settings for one resolution
# =====================================================================
def quick_settings(resolution):
    """
    config.QUICK_SCAN_PRESETS[resolution], plus the ladder length that
    takes xi from config.XI_START up to the pipeline's cap config.XI_MAX
    at this preset's xi_mult.
    """
    if resolution not in config.QUICK_SCAN_PRESETS:
        raise KeyError(f"QUICK_SCAN_RESOLUTION = {resolution!r} is not a key of "
                       f"QUICK_SCAN_PRESETS {sorted(config.QUICK_SCAN_PRESETS)}.")
    settings = dict(config.QUICK_SCAN_PRESETS[resolution])
    settings["max_xi_steps"] = int(math.ceil(
        math.log(config.XI_MAX / config.XI_START) / math.log(settings["xi_mult"]))) + 1
    return settings


def _new_spectra(potential_func, settings):
    """A ScaledSpectra cache with this resolution's level coverage and DVR tolerance."""
    return ScaledSpectra(potential_func, mass=config.MASS, hbar=config.HBAR,
                         thermal_coverage=settings["thermal_coverage"],
                         dvr_tolerance=settings["dvr_tolerance"])


# =====================================================================
# Temperature grid (shared by every variant, as in Section 7)
# =====================================================================
def temperature_grid(base_spectra, n_beta, beta_range=None):
    """
    Log-spaced beta grid: `beta_range` if given, else the pipeline's
    BETA_MIN / BETA_MAX, with BETA_MIN = None resolved exactly as the
    pipeline does (Cv_AutoTune.resolve_beta_min: T_max = 10 (E1 - E0)),
    E1 - E0 taken from a small xi = 1 solve of the base potential.
    """
    if beta_range is not None:
        beta_min, beta_max = (float(b) for b in beta_range)
    else:
        beta_max = config.BETA_MAX
        beta_min = config.BETA_MIN
        if beta_min is None:
            beta_min = resolve_beta_min(None, base_spectra.get(1.0, 1.0 / beta_max))
    return np.geomspace(beta_min, beta_max, n_beta)


# =====================================================================
# One potential: quantum Cv + classical limit
# =====================================================================
def scan_potential(potential_func, beta_arr, settings, spectra=None):
    """
    The quantum and classical-limit Cv(T) of one potential at the given
    resolution, from the pipeline's own functions (module docstring).

    Returns
    -------
    dict with keys:
        cv_quantum, cv_classical, error_estimate : ndarray
            error_estimate is the xi-scan's distance-to-limit estimate
            (NaN where the classical limit did not converge).
        stop_reasons : dict  -- {stop reason: count} over temperatures
        energies : ndarray   -- the xi = 1 spectrum (that of V)
        solves, xi_top, seconds
    """
    t0 = time.time()
    spectra = spectra if spectra is not None else _new_spectra(potential_func, settings)
    solves_before = len(spectra.solves)
    sweep = sweep_temperature_range(
        potential_func, beta_arr,
        config.XI_START, settings["tol_xi"], settings["min_stable_xi"], settings["xi_mult"],
        settings["max_xi_steps"], config.TOL_CV, config.MIN_STABLE_N,
        mass=config.MASS, hbar=config.HBAR, thermal_coverage=settings["thermal_coverage"],
        xi_max=config.XI_MAX, spectra=spectra, verbose=True,
    )
    # The xi = 1 rung, solved for the hottest temperature, is the spectrum of V.
    energies = spectra.get(1.0, 1.0 / float(np.min(beta_arr)))
    cv_quantum = compute_quantum_heat_capacity_curve(energies, beta_arr, xi=1.0)

    stop_reasons = {}
    for xr in sweep["xi_results"]:
        stop_reasons[xr["stop_reason"]] = stop_reasons.get(xr["stop_reason"], 0) + 1
    new_solves = spectra.solves[solves_before:]
    return {
        "cv_quantum": cv_quantum, "cv_classical": sweep["cv_classical"],
        "error_estimate": sweep["error_estimate"], "stop_reasons": stop_reasons,
        "energies": energies, "solves": len(new_solves),
        "xi_top": max([s["xi"] for s in new_solves], default=float("nan")),
        "seconds": time.time() - t0,
    }


# =====================================================================
# Verdict: is the quantum Cv above the classical one, beyond the error?
# =====================================================================
def compare_quantum_classical(T_arr, cv_quantum, cv_classical, error_estimate):
    """
    Classify every temperature by d = Cv_q - Cv_cl against the classical
    value's error estimate eps, with m = ERROR_MARGIN: above (d > m eps),
    below (d < -m eps) or unresolved (|d| <= m eps); see the module
    docstring.

    Returns
    -------
    dict with keys:
        verdict : str  -- "ABOVE", "below" or "unresolved"
        d, eps : ndarray
        above, below, unresolved, missing : ndarray of bool
        hot_tail_only : bool
            True if every unresolved temperature lies in the contiguous
            band at the hot end, where Cv_q -> Cv_cl anyway.
        hot_tail_T : float or None  -- coldest temperature of that band
        lines : list of str  -- the console report
    """
    T_arr = np.asarray(T_arr, dtype=float)
    d = np.asarray(cv_quantum, dtype=float) - np.asarray(cv_classical, dtype=float)
    eps = np.asarray(error_estimate, dtype=float)
    ok = np.isfinite(d) & np.isfinite(eps)
    with np.errstate(invalid="ignore"):          # NaN where the classical limit is missing
        above = ok & (d > ERROR_MARGIN * eps)
        below = ok & (d < -ERROR_MARGIN * eps)
    unresolved = ok & ~above & ~below
    missing = ~ok

    # Unresolved band contiguous with the hottest temperature.
    tail = 0
    for i in np.argsort(T_arr)[::-1]:
        if not unresolved[i]:
            break
        tail += 1
    hot_tail_only = tail > 0 and tail == int(unresolved.sum())
    hot_tail_T = float(np.sort(T_arr)[::-1][tail - 1]) if tail else None

    def _t_range(mask):
        t = T_arr[mask]
        return f"T {t.min():.3g}-{t.max():.3g}"

    lines = []
    if above.any():
        verdict = "ABOVE"
        i = int(np.nanargmax(np.where(above, d / eps, np.nan)))
        lines.append(f"quantum ABOVE classical at {int(above.sum())} temperatures ({_t_range(above)}); "
                     f"clearest excess {d[i]:+.4f} (error estimate {eps[i]:.1e}) at T = {T_arr[i]:.3g}")
    elif below.any() and (not unresolved.any() or hot_tail_only):
        verdict = "below"
        i = int(np.nanargmax(np.where(below, d, np.nan)))
        lines.append(f"quantum below classical at every resolved temperature "
                     f"(closest: {d[i]:+.4f}, error estimate {eps[i]:.1e}, at T = {T_arr[i]:.3g})")
    elif not ok.any():
        verdict = "unresolved"
        lines.append("no temperature with a classical value to compare against")
    else:
        verdict = "unresolved"
        lines.append("no resolved excess, but the comparison is not resolved everywhere (below)")
    if unresolved.any():
        if hot_tail_only:
            lines.append(f"unresolved (|d| <= {ERROR_MARGIN:g} x error) only in the hot tail T >= {hot_tail_T:.3g}, "
                         f"where the quantum Cv approaches the classical one anyway")
        else:
            lines.append(f"unresolved (|d| <= {ERROR_MARGIN:g} x error) at {int(unresolved.sum())} temperatures "
                         f"({_t_range(unresolved)}), not only in the hot tail -- raise "
                         f"QUICK_SCAN_RESOLUTION, or zoom QUICK_SCAN_BETA_RANGE onto that window")
    if missing.any():
        lines.append(f"classical limit missing at {int(missing.sum())} temperatures ({_t_range(missing)}) "
                     f"-- no xi plateau there (stop reasons above)")
    return {"verdict": verdict, "d": d, "eps": eps, "above": above, "below": below,
            "unresolved": unresolved, "missing": missing, "hot_tail_only": hot_tail_only,
            "hot_tail_T": hot_tail_T, "lines": lines}


# =====================================================================
# Orchestrator
# =====================================================================
def run_quick_scan(mode=None, resolution=None, beta_range=None):
    """
    Quick scan of config's potential ("single") or of config.SCAN_PARAM's
    sweep ("sweep"): quantum and classical-limit Cv(T) of every
    potential at the chosen resolution, one figure, a console verdict.

    Parameters
    ----------
    mode : "single", "sweep" or None (None = config.QUICK_SCAN_MODE)
    resolution : int or None (None = config.QUICK_SCAN_RESOLUTION)
    beta_range : (float, float) or None
        Zoom window (beta_min, beta_max); None = config.QUICK_SCAN_BETA_RANGE,
        and if that is None too, the pipeline's BETA_MIN / BETA_MAX.

    Returns
    -------
    dict with keys:
        T_arr, beta_arr : ndarray
        settings : dict  -- the resolution's preset (+ max_xi_steps)
        results : list of dict  -- one per potential that ran: its
            `scan_potential` output plus "params", "value" and "verdict"
            (the `compare_quantum_classical` output)
        failed : list of dict  -- {"params", "value", "error", "diagnosis"}
        figure_path : str
        seconds : float
    """
    t_start = time.time()
    mode = mode or config.QUICK_SCAN_MODE
    resolution = resolution if resolution is not None else config.QUICK_SCAN_RESOLUTION
    beta_range = beta_range if beta_range is not None else config.QUICK_SCAN_BETA_RANGE
    if mode not in ("single", "sweep"):
        raise ValueError(f"QUICK_SCAN_MODE must be 'single' or 'sweep', not {mode!r}.")
    settings = quick_settings(resolution)

    base_params = config.POTENTIAL_PARAMS
    set_figure_context(config.SYSTEM_NAME, base_params)
    if mode == "sweep":
        label_param = config.SCAN_PARAM
        variants = generate_variant_params(base_params, config.SCAN_PARAM, config.SCAN_STEP,
                                           config.SCAN_COUNT, config.SCAN_SYMMETRIC_VALUE)
    else:
        # The legend names the potential by one of its parameters.
        label_param = config.SCAN_PARAM if config.SCAN_PARAM in base_params else next(iter(base_params))
        variants = [dict(base_params)]

    base_func = functools.partial(config.my_potential, p=base_params)
    base_spectra = _new_spectra(base_func, settings)
    beta_arr = temperature_grid(base_spectra, settings["n_beta"], beta_range)
    T_arr = 1.0 / beta_arr

    rule = "=" * 60
    print(f"\n{rule}\n  QUICK SCAN ({mode}, resolution {resolution})  {config.SYSTEM_NAME}\n{rule}")
    print(f"  {len(beta_arr)} temperatures, T = {T_arr.min():.4g} .. {T_arr.max():.4g}; "
          f"tol_xi = {settings['tol_xi']:g}, xi ladder x{settings['xi_mult']:g} from "
          f"{config.XI_START:g} (cap {config.XI_MAX:g}), plateau {settings['min_stable_xi']} steps, "
          f"coverage {settings['thermal_coverage']:g} k_B T, DVR tolerance {settings['dvr_tolerance']:g}")
    if mode == "sweep":
        print(f"  {label_param} over [{', '.join(f'{p[label_param]:g}' for p in variants)}]")

    results, failed = [], []
    for i, params in enumerate(variants):
        value = params[label_param]
        print(f"\n  [{i + 1}/{len(variants)}] {label_param} = {value:g}", flush=True)
        is_base = all(np.isclose(params[k], base_params[k]) for k in base_params)
        potential_func = base_func if is_base else functools.partial(config.my_potential, p=params)
        try:
            r = scan_potential(potential_func, beta_arr, settings, base_spectra if is_base else None)
        except Exception as exc:
            diagnosis = _diagnose_variant_failure(potential_func)
            print(f"  ✗ {label_param} = {value:g} failed: {exc}")
            if diagnosis:
                print(f"    diagnosis: {diagnosis}")
            failed.append({"params": params, "value": value, "error": str(exc), "diagnosis": diagnosis})
            continue
        r["params"], r["value"] = params, value
        r["verdict"] = compare_quantum_classical(T_arr, r["cv_quantum"], r["cv_classical"],
                                                 r["error_estimate"])
        print(f"    {r['solves']} DVR solves (xi up to {r['xi_top']:.4g}), {r['seconds']:.0f}s; "
              f"xi-scan outcomes {r['stop_reasons']}")
        for line in r["verdict"]["lines"]:
            print(f"    -> {line}")
        results.append(r)

    if mode == "sweep":
        title = (f"{config.SYSTEM_NAME} — quick scan (resolution {resolution}): "
                 f"quantum vs classical $C_v(T)$ across {label_param}")
        formula_text = format_potential_formula(config.POTENTIAL_FORMULA, base_params, label_param)
    else:
        title = f"{config.SYSTEM_NAME} — quick scan (resolution {resolution}): quantum vs classical $C_v(T)$"
        formula_text = format_potential_formula(config.POTENTIAL_FORMULA, base_params, None)
    name = f"quick_{mode}_res{resolution}"
    if beta_range is not None:
        name += f"_beta{beta_range[0]:g}_to_{beta_range[1]:g}"
    figure_path = plot_coefficient_sweep(
        T_arr, [r["value"] for r in results], [r["cv_quantum"] for r in results],
        [r["cv_classical"] for r in results], label_param, config.SYSTEM_NAME, formula_text,
        T_units_label=config.T_UNITS_LABEL,
        variant_deltas=[float(r["energies"][1] - r["energies"][0]) for r in results],
        symmetric_value=config.SCAN_SYMMETRIC_VALUE if label_param == config.SCAN_PARAM else None,
        classical_incomplete={r["value"] for r in results if np.isnan(r["cv_classical"]).any()},
        title=title, category="quick_scan", name=name,
    )

    seconds = time.time() - t_start
    print(f"\n{rule}\n  Quick scan summary (resolution {resolution}, {seconds:.0f}s total; "
          f"resolved = |Cv_q - Cv_cl| > {ERROR_MARGIN:g} x the classical error estimate)")
    print(f"  {label_param:>8}  {'verdict':<10}  {'max(Cv_q - Cv_cl)':>17}  {'error est.':>10}  {'at T':>7}  "
          f"{'unresolved':>10}  {'cl. missing':>11}")
    for r in results:
        v = r["verdict"]
        d = np.where(np.isfinite(v["eps"]), v["d"], np.nan)
        if np.isfinite(d).any():
            j = int(np.nanargmax(d))
            d_text, e_text, t_text = f"{d[j]:+.4f}", f"{v['eps'][j]:.1e}", f"{T_arr[j]:.3g}"
        else:
            d_text = e_text = t_text = "-"
        unres = int(v["unresolved"].sum())
        unres_text = f"{unres} (hot tail)" if v["hot_tail_only"] else str(unres)
        print(f"  {r['value']:>8g}  {v['verdict']:<10}  {d_text:>17}  {e_text:>10}  {t_text:>7}  "
              f"{unres_text:>10}  {int(v['missing'].sum()):>11}")
    for f in failed:
        print(f"  {f['value']:>8g}  failed: {f['error'].splitlines()[0]}")
    print(f"  Figure: {figure_path}\n{rule}\n")

    return {"T_arr": T_arr, "beta_arr": beta_arr, "settings": settings, "results": results,
            "failed": failed, "figure_path": figure_path, "seconds": seconds}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Quick, low-resolution quantum vs classical Cv(T) scan "
                    "(Quick_Scan.py's docstring; settings in config.py).")
    parser.add_argument("--mode", choices=["single", "sweep"],
                        help="single potential or coefficient sweep (default: config.QUICK_SCAN_MODE)")
    parser.add_argument("--resolution", type=int, choices=sorted(config.QUICK_SCAN_PRESETS),
                        help="key of config.QUICK_SCAN_PRESETS (default: config.QUICK_SCAN_RESOLUTION)")
    parser.add_argument("--beta-range", type=float, nargs=2, metavar=("BETA_MIN", "BETA_MAX"),
                        help="zoom window (default: config.QUICK_SCAN_BETA_RANGE, else BETA_MIN/BETA_MAX)")
    args = parser.parse_args()
    run_quick_scan(args.mode, args.resolution, args.beta_range)
