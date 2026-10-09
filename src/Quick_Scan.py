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

FASTER WITHOUT TOUCHING THE ENGINE
---------------------------------------------------------------------
The engine sizes every DVR solve with a WKB level count whose bisection
re-samples V on the same 40,001-point windows dozens of times; at
resolution 1 that was about half the run time. The quick scan hands the
engine `RememberedPotential(V)`, which returns stored copies for inputs
it has already seen. The results are bit-identical and the
classical-limit engine (Classical_Limit_Numerical.py) is unchanged.

WHERE QUANTUM AND CLASSICAL MERGE: T_merge
---------------------------------------------------------------------
At high T the quantum Cv approaches the classical one: the leading
quantum correction (Wigner-Kirkwood) is of order (hbar omega / k_B T)^2
and falls off as T^-p with p >= 1 (FINDINGS.md, "Why quantum and
classical Cv merge at high T"). Computing there is both unnecessary and
expensive: every temperature needs all its thermally accessible levels,
whose number grows with T, and a dense DVR solve costs ~grid^3.

So, when the window reaches above T0 = 10 max(E1 - E0, hbar omega_min)
(the pipeline's own "T -> infinity" checkpoint, 10 (E1 - E0), guarded
against tunnelling doublets by the vibrational quantum at the minimum),
`merge_temperature` probes d = Cv_q - Cv_cl at T0, 2 T0, 4 T0, ... with
this resolution's own settings. T_merge is the first T of two
consecutive probes that both have |d| + 2 eps <= tol_xi, and whose d
falls at least as fast as 1/T (or is below the error), i.e. is in the
asymptotic regime. Temperatures above T_merge are not computed; the
console says that quantum = classical within tol_xi there and how many
levels they would need. A window lying wholly above T_merge gets the
verdict "merged".

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
import collections
import functools
import math
import time

import numpy as np

import config
from Classical_Limit_Numerical import ScaledSpectra, sweep_temperature_range, _wkb_level_count
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
# A potential that remembers its recent evaluations (speed only)
# =====================================================================
class RememberedPotential:
    """
    V(x) that remembers its last `size` array evaluations, keyed by the
    exact input array.

    The engine sizes every DVR solve with a WKB level count
    (Classical_Limit_Numerical._wkb_energy_for_count), whose bisection
    re-samples V on the same 40,001-point windows dozens of times per
    solve; at resolution 1 that was about half of the run time. Handing
    the engine this wrapper instead of V removes the repeated
    evaluations without touching the engine. A repeated call returns a
    copy of exactly what V returned for the identical input, so every
    number is unchanged.
    """

    def __init__(self, func, size=64):
        self.func = func
        self.size = size
        self._store = collections.OrderedDict()

    def __call__(self, x):
        if np.ndim(x) == 0:
            return self.func(x)
        x = np.asarray(x)
        key = (x.dtype.str, x.shape, x.tobytes())
        hit = self._store.get(key)
        if hit is not None:
            self._store.move_to_end(key)
            return hit.copy()
        v = self.func(x)
        self._store[key] = np.array(v, copy=True)
        if len(self._store) > self.size:
            self._store.popitem(last=False)
        return v


# =====================================================================
# Where quantum and classical merge: the high-temperature limit T_merge
# =====================================================================
# Most temperature doublings tried before giving up on finding T_merge.
MAX_MERGE_PROBES = 20


def estimated_levels(spectra, T):
    """
    Levels the quantum Cv at temperature T needs (the spectrum of V
    itself, xi = 1), from the engine's own sizing (ScaledSpectra._solve:
    1.15 x the WKB level count up to E_0 + 1.2 x thermal_coverage x
    k_B T, plus 10). The xi-scan needs ~xi times more at a rung xi.
    """
    return 1.15 * _wkb_level_count(spectra.potential_func,
                                   spectra.v_min + 1.2 * spectra.thermal_coverage * T,
                                   spectra.v_min, spectra.mass, spectra.hbar) + 10


def _hbar_omega_at_minimum(spectra):
    """hbar sqrt(V''(x_min) / m) at the minimum the engine located (central
    difference); 0 if the curvature is not positive (e.g. a pure quartic)."""
    x = spectra.x_at_min
    h = 1e-3 * (1.0 + abs(x))
    v = np.asarray(spectra.potential_func(np.array([x - h, x, x + h])), dtype=float)
    curvature = (v[0] - 2.0 * v[1] + v[2]) / h**2
    return spectra.hbar * math.sqrt(curvature / spectra.mass) if curvature > 0 else 0.0


def _difference_at(potential_func, spectra, T, settings):
    """Cv_q - Cv_cl and the classical value's error estimate at one
    temperature, from the same xi-scan and settings as the scan itself."""
    beta = np.array([1.0 / T])
    sweep = sweep_temperature_range(
        potential_func, beta,
        config.XI_START, settings["tol_xi"], settings["min_stable_xi"], settings["xi_mult"],
        settings["max_xi_steps"], config.TOL_CV, config.MIN_STABLE_N,
        mass=config.MASS, hbar=config.HBAR, thermal_coverage=settings["thermal_coverage"],
        xi_max=config.XI_MAX, spectra=spectra, verbose=False,
    )
    cv_quantum = compute_quantum_heat_capacity_curve(spectra.get(1.0, T), beta, xi=1.0)[0]
    return float(cv_quantum - sweep["cv_classical"][0]), float(sweep["error_estimate"][0])


def merge_temperature(potential_func, spectra, settings, T_hot):
    """
    The temperature T_merge above which the quantum and classical Cv of
    this potential agree within the resolution's tol_xi, shown
    numerically (module docstring; the physics is in FINDINGS.md, "Why
    quantum and classical Cv merge at high T").

    d = Cv_q - Cv_cl is probed at T0, 2 T0, 4 T0, ..., with
    T0 = 10 max(E1 - E0, hbar omega_min): the pipeline's own "T -> inf"
    checkpoint 10 (E1 - E0), guarded against a tunnelling doublet (tiny
    E1 - E0) by the vibrational quantum at the minimum. T_merge is the
    colder of the first two consecutive probes that both have
    |d| + ERROR_MARGIN eps <= tol_xi and between which |d| falls at least
    as fast as 1/T (or the hotter |d| is within its error): the
    asymptotic regime, where |d| ~ T^-p with p >= 1 keeps falling.
    Probing stops once past T_hot.

    Returns
    -------
    dict with keys:
        status : str
            "not needed" (T_hot <= T0), "found", "not found" or "failed".
        T_merge : float or None
        T0 : float
        probes : list of (T, d, eps)
        slope : float or None
            d ln|d| / d ln T between the two accepting probes (-inf if
            the colder |d| was 0).
        error : str or None   -- why probing failed
    """
    energies = spectra.get(1.0, 1e-12)        # any xi = 1 spectrum: only E1 - E0 is used
    T0 = 10.0 * max(float(energies[1] - energies[0]), _hbar_omega_at_minimum(spectra))
    out = {"status": "not needed", "T_merge": None, "T0": T0, "probes": [],
           "slope": None, "error": None}
    if T_hot <= T0:
        return out
    tol = settings["tol_xi"]
    T = T0
    try:
        for _ in range(MAX_MERGE_PROBES):
            out["probes"].append((T,) + _difference_at(potential_func, spectra, T, settings))
            if len(out["probes"]) >= 2:
                (Ta, da, ea), (Tb, db, eb) = out["probes"][-2:]
                within = (np.all(np.isfinite([da, ea, db, eb]))
                          and abs(da) + ERROR_MARGIN * ea <= tol
                          and abs(db) + ERROR_MARGIN * eb <= tol)
                if within:
                    slope = (math.log(abs(db) / abs(da)) / math.log(Tb / Ta)
                             if da != 0 and db != 0 else -math.inf)
                    if abs(db) <= ERROR_MARGIN * eb or slope <= -1.0:
                        out.update(status="found", T_merge=Ta, slope=slope)
                        return out
            if T > T_hot:
                break
            T *= 2.0
    except Exception as exc:            # e.g. MemoryError: too many levels to probe
        out.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        return out
    out["status"] = "not found"
    return out


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

    If the window reaches above this potential's merge temperature
    (`merge_temperature`), the temperatures above T_merge are not
    computed: both curves are NaN there, "skipped" marks them, and a
    note says why (quantum = classical within tol_xi, and how many levels
    they would need).

    Returns
    -------
    dict with keys:
        cv_quantum, cv_classical, error_estimate : ndarray
            error_estimate is the xi-scan's distance-to-limit estimate
            (NaN where the classical limit did not converge).
        skipped : ndarray of bool  -- temperatures above T_merge
        merge : dict               -- the `merge_temperature` output
        notes : list of str        -- console lines about the merge limit
        stop_reasons : dict  -- {stop reason: count} over temperatures
        energies : ndarray   -- the xi = 1 spectrum (that of V)
        solves, xi_top, seconds
    """
    t0 = time.time()
    spectra = spectra if spectra is not None else _new_spectra(potential_func, settings)
    solves_before = len(spectra.solves)
    T_arr = 1.0 / np.asarray(beta_arr, dtype=float)
    tol = settings["tol_xi"]
    merge = merge_temperature(potential_func, spectra, settings, float(T_arr.max()))
    keep = np.ones(len(T_arr), bool)
    notes = []
    if merge["status"] == "found":
        keep = T_arr <= merge["T_merge"] * (1.0 + 1e-12)
        if not keep.all():
            (Ta, da, _), (Tb, db, _) = merge["probes"][-2:]
            hot = T_arr[~keep]
            trend = (f"falling as T^{merge['slope']:.1f}" if np.isfinite(merge["slope"])
                     else "the hotter one within its error")
            notes.append(
                f"{int((~keep).sum())} of {len(T_arr)} temperatures (T {hot.min():.3g}-{hot.max():.3g}) "
                f"lie above T_merge = {merge['T_merge']:.3g} and were not computed: there Cv_q = Cv_cl "
                f"within {tol:g}, and computing them would need ~{estimated_levels(spectra, hot.min()):,.0f}"
                f"-{estimated_levels(spectra, hot.max()):,.0f} levels. Shown at T = {Ta:.3g} and {Tb:.3g} "
                f"(Cv_q - Cv_cl = {da:+.1e}, {db:+.1e}, {trend}; FINDINGS.md: why quantum and "
                f"classical Cv merge at high T)")
    elif merge["status"] == "not found":
        notes.append(f"no merge limit up to T = {merge['probes'][-1][0]:.3g} (Cv_q and Cv_cl not yet "
                     f"within {tol:g} there): the whole window is computed")
    elif merge["status"] == "failed":
        notes.append(f"the merge limit could not be established ({merge['error']}): the whole "
                     f"window is computed")

    def _full(values):
        out = np.full(len(T_arr), np.nan)
        out[keep] = values
        return out

    stop_reasons = {}
    if keep.any():
        sweep = sweep_temperature_range(
            potential_func, beta_arr[keep],
            config.XI_START, settings["tol_xi"], settings["min_stable_xi"], settings["xi_mult"],
            settings["max_xi_steps"], config.TOL_CV, config.MIN_STABLE_N,
            mass=config.MASS, hbar=config.HBAR, thermal_coverage=settings["thermal_coverage"],
            xi_max=config.XI_MAX, spectra=spectra, verbose=True,
        )
        for xr in sweep["xi_results"]:
            stop_reasons[xr["stop_reason"]] = stop_reasons.get(xr["stop_reason"], 0) + 1
        # The xi = 1 rung, solved for the hottest temperature, is the spectrum of V.
        energies = spectra.get(1.0, float(T_arr[keep].max()))
        cv_quantum = _full(compute_quantum_heat_capacity_curve(energies, beta_arr[keep], xi=1.0))
        cv_classical, error_estimate = _full(sweep["cv_classical"]), _full(sweep["error_estimate"])
    else:
        energies = spectra.get(1.0, 1e-12)
        cv_quantum = cv_classical = error_estimate = np.full(len(T_arr), np.nan)
    new_solves = spectra.solves[solves_before:]
    return {
        "cv_quantum": cv_quantum, "cv_classical": cv_classical, "error_estimate": error_estimate,
        "skipped": ~keep, "merge": merge, "notes": notes, "stop_reasons": stop_reasons,
        "energies": energies, "solves": len(new_solves),
        "xi_top": max([s["xi"] for s in new_solves], default=float("nan")),
        "seconds": time.time() - t0,
    }


# =====================================================================
# Verdict: is the quantum Cv above the classical one, beyond the error?
# =====================================================================
def compare_quantum_classical(T_arr, cv_quantum, cv_classical, error_estimate, skipped=None):
    """
    Classify every temperature by d = Cv_q - Cv_cl against the classical
    value's error estimate eps, with m = ERROR_MARGIN: above (d > m eps),
    below (d < -m eps) or unresolved (|d| <= m eps); see the module
    docstring. Temperatures in `skipped` (above the potential's merge
    temperature, see `scan_potential`) are left out of every class,
    including "missing"; if every temperature is skipped, the verdict is
    "merged".

    Returns
    -------
    dict with keys:
        verdict : str  -- "ABOVE", "below", "unresolved" or "merged"
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
    skipped = np.zeros(len(T_arr), bool) if skipped is None else np.asarray(skipped, bool)
    missing = ~ok & ~skipped

    # Unresolved band contiguous with the hottest computed temperature.
    tail, hot_tail_T = 0, None
    for i in np.argsort(T_arr)[::-1]:
        if skipped[i]:
            continue
        if not unresolved[i]:
            break
        tail += 1
        hot_tail_T = float(T_arr[i])
    hot_tail_only = tail > 0 and tail == int(unresolved.sum())

    def _t_range(mask):
        t = T_arr[mask]
        return f"T {t.min():.3g}-{t.max():.3g}"

    lines = []
    if skipped.all():
        verdict = "merged"
        lines.append("the whole window lies above T_merge, where Cv_q = Cv_cl within this "
                     "resolution's tolerance (note above): nothing to compute")
    elif above.any():
        verdict = "ABOVE"
        i = int(np.nanargmax(np.where(above, d / eps, np.nan)))
        lines.append(f"quantum ABOVE classical at {int(above.sum())} temperatures ({_t_range(above)}); "
                     f"clearest excess {d[i]:+.2e} (error estimate {eps[i]:.1e}) at T = {T_arr[i]:.3g}")
    elif below.any() and (not unresolved.any() or hot_tail_only):
        verdict = "below"
        i = int(np.nanargmax(np.where(below, d, np.nan)))
        lines.append(f"quantum below classical at every resolved temperature "
                     f"(closest: {d[i]:+.2e}, error estimate {eps[i]:.1e}, at T = {T_arr[i]:.3g})")
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

    base_func = RememberedPotential(functools.partial(config.my_potential, p=base_params))
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
        potential_func = (base_func if is_base else
                          RememberedPotential(functools.partial(config.my_potential, p=params)))
        try:
            r = scan_potential(potential_func, beta_arr, settings, base_spectra if is_base else None)
        except Exception as exc:
            diagnosis = _diagnose_variant_failure(potential_func)
            print(f"  ✗ {label_param} = {value:g} failed: {exc}")
            if isinstance(exc, MemoryError):
                print("    these temperatures need too many levels for a DVR (the matrix does not "
                      "fit in memory); at such temperatures Cv_q = Cv_cl to O((hbar omega / k_B T)^2) "
                      "(FINDINGS.md) -- choose a colder window (QUICK_SCAN_BETA_RANGE)")
            if diagnosis:
                print(f"    diagnosis: {diagnosis}")
            failed.append({"params": params, "value": value, "error": str(exc), "diagnosis": diagnosis})
            continue
        r["params"], r["value"] = params, value
        r["verdict"] = compare_quantum_classical(T_arr, r["cv_quantum"], r["cv_classical"],
                                                 r["error_estimate"], r["skipped"])
        print(f"    {r['solves']} DVR solves (xi up to {r['xi_top']:.4g}), {r['seconds']:.0f}s; "
              f"xi-scan outcomes {r['stop_reasons']}")
        for note in r["notes"]:
            print(f"    ⚠ {note}")
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
    figure_path = None
    plotted = [r for r in results if not r["skipped"].all()]     # a "merged" potential has no curve
    if plotted:
        figure_path = plot_coefficient_sweep(
            T_arr, [r["value"] for r in plotted], [r["cv_quantum"] for r in plotted],
            [r["cv_classical"] for r in plotted], label_param, config.SYSTEM_NAME, formula_text,
            T_units_label=config.T_UNITS_LABEL,
            variant_deltas=[float(r["energies"][1] - r["energies"][0]) for r in plotted],
            symmetric_value=config.SCAN_SYMMETRIC_VALUE if label_param == config.SCAN_PARAM else None,
            classical_incomplete={r["value"] for r in plotted
                                  if np.isnan(r["cv_classical"][~r["skipped"]]).any()},
            title=title, category="quick_scan", name=name,
        )

    seconds = time.time() - t_start
    print(f"\n{rule}\n  Quick scan summary (resolution {resolution}, {seconds:.0f}s total; "
          f"resolved = |Cv_q - Cv_cl| > {ERROR_MARGIN:g} x the classical error estimate)")
    print(f"  {label_param:>8}  {'verdict':<10}  {'max(Cv_q - Cv_cl)':>17}  {'error est.':>10}  {'at T':>7}  "
          f"{'unresolved':>10}  {'cl. missing':>11}  {'T_merge':>8}  {'above it':>8}")
    for r in results:
        v = r["verdict"]
        d = np.where(np.isfinite(v["eps"]), v["d"], np.nan)
        if np.isfinite(d).any():
            j = int(np.nanargmax(d))
            d_text, e_text, t_text = f"{d[j]:+.2e}", f"{v['eps'][j]:.1e}", f"{T_arr[j]:.3g}"
        else:
            d_text = e_text = t_text = "-"
        unres = int(v["unresolved"].sum())
        unres_text = f"{unres} (hot tail)" if v["hot_tail_only"] else str(unres)
        m = r["merge"]
        merge_text = {"found": f"{m['T_merge']:.3g}" if m["T_merge"] else "-",
                      "not needed": "-", "not found": "none", "failed": "failed"}[m["status"]]
        print(f"  {r['value']:>8g}  {v['verdict']:<10}  {d_text:>17}  {e_text:>10}  {t_text:>7}  "
              f"{unres_text:>10}  {int(v['missing'].sum()):>11}  {merge_text:>8}  "
              f"{int(r['skipped'].sum()):>8}")
    for f in failed:
        reason = f["error"].splitlines()[0] if f["error"] else "error"
        print(f"  {f['value']:>8g}  failed: {reason if len(reason) <= 80 else reason[:77] + '...'}")
    print(f"  Figure: {figure_path or 'none (nothing was computed)'}\n{rule}\n")

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
