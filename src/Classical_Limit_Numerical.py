"""
Classical_Limit_Numerical.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
General-purpose (system-agnostic) numerics for the CLASSICAL LIMIT of
the heat capacity Cv(T) of any smooth 1-D potential V(x), using the
xi-scaling of Gelbwaser-Klimovsky et al., "Single-atom heat machines
enabled by energy quantization" (docs/refs; main text, SI-III item 2,
SI-IV): scale the potential AND the temperature by the same factor,

    V -> xi^2 V,    T -> xi^2 T,

and let xi grow. Their Eq. S7,

    E_n(hbar, xi^2 V) = xi^2 E_n(hbar/xi, V),

shows this is exactly hbar -> hbar/xi at fixed V and T, so xi -> inf
is the classical limit (hbar -> 0) of the SAME system at the SAME
temperature. Cv is scanned on a geometric xi ladder and its plateau is
the classical-limit Cv at that temperature.

HOW ONE XI STEP IS EVALUATED
---------------------------------------------------------------------
For each xi on the ladder the DVR is RE-SOLVED for the scaled
potential xi^2 V (DVR_Algorithm, with the same 3-pass convergence check
as the base run; see `ScaledSpectra`), and then

    Cv(T; xi) = compute_cv(E_n(xi^2 V), beta, xi)

-- compute_cv's weights exp(-beta E / xi^2) are the T -> xi^2 T half of
the prescription, the re-solved spectrum is the V -> xi^2 V half. The
spectrum of xi^2 V does not depend on T, so ONE solve per xi serves
every temperature (cached in `ScaledSpectra`).

WHY THE SPECTRUM MUST BE RE-SOLVED (branch correct-classical)
---------------------------------------------------------------------
Until this correction the scan reused the spectrum of the UNSCALED V
at every xi, i.e. applied only T -> xi^2 T: Cv_xi(T) was the plain
quantum Cv at the hotter temperature xi^2 T, whose plateau is the
system's high-temperature Cv at EVERY T. That coincides with the
classical limit only when every level gap scales by one common factor
under V -> xi^2 V (harmonic oscillator, box, |x|^k) -- the only systems
the method had been checked on. For the double well it returned a
nearly flat ~0.73 where the classical Cv runs from 0.71 to 1.45. Full
audit: audit/classical_limit/AUDIT.md.

CONVERGENCE (the plateau) AND ITS ERROR ESTIMATE
---------------------------------------------------------------------
Quantum corrections to the classical Cv start at order hbar^2 (the
Wigner-Kirkwood expansion), i.e. Cv(xi) ~ Cv_cl + a / xi^2 for large
xi. On a ladder xi_k = xi_start * xi_mult^k that gives, for the
distance still left to the limit,

    eps_k = |Cv(xi_k) - Cv(xi_{k-1})| / (xi_mult^2 - 1)

so `tol_xi` is a tolerance on the CLASSICAL VALUE itself, not on the
per-step change. A step counts as stable when eps_k < tol_xi AND
Cv >= 1/2 - tol_xi; convergence needs `min_stable` consecutive stable
steps, and the reported value is the LAST point of that verified streak
(the most classical one computed). The 1/2 guard is a rigorous bound:
for H = p^2/2m + V(x) the classical Cv is 1/2 + beta^2 Var(V) >= 1/2,
so a "plateau" below it can only be the frozen-out quantum regime
(Cv ~ 0 for several steps at very low T), never the classical limit.

If the ladder ends (max_xi_steps, or the hard cap xi_max) or a DVR solve
fails before convergence, the classical value at that temperature is
NaN -- never a stand-in such as the quantum Cv.

NO FINITE-N COLLAPSE, BY CONSTRUCTION
---------------------------------------------------------------------
Each xi's spectrum is solved with enough levels to cover the
temperature it is used at, (E_max - E_0) / xi^2 >= thermal_coverage *
k_B T, the same safety margin (HOT_STATE_SAFETY) the base run uses for
its quantum Cv. Temperatures are processed hot -> cold, so a spectrum
solved for one temperature already covers every colder one. The
companion n-scan (`converge_n`) re-checks this on the converged xi's
spectrum.

This file contains NO plotting. Plotting and system-specific wiring
live in the files that import it (Quantum_Classical_Combined.py,
Cv_Numerical_Benchmark.py). `sweep_temperature_range`'s per-temperature
trace also feeds Cv_AutoTune.diagnose_escalation.
=====================================================================
"""

import contextlib
import io
import time

import numpy as np
from tqdm import tqdm

from DVR.DVR_Algorithm import (auto_configure_dvr, get_fully_converged_energy_levels,
                               colbert_miller_dvr_1d)
from DVR.DVR_Reference_Generator import compute_reference_grid_params


# Fewest levels any scaled spectrum is solved with (the n-scan needs a
# stretch of stable steps to work with).
MIN_SCALED_LEVELS = 30
# Attempts at growing a scaled spectrum's level count before giving up.
MAX_LEVEL_ATTEMPTS = 6


# =====================================================================
# Core Cv formula (also used directly for the "real" quantum Cv(T))
# =====================================================================
def compute_cv(energies, beta, xi=1.0):
    """
    Compute Cv/k_B at a single inverse temperature `beta` from a
    discrete energy spectrum, at the xi-scaled temperature xi^2 T.

    With xi=1.0 this is the literal, physical quantum heat capacity
    for the given (finite) spectrum. With xi != 1.0 it evaluates the
    spectrum at the temperature xi^2 T -- the T -> xi^2 T half of the
    classical-limit scaling. It is the classical-limit building block
    ONLY when `energies` is the spectrum of the scaled potential xi^2 V
    (see the module docstring and `converge_xi`).

    Parameters
    ----------
    energies : array_like
        Energy eigenvalues (any order; need not be pre-sorted).
    beta : float
        Inverse temperature, 1 / (k_B * T).
    xi : float, optional
        Scaling factor (default 1.0).

    Returns
    -------
    Cv : float
        Heat capacity (in units of k_B) at this beta and xi. Returns
        np.nan if the partition function underflows to zero or is
        otherwise non-finite (can happen at extreme beta*xi combos).
    """
    energies = np.asarray(energies, dtype=float)
    a = beta * energies / xi**2
    # Subtract the minimum exponent before exponentiating to avoid
    # overflow/underflow in the Boltzmann weights.
    w = np.exp(-(a - a.min()))
    Z = w.sum()
    if Z == 0 or not np.isfinite(Z):
        return np.nan
    avg_E = np.dot(w, energies) / Z
    avg_E2 = np.dot(w, energies**2) / Z
    return float((beta**2 / xi**4) * (avg_E2 - avg_E**2))


# =====================================================================
# Potential geometry + semiclassical level counting (sizes each solve)
# =====================================================================
def _allowed_region(potential_func, height):
    """
    Sampled window [x_lo, x_hi] holding every x with V(x) - V_min <=
    `height`, plus V_min and where it sits. Grows a symmetric window
    about the origin until V exceeds V_min + height at both ends.

    Raises
    ------
    RuntimeError
        If V never rises that far within |x| < 1e4 (non-confining).
    """
    L = 1.0
    while True:
        x = np.linspace(-L, L, 40001)
        v = np.asarray(potential_func(x), dtype=float)
        if min(v[0], v[-1]) - v.min() > height:
            break
        L *= 1.5
        if L > 1e4:
            raise RuntimeError(
                "Classical limit: V does not rise above V_min + "
                f"{height:.3g} within |x| < 1e4 -- the potential is not confining.")
    keep = np.nonzero(v - v.min() <= height)[0]
    i_min = int(np.argmin(v))
    return (float(x[max(keep[0] - 1, 0)]), float(x[min(keep[-1] + 1, len(x) - 1)]),
            float(v[i_min]), float(x[i_min]))


def _wkb_level_count(potential_func, energy, v_min, mass, hbar):
    """Semiclassical number of levels below `energy`:
    (1/pi hbar) * \\int sqrt(2 m (energy - V)) dx over the allowed region."""
    x_lo, x_hi, _, _ = _allowed_region(potential_func, energy - v_min)
    x = np.linspace(x_lo, x_hi, 20001)
    p = np.sqrt(np.clip(2.0 * mass * (energy - np.asarray(potential_func(x), dtype=float)), 0.0, None))
    return float(p.sum() * (x[1] - x[0]) / (np.pi * hbar))


def _wkb_energy_for_count(potential_func, count, v_min, mass, hbar):
    """Inverse of `_wkb_level_count`: the energy below which ~`count`
    levels lie (bisection; the count is monotonic in energy)."""
    lo, hi = v_min, v_min + 1.0
    while _wkb_level_count(potential_func, hi, v_min, mass, hbar) < count:
        hi = v_min + 2.0 * (hi - v_min)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _wkb_level_count(potential_func, mid, v_min, mass, hbar) < count:
            lo = mid
        else:
            hi = mid
    return hi


# =====================================================================
# The spectra of the scaled potentials xi^2 V, solved once and cached
# =====================================================================
class ScaledSpectra:
    """
    The DVR spectra E_n(hbar, xi^2 V) used by the xi-scan, one per xi,
    solved on demand and cached. A spectrum solved for temperature T
    holds every level up to E_0 + thermal_coverage * k_B T (in units of
    the unscaled V, i.e. after dividing by xi^2), so it also serves every
    colder temperature -- which is why `sweep_temperature_range` walks
    temperatures hot -> cold. Pass the same object to several sweeps
    (e.g. successive auto-tune rounds) to reuse its solves.

    Each solve uses the pipeline's own DVR: `auto_configure_dvr` for the
    grid (given an energy ceiling from WKB, since xi^2 V has level
    spacings ~xi times those of V) and `get_fully_converged_energy_levels`
    for the 3-pass resolution/span check, with both tolerances scaled by
    xi^2 so they mean the same thing (1e-5 in units of V) at every xi.
    With span_factor/dx_factor != 1 (Section 6's numerical reference),
    the auto-configured grid is widened/refined exactly as
    DVR_Reference_Generator does for the base spectrum and solved once.

    Parameters
    ----------
    potential_func : callable
        The UNSCALED V(x) (smooth, confining).
    mass, hbar : float, optional
        Physical constants of the system (default 1.0).
    thermal_coverage : float, optional
        Levels are kept up to E_0 + thermal_coverage * k_B T (default
        20.0, matching config.HOT_STATE_SAFETY).
    span_factor, dx_factor : float, optional
        Grid widening / refinement relative to the auto-configured grid
        (default 1.0, 1.0 = the base grid).
    dvr_tolerance : float, optional
        Convergence tolerance of every DVR solve, in units of V
        (default 1e-5, the base run's own tolerance).
    verbose : bool, optional
        Print the DVR's own per-solve console output (default False).
    """

    def __init__(self, potential_func, mass=1.0, hbar=1.0, thermal_coverage=20.0,
                 span_factor=1.0, dx_factor=1.0, dvr_tolerance=1e-5, verbose=False):
        self.potential_func = potential_func
        self.mass = mass
        self.hbar = hbar
        self.thermal_coverage = thermal_coverage
        self.span_factor = span_factor
        self.dx_factor = dx_factor
        self.dvr_tolerance = dvr_tolerance
        self.verbose = verbose
        self._cache = {}
        self.solves = []        # one record per DVR solve, for diagnostics
        _, _, self.v_min, self.x_at_min = _allowed_region(potential_func, 1.0)

    def get(self, xi, T):
        """Spectrum of xi^2 V covering temperature T (solved if needed)."""
        rec = self._cache.get(xi)
        if rec is None or rec["T_covered"] < T:
            rec = self._solve(xi, T)
            self._cache[xi] = rec
        return rec["energies"]

    def _solve(self, xi, T):
        t0 = time.time()
        need = self.thermal_coverage * T          # energy span to cover, units of V
        hbar_eff = self.hbar / xi                 # Eq. S7: same level count as V at hbar/xi
        n_levels = max(MIN_SCALED_LEVELS, int(np.ceil(
            1.15 * _wkb_level_count(self.potential_func, self.v_min + 1.2 * need,
                                    self.v_min, self.mass, hbar_eff))) + 10)
        V_scaled = lambda x: xi**2 * self.potential_func(x)
        tol = self.dvr_tolerance * xi**2          # the scaled spectrum is xi^2 larger
        for _ in range(MAX_LEVEL_ATTEMPTS):
            # Grid must resolve the top requested level: its energy from WKB, with margin.
            e_top = _wkb_energy_for_count(self.potential_func, 1.1 * n_levels,
                                          self.v_min, self.mass, hbar_eff)
            # The allowed window at e_top is the same for V and xi^2 V. Start
            # from it with a padding of a quarter of its width (never more
            # than the DVR's default 2.0): at large xi the tails are short and
            # the default padding would dominate the grid. Stage 2 of
            # auto_configure_dvr still widens the span until it converges.
            x_left, x_right, _, _ = _allowed_region(self.potential_func, e_top - self.v_min)
            padding = min(2.0, 0.25 * (x_right - x_left))
            with _maybe_quiet(not self.verbose):
                x_min, x_max, n_grid = auto_configure_dvr(
                    V_scaled, n_levels, mass=self.mass, hbar=self.hbar, x0_guess=self.x_at_min,
                    span_tol=tol, energy_ceiling=xi**2 * e_top,
                    turning_points=(x_left, x_right), padding=padding)
                if self.span_factor == 1.0 and self.dx_factor == 1.0:
                    energies = get_fully_converged_energy_levels(
                        V_scaled, n_levels, x_min, x_max, n_grid,
                        mass=self.mass, hbar=self.hbar, tolerance=tol)
                else:
                    g = compute_reference_grid_params(x_min, x_max, n_grid,
                                                      self.span_factor, self.dx_factor)
                    n_grid = g["num_points_ref"]
                    energies = colbert_miller_dvr_1d(
                        V_scaled, n_levels, g["x_min_ref"], g["x_max_ref"], n_grid,
                        self.mass, self.hbar)
            if (energies[-1] - energies[0]) / xi**2 >= need:
                break
            n_levels = int(1.3 * n_levels) + 1
        else:
            raise RuntimeError(
                f"Classical limit: could not cover k_B T = {T:.4g} at xi = {xi:.4g} "
                f"after {MAX_LEVEL_ATTEMPTS} level-count increases (last: {n_levels}).")
        self.solves.append({"xi": xi, "T_covered": T, "num_levels": n_levels,
                            "grid_points": int(n_grid), "seconds": time.time() - t0})
        return {"energies": energies, "T_covered": T}


@contextlib.contextmanager
def _maybe_quiet(quiet):
    """Silence the DVR's per-solve console output when `quiet`."""
    if quiet:
        with contextlib.redirect_stdout(io.StringIO()):
            yield
    else:
        yield


# =====================================================================
# xi-convergence at one temperature: walk the ladder to the plateau
# =====================================================================
def converge_xi(spectra, beta, xi_ladder, tol_xi, min_stable, xi_mult, ladder_capped=False):
    """
    Walk the xi ladder at one temperature, evaluating
    Cv(T; xi) = compute_cv(E_n(xi^2 V), beta, xi) with the spectrum of
    the SCALED potential at every step (from `spectra`), until the
    plateau -- the classical-limit Cv -- is reached.

    A step is stable when its estimated distance to the xi -> inf limit,
    eps = |Cv_k - Cv_{k-1}| / (xi_mult^2 - 1) (hbar^2 law, see module
    docstring), is below `tol_xi` AND Cv_k >= 1/2 - tol_xi (the classical
    Cv can never be below 1/2). `min_stable` consecutive stable steps
    are required; the answer is the last point of that streak.

    Parameters
    ----------
    spectra : ScaledSpectra
        Source of the scaled spectra (shared across temperatures).
    beta : float
        Inverse temperature being probed.
    xi_ladder : ndarray
        The xi values to walk, ascending (geometric, ratio xi_mult).
    tol_xi : float
        Tolerance on the classical value, in units of k_B.
    min_stable : int
        Consecutive stable steps required.
    xi_mult : float
        Ratio between consecutive ladder values.
    ladder_capped : bool, optional
        True if `xi_ladder` was shortened by the hard cap xi_max; only
        changes the stop reason reported when the ladder runs out.

    Returns
    -------
    dict with keys:
        xi_converged, cv_converged : float or None
            The last point of the verified plateau (None if not converged).
        error_estimate : float or None
            eps at that point -- the estimated distance to the limit.
        plateau : (int, int) or None
            First and last index (into xi_values) of the verified plateau.
        converged : bool
        stop_reason : str
            "converged", "max_steps" (ladder ran out), "xi_cap" (ladder
            ran out at the hard cap xi_max) or "dvr_failed".
        xi_values, cv_values, deltas, errors : list
            Full trace of the scan (deltas = |dCv| per step, errors = eps).
    """
    T = 1.0 / beta
    xis, cvs, errs = [], [], []
    streak = 0
    stop_reason = "xi_cap" if ladder_capped else "max_steps"
    for xi in xi_ladder:
        try:
            energies = spectra.get(xi, T)
        except RuntimeError:
            stop_reason = "dvr_failed"
            break
        cv = compute_cv(energies, beta, xi)
        xis.append(float(xi))
        cvs.append(cv)
        if len(cvs) == 1:
            errs.append(None)
            continue
        err = abs(cv - cvs[-2]) / (xi_mult**2 - 1.0)
        errs.append(err)
        stable = np.isfinite(err) and err < tol_xi and cv >= 0.5 - tol_xi
        streak = streak + 1 if stable else 0
        if streak >= min_stable:
            stop_reason = "converged"
            break

    converged = stop_reason == "converged"
    deltas = [None] + [abs(cvs[i] - cvs[i - 1]) for i in range(1, len(cvs))]
    return {
        "xi_converged": xis[-1] if converged else None,
        "cv_converged": cvs[-1] if converged else None,
        "error_estimate": errs[-1] if converged else None,
        "plateau": (len(xis) - 1 - min_stable, len(xis) - 1) if converged else None,
        "converged": converged, "stop_reason": stop_reason,
        "xi_values": xis, "cv_values": cvs, "deltas": deltas, "errors": errs,
    }


# =====================================================================
# n-convergence: how many energy levels are actually needed?
# =====================================================================
def converge_n(energies, beta, xi, tol_cv, min_stable):
    """
    Sweep how many of the (ascending) energy levels are included in
    the partition sum, n = 2, 3, ..., N, tracking Cv(n) at a fixed
    (beta, xi), and find the smallest n after which Cv stays stable
    for the rest of the available spectrum.

    This answers "did we include enough states for this answer to be
    trustworthy at this temperature?" -- independent of (but used
    alongside) the xi-convergence check. For the classical limit it is
    run on the converged xi's own spectrum (that of xi^2 V).

    Parameters
    ----------
    energies : array_like
        Full available energy spectrum (ascending order expected).
    beta : float
        Inverse temperature being probed.
    xi : float
        Scaling factor to use while sweeping n (the converged xi, with
        `energies` the spectrum of xi^2 V; or 1.0 for the physical,
        unscaled quantum Cv).
    tol_cv : float
        Maximum |Cv(n) - Cv(n-1)| for two consecutive n to count as stable.
    min_stable : int
        Minimum run length of stable steps, extending all the way to
        n=N, required to declare convergence.

    Returns
    -------
    dict with keys:
        n_converged, cv_converged : int/float or None
            Smallest converged n and its Cv value (None if not converged).
        converged : bool
        n_values, cv_values, deltas : list
            Full trace, useful for diagnostic plotting.
    """
    N = len(energies)
    n_values = list(range(2, N + 1))

    energies_arr = np.asarray(energies, dtype=float)
    a = beta * energies_arr / (xi**2)
    w = np.exp(-(a - a.min()))

    # Cumulative sums let us compute Cv(n) for every n in one vectorized pass
    # instead of recomputing the partition sum from scratch each time.
    Z_n = np.cumsum(w)
    E_w_n = np.cumsum(w * energies_arr)
    E2_w_n = np.cumsum(w * (energies_arr**2))

    with np.errstate(divide='ignore', invalid='ignore'):
        avg_E = E_w_n / Z_n
        avg_E2 = E2_w_n / Z_n
        cv_array = (beta**2 / xi**4) * (avg_E2 - avg_E**2)

    cv_values = cv_array[1:].tolist()
    deltas = [None] + [abs(cv_values[i] - cv_values[i - 1]) for i in range(1, len(cv_values))]
    stable = [False] + [(d < tol_cv) for d in deltas[1:]]

    # Find runs of consecutive stable steps at least min_stable long.
    runs = []
    i = 0
    while i < len(stable):
        if stable[i]:
            j = i
            while j < len(stable) and stable[j]:
                j += 1
            if j - i >= min_stable:
                runs.append((i, j - 1))
            i = j
        else:
            i += 1

    # Convergence only counts if the stable run extends all the way to
    # the end of the available spectrum (otherwise it might just be a
    # coincidental flat stretch followed by renewed drift).
    last_idx = len(n_values) - 1
    n_converged = cv_conv = None
    for start_idx, end_idx in runs:
        if end_idx == last_idx:
            n_converged = n_values[start_idx]
            cv_conv = cv_values[start_idx]
            break

    return {
        "n_converged": n_converged, "cv_converged": cv_conv,
        "converged": n_converged is not None,
        "n_values": n_values, "cv_values": cv_values, "deltas": deltas,
    }


# =====================================================================
# Sweep both convergence checks across a full temperature range
# =====================================================================
def sweep_temperature_range(potential_func, beta_arr,
                             xi_start, tol_xi, min_stable_xi, xi_mult, max_xi_steps,
                             tol_cv, min_stable_n, mass=1.0, hbar=1.0,
                             thermal_coverage=20.0, xi_max=np.inf, spectra=None,
                             span_factor=1.0, dx_factor=1.0, verbose=True):
    """
    The classical-limit Cv(T) of `potential_func` at every beta in
    `beta_arr`: `converge_xi` on the shared xi ladder at each
    temperature, then `converge_n` on the converged xi's spectrum as a
    truncation check. Temperatures are processed hot -> cold so each
    scaled spectrum is solved once (see `ScaledSpectra`).

    Parameters
    ----------
    potential_func : callable
        The UNSCALED potential V(x) of the system.
    beta_arr : array_like
        Inverse temperatures to sweep over.
    xi_start, tol_xi, min_stable_xi, xi_mult, max_xi_steps :
        The xi ladder (xi_start * xi_mult**k, k < max_xi_steps) and the
        plateau criterion -- see `converge_xi`.
    tol_cv, min_stable_n :
        Passed through to `converge_n`.
    mass, hbar : float, optional
        Physical constants of the system (default 1.0).
    thermal_coverage : float, optional
        See `ScaledSpectra` (default 20.0 = config.HOT_STATE_SAFETY).
    xi_max : float, optional
        Hard cap on xi, independent of max_xi_steps (grid sizes grow
        ~xi, so this bounds memory/time even after auto-tune escalation).
    spectra : ScaledSpectra or None, optional
        Reuse an existing cache (must be for the same potential and grid
        factors); None builds a new one.
    span_factor, dx_factor : float, optional
        Grid factors for a new cache (Section 6's reference uses > 1).
    verbose : bool, optional
        tqdm progress bar + a one-line summary.

    Returns
    -------
    dict with keys:
        cv_classical, xi_conv, n_conv : ndarray, shape (len(beta_arr),)
            Classical-limit Cv, the xi of its plateau point, and the
            converged n (NaN where the respective scan failed). Where
            the xi-scan failed, cv_classical is NaN -- no fallback.
        n_available : ndarray
            Levels in the spectrum the n-scan ran on (NaN if not run).
        error_estimate : ndarray
            Estimated distance of each value from the xi -> inf limit.
        xi_results, n_results : list of dict (n_results entry None where
            the xi-scan failed)
        xi_fail_mask, n_fail_mask : ndarray of bool
        xi_ladder : ndarray
        spectra : ScaledSpectra (its .solves lists every DVR solve)
    """
    beta_arr = np.asarray(beta_arr, dtype=float)
    if spectra is None:
        spectra = ScaledSpectra(potential_func, mass=mass, hbar=hbar,
                                thermal_coverage=thermal_coverage,
                                span_factor=span_factor, dx_factor=dx_factor)
    full_ladder = xi_start * xi_mult ** np.arange(max_xi_steps)
    xi_ladder = full_ladder[full_ladder <= xi_max]
    ladder_capped = len(xi_ladder) < len(full_ladder)

    n_T = len(beta_arr)
    cv_classical = np.full(n_T, np.nan)
    xi_conv_arr = np.full(n_T, np.nan)
    n_conv_arr = np.full(n_T, np.nan)
    n_available = np.full(n_T, np.nan)
    error_estimate = np.full(n_T, np.nan)
    xi_results, n_results = [None] * n_T, [None] * n_T

    solves_before = len(spectra.solves)
    t0 = time.time()
    order = np.argsort(beta_arr)                 # hot (small beta) first
    it = tqdm(order, desc="  Classical limit (xi-scan, xi^2 V re-solved)", unit="T") if verbose else order
    for idx in it:
        beta = beta_arr[idx]
        xr = converge_xi(spectra, beta, xi_ladder, tol_xi, min_stable_xi, xi_mult, ladder_capped)
        xi_results[idx] = xr
        if not xr["converged"]:
            continue
        cv_classical[idx] = xr["cv_converged"]
        xi_conv_arr[idx] = xr["xi_converged"]
        error_estimate[idx] = xr["error_estimate"]
        energies = spectra.get(xr["xi_converged"], 1.0 / beta)
        nr = converge_n(energies, beta, xi=xr["xi_converged"], tol_cv=tol_cv, min_stable=min_stable_n)
        n_results[idx] = nr
        n_available[idx] = len(energies)
        if nr["converged"]:
            n_conv_arr[idx] = nr["n_converged"]

    if verbose:
        n_xi_fail = int(np.isnan(xi_conv_arr).sum())
        n_n_fail = int(np.sum(np.isfinite(xi_conv_arr) & np.isnan(n_conv_arr)))
        new_solves = spectra.solves[solves_before:]
        print(f"  {len(new_solves)} DVR solves of xi^2 V "
              f"(xi up to {max([s['xi'] for s in new_solves], default=float('nan')):.4g}, "
              f"largest grid {max([s['grid_points'] for s in new_solves], default=0)} pts), "
              f"{time.time() - t0:.0f}s total")
        if n_xi_fail:
            reasons = {}
            for xr in xi_results:
                if not xr["converged"]:
                    reasons[xr["stop_reason"]] = reasons.get(xr["stop_reason"], 0) + 1
            print(f"  ⚠  ξ-convergence failed at {n_xi_fail}/{n_T} temperatures {reasons} "
                  f"-- classical limit left as NaN there.")
        if n_n_fail:
            print(f"  ⚠  n-convergence failed at {n_n_fail}/{n_T} temperatures.")
        if not n_xi_fail and not n_n_fail:
            print(f"  ✓  Both ξ and n converged at all {n_T} temperatures.")

    return {
        "cv_classical": cv_classical, "xi_conv": xi_conv_arr, "n_conv": n_conv_arr,
        "n_available": n_available, "error_estimate": error_estimate,
        "xi_results": xi_results, "n_results": n_results,
        "xi_fail_mask": np.isnan(xi_conv_arr), "n_fail_mask": np.isnan(n_conv_arr),
        "xi_ladder": xi_ladder, "spectra": spectra,
    }
