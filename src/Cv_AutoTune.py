"""
Cv_AutoTune.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
System-agnostic helpers that let Quantum_HO_Master.py treat the
temperature range as the only knob that needs hand-tuning, per system.
`NUM_STATES` (energy levels), `XI_START`, and `MAX_XI_STEPS` (the
classical-limit xi-scan) are all escalated automatically when the
pipeline's own diagnostics show a truncation artifact, instead of
requiring the user to notice and fix it by hand.

WHICH KNOB FIXES WHICH FAILURE (after the correct-classical fix)
---------------------------------------------------------------------
The two curves now depend on different things:

    - QUANTUM Cv(T): the base spectrum of NUM_STATES levels. It is
      exact as long as the top level stays thermally inaccessible at
      the hottest temperature, E_max >= HOT_STATE_SAFETY * k_B T_hot.
      Failing that -> grow NUM_STATES.
    - CLASSICAL limit: the xi-scan re-solves xi^2 V at every xi and
      sizes each of those spectra itself (Classical_Limit_Numerical.
      ScaledSpectra, same HOT_STATE_SAFETY margin), so NUM_STATES no
      longer enters it and the old finite-N collapse cannot occur.
      What can fail is the ladder running out before the plateau
      (`converge_xi`'s "max_steps"), typically at the cold end where
      the classical limit needs the largest xi -> grow XI_START /
      MAX_XI_STEPS. A ladder cut short by the hard cap XI_MAX
      ("xi_cap") is NOT escalated -- raising XI_MAX is a deliberate
      cost decision (grid sizes grow ~xi) -- and neither is a failed
      DVR solve ("dvr_failed"); both are reported.

`diagnose_escalation` below reads those signatures from the sweep's
own output (nothing new is computed here) and reports which knob(s)
need to grow. The escalation loop itself lives in Quantum_HO_Master.py.
=====================================================================
"""

import numpy as np


# =====================================================================
# T_max approximation and BETA_MIN auto-fill
# =====================================================================
def estimate_T_max(energies, k_B=1.0):
    """
    Approximate "T -> infinity" as 10x the fundamental gap E1 - E0,
    per the project's convention: we cannot compute an infinite
    temperature or an infinite number of energy levels, but by this
    point the classical-limit plateau should already be reached (not
    a hard limit -- just a practical checkpoint).

    Parameters
    ----------
    energies : array_like
        Ascending energy levels; only the two lowest are used.
    k_B : float, optional
        Boltzmann constant (default 1.0, dimensionless units).

    Returns
    -------
    T_max : float
        10 * (E1 - E0) / k_B.
    """
    delta_E = float(energies[1] - energies[0])
    return 10.0 * delta_E / k_B


def resolve_beta_min(beta_min_config, energies, k_B=1.0):
    """
    Pass through a user-supplied BETA_MIN, or auto-derive it from
    `estimate_T_max` when left as the "auto" sentinel (None).

    Parameters
    ----------
    beta_min_config : float or None
        config.BETA_MIN as written -- None means "auto".
    energies : array_like
        Current round's base energy spectrum (only E0, E1 are used).
    k_B : float, optional

    Returns
    -------
    beta_min : float
    """
    if beta_min_config is not None:
        return float(beta_min_config)

    T_max = estimate_T_max(energies, k_B)
    beta_min = k_B / T_max
    delta_E = float(energies[1] - energies[0])
    print(f"  [Auto-Tune] BETA_MIN auto-filled from T_max = 10*ΔE/k_B: "
          f"ΔE={delta_E:.4g}, T_max={T_max:.4g} -> BETA_MIN={beta_min:.4g}")
    return beta_min


# =====================================================================
# Escalation diagnosis
# =====================================================================
def diagnose_escalation(sweep, beta_arr, num_states, e_max,
                         hot_state_safety=20.0, frac_threshold=0.05,
                         n_margin_frac=0.9):
    """
    Inspect one round's `sweep_temperature_range` output and decide
    whether NUM_STATES (quantum curve) and/or XI_START/MAX_XI_STEPS
    (classical xi ladder) need to grow before the result can be
    trusted, using only diagnostics the sweep already computes (see
    module docstring for which failure maps to which knob).

    Parameters
    ----------
    sweep : dict
        Output of `Classical_Limit_Numerical.sweep_temperature_range`
        (must contain "xi_results", "n_conv" and "n_available").
    beta_arr : array_like
        The inverse-temperature array the sweep was run over.
    num_states : int
        Number of levels in this round's base spectrum (kept in the
        signature for the caller's bookkeeping; the decision only needs
        e_max).
    e_max : float
        Highest level of this round's base spectrum.
    hot_state_safety : float, optional
        Target ratio of e_max to k_B*T_hot (default 20.0, i.e. the top
        level's Boltzmann weight ~exp(-20) is already negligible).
    frac_threshold : float, optional
        Fraction of temperatures allowed to show a failure signature
        before escalation is triggered (default 0.05).
    n_margin_frac : float, optional
        A converged n within this fraction of the scaled spectrum's own
        level count counts as "no safety margin" (default 0.9).
        Reported only -- the scaled spectra are sized by the sweep itself.

    Returns
    -------
    dict with keys:
        need_states, need_xi : bool
        maxsteps_frac, xicap_frac, dvrfail_frac : float
            Fraction of temperatures whose xi-scan stopped for that reason.
        n_marginal_frac : float
            Fraction of converged temperatures whose n-check had no margin
            (or failed) on the converged xi's spectrum.
        direct_hot_coverage_fail : bool
    """
    beta_arr = np.asarray(beta_arr, dtype=float)
    stop_reasons = np.array([xr["stop_reason"] for xr in sweep["xi_results"]])
    maxsteps_frac = float(np.mean(stop_reasons == "max_steps"))
    xicap_frac = float(np.mean(stop_reasons == "xi_cap"))
    dvrfail_frac = float(np.mean(stop_reasons == "dvr_failed"))

    converged = stop_reasons == "converged"
    if converged.any():
        n_conv = np.asarray(sweep["n_conv"], dtype=float)[converged]
        n_avail = np.asarray(sweep["n_available"], dtype=float)[converged]
        n_marginal_frac = float(np.mean(np.isnan(n_conv) | (n_conv >= n_margin_frac * n_avail)))
    else:
        n_marginal_frac = 0.0

    direct_hot_coverage_fail = bool(e_max < hot_state_safety / beta_arr.min())

    return {
        "need_states": direct_hot_coverage_fail,
        "need_xi": maxsteps_frac > frac_threshold,
        "maxsteps_frac": maxsteps_frac,
        "xicap_frac": xicap_frac,
        "dvrfail_frac": dvrfail_frac,
        "n_marginal_frac": n_marginal_frac,
        "direct_hot_coverage_fail": direct_hot_coverage_fail,
    }
