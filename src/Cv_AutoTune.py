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

WHY THIS IS SAFE / GROUNDED IN EXISTING DIAGNOSTICS
---------------------------------------------------------------------
FINDINGS.md documents the classical-limit plateau's only trustworthy
window as

    sqrt(beta * DeltaE)  <<  xi  <<  sqrt(beta * E_max)

DeltaE = the spectrum's level spacing, E_max = the highest computed
level. The window's ABSOLUTE POSITION grows with sqrt(beta), but its
WIDTH (the ratio of the two bounds) does not -- it depends only on
E_max/DeltaE. That asymmetry is exactly why FINDINGS.md's own
Troubleshooting table gives two different remedies for the two ends
of the sweep:
    - cold end (large beta, low T): "Classical limit drops -> increase
      XI_START" -- a fixed XI_START can start below the window once it
      has moved to large xi at cold T, and there may not be enough
      scan steps left to climb into it (`converge_xi`'s "max_steps").
    - hot end (small beta, high T): "Classical limit / benchmark error
      -> increase NUM_STATES" -- the ABSOLUTE thermal energy k_B*T_hot
      needs E_max to stay comfortably above it, or the partition sum
      truncates thermally-accessible states (`converge_xi`'s
      "finite_n", or `converge_n` never stabilising before N).

`diagnose_escalation` below inspects exactly those two existing
failure signatures (already computed by `sweep_temperature_range` in
Classical_Limit_Numerical.py -- nothing new is computed here) on the
matching half of the temperature sweep, and reports which knob(s) need
to grow. The escalation loop itself lives in Quantum_HO_Master.py.
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
    whether NUM_STATES and/or XI_START/MAX_XI_STEPS need to grow
    before the result can be trusted, using only the diagnostics the
    sweep already computes (see module docstring for the physics).

    Parameters
    ----------
    sweep : dict
        Output of `Classical_Limit_Numerical.sweep_temperature_range`
        (must contain "xi_results" and "n_conv").
    beta_arr : array_like
        The inverse-temperature array the sweep was run over.
    num_states : int
        Number of energy levels used for this round.
    e_max : float
        Highest energy level actually computed this round.
    hot_state_safety : float, optional
        Target ratio of e_max to k_B*T_hot (default 20.0, i.e. the top
        level's Boltzmann weight ~exp(-20) is already negligible).
    frac_threshold : float, optional
        Fraction of the relevant sweep half allowed to show a failure
        signature before escalation is triggered (default 0.05).
    n_margin_frac : float, optional
        A converged n within this fraction of num_states is treated as
        "no safety margin" (default 0.9), matching FINDINGS.md's
        n-convergence diagnostic guidance.

    Returns
    -------
    dict with keys:
        need_states, need_xi : bool
        finite_n_hot_frac, n_marginal_hot_frac, maxsteps_cold_frac : float
        direct_hot_coverage_fail : bool
    """
    beta_arr = np.asarray(beta_arr, dtype=float)
    beta_median = np.median(beta_arr)
    hot_mask = beta_arr <= beta_median
    cold_mask = ~hot_mask

    xi_results = sweep["xi_results"]
    n_conv = np.asarray(sweep["n_conv"], dtype=float)
    stop_reasons = np.array([xr["stop_reason"] for xr in xi_results])

    finite_n_hot_frac = float(np.mean(stop_reasons[hot_mask] == "finite_n")) if hot_mask.any() else 0.0
    maxsteps_cold_frac = float(np.mean(stop_reasons[cold_mask] == "max_steps")) if cold_mask.any() else 0.0

    n_conv_hot = n_conv[hot_mask]
    if hot_mask.any():
        marginal = np.isnan(n_conv_hot) | (n_conv_hot >= n_margin_frac * num_states)
        n_marginal_hot_frac = float(np.mean(marginal))
    else:
        n_marginal_hot_frac = 0.0

    direct_hot_coverage_fail = bool(e_max < hot_state_safety / beta_arr.min())

    need_states = (
        finite_n_hot_frac > frac_threshold
        or n_marginal_hot_frac > frac_threshold
        or direct_hot_coverage_fail
    )
    need_xi = maxsteps_cold_frac > frac_threshold

    return {
        "need_states": bool(need_states),
        "need_xi": bool(need_xi),
        "finite_n_hot_frac": finite_n_hot_frac,
        "n_marginal_hot_frac": n_marginal_hot_frac,
        "maxsteps_cold_frac": maxsteps_cold_frac,
        "direct_hot_coverage_fail": direct_hot_coverage_fail,
    }
