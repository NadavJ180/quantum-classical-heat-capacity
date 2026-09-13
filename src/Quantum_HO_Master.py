"""
Quantum_HO_Master.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Fully numerical entry-point driver for the Harmonic Oscillator
(or any smooth potential). No analytical solutions, formulas, or
closed-form references are used anywhere in this file or in the
modules it calls. The only ground truth is a second, higher-quality
DVR computation on a finer/wider grid.

    SECTION 0     -- Configuration              (all parameters, in config.py)
    SECTION 1 & 4 -- DVR base computation + Cv pipeline, auto-tuned
                     (DVR_Algorithm, Quantum_Classical_Combined, Cv_AutoTune)
                     runs in a closed loop: solve the base spectrum, run
                     the quantum Cv(T) + xi/n classical-limit sweep, then
                     escalate NUM_STATES / XI_START / MAX_XI_STEPS and
                     retry if the sweep shows a truncation artifact (see
                     Cv_AutoTune.py); also saves the potential-shape
                     figure (plot_potential.py) once the loop settles,
                     reusing the final grid/energies -- no separate run
                     needed
    SECTION 2     -- Numerical reference        (DVR_Reference_Generator)
    SECTION 3     -- Energy-level accuracy      (error_energylevels)
                     base DVR vs reference DVR, level by level
    SECTION 5     -- DVR limit analysis         (DVR_Limit_Finder)
                     minimum dx and maximum n, both checked vs reference
    SECTION 6     -- Cv numerical benchmark     (Cv_Numerical_Benchmark)
                     quantum Cv and classical limit: base vs reference
    SECTION 7     -- Coefficient sweep          (Cv_Coefficient_Sweep)
                     quantum Cv(T) across several variants of the base
                     potential (one named POTENTIAL_PARAMS coefficient
                     swept) vs. the base run's classical limit

Sections keep their original numbering (matching FINDINGS.md/README.md);
1 and 4 are a merged, auto-tuned loop rather than two fixed, single-shot
steps run back-to-back.

GENERALITY
---------------------------------------------------------------------
Every function called here accepts any smooth V(x). To run this
pipeline on a different potential, change `my_potential` and the
label strings in config.py. Nothing else needs to change.
=====================================================================
"""

import threading
import time
import multiprocessing
import warnings

from DVR.DVR_Algorithm              import (auto_configure_dvr, 
                                            get_fully_converged_energy_levels, 
                                            colbert_miller_dvr_1d)
from error.error_energylevels       import (compute_energy_level_errors,
                                           plot_energy_level_comparison,
                                           plot_energy_level_error,
                                           print_accuracy_summary)
from Quantum_Classical_Combined     import run as run_general_cv_pipeline
from DVR.DVR_Limit_Finder           import run_dvr_limit_analysis
from DVR.DVR_Reference_Generator    import generate_reference_energies
from Cv_Numerical_Benchmark         import run_cv_numerical_benchmark
from Cv_AutoTune                    import resolve_beta_min, diagnose_escalation
from Cv_Coefficient_Sweep           import run_coefficient_sweep
from figures.output_paths           import set_context as set_figure_context
from figures.plot_potential         import plot_potential_with_spectrum
from config                         import (MASS, HBAR, my_potential,
                                            SYSTEM_NAME, T_UNITS_LABEL,
                                            POTENTIAL_PARAMS,
                                            NUM_STATES, BETA_MIN, BETA_MAX,
                                            N_BETA, XI_START, TOL_XI,
                                            MIN_STABLE_XI, XI_MULT,
                                            MAX_XI_STEPS, TOL_CV,
                                            MIN_STABLE_N, LIMIT_TOLERANCE,
                                            INTERACTIVE_REFERENCE_SCALING,
                                            REFERENCE_SPAN_FACTOR,
                                            REFERENCE_DX_FACTOR, ref_label,
                                            AUTO_ESCALATE, MAX_ESCALATION_ROUNDS,
                                            NUM_STATES_GROWTH, NUM_STATES_CAP,
                                            XI_START_GROWTH, MAX_XI_STEPS_GROWTH,
                                            HOT_STATE_SAFETY,
                                            ESCALATION_FRACTION_THRESHOLD,
                                            SCAN_PARAM, SCAN_STEP, SCAN_COUNT)


# =====================================================================
# Lightweight alive-indicator (daemon thread, prints every interval s)
# =====================================================================
class SimpleTimer:
    """
    Daemon thread that prints elapsed time every `interval` seconds
    so it is obvious the pipeline is still running during long solves.

    Parameters
    ----------
    label    : str   -- short description shown in each tick
    interval : float -- seconds between ticks (default 10)
    """
    def __init__(self, label="Running", interval=10):
        self.label    = label
        self.interval = interval
        self._stop    = threading.Event()
        self._thread  = None
        self._start   = None

    def _run(self):
        while not self._stop.wait(timeout=self.interval):
            print(f"  \u23f1  [{self.label}] still running ... "
                  f"{time.time()-self._start:.0f}s", flush=True)

    def __enter__(self):
        self._start = time.time()
        self._stop.clear()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2.0)
        print(f"  \u2713  [{self.label}] done in "
              f"{time.time()-self._start:.1f}s", flush=True)


# =====================================================================
if __name__ == "__main__":
    multiprocessing.freeze_support()

    # =================================================================
    # SECTION 0 -- Configuration
    # Moved to config.py so it can be imported by other scripts (e.g.
    # src/figures/plot_potential.py) without running this entire
    # pipeline as a side effect. Edit config.py to change system,
    # parameters, or reference scaling -- nothing here needs to change.
    #
    # Every figure saved below lands under
    # figures/<system>/<params>/<category>/ -- keyed off SYSTEM_NAME and
    # POTENTIAL_PARAMS, so a later run with different parameters (e.g. a
    # bulk scan over the double well's b) never overwrites this one.
    # =================================================================
    set_figure_context(SYSTEM_NAME, POTENTIAL_PARAMS)

    # =================================================================
    # SECTION 1 & 4 -- DVR base computation + Cv pipeline, auto-tuned
    # The only knob meant to be hand-edited run to run is the
    # temperature range (config.BETA_MAX, and optionally BETA_MIN).
    # Everything else that the Cv(T) curve's correctness depends on --
    # NUM_STATES, XI_START, MAX_XI_STEPS -- is escalated automatically
    # here: each round solves the DVR + runs the full Cv sweep, then
    # `diagnose_escalation` (Cv_AutoTune.py) inspects the sweep's own
    # existing diagnostics for the two known truncation artifacts (see
    # FINDINGS.md's Troubleshooting table) and grows the relevant
    # knob(s) before trying again, up to MAX_ESCALATION_ROUNDS.
    # =================================================================
    print("\n" + "="*60)
    print(f"  SECTION 1 & 4 — DVR base computation + Cv pipeline  ({SYSTEM_NAME})")
    print("="*60)

    num_states_cur   = NUM_STATES
    xi_start_cur     = XI_START
    max_xi_steps_cur = MAX_XI_STEPS
    beta_min_resolved = None

    for escalation_round in range(1, MAX_ESCALATION_ROUNDS + 1):
        print(f"\n  --- Auto-tune round {escalation_round}/{MAX_ESCALATION_ROUNDS} "
              f"(NUM_STATES={num_states_cur}, XI_START={xi_start_cur:.3g}, "
              f"MAX_XI_STEPS={max_xi_steps_cur}) ---")

        x_min, x_max, n_grid = auto_configure_dvr(
            my_potential, num_states_cur, mass=MASS, hbar=HBAR
        )

        with SimpleTimer(f"Round {escalation_round}: DVR 3-pass convergence check"):
            energies_base = get_fully_converged_energy_levels(
                potential_func=my_potential,
                num_levels=num_states_cur,
                x_min=x_min, x_max=x_max, num_points=n_grid,
                mass=MASS, hbar=HBAR,
            )

        if beta_min_resolved is None:
            beta_min_resolved = resolve_beta_min(BETA_MIN, energies_base)

        with SimpleTimer(f"Round {escalation_round}: Cv T-range sweep"):
            base_cv_results = run_general_cv_pipeline(
                energies=energies_base,
                system_name=SYSTEM_NAME,
                beta_min=beta_min_resolved, beta_max=BETA_MAX, n_beta=N_BETA,
                xi_start=xi_start_cur, tol_xi=TOL_XI,
                min_stable_xi=MIN_STABLE_XI,
                xi_multiplier=XI_MULT, max_xi_steps=max_xi_steps_cur,
                tol_cv=TOL_CV, min_stable_n=MIN_STABLE_N,
                cv_analytic=None,          # no analytic overlay
                T_units_label=T_UNITS_LABEL,
            )

        diag = diagnose_escalation(
            base_cv_results["sweep"], base_cv_results["beta_arr"],
            num_states_cur, energies_base[-1],
            hot_state_safety=HOT_STATE_SAFETY,
            frac_threshold=ESCALATION_FRACTION_THRESHOLD,
        )

        if not AUTO_ESCALATE or not (diag["need_states"] or diag["need_xi"]):
            break

        if escalation_round == MAX_ESCALATION_ROUNDS:
            warnings.warn(
                f"[Auto-Tune] Reached MAX_ESCALATION_ROUNDS ({MAX_ESCALATION_ROUNDS}) "
                f"without clearing the truncation diagnostics "
                f"(finite_n_hot_frac={diag['finite_n_hot_frac']:.2f}, "
                f"n_marginal_hot_frac={diag['n_marginal_hot_frac']:.2f}, "
                f"maxsteps_cold_frac={diag['maxsteps_cold_frac']:.2f}). "
                f"Proceeding with the last attempt's results.",
                UserWarning,
            )
            break

        escalate_msgs = []
        if diag["need_states"]:
            num_states_cur = min(int(num_states_cur * NUM_STATES_GROWTH), NUM_STATES_CAP)
            escalate_msgs.append(f"NUM_STATES -> {num_states_cur}")
        if diag["need_xi"]:
            xi_start_cur *= XI_START_GROWTH
            max_xi_steps_cur = int(max_xi_steps_cur * MAX_XI_STEPS_GROWTH)
            escalate_msgs.append(f"XI_START -> {xi_start_cur:.3g}, MAX_XI_STEPS -> {max_xi_steps_cur}")
        print(f"  [Auto-Tune] Escalating: {'; '.join(escalate_msgs)}")

    NUM_STATES   = num_states_cur
    XI_START     = xi_start_cur
    MAX_XI_STEPS = max_xi_steps_cur

    # Potential-shape figures (V(x) with the computed spectrum overlaid;
    # one full-spectrum overview, one zoomed on the well's own minima).
    # Reuses the final round's grid/energies -- no extra DVR solve --
    # so this never needs to be run separately via plot_potential.py.
    potential_fig_paths = plot_potential_with_spectrum(
        x_min, x_max, my_potential, energies_base, SYSTEM_NAME,
        levels_to_draw=NUM_STATES,
    )
    print(f"  Potential-shape figures saved: {potential_fig_paths['full_spectrum']}, "
          f"{potential_fig_paths['zoomed']}")

    # =================================================================
    # SECTION 2 -- Numerical reference generation
    # Run the DVR on a finer/wider grid to produce the high-precision
    # reference spectrum used as ground truth throughout Sections 3-6.
    # The reference is computed ONCE here and reused everywhere.
    #
    # NOTE on level count for Section 5 (DVR limit analysis):
    # The level-count search can test up to min(n_grid-2, NUM_STATES)
    # levels. If you want to probe beyond NUM_STATES in the limit search,
    # increase NUM_STATES or generate a separate reference with more
    # levels specifically for Section 5.
    # =================================================================
    print("\n" + "="*60)
    print("  SECTION 2 — Numerical reference generation")
    print("="*60)

    with SimpleTimer("Section 2: reference DVR solve"):
        reference_result = generate_reference_energies(
            my_potential, NUM_STATES,
            x_min, x_max, n_grid,
            span_factor=REFERENCE_SPAN_FACTOR,
            dx_factor=REFERENCE_DX_FACTOR,
            mass=MASS, hbar=HBAR,
            interactive=INTERACTIVE_REFERENCE_SCALING,
            verbose=True,
        )

    energies_ref = reference_result["energies"]

    # =================================================================
    # SECTION 3 -- Energy-level accuracy: base DVR vs numerical reference
    # Compares the two spectra level-by-level and plots the error.
    # Uses the generic functions from error_energylevels.py, which
    # only ever compare two plain arrays -- no system-specific logic.
    # =================================================================
    print("\n" + "="*60)
    print("  SECTION 3 — Energy-level accuracy (base DVR vs reference)")
    print("="*60)

    energy_error = compute_energy_level_errors(
        energies_base,
        energies_ref[:NUM_STATES],
    )
    print_accuracy_summary(
        energy_error, NUM_STATES,
        system_name=f"{SYSTEM_NAME}  [{ref_label}]",
    )
    # Plot 1: full spectrum + zoom on worst-error state
    plot_energy_level_comparison(
        energies_base, energies_ref[:NUM_STATES],
        error_result=energy_error, zoom=True,
        system_name=f"{SYSTEM_NAME} [{ref_label}]",
    )
    # Plot 2: relative error vs state index n (log y-axis)
    plot_energy_level_error(
        energy_error,
        system_name=f"{SYSTEM_NAME} — base DVR vs {ref_label}",
    )

    # =================================================================
    # SECTION 5 -- DVR limit analysis (numerical reference as truth)
    # Search A: sweep dx at fixed n → maximum safe Δx for NUM_STATES levels
    # Search B: sweep n at fixed grid → maximum trustworthy n for this Δx
    # Both searches compare against the numerical reference (Section 2).
    # =================================================================
    print("\n" + "="*60)
    print("  SECTION 5 — DVR limit analysis (vs numerical reference)")
    print("="*60)

    extend_fn = lambda n_needed: colbert_miller_dvr_1d(
        my_potential, n_needed,
        reference_result["grid"]["x_min_ref"],
        reference_result["grid"]["x_max_ref"],
        reference_result["grid"]["num_points_ref"],
        MASS, HBAR,
    )
    
    with SimpleTimer("Section 5: DVR limit searches"):
        limit_results = run_dvr_limit_analysis(
            potential_func=my_potential,
            system_name=SYSTEM_NAME,
            reference_energies=energies_ref,
            num_levels_for_grid_search=NUM_STATES,
            x_min=x_min, x_max=x_max,
            num_points_for_level_search=n_grid,
            tolerance=LIMIT_TOLERANCE,
            metric="max_abs",
            mass=MASS, hbar=HBAR,
            level_search_kwargs={"extend_reference_func": extend_fn},
        )

    # =================================================================
    # SECTION 6 -- Cv numerical benchmark
    # Runs the FULL Cv pipeline (quantum Cv + classical limit sweep)
    # on the reference energies, then compares both curves against
    # the base results from Section 4. Produces:
    #   Figure 1: quantum Cv(T) base vs reference + error panel
    #   Figure 2: classical limit Cv(T) base vs reference + error panel
    # =================================================================
    print("\n" + "="*60)
    print("  SECTION 6 — Cv numerical benchmark (base vs reference)")
    print("="*60)

    with SimpleTimer("Section 6: reference Cv sweep"):
        cv_benchmark_results = run_cv_numerical_benchmark(
            base_cv_results=base_cv_results,
            reference_energies=energies_ref,
            beta_arr=base_cv_results["beta_arr"],
            system_name=SYSTEM_NAME,
            reference_label=ref_label,
            xi_start=XI_START, tol_xi=TOL_XI,
            min_stable_xi=MIN_STABLE_XI,
            xi_multiplier=XI_MULT, max_xi_steps=MAX_XI_STEPS,
            tol_cv=TOL_CV, min_stable_n=MIN_STABLE_N,
            T_units_label=T_UNITS_LABEL,
        )

    # =================================================================
    # SECTION 7 -- Coefficient sweep comparison plot
    # Compares the quantum Cv(T) curves of several variants of the
    # base potential that differ only in one named POTENTIAL_PARAMS
    # coefficient (config.SCAN_PARAM), against the single classical-
    # limit curve already computed above for the base potential. Does
    # NOT touch the base potential's own pipeline/results above --
    # this section only ever adds this one new figure.
    # =================================================================
    print("\n" + "="*60)
    print(f"  SECTION 7 — Coefficient sweep ({SCAN_PARAM})")
    print("="*60)

    with SimpleTimer("Section 7: coefficient sweep"):
        coefficient_sweep_results = run_coefficient_sweep(
            base_cv_results=base_cv_results,
            my_potential=my_potential,
            base_params=POTENTIAL_PARAMS,
            scan_param=SCAN_PARAM, scan_step=SCAN_STEP, scan_count=SCAN_COUNT,
            num_states=NUM_STATES, mass=MASS, hbar=HBAR,
            system_name=SYSTEM_NAME, T_units_label=T_UNITS_LABEL,
        )
