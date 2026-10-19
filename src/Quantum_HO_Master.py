"""
Quantum_HO_Master.py
Full pipeline for the potential in config.py, with the settings config.FULL_RUN.
It runs as long as it needs (minutes to hours); nothing is changed automatically.

    SECTION 1 & 4 -- quantum Cv, classical Cv (hbar -> hbar/xi) and trust checks
    SECTION 2     -- reference spectrum (wider and finer grid)
    SECTION 3     -- energy levels: base vs reference
    SECTION 5     -- DVR limits: how many levels are correct for each dx / grid width
    SECTION 6     -- Cv: base vs reference
    SECTION 7     -- coefficient sweep (config.SCAN_PARAM over config.SCAN_VALUES)
"""

import threading
import time

import config
import potentials
from analytical.HO_Benchmark import run_ho_benchmark
from DVR.DVR_Algorithm import suggest_grid
from DVR.DVR_Reference_Generator import generate_reference_energies
from DVR.DVR_Limit_Finder import run_dvr_limit_analysis
from error.error_energylevels import (compute_energy_level_errors, plot_energy_level_comparison,
                                      print_accuracy_summary)
from Quantum_Classical_Combined import (run, print_report, print_grid_suggestion, plot_cv_summary,
                                        plot_xi_convergence, plot_n_convergence, potential_label)
from Cv_Numerical_Benchmark import run_cv_numerical_benchmark
from Cv_Coefficient_Sweep import run_coefficient_sweep
from figures.output_paths import set_context
from figures.plot_potential import plot_potential_with_spectrum


class SimpleTimer:
    """Prints the elapsed time every `interval` seconds, so long solves show they are alive."""
    def __init__(self, label="Running", interval=30):
        self.label = label
        self.interval = interval
        self._stop = threading.Event()

    def _run(self):
        while not self._stop.wait(timeout=self.interval):
            print(f"  [{self.label}] still running ... {time.time() - self._start:.0f} s", flush=True)

    def __enter__(self):
        self._start = time.time()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_):
        self._stop.set()
        self._thread.join(timeout=2.0)
        print(f"  [{self.label}] done in {time.time() - self._start:.1f} s", flush=True)


def section(title):
    print("\n" + "=" * 70 + f"\n  {title}\n" + "=" * 70, flush=True)


if __name__ == "__main__":
    t_start = time.time()
    s = config.FULL_RUN
    V = config.my_potential
    name = config.SYSTEM_NAME
    set_context(name, config.POTENTIAL_PARAMS)

    # =================================================================
    # SECTION 1 & 4 -- quantum Cv, classical Cv and trust checks
    # =================================================================
    section(f"SECTION 1 & 4 -- Cv(T): {name}, {potential_label(config.POTENTIAL_PARAMS)}")
    print(f"  settings: {s}")
    print_grid_suggestion(suggest_grid(V, config.T_MAX, config.X_MIN, config.X_MAX,
                                       config.MASS, config.HBAR, config.THERMAL_CUTOFF),
                          s["num_states"], s["dx"])
    with SimpleTimer("Section 1 & 4"):
        result = run(V, s["dx"], s["num_states"], s["xi_values"], s["n_T"], s["cv_tol"])
    print_report(result)
    title = f"{name}\n{potential_label(config.POTENTIAL_PARAMS)}, dx = {s['dx']:g}, {s['num_states']} levels"
    plot_cv_summary(result, title, "cv", "cv_summary")
    plot_xi_convergence(result)
    plot_n_convergence(result)
    plot_potential_with_spectrum(config.X_MIN, config.X_MAX, V, result["energies"], name)
    if config.POTENTIAL is potentials.harmonic:      # exact answers exist: compare with them
        run_ho_benchmark(result, config.HBAR, (config.POTENTIAL_PARAMS.get("k", 1.0) / config.MASS)**0.5)

    # =================================================================
    # SECTION 2 -- reference spectrum: wider and finer grid
    # =================================================================
    section("SECTION 2 -- Reference spectrum")
    with SimpleTimer("Section 2"):
        reference = generate_reference_energies(
            V, s["num_states"], config.X_MIN, config.X_MAX, s["dx"],
            config.REFERENCE_SPAN_FACTOR, config.REFERENCE_DX_FACTOR,
            config.MASS, config.HBAR, config.WALLS)

    # =================================================================
    # SECTION 3 -- energy levels: base vs reference
    # =================================================================
    section("SECTION 3 -- Energy levels: base vs reference")
    level_errors = compute_energy_level_errors(result["energies"], reference["energies"])
    print_accuracy_summary(level_errors, name)
    plot_energy_level_comparison(result["energies"], reference["energies"], level_errors, name,
                                 config.LIMIT_TOLERANCE)

    # =================================================================
    # SECTION 5 -- DVR limits: correct levels for each dx and grid width
    # =================================================================
    section("SECTION 5 -- DVR limit analysis")
    with SimpleTimer("Section 5"):
        run_dvr_limit_analysis(V, name, reference["energies"], s["num_states"], config.X_MIN, config.X_MAX,
                               s["dx"], config.LIMIT_DX_FACTORS, config.LIMIT_SPAN_FACTORS,
                               config.LIMIT_TOLERANCE, config.MASS, config.HBAR, config.WALLS)

    # =================================================================
    # SECTION 6 -- Cv: base vs reference
    # =================================================================
    section("SECTION 6 -- Cv benchmark: base vs reference")
    with SimpleTimer("Section 6"):
        run_cv_numerical_benchmark(result, reference, V, name, config.RUN_CLASSICAL_REFERENCE)

    # =================================================================
    # SECTION 7 -- coefficient sweep
    # =================================================================
    if config.SCAN_VALUES:
        section(f"SECTION 7 -- Coefficient sweep: {config.SCAN_PARAM} = {config.SCAN_VALUES}")
        with SimpleTimer("Section 7"):
            run_coefficient_sweep(s["dx"], s["num_states"], s["xi_values"], s["n_T"], s["cv_tol"])

    print(f"\n  Full pipeline finished in {(time.time() - t_start) / 60:.1f} min")
