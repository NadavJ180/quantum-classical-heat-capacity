"""
Quick_Scan.py
A fast first look at a potential (a few minutes at most): the quantum and
classical Cv(T) with the coarse settings config.QUICK_SCAN, where each curve
can be trusted, and where the quantum Cv might be above the classical one.

    python Quick_Scan.py                    # config.QUICK_SCAN_MODE
    python Quick_Scan.py --mode sweep       # every value in config.SCAN_VALUES
"""

import argparse
import time

import config
from DVR.DVR_Algorithm import suggest_grid
from Quantum_Classical_Combined import (run, print_report, print_grid_suggestion,
                                        plot_cv_summary, potential_label)
from Cv_Coefficient_Sweep import run_coefficient_sweep
from figures.output_paths import set_context


def quick_scan(mode):
    t_start = time.time()
    s = config.QUICK_SCAN
    set_context(config.SYSTEM_NAME, config.POTENTIAL_PARAMS)
    print("\n" + "=" * 70)
    print(f"  QUICK SCAN ({mode}): {config.SYSTEM_NAME}, {potential_label(config.POTENTIAL_PARAMS)}")
    print(f"  T = {config.T_MIN:g} .. {config.T_MAX:g} ({s['n_T']} points), settings {s}")
    print("=" * 70)
    print_grid_suggestion(suggest_grid(config.my_potential, config.T_MAX, config.X_MIN, config.X_MAX,
                                       config.MASS, config.HBAR, config.THERMAL_CUTOFF),
                          s["num_states"], s["dx"])

    if mode == "sweep":
        run_coefficient_sweep(s["dx"], s["num_states"], s["xi_values"], s["n_T"], s["cv_tol"],
                              category="quick_scan", name="quick_sweep")
    else:
        result = run(config.my_potential, s["dx"], s["num_states"], s["xi_values"], s["n_T"], s["cv_tol"])
        print_report(result)
        title = (f"{config.SYSTEM_NAME}: quick scan\n{potential_label(config.POTENTIAL_PARAMS)}, "
                 f"dx = {s['dx']:g}, {s['num_states']} levels")
        path = plot_cv_summary(result, title, "quick_scan", "quick_single")
        print(f"\n  Figure: {path}")
    print(f"  Quick scan finished in {time.time() - t_start:.0f} s\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quick quantum vs classical Cv(T) scan (settings in config.py)")
    parser.add_argument("--mode", choices=["single", "sweep"], default=config.QUICK_SCAN_MODE)
    quick_scan(parser.parse_args().mode)
