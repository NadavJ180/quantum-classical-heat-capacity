"""
verify_corrected_engine.py
=====================================================================
Checks the CORRECTED classical-limit engine on branch correct-classical
(src/Classical_Limit_Numerical.py: xi^2 V re-solved at every xi, fixed
plateau picker, no fallback) against exact classical answers, with the
pipeline's own config.py settings:

    HO, V = x^2/2         exact Cv_cl = 1      (homogeneous)
    quartic, V = x^4      exact Cv_cl = 3/4    (homogeneous, anharmonic)
    the double well       exact Cv_cl from the phase-space integral
                          1/2 + beta^2 Var(V)  (inhomogeneous -- the case
                          the old engine got wrong)

For every converged temperature it reports the true error and the
engine's own error estimate, which should agree (hbar^2 law).

Run:  python -u verify_corrected_engine.py     (~15-20 minutes)
Writes VERIFY_CORRECTED.txt next to this script.
=====================================================================
"""
import functools
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import numpy as np

import config
from Classical_Limit_Numerical import sweep_temperature_range
from audit_classical_limit import exact_classical_R

N_T = 40
OUT = []


def log(line=""):
    print(line, flush=True)
    OUT.append(line)


def check(name, V, beta_arr, exact):
    t0 = time.time()
    sw = sweep_temperature_range(
        V, beta_arr, config.XI_START, config.TOL_XI, config.MIN_STABLE_XI, config.XI_MULT,
        config.MAX_XI_STEPS, config.TOL_CV, config.MIN_STABLE_N, mass=config.MASS, hbar=config.HBAR,
        thermal_coverage=config.HOT_STATE_SAFETY, xi_max=config.XI_MAX, verbose=False)
    T = 1.0 / beta_arr
    ok = ~sw["xi_fail_mask"]
    err = sw["cv_classical"] - exact
    est = sw["error_estimate"]
    solves = sw["spectra"].solves
    stops = {}
    for xr in sw["xi_results"]:
        stops[xr["stop_reason"]] = stops.get(xr["stop_reason"], 0) + 1
    log(f"\n{name}: T in [{T.min():.3g}, {T.max():.3g}], {len(T)} log-spaced points")
    log(f"  {time.time() - t0:.0f}s; {len(solves)} DVR solves of xi^2 V; levels "
        f"{min(s['num_levels'] for s in solves)}..{max(s['num_levels'] for s in solves)}; "
        f"grid {min(s['grid_points'] for s in solves)}..{max(s['grid_points'] for s in solves)} pts")
    log(f"  xi-scan outcomes: {stops};  n-check failed at {int(np.sum(ok & sw['n_fail_mask']))}")
    if ok.any():
        log(f"  max |Cv - exact| over converged T: {np.max(np.abs(err[ok])):.2e};  "
            f"max |true error - estimate|: {np.max(np.abs(np.abs(err[ok]) - est[ok])):.1e};  "
            f"xi_conv {np.nanmin(sw['xi_conv']):.3g}..{np.nanmax(sw['xi_conv']):.4g}")
    log(f"  {'T':>8} {'pipeline':>9} {'exact':>9} {'error':>9} {'estimate':>9} {'xi_conv':>8}")
    for k in range(0, len(T), 4):
        log(f"  {T[k]:8.4f} {sw['cv_classical'][k]:9.5f} {exact[k]:9.5f} {err[k]:+9.1e} "
            f"{est[k]:9.1e} {sw['xi_conv'][k]:8.4g}")


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    log(f"Corrected classical-limit engine, {time.strftime('%Y-%m-%d %H:%M')}; config.py settings: "
        f"XI_START={config.XI_START}, XI_MULT={config.XI_MULT}, TOL_XI={config.TOL_XI}, "
        f"MIN_STABLE_XI={config.MIN_STABLE_XI}, MAX_XI_STEPS={config.MAX_XI_STEPS}, XI_MAX={config.XI_MAX}, "
        f"HOT_STATE_SAFETY={config.HOT_STATE_SAFETY}")

    beta_ho = np.geomspace(1 / 10.0, 1 / 0.02, N_T)
    check("HO, V = x^2/2 (exact 1)", lambda x: 0.5 * x**2, beta_ho, np.ones(N_T))

    beta_q = np.geomspace(1 / 20.0, 1 / 0.05, N_T)
    check("quartic, V = x^4 (exact 3/4)", lambda x: x**4, beta_q, np.full(N_T, 0.75))

    V = functools.partial(config.my_potential, p=config.POTENTIAL_PARAMS)
    beta_dw = np.geomspace(1 / 13.3, config.BETA_MAX, N_T)      # the pipeline's T range
    check(f"double well {config.POTENTIAL_PARAMS} (exact: phase-space integral)",
          V, beta_dw, exact_classical_R(V, beta_dw)[0])

    with open(os.path.join(HERE, "VERIFY_CORRECTED.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
