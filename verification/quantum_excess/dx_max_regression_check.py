"""
dx_max_regression_check.py
=====================================================================
Regression check for the DVR_DX_MAX change (SUMMARY.md, section 5):
with the cap off (None), or on but looser than the momentum rule's own
dx, every grid-building path of the pipeline must give BIT-IDENTICAL
results to the code before the change.

Paths exercised, on the committed double well 0.25x^4 - 0.5x^3 - 0.5x^2
(pinned here, independent of the local config.py):
  1. auto_configure_dvr + get_fully_converged_energy_levels with
     NUM_STATES levels (Section 1 of the master, and
     Cv_Coefficient_Sweep.solve_variant_spectrum for Section 7);
  2. sweep_temperature_range building its own ScaledSpectra (Section 4
     via Quantum_Classical_Combined.run, and Section 7's classical curves);
  3. the same with span/dx factors 1.5 (Section 6's classical reference);
  4. Quick_Scan.scan_potential at resolution 1 (its _new_spectra cache).

    python dx_max_regression_check.py --baseline   # BEFORE the change: writes the .npz
    python dx_max_regression_check.py              # AFTER: compares, writes DX_MAX_REGRESSION.md
=====================================================================
"""
import contextlib
import functools
import inspect
import io
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))

import config  # noqa: E402
from DVR.DVR_Algorithm import auto_configure_dvr, get_fully_converged_energy_levels  # noqa: E402
from Classical_Limit_Numerical import sweep_temperature_range  # noqa: E402
import Quick_Scan  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

BASELINE = os.path.join(HERE, "dx_max_regression_baseline.npz")
PARAMS = {"a": 0.25, "b": -0.5, "c": -0.5, "d": 0.0}


def V(x):
    return PARAMS["a"] * x**4 + PARAMS["b"] * x**3 + PARAMS["c"] * x**2 + PARAMS["d"] * x


def _kw(func, dx_max):
    """{'dx_max': ...} if this version of `func` accepts it and a cap is requested."""
    if dx_max is None or "dx_max" not in inspect.signature(func).parameters:
        return {}
    return {"dx_max": dx_max}


def run_all(dx_max=None):
    out = {}
    quiet = contextlib.redirect_stdout(io.StringIO())
    with quiet:
        # 1. base spectrum, as Section 1 / solve_variant_spectrum
        x_min, x_max, n_grid = auto_configure_dvr(V, config.NUM_STATES, mass=config.MASS,
                                                  hbar=config.HBAR, **_kw(auto_configure_dvr, dx_max))
        out["grid"] = np.array([x_min, x_max, n_grid], dtype=float)
        out["energies"] = get_fully_converged_energy_levels(V, config.NUM_STATES, x_min, x_max, n_grid,
                                                            mass=config.MASS, hbar=config.HBAR)
        # 2./3. classical-limit sweeps that build their own ScaledSpectra
        beta = 1.0 / np.array([8.0, 2.0, 0.5, 0.15])
        for tag, f in (("sweep", 1.0), ("sweep_ref15", 1.5)):
            s = sweep_temperature_range(V, beta, config.XI_START, 1e-2, 2, 1.5, 25, config.TOL_CV,
                                        config.MIN_STABLE_N, mass=config.MASS, hbar=config.HBAR,
                                        thermal_coverage=12.0, xi_max=config.XI_MAX,
                                        span_factor=f, dx_factor=f, verbose=False,
                                        **_kw(sweep_temperature_range, dx_max))
            for k in ("cv_classical", "error_estimate", "xi_conv"):
                out["%s_%s" % (tag, k)] = s[k]
        # 4. the quick scan's own path (its _new_spectra reads config.DVR_DX_MAX)
        saved = getattr(config, "DVR_DX_MAX", None)
        if hasattr(config, "DVR_DX_MAX"):
            config.DVR_DX_MAX = dx_max
        try:
            settings = Quick_Scan.quick_settings(1)
            r = Quick_Scan.scan_potential(Quick_Scan.RememberedPotential(V),
                                          1.0 / np.geomspace(10.0, 0.2, 6), settings)
        finally:
            if hasattr(config, "DVR_DX_MAX"):
                config.DVR_DX_MAX = saved
        for k in ("cv_quantum", "cv_classical", "error_estimate"):
            out["quick_" + k] = r[k]
    return out


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and bool(np.all((a == b) | (np.isnan(a) & np.isnan(b))))


if __name__ == "__main__":
    t0 = time.time()
    if "--baseline" in sys.argv:
        np.savez(BASELINE, **run_all())
        print("baseline written: %s (%.0f s)" % (BASELINE, time.time() - t0))
        sys.exit(0)
    base = dict(np.load(BASELINE))
    lines = ["# DVR_DX_MAX regression check", "",
             "Written by `dx_max_regression_check.py`. Baseline: the code before the change "
             "(`dx_max_regression_baseline.npz`). Bit-identical = every element equal (NaN = NaN).", "",
             "| run | quantity | bit-identical to baseline |", "|---|---|---|"]
    ok_all = True
    for label, cap in (("DVR_DX_MAX = None", None), ("DVR_DX_MAX = 10 (looser than the momentum rule)", 10.0)):
        res = run_all(cap)
        for k in sorted(base):
            ok = same(base[k], res[k])
            ok_all &= ok
            lines.append("| %s | %s | %s |" % (label, k, "yes" if ok else "**NO**"))
    lines += ["", "All identical: %s. Run time %.0f s." % ("yes" if ok_all else "**NO**", time.time() - t0)]
    with open(os.path.join(HERE, "DX_MAX_REGRESSION.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
