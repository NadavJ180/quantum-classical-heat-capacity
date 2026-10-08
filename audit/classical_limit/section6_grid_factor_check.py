"""
section6_grid_factor_check.py
=====================================================================
Can Section 6's classical-limit reference use lighter grid factors than
Section 2's (span x2, dx /2) and still do its job?

Section 6 re-runs the corrected xi-scan with every xi^2 V solve on a grid
widened by `span_factor` and refined by `dx_factor`, and compares the
result with the base (Section 4) curve. With (2, 2) that reference took
2176 s in the full run (corrected_pipeline_run.txt) and agreed with the
base to a relative 5.5e-12. This script runs the base curve and the
reference at lighter factors on the pipeline's own 1000-point
temperature grid (config.py settings), and reports for each:

    runtime, largest grid, max / mean relative difference to the base,
    and the max relative difference to the (2, 2) value's scale.

The base curve is also compared with the exact phase-space classical Cv
(1/2 + beta^2 Var(V)) as an independent anchor.

Writes SECTION6_GRID_FACTORS.txt and figures/section6_grid_factors.png.
Run:  python -u section6_grid_factor_check.py     (~25-30 minutes)
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from Classical_Limit_Numerical import sweep_temperature_range
from DVR.DVR_Algorithm import auto_configure_dvr, get_fully_converged_energy_levels
from Cv_AutoTune import resolve_beta_min
from audit_classical_limit import exact_classical_R

FACTORS = [(1.25, 1.25), (1.5, 1.5)]
# Measured in the full pipeline run (corrected_pipeline_run.txt), not re-run here.
REFERENCE_2X2 = {"seconds": 2176.0, "max_rel": 5.458e-12, "mean_rel": 2.402e-13, "largest_grid": 12785}

OUT = []


def log(line=""):
    print(line, flush=True)
    OUT.append(line)


def sweep(V, beta_arr, span_factor, dx_factor):
    t0 = time.time()
    sw = sweep_temperature_range(
        V, beta_arr, config.XI_START, config.TOL_XI, config.MIN_STABLE_XI, config.XI_MULT,
        config.MAX_XI_STEPS, config.TOL_CV, config.MIN_STABLE_N, mass=config.MASS, hbar=config.HBAR,
        thermal_coverage=config.HOT_STATE_SAFETY, xi_max=config.XI_MAX,
        span_factor=span_factor, dx_factor=dx_factor, verbose=False)
    return sw, time.time() - t0


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    V = functools.partial(config.my_potential, p=config.POTENTIAL_PARAMS)

    # The pipeline's own temperature grid (BETA_MIN auto from the base E1 - E0).
    x_min, x_max, n_grid = auto_configure_dvr(V, config.NUM_STATES, mass=config.MASS, hbar=config.HBAR)
    E = get_fully_converged_energy_levels(V, config.NUM_STATES, x_min, x_max, n_grid,
                                          mass=config.MASS, hbar=config.HBAR)
    beta_arr = np.geomspace(resolve_beta_min(config.BETA_MIN, E), config.BETA_MAX, config.N_BETA)
    T = 1.0 / beta_arr
    log(f"Section 6 grid-factor check, {time.strftime('%Y-%m-%d %H:%M')}: potential {config.POTENTIAL_PARAMS}, "
        f"{len(T)} temperatures in [{T.min():.4g}, {T.max():.4g}]")

    base, t_base = sweep(V, beta_arr, 1.0, 1.0)
    exact = exact_classical_R(V, beta_arr)[0]
    ok = np.isfinite(base["cv_classical"])
    log(f"\nbase (Section 4 grids): {t_base:.0f}s, {len(base['spectra'].solves)} solves, largest grid "
        f"{max(s['grid_points'] for s in base['spectra'].solves)} pts, converged {ok.sum()}/{len(T)}")
    log(f"  vs exact phase-space Cv: max |diff| {np.max(np.abs(base['cv_classical'][ok] - exact[ok])):.2e} "
        f"(the xi-scan's own convergence error, set by TOL_XI = {config.TOL_XI})")

    rows = [("2 / 2 (current; full run)", REFERENCE_2X2["seconds"], REFERENCE_2X2["largest_grid"],
             REFERENCE_2X2["max_rel"], REFERENCE_2X2["mean_rel"], None)]
    curves = {}
    for span_f, dx_f in FACTORS:
        ref, t_ref = sweep(V, beta_arr, span_f, dx_f)
        both = ok & np.isfinite(ref["cv_classical"])
        rel = np.abs(base["cv_classical"] - ref["cv_classical"]) / np.abs(ref["cv_classical"])
        curves[(span_f, dx_f)] = rel
        rows.append((f"{span_f:g} / {dx_f:g}", t_ref, max(s["grid_points"] for s in ref["spectra"].solves),
                     float(np.nanmax(rel[both])), float(np.nanmean(rel[both])), int(both.sum())))
        log(f"\nreference span x{span_f:g}, dx /{dx_f:g}: {t_ref:.0f}s, converged "
            f"{int(np.isfinite(ref['cv_classical']).sum())}/{len(T)}")

    log("\n| span / dx factor | reference sweep time | largest grid | max rel. diff to base | mean rel. diff |")
    log("|---|---|---|---|---|")
    for name, secs, grid, mx, mn, _ in rows:
        log(f"| {name} | {secs / 60:.1f} min | {grid} pts | {mx:.1e} | {mn:.1e} |")

    fig, ax = plt.subplots(figsize=(9, 4.8))
    for (span_f, dx_f), rel in curves.items():
        ax.plot(T, rel, linewidth=1.2, label=f"span x{span_f:g}, dx /{dx_f:g}")
    ax.axhline(REFERENCE_2X2["max_rel"], color="k", linestyle=":", linewidth=1,
               label="max with span x2, dx /2 (full run)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(config.T_UNITS_LABEL)
    ax.set_ylabel(r"$|C_v^{base} - C_v^{ref}| \,/\, C_v^{ref}$ (classical limit)")
    ax.set_title("Section 6 classical reference: lighter grid factors vs the base curve")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "section6_grid_factors.png"), dpi=150)
    with open(os.path.join(HERE, "SECTION6_GRID_FACTORS.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
