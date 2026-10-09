"""
validate_quick_scan.py
=====================================================================
Does src/Quick_Scan.py give the same answer as the full pipeline, and
is its error estimate honest?

Runs Quick_Scan.run_quick_scan exactly as `python Quick_Scan.py` would,
on the double well of config.py (b swept over -0.9 .. 0, as Section 7),
at resolutions 1 and 2 (sweep) and 3 (single potential), and compares,
on each quick run's own temperature grid:

  quantum Cv     with the full pipeline's own quantum curve: the
                 NUM_STATES-level base DVR solve of each variant
                 (Cv_Coefficient_Sweep.solve_variant_spectrum, as
                 Section 7), evaluated on the same temperatures;
  classical Cv   with the exact phase-space classical Cv,
                 1/2 + beta^2 Var(V) (audit_classical_limit.
                 exact_classical_R) -- an independent check, used here
                 only to validate, not in any pipeline -- and with the
                 full pipeline's own xi-scan curve
                 (full_resolution_section7.npz, interpolated in log T);
  verdict        every temperature the quick scan calls resolved
                 ("above" or "below") must have the same sign as
                 d_exact = Cv_q(full) - Cv_cl(exact).

full_resolution_section7.npz holds the per-variant quantum and
classical curves of the full-resolution Section 7 run (Sections 1 & 4
+ 7 as the master calls them, config.py settings, 1000 temperatures,
2026-10-09; summarized in ../classical_limit/SECTION7_RESULTS.txt).

Writes QUICK_SCAN_VALIDATION.txt and the quick-scan figures under
figures/. Run from this folder:  python -u validate_quick_scan.py
(~35 minutes).
=====================================================================
"""
import functools
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "src")
sys.path.insert(0, SRC)
sys.path.insert(0, os.path.join(HERE, "..", "classical_limit"))

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import contextlib
import io

import numpy as np

import config
import figures.output_paths as output_paths
from Quick_Scan import run_quick_scan
from Cv_Coefficient_Sweep import solve_variant_spectrum
from Quantum_Classical_Combined import compute_quantum_heat_capacity_curve
from audit_classical_limit import exact_classical_R

OUT = []


def log(line=""):
    print(line, flush=True)
    OUT.append(line)


def full_pipeline_quantum(params, beta_arr):
    """Section 7's quantum Cv of one variant, on beta_arr."""
    with contextlib.redirect_stdout(io.StringIO()):
        energies, _, _ = solve_variant_spectrum(config.my_potential, params, config.NUM_STATES,
                                                config.MASS, config.HBAR)
    return compute_quantum_heat_capacity_curve(energies, beta_arr)


def check(run, full):
    """Per-variant comparison table for one quick run."""
    T, beta = run["T_arr"], run["beta_arr"]
    log(f"  {run['seconds']:.0f}s total; {len(T)} temperatures, T = {T.min():.4g} .. {T.max():.4g}")
    log("")
    log("  | b | time | DVR solves (xi max) | quantum: max abs diff to full pipeline | classical: max abs error vs exact "
        "| max error / estimate | classical: max abs diff to full-resolution xi-scan | verdict | wrong resolved verdicts |")
    log("  |---|---|---|---|---|---|---|---|---|")
    worst = {"q": 0.0, "cl": 0.0, "ratio": 0.0, "wrong": 0}
    for r in run["results"]:
        V = functools.partial(config.my_potential, p=r["params"])
        q_full = full_pipeline_quantum(r["params"], beta)
        cl_exact = exact_classical_R(V, beta)[0]
        q_diff = float(np.max(np.abs(r["cv_quantum"] - q_full)))
        ok = np.isfinite(r["cv_classical"])
        err = np.abs(r["cv_classical"] - cl_exact)
        ratio = float(np.nanmax(err[ok] / r["error_estimate"][ok])) if ok.any() else float("nan")
        k = int(np.argmin(np.abs(full["b"] - r["value"])))
        if np.isclose(full["b"][k], r["value"]):
            order = np.argsort(full["T"])           # np.interp needs ascending T
            cl_full = np.interp(np.log(T), np.log(full["T"][order]), full["classical"][k][order])
            full_text = f"{np.nanmax(np.abs(r['cv_classical'] - cl_full)):.1e}"
        else:
            full_text = "-"
        v = r["verdict"]
        d_exact = q_full - cl_exact
        wrong = int(np.sum(v["above"] & (d_exact <= 0)) + np.sum(v["below"] & (d_exact >= 0)))
        log(f"  | {r['value']:g} | {r['seconds']:.0f} s | {r['solves']} ({r['xi_top']:.4g}) | {q_diff:.1e} | "
            f"{np.nanmax(err[ok]):.1e} | {ratio:.2f} | {full_text} | {v['verdict']}"
            f"{' (hot tail unresolved)' if v['hot_tail_only'] else ''} | {wrong} |")
        worst["q"] = max(worst["q"], q_diff)
        worst["cl"] = max(worst["cl"], float(np.nanmax(err[ok])))
        worst["ratio"] = max(worst["ratio"], ratio)
        worst["wrong"] += wrong
    log("")
    log(f"  worst: quantum {worst['q']:.1e}; classical {worst['cl']:.1e} "
        f"(at most {worst['ratio']:.2f} x its own estimate); wrong resolved verdicts: {worst['wrong']}")
    return worst


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    output_paths.FIGURES_ROOT = os.path.join(HERE, "figures")
    full = np.load(os.path.join(HERE, "full_resolution_section7.npz"))

    log(f"Quick-scan validation, {time.strftime('%Y-%m-%d %H:%M')}: potential {config.POTENTIAL_PARAMS}, "
        f"sweep of {config.SCAN_PARAM} (step {config.SCAN_STEP}, count {config.SCAN_COUNT}, "
        f"symmetric value {config.SCAN_SYMMETRIC_VALUE})")
    log("Full-resolution reference (Sections 1 & 4 + 7, config.py settings, 1000 temperatures): ~44 min; "
        "max(Cv_q - Cv_cl) per b: " + ", ".join(
            f"{b:g}: {np.nanmax(q - c):+.4f}" for b, q, c in zip(full["b"], full["quantum"], full["classical"])))

    for mode, res in (("sweep", 1), ("sweep", 2), ("single", 3)):
        log("")
        log(f"## {mode}, resolution {res}: {config.QUICK_SCAN_PRESETS[res]}")
        run = run_quick_scan(mode, res)
        check(run, full)
        log(f"  figure: {os.path.relpath(run['figure_path'], HERE)}")

    with open(os.path.join(HERE, "QUICK_SCAN_VALIDATION.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
