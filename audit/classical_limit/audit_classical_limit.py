"""
audit_classical_limit.py
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Audit of how the pipeline computes the CLASSICAL heat capacity
(branch correct-classical). Nothing in src/ is modified: every test
imports and calls the pipeline's own functions, unchanged, and checks
them against references that do not share any code with them.

Question 1: is the classical limit computed wrong?
Question 2: if so, why did the harmonic-oscillator (HO) validation not
            catch it?

THE TWO PRESCRIPTIONS BEING COMPARED
---------------------------------------------------------------------
Gelbwaser-Klimovsky et al., "Single-atom heat machines enabled by
energy quantization" (docs/refs), SI-III item 2 + SI-IV: scale the
potential AND the temperature, V -> xi^2 V and T -> xi^2 T. Eq. S7,

    E_n(hbar, xi^2 V) = xi^2 E_n(hbar/xi, V),

makes that hbar -> hbar/xi at fixed V and T; xi -> infinity is the
classical limit. A Boltzmann sum under this prescription therefore
needs E_n(hbar, xi^2 V): the spectrum of the RESCALED potential.

The pipeline (Classical_Limit_Numerical.converge_xi) evaluates
compute_cv(energies, beta, xi) = beta^2/xi^4 Var(E) with weights
exp(-beta E/xi^2), passing the spectrum of the UNSCALED V at every xi.

TESTS (all numbers go to RESULTS.txt, figures to figures/)
---------------------------------------------------------------------
T1  Runtime trace: which spectrum does compute_cv actually receive
    during a real sweep_temperature_range call?
T2  Algebra: compute_cv(E, beta, xi) == compute_cv(E, beta/xi^2, 1),
    i.e. the pipeline's scan is the quantum Cv at xi^2 T. Checked on
    random spectra and against the repo's own analytic Einstein Cv.
T3  Eq. S7 on the repo's DVR: E_n(hbar, xi^2 V)/xi^2 vs E_n(hbar/xi, V),
    and the pipeline's substitute E_n(hbar, V)/xi^2.
T4  Level-scaling homogeneity, gap ratio (E_n-E_0)(xi^2 V)/(E_n-E_0)(V):
    HO, box, x^4, Poschl-Teller, the double well.
T5  The HO validation, reproduced: every check the repo ran on the HO,
    plus the one it never ran, plus proof that the pipeline's xi-scan
    is IDENTICAL to the correct one for the HO (xi relabelled to xi^2).
T6  Infinite square well (exact spectrum): the case the method was first
    written for (commit e9fa799).
T7  Pure quartic x^4 (DVR): homogeneous but anharmonic.
T8  Poschl-Teller V0 tan^2 x: exact spectrum for every xi^2 V AND an
    analytic, temperature-dependent classical Cv -- the decisive test,
    no DVR involved.
T9  The project's double well (config.py), pipeline vs exact classical
    Cv and the correct xi-scaling.

The exact classical Cv of any 1-D H = p^2/2m + V(x) is
    Cv_cl/k_B = 1/2 + beta^2 Var_beta(V),   weights exp(-beta V(x)),
(the momentum integral is Gaussian; hbar cancels). It is validated here
against HO (1), x^4 (3/4) and the analytic Poschl-Teller formula.

Run:  python -u audit_classical_limit.py     (about 10 minutes)
=====================================================================
"""
import contextlib
import functools
import io
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
from scipy.special import erfcx

import config

# The audit tests the classical-limit engine AS IT WAS (main, commit
# 1ccfb7d), not the corrected one on branch correct-classical, so it loads
# that version of Classical_Limit_Numerical.py straight from git.
import subprocess
import types

AUDITED_COMMIT = "1ccfb7d"
CLN = types.ModuleType("Classical_Limit_Numerical_audited")
exec(subprocess.check_output(["git", "show", f"{AUDITED_COMMIT}:src/Classical_Limit_Numerical.py"],
                             cwd=os.path.join(HERE, "..", "..")).decode("utf-8"), CLN.__dict__)
compute_cv, sweep_temperature_range = CLN.compute_cv, CLN.sweep_temperature_range
# ...and with the xi-scan settings config.py had at that commit.
AUDITED = dict(XI_START=3.0, TOL_XI=5e-3, MIN_STABLE_XI=3, XI_MULT=1.1, MAX_XI_STEPS=80,
               TOL_CV=1e-4, MIN_STABLE_N=3)
from Quantum_Classical_Combined import compute_quantum_heat_capacity_curve
from DVR.DVR_Algorithm import auto_configure_dvr, get_fully_converged_energy_levels
from DVR.DVR_Reference_Generator import generate_reference_energies
from analytical.HO_Analytical import analytic_energy_levels_HO, analytic_cv_HO_quantum
from analytical.HO_Benchmark import compute_cv_benchmark_error
from Cv_Numerical_Benchmark import compute_cv_comparison_error
from Cv_AutoTune import resolve_beta_min

FIG_DIR = os.path.join(HERE, "figures")
RESULTS_PATH = os.path.join(HERE, "RESULTS.txt")
OUT = []

N_PIPE = config.NUM_STATES            # levels the pipeline itself uses (500)
XI_ANALYTIC = [1, 2, 4, 8, 16, 32, 64]
XI_DVR = [1, 2, 4, 8]
PT_V0 = 1.0


def log(line=""):
    print(line, flush=True)
    OUT.append(line)


def quiet(fn, *args, **kwargs):
    """Run a pipeline function with its console chatter suppressed."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


# =====================================================================
# Shared pieces
# =====================================================================
def pipeline_sweep(energies, beta_arr):
    """The pipeline's classical-limit engine with config.py's settings --
    exactly what Quantum_Classical_Combined.run and Cv_Numerical_Benchmark call."""
    return sweep_temperature_range(
        energies, beta_arr, AUDITED["XI_START"], AUDITED["TOL_XI"], AUDITED["MIN_STABLE_XI"],
        AUDITED["XI_MULT"], AUDITED["MAX_XI_STEPS"], AUDITED["TOL_CV"], AUDITED["MIN_STABLE_N"],
        verbose=False)


def s7_curve(levels_scaled, beta_arr, xi):
    """The paper's prescription with the pipeline's own compute_cv: the spectrum
    of xi^2 V, weighted at xi^2 T."""
    return np.array([compute_cv(levels_scaled, b, xi) for b in beta_arr])


def dvr_levels(V, n, hbar=1.0):
    x_min, x_max, n_grid = quiet(auto_configure_dvr, V, n, mass=1.0, hbar=hbar)
    return quiet(get_fully_converged_energy_levels, V, n, x_min, x_max, n_grid, mass=1.0, hbar=hbar)


def classical_cv_grid(V, beta_arr, x_lo, x_hi, tol=1e-10):
    """Cv_cl = 1/2 + beta^2 Var_beta(V) on a uniform grid over (x_lo, x_hi)
    (endpoints excluded), doubled until two grids agree to `tol`.
    Returns (cv, last change, points)."""
    n, cv, change = 4001, None, np.inf
    while n <= 1_000_001:
        x = np.linspace(x_lo, x_hi, n)[1:-1]
        u = V(x)
        u = u - u.min()
        new = np.empty(len(beta_arr))
        for i, b in enumerate(beta_arr):
            w = np.exp(-b * u)
            Z = w.sum()
            m = (w * u).sum() / Z
            new[i] = 0.5 + b * b * (w * (u - m) ** 2).sum() / Z
        if cv is not None:
            change = float(np.max(np.abs(new - cv)))
        cv = new
        if change < tol:
            break
        n = 2 * n - 1
    return cv, change, n


def window_on_R(V, height):
    """[x_lo, x_hi] containing every x with V - V_min <= height (V confining on R)."""
    L = 1.0
    while True:
        x = np.linspace(-L, L, 40001)
        v = V(x)
        if min(v[0], v[-1]) - v.min() > height:
            break
        L *= 1.5
    keep = np.nonzero(v - v.min() <= height)[0]
    return x[max(keep[0] - 1, 0)], x[min(keep[-1] + 1, len(x) - 1)]


def exact_classical_R(V, beta_arr):
    return classical_cv_grid(V, beta_arr, *window_on_R(V, 60.0 / beta_arr.min()))


def gap_ratio(levels_scaled, levels_unscaled, n_max=20):
    return ((levels_scaled[1:n_max + 1] - levels_scaled[0]) /
            (levels_unscaled[1:n_max + 1] - levels_unscaled[0]))


COUNTS = {}


def classify(sw):
    """
    How the pipeline produced each temperature's 'classical limit':
      plateau   converge_xi converged and the value lies on a real plateau
      tail      converge_xi reported 'converged', but the value was taken from
                the collapsing finite-N tail (< half the scan's peak) -- defect S1
      fallback  converge_xi failed, so sweep_temperature_range used converge_n
                at xi = 1, i.e. the quantum Cv at T itself -- defect S2
      none      no value at all (NaN)
    """
    kind = []
    for xr, cv in zip(sw["xi_results"], sw["cv_classical"]):
        if xr["converged"]:
            kind.append("tail" if xr["cv_converged"] < 0.5 * max(xr["cv_values"]) else "plateau")
        elif np.isfinite(cv):
            kind.append("fallback")
        else:
            kind.append("none")
    return np.array(kind)


def report_kinds(name, T, curve, ref, kind):
    """Log the main-defect error on genuine plateau points and, separately,
    how many points each secondary defect produced and how wrong they are."""
    counts = {k: int(np.sum(kind == k)) for k in ("plateau", "tail", "fallback", "none")}
    COUNTS[name] = counts
    log(f"  per-temperature origin of the pipeline's value: plateau {counts['plateau']}, "
        f"tail-pick (S1) {counts['tail']}, n-scan fallback (S2) {counts['fallback']}, NaN {counts['none']}")
    out = {}
    for k in ("plateau", "tail", "fallback"):
        if counts[k]:
            d, Td = stats_vs(T, curve, ref, mask=(kind == k))
            out[k] = d
            log(f"    max |pipeline - exact| on {k:8s} points: {d:.3g} at T = {Td:.3g}")
    return out


def stats_vs(T, curve, ref, mask=None):
    """max |curve - ref| and its T, over finite points (and `mask`)."""
    ok = np.isfinite(curve) & (np.ones_like(T, bool) if mask is None else mask)
    d = np.abs(curve - ref)[ok]
    i = int(np.argmax(d))
    return float(d[i]), float(T[ok][i])


# ---- Poschl-Teller V0 tan^2 x on (-pi/2, pi/2), hbar = m = 1 ----
def pt_levels(xi, n):
    """Exact spectrum of xi^2 V0 tan^2 x: E_k = (k + lam)^2 / 2 - xi^2 V0,
    lam = 1/2 + sqrt(1/4 + 2 xi^2 V0)  (tan^2 = sec^2 - 1)."""
    lam = 0.5 + np.sqrt(0.25 + 2.0 * xi**2 * PT_V0)
    k = np.arange(n, dtype=float)
    return 0.5 * (k + lam) ** 2 - xi**2 * PT_V0


def pt_classical_analytic(beta_arr):
    """\\int exp(-a tan^2 x) dx = pi e^a erfc(sqrt a), a = beta V0, so
    Cv_cl = 1/2 + a^2 d^2/da^2 ln erfc(sqrt a), written with erfcx for stability."""
    a = beta_arr * PT_V0
    s = np.sqrt(a)
    G = np.sqrt(np.pi) * s * erfcx(s)                 # = sqrt(pi a) e^a erfc(sqrt a)
    # With g = erfc(sqrt a): g'/g = -1/G and g''/g = (1 + 1/(2a))/G, so
    # (ln g)'' = g''/g - (g'/g)^2.
    f2 = (1.0 + 1.0 / (2.0 * a)) / G - 1.0 / G**2
    return 0.5 + a**2 * f2


# =====================================================================
# T1 -- runtime trace of what compute_cv receives
# =====================================================================
def test_T1():
    log("=" * 78)
    log("T1  Runtime trace: which spectrum does compute_cv receive inside the pipeline?")
    log("=" * 78)
    E = analytic_energy_levels_HO(N_PIPE)
    beta_arr = np.geomspace(0.1, 50.0, 5)
    calls = []
    original = CLN.compute_cv

    def recorder(energies, beta, xi=1.0):
        calls.append((np.asarray(energies, float).tobytes(), float(beta), float(xi)))
        return original(energies, beta, xi)

    CLN.compute_cv = recorder
    try:
        pipeline_sweep(E, beta_arr)
    finally:
        CLN.compute_cv = original
    spectra = {c[0] for c in calls}
    xis = sorted({c[2] for c in calls})
    log(f"  sweep_temperature_range over {len(beta_arr)} temperatures -> {len(calls)} compute_cv calls")
    log(f"  distinct xi values used: {len(xis)} (from {xis[0]:.3g} to {xis[-1]:.3g})")
    log(f"  distinct spectra passed: {len(spectra)}  "
        f"(identical to the input spectrum: {spectra == {E.tobytes()}})")
    log("  => every xi step is evaluated on the SAME, unscaled spectrum; no spectrum of xi^2 V")
    log("     is ever computed. (converge_n, the other half of the sweep, also reuses it:")
    log("     Classical_Limit_Numerical.py:238.)")


# =====================================================================
# T2 -- the algebraic identity
# =====================================================================
def test_T2():
    log("\n" + "=" * 78)
    log("T2  compute_cv(E, beta, xi) == compute_cv(E, beta/xi^2, 1): the scan is Cv_quantum(xi^2 T)")
    log("=" * 78)
    rng = np.random.RandomState(1)
    worst = 0.0
    for _ in range(20):
        E = np.sort(rng.uniform(-5, 50, 300)) + rng.normal(0, 1)
        for _ in range(50):
            b = 10 ** rng.uniform(-2, 2)
            xi = 10 ** rng.uniform(0, 2.5)
            lhs, rhs = compute_cv(E, b, xi), compute_cv(E, b / xi**2, 1.0)
            if np.isfinite(lhs) and np.isfinite(rhs):
                worst = max(worst, abs(lhs - rhs))
    log(f"  1000 random (spectrum, beta, xi) triples: max |difference| = {worst:.1e}")
    E = analytic_energy_levels_HO(200000)
    worst, T_list = 0.0, np.geomspace(0.02, 10, 40)
    for xi in (2.0, 5.0, 20.0):
        for T in T_list:
            worst = max(worst, abs(compute_cv(E, 1 / T, xi) - float(analytic_cv_HO_quantum(xi**2 * T))))
    log(f"  HO, exact levels: max |compute_cv(E, 1/T, xi) - Einstein Cv(xi^2 T)| = {worst:.1e}  "
        f"(xi = 2, 5, 20; 40 temperatures; repo's own analytic_cv_HO_quantum)")
    log("  => at every xi step the pipeline evaluates the ordinary quantum Cv of V at the")
    log("     temperature xi^2 T. The scan heats the system; it never rescales V.")


# =====================================================================
# T3 -- Eq. S7 on the repo's DVR
# =====================================================================
def test_T3():
    log("\n" + "=" * 78)
    log("T3  Eq. S7 with the repo's DVR: E_n(hbar, xi^2 V)/xi^2 vs E_n(hbar/xi, V) vs E_n(hbar, V)/xi^2")
    log("=" * 78)
    n = 150
    systems = [("double well (config.py)", functools.partial(config.my_potential, p=config.POTENTIAL_PARAMS)),
               ("HO, V = x^2/2", lambda x: 0.5 * x**2)]
    for name, V in systems:
        E_V = dvr_levels(V, n)
        for xi in (2.0, 4.0):
            E_s7 = dvr_levels(lambda x, xi=xi: xi**2 * V(x), n) / xi**2
            E_hb = dvr_levels(V, n, hbar=1.0 / xi)
            E_pipe = E_V / xi**2
            log(f"  {name}, xi = {xi:g}:  max|S7 - hbar/xi| = {np.max(np.abs(E_s7 - E_hb)):.1e};  "
                f"max|pipeline - hbar/xi| = {np.max(np.abs(E_pipe - E_hb)):.3g};  "
                f"E_0: S7 {E_s7[0]:.4f}, hbar/xi {E_hb[0]:.4f}, pipeline {E_pipe[0]:.4f}")
    log("  => Eq. S7 holds for the solver to ~1e-12, so 'scale V and T by xi^2' and 'hbar -> hbar/xi'")
    log("     are the same calculation. The pipeline's E_n(V)/xi^2 is neither, for the double well")
    log("     AND for the HO (whose levels it shrinks as hbar -> hbar/xi^2 instead; see T5).")


# =====================================================================
# T4 -- level-scaling homogeneity
# =====================================================================
def test_T4():
    log("\n" + "=" * 78)
    log("T4  Gap ratio (E_n - E_0)(xi^2 V) / (E_n - E_0)(V), n = 1..20, xi = 2")
    log("=" * 78)
    log("  Theorem (proved in AUDIT.md, checked numerically in T5): if this ratio is the same q for")
    log("  every n, the pipeline's xi-scan at xi_p equals the correct xi-scan at xi = xi_p sqrt(q(xi)),")
    log("  so both reach the same plateau. If it depends on n, no relabelling exists.")
    xi = 2.0
    dw = functools.partial(config.my_potential, p=config.POTENTIAL_PARAMS)
    rows = [
        ("HO (exact)", (xi * (np.arange(21) + 0.5)), np.arange(21) + 0.5, "xi = 2"),
        ("box (exact)", (np.pi**2 / 2) * np.arange(1, 22)**2, (np.pi**2 / 2) * np.arange(1, 22)**2, "1"),
        ("x^4 (DVR)", dvr_levels(lambda x: xi**2 * x**4, 60), dvr_levels(lambda x: x**4, 60),
         f"xi^(2/3) = {xi**(2/3):.4f}"),
        ("Poschl-Teller (exact)", pt_levels(xi, 21), pt_levels(1.0, 21), "none (inhomogeneous)"),
        ("double well (DVR)", dvr_levels(lambda x: xi**2 * dw(x), 60), dvr_levels(dw, 60), "none (inhomogeneous)"),
    ]
    ratios = {}
    for name, Es, Eu, expected in rows:
        r = gap_ratio(Es, Eu)
        ratios[name] = r
        log(f"  {name:24s} min {r.min():.6f}  max {r.max():.6f}  spread {r.max() - r.min():.1e}   "
            f"homogeneous value: {expected}")
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, r in ratios.items():
        ax.plot(np.arange(1, 21), r, marker="o", markersize=4, label=name)
    ax.set_xlabel("level n")
    ax.set_ylabel(r"$(E_n-E_0)[\xi^2 V] \;/\; (E_n-E_0)[V]$,  $\xi = 2$")
    ax.set_title("Level-scaling homogeneity: flat = the pipeline's xi-scan is exact")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(fontsize=9)
    fig.tight_layout()
    _save(fig, "T4_gap_ratios.png")


# =====================================================================
# T5 -- the HO validation, reproduced
# =====================================================================
def test_T5():
    log("\n" + "=" * 78)
    log("T5  The HO validation, reproduced with the pipeline's own functions")
    log("=" * 78)
    V = lambda x: 0.5 * x**2
    beta_arr = np.geomspace(1 / 10.0, 1 / 0.02, 500)     # the HO run's grid: T in [0.02, 10], 500 points
    T = 1 / beta_arr
    x_min, x_max, n_grid = quiet(auto_configure_dvr, V, N_PIPE)
    E_base = quiet(get_fully_converged_energy_levels, V, N_PIPE, x_min, x_max, n_grid)
    E_ref = quiet(generate_reference_energies, V, N_PIPE, x_min, x_max, n_grid,
                  span_factor=2.0, dx_factor=2.0, interactive=False, verbose=False)["energies"]
    sw_b, sw_r = pipeline_sweep(E_base, beta_arr), pipeline_sweep(E_ref, beta_arr)
    n_q = int(np.nanmax(sw_b["n_conv"]))
    cv_q = compute_quantum_heat_capacity_curve(E_base[:n_q], beta_arr)

    e_q = compute_cv_benchmark_error(cv_q, analytic_cv_HO_quantum(T))
    e_c = compute_cv_comparison_error(sw_b["cv_classical"], sw_r["cv_classical"])
    log("  Checks the repo actually ran on the HO:")
    log(f"   (i)  quantum Cv, DVR vs Einstein formula (HO_Benchmark):    max abs error {e_q['max_abs_error']:.1e}")
    log(f"   (ii) classical limit, base DVR vs reference DVR (Section 6): max rel error {e_c['max_rel']:.1e}"
        f"   <- the 'machine precision' classical figure")
    log("  Check the repo never quantified (HO_Benchmark only draws k_B next to the curve):")
    log("   (iii) classical limit vs the exact value 1:")
    report_kinds("HO (DVR, pipeline run)", T, sw_b["cv_classical"], 1.0, classify(sw_b))

    # Relabelling identity: pipeline at xi_p == correct S7 at xi = xi_p^2 for the HO.
    E_ex = analytic_energy_levels_HO(N_PIPE)
    worst = 0.0
    for xi_p in (1.5, 3.0, 7.0, 20.0):
        correct = s7_curve(xi_p**2 * E_ex, beta_arr, xi_p**2)   # spectrum of (xi_p^2)^2 V is xi_p^2 (n + 1/2)
        pipe = np.array([compute_cv(E_ex, b, xi_p) for b in beta_arr])
        worst = max(worst, float(np.max(np.abs(correct - pipe))))
    log(f"  Relabelling identity (HO): max |pipeline(xi) - correct S7(xi^2)| = {worst:.1e} over 4 xi x 500 T")
    log("  => for the HO the pipeline's scan IS the correct scan (with xi read as xi^2), and the")
    log("     correct answer is T-independent. Checks (i)-(iii) are therefore blind to the defect:")
    log("     (i) never touches the classical code, (ii) compares the code with itself, (iii) passes")
    log("     for any method that returns the high-T value, which is what the defect produces.")

    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(T, cv_q, color="#1f77b4", linewidth=2, label="quantum $C_v$ (DVR)")
    ax.plot(T, analytic_cv_HO_quantum(T), color="#ff7f0e", linewidth=1.2, linestyle="--", label="Einstein formula")
    ax.plot(T, sw_b["cv_classical"], color="#2ca02c", linewidth=2, linestyle="--", label="pipeline $\\xi$-scan, base")
    ax.plot(T, sw_r["cv_classical"], color="#d62728", linewidth=1.2, linestyle=":", label="pipeline $\\xi$-scan, reference")
    ax.axhline(1.0, color="k", linewidth=1, label="exact classical $C_v$ = 1")
    ax.set_xscale("log")
    ax.set_xlabel(r"$k_BT/\hbar\omega$")
    ax.set_ylabel(r"$C_v/k_B$")
    ax.set_title("HO: every check passes -- and would pass for the defect (T5)")
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    _save(fig, "T5_ho_validation.png")
    return {"T": T, "cv_q": cv_q, "cv_xi": sw_b["cv_classical"]}


# =====================================================================
# T6 / T7 / T8 -- box, x^4, Poschl-Teller
# =====================================================================
def run_analytic_case(name, levels_fn, beta_arr, cv_exact):
    """Pipeline xi-scan on the unscaled spectrum (N_PIPE levels, as the pipeline
    would hold it) vs the correct S7 scan (spectrum of xi^2 V) vs exact Cv_cl."""
    T = 1 / beta_arr
    sw = pipeline_sweep(levels_fn(1.0, N_PIPE), beta_arr)
    log("  pipeline xi-scan (as-is, N_PIPE exact levels) vs exact classical Cv:")
    report_kinds(name, T, sw["cv_classical"], cv_exact, classify(sw))
    s7 = {}
    for xi in XI_ANALYTIC:
        need = 20.0 / beta_arr.min() * xi**2
        n = 2000
        while True:
            E = levels_fn(float(xi), n)
            if E[-1] - E[0] >= need:
                break
            n *= 2
        s7[xi] = s7_curve(E, beta_arr, float(xi))
    hot = T >= 1.0
    devs = [float(np.max(np.abs(s7[xi] - cv_exact)[hot])) for xi in XI_ANALYTIC]
    log("  correct S7 xi-scan, max |Cv - Cv_cl| over T >= 1:  " +
        ",  ".join(f"xi={xi}: {dv:.1e}" for xi, dv in zip(XI_ANALYTIC, devs)))
    return sw["cv_classical"], s7


def test_T6():
    log("\n" + "=" * 78)
    log("T6  Infinite square well, L = 1 (exact levels pi^2 n^2 / 2; xi^2 V = V inside, walls stay infinite)")
    log("=" * 78)
    beta_arr = np.geomspace(1 / 50.0, 1 / 0.05, 300)
    levels = lambda xi, n: 0.5 * np.pi**2 * np.arange(1, n + 1, dtype=float) ** 2
    cv_xi, s7 = run_analytic_case("box", levels, beta_arr, 0.5)
    log("  => for the box, 'spectrum of xi^2 V' IS 'spectrum of V', so the pipeline is the paper's")
    log("     prescription exactly -- the method was first written for this system (commit e9fa799:")
    log("     `E_n = n**2 * E_g / xi**2`, the box spectrum at hbar/xi).")
    return beta_arr, cv_xi, s7, np.full(len(beta_arr), 0.5)


def test_T7():
    log("\n" + "=" * 78)
    log("T7  Pure quartic V = x^4 (DVR; homogeneous but anharmonic, exact Cv_cl = 3/4)")
    log("=" * 78)
    beta_arr = np.geomspace(1 / 20.0, 1 / 0.05, 300)
    T = 1 / beta_arr
    sw = pipeline_sweep(dvr_levels(lambda x: x**4, N_PIPE), beta_arr)
    log("  pipeline xi-scan vs 3/4:")
    report_kinds("x^4 (DVR)", T, sw["cv_classical"], 0.75, classify(sw))
    log("  => passes: what the pipeline needs is homogeneous level scaling (T4), not a harmonic well.")
    return beta_arr, sw["cv_classical"]


def test_T8():
    log("\n" + "=" * 78)
    log(f"T8  Poschl-Teller V = {PT_V0:g} tan^2 x on (-pi/2, pi/2): exact levels for every xi^2 V,")
    log("    analytic T-dependent classical Cv (1 at low T -> 1/2 at high T). No DVR anywhere.")
    log("=" * 78)
    # (a) the level formula, against an independent finite-difference solve
    for xi in (1.0, 4.0):
        N = 6000
        x = np.linspace(-np.pi / 2, np.pi / 2, N + 2)[1:-1]
        h = x[1] - x[0]
        Ev = eigh_tridiagonal(1 / h**2 + xi**2 * PT_V0 * np.tan(x) ** 2, -0.5 / h**2 * np.ones(N - 1),
                              select="i", select_range=(0, 9))[0]
        rel = np.max(np.abs(Ev - pt_levels(xi, 10)) / np.abs(pt_levels(xi, 10)))
        log(f"  (a) level formula vs finite differences (6000 pts), xi = {xi:g}, lowest 10 levels: "
            f"max rel diff {rel:.1e}")
    # (b) the analytic classical Cv, against the numerical phase-space integral
    beta_arr = np.geomspace(1 / 20.0, 1 / 0.02, 300)
    T = 1 / beta_arr
    cv_cl = pt_classical_analytic(beta_arr)
    cv_num, ch, npts = classical_cv_grid(lambda x: PT_V0 * np.tan(x) ** 2, beta_arr, -np.pi / 2, np.pi / 2)
    log(f"  (b) analytic Cv_cl vs numerical integral: max abs diff {np.max(np.abs(cv_cl - cv_num)):.1e} "
        f"(grid change {ch:.0e}, {npts} pts); Cv_cl runs {cv_cl[-1]:.4f} (T={T[-1]:.3g}) -> "
        f"{cv_cl[0]:.4f} (T={T[0]:.3g})")
    # (c) pipeline vs correct
    log("  (c)")
    cv_xi, s7 = run_analytic_case("PT", pt_levels, beta_arr, cv_cl)
    log("  => the pipeline's xi-scan returns the high-T value at every T; the paper's prescription,")
    log("     evaluated with the pipeline's own compute_cv, converges onto the analytic answer.")
    return beta_arr, cv_xi, s7, cv_cl


# =====================================================================
# T9 -- the project's double well
# =====================================================================
def test_T9():
    log("\n" + "=" * 78)
    log(f"T9  The project's double well, {config.POTENTIAL_PARAMS}")
    log("=" * 78)
    V = functools.partial(config.my_potential, p=config.POTENTIAL_PARAMS)
    E = dvr_levels(V, N_PIPE)
    beta_min = quiet(resolve_beta_min, config.BETA_MIN, E)
    beta_arr = np.geomspace(beta_min, config.BETA_MAX, config.N_BETA)    # the pipeline's own grid
    T = 1 / beta_arr
    sw = pipeline_sweep(E, beta_arr)
    cv_q = compute_quantum_heat_capacity_curve(E, beta_arr)
    cv_cl, ch, npts = exact_classical_R(V, beta_arr)
    log(f"  grid: T in [{T.min():.4g}, {T.max():.4g}], {len(T)} points (BETA_MIN auto, as in Quantum_HO_Master)")
    log(f"  exact classical integral: grid change {ch:.0e} ({npts} pts)")
    log(f"  pipeline range {np.nanmin(sw['cv_classical']):.4f}..{np.nanmax(sw['cv_classical']):.4f},  "
        f"exact classical range {cv_cl.min():.4f}..{cv_cl.max():.4f}")
    log("  pipeline xi-scan vs exact classical:")
    report_kinds("double well (DVR, pipeline run)", T, sw["cv_classical"], cv_cl, classify(sw))
    s7 = {}
    for xi in XI_DVR:
        need = 20.0 / beta_arr.min()
        n = int(N_PIPE * xi / 4) + 100
        while True:
            Es = dvr_levels(lambda x, xi=xi: xi**2 * V(x), n)
            if (Es[-1] - Es[0]) / xi**2 >= need:
                break
            n = int(1.3 * n)
        s7[xi] = s7_curve(Es, beta_arr, float(xi))
        T_agree = _agree_above(T, s7[xi], cv_cl)
        log(f"  correct S7 xi = {xi}: {n} levels of xi^2 V;  max|Cv - Cv_cl| over T >= 1: "
            f"{np.max(np.abs(s7[xi] - cv_cl)[T >= 1]):.1e};  within 0.01 of Cv_cl for T >= {T_agree:.3g}")
    exceed = cv_q - cv_cl
    i = int(np.argmax(exceed))
    log(f"  quantum Cv - exact classical Cv: max {exceed[i]:+.4f} at T = {T[i]:.3g}  "
        f"(vs the pipeline's curve: max {np.nanmax(cv_q - sw['cv_classical']):+.4f})")
    return beta_arr, cv_q, sw["cv_classical"], s7, cv_cl


def _agree_above(T, cv, ref, tol=0.01):
    order = np.argsort(T)[::-1]
    bad = np.nonzero(np.abs(cv - ref)[order] >= tol)[0]
    if len(bad) == 0:
        return float(T.min())
    return float(T[order][bad[0] - 1]) if bad[0] > 0 else np.nan


# =====================================================================
# Figures for T6-T9
# =====================================================================
def plot_cases(cases):
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2 * len(cases), 4.6), sharey=False)
    for ax, (title, beta_arr, cv_xi, s7, cv_cl, cv_q) in zip(axes, cases):
        T = 1 / beta_arr
        xis = sorted(s7)
        colors = plt.cm.viridis(np.linspace(0, 0.85, len(xis)))
        for xi, c in zip(xis, colors):
            ax.plot(T, s7[xi], color=c, linewidth=1.2, label=f"paper's scaling, $\\xi$={xi}")
        if cv_q is not None:
            ax.plot(T, cv_q, color="#1f77b4", linewidth=2.2, label="quantum $C_v$")
        ax.plot(T, cv_cl, color="k", linewidth=2.2, label="exact classical $C_v$")
        ax.plot(T, cv_xi, color="#d62728", linewidth=2.0, linestyle="--", label="pipeline $\\xi$-scan")
        ax.set_xscale("log")
        ax.set_ylim(0, 1.15 * max(np.nanmax(cv_cl), np.nanmax(cv_xi), 1.0))
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(r"$k_BT$")
        ax.grid(True, linestyle="--", alpha=0.4)
    axes[0].set_ylabel(r"$C_v/k_B$")
    axes[-1].legend(fontsize=7.5, loc="lower right")
    fig.suptitle("Pipeline $\\xi$-scan vs the paper's $\\xi$-scaling ($V\\to\\xi^2V$, $T\\to\\xi^2T$, "
                 "same compute_cv) vs exact classical $C_v$", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    _save(fig, "T6-T9_pipeline_vs_correct.png")


def _save(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    fig.savefig(os.path.join(FIG_DIR, name), dpi=150)
    plt.close(fig)
    log(f"  [figure] figures/{name}")


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    t0 = time.time()
    log(f"Classical-limit audit, branch correct-classical, {time.strftime('%Y-%m-%d %H:%M')}")
    log(f"Audited engine: src/Classical_Limit_Numerical.py at commit {AUDITED_COMMIT}, with that "
        f"commit's settings NUM_STATES={N_PIPE}, " + ", ".join(f"{k}={v}" for k, v in AUDITED.items()) + "\n")
    test_T1()
    test_T2()
    test_T3()
    test_T4()
    test_T5()
    b6, xi6, s6, cl6 = test_T6()
    b7, xi7 = test_T7()
    b8, xi8, s8, cl8 = test_T8()
    b9, q9, xi9, s9, cl9 = test_T9()
    plot_cases([
        ("T6  box (exact levels)", b6, xi6, s6, cl6, None),
        ("T8  Poschl-Teller (exact levels)", b8, xi8, s8, cl8, None),
        ("T9  double well, config.py (DVR)", b9, xi9, s9, cl9, q9),
    ])
    log("\n" + "=" * 78)
    log("T10 Secondary defects in the same engine, counted over every sweep above")
    log("=" * 78)
    log("  S1 tail-pick: converge_xi (Classical_Limit_Numerical.py:168-172) takes the LAST step with")
    log("     |dCv| < TOL_XI (absolute), which can sit in the collapsing finite-N tail.")
    log("  S2 fallback: sweep_temperature_range (Classical_Limit_Numerical.py:343-353) stores")
    log("     converge_n at xi = 1 -- the quantum Cv at T -- as the classical limit when xi fails.")
    for name, c in COUNTS.items():
        log(f"  {name:34s} plateau {c['plateau']:4d}   S1 tail {c['tail']:4d}   S2 fallback {c['fallback']:4d}   "
            f"NaN {c['none']:4d}")
    log(f"\nTotal time {(time.time() - t0) / 60:.1f} min")
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
