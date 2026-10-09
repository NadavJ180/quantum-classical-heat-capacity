"""
merge_check.py
=====================================================================
Checks the high-temperature argument behind the quick scan's merge
limit (FINDINGS.md, "Why quantum and classical Cv merge at high T"):

    Cv_q - Cv_cl  ~  dCv_WK = -k_B beta^2 g''(beta),
    g(beta) = (hbar^2 / 24 m) beta^2 <V''>_cl(beta)

(the leading Wigner-Kirkwood correction), which for V ~ |x|^k at large
|x| falls off as T^-(1 + 2/k): T^-2 for a harmonic well, T^-1.5 for a
quartic one.

For each potential and temperature it compares
    d     = Cv_q (DVR spectrum, the pipeline's own solver via
            ScaledSpectra at xi = 1) - Cv_cl (exact phase-space value,
            1/2 + beta^2 Var(V), by quadrature)
    d_WK  = the formula above (<V''> by quadrature)
and reports their ratio, which should approach 1 as T grows, and the
log-log slope of |d| between the two hottest temperatures. The quadrature
values are used here only as an independent check; no pipeline uses them.

Writes MERGE_CHECK.txt. Run from this folder:  python -u merge_check.py
(a few minutes).
=====================================================================
"""
import contextlib
import functools
import io
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import numpy as np

import config
from Classical_Limit_Numerical import ScaledSpectra, compute_cv

OUT = []


def log(line=""):
    print(line, flush=True)
    OUT.append(line)


def _window(V, beta, n=200001):
    """x grid holding every x with beta (V - V_min) < 60, plus a margin."""
    L = 1.0
    while True:
        x = np.linspace(-L, L, 40001)
        v = V(x)
        if min(v[0], v[-1]) - v.min() > 60.0 / beta:
            break
        L *= 1.5
    keep = np.nonzero(beta * (v - v.min()) < 60.0)[0]
    lo, hi = x[max(keep[0] - 1, 0)], x[min(keep[-1] + 1, len(x) - 1)]
    return np.linspace(lo, hi, n)


def classical_and_wk(V, beta, mass=1.0, hbar=1.0):
    """Exact classical Cv/k_B and the leading Wigner-Kirkwood correction at beta."""
    x = _window(V, beta)
    v = V(x)
    dx = x[1] - x[0]
    vpp = np.gradient(np.gradient(v, dx), dx)

    def avg(f, b):
        w = np.exp(-b * (v - v.min()))
        return float(np.sum(f * w) / np.sum(w))

    cv_cl = 0.5 + beta**2 * (avg(v**2, beta) - avg(v, beta)**2)

    def g(b):
        return hbar**2 / (24.0 * mass) * b**2 * avg(vpp, b)

    h = 1e-3 * beta
    g2 = (g(beta + h) - 2.0 * g(beta) + g(beta - h)) / h**2
    return cv_cl, -beta**2 * g2


def check(name, V, temperatures, slope_expected):
    t0 = time.time()
    spectra = ScaledSpectra(V, mass=config.MASS, hbar=config.HBAR, thermal_coverage=30.0)
    with contextlib.redirect_stdout(io.StringIO()):
        energies = spectra.get(1.0, max(temperatures))
    log(f"\n{name}: {len(energies)} DVR levels (covering 30 k_B T at T = {max(temperatures):g}), "
        f"{time.time() - t0:.0f}s")
    log("  |  T | Cv_q - Cv_cl | WK prediction | ratio |")
    log("  |---|---|---|---|")
    ds = []
    for T in temperatures:
        cv_q = compute_cv(energies, 1.0 / T)
        cv_cl, d_wk = classical_and_wk(V, 1.0 / T)
        d = cv_q - cv_cl
        ds.append(d)
        log(f"  | {T:g} | {d:+.3e} | {d_wk:+.3e} | {d / d_wk:.3f} |")
    slope = math.log(abs(ds[-1]) / abs(ds[-2])) / math.log(temperatures[-1] / temperatures[-2])
    log(f"  slope of |Cv_q - Cv_cl| between the two hottest T: {slope:.2f} (expected {slope_expected})")


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    log(f"High-temperature merge check, {time.strftime('%Y-%m-%d %H:%M')}: "
        f"Cv_q - Cv_cl against the leading Wigner-Kirkwood term")
    dw = functools.partial(config.my_potential, p={"a": 0.25, "b": -0.5, "c": -0.5, "d": 0.0})
    deep = functools.partial(config.my_potential, p={"a": 0.25, "b": -10.0, "c": -0.5, "d": 0.0})
    check("HO, V = x^2/2", lambda x: 0.5 * x**2, [1, 2, 5, 10, 20], "-2")
    check("quartic, V = x^4", lambda x: x**4, [1, 2, 5, 10, 20, 40], "-1.5")
    check("double well b = -0.5", dw, [2, 5, 10, 20, 40], "-> -1.5 (quartic at large |x|)")
    check("double well b = -10", deep, [100, 200, 400, 800], "-2 (harmonic deep well, hbar omega = 30)")
    with open(os.path.join(HERE, "MERGE_CHECK.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
