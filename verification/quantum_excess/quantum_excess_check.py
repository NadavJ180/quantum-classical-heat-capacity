"""
quantum_excess_check.py
=====================================================================
Verification for verification/quantum_excess/SUMMARY.md:

  A. Candidates: where (if anywhere) Cv_q(T) > Cv_cl(T), for power laws
     |x|^n, the quartic double-well family, and narrow dimples in wide
     traps.
  B. The general criterion and its consequences:
     B1  Cv_q - Cv_cl = d(U_q - U_cl)/dT, on the best candidate;
     B1b the Otto-cycle work ratio W_q / W_cl over windows of that candidate;
     B2  the sum rule  int_0^inf (Cv_q - Cv_cl) dT = -(E_0 - V_min)
         (HO, exact) and its premise U_q - U_cl -> 0 (x^4 vs Wigner-Kirkwood);
     B3  the hard-wall (box) asymptote  Cv_q - 1/2 ~ 1/(8 Z_cl);
     B4  the two-region (bound state + reservoir) model of the dimples.
  C. Homogeneous eigen-energy scaling E_n(V_2) - E_0(V_2) = q (E_n(V_1) - E_0(V_1))
     under three kinds of stroke: amplitude V -> kappa V, dilation
     V(x) -> V(x / gamma), and scale-invariant V(x) -> V(x / gamma) / gamma^2.

Units hbar = m = k_B = 1. Quantum Cv: the pipeline's own Colbert-Miller DVR
(src/DVR/DVR_Algorithm.colbert_miller_dvr_1d) on grids chosen here, fine
enough for the narrow dimples (dx <= s/3). Every spectrum is checked against
a 25% wider, 1.4x finer grid; only levels agreeing to 1e-6 (relative) are
kept, and only temperatures those levels cover (E_top - E_0 >= 25 k_B T).
Classical Cv: the exact 1/2 + beta^2 Var(V) by quadrature, which is the
xi -> inf limit of the pipeline's xi-scan (CHANGES_SUMMARY.md, section F).

Writes RESULTS.md and figures/*.png next to this file. About 6 minutes.
    python verification/quantum_excess/quantum_excess_check.py
=====================================================================
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from DVR.DVR_Algorithm import colbert_miller_dvr_1d  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

FIG_DIR = os.path.join(HERE, "figures")
COVERAGE = 25.0          # keep T with E_top - E_0 >= COVERAGE * T
REPORT = []              # lines of RESULTS.md


def out(line=""):
    print(line, flush=True)
    REPORT.append(line)


# =====================================================================
# Engine
# =====================================================================
def levels(V, x_lo, x_hi, dx, nlev):
    n = int(round((x_hi - x_lo) / dx)) + 1
    return colbert_miller_dvr_1d(V, min(nlev, n - 1), x_lo, x_hi, n)


def converged_levels(V, x_lo, x_hi, dx, nlev, rtol=1e-6):
    """Levels of the base grid that a 25% wider, 1.4x finer grid confirms."""
    E1 = levels(V, x_lo, x_hi, dx, nlev)
    c, h = 0.5 * (x_lo + x_hi), 0.625 * (x_hi - x_lo)
    E2 = levels(V, c - h, c + h, dx / 1.4, nlev)
    k = min(len(E1), len(E2))
    bad = np.nonzero(np.abs(E1[:k] - E2[:k]) > rtol * np.maximum(1.0, np.abs(E2[:k])))[0]
    return E2[:(bad[0] if len(bad) else k)]


def quantum_thermo(E, T):
    """U and Cv of a discrete spectrum at each T."""
    U, C = np.empty(len(T)), np.empty(len(T))
    for i, t in enumerate(T):
        w = np.exp(-(E - E[0]) / t)
        w /= w.sum()
        U[i] = np.dot(w, E)
        C[i] = np.dot(w, (E - U[i]) ** 2) / t**2
    return U, C


def classical_thermo(V, x_lo, x_hi, T, n=400001):
    """U = T/2 + <V> and Cv = 1/2 + beta^2 Var(V), weight exp(-V/T)."""
    x = np.linspace(x_lo, x_hi, n)
    v = V(x)
    vmin = v.min()
    U, C = np.empty(len(T)), np.empty(len(T))
    for i, t in enumerate(T):
        w = np.exp(-(v - vmin) / t)
        w /= w.sum()
        m = np.dot(w, v)
        U[i] = 0.5 * t + m
        C[i] = 0.5 + np.dot(w, (v - m) ** 2) / t**2
    return U, C


def windows(T, d):
    """Temperature ranges (first, last grid T) where d > 0."""
    res, start = [], None
    for i, pos in enumerate(d > 0):
        if pos and start is None:
            start = i
        if start is not None and (not pos or i == len(d) - 1):
            res.append((T[start], T[i if pos else i - 1]))
            start = None
    return res


def fmt_windows(ws):
    return "none" if not ws else ", ".join("%.3g-%.3g" % w for w in ws)


def compare(name, V, x_lo, x_hi, dx, nlev, T):
    t0 = time.time()
    E = converged_levels(V, x_lo, x_hi, dx, nlev)
    T = T[(E[-1] - E[0]) >= COVERAGE * T]
    Uq, Cq = quantum_thermo(E, T)
    Ucl, Ccl = classical_thermo(V, x_lo, x_hi, T)
    d = Cq - Ccl
    i = int(np.argmax(d))
    rec = dict(name=name, V=V, span=(x_lo, x_hi), E=E, T=T, Uq=Uq, Cq=Cq, Ucl=Ucl, Ccl=Ccl,
               d=d, i=i, windows=windows(T, d))
    print("  %-46s %4d levels  T %.3g..%.3g  max(Cq-Ccl) = %+.4f at T = %.3g  (%.0fs)"
          % (name, len(E), T[0], T[-1], d[i], T[i], time.time() - t0), flush=True)
    return rec


def gauss(D, s):
    return lambda x: -D * np.exp(-x**2 / (2.0 * s**2))


def sech2(D, s):
    """-D sech^2(x/s), written without cosh overflow."""
    def f(x):
        e = np.exp(-2.0 * np.abs(x / s))
        return -D * 4.0 * e / (1.0 + e) ** 2
    return f


def wall(L, k):
    return lambda x: np.minimum((x / L) ** k, 1e5)


def add(*fs):
    return lambda x: sum(f(x) for f in fs)


def lam_of(D, s):
    """Poschl-Teller lam with lam (lam + 1) = 2 D s^2."""
    return 0.5 * (-1.0 + math.sqrt(1.0 + 8.0 * D * s**2))


# =====================================================================
# A. Candidates
# =====================================================================
def part_a():
    out("## A. Candidates: max(Cv_q - Cv_cl) and the windows where it is positive")
    out()
    out("### A1. Power laws |x|^n (classical Cv = 1/2 + 1/n at every T)")
    out()
    out("| n | levels | T range | max(Cv_q - Cv_cl) | at T | Cv_q / Cv_cl there | window(s) with Cv_q > Cv_cl |")
    out("|---|---|---|---|---|---|---|")
    T = np.geomspace(0.05, 40, 160)
    for n in (2, 4, 6, 8, 10, 12, 13, 14, 16, 20, 30, 50):
        xm = 1.15 * 1000.0 ** (1.0 / n)
        r = compare("|x|^%d" % n, lambda x, n=n: np.minimum(np.abs(x) ** n, 1e5),
                    -xm, xm, 0.03, 900, T)
        i = r["i"]
        out("| %d | %d | %.3g-%.3g | %+.4f | %.3g | %.3f / %.3f | %s |" % (
            n, len(r["E"]), r["T"][0], r["T"][-1], r["d"][i], r["T"][i], r["Cq"][i], r["Ccl"][i],
            fmt_windows(r["windows"])))
    out()

    out("### A2. Quartic double wells V = A [x^4 - 4/3 (1 + r) x^3 + 2 r x^2]")
    out()
    out("Minima at x = 0 and x = 1, barrier at x = r. For r > 1/2 the well at 0 is the deeper "
        "AND the stiffer one (omega_0^2 / omega_1^2 = r / (1 - r)): the dimple-like topology. "
        "A sets the depth against hbar^2/m (x in units of the well separation). Entry: "
        "max(Cv_q - Cv_cl) over T, with ZPE/gap = (hbar omega_0 / 2) / (V(1) - V(0)) in brackets.")
    out()
    As = (3, 10, 30, 100, 300, 1000, 3000)
    out("| r \\ A | " + " | ".join(str(a) for a in As) + " |")
    out("|---|" + "---|" * len(As))
    worst = -np.inf
    for r in (0.6, 0.75, 0.9, 0.97, 0.99):
        cells = []
        for A in As:
            V = lambda x, r=r, A=A: A * (x**4 - 4.0 / 3 * (1 + r) * x**3 + 2 * r * x**2)
            hw0 = math.sqrt(4 * A * r)
            gap = A * (2 * r - 1) / 3
            Ecap = 40 * max(gap, hw0)
            x = np.linspace(-20, 20, 200001)
            ok = x[V(x) < Ecap]
            xl, xr = ok[0] - 0.3, ok[-1] + 0.3
            dx = max(min(0.5 * np.pi / math.sqrt(2 * Ecap), (xr - xl) / 200), (xr - xl) / 2200)
            E = converged_levels(V, xl, xr, dx, 600, rtol=1e-5)
            T = np.geomspace((E[1] - E[0]) / 30, (E[-1] - E[0]) / COVERAGE, 80)
            d = quantum_thermo(E, T)[1] - classical_thermo(V, xl - 1, xr + 1, T)[1]
            worst = max(worst, d.max())
            cells.append("%+.4f (%.2g)" % (d.max(), 0.5 * hw0 / gap))
        out("| %.2f | " % r + " | ".join(cells) + " |")
        print("  quartic r = %.2f done" % r, flush=True)
    out()
    out("Largest value over all 35 shapes: %+.4f." % worst)
    out()

    out("### A3. Narrow dimples in wide traps")
    out()
    out("ZPE = E_0 - V_min (V_min = -D, the trap is 0 at x = 0). lam is the Poschl-Teller "
        "parameter of the sech^2 dimples, lam (lam + 1) = 2 D s^2; their exact ground energy "
        "-lam^2 / (2 s^2) is shown next to the DVR value.")
    out()
    out("| potential | bound | E_0 (DVR) | ZPE / D | max(Cv_q - Cv_cl) | at T | Cv_q / Cv_cl there | window(s) with Cv_q > Cv_cl |")
    out("|---|---|---|---|---|---|---|---|")
    T = np.geomspace(0.3, 15, 140)
    cases = [
        ("Gauss D=60 s=0.1 in x^2/2", 60, add(lambda x: 0.5 * x**2, gauss(60, 0.1)), 22, 0.025, None),
        ("Gauss D=30 s=0.1 in (x/5)^20", 30, add(wall(5, 20), gauss(30, 0.1)), 8, 0.025, None),
        ("Gauss D=30 s=0.1 in (x/10)^20", 30, add(wall(10, 20), gauss(30, 0.1)), 15, 0.025, None),
        ("Gauss D=30 s=0.1 in (x/20)^20", 30, add(wall(20, 20), gauss(30, 0.1)), 27, 0.025, None),
        ("Gauss D=30 s=0.1 in (x/10)^6", 30, add(wall(10, 6), gauss(30, 0.1)), 27, 0.025, None),
        ("Gauss D=80 s=0.1 in (x/10)^20", 80, add(wall(10, 20), gauss(80, 0.1)), 15, 0.025, None),
        ("Gauss D=60 s=0.15 in (x/10)^20", 60, add(wall(10, 20), gauss(60, 0.15)), 15, 0.03, None),
    ]
    for lam in (0.4, 0.5, 0.7, 1.0):
        D = lam * (lam + 1) / (2 * 0.1**2)
        cases.append(("sech2 lam=%.1f (D=%.4g) s=0.1 in (x/10)^6" % (lam, D), D,
                      add(wall(10, 6), sech2(D, 0.1)), 27, 0.025, -lam**2 / (2 * 0.1**2)))
    recs = {}
    for name, D, V, W, dx, e0_exact in cases:
        r = compare(name, V, -W, W, dx, 1500, T)
        r["D"] = D
        E, i = r["E"], r["i"]
        nb = int(np.sum(E < 0))
        e0 = "%.4f" % E[0] + ("" if e0_exact is None else " (exact %.4f)" % e0_exact)
        out("| %s | %d | %s | %.2f | %+.3f | %.3g | %.3f / %.3f | %s |" % (
            name, nb, e0, (E[0] + D) / D, r["d"][i], r["T"][i], r["Cq"][i], r["Ccl"][i],
            fmt_windows(r["windows"])))
        recs[name] = r
    out()
    return recs


# =====================================================================
# B. The criterion and its consequences
# =====================================================================
def part_b(best_name, best_V, D, s, W):
    out("## B. The criterion Cv_q - Cv_cl = d(U_q - U_cl)/dT and its consequences")
    out()
    # B1 -------------------------------------------------------------
    T = np.geomspace(0.3, 12, 400)
    r = compare(best_name + " (fine T grid)", best_V, -W, W, 0.025, 1500, T)
    T, dU, d = r["T"], r["Uq"] - r["Ucl"], r["d"]
    ddU = np.gradient(dU, T)
    zpe = r["E"][0] + D
    inner = slice(2, -2)
    agree = np.mean(np.sign(ddU[inner]) == np.sign(d[inner]))
    rel = np.max(np.abs(ddU[inner] - d[inner])) / np.max(np.abs(d))
    out("### B1. Best candidate, %s" % best_name)
    out()
    out("- ZPE = E_0 - V_min = %.4f; Delta U = U_q - U_cl at the coldest T (%.3g): %.4f "
        "(-> ZPE as T -> 0; here ZPE - Delta U = %.3f, the classical thermal energy of the dimple "
        "the quantum one does not yet have)." % (zpe, T[0], dU[0], zpe - dU[0]))
    out("- Finite-difference d(Delta U)/dT vs Cv_q - Cv_cl: same sign at %.1f%% of the %d interior "
        "temperatures; largest difference %.2g of max|Cv_q - Cv_cl| (discretisation of the T grid)."
        % (100 * agree, len(T) - 4, rel))
    out("- Windows where Delta U rises: %s; windows where Cv_q > Cv_cl: %s."
        % (fmt_windows(windows(T[inner], ddU[inner])), fmt_windows(windows(T, d))))
    out("- Area: int (Cv_q - Cv_cl) dT over the excess window(s) = %+.3f; over the whole range "
        "%.3g-%.3g = %+.3f (= Delta U(T_max) - Delta U(T_min) = %+.3f)." % (
            sum(trap(T[(T >= a) & (T <= b)], d[(T >= a) & (T <= b)]) for a, b in windows(T, d)),
            T[0], T[-1], trap(T, d), dU[-1] - dU[0]))
    out()
    out("### B1b. What the window is worth in an Otto cycle (homogeneous stroke q, SUMMARY.md section 4)")
    out()
    out("W = (q - 1) int_{T_c}^{T_h/q} Cv dT for both the quantum and the classical working medium "
        "(same q, so the same efficiency 1 - 1/q): W_q / W_cl = int Cv_q / int Cv_cl over [T_c, T_h/q].")
    out()
    out("| [T_c, T_h/q] | int Cv_q dT | int Cv_cl dT | W_q / W_cl |")
    out("|---|---|---|---|")
    w0 = windows(T, d)[0]
    i = r["i"]
    for a, b in (w0, (T[i] / 1.15, T[i] * 1.15), (1.0, 3.0), (0.5, 10.0)):
        m = (T >= a) & (T <= b)
        iq, ic = trap(T[m], r["Cq"][m]), trap(T[m], r["Ccl"][m])
        out("| %.3g-%.3g | %.3f | %.3f | %.3f |" % (a, b, iq, ic, iq / ic))
    out()
    plot_best(r, zpe)
    plot_delta_u(r, zpe)

    # B2 -------------------------------------------------------------
    out("### B2. Sum rule int_0^inf (Cv_q - Cv_cl) dT = -(E_0 - V_min)")
    out()
    Tg = np.geomspace(1e-4, 1e5, 200001)
    x = 1.0 / Tg
    with np.errstate(over="ignore"):
        einstein = (0.5 * x / np.sinh(0.5 * x)) ** 2
    einstein[~np.isfinite(einstein)] = 0.0
    area = trap(Tg, einstein - 1.0) - 1e-4 - 1.0 / (12 * 1e5)     # + [0, 1e-4] and the 1/(12 T^2) tail
    out("- HO (hbar omega = 1): int_0^inf (C_Einstein - 1) dT = %.6f, against -ZPE = -0.5." % area)
    T4 = np.geomspace(1.0, 40, 30)
    r4 = compare("|x|^4 (Delta U at high T)", lambda x: np.minimum(np.abs(x) ** 4, 1e5),
                 -6.5, 6.5, 0.03, 900, T4)
    c = math.gamma(0.75) / math.gamma(0.25)
    out("- The premise U_q - U_cl -> 0 at high T, on x^4 (Wigner-Kirkwood: Delta U = (hbar^2/24m) "
        "d/dbeta [beta^2 <V''>] = 0.75 (Gamma(3/4)/Gamma(1/4)) T^-1/2 = %.4f T^-1/2):" % (0.75 * c))
    out()
    out("| T | Delta U (DVR - exact classical) | Wigner-Kirkwood | ratio |")
    out("|---|---|---|---|")
    for j in np.linspace(0, len(r4["T"]) - 1, 6).astype(int):
        t = r4["T"][j]
        du = r4["Uq"][j] - r4["Ucl"][j]
        wk = 0.75 * c / math.sqrt(t)
        out("| %.3g | %.5f | %.5f | %.3f |" % (t, du, wk, du / wk))
    out()

    # B3 -------------------------------------------------------------
    out("### B3. Hard walls: the box, E_n = n^2 (units of E_1), Cv_cl = 1/2")
    out()
    E = np.arange(1, 20001, dtype=float) ** 2
    Tb = np.geomspace(0.3, 3e5, 400)
    Cb = quantum_thermo(E, Tb)[1]
    j = int(np.argmax(Cb))
    out("- Peak Cv_q = %.4f at k_B T = %.3g E_1 (Cv_q > 1/2 for every T above %.3g E_1)."
        % (Cb[j], Tb[j], Tb[np.argmax(Cb > 0.5)]))
    out("- Leading O(hbar) prediction Cv_q - 1/2 = 1/(8 Z_cl), Z_cl = (1/2) sqrt(pi T / E_1):")
    out()
    out("| T / E_1 | Cv_q - 1/2 | 1/(8 Z_cl) | ratio |")
    out("|---|---|---|---|")
    for t in (10, 100, 1e3, 1e4, 1e5):
        k = int(np.argmin(np.abs(Tb - t)))
        pred = 1.0 / (8 * 0.5 * math.sqrt(math.pi * Tb[k]))
        out("| %.3g | %.3e | %.3e | %.4f |" % (Tb[k], Cb[k] - 0.5, pred, (Cb[k] - 0.5) / pred))
    out()
    return r


def part_b4(recs):
    out("### B4. Two-region model of the dimples: Z_q ~ sum over bound states exp(-E_n/T) + G(T)")
    out()
    out("G(T) = sqrt(T/2 pi) int exp(-V_trap/T) dx is the trap's own (classical) partition function, "
        "energies from the trap bottom. T_q = quantum escape temperature, exp(E_b/T_q) = G(T_q) with "
        "E_b = -E_0. Design rules: SUMMARY.md section 2.3 (c), R2 (ZPE / T_q, the classical particle is "
        "still trapped) and R3 (ln G(T_q), a large reservoir).")
    out()
    out("| potential | E_b | ZPE | T_q | ZPE / T_q | ln G(T_q) | window (DVR) | model max|Cv - Cv_DVR| in window |")
    out("|---|---|---|---|---|---|---|---|")
    traps = {"x^2/2": lambda x: 0.5 * x**2, "(x/5)^20": wall(5, 20), "(x/10)^20": wall(10, 20),
             "(x/20)^20": wall(20, 20), "(x/10)^6": wall(10, 6)}
    for name, r in recs.items():
        trap_V = [v for k, v in traps.items() if name.endswith(k)][0]
        x = np.linspace(r["span"][0], r["span"][1], 200001)
        vt = trap_V(x)
        dxg = x[1] - x[0]
        G = lambda t: math.sqrt(t / (2 * math.pi)) * np.exp(-vt / t).sum() * dxg
        bound = r["E"][r["E"] < 0]
        eb, zpe = -bound[0], bound[0] + r["D"]
        lo, hi = 0.05, 50.0
        for _ in range(80):
            mid = math.sqrt(lo * hi)
            if eb / mid > math.log(G(mid)):
                lo = mid
            else:
                hi = mid
        Tq = math.sqrt(lo * hi)

        def lnZ(beta):
            return math.log(np.exp(-beta * bound).sum() + G(1.0 / beta))

        Cm = []
        for t in r["T"]:
            b = 1.0 / t
            h = 1e-3 * b
            Cm.append(b**2 * (lnZ(b + h) - 2 * lnZ(b) + lnZ(b - h)) / h**2)
        Cm = np.array(Cm)
        mask = r["d"] > 0
        err = np.max(np.abs(Cm - r["Cq"])[mask]) if mask.any() else float("nan")
        out("| %s | %.2f | %.2f | %.3g | %.1f | %.2f | %s | %.3f |" % (
            name, eb, zpe, Tq, zpe / Tq, math.log(G(Tq)), fmt_windows(r["windows"]), err))
    out()


def trap(x, y):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


# =====================================================================
# C. Homogeneous eigen-energy scaling under three strokes
# =====================================================================
def part_c(best_rec, best_V, W_best):
    out("## C. Homogeneous eigen-energy scaling (needed for Q = int Cv dT in an Otto cycle)")
    out()
    out("Gap ratios q_n = (E_n - E_0)[V_2] / (E_n - E_0)[V_1], n = 1..39. Homogeneous means every "
        "q_n is the same. Spread = max q_n / min q_n - 1. Strokes: amplitude V_2 = 2 V_1; dilation "
        "V_2(x) = V_1(x / 1.5); scale-invariant V_2(x) = V_1(x / 1.5) / 1.5^2 (expected q = 1/1.5^2 = 0.4444).")
    out()
    out("| potential | stroke | q_1 | min q_n | max q_n | spread | homogeneous? |")
    out("|---|---|---|---|---|---|---|")
    g = 1.5
    cands = [
        ("HO x^2/2", lambda x: 0.5 * x**2, 12.0, 0.05),
        ("|x|^20", lambda x: np.minimum(np.abs(x) ** 20, 1e5), 1.7, 0.012),
        ("double well 0.25x^4 - 0.5x^3 - 0.5x^2", lambda x: 0.25 * x**4 - 0.5 * x**3 - 0.5 * x**2, 7.5, 0.03),
        ("Gauss D=30 s=0.1 in (x/10)^20", add(wall(10, 20), gauss(30, 0.1)), 13.0, 0.025),
        ("sech2 lam=0.5 s=0.1 in (x/10)^6", best_V, 24.0, 0.025),
    ]
    strokes = [("amplitude x2", lambda V: (lambda x: 2.0 * V(x)), 1.0),
               ("dilation x1.5", lambda V: (lambda x: V(x / g)), g),
               ("scale-invariant x1.5", lambda V: (lambda x: V(x / g) / g**2), g)]
    curves = {}
    for name, V, W, dx in cands:
        E1 = converged_levels(V, -W, W, dx, 60, rtol=1e-8)
        for sname, make, wf in strokes:
            E2 = converged_levels(make(V), -W * wf, W * wf, dx, 60, rtol=1e-8)
            k = min(40, len(E1), len(E2))
            q = (E2[1:k] - E2[0]) / (E1[1:k] - E1[0])
            spread = q.max() / q.min() - 1
            out("| %s | %s | %.6f | %.6f | %.6f | %.1e | %s |" % (
                name, sname, q[0], q.min(), q.max(), spread, "yes" if spread < 1e-6 else "NO"))
            if name.startswith("sech2"):
                curves[sname] = q
        print("  homogeneity: %s done" % name, flush=True)
    out()

    # The whole Cv(T) maps under the scale-invariant stroke: Cv_2(T) = Cv_1(gamma^2 T)
    T1 = np.geomspace(0.5, 8, 40)
    V2 = (lambda x: best_V(x / g) / g**2)
    E1 = converged_levels(best_V, -W_best, W_best, 0.025, 1500)
    E2 = converged_levels(V2, -W_best * g, W_best * g, 0.025 * g, 1500)
    T1 = T1[((E1[-1] - E1[0]) >= COVERAGE * T1) & ((E2[-1] - E2[0]) >= COVERAGE * T1 / g**2)]
    cq1, cq2 = quantum_thermo(E1, T1)[1], quantum_thermo(E2, T1 / g**2)[1]
    cc1 = classical_thermo(best_V, -W_best, W_best, T1)[1]
    cc2 = classical_thermo(V2, -W_best * g, W_best * g, T1 / g**2)[1]
    out("Scale-invariant stroke on the best candidate, whole curves: max |Cv_q[V_2](T / 1.5^2) - "
        "Cv_q[V_1](T)| = %.1e and max |Cv_cl[V_2](T / 1.5^2) - Cv_cl[V_1](T)| = %.1e over "
        "T = %.3g..%.3g (%d temperatures)." % (np.max(np.abs(cq2 - cq1)), np.max(np.abs(cc2 - cc1)),
                                              T1[0], T1[-1], len(T1)))
    out()
    plot_homogeneity(curves)


# =====================================================================
# Figures (reference palette: slot 1 blue, 2 orange, 3 aqua)
# =====================================================================
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def _axes(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)
    ax.tick_params(colors=INK2)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _shade(ax, r):
    for a, b in r["windows"]:
        ax.axvspan(a, b, color=BLUE, alpha=0.10, linewidth=0)


def plot_best(r, zpe):
    fig, ax = plt.subplots(figsize=(8.5, 4.8), facecolor=SURFACE)
    _axes(ax)
    _shade(ax, r)
    ax.plot(r["T"], r["Cq"], color=BLUE, linewidth=2, label="quantum $C_v$")
    ax.plot(r["T"], r["Ccl"], color=ORANGE, linewidth=2, label="classical $C_v$ (exact $\\hbar\\to0$)")
    ax.set_xscale("log")
    ax.set_xlabel("$k_B T$   ($\\hbar = m = 1$)", color=INK)
    ax.set_ylabel("$C_v / k_B$", color=INK)
    i = r["i"]
    ax.annotate("max excess %+.2f at $k_BT$ = %.2f" % (r["d"][i], r["T"][i]),
                xy=(r["T"][i], r["Cq"][i]), xytext=(0.35, r["Cq"][i] + 0.6), color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, linewidth=0.8))
    ax.set_title("$V = (x/10)^6 - 37.5\\,\\mathrm{sech}^2(x/0.1)$: quantum above classical in the shaded window",
                 color=INK, fontsize=11)
    ax.legend(frameon=False, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "best_candidate_cv.png"), dpi=150)
    plt.close(fig)


def plot_delta_u(r, zpe):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True, facecolor=SURFACE)
    for ax in (a1, a2):
        _axes(ax)
        _shade(ax, r)
    a1.plot(r["T"], r["Uq"] - r["Ucl"], color=BLUE, linewidth=2)
    a1.axhline(zpe, color=INK2, linewidth=1, linestyle="--")
    a1.text(r["T"][0], zpe * 1.02, "zero-point energy $E_0 - V_{min}$ = %.1f" % zpe, color=INK2,
            va="bottom", fontsize=9)
    a1.set_ylabel("$\\Delta U = U_q - U_{cl}$", color=INK)
    a1.set_title("$C_v^{q} > C_v^{cl}$ exactly where $\\Delta U$ rises (shaded)", color=INK, fontsize=11)
    a2.plot(r["T"], r["d"], color=BLUE, linewidth=2)
    a2.axhline(0, color=INK2, linewidth=1)
    a2.set_ylabel("$\\Delta C_v = d\\Delta U / dT$", color=INK)
    a2.set_xscale("log")
    a2.set_xlabel("$k_B T$   ($\\hbar = m = 1$)", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "delta_u_criterion.png"), dpi=150)
    plt.close(fig)


def plot_homogeneity(curves):
    fig, ax = plt.subplots(figsize=(8.5, 4.6), facecolor=SURFACE)
    _axes(ax)
    for (sname, q), color in zip(curves.items(), (BLUE, ORANGE, AQUA)):
        n = np.arange(1, len(q) + 1)
        ax.plot(n, q, color=color, linewidth=2, marker="o", markersize=4, label=sname)
    ax.set_xlabel("level n", color=INK)
    ax.set_ylabel("$q_n = \\Delta E_n[V_2] \\,/\\, \\Delta E_n[V_1]$", color=INK)
    ax.set_title("sech$^2$ dimple in $(x/10)^6$: gap ratios under three strokes "
                 "(flat = homogeneous scaling)", color=INK, fontsize=11)
    ax.legend(frameon=False, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "homogeneity_gap_ratios.png"), dpi=150)
    plt.close(fig)


# =====================================================================
if __name__ == "__main__":
    os.makedirs(FIG_DIR, exist_ok=True)
    t_start = time.time()
    REPORT.extend(["# Quantum excess check: results", "",
                   "Written by `quantum_excess_check.py` (units hbar = m = k_B = 1). Discussion and "
                   "derivations: [SUMMARY.md](SUMMARY.md).", ""])
    recs = part_a()
    lam = 0.5
    D = lam * (lam + 1) / (2 * 0.1**2)
    best_V = add(wall(10, 6), sech2(D, 0.1))
    best = part_b("sech2 lam=0.5 (D=37.5) s=0.1 in (x/10)^6", best_V, D, 0.1, 27.0)
    part_b4(recs)
    part_c(best, best_V, 27.0)
    out("Figures: [best_candidate_cv.png](figures/best_candidate_cv.png), "
        "[delta_u_criterion.png](figures/delta_u_criterion.png), "
        "[homogeneity_gap_ratios.png](figures/homogeneity_gap_ratios.png). "
        "Total run time %.1f min." % ((time.time() - t_start) / 60))
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(REPORT) + "\n")
