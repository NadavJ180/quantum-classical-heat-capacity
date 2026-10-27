"""
check_formulas.py -- numerical checks of every formula and number quoted in
PT_Box_Summary.tex and PT_Box_Derivation.tex.   python check_formulas.py   (~30 s)

Units: hbar = m = alpha = k_B = 1, so hbar^2 alpha^2 / 2m = 1/2 and
theta = k_B T / (hbar^2 alpha^2 / 2m) = 2 T.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from cv_plot import cv_quantum, cv_classical, cv_quantum_box_only

EPS = 0.5


def exact_levels(L, n_box, alpha=1.0):
    """kappa = c tanh(kappa L/2),  k_j L = j pi + 2 arctan(k_j/c),  c = alpha tanh(alpha L/2)."""
    c = alpha * np.tanh(alpha * L / 2)
    kappa = alpha
    for _ in range(300):
        kappa = c * np.tanh(kappa * L / 2)
    j = np.arange(1, n_box + 1)
    k = j * np.pi / L
    for _ in range(300):
        k = (j * np.pi + 2 * np.arctan(k / c)) / L
    resid = np.max(np.abs(k * L - j * np.pi - 2 * np.arctan(k / c)))
    return np.concatenate([[-kappa**2 / 2], k**2 / 2]), resid


def cv_level_sum(T, E):
    cv = np.zeros(len(T))
    for i in range(len(T)):
        beta = 1 / T[i]
        p = np.exp(-beta * (E - E[0]))
        p = p / np.sum(p)
        cv[i] = beta**2 * (np.sum(p * E**2) - np.sum(p * E)**2)
    return cv


def cv_classical_quad(T, L):
    """Cv_cl = 1/2 + beta^2 Var(V) over P(x) ~ exp(-beta V), straight from the phase-space integral."""
    cv = np.zeros(len(T))
    X = 40.0                                     # V < 1e-34 beyond |x| = 40
    for i in range(len(T)):
        beta = 1 / T[i]
        V = lambda x: -1 / np.cosh(x)**2
        w = lambda x: np.exp(-beta * (V(x) + 1))     # shifted by the floor V = -1
        Z = quad(w, -X, X, limit=400, points=[0])[0] + (L - 2 * X) * np.exp(-beta)
        V1 = quad(lambda x: V(x) * w(x), -X, X, limit=400, points=[0])[0] / Z
        V2 = quad(lambda x: V(x)**2 * w(x), -X, X, limit=400, points=[0])[0] / Z
        cv[i] = 0.5 + beta**2 * (V2 - V1**2)
    return cv


print("1. Exact levels")
for L in [20.0, 60.0]:
    E, resid = exact_levels(L, 8)
    print("   alpha L = %g: E_b/(hbar^2 alpha^2/2m) = %.15f, quantization residual %.1e" % (L, E[0] / EPS, resid))
    print("      first box levels E_j/(hbar^2 alpha^2/2m):", np.round(E[1:6] / EPS, 6))
    print("      empty box (j pi/alpha L)^2:              ", np.round((np.arange(1, 6) * np.pi / L)**2, 6))

print("\n2. Closed-form Cv_q vs Var(E) over the exact levels  (theta in [0.02, 5])")
for L in [60.0, 100.0, 1000.0]:
    T = np.linspace(0.02, 5, 400) * EPS
    E, _ = exact_levels(L, int(L / np.pi * np.sqrt(2 * 60 * T.max())) + 50)
    print("   alpha L = %6g: max |difference| = %.1e k_B  (%d levels)" % (L, np.max(np.abs(cv_quantum(T, 1.0, L) - cv_level_sum(T, E))), len(E)))

print("\n3. Closed-form Cv_cl vs direct phase-space quadrature  (theta in [0.02, 5])")
for L in [100.0, 1000.0]:
    T = np.linspace(0.02, 5, 40) * EPS
    print("   alpha L = %6g: max |difference| = %.1e k_B" % (L, np.max(np.abs(cv_classical(T, 1.0, L) - cv_classical_quad(T, L)))))
a = np.array([0.1, 1.0, 5.0, 20.0, 50.0])
series = np.zeros(len(a))
term = 2 * a
for l in range(1, 400):
    if l > 1:
        term = term * 2 * a / (2 * l - 1)
    series = series + term / l
t, w = np.polynomial.legendre.leggauss(1000)
gl = np.array([np.sum(w * (np.exp(ai * (1 - t**2)) - 1) / (1 - t**2)) for ai in a])
print("   alpha*I(a): series sum (2a)^l/(l (2l-1)!!) vs Gauss-Legendre integral, max rel. diff = %.1e" % np.max(np.abs(series / gl - 1)))
print("   low T: Cv_cl(theta = 0.002) = %.4f,  Cv_q = %.1e;  high T (theta = 1000): Cv_cl = %.5f, Cv_q = %.5f"
      % (cv_classical(np.array([0.001]), 1.0, 100.0)[0], cv_quantum(np.array([0.001]), 1.0, 100.0)[0],
         cv_classical(np.array([500.0]), 1.0, 100.0)[0], cv_quantum(np.array([500.0]), 1.0, 100.0)[0]))
print("   wall term lambda_th/(8L) at theta = 1000, alpha L = 100: %.5f" % (np.sqrt(2 * np.pi / 500.0) / 800))

print("\n4. The anomaly  (theta = k_B T/(hbar^2 alpha^2/2m))")
print("   alpha L | max(Cv_q-Cv_cl) | theta  | Cv_q   | Cv_cl  | window(s) theta < 1     | Cv_q peak (theta) | Cv_cl peak (theta) | min(Cv_q-Cv_cl) (theta)")
for L in [30.0, 100.0, 400.0, 1000.0]:
    th = np.linspace(0.02, 1.0, 20000)
    cq, ccl = cv_quantum(th * EPS, 1.0, L), cv_classical(th * EPS, 1.0, L)
    d = cq - ccl
    i, iq, ic, im = np.argmax(d), np.argmax(cq), np.argmax(ccl), np.argmin(d)
    pos = d > 0
    edges = np.where(np.diff(pos.astype(int)) != 0)[0]
    win = ", ".join("%.3f" % th[e] for e in edges)
    print("   %7g | %+15.3f | %6.3f | %6.3f | %6.3f | %-23s | %5.2f (%.3f)      | %5.2f (%.3f)       | %+.2f (%.3f)"
          % (L, d[i], th[i], cq[i], ccl[i], win, cq[iq], th[iq], ccl[ic], th[ic], d[im], th[im]))

print("\n5. Equal-weight (desorption) estimates, theta_q and theta_cl")
for L in [100.0, 400.0, 1000.0]:
    th_q = brentq(lambda th: 1 / th - np.log(L * np.sqrt(th / (4 * np.pi))), 0.05, 2)
    a_cl = brentq(lambda a: a - 0.5 * np.log(a) - np.log(L / np.sqrt(np.pi)), 1.0, 50)
    print("   alpha L = %6g: theta_q = %.3f, theta_cl = %.3f" % (L, th_q, 2 / a_cl))

print("\n6. Empty box (L = 1): Cv_q - 1/2 vs r(1-r)/(2(2-r)^2),  r = lambda_th / L")
E1 = np.pi**2 / 2
T = E1 * np.logspace(-0.5, 3, 3000)
c = cv_quantum_box_only(T, 1.0)
r = np.sqrt(2 * np.pi / T)
for r_max in [1.0, 2.0 / 3.0, 0.5]:
    m = r < r_max
    print("   r < %.3f: max |difference| = %.1e" % (r_max, np.max(np.abs(c[m] - 0.5 - r[m] * (1 - r[m]) / (2 * (2 - r[m])**2)))))
i = np.argmax(c)
print("   max Cv_q = %.6f at k_B T/E_1 = %.3f   (prediction 1/2 + 1/16 = 0.5625 at 9/pi = %.3f)" % (c[i], T[i] / E1, 9 / np.pi))
print("   Cv_q > 1/2 from k_B T/E_1 = %.3f   (prediction L = lambda_th: 4/pi = %.3f)" % (T[np.argmax(c > 0.5)] / E1, 4 / np.pi))

print("\n7. Scale-invariant stroke: Cv(T; alpha/gamma, gamma L) vs Cv(gamma^2 T; alpha, L)")
T = np.linspace(0.01, 0.6, 300)
for gamma in [0.7, 1.5, 2.0]:
    dq = np.max(np.abs(cv_quantum(T, 1 / gamma, 400 * gamma) - cv_quantum(gamma**2 * T, 1.0, 400.0)))
    dc = np.max(np.abs(cv_classical(T, 1 / gamma, 400 * gamma) - cv_classical(gamma**2 * T, 1.0, 400.0)))
    E_after, _ = exact_levels(400 * gamma, 50, alpha=1 / gamma)
    E_before, _ = exact_levels(400.0, 50)
    dE = np.max(np.abs(E_after * gamma**2 - E_before))
    print("   gamma = %.1f: max |gamma^2 E_n(after) - E_n(before)| = %.1e,  max |dCv_q| = %.1e, max |dCv_cl| = %.1e"
          % (gamma, dE, dq, dc))
