"""
check_explained.py -- symbolic and numerical checks of every claim in PT_Box_Explained.tex.
Each block prints its check id C1, C2, ...; the text cites them as "(C<n>)".
At the end it draws figures/fig_explained_levels.pdf (and .png).

    python check_explained.py        (about 1 min; numpy, scipy, sympy, mpmath, matplotlib)

Numerical units: hbar = m = alpha = k_B = 1 unless stated.  Then a = beta hbar^2 alpha^2/m = 1/T,
hbar^2 alpha^2/2m = 1/2, and k_B T/(hbar^2 alpha^2/2m) = 2/a.
"""

import sys
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erf, erfc

from cv_plot import cv_quantum, cv_classical, cv_quantum_box_only

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass


def zero(expr, x, extra_subs=None):
    """True if expr == 0 at 8 random points (30-digit evaluation; parameters from extra_subs).
    The expressions are analytic, so vanishing at random points means the identity holds."""
    subs = dict(extra_subs or {})
    rng = np.random.RandomState(1)
    for _ in range(8):
        pt = dict(subs)
        pt[x] = sp.Rational(int(rng.uniform(-2000, 2000)), 1000)   # exact input: no 15-digit rounding
        val = complex(sp.N(expr.subs(pt), 40))
        if abs(val) > 1e-25:
            return False
    return True


def exact_levels(L, n_box, alpha=1.0):
    """Box levels k_j L = j pi + 2 arctan(k_j/c) and bound state kappa = c tanh(kappa L/2), c = alpha tanh(alpha L/2)."""
    c = alpha * np.tanh(alpha * L / 2)
    kappa = alpha
    for _ in range(300):
        kappa = c * np.tanh(kappa * L / 2)
    j = np.arange(1, n_box + 1)
    k = j * np.pi / L
    for _ in range(300):
        k = (j * np.pi + 2 * np.arctan(k / c)) / L
    return c, kappa, k


def cv_level_sum(T, E):
    cv = np.zeros(len(T))
    for i in range(len(T)):
        beta = 1 / T[i]
        p = np.exp(-beta * (E - E[0]))
        p = p / np.sum(p)
        cv[i] = beta**2 * (np.sum(p * E**2) - np.sum(p * E)**2)
    return cv


print("=" * 100)
print("Symbolic checks (sympy)")
x, k, kap, al, hb, m, a, s, c = sp.symbols('x k kappa alpha hbar m a s c', positive=True)
f = sp.Function('f')(x)
W = hb * al / sp.sqrt(2 * m) * sp.tanh(al * x)
A = lambda g: hb / sp.sqrt(2 * m) * sp.diff(g, x) + W * g
Ad = lambda g: -hb / sp.sqrt(2 * m) * sp.diff(g, x) + W * g
V = -hb**2 * al**2 / m / sp.cosh(al * x)**2
H_PT = lambda g: -hb**2 / (2 * m) * sp.diff(g, x, 2) + V * g
H_0 = lambda g: -hb**2 / (2 * m) * sp.diff(g, x, 2)
E0 = hb**2 * al**2 / (2 * m)
num = {hb: sp.Rational(11, 10), m: sp.Rational(7, 10), al: sp.Rational(13, 10), k: sp.Rational(9, 10), kap: sp.Rational(3, 5)}
ftest = sp.sin(sp.Rational(17, 10) * x) + x**3 * sp.exp(-x**2)        # exact coefficients (see zero())

print("C1  hyperbolic identities:",
      zero(sp.cosh(x)**2 - sp.sinh(x)**2 - 1, x),
      zero(1 - sp.tanh(x)**2 - 1 / sp.cosh(x)**2, x),
      zero(sp.diff(sp.tanh(x), x) - 1 / sp.cosh(x)**2, x),
      zero(sp.cosh(x) - sp.tanh(x) * sp.sinh(x) - 1 / sp.cosh(x), x),
      zero(sp.sinh(2 * x) - 2 * sp.sinh(x) * sp.cosh(x), x))

lhs1 = Ad(A(ftest)) - (H_PT(ftest) + E0 * ftest)
lhs2 = A(Ad(ftest)) - (H_0(ftest) + E0 * ftest)
print("C2  A^dag A = H_PT + hbar^2 alpha^2/2m:", zero(lhs1, x, num), "  A A^dag = H_0 + hbar^2 alpha^2/2m:", zero(lhs2, x, num))
print("    generic f: A^dag A f = -hbar^2 f''/2m + (W^2 - hbar W'/sqrt(2m)) f:",
      sp.simplify(sp.expand(Ad(A(f)) - (-hb**2 / (2 * m) * sp.diff(f, x, 2) + (W**2 - hb / sp.sqrt(2 * m) * sp.diff(W, x)) * f))) == 0,
      "  A A^dag f = -hbar^2 f''/2m + (W^2 + hbar W'/sqrt(2m)) f:",
      sp.simplify(sp.expand(A(Ad(f)) - (-hb**2 / (2 * m) * sp.diff(f, x, 2) + (W**2 + hb / sp.sqrt(2 * m) * sp.diff(W, x)) * f))) == 0)
print("C3  H_PT A^dag = A^dag H_0:", zero(H_PT(Ad(ftest)) - Ad(H_0(ftest)), x, num),
      "   A H_PT = H_0 A:", zero(A(H_PT(ftest)) - H_0(A(ftest)), x, num),
      "   generic f (symbolic):", sp.simplify(sp.expand((H_PT(Ad(f)) - Ad(H_0(f))).rewrite(sp.exp))) == 0,
      sp.simplify(sp.expand((A(H_PT(f)) - H_0(A(f))).rewrite(sp.exp))) == 0)

psi_e = k * sp.cos(k * x) - al * sp.tanh(al * x) * sp.sin(k * x)
psi_o = al * sp.tanh(al * x) * sp.cos(k * x) + k * sp.sin(k * x)
chi_e = kap * sp.cosh(kap * x) - al * sp.tanh(al * x) * sp.sinh(kap * x)
chi_o = al * sp.tanh(al * x) * sp.cosh(kap * x) - kap * sp.sinh(kap * x)
pre = hb / sp.sqrt(2 * m)
print("C4  A^dag cos = pre*psi_o:", zero(Ad(sp.cos(k * x)) - pre * psi_o, x, num),
      " A^dag sin = -pre*psi_e:", zero(Ad(sp.sin(k * x)) + pre * psi_e, x, num),
      " A^dag cosh = pre*chi_o:", zero(Ad(sp.cosh(kap * x)) - pre * chi_o, x, num),
      " A^dag sinh = -pre*chi_e:", zero(Ad(sp.sinh(kap * x)) + pre * chi_e, x, num))
print("    H_PT psi_e = E psi_e:", zero(H_PT(psi_e) - hb**2 * k**2 / (2 * m) * psi_e, x, num),
      " psi_o:", zero(H_PT(psi_o) - hb**2 * k**2 / (2 * m) * psi_o, x, num),
      " chi_e (E=-hb^2 kap^2/2m):", zero(H_PT(chi_e) + hb**2 * kap**2 / (2 * m) * chi_e, x, num),
      " chi_o:", zero(H_PT(chi_o) + hb**2 * kap**2 / (2 * m) * chi_o, x, num))
z_odd = al * sp.tanh(al * x)
z_even = al * x * sp.tanh(al * x) - 1
print("    E = 0: A^dag 1 = pre*alpha tanh:", zero(Ad(sp.Integer(1) + 0 * x) - pre * z_odd, x, num),
      " A^dag x = pre*(alpha x tanh - 1):", zero(Ad(x) - pre * z_even, x, num),
      " both solve H_PT psi = 0:", zero(H_PT(z_odd), x, num), zero(H_PT(z_even), x, num))
print("    psi_e(0) = k:", sp.simplify(psi_e.subs(x, 0)) == k, "  psi_o'(0) = alpha^2 + k^2:",
      sp.simplify(sp.diff(psi_o, x).subs(x, 0) - al**2 - k**2) == 0,
      "  chi_o'(0) = alpha^2 - kappa^2:", sp.simplify(sp.diff(chi_o, x).subs(x, 0) - al**2 + kap**2) == 0)

y1 = 1 / sp.cosh(al * x)
y2 = (x / 2 + sp.sinh(2 * al * x) / (4 * al)) / sp.cosh(al * x)
print("C5  E = -hbar^2 alpha^2/2m:  A y1 = 0:", zero(A(y1), x, num),
      "  H_PT y1 = E y1:", zero(H_PT(y1) + E0 * y1, x, num),
      "  H_PT y2 = E y2:", zero(H_PT(y2) + E0 * y2, x, num),
      "  A^dag cosh(alpha x) = 0:", zero(Ad(sp.cosh(al * x)), x, num),
      "  chi_o at kappa = alpha is 0:", zero(chi_o.subs(kap, al), x, num),
      "  (psi cosh)' = C cosh^2 integral:", zero(sp.diff(x / 2 + sp.sinh(2 * al * x) / (4 * al), x) - sp.cosh(al * x)**2, x, num))

B = sp.Symbol('Lal', positive=True) / sp.sqrt(2 * sp.pi * a)
G = sp.exp(a / 2) * sp.erf(sp.sqrt(a / 2))
print("C12 B' = -B/2a:", sp.simplify(sp.diff(B, a) + B / (2 * a)) == 0,
      " B'' = 3B/4a^2:", sp.simplify(sp.diff(B, a, 2) - 3 * B / (4 * a**2)) == 0,
      " G' = G/2 + 1/sqrt(2 pi a):", sp.simplify(sp.diff(G, a) - G / 2 - 1 / sp.sqrt(2 * sp.pi * a)) == 0,
      " G'' = G'/2 - a^(-3/2)/(2 sqrt(2 pi)):",
      sp.simplify(sp.diff(G, a, 2) - sp.diff(G, a) / 2 + a**sp.Rational(-3, 2) / (2 * sp.sqrt(2 * sp.pi))) == 0)
aJ1 = sp.sqrt(sp.pi) * sp.exp(a) * sp.erf(sp.sqrt(a)) / sp.sqrt(a)
print("    alpha J2 = alpha J1 (1 - 1/2a) + 1/a:", sp.simplify(sp.diff(aJ1, a) - aJ1 * (1 - 1 / (2 * a)) - 1 / a) == 0,
      "   a^2 (d/da)^2 (-ln(a)/2) = 1/2:", sp.simplify(a**2 * sp.diff(-sp.log(a) / 2, a, 2)) == sp.Rational(1, 2))
ser_G = sp.series(G / sp.sqrt(2 * a / sp.pi), a, 0, 4).removeO()
print("C15a G/sqrt(2a/pi) =", sp.expand(ser_G), " (expect 1 + a/3 + a^2/15 + a^3/105)")
n = sp.Symbol('n', integer=True, nonnegative=True)
print("     gap term coefficients (2^n/(n+1) - 1)/(2n+1)!!, n = 0..5:",
      [sp.Rational(2**i, i + 1) - 1 for i in range(6)], "-> zero for n = 0, 1; positive for n >= 2")
yv = sp.Symbol('y', positive=True)
Zb = yv - sp.Rational(1, 2)
bt = sp.Symbol('beta', positive=True)
yb = sp.Symbol('C0', positive=True) * bt**sp.Rational(-1, 2)
Cbox = sp.simplify(bt**2 * sp.diff(sp.log(yb - sp.Rational(1, 2)), bt, 2))
print("C16a empty box: beta^2 (ln Z)'' - 1/2 - (y-1)/(8(y-1/2)^2) =",
      sp.simplify(Cbox.subs(sp.Symbol('C0', positive=True), yv * sp.sqrt(bt)) - sp.Rational(1, 2) - (yv - 1) / (8 * (yv - sp.Rational(1, 2))**2)))
rr = sp.Symbol('r', positive=True)
exc = rr * (1 - rr) / (2 * (2 - rr)**2)
print("     r(1-r)/(2(2-r)^2) == (y-1)/(8(y-1/2)^2) at y = 1/r:", sp.simplify(exc - ((1 / rr - 1) / (8 * (1 / rr - sp.Rational(1, 2))**2))) == 0,
      "  d/dr = 0 at r =", sp.solve(sp.diff(exc, rr), rr), "  max =", exc.subs(rr, sp.Rational(2, 3)))
print("C18a sign rule: beta^2 (d/dbeta)^2 (-C beta^p) =", sp.factor(bt**2 * sp.diff(-sp.Symbol('C') * bt**sp.Symbol('p'), bt, 2)))

print("=" * 100)
print("C6  Box levels, alpha L = 20 (units hbar^2 alpha^2/2m); exactly one root per j in (j pi, (j+1) pi)")
L = 20.0
cc, kb, kj = exact_levels(L, 8)
print("    c/alpha = %.15f  (1 - c/alpha = %.3e,  2e^{-alpha L} = %.3e)" % (cc, 1 - cc, 2 * np.exp(-L)))
print("    cL - 2 = %.3f > 0 ;  min f'(k) = L - 2/c = %.4f" % (cc * L - 2, L - 2 / cc))
print("    E_j:", np.round(kj[:5]**2, 6), " empty box (j pi/alpha L)^2:", np.round((np.arange(1, 6) * np.pi / L)**2, 6))
print("    k_j L/pi:", np.round(kj[:6] * L / np.pi, 4), "  residual max:", np.max(np.abs(kj * L - np.arange(1, 9) * np.pi - 2 * np.arctan(kj / cc))))
print("    slide fraction 2 arctan(k_j/alpha)/pi:", np.round(2 * np.arctan(kj[:6]) / np.pi, 4))
for alL in [2.0, 2.4, 2.5, 3.0]:
    print("    alpha L = %.1f: c L = alpha L tanh(alpha L/2) = %.4f" % (alL, alL * np.tanh(alL / 2)))
print("    threshold alpha L tanh(alpha L/2) = 2 at alpha L = %.4f" % brentq(lambda q: q * np.tanh(q / 2) - 2, 1, 4))

print("C7  Bound state kappa_b = c tanh(kappa_b L/2)  (mpmath, 40 digits)")
mp.mp.dps = 40
for alL in [10, 20, 40]:
    Lm = mp.mpf(alL)
    cm = mp.tanh(Lm / 2)
    km = mp.mpf(1)
    for _ in range(3000):
        km = cm * mp.tanh(km * Lm / 2)
    print("    alpha L = %2d: kappa_b/alpha = %s,  1 - 4e^{-aL} = %s;  E_b/(hb^2 al^2/2m) = -%s,  1 - 8e^{-aL} = %s;  diff/(aL e^{-2aL}) = %s"
          % (alL, mp.nstr(km, 18), mp.nstr(1 - 4 * mp.e**(-Lm), 18), mp.nstr(km**2, 14), mp.nstr(1 - 8 * mp.e**(-Lm), 14),
             mp.nstr((km**2 - (1 - 8 * mp.e**(-Lm))) / (Lm * mp.e**(-2 * Lm)), 4)))
gk = lambda q: cc * np.tanh(q * L / 2) - q
print("    g(kappa) = c tanh(kappa L/2) - kappa at alpha L = 20: g'(0) = cL/2 - 1 = %.3f, root in (0, c): %s" % (cc * L / 2 - 1, brentq(gk, 0.5, cc)))
print("    odd condition kappa tanh(kappa L/2) = c: root = %.15f (= alpha, spurious)" % brentq(lambda q: q * np.tanh(q * L / 2) - cc, 0.5, 2))
print("    y1(L/2) = %.3e > 0,  y2(L/2) = %.3e > 0  (alpha L = 20)" % (1 / np.cosh(L / 2), (L / 4 + np.sinh(L) / 4) / np.cosh(L / 2)))

print("C8  Far field: psi_e(x) vs sqrt(k^2+alpha^2) cos(kx + arctan(alpha/k)), k = 0.5 alpha")
kk = 0.5
for xx in [2.0, 5.0, 10.0]:
    pe = kk * np.cos(kk * xx) - np.tanh(xx) * np.sin(kk * xx)
    po = np.tanh(xx) * np.cos(kk * xx) + kk * np.sin(kk * xx)
    th = np.arctan(1 / kk)
    print("    alpha x = %4.1f: |psi_e - far| = %.2e,  |psi_o - far| = %.2e,  2 e^{-2 alpha x} = %.2e"
          % (xx, abs(pe - np.sqrt(kk**2 + 1) * np.cos(kk * xx + th)), abs(po - np.sqrt(kk**2 + 1) * np.sin(kk * xx + th)), 2 * np.exp(-2 * xx)))
print("    removed-state density: integral_0^inf (2 alpha/pi)/(k^2+alpha^2) dk = %.15f" % quad(lambda q: 2 / np.pi / (q**2 + 1), 0, np.inf)[0])

print("=" * 100)
print("C9  Sum -> integral remainder R = sum_{j>=1} F(j) - [int_0^inf F dj - 1/2]")
mp.mp.dps = 80
print("    (a) empty box, F = e^{-omega y^2}: R exact = sqrt(pi/omega) sum_l e^{-pi^2 l^2/omega} = 2(L/lam) sum_l e^{-4 pi l^2 (L/lam)^2}")
for Llam in [0.5, 1.0, 1.5, 2.0, 3.0]:
    om = mp.pi / (4 * mp.mpf(Llam)**2)                     # omega = beta E_1 = pi lam^2 / 4L^2
    direct = mp.nsum(lambda q: mp.e**(-om * q**2), [1, mp.inf])
    Rd = direct - (mp.sqrt(mp.pi / om) / 2 - mp.mpf(1) / 2)
    Rf = 2 * Llam * mp.fsum(mp.e**(-4 * mp.pi * l**2 * mp.mpf(Llam)**2) for l in range(1, 8))
    print("       L/lam = %.1f: R direct = %s,  formula = %s" % (Llam, mp.nstr(Rd, 6), mp.nstr(Rf, 6)))


def mp_levels(Lm, cm, beta, cutoff=220):
    ks = []
    j = 1
    while True:
        q = j * mp.pi / Lm
        for _ in range(300):
            q = (j * mp.pi + 2 * mp.atan(q / cm)) / Lm
        ks.append(q)
        if beta * q**2 / 2 > cutoff:
            break
        j += 1
    return ks


print("    (b) PT box (exact c and levels, 80 digits): R and R / e^{a/2}  (e^{a/2} = bound-state weight)")
for alL in [20, 30, 60]:
    for Llam in [1.0, 2.0, 3.0]:
        Lm = mp.mpf(alL)
        cm = mp.tanh(Lm / 2)
        lam = Lm / Llam
        beta = lam**2 / (2 * mp.pi)
        ks = mp_levels(Lm, cm, beta)
        S = mp.fsum(mp.e**(-beta * q**2 / 2) for q in ks)
        R = S - (Lm / lam - mp.mpf(1) / 2 - mp.e**(beta * cm**2 / 2) * mp.erfc(cm * mp.sqrt(beta / 2)))
        Rb = 2 * Llam * mp.fsum(mp.e**(-4 * mp.pi * l**2 * mp.mpf(Llam)**2) for l in range(1, 8))
        print("       alpha L = %2d, L/lam = %.0f (a = %7.2f): R = %10s, empty-box part %10s, (R - empty)/e^{a/2} = %10s, (alpha L)^2 e^{-2 alpha L} = %s"
              % (alL, Llam, float(beta), mp.nstr(R, 4), mp.nstr(Rb, 4), mp.nstr((R - Rb) / mp.e**(beta / 2), 4), mp.nstr(Lm**2 * mp.e**(-2 * Lm), 3)))
print("    (c) |B_2m|/(2m)! <= 4/(2 pi)^{2m}:", all(abs(mp.bernoulli(2 * q)) / mp.factorial(2 * q) <= 4 / (2 * mp.pi)**(2 * q) for q in range(1, 30)),
      "  max|P_2| = 1/12:", abs(max(abs((t**2 - t + 1 / 6) / 2) for t in np.linspace(0, 1, 100001)) - 1 / 12) < 1e-12)

print("C10 Gaussian-Lorentzian integral int_0^inf e^{-s u^2}/(u^2+c^2) du = (pi/2c) e^{s c^2} erfc(c sqrt s)")
for sv, cv_ in [(0.1, 1.0), (1.0, 1.0), (3.0, 0.5), (10.0, 2.0)]:
    qv = quad(lambda u: np.exp(-sv * u**2) / (u**2 + cv_**2), 0, np.inf, epsabs=1e-14, epsrel=1e-13)[0]
    cf = np.pi / (2 * cv_) * np.exp(sv * cv_**2) * erfc(cv_ * np.sqrt(sv))
    print("    s = %5.1f, c = %.1f: quad = %.15f, closed form = %.15f, rel diff %.1e" % (sv, cv_, qv, cf, abs(qv / cf - 1)))
print("    Gaussian integral int_0^inf e^{-s u^2} du = sqrt(pi/s)/2 at s = 2.3: %.1e" % abs(quad(lambda u: np.exp(-2.3 * u**2), 0, np.inf)[0] - 0.5 * np.sqrt(np.pi / 2.3)))
print("    W_n = int_{-1}^{1} (1-t^2)^n dt vs 2^{n+1} n!/(2n+1)!!:",
      max(abs(quad(lambda t: (1 - t**2)**nn, -1, 1)[0] / (2**(nn + 1) * float(sp.factorial(nn)) / float(sp.factorial2(2 * nn + 1))) - 1) for nn in range(0, 12)))

print("C11 Z_q closed form vs direct level sum, alpha L = 100")
L = 100.0
cc, kb, kj = exact_levels(L, 4000)
for T in [0.02, 0.1, 0.5, 2.0]:
    beta = 1 / T
    direct = np.exp(beta * kb**2 / 2) + np.sum(np.exp(-beta * kj**2 / 2))
    lam = np.sqrt(2 * np.pi * beta)
    exact_c = L / lam - 0.5 + np.exp(beta * kb**2 / 2) - np.exp(beta * cc**2 / 2) * erfc(cc * np.sqrt(beta / 2))
    a_ = beta
    closed = L / lam - 0.5 + np.exp(a_ / 2) * erf(np.sqrt(a_ / 2))
    print("    k_B T/(hb^2 al^2/2m) = %4.2f: direct = %.15e, with exact c,kappa_b: rel %.1e, closed form: rel %.1e"
          % (2 * T, direct, abs(exact_c / direct - 1), abs(closed / direct - 1)))

print("C13 closed-form Cv_q vs k_B beta^2 Var(E) over the exact levels, k_B T in [0.02, 5] (hb^2 al^2/2m)")
for L in [60.0, 100.0, 1000.0]:
    T = np.linspace(0.02, 5, 400) * 0.5
    nb = int(L / np.pi * np.sqrt(2 * 60 * T.max())) + 50
    cc, kb, kj = exact_levels(L, nb)
    E = np.concatenate([[-kb**2 / 2], kj**2 / 2])
    print("    alpha L = %6g: max |diff| = %.1e k_B  (%d levels)" % (L, np.max(np.abs(cv_quantum(T, 1.0, L) - cv_level_sum(T, E))), len(E)))
print("    low T (mpmath): Cv_q vs Schottky tail (a/2)^2 Z_c e^{-a/2}, Z_c = L/lam - 1/2 - e^{a/2} erfc(sqrt(a/2)) (continuum part)")
mp.mp.dps = 50
for alL in [100, 1000]:
    for a_ in [20, 40, 80, 160]:
        am = mp.mpf(a_)
        Zf = lambda q: alL / mp.sqrt(2 * mp.pi * q) - mp.mpf(1) / 2 + mp.e**(q / 2) * mp.erf(mp.sqrt(q / 2))
        cvq = am**2 * mp.diff(lambda q: mp.log(Zf(q)), am, 2)
        Zc = alL / mp.sqrt(2 * mp.pi * am) - mp.mpf(1) / 2 - mp.e**(am / 2) * mp.erfc(mp.sqrt(am / 2))
        tail = (am / 2)**2 * Zc * mp.e**(-am / 2)
        print("       alpha L = %4d, a = %4d (k_B T = %.4f hb^2 al^2/2m): Cv_q = %s, tail = %s, ratio = %s, 1 + 2/a = %.4f"
              % (alL, a_, 2.0 / a_, mp.nstr(cvq, 5), mp.nstr(tail, 5), mp.nstr(cvq / tail, 5), 1 + 2.0 / a_))

print("C14 classical side")
aa = np.array([0.1, 1.0, 5.0, 20.0, 50.0])
series = np.zeros(len(aa))
term = 2 * aa
for l in range(1, 400):
    if l > 1:
        term = term * 2 * aa / (2 * l - 1)
    series = series + term / l
quadI = np.array([quad(lambda t: (np.exp(ai * (1 - t**2)) - 1) / (1 - t**2) if abs(t) < 1 else ai, -1, 1, epsrel=1e-13)[0] for ai in aa])
quadI_x = np.array([quad(lambda xx: np.exp(ai / np.cosh(xx)**2) - 1, -40, 40, points=[0], limit=400, epsrel=1e-13)[0] for ai in aa])
print("    alpha I: series vs t-integral max rel diff %.1e;  t-integral vs x-integral max rel diff %.1e" % (np.max(np.abs(series / quadI - 1)), np.max(np.abs(quadI_x / quadI - 1))))
J1q = np.array([quad(lambda xx: np.exp(ai / np.cosh(xx)**2) / np.cosh(xx)**2, -40, 40, points=[0], limit=400, epsrel=1e-13)[0] for ai in aa])
J2q = np.array([quad(lambda xx: np.exp(ai / np.cosh(xx)**2) / np.cosh(xx)**4, -40, 40, points=[0], limit=400, epsrel=1e-13)[0] for ai in aa])
J1c = np.sqrt(np.pi) * np.exp(aa) * erf(np.sqrt(aa)) / np.sqrt(aa)
J2c = J1c * (1 - 1 / (2 * aa)) + 1 / aa
print("    J1 closed vs x-quadrature max rel %.1e;  J2 closed vs x-quadrature max rel %.1e" % (np.max(np.abs(J1c / J1q - 1)), np.max(np.abs(J2c / J2q - 1))))
J1ser = 2 * np.array([sum((2 * ai)**q / float(sp.factorial2(2 * q + 1)) for q in range(0, 200)) for ai in aa[:3]])
print("    J1 = 2 sum_n (2a)^n/(2n+1)!!  (DLMF 7.6.2) vs closed form, a = 0.1, 1, 5: max rel %.1e" % np.max(np.abs(J1ser / J1c[:3] - 1)))


def cv_classical_quad(T, L):
    out = np.zeros(len(T))
    X = 40.0
    for i in range(len(T)):
        beta = 1 / T[i]
        w = lambda xx: np.exp(beta * (1 / np.cosh(xx)**2 - 1))
        Z = quad(w, -X, X, limit=400, points=[0])[0] + (L - 2 * X) * np.exp(-beta)
        V1 = quad(lambda xx: -1 / np.cosh(xx)**2 * w(xx), -X, X, limit=400, points=[0])[0] / Z
        V2 = quad(lambda xx: 1 / np.cosh(xx)**4 * w(xx), -X, X, limit=400, points=[0])[0] / Z
        out[i] = 0.5 + beta**2 * (V2 - V1**2)
    return out


for L in [100.0, 1000.0]:
    T = np.linspace(0.02, 5, 40) * 0.5
    print("    alpha L = %6g: Cv_cl closed form vs 1/2 + beta^2 Var(V) by quadrature, max |diff| = %.1e k_B" % (L, np.max(np.abs(cv_classical(T, 1.0, L) - cv_classical_quad(T, L)))))
for a_ in [100.0, 1000.0]:
    print("    low T: a = %6.0f: Cv_cl = %.6f,  1 + 1/a = %.6f" % (a_, cv_classical(np.array([1 / a_]), 1.0, 100.0)[0], 1 + 1 / a_))
print("    extension error of I to all x, alpha L = 100, a = 20: (4a/alpha) e^{-alpha L} = %.1e" % (4 * 20 * np.exp(-100.0)))
print("    low T: alpha I vs sqrt(pi/a) e^a (1 + 1/2a), a = 20, 50: %s" % [float(quadI[i] / (np.sqrt(np.pi / aa[i]) * np.exp(aa[i]) * (1 + 1 / (2 * aa[i])))) for i in [3, 4]])

print("C15 high T: G vs I/lambda_th and Cv_q - Cv_cl vs lambda_th/8L - (4/15) a^3/(alpha L)")
for a_ in [0.001, 0.01, 0.05]:
    Gv = np.exp(a_ / 2) * erf(np.sqrt(a_ / 2))
    Iv = quad(lambda t: (np.exp(a_ * (1 - t**2)) - 1) / (1 - t**2) if abs(t) < 1 else a_, -1, 1, epsrel=1e-14)[0] / np.sqrt(2 * np.pi * a_)
    D = Iv - Gv
    pred = np.sqrt(2 * a_ / np.pi) * (a_**2 / 45 + a_**3 / 105)
    T = np.array([1 / a_])
    dcv = cv_quantum(T, 1.0, 100.0)[0] - cv_classical(T, 1.0, 100.0)[0]
    wall = np.sqrt(2 * np.pi * a_) / 800
    print("    a = %.3f: I/lam - G = %.4e (series a^2/45 + a^3/105 terms: %.4e);  dCv = %.6e, wall term = %.6e, wall + well = %.6e"
          % (a_, D, pred, dcv, wall, wall - 4 / 15 * a_**3 / 100))
print("    well term alone (alpha L = 1e4, wall and wall x well terms small): dCv - lambda_th/8L vs -(4/15) a^3/(alpha L)")
for a_ in [0.05, 0.1]:
    T = np.array([1 / a_])
    dcv = cv_quantum(T, 1.0, 1e4)[0] - cv_classical(T, 1.0, 1e4)[0]
    wall = np.sqrt(2 * np.pi * a_) / 8e4
    cross = np.sqrt(2 * np.pi) * 0.75 * a_**1.5 / 1e8
    print("       a = %.2f: dCv - wall = %.3e, -(4/15)a^3/(alpha L) = %.3e, wall x well cross term +sqrt(2 pi)(3/4)a^1.5/(alpha L)^2 = %.1e"
          % (a_, dcv - wall, -4 / 15 * a_**3 / 1e4, cross))
aa2 = np.logspace(-2, 2, 30)
Dpos = [quad(lambda t: (np.exp(ai * (1 - t**2)) - 1) / (1 - t**2) if abs(t) < 1 else ai, -1, 1, epsrel=1e-13)[0] / np.sqrt(2 * np.pi * ai) - np.exp(ai / 2) * erf(np.sqrt(ai / 2)) for ai in aa2]
print("    I/lam - G > 0 for all 30 a in [0.01, 100]:", all(d > 0 for d in Dpos))

print("C16 empty box (L = 1): Cv_q - 1/2 vs r(1-r)/(2(2-r)^2),  r = lambda_th/L")
E1 = np.pi**2 / 2
T = E1 * np.logspace(-0.5, 3, 3000)
cb = cv_quantum_box_only(T, 1.0)
r = np.sqrt(2 * np.pi / T)
for r_max in [1.0, 2.0 / 3.0, 0.5]:
    msk = r < r_max
    print("    r < %.3f: max |diff| = %.1e" % (r_max, np.max(np.abs(cb[msk] - 0.5 - r[msk] * (1 - r[msk]) / (2 * (2 - r[msk])**2)))))
i = np.argmax(cb)
print("    max Cv_q = %.6f at k_B T/E_1 = %.3f  (prediction 1/2 + 1/16 at 9/pi = %.3f);  Cv_q > 1/2 from k_B T/E_1 = %.3f (prediction 4/pi = %.3f)"
      % (cb[i], T[i] / E1, 9 / np.pi, T[np.argmax(cb > 0.5)] / E1, 4 / np.pi))

print("C17 anomaly table (k_B T in units hb^2 al^2/2m, grid up to 1)")
for L in [30.0, 70.0, 80.0, 100.0, 400.0, 1000.0]:
    th = np.linspace(0.02, 1.0, 20000)
    cq, ccl = cv_quantum(th * 0.5, 1.0, L), cv_classical(th * 0.5, 1.0, L)
    d = cq - ccl
    i, iq, ic, im = np.argmax(d), np.argmax(cq), np.argmax(ccl), np.argmin(d)
    edges = np.where(np.diff((d > 0).astype(int)) != 0)[0]
    print("    alpha L = %6g: max dCv = %+.3f at %.3f (Cv_q %.3f, Cv_cl %.3f); window edges %s; Cv_q peak %.2f (%.3f); Cv_cl peak %.2f (%.3f); min dCv %+.2f (%.3f)"
          % (L, d[i], th[i], cq[i], ccl[i], ", ".join("%.3f" % th[e] for e in edges), cq[iq], th[iq], ccl[ic], th[ic], d[im], th[im]))
    print("        L/lambda_th at the max: %.1f" % (L / np.sqrt(2 * np.pi * 2 / th[i])))

print("C18 equal-weight (desorption) estimates")
for L in [100.0, 400.0, 1000.0]:
    a_q = brentq(lambda q: q / 2 - np.log(L / np.sqrt(2 * np.pi * q)), 0.5, 200)
    a_c = brentq(lambda q: q - 0.5 * np.log(q) - np.log(L / np.sqrt(np.pi)), 1.0, 200)
    print("    alpha L = %6g: k_B T_q = %.3f, k_B T_cl = %.3f (units hb^2 al^2/2m);  classical/quantum well weight at a = 10: e^{a/2}/(sqrt2 a) = %.1f"
          % (L, 2 / a_q, 2 / a_c, np.exp(5) / (np.sqrt(2) * 10)))
    th = np.linspace(0.02, 1.0, 20000)
    cq, ccl = cv_quantum(th * 0.5, 1.0, L), cv_classical(th * 0.5, 1.0, L)
    tq, tc = th[np.argmax(cq)], th[np.argmax(ccl)]
    sq = 0.25 * (1 / tq + 0.5)**2
    scl = 0.25 * (2 / tc - 0.5)**2
    print("        Schottky peak estimates (beta Delta)^2/4: quantum %.2f (computed %.2f), classical %.2f (computed %.2f)" % (sq, cq.max(), scl, ccl.max()))

print("C19 the associated-Legendre route (Sec. 4.7)")
tt, mu_ = sp.symbols('t mu')
Wa, Wb, Wc = sp.symbols('W0 W1 W2')
schr_t = -hb**2 / (2 * m) * al**2 * (1 - tt**2) * ((1 - tt**2) * Wc - 2 * tt * Wb) - hb**2 * al**2 / m * (1 - tt**2) * Wa - sp.Symbol('E') * Wa
leg_t = (1 - tt**2) * Wc - 2 * tt * Wb + (2 + (2 * m * sp.Symbol('E') / (hb**2 * al**2)) / (1 - tt**2)) * Wa
print("    Schrodinger in t = tanh(alpha x) equals -(hb^2 al^2/2m)(1-t^2) x Legendre(nu = 1, mu^2 = -2mE/hb^2 al^2):",
      sp.simplify(schr_t + hb**2 * al**2 / (2 * m) * (1 - tt**2) * leg_t) == 0)
ft = sp.sin(3 * tt) + tt**3
d2_direct = sp.diff(ft.subs(tt, sp.tanh(al * x)), x, 2)
d2_chain = (al**2 * (1 - tt**2) * ((1 - tt**2) * sp.diff(ft, tt, 2) - 2 * tt * sp.diff(ft, tt))).subs(tt, sp.tanh(al * x))
print("    chain rule psi'' = alpha^2 (1-t^2)[(1-t^2) w'' - 2t w']:", zero(d2_direct - d2_chain, x, {al: sp.Rational(13, 10)}))
P1m = ((1 - tt) / (1 + tt))**(mu_ / 2) * (tt + mu_) / sp.gamma(2 + mu_)
legeq = (1 - tt**2) * sp.diff(P1m, tt, 2) - 2 * tt * sp.diff(P1m, tt) + (2 - mu_**2 / (1 - tt**2)) * P1m
print("    P_1^{-mu}(t) = ((1-t)/(1+t))^{mu/2} (t+mu)/Gamma(2+mu) solves the nu = 1 equation:", sp.simplify(legeq) == 0)
zz, cc_ = sp.symbols('z c')
print("    terminating series F(2,-1;c;z) = 1 - 2z/c:",
      sp.simplify(sum(sp.rf(2, s_) * sp.rf(-1, s_) / (sp.rf(cc_, s_) * sp.factorial(s_)) * zz**s_ for s_ in range(0, 6)) - (1 - 2 * zz / cc_)) == 0)
mp.mp.dps = 30
worst = 0
for muv in [mp.mpf('0.3'), mp.mpf('1'), mp.mpf('2.5'), mp.mpc(0, '0.7'), mp.mpc(0, '3')]:
    for tv in [mp.mpf('-0.8'), mp.mpf('0.1'), mp.mpf('0.95')]:
        mine = ((1 - tv) / (1 + tv))**(muv / 2) * (tv + muv) / mp.gamma(2 + muv)
        ref = mp.legenp(1, -muv, tv, type=2)
        mine_p = ((1 + tv) / (1 - tv))**(muv / 2) * (tv - muv) / mp.gamma(2 - muv)
        ref_p = mp.legenp(1, muv, tv, type=2)
        worst = max(worst, abs(mine - ref), abs(mine_p - ref_p))
print("    closed forms of P_1^{-mu} and P_1^{+mu} vs mpmath Ferrers functions (real and imaginary mu): max |diff| = %s" % mp.nstr(worst, 3))
worst = 0
for alv, kv, xv in [(1.3, 0.6, -0.9), (1.3, 0.6, 0.4), (0.8, 2.0, 1.7)]:
    alv, kv, xv = mp.mpf(alv), mp.mpf(kv), mp.mpf(xv)
    th_ = mp.tanh(alv * xv)
    ce = kv * mp.cosh(kv * xv) - alv * th_ * mp.sinh(kv * xv)
    co = alv * th_ * mp.cosh(kv * xv) - kv * mp.sinh(kv * xv)
    g = alv * mp.gamma(2 + kv / alv)
    worst = max(worst, abs(g * mp.legenp(1, -kv / alv, th_, type=2) - (ce + co)), abs(g * mp.legenp(1, -kv / alv, -th_, type=2) - (ce - co)))
    pe = kv * mp.cos(kv * xv) - alv * th_ * mp.sin(kv * xv)
    po = alv * th_ * mp.cos(kv * xv) + kv * mp.sin(kv * xv)
    gi = alv * mp.gamma(2 + 1j * kv / alv)
    worst = max(worst, abs(gi * mp.legenp(1, -1j * kv / alv, th_, type=2) - (po + 1j * pe)))
print("    alpha Gamma(2+mu) P_1^{-mu}(+-tanh alpha x) = chi_e +- chi_o (mu = kappa/alpha);  = psi_o + i psi_e (mu = ik/alpha): max |diff| = %s" % mp.nstr(worst, 3))
tv = mp.mpf('0.37')
xv = mp.atanh(tv)            # alpha = 1
print("    mu = 0: P_1(t) = t, Q_1(t) = alpha x tanh(alpha x) - 1:", mp.nstr(abs(mp.legenp(1, 0, tv, type=2) - tv), 3),
      mp.nstr(abs(mp.legenq(1, 0, tv, type=2) - (xv * tv - 1)), 3))
print("    mu = 1: P_1^{-1}(t) = 1/(2 cosh(alpha x)):", mp.nstr(abs(mp.legenp(1, -1, tv, type=2) - 1 / (2 * mp.cosh(xv))), 3))
for muv in [mp.mpf('0.4'), mp.mpf('1'), mp.mpf('1.7'), mp.mpc(0, '0.9'), mp.mpf('0')]:
    f1 = lambda q: mp.legenp(1, -muv, q, type=2)
    f2 = lambda q: mp.legenp(1, -muv, -q, type=2)
    Wn = f1(tv) * mp.diff(f2, tv) - mp.diff(f1, tv) * f2(tv)
    Wf = 2 / (mp.gamma(muv - 1) * mp.gamma(muv + 2) * (1 - tv**2)) if muv not in (0, 1) else 0
    print("    Wronskian W{P_1^{-mu}(t), P_1^{-mu}(-t)} at mu = %s: numeric %s, DLMF 14.2.3: %s" % (mp.nstr(muv, 3), mp.nstr(Wn, 6), mp.nstr(Wf, 6)))
print("    open line: P_1^{-mu}(t) near t = -1 (should blow up unless mu = 1):",
      ", ".join("mu=%s: %s" % (mp.nstr(muv, 2), mp.nstr(mp.legenp(1, -muv, mp.mpf('-0.999999'), type=2), 4)) for muv in [mp.mpf('0.5'), mp.mpf('1'), mp.mpf('1.5')]))
print("    general lambda (Sec. 4.8): bound states E_n = -(hb^2 al^2/2m)(lambda-1-n)^2 vs finite differences (hbar = m = alpha = 1)")
Xg = np.linspace(-20, 20, 2401)[1:-1]
hg = Xg[1] - Xg[0]
for lam_ in [2.0, 2.5, 3.0]:
    Vg = -0.5 * lam_ * (lam_ - 1) / np.cosh(Xg)**2
    Hg = np.diag(1 / hg**2 + Vg) + np.diag(-0.5 / hg**2 * np.ones(len(Xg) - 1), 1) + np.diag(-0.5 / hg**2 * np.ones(len(Xg) - 1), -1)
    ev = np.linalg.eigvalsh(Hg)
    neg = ev[ev < -1e-6] / 0.5
    pred = [-(lam_ - 1 - nn)**2 for nn in range(int(np.ceil(lam_ - 1))) if lam_ - 1 - nn > 0]
    print("       lambda = %.1f: E_n/(hb^2 al^2/2m) numeric %s, predicted %s" % (lam_, np.round(neg, 4), pred))

# ------------------------------------------------------------------------------------------------
# Figure for Sec. 4: phase shift, level slide, removed-state density
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 3, figsize=(13, 3.9))
ax = axs[0]
xx = np.linspace(-12, 12, 2000)
kk = 0.5
th = np.arctan(1 / kk)
pe = (kk * np.cos(kk * xx) - np.tanh(xx) * np.sin(kk * xx)) / np.sqrt(kk**2 + 1)
ax.fill_between(xx, -1.6, -1.6 + 0.9 * (1 - 1 / np.cosh(xx)**2), color="0.88", lw=0)
ax.plot(xx, np.cos(kk * xx), color="0.55", lw=1.2, ls=":", label=r"free even wave $\cos kx$")
ax.plot(xx, pe, color="#2a78d6", lw=3.0, alpha=0.55, label=r"$\psi_e/\sqrt{k^2+\alpha^2}$ (with the well)")
ax.plot(xx, np.cos(kk * np.abs(xx) + th), color="#eb6834", lw=1.3, ls="--", zorder=3,
        label=r"$\cos(k|x|+\delta(k))$, $\delta(k)=\arctan(\alpha/k)$")
ax.set_xlabel(r"$\alpha x$")
ax.set_title(r"(a) outside the well: a free wave, shifted by $\delta(k)$" "\n" r"($k=\alpha/2$; grey: shape of $V$)", fontsize=10)
ax.set_ylim(-1.65, 2.35)
ax.legend(fontsize=7.5, loc="upper center", ncol=1, framealpha=0.95)

ax = axs[1]
L = 20.0
cc, kb, kj = exact_levels(L, 7)
nmax = 6
for nn in range(1, nmax + 1):
    En = (nn * np.pi / L)**2
    ax.plot([0, 1], [En, En], color="0.4", lw=1.4)
ax.plot([2, 3], [-kb**2, -kb**2], color="#2a78d6", lw=1.8)
for jj in range(1, nmax):
    Ej = kj[jj - 1]**2
    ax.plot([2, 3], [Ej, Ej], color="#2a78d6", lw=1.4)
    ax.annotate("", xy=(2, Ej), xytext=(1, ((jj + 1) * np.pi / L)**2), arrowprops={"arrowstyle": "->", "color": "#eb6834", "lw": 0.9})
ax.annotate("", xy=(2, -kb**2), xytext=(1, (np.pi / L)**2), arrowprops={"arrowstyle": "->", "color": "#eb6834", "lw": 0.9})
ax.text(0.5, 1.02, "empty box\n" r"$n=1,2,\dots$", ha="center", va="bottom", fontsize=8.5)
ax.text(2.5, 1.02, "well + box\n" r"$j=1,2,\dots$", ha="center", va="bottom", fontsize=8.5)
ax.text(2.5, -kb**2 + 0.05, r"bound state $E_b$", ha="center", va="bottom", fontsize=8.5, color="#2a78d6")
ax.set_xlim(-0.3, 3.3)
ax.set_ylim(-1.15, 1.32)
ax.set_xticks([])
ax.set_ylabel(r"$E\,/\,(\hbar^2\alpha^2/2m)$")
ax.set_title(r"(b) each level slides down; $n=1$ becomes $E_b$" "\n" r"($\alpha L=20$; arrows: $n=j+1\to j$)", fontsize=10)

ax = axs[2]
kgrid = np.linspace(0, 6, 600)
dens = 2 / np.pi / (kgrid**2 + 1)
ax.plot(kgrid, dens, color="#2a78d6", lw=1.8, label=r"removed density $\frac{1}{\pi}\frac{2\alpha}{k^2+\alpha^2}$ (area 1)")
for a_, ls in [(0.5, "--"), (5.0, ":")]:
    ax.plot(kgrid, dens * np.exp(-a_ * kgrid**2 / 2), color="#eb6834", lw=1.2, ls=ls,
            label=r"$\times\,e^{-\beta\hbar^2k^2/2m}$, $k_BT=%g\,\hbar^2\alpha^2/2m$" % (2 / a_))
ax.set_xlabel(r"$k/\alpha$")
ax.set_title(r"(c) the one state the well removes from the" "\n" r"continuum, and its Boltzmann weight", fontsize=10)
ax.legend(fontsize=7.5)
ax.set_ylim(0, 0.75)
for a_ in axs:
    a_.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/fig_explained_levels.pdf")
fig.savefig("figures/fig_explained_levels.png", dpi=150)
print("figure written: figures/fig_explained_levels.pdf")
for a_ in [0.5, 5.0]:
    w = quad(lambda q: 2 / np.pi / (q**2 + 1) * np.exp(-a_ / 2 * q**2), 0, np.inf)[0]
    print("    removed weight at a = %.1f: %.6f, e^{a/2} erfc(sqrt(a/2)) = %.6f" % (a_, w, np.exp(a_ / 2) * erfc(np.sqrt(a_ / 2))))
