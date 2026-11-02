# PT (lambda = 2) in a box -- verification results

Written by `pt_box_check.py`. Units: hbar = m = k_B = alpha = 1, energy unit eps = 1/2;
theta = k_B T / eps, Lambda = alpha L. All quoted in PT_Box_Cv.tex.

## 1. Exact levels (Lambda = 60)
- quantization residual max |k L - j pi - 2 arctan(k/tau)| = 1.4e-12
- bound state E_b/eps = -1.000000000000000   (L -> infinity: -1)
- first box levels E_j/eps vs the empty box (j pi / L)^2 / 2 / eps:

| j | PT + box | empty box | shift |
|---|---|---|---|
| 1 | 2.93369203e-03 | 2.74155678e-03 | +1.921e-04 |
| 2 | 1.17324157e-02 | 1.09662271e-02 | +7.662e-04 |
| 3 | 2.63892433e-02 | 2.46740110e-02 | +1.715e-03 |
| 4 | 4.68930453e-02 | 4.38649084e-02 | +3.028e-03 |
| 5 | 7.32290633e-02 | 6.85389195e-02 | +4.690e-03 |

## 2. Exact levels vs DVR (Lambda = 60, 120 levels)

Absolute errors |E_DVR - E_exact| / eps.

| alpha dx | repo sinc-DVR: bound | repo sinc-DVR: levels 1-20 | repo sinc-DVR: levels 1-119 | sine-DVR: all 120 |
|---|---|---|---|---|
| 0.2 | 3.6e-15 | 1.9e-03 | 6.4e-02 | 2.3e-12 |
| 0.1 | 3.8e-14 | 9.4e-04 | 3.2e-02 | 1.6e-13 |
| 0.05 | 1.2e-13 | 4.7e-04 | 1.6e-02 | 4.8e-13 |
| 0.025 | 7.2e-13 | 2.3e-04 | 7.9e-03 | 1.8e-12 |

## 3. Integer lambda > 2 (asymptotic quantization) vs sine-DVR (Lambda = 40, 80 levels)

| lambda | bound states E/eps | max |E_DVR - E_asym| / eps, 80 levels |
|---|---|---|
| 2 | -1 | 2.4e-12 |
| 3 | -4, -1 | 2.4e-12 |
| 4 | -9, -4, -1 | 2.9e-12 |
| 6 | -25, -16, -9, -4, -1 | 2.8e-12 |

## 4. Quantum Cv: Var(E) over the exact levels vs the closed-form Z_q

| Lambda | max rel. diff ln Z | max abs diff Cv_q | levels used |
|---|---|---|---|
| 60 | 2.2e-16 | 7.3e-13 | 307 |
| 100 | 2.2e-16 | 4.9e-13 | 498 |
| 1000 | 2.2e-16 | 1.1e-12 | 4795 |

## 5. Classical closed forms vs quadrature / finite differences

- I(a): power series (2F2) vs quadrature, max rel. diff over a in [0.05, 50]: 1.3e-14
- J1(a) erf form vs quadrature: 1.8e-15;  J2(a): 1.6e-15
- Cv_cl closed form vs finite-difference beta^2 d^2 lnZ_cl (Lambda = 100): max abs diff 5.1e-06 (FD truncation)

## 6. xi formula along the integer ladder xi_n = sqrt(n(n-1)/2)  (Lambda = 100)

| theta | Cv_cl | n=2 (xi=1) | n=3 (xi=1.73) | n=5 (xi=3.16) | n=10 (xi=6.71) | n=20 (xi=13.8) | n=50 (xi=35) | n=100 (xi=70.4) | n=300 (xi=212) | n=1000 (xi=707) | n=3000 (xi=2.12e+03) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 1.822012 | +3.1e-01 | -3.7e-01 | -5.1e-01 | -2.4e-01 | -7.2e-02 | -1.2e-02 | -3.2e-03 | -4.3e-04 | -6.1e-05 | -1.4e-05 |
| 0.4 | 5.507586 | -3.0e+00 | -1.7e+00 | -7.0e-01 | -1.8e-01 | -4.4e-02 | -6.6e-03 | -1.5e-03 | -1.1e-04 | +6.4e-06 | +6.0e-06 |
| 1.0 | 0.735324 | -8.5e-02 | -3.8e-02 | -1.2e-02 | -2.1e-03 | -1.3e-04 | +1.6e-04 | +1.1e-04 | +4.4e-05 | +1.4e-05 | +4.7e-06 |

(entries: Cv^(xi_n) - Cv_cl; n = 2 is the physical system itself, xi = 1)

## 7. Cv_q vs Cv_cl (exact levels for Cv_q, closed form for Cv_cl)

| Lambda | max (Cv_q - Cv_cl) | at theta | Cv_q there | Cv_cl there | windows Cv_q > Cv_cl (theta) | Cv_q peak (theta) | Cv_cl peak (theta) | min (Cv_q - Cv_cl) (theta) | wall term lth/8L at max |
|---|---|---|---|---|---|---|---|---|---|
| 30 | +0.007 | 3.668 | 0.524 | 0.517 | [2.11, 5] | 2.00 (0.392) | 3.74 (0.410) | -1.75 (0.419) | 7.7e-03 |
| 100 | +0.343 | 0.209 | 2.400 | 2.057 | [0.177, 0.238]; [2.69, 5] | 3.56 (0.285) | 6.21 (0.344) | -3.18 (0.367) | 9.7e-03 |
| 400 | +2.590 | 0.183 | 5.048 | 2.458 | [0.13, 0.223] | 6.05 (0.216) | 10.02 (0.285) | -6.20 (0.298) | 2.6e-03 |
| 1000 | +4.801 | 0.167 | 7.446 | 2.645 | [0.112, 0.202] | 8.13 (0.185) | 13.13 (0.252) | -9.10 (0.263) | 1.1e-03 |

(the upper window edge 5 is the end of the theta grid; above the dimple window the
 remaining small excess is the wall term lambda_th/(8L), see section 8)

## 8. The empty box (no dimple)
- Cv_q max = 0.5625 at T/E_g = 2.86; Cv_q > 1/2 for all T/E_g > 1.26
- at T/E_g = 100: Cv_q - 1/2 = 0.0140; wall term lambda_th/(8L) = 0.0141

