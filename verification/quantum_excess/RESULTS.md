# Quantum excess check: results

Written by `quantum_excess_check.py` (units hbar = m = k_B = 1). Discussion and derivations: [SUMMARY.md](SUMMARY.md).

## A. Candidates: max(Cv_q - Cv_cl) and the windows where it is positive

### A1. Power laws |x|^n (classical Cv = 1/2 + 1/n at every T)

| n | levels | T range | max(Cv_q - Cv_cl) | at T | Cv_q / Cv_cl there | window(s) with Cv_q > Cv_cl |
|---|---|---|---|---|---|---|
| 2 | 900 | 0.05-40 | -0.0001 | 40 | 1.000 / 1.000 | none |
| 4 | 191 | 0.05-40 | -0.0005 | 40 | 0.749 / 0.750 | none |
| 6 | 121 | 0.05-40 | -0.0008 | 40 | 0.666 / 0.667 | none |
| 8 | 102 | 0.05-40 | -0.0010 | 40 | 0.624 / 0.625 | none |
| 10 | 99 | 0.05-40 | -0.0011 | 40 | 0.599 / 0.600 | none |
| 12 | 102 | 0.05-40 | -0.0011 | 40 | 0.582 / 0.583 | none |
| 13 | 107 | 0.05-40 | +0.0015 | 2.83 | 0.578 / 0.577 | 2.2-4.69 |
| 14 | 102 | 0.05-40 | +0.0044 | 2.71 | 0.576 / 0.571 | 1.86-6.84 |
| 16 | 92 | 0.05-40 | +0.0096 | 2.6 | 0.572 / 0.562 | 1.64-12.3 |
| 20 | 80 | 0.05-40 | +0.0180 | 2.6 | 0.568 / 0.550 | 1.44-28.6 |
| 30 | 55 | 0.05-40 | +0.0311 | 2.6 | 0.564 / 0.533 | 1.38-40 |
| 50 | 25 | 0.05-23.2 | +0.0430 | 2.83 | 0.563 / 0.520 | 1.38-23.2 |

### A2. Quartic double wells V = A [x^4 - 4/3 (1 + r) x^3 + 2 r x^2]

Minima at x = 0 and x = 1, barrier at x = r. For r > 1/2 the well at 0 is the deeper AND the stiffer one (omega_0^2 / omega_1^2 = r / (1 - r)): the dimple-like topology. A sets the depth against hbar^2/m (x in units of the well separation). Entry: max(Cv_q - Cv_cl) over T, with ZPE/gap = (hbar omega_0 / 2) / (V(1) - V(0)) in brackets.

| r \ A | 3 | 10 | 30 | 100 | 300 | 1000 | 3000 |
|---|---|---|---|---|---|---|---|
| 0.60 | -0.0325 (6.7) | -0.0186 (3.7) | -0.0129 (2.1) | -0.0086 (1.2) | -0.0061 (0.67) | -0.0024 (0.37) | -0.0007 (0.21) |
| 0.75 | -0.0277 (3) | -0.0172 (1.6) | -0.0121 (0.95) | -0.0086 (0.52) | -0.0023 (0.3) | -0.0006 (0.16) | -0.0002 (0.095) |
| 0.90 | -0.0244 (2.1) | -0.0153 (1.1) | -0.0119 (0.65) | -0.0043 (0.36) | -0.0011 (0.21) | -0.0003 (0.11) | -0.0001 (0.065) |
| 0.97 | -0.0230 (1.8) | -0.0150 (0.99) | -0.0123 (0.57) | -0.0034 (0.31) | -0.0009 (0.18) | -0.0002 (0.099) | -0.0001 (0.057) |
| 0.99 | -0.0232 (1.8) | -0.0153 (0.96) | -0.0117 (0.56) | -0.0032 (0.3) | -0.0009 (0.18) | -0.0002 (0.096) | -0.0001 (0.056) |

Largest value over all 35 shapes: -0.0001.

### A3. Narrow dimples in wide traps

ZPE = E_0 - V_min (V_min = -D, the trap is 0 at x = 0). lam is the Poschl-Teller parameter of the sech^2 dimples, lam (lam + 1) = 2 D s^2; their exact ground energy -lam^2 / (2 s^2) is shown next to the DVR value.

| potential | bound | E_0 (DVR) | ZPE / D | max(Cv_q - Cv_cl) | at T | Cv_q / Cv_cl there | window(s) with Cv_q > Cv_cl |
|---|---|---|---|---|---|---|---|
| Gauss D=60 s=0.1 in x^2/2 | 1 | -30.5470 | 0.49 | +0.104 | 7.22 | 2.385 / 2.281 | 6.63-7.63 |
| Gauss D=30 s=0.1 in (x/5)^20 | 1 | -11.5892 | 0.61 | +0.562 | 3.02 | 2.188 / 1.627 | 2.34-3.57 |
| Gauss D=30 s=0.1 in (x/10)^20 | 1 | -11.5892 | 0.61 | +1.578 | 2.77 | 3.220 / 1.642 | 1.98-3.47 |
| Gauss D=30 s=0.1 in (x/20)^20 | 1 | -11.5892 | 0.61 | +2.863 | 2.55 | 4.476 / 1.613 | 1.67-3.19 |
| Gauss D=30 s=0.1 in (x/10)^6 | 1 | -11.5892 | 0.61 | +1.876 | 2.77 | 3.571 / 1.695 | 1.92-3.47 |
| Gauss D=80 s=0.1 in (x/10)^20 | 2 | -44.6293 | 0.44 | +0.846 | 8.31 | 3.300 / 2.455 | 6.63-9.83 |
| Gauss D=60 s=0.15 in (x/10)^20 | 2 | -38.3651 | 0.36 | +0.031 | 6.82 | 2.718 / 2.686 | 6.63-7.02 |
| sech2 lam=0.4 (D=28) s=0.1 in (x/10)^6 | 1 | -8.0000 (exact -8.0000) | 0.71 | +2.277 | 2.21 | 3.590 / 1.313 | 1.41-2.93 |
| sech2 lam=0.5 (D=37.5) s=0.1 in (x/10)^6 | 1 | -12.5000 (exact -12.5000) | 0.67 | +2.427 | 3.19 | 3.975 / 1.548 | 2.03-4 |
| sech2 lam=0.7 (D=59.5) s=0.1 in (x/10)^6 | 1 | -24.5000 (exact -24.5000) | 0.59 | +2.422 | 5.45 | 4.418 / 1.996 | 3.67-6.63 |
| sech2 lam=1.0 (D=100) s=0.1 in (x/10)^6 | 1 | -50.0000 (exact -50.0000) | 0.50 | +2.088 | 9.56 | 4.552 / 2.464 | 7.02-11.6 |

## B. The criterion Cv_q - Cv_cl = d(U_q - U_cl)/dT and its consequences

### B1. Best candidate, sech2 lam=0.5 (D=37.5) s=0.1 in (x/10)^6

- ZPE = E_0 - V_min = 25.0000; Delta U = U_q - U_cl at the coldest T (0.3): 24.6988 (-> ZPE as T -> 0; here ZPE - Delta U = 0.301, the classical thermal energy of the dimple the quantum one does not yet have).
- Finite-difference d(Delta U)/dT vs Cv_q - Cv_cl: same sign at 100.0% of the 396 interior temperatures; largest difference 0.00034 of max|Cv_q - Cv_cl| (discretisation of the T grid).
- Windows where Delta U rises: 2.03-4.03; windows where Cv_q > Cv_cl: 2.03-4.03.
- Area: int (Cv_q - Cv_cl) dT over the excess window(s) = +3.101; over the whole range 0.3-12 = -23.079 (= Delta U(T_max) - Delta U(T_min) = -23.079).

### B1b. What the window is worth in an Otto cycle (homogeneous stroke q, SUMMARY.md section 4)

W = (q - 1) int_{T_c}^{T_h/q} Cv dT for both the quantum and the classical working medium (same q, so the same efficiency 1 - 1/q): W_q / W_cl = int Cv_q / int Cv_cl over [T_c, T_h/q].

| [T_c, T_h/q] | int Cv_q dT | int Cv_cl dT | W_q / W_cl |
|---|---|---|---|
| 2.03-4.03 | 6.534 | 3.433 | 1.903 |
| 2.76-3.64 | 3.415 | 1.449 | 2.356 |
| 1-3 | 2.715 | 2.190 | 1.240 |
| 0.5-10 | 17.445 | 38.789 | 0.450 |

### B2. Sum rule int_0^inf (Cv_q - Cv_cl) dT = -(E_0 - V_min)

- HO (hbar omega = 1): int_0^inf (C_Einstein - 1) dT = -0.500000, against -ZPE = -0.5.
- The premise U_q - U_cl -> 0 at high T, on x^4 (Wigner-Kirkwood: Delta U = (hbar^2/24m) d/dbeta [beta^2 <V''>] = 0.75 (Gamma(3/4)/Gamma(1/4)) T^-1/2 = 0.2535 T^-1/2):

| T | Delta U (DVR - exact classical) | Wigner-Kirkwood | ratio |
|---|---|---|---|
| 1 | 0.24219 | 0.25349 | 0.955 |
| 1.89 | 0.17918 | 0.18444 | 0.972 |
| 4.05 | 0.12451 | 0.12593 | 0.989 |
| 8.69 | 0.08565 | 0.08598 | 0.996 |
| 18.6 | 0.05863 | 0.05870 | 0.999 |
| 40 | 0.04006 | 0.04008 | 1.000 |

### B3. Hard walls: the box, E_n = n^2 (units of E_1), Cv_cl = 1/2

- Peak Cv_q = 0.5625 at k_B T = 2.85 E_1 (Cv_q > 1/2 for every T above 1.28 E_1).
- Leading O(hbar) prediction Cv_q - 1/2 = 1/(8 Z_cl), Z_cl = (1/2) sqrt(pi T / E_1):

| T / E_1 | Cv_q - 1/2 | 1/(8 Z_cl) | ratio |
|---|---|---|---|
| 9.91 | 4.268e-02 | 4.481e-02 | 0.9523 |
| 101 | 1.400e-02 | 1.405e-02 | 0.9965 |
| 991 | 4.480e-03 | 4.481e-03 | 0.9997 |
| 1.01e+04 | 1.405e-03 | 1.405e-03 | 1.0000 |
| 9.91e+04 | 4.481e-04 | 4.481e-04 | 1.0000 |

### B4. Two-region model of the dimples: Z_q ~ sum over bound states exp(-E_n/T) + G(T)

G(T) = sqrt(T/2 pi) int exp(-V_trap/T) dx is the trap's own (classical) partition function, energies from the trap bottom. T_q = quantum escape temperature, exp(E_b/T_q) = G(T_q) with E_b = -E_0. Design rules: SUMMARY.md section 2.3 (c), R2 (ZPE / T_q, the classical particle is still trapped) and R3 (ln G(T_q), a large reservoir).

| potential | E_b | ZPE | T_q | ZPE / T_q | ln G(T_q) | window (DVR) | model max|Cv - Cv_DVR| in window |
|---|---|---|---|---|---|---|---|
| Gauss D=60 s=0.1 in x^2/2 | 30.55 | 29.45 | 12.2 | 2.4 | 2.50 | 6.63-7.63 | 0.059 |
| Gauss D=30 s=0.1 in (x/5)^20 | 11.59 | 18.41 | 5.14 | 3.6 | 2.26 | 2.34-3.57 | 0.108 |
| Gauss D=30 s=0.1 in (x/10)^20 | 11.59 | 18.41 | 4.1 | 4.5 | 2.83 | 1.98-3.47 | 0.074 |
| Gauss D=30 s=0.1 in (x/20)^20 | 11.59 | 18.41 | 3.39 | 5.4 | 3.42 | 1.67-3.19 | 0.049 |
| Gauss D=30 s=0.1 in (x/10)^6 | 11.59 | 18.41 | 3.97 | 4.6 | 2.92 | 1.92-3.47 | 0.068 |
| Gauss D=80 s=0.1 in (x/10)^20 | 44.63 | 35.37 | 12.9 | 2.7 | 3.46 | 6.63-9.83 | 0.096 |
| Gauss D=60 s=0.15 in (x/10)^20 | 38.37 | 21.63 | 11.3 | 1.9 | 3.39 | 6.63-7.02 | 0.081 |
| sech2 lam=0.4 (D=28) s=0.1 in (x/10)^6 | 8.00 | 20.00 | 2.94 | 6.8 | 2.72 | 1.41-2.93 | 0.083 |
| sech2 lam=0.5 (D=37.5) s=0.1 in (x/10)^6 | 12.50 | 25.00 | 4.22 | 5.9 | 2.96 | 2.03-4 | 0.068 |
| sech2 lam=0.7 (D=59.5) s=0.1 in (x/10)^6 | 24.50 | 35.00 | 7.35 | 4.8 | 3.33 | 3.67-6.63 | 0.046 |
| sech2 lam=1.0 (D=100) s=0.1 in (x/10)^6 | 50.00 | 50.00 | 13.4 | 3.7 | 3.73 | 7.02-11.6 | 0.019 |

## C. Homogeneous eigen-energy scaling (needed for Q = int Cv dT in an Otto cycle)

Gap ratios q_n = (E_n - E_0)[V_2] / (E_n - E_0)[V_1], n = 1..39. Homogeneous means every q_n is the same. Spread = max q_n / min q_n - 1. Strokes: amplitude V_2 = 2 V_1; dilation V_2(x) = V_1(x / 1.5); scale-invariant V_2(x) = V_1(x / 1.5) / 1.5^2 (expected q = 1/1.5^2 = 0.4444).

| potential | stroke | q_1 | min q_n | max q_n | spread | homogeneous? |
|---|---|---|---|---|---|---|
| HO x^2/2 | amplitude x2 | 1.414214 | 1.414214 | 1.414214 | 4.3e-13 | yes |
| HO x^2/2 | dilation x1.5 | 0.666667 | 0.666667 | 0.666667 | 2.0e-13 | yes |
| HO x^2/2 | scale-invariant x1.5 | 0.444444 | 0.444444 | 0.444444 | 9.0e-13 | yes |
| |x|^20 | amplitude x2 | 1.065041 | 1.065041 | 1.065041 | 2.9e-12 | yes |
| |x|^20 | dilation x1.5 | 0.478447 | 0.478447 | 0.478447 | 7.7e-12 | yes |
| |x|^20 | scale-invariant x1.5 | 0.444444 | 0.444444 | 0.444444 | 4.7e-12 | yes |
| double well 0.25x^4 - 0.5x^3 - 0.5x^2 | amplitude x2 | 1.777505 | 1.241566 | 1.777505 | 4.3e-01 | NO |
| double well 0.25x^4 - 0.5x^3 - 0.5x^2 | dilation x1.5 | 0.863547 | 0.572579 | 0.863547 | 5.1e-01 | NO |
| double well 0.25x^4 - 0.5x^3 - 0.5x^2 | scale-invariant x1.5 | 0.444444 | 0.444444 | 0.444444 | 1.1e-12 | yes |
| Gauss D=30 s=0.1 in (x/10)^20 | amplitude x2 | 2.628861 | 1.731642 | 2.629172 | 5.2e-01 | NO |
| Gauss D=30 s=0.1 in (x/10)^20 | dilation x1.5 | 1.365027 | 0.856935 | 1.366349 | 5.9e-01 | NO |
| Gauss D=30 s=0.1 in (x/10)^20 | scale-invariant x1.5 | 0.444444 | 0.444444 | 0.444444 | 7.4e-13 | yes |
| sech2 lam=0.5 s=0.1 in (x/10)^6 | amplitude x2 | 2.698199 | 2.008870 | 2.698199 | 3.4e-01 | NO |
| sech2 lam=0.5 s=0.1 in (x/10)^6 | dilation x1.5 | 1.408252 | 1.013617 | 1.408305 | 3.9e-01 | NO |
| sech2 lam=0.5 s=0.1 in (x/10)^6 | scale-invariant x1.5 | 0.444444 | 0.444444 | 0.444444 | 5.6e-13 | yes |

Scale-invariant stroke on the best candidate, whole curves: max |Cv_q[V_2](T / 1.5^2) - Cv_q[V_1](T)| = 5.7e-13 and max |Cv_cl[V_2](T / 1.5^2) - Cv_cl[V_1](T)| = 2.4e-14 over T = 0.5..8 (40 temperatures).

Figures: [best_candidate_cv.png](figures/best_candidate_cv.png), [delta_u_criterion.png](figures/delta_u_criterion.png), [homogeneity_gap_ratios.png](figures/homogeneity_gap_ratios.png). Total run time 5.6 min.
