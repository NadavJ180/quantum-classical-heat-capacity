# Option A: the exact classical heat capacity, and the physics behind it

This note explains the idea behind option A from `AUDIT.md`: computing the classical Cv directly from the classical partition function instead of reaching it through a ξ (or ħ) scan. It covers what the quantity is, where the formula comes from, what each term means physically, and when it can and cannot be used. Nothing here is implemented in `src/`; it is background for deciding whether to use it later.

---

## 1. What "the classical limit" means

A quantum system and its classical counterpart share the same Hamiltonian,

$$H = \frac{p^2}{2m} + V(x).$$

They differ in what the state space is:
- **Quantum:** discrete energy levels $E_n$, which depend on ħ.
- **Classical:** a continuum of points $(x, p)$ in phase space.

**The classical limit is ħ → 0 with everything else held fixed: m, V(x) and T.** In that limit:
- level spacings shrink relative to $k_BT$;
- the sum over levels becomes an integral over phase space;
- every quantum prediction turns into its classical counterpart (the correspondence principle).

The temperature stays where it is. That is the essential difference from "going to high temperature". The two coincide only in special cases (section 6).

The paper's ξ-scaling (V → ξ²V, T → ξ²T, Eq. S7) is a way of realizing ħ → ħ/ξ by changing things an experiment can change. It reaches the same limit as ξ → ∞, and it is what the pipeline now implements. Option A computes the destination of that limit directly.

## 2. The classical partition function

In classical statistical mechanics, the canonical partition function is an integral over phase space:

$$Z_{cl}(\beta) = \frac{1}{2\pi\hbar}\int_{-\infty}^{\infty}\!dx\int_{-\infty}^{\infty}\!dp\; e^{-\beta\left[\frac{p^2}{2m}+V(x)\right]},\qquad \beta=\frac{1}{k_BT}.$$

**Why the 1/(2πħ)?** Semiclassically, each quantum state occupies an area 2πħ of phase space (Bohr–Sommerfeld / Weyl counting). Dividing by it makes $Z_{cl}$ dimensionless and equal to the high-temperature limit of the quantum $Z=\sum_n e^{-\beta E_n}$. It is a constant factor, so it drops out of every thermodynamic derivative, including Cv.

The exponent splits into a p-part and an x-part, so the double integral factorizes. The momentum integral is a Gaussian:

$$\int_{-\infty}^{\infty} e^{-\beta p^2/2m}\,dp = \sqrt{\frac{2\pi m}{\beta}}.$$

This leaves **one** integral, over position only:

$$Z_{cl}(\beta) = \sqrt{\frac{m}{2\pi\beta\hbar^2}}\;\int_{-\infty}^{\infty} e^{-\beta V(x)}\,dx .$$

That is the "1-D integral". The whole classical thermodynamics of the system is contained in $I(\beta)=\int e^{-\beta V(x)}dx$, the Boltzmann-weighted "volume" of configuration space.

## 3. From Z to the heat capacity

The standard canonical relations are

$$U = -\frac{\partial \ln Z}{\partial\beta},\qquad C_v = \frac{\partial U}{\partial T} = k_B\,\beta^2\,\frac{\partial^2 \ln Z}{\partial\beta^2}.$$

With $\ln Z_{cl} = \text{const} - \tfrac12\ln\beta + \ln I(\beta)$:

$$U = \frac{1}{2\beta} + \langle V\rangle,\qquad \langle V^k\rangle \equiv \frac{\int V^k e^{-\beta V}dx}{\int e^{-\beta V}dx},$$

$$\boxed{\;\frac{C_v^{cl}}{k_B} = \frac12 + \beta^2\Big(\langle V^2\rangle - \langle V\rangle^2\Big)\;}$$

The averages are taken over positions with Boltzmann weight $e^{-\beta V(x)}$, i.e. over where a classical particle at temperature T actually spends its time.

**ħ has disappeared.** It entered only through the constant prefactor. This is what "classical" means: the answer depends on m, V and T, and on no quantum scale.

## 4. What each term means

- **½ (kinetic energy).** The momentum is always a free Gaussian variable, whatever V is. It contributes $\tfrac12 k_BT$ to U, so ½ k_B to Cv: equipartition for the single quadratic term $p^2/2m$. It never depends on T or on the potential.
- **β²·Var(V) (potential energy).** This is the fluctuation of the potential energy, measured in units of $(k_BT)^2$. It is the classical version of the identity $C_v = k_B\beta^2\,\mathrm{Var}(E)$ (heat capacity = energy fluctuations) that the quantum Cv formula in this project is built on. A large Cv means the energy the particle holds changes a lot as T changes, i.e. it fluctuates strongly at fixed T.
- **A strict lower bound.** A variance is never negative, so $C_v^{cl} \ge \tfrac12 k_B$ for every potential and every T. The corrected pipeline uses this bound to reject false plateaus.

## 5. Worked examples

| potential | classical Cv | why |
|---|---|---|
| HO, $\tfrac12 m\omega^2x^2$ | 1, at every T | V is quadratic: ⟨V⟩ = T/2 and Var(V) = T²/2, so β²Var = ½ |
| power law, $c\lvert x\rvert^k$ | ½ + 1/k, at every T (x⁴: ¾) | generalized equipartition: ⟨x V′⟩ = k_BT and x V′ = kV give ⟨V⟩ = T/k |
| box (V = 0 inside hard walls) | ½, at every T | no potential energy to fluctuate (the k → ∞ power law) |
| Pöschl–Teller, $V_0\tan^2x$ | 1 at low T → ½ at high T | near the bottom it is a HO; at high T the particle feels only the walls (a box) |
| your double well, b = −0.5 | ≈1.00 → 1.45 (T ≈ 0.55) → ≈0.71 | see below |

**The double well is the important case.** Its classical Cv depends on T because V is not a single power law:
- **Low T:** the particle sits at the bottom of the deep well at x = 2, which is locally harmonic, so Cv ≈ 1.
- **Intermediate T:** the particle gains access to a second region of configuration space: the shoulder and shallow well around x ≈ −0.5, about 2 energy units higher. Whether it sits "low" or "high" is a large fluctuation of V. This is the classical analogue of a two-level Schottky anomaly, with the energy gap living in configuration space instead of in the spectrum. Var(V) peaks, and so does Cv (≈ 1.45).
- **High T:** the x⁴ term dominates, so Cv → ¾ from below (0.71 at the hottest grid T).

**Consequence for the project.** A peak in Cv(T) is not by itself a quantum effect, because the classical curve of this potential has one too. A quantum enhancement means $C_v^{q}(T) > C_v^{cl}(T)$ at the same T. For the base potential the audit found that this never happens (largest difference −0.0014).

## 6. Why a "T → ∞" value is not the classical limit

The old ξ-scan effectively returned $\lim_{T'\to\infty}C_v^{q}(T')$ at every T. That is the classical Cv at high T, a single number. It equals the classical Cv at T only if $C_v^{cl}$ does not depend on T. By section 5, that happens exactly for the HO, the box and power laws, the systems where every level gap scales by one factor (homogeneous level scaling). For any potential with more than one energy scale, such as the double well or Pöschl–Teller, the two are different curves.

## 7. How the ξ-scaling relates to option A

As ξ → ∞, the paper's scaling (now in the pipeline) converges to exactly this integral. The size of the remaining quantum correction is known from the Wigner–Kirkwood expansion of the free energy:

$$F_q = F_{cl} + \frac{\hbar^2\beta}{24m}\,\langle V''(x)\rangle_{cl} + O(\hbar^4).$$

For the HO this gives $C_v^q/k_B = 1 - (\beta\hbar\omega)^2/12 + \dots$, the Taylor expansion of the Einstein formula.

Replacing ħ by ħ/ξ, the correction falls off as **1/ξ²**. That is the law the corrected pipeline uses to estimate how far each ξ step still is from the limit. In the tests the estimate matched the true distance (measured against option A) to within ~1e-5.

They give the same curve. Which to use depends on the purpose:

| | ξ-scaling (now in `src/`) | option A (phase-space integral) |
|---|---|---|
| what it computes | quantum Cv of a sequence of rescaled systems, extrapolated by plateau | the limit itself |
| accuracy | ~1e-3 at the default tolerance, with an error estimate | ~1e-12 (grid-converged) |
| cost | ~35 DVR solves per temperature sweep (minutes) | one 1-D integral per T (about a second for a 1000-point sweep) |
| truncation / plateau logic | needed | none |
| connection to the paper | direct: the scaling an experiment would perform | the ħ → 0 end point the paper's scaling approaches |

A sensible combination for later is to keep the ξ-scaling as the method that mirrors the paper, and option A as an independent check of it. It is cheap enough to run for every coefficient-sweep variant, which is exactly where each potential needs its own classical curve.

## 8. When option A applies, and when it does not

The formula $C_v/k_B = \tfrac12 + \beta^2\mathrm{Var}(V)$ assumes all of the following:

1. **Canonical equilibrium** at temperature T. This is the same assumption as everywhere in the project.
2. **One particle in one dimension, with $H = p^2/2m + V(x)$** and a constant mass.
   - With position-dependent mass or velocity-dependent forces (magnetic fields), the momentum integral no longer factorizes so simply.
   - In d dimensions the kinetic part gives d/2 and the configuration integral becomes d-dimensional.
3. **A confining potential**, so that $\int e^{-\beta V}dx$ is finite. Hard walls are fine as a limit; they simply cut the integration range (the box and Pöschl–Teller cases above).
4. **A system that has a classical counterpart.**
   - Spin systems and other purely discrete degrees of freedom do not.
   - Many identical particles need a 1/N! and, at low T, quantum statistics. This does not affect a single particle.

Numerically it is the easiest quantity in the project:
- the integrand is smooth and decays like $e^{-\beta V}$;
- a uniform grid with the trapezoid rule converges exponentially fast;
- the only care needed is a window wide enough to contain $e^{-\beta V}$ down to ~e⁻⁶⁰ at the hottest T, and a spacing fine enough to resolve the narrowest well at the coldest T.

The audit's implementation (`classical_cv_grid` in `audit_classical_limit.py`) does exactly that, with grid doubling until two grids agree to 1e-10. It reproduces the HO (1) and x⁴ (¾) exactly, and the analytic Pöschl–Teller result to 7e-13.
