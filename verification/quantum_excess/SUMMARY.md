# When is the quantum heat capacity above the classical one?

Findings, the general criterion, candidate potentials, the homogeneous-scaling requirement of the Otto-cycle relations, and a verification plan. Written 2026-10-13, branch `correct-classical`. Units: ħ = m = k_B = 1 unless written out.

**Status.** Every number below comes from [`quantum_excess_check.py`](quantum_excess_check.py) (output in [`RESULTS.md`](RESULTS.md), figures in [`figures/`](figures/)):
- **Quantum Cv:** the pipeline's own Colbert–Miller DVR on grids chosen by the script, each spectrum confirmed on a 25% wider, 1.4× finer grid.
- **Classical Cv:** the exact ½ + β²Var(V), which the ξ-scan converges to ([CHANGES_SUMMARY.md, section F](../../CHANGES_SUMMARY.md)).
- **Through the pipeline: partly.** After the `DVR_DX_MAX` fix (section 5.2) the pipeline reproduces the candidate's quantum side. Its classical ξ-scan in the excess window is still too expensive to finish (section 5.3), so the classical numbers here come from the exact integral.

---

## 0. Summary

1. **A Schottky-like bump in the classical Cv is not a quantum effect.** It appears whenever heating opens a second region of configuration space (section 1).
2. **Exact criterion.** Cv_q(T) > Cv_cl(T) exactly where **ΔU = U_q − U_cl increases with T**. ΔU starts at the zero-point energy and, for smooth V, ends at 0, so ∫₀^∞ (Cv_q − Cv_cl) dT = −(E₀ − V_min). An excess is always a *window*, and it is paid for by a larger deficit elsewhere (section 2).
3. **What works: a narrow, deep trap with one bound state, inside a wide, soft reservoir.** The trap's zero-point energy lowers the quantum escape energy far below the classical one, so the quantum particle escapes, and gains the reservoir's entropy, at a lower T. Best case: V = (x/10)⁶ − 37.5·sech²(x/0.1), with Cv_q = 3.98 against Cv_cl = 1.55 at k_BT ≈ 3.2, and quantum above classical for k_BT ≈ 2.0–4.0 (section 3).
4. **Weak: steep walls.** |x|ⁿ with n ≥ 13, and the box, exceed classical by at most 0.04–0.06.
5. **None: the quartic family a x⁴ + b x³ + c x² (+ d x).** 35 shapes were checked and every one stays below classical.
6. **Homogeneous scaling (needed for Q = ∫Cv dT).** The candidates scale homogeneously **only under scale-invariant strokes**, V → γ⁻²V(x/γ). That means changing the trap and dimple width and depth together, (D, σ, L) → (D/γ², γσ, γ^{4/3}L). Under such a stroke the excess window just slides to T/γ². Amplitude-only or width-only strokes are not homogeneous for these potentials (section 4). With scale-invariant strokes, an Otto cycle confined to the window gives **W_q/W_cl ≈ 1.9–2.4 at the same efficiency**. A wide cycle (0.5–10) gives 0.45.
7. **The pipeline and the candidate.** The DVR used to pick dx from the wavefunction momentum alone, which under-samples the narrow dimple. The fix, `DVR_DX_MAX`, is now implemented and bit-identical on existing potentials (section 5.2).
   - **Quantum side:** with the fix, the pipeline reproduces the candidate exactly.
   - **Classical side:** the ξ-scan in the excess window is still too expensive for a dense uniform DVR (section 5.3).

---

## 1. Why the classical Cv can have a Schottky-like bump

Classically, Cv_cl = ½ + β²Var(V). Split x at a barrier into basins L and R, with populations p and 1 − p. The law of total variance then gives an exact identity:

  Cv_cl = ½ + β²[p·Var_L(V) + (1−p)·Var_R(V)] + β²·p(1−p)·(⟨V⟩_R − ⟨V⟩_L)²

The last term is the two-level Schottky formula (βΔ)²p(1−p), with the gap Δ equal to the difference of the basins' mean potential energies. Each basin's phase-space volume plays the role of the degeneracy. A Schottky peak needs a gap and a saturating weight ratio, not ħ.
- **Asymmetric double wells** (0.25x⁴ + b x³ − 0.5x²): b = −1.6 peaks at 1.415 (T ≈ 16), and b = −2.0 at 1.405 (T ≈ 31). These are exactly the dashed curves of `quick_sweep_res1.png`.
- **The symmetric well** (b = 0) has no between-basin term. Its bump (1.14 at T ≈ 0.05) comes from the barrier top becoming reachable.

**What is quantum** is the low-T side: Cv_q → 0, while Cv_cl ≥ ½. The ξ-scan also gives an operational test. A bump that survives ξ → ∞ at the same T is set by the geometry of V and is classical. A quantum bump is set by level spacings and moves as ħ/ξ → 0.

---

## 2. The general criterion

### 2.1 Exact, equivalent forms

For H = p²/2m + V(x), write ΔX = X_q − X_cl at the same T. Every thermodynamic quantity is a linear functional of ln Z, so the standard identities hold for the differences too:

  ΔCv = dΔU/dT = T·dΔS/dT = −T·d²ΔF/dT² = k_B β²·∂²_β ln(Z_q/Z_cl) = k_B β²·[Var_q(H) − Var_cl(H)]

Hence **Cv_q(T) > Cv_cl(T)** ⇔ ΔU rises with T ⇔ ΔS rises with T ⇔ ΔF is concave in T ⇔ ln(Z_q/Z_cl) is convex in β.

### 2.2 Boundary values: why an excess is the exception

- **Z_q ≤ Z_cl at every T (Golden–Thompson).** Tr e^(A+B) ≤ Tr(e^A e^B) with A = −βp̂²/2m and B = −βV(x̂). The right-hand side is exactly Z_cl (CHANGES_SUMMARY F, Step 7). So ΔF ≥ 0 always.
- **T → 0:** ΔU → E₀ − V_min = ZPE > 0. The quantum particle keeps its zero-point energy; the classical one rests at the bottom.
- **T → ∞, smooth V:** ΔU → 0. By Wigner–Kirkwood, ΔU = (ħ²/24m)·∂_β[β²⟨V″⟩_cl] ∝ T^(−2/k) for a |x|^k tail. Checked on x⁴ in [RESULTS B2](RESULTS.md).
- **Sum rule.** It follows that

  **∫₀^∞ (Cv_q − Cv_cl) dT = −(E₀ − V_min)**.

  Exact for the HO: ∫(C_Einstein − 1)dT = −ħω/2, reproduced numerically in RESULTS B2.
- **The generic case:** ΔU falls monotonically from the ZPE to 0, so Cv_q < Cv_cl at every T. This holds for the HO, for |x|ⁿ with n ≤ 12, and for every quartic tested. An excess window needs ΔU to *rise*: in that window the quantum system takes up energy faster than the classical one, which means it reaches new, entropy-rich states at a lower T. The area under the window is more than repaid elsewhere, because the total area is −ZPE.
- **Hard walls break the premise ΔU(∞) = 0.** For the box, ΔU grows like √T. That is why the box stays above classical at every high T (2.3 b).

See [`figures/delta_u_criterion.png`](figures/delta_u_criterion.png): on the best candidate ΔCv > 0 exactly where ΔU rises.

### 2.3 Three regimes and what each predicts

**(a) High T, smooth V (Wigner–Kirkwood, O(ħ²)).** With ΔF ≈ ħ²⟨V″⟩_T/(24 m T):

  ΔCv ≈ −(ħ²/24m)·T·∂²_T[⟨V″⟩_T / T]

- **For every power law |x|^k this is negative**, ∝ −T^(−(1+2/k)) (FINDINGS.md; audit/quick_scan/MERGE_CHECK.txt, ratio 0.998–1.000).
- **It turns positive** only where ⟨V″⟩_T/T has a concave stretch, i.e. where the mean curvature the particle samples changes abruptly with T. A potential with a single length scale cannot produce a high-T excess at this order.

**(b) Hard or very steep walls (O(ħ), the Weyl boundary term).** Each Dirichlet wall removes ¼ of a state, so Z_q ≈ Z_cl − n_w/4. Write Z_cl ∝ T^α, where α is the high-T classical Cv. Then ln Z_q ≈ ln Z_cl + δ/Z_cl with δ = −n_w/4, and

  ΔCv ≈ α(1−α)·(n_w/4)/Z_cl,

which is positive iff 0 < α < 1, i.e. for confinement steeper than harmonic.
- **Box** (n_w = 2, α = ½): ΔCv ≈ 1/(8Z_cl). Confirmed in RESULTS B3, where the ratio tends to 1. The peak is 0.5625 at T ≈ 2.9 E₁.
- **Half-oscillator** (α = 1): no effect at this order, and its exact Cv is an Einstein function, below 1 everywhere.
- **Smooth |x|ⁿ** behaves like a hard wall only while the thermal wavelength exceeds the wall's softness length. In practice that needs n ≥ 13 (RESULTS A1). The excess is at most +0.043 (n = 50) and +0.0625 (box).

**(c) Deep quantum region: a bound state plus a reservoir (non-perturbative).** This is the mechanism that gives a large excess. Model the system as a trap of depth D with a single bound state at −E_b (E_b = D − ZPE), next to a wide reservoir with partition function G(T) (energies from the reservoir bottom):

  Z_q ≈ e^(βE_b) + G(T),  Z_cl ≈ (k_BT/ħω_d)·e^(βD) + G(T)

- **Escape temperatures.** The quantum particle escapes at T_q, where G(T_q) = e^(E_b/k_BT_q). The classical particle escapes only at T_cl, where G(T_cl) = (k_BT_cl/ħω_d)·e^(D/k_BT_cl).
- **The cold flank of the escape peaks.** The ratio of the two escape (Schottky) terms is (E_b/D)²·(k_BT/ħω_d)·e^(ZPE/k_BT). This is exponentially in favour of the quantum system when ZPE ≫ k_BT.
- **The quantum peak height.** It is ≈ (ln G)²/4, against Cv_cl ≈ 1 for the classical particle still in the trap.

This gives four design rules:
- **R1.** One bound state, with ZPE ≳ D/2, so that E_b ≪ D. For sech²: λ ≤ 1, where λ(λ+1) = 2mDσ²/ħ² and the ZPE fraction is 1/(λ+1).
- **R2.** ZPE/k_BT_q ≫ 1: the classical particle is still trapped when the quantum one escapes. Empirically ≳ 3.7 for a large excess; below ≈ 2.5 the excess is marginal (section 3.2).
- **R3.** ln G(T_q) ≳ 2: the reservoir has ≳ 7 thermally accessible states per bound state, so (ln G)²/4 beats ≈ 1. Empirically ≳ 2.8 for a large excess.
- **R4** (implied). The trap is narrower than the thermal wavelength (quantum) while the reservoir is much wider than it (classical): **σ ≪ λ_T ≪ L**. Two very different length scales are essential.

The model is tested against the DVR in [RESULTS B4](RESULTS.md).

**Why the quartic family fails** (RESULTS A2). A quartic has one length scale. When ZPE/gap ≳ 1 (R1), the second well is just as quantum, so R3 fails. When the second well is classical, ZPE/gap ≪ 1, so R1 fails. Your asymmetric wells are in the second case: at b = −1.6, ħω/2 ≈ 2.5 against a gap of 56.

---

## 3. Candidates

### 3.1 Recommended candidate

**V(x) = (x/10)⁶ − 37.5·sech²(x/0.1)**, a sech² dimple with λ = 0.5 in a soft x⁶ trap.
- **Spectrum check.** The DVR gives E₀ = −12.5000, the exact Pöschl–Teller value. The ZPE is 25, ⅔ of the depth.
- **Result.** Cv_q > Cv_cl for k_BT ≈ 2.0–4.0, with the largest excess **+2.43 at k_BT = 3.2 (Cv_q = 3.98 against Cv_cl = 1.55)** ([figure](figures/best_candidate_cv.png)).
- **What happens.** The quantum particle escapes the dimple around T ≈ 3 and the classical one only around T ≈ 5. The classical peak is then taller (9.1), because the classical particle absorbs the full depth D. The total area is −ZPE (RESULTS B1).

Why sech² rather than a Gaussian:
- its bound states are known exactly, which gives a built-in DVR check;
- λ alone sets the number of bound states and the ZPE fraction;
- λ is invariant under the scale-invariant stroke (section 4).

It is in `src/config.py` as an inactive block.

**How it was found.** It was not taken from a paper; it came from the criterion.
1. **What an excess requires.** By section 2.2 the quantum system must gain entropy *earlier* than the classical one.
2. **The generic quantum effect that does this** is the zero-point energy. It lowers the energy needed to escape a narrow trap, which pointed to a trap with about one bound state next to a large, classical reservoir.
3. **First test.** A Gaussian dimple in a harmonic trap with ħω = 1 gave only +0.10: the reservoir was too small.
4. **Box-like traps.** In flat-bottomed traps (x/L)^20 and (x/10)^6 the same dimple gave +0.6 to +2.9.
5. **The switch to sech².** The sech² (Pöschl–Teller) dimple came last, chosen for its exactly known bound states and its single knob λ.

The closest relative in the literature is the infinite well with an attractive δ at its centre in Gelbwaser-Klimovsky et al.; see section 4.4.

### 3.2 All candidates

This merges RESULTS A1–A3 and B4. Here T_q is the quantum escape temperature of the two-region model and G its reservoir partition function.

| potential | bound states | ZPE / D | ZPE / k_BT_q | ln G(T_q) | max(Cv_q − Cv_cl) | window (k_BT) |
|---|---|---|---|---|---|---|
| \|x\|ⁿ, n = 2…12 | – | – | – | – | < 0 everywhere | none |
| \|x\|ⁿ, n = 13 / 20 / 50 | – | – | – | – | +0.0015 / +0.018 / +0.043 | ≈ 2.6–2.8 (units of the \|x\|ⁿ scale) |
| box (hard walls) | – | – | – | – | +0.0625 (peak 0.5625) | every T > 1.28 E₁ |
| quartic double wells, 35 shapes | – | – | – | – | < 0 everywhere (largest −0.0001) | none |
| Gauss D=60 σ=0.1 in ½x² | 1 | 0.49 | 2.4 | 2.50 | +0.10 | 6.6–7.6 |
| Gauss D=60 σ=0.15 in (x/10)²⁰ | 2 | 0.36 | 1.9 | 3.39 | +0.03 | 6.6–7.0 |
| Gauss D=30 σ=0.1 in (x/5)²⁰ | 1 | 0.61 | 3.6 | 2.26 | +0.56 | 2.3–3.6 |
| Gauss D=80 σ=0.1 in (x/10)²⁰ | 2 | 0.44 | 2.7 | 3.46 | +0.85 | 6.6–9.8 |
| Gauss D=30 σ=0.1 in (x/10)²⁰ | 1 | 0.61 | 4.5 | 2.83 | +1.58 | 2.0–3.5 |
| Gauss D=30 σ=0.1 in (x/10)⁶ | 1 | 0.61 | 4.6 | 2.92 | +1.88 | 1.9–3.5 |
| Gauss D=30 σ=0.1 in (x/20)²⁰ | 1 | 0.61 | 5.4 | 3.42 | +2.86 | 1.7–3.2 |
| sech² λ=1.0 (D=100) in (x/10)⁶ | 1 | 0.50 | 3.7 | 3.73 | +2.09 | 7.0–11.6 |
| sech² λ=0.7 (D=59.5) in (x/10)⁶ | 1 | 0.59 | 4.8 | 3.33 | +2.42 | 3.7–6.6 |
| sech² λ=0.4 (D=28) in (x/10)⁶ | 1 | 0.71 | 6.8 | 2.72 | +2.28 | 1.4–2.9 |
| **sech² λ=0.5 (D=37.5) in (x/10)⁶** | 1 | 0.67 | 5.9 | 2.96 | **+2.43** | **2.0–4.0** |

What the table shows about the design rules of section 2.3 (c):
- **A large excess (> 1.5)** appears when R1 holds (one bound state, ZPE ≥ D/2), **ZPE/k_BT_q ≳ 3.7** (R2) and **ln G(T_q) ≳ 2.8** (R3).
- **Fails R2, marginal excess.** Below ZPE/k_BT_q ≈ 2.5 the excess stays at or below 0.1, however large the reservoir: the harmonic trap with ħω = 1, and the Gaussian with two bound states.
- **Fails R3, limited excess.** A small reservoir caps the excess, as for the (x/5)²⁰ box: +0.56 with ln G = 2.26.
- **A bigger reservoir helps.** Widening the box from L = 5 to 10 to 20 raises the excess from 0.56 to 1.58 to 2.86.
- **The model holds.** The two-region model reproduces the DVR quantum Cv inside every window to within 0.02–0.11 (RESULTS B4). The mechanism is the one described in section 2.3 (c).

---

## 4. Homogeneous eigen-energy scaling and the Otto-cycle relations

### 4.1 Why it is needed

This follows the meeting notes (Meetings_Summary.tex): start from Q_h = ⟨H_h⟩_A − ⟨H_h⟩_D for an Otto cycle between H_c and H_h.
- **The adiabatic strokes.** A quantum-adiabatic stroke keeps the populations (1-D bound spectra have no crossings). The state after it is a Gibbs state of the new Hamiltonian **iff every gap scales by one factor**: E_n^h − E_0^h = q·(E_n^c − E_0^c) for all n. It is then Gibbs at temperature qT.
- **Heat and work.** The heats are then integrals of Cv. With C_h(T) = C_c(T/q):

  Q_h = U_h(T_h) − U_h(qT_c) = ∫_{qT_c}^{T_h} C_h dT = q·∫_{T_c}^{T_h/q} C_c dT
  Q_c = −∫_{T_c}^{T_h/q} C_c dT,  W = Q_h + Q_c = (q − 1)·∫_{T_c}^{T_h/q} C_c(T) dT,  η = 1 − 1/q

- **The comparison.** At equal q, the efficiency is the same for the quantum and the classical working medium, and **W_q/W_cl = ∫Cv_q/∫Cv_cl over [T_c, T_h/q]**. The classical system obeys the same relations under the same stroke: a scale-invariant stroke maps classical Gibbs states to Gibbs states too (Deffner, Jarzynski, del Campo 2014, in `docs/refs`).
- **Without homogeneous scaling** the post-stroke state is not thermal, and Q ≠ ∫Cv dT.

### 4.2 Which strokes scale homogeneously

**Exact statement.** For **any** V, the scale-invariant family V_γ(x) = γ⁻²V(x/γ) has E_n(γ) = E_n(1)/γ². Proof: substitute y = x/γ; then −(ħ²/2m)d²/dx² = γ⁻²·[−(ħ²/2m)d²/dy²], so H_γ = γ⁻²·U H₁ U† with the unitary dilation U.

By contrast, the amplitude stroke V → κV and the pure dilation V(x) → V(x/γ) are homogeneous only for homogeneous potentials: the HO, |x|^k and the box.

**Checked numerically** in [RESULTS C](RESULTS.md) and [the figure](figures/homogeneity_gap_ratios.png). Each entry is the spread of the first 39 gap ratios q_n (homogeneous: q_n all equal):

| potential | amplitude V → 2V | dilation V(x/1.5) | scale-invariant V(x/1.5)/1.5² |
|---|---|---|---|
| HO ½x² | yes, q = √2 (spread 4e-13) | yes, q = 2/3 | yes, q = 0.4444 |
| \|x\|²⁰ | yes, q = 2^(1/11) = 1.0650 | yes, q = 1.5^(−20/11) = 0.4784 | yes, q = 0.4444 |
| double well 0.25x⁴ − 0.5x³ − 0.5x² | **no** (spread 0.43) | **no** (0.51) | yes (1e-12) |
| Gauss D=30 σ=0.1 in (x/10)²⁰ | **no** (0.52) | **no** (0.59) | yes (7e-13) |
| sech² λ=0.5 in (x/10)⁶ (the candidate) | **no** (0.34) | **no** (0.39) | yes (6e-13) |

Under the scale-invariant stroke the whole Cv curves map exactly, Cv[V₂](T/1.5²) = Cv[V₁](T): 5.7e-13 for the quantum curve and 2.4e-14 for the classical one. The simple relations Q, W ∝ ∫Cv dT therefore hold for the candidates, **provided the cycle's strokes are scale-invariant**. If a stroke changes only the trap, only the dimple, or only an amplitude, they do not hold.

**What the scale-invariant stroke means for each family:**
- **sech² dimple in (x/L)^k:** (D, σ, L) → (D/γ², γσ, γ^(1+2/k)L), so L → γ^(4/3)L for k = 6. λ (2mDσ²/ħ²) is invariant, so the number of bound states and the ZPE fraction stay fixed along the stroke, and the excess window slides to T/γ².
- **Gaussian dimple:** the same rule.
- **Quartic a x⁴ + b x³ + c x² + d x:** (a, b, c, d) → (a/γ⁶, b/γ⁵, c/γ⁴, d/γ³). This answers the open meeting question about constraints on (a, b, c): the two engine edges must share the invariants b⁶/a⁵, c³/a² and d²/a.
- **The ξ-scaling of the classical-limit method** (ξ²V at fixed ħ) is not such a stroke. It changes λ, i.e. the effective ħ. It is a mathematical device for reaching ħ → 0, not an engine stroke.

### 4.3 What the window is worth

For the recommended candidate, take a scale-invariant stroke with factor q, so both media have η = 1 − 1/q. From [RESULTS B1b](RESULTS.md):

| [T_c, T_h/q] | ∫Cv_q dT | ∫Cv_cl dT | W_q / W_cl |
|---|---|---|---|
| 2.03–4.03 (the whole window) | 6.53 | 3.43 | **1.90** |
| 2.76–3.64 (around the maximum) | 3.42 | 1.45 | **2.36** |
| 1–3 | 2.72 | 2.19 | 1.24 |
| 0.5–10 (wide) | 17.4 | 38.8 | 0.45 |

- **Inside the window** the quantum engine delivers about twice the classical work per cycle, at the same efficiency.
- **Over a wide range** the sum rule takes over and the classical engine wins. The cycle must be confined to the window, i.e. T_c and T_h/q must both lie inside it.
- **This is the plain Otto bookkeeping** with quasi-static strokes. Power, finite-time strokes and the cost of a scale-invariant (or shortcut-to-adiabaticity) drive are not included.

### 4.4 Relation to Gelbwaser-Klimovsky et al. (PRL 120, 170601, 2018)

**What the paper does.** Read from `docs/refs` (text extracted with `pdftotext`):
- **Main example.** An ideal-gas particle in a 1-D infinite well with a δ-function of strength g at its centre: a barrier for g > 0 or a well for g < 0 (both appear in their Fig. 2A).
  - The strokes switch the δ on and/or change the length L.
  - The δ shifts only the even levels, so the level scaling is **inhomogeneous**.
  - Classically the δ changes nothing. At fixed L (r = 1) the classical engine therefore extracts no work, while the quantum one does.
  - Their advantage is in **efficiency and operation mode**.
- **Experimental proposal (their Eq. 4).** An ion in a harmonic Paul trap plus an optical-lattice term that puts a barrier at the centre (flat-bottomed at κ = 1, a double well for κ > 1). The strokes tune κ and ω, again a change of shape, so again inhomogeneous.
- **Homogeneous scaling (their Eq. 2 and SI Eq. S3).** For E_n,h = q·E_n,c, W = ((1 − q)/q)·∫_{qT_c}^{T_h} Cv dT, and the same formula holds for a classical machine with the classical Cv. Quantum and classical machines then share the efficiency and the operation mode, which is the paper's sense of "behave like classical machines". **The work and the heats are still integrals of Cv**, so they differ exactly where Cv_q ≠ Cv_cl. For the HO, for instance, W_q < W_cl because Cv_q < 1.

**How our candidate relates to it:**
- **Same geometric family.** (x/10)⁶ − D·sech²(x/σ) is a smooth version of their infinite well with an *attractive* δ (g < 0): a wide, flat-bottomed trap with a narrow central well holding one bound state (λ ≤ 1). This is no coincidence. Both need a feature narrower than the thermal wavelength, which the quantum levels feel and classical mechanics feels differently or not at all.
- **Different stroke, different mechanism.**
  - **Theirs:** the stroke changes the narrow feature, which makes the scaling inhomogeneous and yields the efficiency/mode advantage. Q ≠ ∫Cv dT there.
  - **Ours:** the feature is kept and the whole potential is rescaled, (D, σ, L) → (D/γ², γσ, γ^{4/3}L). That is homogeneous (section 4.2, verified to 1e-12), so the efficiency is the same and the gain is in the work, through Cv_q > Cv_cl. This is exactly the homogeneous regime the paper sets aside, and its Eq. 2 is the formula used here.
- **A different classical side.** Their δ is invisible classically by construction. Our dimple has a finite depth and width: it traps a classical particle at low T and gives a large classical Cv peak (9.1 at T ≈ 5). Our excess comes from the zero-point energy moving the quantum escape to a lower T, not from classical mechanics ignoring the feature.
- **Homogeneity belongs to the stroke, not to the potential.** Even their δ-in-a-box has a homogeneous stroke: γ⁻²·g·δ(x/γ) = (g/γ)·δ(x), so L → γL together with g → g/γ scales every level by 1/γ². The dimensionless strength mgL/ħ² stays fixed, just as λ does for sech².

**So the candidate is relevant to our Cv-based check, provided the cycle's strokes are scale-invariant.** With strokes of the paper's kind (changing D, σ or L alone), the ∫Cv relations do not hold. The cycle must then be evaluated from the level populations directly (their Eqs. S1–S2), which is plan item V9.

**A practical caveat.** A scale-invariant stroke moves three knobs together. For an optical dimple σ is the beam waist, which is hard to change during a cycle (plan item V10).

---

## 5. DVR grid spacing: the problem and the `DVR_DX_MAX` fix

### 5.1 The problem

I added the sech² candidate to [`src/config.py`](../../src/config.py) as an inactive block. Switched on in a scratch copy of `src/`, the quick scan (resolution 1, k_BT = 0.5–10) failed at its first solve with `DVR CONVERGENCE ERROR` (resolution error 1.53).

**Diagnosis.** It was made with the pipeline's own `auto_configure_dvr`, `_wkb_*` and `colbert_miller_dvr_1d` (the sech² ground state is exactly E₀ = −12.5):

| solve | grid (auto) | dx / σ | E₀ on that grid | finer grids |
|---|---|---|---|---|
| T = 10 (223 levels) | ±42, 1077 pts | 0.78 | −12.514 | −12.5000 at dx/σ ≤ 0.26 |
| T = 3 (41 levels) | ±232, 2962 pts | 1.57 | −10.744 | −12.4999 at dx/σ = 0.52, −12.5000 at 0.39 |

**Root cause.** `auto_configure_dvr` sets dx = π/(2k_max) from the *wavefunction* momentum at the energy ceiling. That rule does not see that V itself has a feature of width σ. Two consequences:
- The bound state, the very state the effect is about, is under-sampled.
- Moving the grid points changes it, so the span-expansion loop also misfires: the T = 3 grid grew to ±232.

**No parameter choice avoids this.** dx/σ is invariant under the scale-invariant transformation. And the effect *requires* σ to be below the thermal wavelength (R4), which is exactly where the momentum rule gives dx ≳ σ.

### 5.2 The fix (implemented 2026-10-13, with your approval)

**The idea.** The DVR has two separate resolution needs:
- **the wavefunctions.** Their shortest wavelength is set by the highest energy solved for, and the momentum rule dx = π/2k_max handles it.
- **The potential itself.** The DVR samples V only at the grid points, so the grid must also resolve V's narrowest feature.

For every potential the project has used so far (HO, quartic wells), the first need is always the stricter one, so only it was coded. A narrow dimple reverses this. The fix lets a potential declare its own feature scale: **dx = min(π/2k_max, DVR_DX_MAX)**.
- **`DVR_DX_MAX` belongs to the potential**, like its parameters. It defaults to `None`, meaning no cap, and the sech² block sets it to σ/4.
- **Every pipeline solve applies it**: base, ξ²V, references, quick scan and coefficient sweep. ξ²V has the same feature widths as V, so one value serves every ξ. At large ξ the momentum rule is finer anyway and the cap is inactive.
- **It is threaded explicitly**, as a keyword argument, through every function that builds a grid, with default `None`. This is the same pattern as the earlier `energy_ceiling`/`padding` additions: no global state, and with `None` every code path is unchanged.

**The exact code changes** (`git diff` against 313e442):

| file | change |
|---|---|
| `src/DVR/DVR_Algorithm.py` | `auto_configure_dvr(..., dx_max=None)`: documented. After `dx_target = np.pi / (2.0 * k_max)` it adds `if dx_max is not None: dx_target = min(dx_target, float(dx_max))`. Nothing else in the DVR changes; the 3-pass check and the reference grids scale from this dx. |
| `src/Classical_Limit_Numerical.py` | `ScaledSpectra(..., dx_max=None)` stores `self.dx_max` and passes it to `auto_configure_dvr` in `_solve`. `sweep_temperature_range(..., dx_max=None)` passes it to the `ScaledSpectra` it builds when none is given. |
| `src/Quick_Scan.py` | `_new_spectra` passes `dx_max=config.DVR_DX_MAX`. |
| `src/Quantum_Classical_Combined.py` | `run(..., dx_max=None)` passes it to `sweep_temperature_range` (Section 4). |
| `src/Cv_Numerical_Benchmark.py` | `_run_cv_pipeline(..., dx_max=None)` and `run_cv_numerical_benchmark(..., dx_max=None)` pass it through (Section 6). |
| `src/Cv_Coefficient_Sweep.py` | `solve_variant_spectrum`, `solve_variant_with_hot_coverage`, `_variant_classical_limit` and `run_coefficient_sweep` each take `dx_max=None` and pass it down to `auto_configure_dvr` / `sweep_temperature_range` (Section 7). |
| `src/Quantum_HO_Master.py` | Imports `DVR_DX_MAX` from config and passes it to `auto_configure_dvr` (Section 1), `run_general_cv_pipeline` (Section 4), `run_cv_numerical_benchmark` (Section 6) and `run_coefficient_sweep` (Section 7). |
| `src/figures/plot_potential.py` | The standalone `main()` passes `dx_max=DVR_DX_MAX`. |
| `src/config.py` | `DVR_DX_MAX = None` after the physical constants, with its rationale. The sech² block sets `DVR_DX_MAX = POTENTIAL_PARAMS["s"] / 4`. |

**Not changed:**
- `colbert_miller_dvr_1d` and `get_fully_converged_energy_levels`;
- `DVR_Reference_Generator` (it refines the already-capped base grid);
- `DVR_Limit_Finder` (it scans dx itself).

**Regression** ([`dx_max_regression_check.py`](dx_max_regression_check.py), output [`DX_MAX_REGRESSION.md`](DX_MAX_REGRESSION.md)). Before any code change, a baseline was recorded on the committed double well 0.25x⁴ − 0.5x³ − 0.5x². It covers four paths:
- the base spectrum (NUM_STATES = 500 levels, 3-pass check);
- `sweep_temperature_range` building its own cache;
- the same with span/dx factors 1.5;
- `Quick_Scan.scan_potential` at resolution 1.

After the change, every quantity (energies, grid, classical Cv, error estimates, ξ_conv, quick-scan curves) is **bit-identical**, both with `DVR_DX_MAX = None` and with a cap looser than the momentum rule (10).

### 5.3 The candidate through the pipeline after the fix (plan item V3, partial)

These were run on a scratch copy of the new `src/` with the sech² block switched on:
- **The quantum side works.** `ScaledSpectra` (the quick scan's cache, coverage 12, tolerance 1e-3) gives E₀ = −12.500000 at ξ = 1. The quantum Cv matches the exact one within 1.3e-3 over k_BT = 1.5–5. The same call without the cap still fails with `DVR CONVERGENCE ERROR` (resolution error 0.28).
- **The classical ξ-scan does not finish in the window on this machine (15.8 GB RAM).**
  - Over k_BT = 0.5–10 it spent more than 15 minutes on the hottest temperature, growing past 2.3 GB.
  - Over k_BT = 1.5–5 it passed 3.6 GB at T = 5, still before the plateau, and I stopped both runs.
  - The reason is built into the candidate. The dimple is deeply quantum in the window (ħω_d/k_BT ≈ 17–58), so the ħ² law needs ξ ≈ 30–80 for the plateau. A *uniform* grid that resolves the dimple at that ξ across the whole trap needs 10⁴–10⁵ points, beyond a dense diagonalisation.
- **Thermal coverage leaves a small bias.** The scan keeps levels up to ≈ V_min + 1.2·c·k_BT (c = the preset's coverage). For the dimple this cuts into the trap's states at low T. Computed exactly with the cut-off phase-space integral, the classical Cv comes out low by up to **0.030 at c = 12 (resolution 1)**, ≤ 0.002 at c = 15 and ≤ 0.001 at c = 20. That is negligible next to the excess of 2.4, but at resolution 1 it is 3× tol_xi, and the n-scan cannot see it.

**Options for verifying the classical side** (not implemented; your choice):
1. **The exact integral ½ + β²Var(V)** (the open "option A"), as a check on this potential only. It costs a second.
2. **Richardson extrapolation in ξ** (open item 2). Fitting Cv(ξ) = Cv_cl + a/ξ² + b/ξ⁴ on moderate ξ (≲ 15) would cut the ξ needed several-fold.
3. **A mapped, non-uniform grid** (the adaptive-mesh proposal in Meetings_Summary.tex). Fine where the dimple is, coarse in the flat trap, it attacks the two-length-scale cost directly.

---

## 6. Verification plan

Run in order; each step has a pass criterion.

| # | Check | How | Pass if |
|---|---|---|---|
| V1 | Reproduce this folder | `python verification/quantum_excess/quantum_excess_check.py` (~5 min) | RESULTS.md unchanged |
| V2 | DVR grid fix | **done 2026-10-13:** `DVR_DX_MAX` added (section 5.2); `dx_max_regression_check.py` against a pre-change baseline | bit-identical results: **passed** on all four paths, with the cap off and with a loose cap |
| V3 | Quick scan of the sech² candidate | activate the config block, resolution 1, then 2 | **partial (section 5.3).** Quantum side **passed** (E₀ = −12.500000, Cv_q within 1.3e-3). Classical ξ-scan **did not finish** in the window (memory). Still open: classical within 2ε of ½ + β²Var(V), and verdict "quantum above" with \|d\| > 2ε |
| V4 | Full pipeline on the candidate | `Quantum_HO_Master.py` (Sections 1–6), with `SCAN_PARAM = "D"` | Section 2 and Section 6 references agree with the base run; same window as V3 |
| V5 | Sum rule on the pipeline's own curves | ∫(Cv_q − Cv_cl)dT over the sweep against ΔU(T_max) − ΔU(T_min) (U_q from the spectrum, U_cl = T/2 + ⟨V⟩ as a check only) | agreement within the ξ-scan's integrated error estimate |
| V6 | Map the design rules | scan λ ∈ [0.3, 1.3], L/σ ∈ [25, 400], trap exponent k ∈ {2, 4, 6, 20} | excess appears exactly where R1–R3 hold; boundary located |
| V7 | Shape independence | sech² vs Gaussian vs Lorentzian dimple at equal λ-like ZPE fraction | same window to within the ZPE/E_b differences |
| V8 | Homogeneity of the actual cycle | for the chosen (V_c, V_h) pair: gap ratios (as in RESULTS C) and Cv_h(T) = Cv_c(T/q) | spread < 1e-6; curves equal to DVR precision |
| V9 | Otto cycle, directly | propagate populations through the four strokes (adiabatic = same n, isochores = Gibbs) and compute Q_h, W | equals (q − 1)∫Cv dT; W_q/W_cl as in section 4.3 |
| V10 | Physical units | map (D, σ, L) to an atom in an optical trap with a dimple (recoil units) | parameters feasible |

---

## 7. Files

- `SUMMARY.md`: this file.
- [`quantum_excess_check.py`](quantum_excess_check.py): every number above. It imports the pipeline's DVR from `src/`.
- [`RESULTS.md`](RESULTS.md): its output.
- [`figures/`](figures/): `best_candidate_cv.png`, `delta_u_criterion.png`, `homogeneity_gap_ratios.png`.
- [`src/config.py`](../../src/config.py): the inactive sech² block, with the swap instructions and the known limitation of section 5.

**Related:**
- [FINDINGS.md](../../FINDINGS.md): the Wigner–Kirkwood merge and the box overshoot.
- [CHANGES_SUMMARY.md](../../CHANGES_SUMMARY.md) F: the exact classical Cv and the ħ → 0 limit.
- [audit/quick_scan/MERGE_CHECK.txt](../../audit/quick_scan/MERGE_CHECK.txt).
