# Proof sketch: why Cv_q > Cv_cl needs a low-entropy trap inside an entropy-rich reservoir

**Status: sketch only (2026-10-15).** This document sets out a proof strategy, the lemmas it needs, and which steps are already rigorous, which are plausible, and which are open. The worked example is the reflectionless Pöschl–Teller well in a box ([PT_Box_Cv.tex](PT_Box_Cv.tex)), because every object below has a closed form there.

The numbers marked *(scratch)* come from exploratory checks run on 2026-10-15. They are not yet in [`pt_box_check.py`](pt_box_check.py), and should be added there before they are cited.

Units are ħ = m = k_B = 1 unless written out. The exact anchors of Section 1 are already derived in [`../quantum_excess/SUMMARY.md`](../quantum_excess/SUMMARY.md), section 2.

---

## 0. The claim, informally

Take a 1-D system H = p²/2m + V(x) made of two parts:

- a **trap**: a narrow region whose quantum spectrum below some threshold E_th is one level (or a few, widely spaced);
- a **reservoir**: a wide, flat region above E_th with a dense spectrum and large entropy.

Then:

1. **Sufficiency.** If the reservoir's entropy is large enough, Cv_q > Cv_cl in a temperature window. The window lies on the low-T flank of the quantum "desorption" peak.
2. **Asymptotic version (the cleanest target).** For any trap, max_T (Cv_q − Cv_cl) → +∞ as the reservoir size L → ∞, roughly like ¼(ln L)².
3. **Necessity (conjecture).** Without the two-phase structure there is no excess, apart from a weak hard-wall effect of order λ_th/L. "Two-phase" means a trap phase whose quantum free energy is shifted up strongly relative to its classical one, plus a reservoir phase where the shift is negligible.

---

## 1. Exact anchors (known)

- **(A1)** ΔCv = Cv_q − Cv_cl = −T·∂²ΔF/∂T², with ΔF = F_q − F_cl. So **an excess at T ⇔ ΔF is concave at T**.
- **(A2) Golden–Thompson:** Z_q ≤ Z_cl, i.e. ΔF ≥ 0 at every T.
  - This also holds with Dirichlet walls, because the Dirichlet heat kernel is bounded by the free one.
  - Checked on PT-in-box *(scratch)*: min over θ of Z_cl/Z_q is 1.026, 1.008, 1.002 and 1.0008 for αL = 30, 100, 400, 1000.
- **(A3) Limits.** ΔF(0) = E₀ − V_min ≡ ZPE > 0. For smooth confining V, ΔF(∞) = 0 (Wigner–Kirkwood). Hence the **sum rule** ∫₀^∞ ΔCv dT = −ZPE.
- **(A4) Hard walls are the exception.** A Dirichlet wall adds −½ to Z_q, so ΔF ⊃ +k_BT·λ_th/(2L) ∝ √T. That term is concave, so walls alone give ΔCv → λ_th/(8L) > 0 at high T, and the sum rule fails (PT_Box_Cv §6.3). Any necessity statement must be made modulo this term, or for smooth V.

---

## 2. Definitions: two phases, split at the threshold

Choose the threshold E_th: the reservoir floor (0 in the example). Split both systems **by energy**, the same way on both sides:

| | Quantum | Classical |
|---|---|---|
| Trap phase **B** | eigenstates with E_n < E_th | phase-space region H(x,p) < E_th |
| Reservoir phase **R** | eigenstates with E_n > E_th | phase-space region H > E_th |

For each phase X ∈ {B, R} define Z_X, the weight p_X = Z_X/Z, the mean energy U_X, and the phase heat capacity c_X = β²Var_X(E).

Three derived quantities control everything:

- **Gap:** Δ(T) = U_R − U_B.
- **Trap ratio:** r(T) = Z_B^cl / Z_B^q.
- **Reservoir ratio:** ρ(T) = Z_R^q / Z_R^cl. It tends to 1 for a dense, quasi-classical reservoir.

The quantum split is exact, because eigenstates do not mix in equilibrium. There is no tunnelling issue.

**Example (closed forms).** With λ_th = √(2πβ):

| Quantity | Closed form |
|---|---|
| Z_B^q | e^{βε} (one level, so c_B^q = 0 and S_B^q = 0 exactly) |
| Z_R^q | L/λ_th − ½ − e^{βε}·erfc√(βε) |
| Z_B^cl | (1/λ_th) ∫ e^{−βV} erf√(−βV) dx (momenta restricted to H < 0) |
| Z_R^cl | Z_cl − Z_B^cl, with Z_cl given by the ₂F₂ form of PT_Box_Cv §5 |

---

## 3. Lemma 1: two-phase decomposition (exact)

By the law of total variance, for **any** partition of states into B and R:

  **Cv/k_B = p_B·c_B + p_R·c_R + β²·p_B·p_R·Δ²**

All three terms are non-negative. The last one is the "transition" (Schottky) term.

At the phase crossing (p_B = p_R = ½, so F_B = F_R and Δ = T·ΔS), the transition term equals **(ΔS)²/4**. The height of a desorption peak is therefore set by the entropy jump between the phases.

**Check on the example** *(scratch)*: the decomposition reproduces Cv_q to every printed digit. At the quantum crossing:

| αL | S_R | S_R²/4 | + p_R·c_R | = Cv_q |
|---|---|---|---|---|
| 1000 | 5.35 | 7.15 | 0.25 | 7.40 |
| 100 | 3.30 | 2.72 | 0.26 | 2.97 |

**Status:** rigorous.

---

## 4. Lemma 2: the quantum trap is less stable than the classical one (r ≥ 1, and r grows as T → 0)

**Claim.** r(T) ≥ 1, and for a trap with a harmonic bottom of frequency ω_d,

  r(T) ≈ (k_BT/ħω_d)·e^{β·ZPE},

which diverges as T → 0.

**Route.**
1. Golden–Thompson applied to a *trap Hamiltonian*.
2. The energy-split phase B is not itself a Hamiltonian, so it needs to be made rigorous with **Dirichlet–Neumann bracketing**: insert walls at the trap edges. That bounds Z_B^q above and below by trap-in-a-box spectra, to which Golden–Thompson applies.
3. The low-T asymptotics follow from Laplace's method (classical side) and the ground-state term (quantum side).

**Check on the example** *(scratch)*: r is independent of L, as it should be.

| θ | r measured | harmonic estimate |
|---|---|---|
| 0.11 | 355.5 | 345 |
| 0.25 | 5.2 | 4.8 |

The classical trap heat capacity is c_B^cl = 1.07–1.27 over this range: it approaches 1 by equipartition, plus anharmonic softening. The quantum value is c_B^q = 0.

**Status:** the inequality is plausible-to-rigorous via bracketing. The asymptotics are standard.

---

## 5. Lemma 3: window criterion (activation regime)

**The excess is not at the crossing.** At the quantum crossing (p_R^q = ½) the classical activation term has already caught up. *(scratch)*: ΔCv(θ*_q) is negative for every αL tested:

| αL | 30 | 100 | 400 | 1000 |
|---|---|---|---|---|
| ΔCv(θ*_q) | −1.14 | −3.17 | −3.19 | −0.68 |

The window is on the **rising flank**, where p_R^q is noticeable but p_R^cl is still tiny.

At αL = 1000 *(scratch)*:

| θ | p_R^q | p_R^cl |
|---|---|---|
| 0.15 | 0.12 | 0.003 |
| 0.19 | 0.39 | 0.044 |

**Flank form.** Since the odds satisfy p_R/p_B = Z_R/Z_B, we have p_R^cl ≈ p_R^q/(r·ρ). For p_R ≪ 1, Lemma 1 gives

  **ΔCv ≈ p_R^q·[(βΔ_q)² − (βΔ_cl)²/(r·ρ)] − c_B^cl + O(p_R·c_R)**,

using c_B^q = 0.

**Criterion.** A window exists if, at some T,

1. **r·ρ > (Δ_cl/Δ_q)²**: the quantum trap activates first.
   - In the example Δ_cl/Δ_q ≈ (2 − θ/2)/(1 + θ/2), so this needs r ≳ 3. That happens for θ ≲ 0.3, which is why every window ends near θ ≈ 0.20–0.24.
2. **p_R^q·(βΔ_q)² exceeds c_B^cl ≈ 1** while condition 1 still holds: the quantum activation is strong enough to beat the classical trap's own equipartition heat capacity.

The **three ingredients of the claim** map onto these terms:

| Ingredient | Where it enters |
|---|---|
| Low-entropy trap | c_B^q = 0 and S_B^q = 0: no internal heat capacity, and the least stable possible bound phase |
| Zero-point shift | r ≫ 1, a property of the trap alone |
| Entropy-rich reservoir | a large Z_R pushes activation to low T, where βΔ_q is large, so condition 2 is reached while r is still big |

**Status:** heuristic as written. To make it rigorous:
- keep the (1 − p) factors exactly instead of expanding;
- bound ρ − 1 and its β-derivatives (Beth–Uhlenbeck or Levinson for the reservoir: corrections of order λ_th/L);
- bound c_B^cl from above.

---

## 6. Theorem target: the large-reservoir limit

**Statement (target).** Let the trap have a non-degenerate ground state E₀ > V_min (always true) and a threshold gap Δ_q = E_th − E₀ > 0. Embed it in a flat reservoir of length L: a box, or a flat floor much deeper than the relevant T. Then

  **sup_T [Cv_q(T) − Cv_cl(T)] → +∞ as L → ∞**, growing like ¼·[ln(L/λ_th)]².

**Proof route.**

1. **Quantum activation temperature.** e^{βΔ_q} ≈ Z_R ≈ L/λ_th gives T*_q ≈ Δ_q/ln(L/λ_th) → 0.
2. **Low entropy is automatic.** At T*_q → 0, the trap's excited levels are frozen out (T*_q ≪ E₁ − E₀). So S_B^q → 0 and c_B^q → 0, even if the trap has several bound states.
3. **Classical activation comes later by a fixed ratio.** T*_cl ≈ Δ_cl/ln(L/λ_th) up to log-log corrections, with Δ_cl → E_th − V_min. The ratio T*_cl/T*_q → (E_th − V_min)/(E_th − E₀) = 1 + ZPE/Δ_q > 1, while both transitions sharpen: their relative width ~1/ln L. So at T*_q the classical reservoir weight p_R^cl → 0 faster than any power of 1/ln L.
4. **Classical side stays bounded.** Cv_cl(T*_q) → c_B^cl(T*_q) → 1 (harmonic bottom, equipartition).
5. **Quantum side grows.** Cv_q near T*_q is about S_R²/4 → ∞ by Lemma 1, with S_R ≈ ln(L/λ_th) + ½.
6. Therefore ΔCv → +∞. To finish: a uniform Laplace-method estimate for steps 3–4, and the correction from evaluating slightly below T*_q (Lemma 3).

**Corollaries.**

- **(a) Every trap works for large enough L,** including traps with several bound states.
- **(b) The required size explodes in the classical limit.** The relative separation ZPE/Δ_q must beat the relative width ~c/ln L, so ln L_min ~ Δ_q/ZPE. As ZPE/depth → 0 (for example ξ → ∞), L_min → ∞. This is consistent with the ξ formula recovering Cv_cl at fixed L.
- **(c) Peak temperature.** θ at the maximum drifts like ~1/ln Λ.

**Numerical support** *(scratch, closed forms, PT trap of parameter λ in a box)*. Entries are max_T ΔCv:

| λ (ZPE/depth) | αL = 10² | 10³ | 10⁴ | 10⁵ | 10⁶ | 10⁸ |
|---|---|---|---|---|---|---|
| 2 (0.50) | +0.34 | +4.80 | +12.65 | +22.82 | +34.98 | +66.34 |
| 3 (0.33) | 0 | +2.42 | +8.58 | +19.00 | +33.25 | +69.63 |
| 4 (0.25) | 0 | +1.19 | +5.24 | +12.96 | +25.06 | +61.97 |
| 6 (0.17) | 0 | +0.34 | +2.24 | +6.28 | +13.31 | +39.13 |

- For λ = 2 and αL ≥ 10⁴ the data fit ¼(ln Λ − b)² with b ≈ 2.0–2.1.
- θ at the maximum times ln Λ stays at about 1.15–1.21.
- Smaller ZPE/depth needs a larger reservoir, as corollary (b) predicts.

**Status:** each step is standard asymptotics. The main work is a uniform error control in step 3.

---

## 7. Finite-size criterion: what "entropy-rich and shallow" must mean quantitatively

From Lemma 3, the reservoir must supply roughly **S_R ≳ 2** at a temperature where r(T) ≫ (Δ_cl/Δ_q)². In the example this is marginal at αL = 30 (no excess) and enough at αL = 100.

How fast different reservoirs gain entropy:

| Reservoir | S_R(T) | Verdict (scratch, envelope ε(2x/L)^{2k} + PT, αL = 100) |
|---|---|---|
| box / flat floor of length L | ln(L/λ_th) + ½ | best: +0.34 |
| power-law envelope x^{2k} | ln(L/λ_th) + (1/2k)·ln(k_BT/ε) + const | +0.09 (x⁴), +0.17 (x⁶), +0.25 (x¹⁰), +0.32 (x²⁰) |
| harmonic envelope (k = 1) | 1 + ln(k_BT/ħω_env) | no excess unless ħω_env ≪ T* |

"Shallow" means two things:
- the reservoir floor lies **above** the trap's quantum ground state (Δ_q > 0);
- the reservoir is **dense**: level spacing ≪ k_BT*, so that ρ ≈ 1 and the reservoir carries no quantum deficit of its own.

Flatness matters because S_R grows fastest when the accessible length does not shrink with T.

---

## 8. Necessity (conjectures and routes)

- **(N1) One phase only ⇒ no excess (smooth V).**
  - *Conjecture:* if V has a single basin and the quantum correction to F varies slowly across the thermally visited region, ΔF is convex at every T.
  - *Evidence:* the HO (exact); |x|ⁿ for n ≤ 12; the whole quartic family; single PT and Morse wells in bound-only models; the harmonic envelope above.
  - *Routes:*
    - the sandwich F_cl ≤ F_q ≤ F_FK: Golden–Thompson below, the Feynman–Kleinert variational bound (Gibbs–Bogoliubov) above;
    - Wigner–Kirkwood at high T (ΔF ≈ ħ²β⟨V″⟩/24m, convex for power laws);
    - exponential approach to ZPE at low T.
  - The hard part is controlling the convexity of ΔF itself, not just bounds on ΔF.
- **(N2) The trap must be low-entropy at T*.** If T* is much larger than the trap's level spacing, Wigner–Kirkwood applies inside the trap, r → 1 + O(ħ²), and the two activations merge. This is the ξ → ∞ statement.
- **(N3) The reservoir must carry entropy.** If Z_R is bounded, as in bound-only models without a box, the quantum activation term cannot exceed c_B^cl before the classical one catches up. Single plateau wells in the bound-only model give no excess (numerics, 2026-10-15).
- **(N4) Exception: hard or very steep walls.** These give a trap-free excess ≤ λ_th/(8L), and up to 0.04–0.06 for |x|ⁿ with n ≥ 13 (SUMMARY.md). A clean necessity theorem must either exclude it or state the excess threshold above it.

---

## 9. Toolbox: where each rigorous tool enters

| Tool | Used for |
|---|---|
| Law of total variance | Lemma 1 (exact) |
| Golden–Thompson, including Dirichlet walls | ΔF ≥ 0; r ≥ 1 |
| Dirichlet–Neumann bracketing | making the trap/reservoir split into Hamiltonians (Lemma 2) |
| Wigner–Kirkwood expansion | high-T behaviour; N2; ρ → 1 |
| Feynman–Kleinert / Gibbs–Bogoliubov | upper bound on F_q (N1) |
| Beth–Uhlenbeck / Levinson | reservoir correction ρ − 1 = O(λ_th/L) (Lemma 3) |
| Laplace's method | T → 0 and L → ∞ asymptotics (Lemma 2, Theorem) |

---

## 10. Plan: validate each step on the PT-in-box example first

Every quantity is explicit (Section 2), so each lemma can be checked before it is proved generally:

1. Lemma 1 at all θ (exact identity; a regression test).
2. Lemma 2: r(θ) against the harmonic estimate and against the Golden–Thompson bound; c_B^cl(θ).
3. Lemma 3: predict the window edges from the flank inequality and compare with the exact window (PT_Box_Cv, Table 5).
4. Theorem: max ΔCv against ln Λ up to 10⁸ for λ = 2, 3, 4, 6. Fit the ¼(ln Λ − b)² law and θ_max·ln Λ.
5. Corollary (b): along the ξ ladder (λ = n), find Λ_min(n) at which an excess first appears, and test ln Λ_min ∝ Δ_q/ZPE.

Items 1–5 would extend `pt_box_check.py` (new sections plus a RESULTS.md table) and become the numerical appendix of the eventual proof.

---

## 11. Open questions

- Uniform error control in the Theorem (step 3): how large must L be, explicitly, as a function of ZPE/Δ_q?
- Is ¼(ln L/λ_th)² the exact leading term? What is the subleading constant (b ≈ 2 in the example)?
- N1 for smooth single-basin potentials: is there a convexity proof of ΔF, or a counterexample short of steep walls?
- Several traps, or a trap with nearly degenerate levels: does a low-entropy *manifold* (S_B = ln g) still work? Presumably yes, with S_R − S_B replacing S_R.
- Effect on the engine: the window is on the flank and is followed by a larger deficit. Which Otto-cycle temperatures maximize W_q/W_cl under scale-invariant strokes?
