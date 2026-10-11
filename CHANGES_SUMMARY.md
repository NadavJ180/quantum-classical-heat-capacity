# Changes Since the Scaling Error: What, Where, and Why

A short, physics-first list of every change made on branch `correct-classical` after we found that the classical-limit ξ-scan never scaled the potential. For each change it gives:
- what changed;
- where in the code (line numbers as of commit `9e8621b`);
- the physical principle it relies on.

"Numerics only" marks changes with no physics behind them. The detailed, file-by-file change log is [`audit/classical_limit/CHANGES.md`](audit/classical_limit/CHANGES.md), and the narrative is in [`HISTORY.md`](HISTORY.md).

---

## A. Fixing the classical limit (commit `f747e61`)

**A1. Solve ξ²V at every ξ, not V.**
- **Where:** `ScaledSpectra._solve`, [`src/Classical_Limit_Numerical.py`](src/Classical_Limit_Numerical.py#L275) line 275.
- **Principle:** the classical limit is ħ → 0 with V and T fixed. Factoring ξ² out of the Hamiltonian,

  −(ħ²/2m) d²/dx² + ξ²V = ξ² · [−((ħ/ξ)²/2m) d²/dx² + V],

  so the spectrum of ξ²V is ξ² times the spectrum of V with a smaller ħ (Gelbwaser Eq. S7). The old code reused V's spectrum, which only rescaled T.
- **Not changed:** the division by ξ² ([lines 138, 147](src/Classical_Limit_Numerical.py#L138)) was already in the code. Boltzmann factors depend only on E/k_BT, so dividing E by ξ² is the same as raising T to ξ²T. It now finally receives the right spectrum.

**A2. Keep every level up to E₀ + 20·k_BT, counted in advance.**
- **Where:** [lines 272](src/Classical_Limit_Numerical.py#L272) and [304](src/Classical_Limit_Numerical.py#L304), using `_wkb_level_count` at [line 181](src/Classical_Limit_Numerical.py#L181).
- **Principle 1 (Boltzmann factor):** a level 20·k_BT above the ground state has weight e⁻²⁰ ≈ 2·10⁻⁹, which is negligible.
- **Principle 2 (WKB, Bohr–Sommerfeld):** the number of levels below E ≈ (area of the classical orbit in phase space)/(2πħ) = (1/πħ)∮p dx. That gives how many levels to ask the DVR for.

**A3. One DVR solve per ξ, with temperatures processed hot → cold.**
- **Where:** [line 583](src/Classical_Limit_Numerical.py#L583).
- **Principle:** a spectrum complete enough for a hot temperature also contains everything a colder one needs, by the Boltzmann factor again. This is purely about speed.

**A4. An error estimate on the classical value: ε = |ΔCv|/(m² − 1).**
- **Where:** [line 392](src/Classical_Limit_Numerical.py#L392).
- **Principle:** quantum corrections to thermodynamics start at order ħ² (Wigner–Kirkwood). Since ħ → ħ/ξ, Cv(ξ) ≈ Cv_cl + a/ξ². On a ladder ξ_k = m·ξ_{k−1}, one step changes Cv by (a/ξ_k²)·(m² − 1), so the distance still left to the limit is |ΔCv|/(m² − 1). This is the same idea as Richardson extrapolation.

**A5. A plateau only counts if Cv ≥ ½ − tol.**
- **Where:** [line 394](src/Classical_Limit_Numerical.py#L394).
- **Principle:** Cv_cl = ½ + β²·Var(V).
  - The ½ is the kinetic energy's share: the momentum integral is Gaussian (equipartition).
  - A variance is never negative, so the classical Cv can never be below ½.
- **What it fixes:** it rejects the frozen-out "plateau" at Cv ≈ 0 at low T.
- **OPEN: still to be verified and proven.** This guard is the one place where the pipeline relies on classical theory (Cv_cl ≥ ½, section F) rather than on the ξ-scan alone. It is marked with a comment in the code. Either prove that keeping it is safe here, or replace it with a criterion that uses only the ξ-scan. The options proposed (not implemented):
  1. **A relative stability test, |ΔCv|/(m² − 1) < tol·|Cv|.** In the frozen regime Cv grows exponentially from rung to rung, so its relative change is huge and the step is rejected. `TOL_XI` would become a relative tolerance.
  2. **A check that successive changes shrink by about 1/m² per step:** the ħ² law, which the error estimate A4 already assumes. In the frozen regime the changes grow instead, so the step is rejected. Needs a small floor against round-off. *(Recommended, possibly combined with 1.)*
  3. **A minimum number of thermally populated levels,** N_eff = 1/Σp_n² ≫ 1 (correspondence principle). Adds a threshold knob, and needs care with tunnelling doublets.

  Whichever is chosen, the engine verification against the exact Cv (HO, x⁴, double well) and the quick-scan validation should be re-run to show that nothing else changes.

**A6. Report the last plateau point, and NaN when no plateau is reached.** *Numerics only.*
- **Where:** [lines 403](src/Classical_Limit_Numerical.py#L403) and [589](src/Classical_Limit_Numerical.py#L589).
- **Why:** the last point is the most classical value computed. The quantum Cv used to be stored when the scan failed, and it is not a classical value.

**A7. DVR grid options for ξ²V: energy ceiling, turning points, padding.**
- **Where:** [`src/DVR/DVR_Algorithm.py`](src/DVR/DVR_Algorithm.py#L145), lines 36 and 145–172.
- **Principle (sampling):** the grid must resolve the shortest de Broglie wavelength, so Δx ≤ π/(2k_max) with k_max = √(2m(E_max − V_min))/ħ.
  - ξ²V has ξ times more levels per unit energy, so E_max now comes from WKB rather than the harmonic-oscillator guess.
  - The window comes from the classical turning points, where V(x) = E.

**A8. New ladder settings:** `XI_START` = 1, `XI_MULT` = 1.25, `TOL_XI` = 2e-3, `MAX_XI_STEPS` = 35, `XI_MAX` = 2000.
- **Where:** [`src/config.py`](src/config.py#L148), lines 148–153.
- **Why:** ξ = 1 is the physical system itself. The cap exists because the grid grows ∝ ξ and a dense diagonalization costs ∝ (grid)³.

**A9. Callers adapted.**
- **Quantum Cv:** taken from the base spectrum only ([`src/Quantum_Classical_Combined.py`](src/Quantum_Classical_Combined.py#L336), line 336).
- **Auto-tune escalation** ([`src/Cv_AutoTune.py`](src/Cv_AutoTune.py#L162), lines 162–163): more levels if the quantum curve is truncated, and a longer ladder if the ξ-scan runs out.
- **Section 6:** re-solves ξ²V on refined grids.
- **HO docstring:** its classical Cv = 1 because V is quadratic (two quadratic terms × ½ each, by equipartition), not by definition.

## B. Section 7 compares each variant with its own classical limit (commit `f3146ae`)

- **Where:** [`src/Cv_Coefficient_Sweep.py`](src/Cv_Coefficient_Sweep.py#L979), lines 920 and 979; plot at [line 442](src/Cv_Coefficient_Sweep.py#L442).
- **Principle:** the classical Cv depends on the potential (½ + β²·Var(V)). "Quantum above classical" only means something for the same Hamiltonian.

## C. Section 6 uses a lighter grid for its classical reference (commit `59dc07d`)

- **Where:** [`src/config.py`](src/config.py#L182), lines 182–184, and [`src/Quantum_HO_Master.py`](src/Quantum_HO_Master.py#L377), lines 377–400.
- **Principle (grid-convergence test):** a converged DVR result doesn't change when you change the grid. Sinc-DVR errors shrink exponentially as the grid spacing and span are refined, so a 50% change (1.5/1.5) exposes an unconverged grid just as well as 2/2, at a quarter of the cost. The 2/2 grid remains available via `CLASSICAL_REFERENCE_PRECISE`.

## D. The quick scan (commits `0cad903`, `4707a51`, `9e8621b`)

**D1. Same ξ method at coarser settings.**
- **Where:** [`src/Quick_Scan.py`](src/Quick_Scan.py), presets at [`src/config.py`](src/config.py#L219) line 219.
- **Principle:** by the ħ² law, the ξ needed grows ∝ 1/√tol. The grid grows ∝ ξ, and a dense diagonalization costs ∝ (grid)³, so the cost ∝ tol^(−3/2). A looser tolerance is much cheaper.

**D2. The quantum Cv comes from the ξ = 1 solve.**
- **Where:** [`src/Quick_Scan.py`](src/Quick_Scan.py#L404), line 404.
- **Principle:** at ξ = 1, ξ²V = V: it is the physical spectrum, already computed by the scan.

**D3. A verdict uses a 2ε margin.**
- **Where:** [`src/Quick_Scan.py`](src/Quick_Scan.py#L423), lines 138 and 423.
- **Principle:** by the ħ² law, a coarse classical value sits *below* its limit by about ε. That pushes the quantum − classical difference up by about ε, so requiring |difference| > 2ε prevents false "quantum above classical" calls.

**D4. `RememberedPotential`.** *Numerics only.*
- **Where:** [`src/Quick_Scan.py`](src/Quick_Scan.py#L169), line 169.
- **What:** it caches repeated evaluations of V. Results are bit-identical, about 2.5× faster, and the engine is untouched.

**D5. T_merge: stop computing where quantum = classical.**
- **Where:** [`src/Quick_Scan.py`](src/Quick_Scan.py#L249), line 249. The starting temperature T0 is set at line 280, the acceptance test at line 298.
- **Principle (Wigner–Kirkwood):** Cv_q − Cv_cl = −β²·g″(β), where g = (ħ²/24m)·β²·⟨V″⟩.
  - For the harmonic oscillator this is −(ħω/k_BT)²/12, the first correction in Einstein's Cv.
  - For a V ~ |x|^k tail it falls as T^−(1+2/k), always faster than 1/T, so once that regime is reached the difference only shrinks.
- **The starting point:** T0 = 10·max(E₁ − E₀, ħω at the minimum), because quantum effects are controlled by ħω/k_BT.
- **More:** the sketch and its numerical check are in [`FINDINGS.md`](FINDINGS.md), "Why Quantum and Classical Cv Merge at High T".

## E. Checks added (results unchanged)

- **What:** [`audit/classical_limit/`](audit/classical_limit/) and [`audit/quick_scan/`](audit/quick_scan/) (the old-engine audit, verification of the corrected engine, Section 6 grid factors, quick-scan validation, the merge check).
- **Principle:** they compare against independent references:
  - the exact classical Cv, ½ + β²·Var(V) by numerical integration (option A, used only as a check);
  - the Pöschl–Teller potential, whose spectrum is exactly known, as a test that needs no DVR.

## F. The exact classical Cv, ½ + β²·Var(V): NOT used by the pipeline

**Status.** No code in `src/` computes this formula. The pipeline's classical curve always comes from the ξ-scan, i.e. from DVR spectra of ξ²V. The pipeline does rely on **one consequence** of it: the plateau test of A5 accepts a ξ step only if Cv ≥ ½ − tol (Step 5 below). That reliance is an open item, still to be verified and proven (see A5). The formula appears in `src/` only in two docstrings, which state it to justify that rule:
- [`src/Classical_Limit_Numerical.py`](src/Classical_Limit_Numerical.py#L62), line 62;
- [`src/Cv_Numerical_Benchmark.py`](src/Cv_Numerical_Benchmark.py#L292), line 292.

It is computed (by numerical integration) only in [`audit/`](audit/), as an independent check of the ξ-scan. Using it in the pipeline itself is "option A", which is still undecided.

### Part 1: Deriving Cv_cl = ½ + β²·Var(V)

**Setup.** One particle in 1-D, H = p²/2m + V(x), in contact with a bath at temperature T (canonical ensemble), with β = 1/(k_B T).

**Step 1: the classical partition function.**

  Z_cl = (1/2πħ) ∫∫ e^(−βH(x,p)) dx dp

The 1/2πħ counts one quantum state per phase-space area 2πħ (Step 7 shows where it comes from). It is a constant, and Step 4 shows it drops out.

**Step 2: it factorizes.** Because e^(−β(p²/2m + V)) = e^(−βp²/2m)·e^(−βV(x)), the two integrals separate:

  Z_cl = (1/2πħ) · [∫ e^(−βp²/2m) dp] · [∫ e^(−βV(x)) dx]

The p-integral is a Gaussian: ∫ e^(−βp²/2m) dp = √(2πm/β). Call the x-integral Q(β) = ∫ e^(−βV(x)) dx. Then

  Z_cl = (1/2πħ)·√(2πm)·β^(−1/2)·Q(β),  ln Z_cl = const − ½ ln β + ln Q(β)

**Step 3: the mean energy.** U = −∂(ln Z)/∂β. The constant drops out, −½ ln β gives 1/(2β), and

  −Q′/Q = ∫ V e^(−βV) dx / ∫ e^(−βV) dx = ⟨V⟩,

the average of V weighted by e^(−βV). So

  U = 1/(2β) + ⟨V⟩ = ½ k_B T + ⟨V⟩

The ½ k_B T is the kinetic energy (equipartition).

**Step 4: the heat capacity.** Cv = dU/dT, and since dβ/dT = −k_B β², we have Cv = −k_B β²·dU/dβ.
- The kinetic term: d/dβ [1/(2β)] = −1/(2β²).
- The potential term, differentiating the ratio:

  d⟨V⟩/dβ = −∫V² e^(−βV)/Q + (∫V e^(−βV)/Q)² = −(⟨V²⟩ − ⟨V⟩²) = −Var(V)

Putting them together:

  **Cv_cl/k_B = −β²·[−1/(2β²) − Var(V)] = ½ + β²·Var(V)**

No ħ is left: the 1/2πħ was a constant inside ln Z, and U and Cv only involve β-derivatives of ln Z.

**Step 5: consequences and checks.**
- **The bound:** a variance is never negative, so Cv_cl ≥ ½. This is the rule the ξ-scan uses (A5).
- **HO**, V = ½mω²x²: Q ∝ β^(−1/2), so ⟨V⟩ = 1/(2β) and Var(V) = 1/(2β²). Then Cv = ½ + ½ = 1, the equipartition value.
- **|x|^k** in general: Q ∝ β^(−1/k), so Var(V) = 1/(kβ²) and Cv = ½ + 1/k. That gives ¾ for x⁴, and ½ for the box (k → ∞, where V = 0 inside and Var = 0).
- **A handy identity:** Var(V) = ∂²(ln Q)/∂β². This is the "fluctuations = second derivative of ln Z" rule, applied to the x-part only.

### Part 2: Why this is the classical limit of the quantum Cv

**Step 6: the quantum partition function.**

  Z_q = Σ_n e^(−βE_n) = Tr e^(−βĤ) = ∫ dx ⟨x|e^(−βĤ)|x⟩

The pipeline's quantum Cv (k_B β²·Var(E)) comes from this Z_q.

**Step 7: ignore, for a moment, that x̂ and p̂ don't commute.** If e^(−β(T̂+V̂)) were equal to e^(−βV(x̂))·e^(−βp̂²/2m), then by inserting momentum states ⟨x|p⟩ = e^(ipx/ħ)/√(2πħ):

  ⟨x|e^(−βV)e^(−βp²/2m)|x⟩ = e^(−βV(x)) ∫ dp |⟨x|p⟩|² e^(−βp²/2m) = e^(−βV(x)) · (1/2πħ) ∫ dp e^(−βp²/2m)

Integrating over x gives exactly Z_cl from Step 1. This is where the 1/2πħ comes from: |⟨x|p⟩|² = 1/2πħ.

**Step 8: the error from that assumption vanishes as ħ → 0.** e^(−β(T̂+V̂)) and e^(−βV̂)e^(−βT̂) differ by terms built from the commutator [T̂, V̂], which is proportional to ħ (from [x̂, p̂] = iħ). Worked out properly (Wigner 1932, Kirkwood 1933), the result is

  Z_q = Z_cl·[1 − (ħ²β²/24m)·⟨V″⟩ + O(ħ⁴)]

So Z_q → Z_cl as ħ → 0, and with it Cv_q → Cv_cl. Physically, the correction is small when the thermal de Broglie wavelength ħ√(β/m) is small compared with the distance over which V changes. Equivalently, k_BT is much larger than the level spacing ħω.

**Step 9: the link to the ξ-scan.** By Eq. S7, the ξ-scan computes the quantum Cv of the same V at the same T, but with ħ/ξ. As ξ → ∞ that is ħ → 0, so by Step 8 the scan converges to ½ + β²·Var(V). The Wigner–Kirkwood term also shows that the remaining gap shrinks as ħ², i.e. as 1/ξ², which is the law behind the error estimate ε (A4).

**Conditions this needs (all true in this problem):**
- a single distinguishable particle (no Bose/Fermi statistics);
- a Hamiltonian of the form p²/2m + V(x), with the kinetic term a pure p², so the p-integral is Gaussian;
- the canonical ensemble;
- a confining V, so that Q(β) is finite;
- a smooth V. Hard walls change the correction to O(ħ) but not the limit.

**Numerical confirmation:**
- The ξ-scan converges to this formula for the HO, x⁴ and the double well, within 8.2·10⁻⁴ ([`audit/classical_limit/VERIFY_CORRECTED.txt`](audit/classical_limit/VERIFY_CORRECTED.txt)).
- The Wigner–Kirkwood term of Step 8 matches the DVR's Cv_q − Cv_cl to a ratio of 0.998–1.000 at high T ([`audit/quick_scan/MERGE_CHECK.txt`](audit/quick_scan/MERGE_CHECK.txt)).

**References for section F:**
- Steps 1–5:
  - Schroeder, *An Introduction to Thermal Physics*, ch. 6 (partition function, equipartition);
  - Pathria & Beale, *Statistical Mechanics*, ch. 3 (classical canonical ensemble, energy fluctuations);
  - [`audit/classical_limit/OPTION_A_physics.md`](audit/classical_limit/OPTION_A_physics.md).
- Steps 6–8:
  - Feynman, *Statistical Mechanics: A Set of Lectures*, ch. 2 (the Wigner correction);
  - Landau & Lifshitz, *Statistical Physics* Part 1, ch. III, "Expansion in powers of ħ";
  - Pathria & Beale, ch. 5 (density matrix and the classical limit);
  - E. Wigner, Phys. Rev. 40, 749 (1932); J. G. Kirkwood, Phys. Rev. 44, 31 (1933);
  - for the non-commuting exponentials: the Baker–Campbell–Hausdorff formula.

---

## Further reading

- **The ξ-scaling and classical limit:** Gelbwaser-Klimovsky et al., "Single-atom heat machines enabled by energy quantization", SI sections III–IV (in `docs/refs`).
- **Boltzmann factor, partition function, Cv = k_B·β²·Var(E), equipartition:** Schroeder, *An Introduction to Thermal Physics*, ch. 6; Pathria & Beale, *Statistical Mechanics*, ch. 3.
- **WKB and Bohr–Sommerfeld level counting:** Griffiths & Schroeter, *Introduction to Quantum Mechanics*, the WKB chapter.
- **Quantum corrections as a series in ħ (Wigner–Kirkwood):**
  - Landau & Lifshitz, *Statistical Physics* Part 1, ch. III, "Expansion in powers of ħ";
  - E. Wigner, Phys. Rev. 40, 749 (1932); J. G. Kirkwood, Phys. Rev. 44, 31 (1933);
  - a good undergrad exercise: expand Einstein's oscillator Cv at high T and recover −(βħω)²/12.
- **The error estimate:** "Richardson extrapolation" in any numerical-analysis text, e.g. Burden & Faires.
- **The DVR and grid resolution:** Colbert & Miller, J. Chem. Phys. 96, 1982 (1992); the Nyquist sampling theorem.
- **The exact classical Cv integral:** section F above (full derivation, not used by the pipeline), and [`audit/classical_limit/OPTION_A_physics.md`](audit/classical_limit/OPTION_A_physics.md).
- **The quantum → classical limit of the partition function:** Feynman, *Statistical Mechanics: A Set of Lectures*, ch. 2; Pathria & Beale, ch. 5; the Baker–Campbell–Hausdorff formula.
