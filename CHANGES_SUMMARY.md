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
- **The exact classical Cv integral:** [`audit/classical_limit/OPTION_A_physics.md`](audit/classical_limit/OPTION_A_physics.md).
