# Audit: the pipeline's classical Cv limit

Branch `correct-classical`, 2026-10-08. This audits the code as it was on `main` at commit `1ccfb7d`. Every number below comes from `audit_classical_limit.py` in this folder (full output in `RESULTS.txt`, figures in `figures/`). The script loads that commit's `Classical_Limit_Numerical.py` from git, with that commit's ξ settings, so it keeps testing the old engine even on this branch.

**Status:** the defect has since been fixed on this branch. See [CHANGES.md](CHANGES.md) for what changed and the evidence, and [OPTION_A_physics.md](OPTION_A_physics.md) for the physics of option A.

## Verdict

**The classical limit is computed wrong for any potential whose energy levels do not scale homogeneously, which includes the project's double well.**

- The pipeline applies only the temperature half of the ξ-scaling prescription, T → ξ²T, and never rescales the potential.
- Its "classical limit" is therefore the quantum Cv of the same system at a much higher temperature, so it is nearly the same value at every T.
- For the harmonic oscillator, the infinite well and any |x|^k potential, this shortcut is mathematically identical to the correct prescription. Those are the only systems it was ever checked on, which is why no check ever failed.

Two smaller, independent defects in the same code are listed in section 6.

## 1. What the code does

There is a single classical-limit engine, `Classical_Limit_Numerical.sweep_temperature_range`. Section 4 (`Quantum_Classical_Combined.run`), Section 6 (`Cv_Numerical_Benchmark`) and Section 7 (`Cv_Coefficient_Sweep`, which reuses the base curve) all consume its output.

- **The scan:** at each T, `converge_xi` loops over ξ and calls `compute_cv(energies, beta, xi)` (`Classical_Limit_Numerical.py:133`). It always passes the same `energies`: the spectrum of the unscaled potential V.
- **The function:** `compute_cv` weights each level by exp(−βE/ξ²) (line 69) and multiplies by β²/ξ⁴ (line 78).
- **Together:** compute_cv(E, β, ξ) = compute_cv(E, β/ξ², 1), the ordinary quantum Cv of V at temperature ξ²T.

The scan heats the system. It never computes the spectrum of a rescaled potential.

## 2. What the physics requires

Gelbwaser-Klimovsky et al. (docs/refs, *Single-atom heat machines enabled by energy quantization*) give the prescription in the main text, in SI-III item 2 and in SI-IV: "increase the potentials and temperatures by the same multiplicative factor ξ²", that is, V → ξ²V and T → ξ²T. Their Eq. S7,

  E_n(ħ, ξ²V) = ξ² E_n(ħ/ξ, V),

shows that this is ħ → ħ/ξ at fixed V and T. Its left-hand side is the spectrum of the **re-solved, rescaled potential** ξ²V.

The pipeline uses E_n(ħ, V) in place of E_n(ħ, ξ²V). Equivalently, it assumes E_n(ħ/ξ, V) = E_n(ħ, V)/ξ². That is the infinite-square-well rule E_n ∝ ħ², where V = 0 inside the well and rescaling V changes nothing. The method was first written for exactly that system: commit `e9fa799`, "classical limit of the box potential", `E_n = n**2 * E_g / xi**2`. Commit `cb09603` (2026-06-13) later generalized it into today's `compute_cv(energies, beta, xi)` / `converge_xi`, which divide any spectrum by ξ².

The report states the method the same way (`docs/summaries/IEEE_Summary.tex:77`: "realized numerically by scaling both temperature and the spectrum by ξ²"). Scaling the spectrum of V is not the same as using the spectrum of ξ²V.

## 3. When the shortcut is exact: the condition

**Claim.** Suppose every gap scales by one common factor, (E_n − E₀)[ξ²V] = q(ξ)·(E_n − E₀)[V] for all n.

**Consequence.** The correct ξ-scan at ξ then equals the pipeline's scan at ξ_p = ξ/√q(ξ).

**Proof.** Cv only depends on gaps. Substituting q·E[V] for E[ξ²V] in compute_cv(E[ξ²V], β, ξ) gives
- weights exp(−βqE/ξ²), and
- prefactor β²q²/ξ⁴ multiplying Var(E[V]).

That is exactly compute_cv(E[V], β, ξ_p) with ξ_p² = ξ²/q. ∎

**The homogeneous cases** (with E_n ∝ ħ^(2k/(k+2)) for |x|^k):

| system | q(ξ) | pipeline's ξ_p |
|---|---|---|
| HO | ξ | √ξ |
| box | 1 | ξ |
| \|x\|^k | ξ^(4/(k+2)) | ξ^(k/(k+2)) |

In every case ξ_p → ∞ as ξ → ∞, so the shortcut reaches the correct plateau.

**When q depends on n**, no relabelling exists, and the shortcut converges instead to the high-T limit of the quantum Cv at every T. This is the paper's own distinction between homogeneous and inhomogeneous level scaling. It is also what `homog/RESULTS.md` found for this potential family: only (A/g⁶, B/g⁵, C/g⁴) rescalings scale homogeneously, and ξ²·(A, B, C) is not of that form.

## 4. Evidence

| test | what it checks | result |
|---|---|---|
| T1 | runtime trace of `compute_cv` inside a real `sweep_temperature_range` | 164 calls, 49 distinct ξ, **1 distinct spectrum** (the unscaled input) |
| T2 | compute_cv(E, β, ξ) = compute_cv(E, β/ξ², 1) | 5e-12 on 1,000 random cases; equals the repo's Einstein Cv at ξ²T to 5e-13 |
| T3 | Eq. S7 on the repo's DVR | E(ξ²V)/ξ² vs E(ħ/ξ, V): 1e-12. Pipeline's E(V)/ξ² is off by up to 93 (double well, ξ = 2: E₀ = −0.25 vs −1.46) |
| T4 | gap ratio at ξ = 2, n = 1..20 | HO 2.000000; box 1.000000; x⁴ 1.587401 (= 2^(2/3)); Pöschl–Teller 1.11–1.55; **double well 1.55–2.86** |
| T5 | HO relabelling identity | pipeline(ξ) vs correct(ξ²): 7e-14 |
| T6 | box, exact levels | pipeline ≡ paper's prescription by construction; error on genuine plateaus 0.008 |
| T7 | x⁴ (homogeneous, anharmonic) | passes, 4e-4 vs 3/4 |
| T8 | Pöschl–Teller V₀ tan²x: exact levels for every ξ²V; analytic, T-dependent Cv_cl (0.98 → 0.56) | **pipeline fails by 0.47 on genuine plateau points** (flat ≈ 0.51). Correct scan, same `compute_cv`: 0.21, 0.043, 0.010, 2.5e-3, 6.1e-4, 1.5e-4, 3.8e-5 for ξ = 1…64 (4× per doubling, the ħ² law) |
| T9 | the project's double well, pipeline's own T grid | pipeline flat 0.728–0.740; exact Cv_cl 0.713–1.449. **Error up to 0.72**, on 1000/1000 genuine plateau points. Correct scan converges: 0.30, 0.098, 0.026, 0.0067 over T ≥ 1 for ξ = 1, 2, 4, 8 |

Independent cross-checks behind the reference values:

- **Pöschl–Teller spectrum formula:** checked against a finite-difference solve, to 3e-6.
- **Pöschl–Teller analytic Cv_cl:** checked against direct numerical quadrature, to 7e-13.
- **Exact classical integral used for the double well:** reproduces the HO (1) and x⁴ (3/4) limits to 1e-16, and is grid-converged to 1e-15.

## 5. Why the HO validation could not catch it

1. **For the HO the pipeline's algorithm is not approximately right but exactly right.** By section 3 it is the correct ξ-scan with ξ read as ξ². T5 confirms this numerically to 7e-14, so no HO test of any precision could expose the defect.
2. **The HO's classical Cv is T-independent** (= 1). Any method that returns the high-temperature value passes, and that is exactly what the defect produces for every potential.
3. **The "machine precision" for the classical limit is a self-comparison.** `fig_cv_benchmark_classical.png` (rel. error 6e-14; reproduced here as 2e-13) compares the base DVR and the reference DVR running the same ξ algorithm. It measures reproducibility, not correctness. The report itself notes this limitation for shared flaws (`IEEE_Summary.tex:146`), and the analytic check meant to close the gap only covered the quantum curve.
4. **The analytic HO benchmark never quantified the classical limit.** `HO_Benchmark.py` computes errors for the quantum Cv only and draws k_B as a reference line (`HO_Benchmark.py:126`). It has not been called since commit `705c1f9`; `Quantum_HO_Master.py:195` passes `cv_analytic=None`.
5. **Every system the method was ever checked on is in the blind class:** box (where it was written), HO, and implicitly power laws. The early Morse script drew its "classical limit" as the constant 1 and never ran it through the ξ method. A temperature-dependent classical reference was never computed for any system.
6. **The double-well runs raise no internal warning.** All 1,000 temperatures land on a genuine plateau (it is a real plateau, of the quantum Cv at high T); the n-scan converges; base and reference agree; the auto-tune loop ends clean. Every diagnostic tests numerical consistency of the algorithm, none tests its physics.
7. **Extra finding.** In the reproduced HO run, 90 of the 500 "classical limit" points (the hot end, T ≳ 2.5) are not ξ-scan results at all. They are the quantum Cv itself, via the fallback in section 6 (S2), which is ≈ 0.99 there and therefore indistinguishable from 1.

## 6. Two secondary defects in the same engine

Neither is triggered in the project's double-well run (1000/1000 genuine plateaus), but both corrupt results elsewhere.

- **S1, plateau picker (`Classical_Limit_Numerical.py:168-172`).** After a finite-N collapse, the search takes the *last* step with |ΔCv| < `TOL_XI`. Because that tolerance is absolute (0.005), the step can sit in the decaying collapse tail and be reported as "converged". Example: box, T = 12.8, reported Cv = 0.011. This happened at 30/300 points for the box and 16/300 for Pöschl–Teller.
- **S2, fallback (`Classical_Limit_Numerical.py:343-353`).** When the ξ-scan fails, `converge_n` is run at ξ = 1 and its value is stored as the classical limit. That value is the quantum Cv at T. This happened at 90/500 points for the HO and 50/300 for the box (cold end, reported as 0).

## 7. What is and is not affected

**Affected:**
- Every "numerical classical limit" curve for non-homogeneous potentials:
  - Section 4: `cv_summary`.
  - Section 6: the classical benchmark, which compares the defect with itself.
  - Section 7: the shared reference curve in `cv_coefficient_sweep.png`.
- `Cv_AutoTune`'s ξ-related escalation (`need_xi`, `finite_n` fractions), since it tunes this scan.
- Report and documentation claims:
  - `IEEE_Summary.tex:77` (method), `:182` (the plateau as an equipartition value), `:198` ("machine precision ... for both the quantum curve and its classical limit");
  - the FINDINGS.md section "The Classical Limit and the ξ-Scaling Trick";
  - any conclusion drawn from "quantum Cv above the classical limit" for the double well.

**Not affected:**
- DVR spectra (3-pass converged; Eq. S7 holds on them to 1e-12).
- Quantum Cv(T) curves (HO vs Einstein: 3e-13).
- `compute_cv` itself, which is correct when handed the spectrum of ξ²V (T8, T9).
- `homog/` (quantum only).

**Consequence for the toy model (base potential b = −0.5 only).** Against its exact classical Cv, the quantum Cv **never exceeds** the classical curve: max(Cv_q − Cv_cl) = −0.0014, at the hottest T. The apparent excess of up to +0.13 comes entirely from the wrong reference. The other b values have not been audited here.

## 8. Options for a fix (proposals only, nothing implemented)

- **A. Exact phase-space integral (recommended as the reference).** Cv_cl = ½ + β²·Var_β(V) with weights e^(−βV(x)). It is exact for any 1-D H = p²/2m + V(x), uses one 1-D integral per T, and involves no spectrum, truncation or plateau detection.
- **B. The paper's ξ-scaling, done as written.** Re-solve the DVR for ξ²V at each ξ; each solve is reused for every T. Levels needed grow ∝ ξ: about 225, 350, 600 and 1,100 levels for ξ = 1, 2, 4, 8 here. Keep plateau detection over ξ. This keeps the ξ narrative of the report (an experimentally realizable route to ħ → 0) and is a valuable cross-check of A.
- **C. Independently of A and B:**
  - make the plateau picker use the streak that `converge_xi` already verified (or a relative tolerance), fixing S1;
  - return NaN instead of the quantum Cv when the scan fails, fixing S2.

## Reproduce

```bash
cd audit/classical_limit
python -u audit_classical_limit.py
```

Runtime is about 5 minutes. It writes `RESULTS.txt` and `figures/`, and needs the Anaconda Python 3.7 environment the pipeline uses.
