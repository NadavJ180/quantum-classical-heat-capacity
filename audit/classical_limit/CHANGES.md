# Correcting the classical limit: what changed and why

Branch `correct-classical`. This document is written to be followed step by step:
- section 1 restates the check you asked for;
- section 2 explains the new method in words;
- section 3 walks through every changed file;
- section 4 lists what you will see differently;
- section 5 gives the evidence it works;
- section 6 lists what was deliberately left untouched.

The physics of option A (the exact classical integral) is explained separately in [OPTION_A_physics.md](OPTION_A_physics.md); the original diagnosis is in [AUDIT.md](AUDIT.md).

---

## 1. The double-check: is the potential ever scaled by ξ?

**No, in no file and at no point in the history.** Checked four ways on the code as it was (commit `1ccfb7d`):

1. **Arithmetic involving ξ.** Every such line in every `.py` file (src, homog, scripts) is one of five lines in `Classical_Limit_Numerical.py`: 69 and 78 (`compute_cv`), 159 (`xi *= xi_mult`), 238 and 250 (`converge_n`). All five act on the energies and β.
2. **Potential evaluations.** Every evaluation of a potential, including the one place the DVR builds V on its grid (`DVR_Algorithm.py`, `v_diag = potential_func(x)`), receives the unscaled `my_potential`, or the coefficient-sweep variant built with `functools.partial`. No expression multiplies or divides a potential by anything.
3. **ħ.** No call passes an ħ other than 1.0.
4. **Git history, all branches.** The only past `* xi**2` is the early box code's `beta = 1.0 / (T * kB * xi**2)` (commits `bc27e6c`, `c55b891`), which again scales temperature only.

So the classical limit has always been computed with T → ξ²T alone, which is the defect described in `AUDIT.md`.

## 2. The new method, in words

The pipeline now follows the paper's prescription literally: scale **both** the potential and the temperature by ξ² (Gelbwaser et al., main text, SI-III item 2, SI-IV; Eq. S7 makes this ħ → ħ/ξ). At every temperature T:

1. **Walk a ladder of ξ values,** 1, 1.25, 1.56, … (`XI_START · XI_MULT^k`).
2. **At each ξ, solve the Schrödinger equation again** with the DVR for the scaled potential ξ²V. This is the step that was missing. Then evaluate `compute_cv(E_n(ξ²V), β, ξ)`:
   - `compute_cv`'s weights exp(−βE/ξ²) are the T → ξ²T half;
   - the re-solved spectrum is the V → ξ²V half.
   - By Eq. S7, the result is the quantum Cv of the original system with ħ/ξ at the original T.
3. **Decide when the plateau is reached, with an error estimate.**
   - Quantum corrections to the classical Cv shrink as ħ², i.e. as 1/ξ².
   - On a ladder with ratio m = 1.25, the remaining distance to the ξ → ∞ limit is therefore |Cv(ξ_k) − Cv(ξ_{k−1})| / (m² − 1).
   - A step is *stable* when this estimated distance is below `TOL_XI` **and** Cv ≥ ½ − `TOL_XI`. The second condition uses the rigorous bound C_v^cl = ½ + β²Var(V) ≥ ½. It stops the frozen-out regime, where Cv ≈ 0 for several steps at very low T, from being mistaken for a plateau.
   - Three consecutive stable steps (`MIN_STABLE_XI`) make a verified plateau.
4. **Report the last point of the verified plateau,** the most classical value computed, together with its estimated distance to the limit.
5. **If the scan fails, report NaN.** This happens if the ladder runs out (`max_steps`), hits the hard cap `XI_MAX` (`xi_cap`) or a DVR solve fails (`dvr_failed`). It is never a substitute value.

**Making it affordable.**
- The spectrum of ξ²V does not depend on T, so one solve per ξ serves every temperature. Solves are cached.
- Each solve keeps only the levels needed: up to E₀ + `HOT_STATE_SAFETY`·k_BT in units of V, for the hottest temperature that uses it. Temperatures are processed hot → cold, so a solve made for one temperature already covers every colder one.
- Because every solve keeps the levels the temperature needs, the old "finite-N collapse" cannot happen any more. The companion n-scan still re-checks this.

## 3. File by file

### `src/Classical_Limit_Numerical.py` (rewritten)

- **`compute_cv`:** unchanged. Only its docstring now says it is the classical-limit building block *only* when handed the spectrum of ξ²V.
- **`ScaledSpectra` (new):** the cache of DVR solves of ξ²V. For each (ξ, T) it:
  1. estimates the level count with a semiclassical (WKB) count at ħ/ξ, which by Eq. S7 equals the count for ξ²V;
  2. calls the pipeline's own `auto_configure_dvr` and `get_fully_converged_energy_levels` (3-pass resolution/span check);
  3. checks the coverage and grows the level count if it falls short.
  - Both DVR tolerances are multiplied by ξ², so they mean the same thing, 1e-5 in units of V, at every ξ.
  - For Section 6's reference it instead widens/refines each auto-configured grid by `span_factor`/`dx_factor`, as `DVR_Reference_Generator` does, and solves once.
- **Helpers `_allowed_region`, `_wkb_level_count`, `_wkb_energy_for_count` (new):** the potential's classically allowed window and WKB level counting, used only to size each solve.
- **`converge_xi`:** rewritten as in section 2. Now takes the cache and the ladder instead of a fixed spectrum. Returns the plateau's indices and the error estimate. **S1 is fixed by construction:** the value can only come from the verified streak, and there is no collapse tail to pick from.
- **`converge_n`:** unchanged. It now runs on the converged ξ's own spectrum.
- **`sweep_temperature_range`:** new signature, taking the potential instead of energies plus `mass`, `hbar`, `thermal_coverage`, `xi_max`, `spectra`, `span_factor`, `dx_factor`. **S2 is fixed:** no fallback to `converge_n` at ξ = 1 when the ξ-scan fails. Also returns `n_available`, `error_estimate`, `xi_ladder`, `spectra`.

### `src/DVR/DVR_Algorithm.py` (three optional arguments, defaults = old behavior)

`auto_configure_dvr` gained `energy_ceiling=None`, `turning_points=None` and `padding=2.0`.

- **`energy_ceiling`:** its built-in ceiling `v_min + 1.5·num_levels` assumes level spacings of order 1. For ξ²V they are ~ξ times larger, so the scaled solves pass a WKB-based ceiling.
- **`turning_points`, `padding`:** at large ξ the wavefunctions occupy a far narrower region than the fixed 2.0-unit padding allows for. Starting there made the grid ~3× larger than needed (5,094 instead of ~2,000 points at ξ = 517), and the cost grows as the cube.
- **Safety is unchanged:** stage 2 of `auto_configure_dvr` still widens the span until the energies stop changing, and `get_fully_converged_energy_levels` still validates every solve.
- Every existing caller passes none of the new arguments, so its behavior is unchanged.

### `src/Quantum_Classical_Combined.py`

- **`run(energies, potential_func, system_name, …)`:** the classical limit needs the potential. New arguments `mass`, `hbar`, `thermal_coverage`, `xi_max`, `spectra`; the result also returns `spectra` for reuse.
- **The quantum Cv now uses the whole base spectrum.** It used to be cut at the largest n the classical n-scan needed, but that n now refers to the scaled spectra. Using the whole spectrum is exact while the top level is thermally inaccessible, which the auto-tune loop enforces.
- **ξ-diagnostic plot:** x-axis now logarithmic. Points are colored *verified plateau* / *stable step* / *still moving* instead of the obsolete "falling (finite-N collapse)", and the annotation shows the estimated distance to the limit.

### `src/Cv_Numerical_Benchmark.py` (Section 6)

- **The classical "reference" is now meaningful.** It is the same ξ-scan with every ξ²V solve on a grid widened/refined by Section 2's factors (span×2, dx÷2). It previously re-ran the identical algorithm on the reference spectrum.
- **Signature:** `run_cv_numerical_benchmark(…)` gained `potential_func`, `mass`, `hbar`, `thermal_coverage`, `xi_max`, `span_factor`, `dx_factor`.
- **Classical plot:** the fixed y-limit of 1.1 would clip the double well's classical curve, which peaks at 1.45, so it now scales to the data. A docstring that said the classical curve "falls toward zero at cold T" was corrected; it never does, since it is ≥ ½.

### `src/Cv_AutoTune.py`

`diagnose_escalation` now maps failures to the knob that actually controls them:
- **`need_states`:** only when the *base* spectrum's top level is thermally accessible at the hottest T. NUM_STATES now only affects the quantum curve.
- **`need_xi`:** when more than 5% of temperatures ran out of ladder (`max_steps`).
- **Reported, not escalated:** `xi_cap` and `dvr_failed` (with a warning from the master script), and `n_marginal_frac`.

The old hot-half/cold-half split and the `finite_n` signal no longer apply.

### `src/Quantum_HO_Master.py`

- Passes the potential and the new settings to `run()` and to Section 6.
- Reuses the scaled-spectrum cache across auto-tune rounds, so a round that only grows NUM_STATES repeats no classical solve.
- New warning if any temperature stopped at `XI_MAX` or on a failed DVR solve.

### `src/config.py`

| setting | was | now | why |
|---|---|---|---|
| `XI_START` | 3.0 | 1.0 | start from the physical system; at high T the plateau is reached by ξ ≈ 2, and starting at 3 tripled the levels needed there |
| `XI_MULT` | 1.1 | 1.25 | every rung is now a DVR solve; ×1.25 reaches ξ ≈ 1,300 in ~33 solves |
| `TOL_XI` | 5e-3 | 2e-3 | **meaning changed:** now the estimated distance of the reported value from the ξ → ∞ limit, not the per-step change |
| `MIN_STABLE_XI` | 3 | 3 | — |
| `MAX_XI_STEPS` | 80 | 35 | ladder up to 1.25³⁴ ≈ 1,970 |
| `XI_MAX` | — | 2000 (new) | hard cap that auto-tune never raises (grid size grows ~ξ) |

### `src/analytical/HO_Analytical.py`

Docstring only. It said the classical Cv is constant "by definition of being the classical limit — what the quantum Cv approaches as T → ∞". That is true for the HO but false in general, and it is exactly the assumption the old scan relied on.

### `audit/` (not part of the pipeline)

- `audit_classical_limit.py` now loads the **pre-fix** engine from git (commit `1ccfb7d`) with that commit's settings, so it keeps auditing the old code on this branch.
- `verify_corrected_engine.py` (new) checks the corrected engine against exact answers (section 5).
- `OPTION_A_physics.md` (new) and this file.

## 4. What you will see differently

- **The classical curve is the classical curve.** For the base double well it runs 1.00 → 1.45 (T ≈ 0.55) → 0.71, instead of a flat ~0.73.
- **No fake classical values.** A temperature where the scan cannot converge shows a gap (NaN) and a console warning, instead of the quantum Cv silently standing in.
- **`ξ_conv` now means something physical:** the effective ħ at the plateau was ħ/ξ_conv. It grows like 1/T, reaching ~1,262 at T = 0.022 for the double well.
- **Runtime.** Sections 1 & 4 take ~7 min for the double well (33 DVR solves of ξ²V). Section 6's reference sweep takes ~36 min (section 5). The rest of the pipeline is unchanged; the whole run took ~52 min.

## 5. Evidence

### 5a. The corrected engine against exact answers (`verify_corrected_engine.py` → `VERIFY_CORRECTED.txt`)

This uses the pipeline's own `sweep_temperature_range` with the `config.py` settings, 40 log-spaced temperatures per case:

| system | exact classical Cv | converged | max error | error estimate vs true error | ξ_conv range | time |
|---|---|---|---|---|---|---|
| HO, x²/2 | 1 (all T) | 40/40 | 7.8e-4 | agree to 6e-7 | 1.95 – 517 | 191 s |
| quartic, x⁴ | ¾ (all T) | 40/40 | 8.2e-4 | agree to 3e-6 | 1.95 – 136 | 36 s |
| **your double well** | phase-space integral (1.00 → 1.45 → 0.71) | 40/40 | **8.2e-4** | agree to 6e-6 | 1.95 – 1262 | 345 s |

Every error is below `TOL_XI` = 2e-3, and every value approaches the limit from below by almost exactly its own error estimate. For comparison, the old engine missed the double well by up to **0.72** (AUDIT.md, T9).

### 5b. The full pipeline, end to end (`corrected_pipeline_run.txt`, figures in `figures/corrected_run/`)

`Quantum_HO_Master.py` was run unmodified on the double well, with figures redirected to a scratch folder and copies kept here. Exit code 0, no warnings.

- **Sections 1 & 4:** one auto-tune round (no escalation). ξ and n converged at **1000/1000** temperatures, using 33 DVR solves of ξ²V (ξ up to 1,262 at T = 0.022, largest grid 3,197 points). 399 s.
- **Sections 2, 3, 5:** unchanged by this work, and still passing.
  - Section 3: base vs reference levels agree to 1.4e-9.
  - Section 5: maximum dx and trustworthy n found as before.
- **Section 6:**
  - Classical limit, base vs reference grids (span×2, dx÷2): max relative error **5.5e-12** over 1000/1000 temperatures, so the new classical curve is fully converged in the DVR grid.
  - Quantum Cv: max absolute error 8.5e-12.
  - Cost: **36 minutes**, because the reference's 33 solves reach 12,785 grid points. This is now the most expensive part of the run (see section 6 of this file).
- **Section 7:** ran as before in 67 s. See the caveat in section 6 of this file.

What the figures show:
- **`cv_summary.png`:** the classical curve is now 1.00 → 1.45 → 0.71, and the quantum curve lies below it at every T.
- **`xi_convergence.png`:** the scan at the coldest temperature. Frozen out (Cv ≈ 0, correctly rejected by the ½ bound), then rising, then a verified plateau at 1.0043.

### 5c. The audit still audits the old code

`audit_classical_limit.py` now loads the pre-fix engine from commit `1ccfb7d`. Re-run on this branch, its `RESULTS.txt` is **identical** to the original audit's, apart from the header and timing. This also shows that the three optional `auto_configure_dvr` arguments leave existing behavior unchanged: tests T3 and T5 call it without them and reproduce their numbers exactly.

### 5d. One possible refinement (not implemented)

Because the reported value sits below the limit by almost exactly its estimate, adding the estimate back (Richardson extrapolation) would give the classical Cv to ~1e-5 at the same cost. It would also allow a looser `TOL_XI` and therefore fewer DVR solves. It changes the reported number from "computed" to "extrapolated", so it is your call.

## 6. Deliberately not changed, and needs your decision

- **Section 7 (coefficient sweep) is now misleading in a new way.**
  - It still plots every variant against the **base** potential's classical curve, as its title says ("shared classical limit reused from the base run"). That curve is now correct for b = −0.5.
  - The b = −0.9 and −0.7 quantum curves therefore appear to rise above "the classical limit" at T ≈ 1.3–10 (`figures/corrected_run/cv_coefficient_sweep.png`). That is a comparison between different potentials. Each b has its own, quite different, classical curve. For example, the earlier quick calculation gave b = −0.9's own classical Cv as ≈ 1.36 at T = 3, above its quantum peak of ≈ 1.28.
  - **Do not read that plot as a quantum excess.**
  - Fixing it means one classical curve per variant: about 6 minutes each with the ξ-scaling, or about a second each with option A. I have not changed it.
- **Section 6 cost.** The classical reference sweep takes ~36 minutes for the double well. Options:
  - keep it;
  - lighten it (e.g. reference factors 1.5 instead of 2 for the classical part only);
  - make it optional with a config flag, relying on each ξ solve's own 3-pass check (or on option A) instead.
- **Documentation outside `src/`:**
  - `README.md`, `FINDINGS.md` (section "The Classical Limit and the ξ-Scaling Trick", the "ξ-convergence diagnostic" row of the diagnostics table, and the Troubleshooting entries about `XI_START` and finite-N collapse) and `HISTORY.md` describe the old method.
  - `docs/summaries/IEEE_Summary.tex` lines 77, 182 and 198 do too. It is synced with Overleaf, and you write it yourself.
- **`verification/classical_limit/`** (the paused toy-model check) was written against the old API and has not been updated.
- **Figures** under `figures/` were not regenerated. The verification run wrote its figures to a scratch folder so your committed figures stay as they are.
