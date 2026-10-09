# Correcting the classical limit: what changed and why

Branch `correct-classical`. This document is written to be followed step by step:
- section 1 restates the check you asked for;
- section 2 explains the new method in words;
- section 3 walks through every changed file;
- section 4 lists what you will see differently;
- section 5 gives the evidence it works;
- section 6 lists what was deliberately left untouched;
- sections 7 and 8 are the follow-up rounds: per-variant classical curves and the documentation (7); Section 6's lighter classical reference and the quick initial scan (8).

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
- **Section 7:** ran in 67 s, at that point still against the base potential's classical curve (`figures/corrected_run/cv_coefficient_sweep.png`, kept as a record). Superseded by per-variant classical curves (section 7a).

What the figures show:
- **`cv_summary.png`:** the classical curve is now 1.00 → 1.45 → 0.71, and the quantum curve lies below it at every T.
- **`xi_convergence.png`:** the scan at the coldest temperature. Frozen out (Cv ≈ 0, correctly rejected by the ½ bound), then rising, then a verified plateau at 1.0043.

### 5c. The audit still audits the old code

`audit_classical_limit.py` now loads the pre-fix engine from commit `1ccfb7d`. Re-run on this branch, its `RESULTS.txt` is **identical** to the original audit's, apart from the header and timing. This also shows that the three optional `auto_configure_dvr` arguments leave existing behavior unchanged: tests T3 and T5 call it without them and reproduce their numbers exactly.

### 5d. One possible refinement (not implemented)

Because the reported value sits below the limit by almost exactly its estimate, adding the estimate back (Richardson extrapolation) would give the classical Cv to ~1e-5 at the same cost. It would also allow a looser `TOL_XI` and therefore fewer DVR solves. It changes the reported number from "computed" to "extrapolated", so it is your call.

## 6. Open items

- **Section 6 cost.** Measured in section 7c, applied in section 8a.
- **`docs/summaries/IEEE_Summary.tex`** lines 77, 182 and 198 still describe the old method. It is synced with Overleaf, and you write it yourself.
- **`verification/classical_limit/`** (the paused toy-model check) was written against the old API and has not been updated.
- **Figures** under `figures/` were not regenerated. The runs wrote their figures to a scratch folder, and copies are kept in `figures/corrected_run/` here.

## 7. Follow-up: per-variant classical curves, documentation, Section 6 factors

### 7a. Section 7 compares each variant with its own classical limit

**`src/Cv_Coefficient_Sweep.py`**
- **Computation.** `run_coefficient_sweep` now computes every variant's classical limit with the same corrected ξ-scan and the base run's settings (new helper `_variant_classical_limit`). It takes new arguments `xi_start`, `tol_xi`, `min_stable_xi`, `xi_multiplier`, `max_xi_steps`, `xi_max`, `tol_cv`, `min_stable_n`.
  - The base variant reuses the base run's curve, since it is the identical computation.
  - A variant whose classical scan fails is kept, with gaps in its curve and a † in the legend, instead of aborting the sweep.
  - The result also returns `variant_classical` and `classical_incomplete`.
- **Plot.** `plot_coefficient_sweep` draws each variant's quantum $C_v$ solid and its own classical limit dashed, in the same color. The legend sits outside the axes, with one entry per variant plus two neutral entries for the line styles. The title states the convention.
- **Docstring.** The module docstring's "why no ξ/n search here" section became "each variant gets its own classical limit".

**`src/Quantum_HO_Master.py`** passes the base run's final ξ settings to Section 7 and updates the Section 7 comments. All four of its calls into the changed functions were checked to bind to the new signatures.

**Run** ([`SECTION7_RESULTS.txt`](SECTION7_RESULTS.txt), [`figures/corrected_run/cv_coefficient_sweep_own_classical.png`](figures/corrected_run/cv_coefficient_sweep_own_classical.png)): Sections 1 & 4 + Section 7 exactly as the master calls them, 1000 temperatures, ~44 min.
- All six variants converged at all temperatures.
- **For every b the quantum $C_v$ stays below its own classical $C_v$ at every temperature** (closest at the hottest T, by about 1e-3).
- The quantum peaks of b = −0.9 and −0.7 (1.28 and 1.13) sit below classical peaks of 1.42 and 1.43.
- b = −0.9 needed ξ up to ≈1972, the last rung of the default ladder. It converged, but with no margin.

### 7b. Documentation

- **`README.md`:** rewritten to describe only the current pipeline: the classical-limit method, the new auto-tuning logic, Section 6's classical reference, Section 7's per-variant curves, runtimes, and running a new potential.
- **`FINDINGS.md`:** updated physics, diagnostics, parameters and troubleshooting.
  - The classical limit is now explained from its definition and the exact phase-space formula; the ξ-scaling, the ħ² error estimate and the ½ bound are covered.
  - Results now include the classical-limit verification, the base double well, the per-variant table and the Section 6 factor measurements.
  - The ξ-diagnostic and coefficient-sweep plot descriptions match the new figures.
  - The troubleshooting covers the NaN stop reasons (`max_steps`, `xi_cap`, `dvr_failed`).
- **`HISTORY.md`:** a new chapter tells the October 2026 process: the toy-model verification plan, the first check, the audit and why the HO validation could not catch the defect, the correction, its verification, the per-variant Section 7 and its result, and the Section 6 factor check. Earlier entries are kept as history, annotated where the correction made them obsolete. The version table has new rows for every changed module.
- **Not touched:** the `.tex` files.

### 7c. Can Section 6 use lighter grid factors?

**Yes.** Measured on the full 1000-temperature grid ([`section6_grid_factor_check.py`](section6_grid_factor_check.py), [`SECTION6_GRID_FACTORS.txt`](SECTION6_GRID_FACTORS.txt), [`figures/section6_grid_factors.png`](figures/section6_grid_factors.png)):

| span / dx factor | reference sweep | largest grid | max relative difference to base | mean |
|---|---|---|---|---|
| 2 / 2 (current) | 36.3 min | 12,785 pts | 5.5e-12 | 2.4e-13 |
| 1.5 / 1.5 | 9.8 min | 7,192 pts | 5.4e-12 | 2.0e-13 |
| 1.25 / 1.25 | 5.7 min | 4,995 pts | 6.9e-12 | 1.9e-13 |

**Interpretation.**
- **All three show the same thing.** Every reference agrees with the base curve at the round-off level, ~1e-12 at the cold end, where energies of ξ²V are ~10⁶ times larger, and ~1e-13 at the hot end. Every factor shows the same fact: the base classical curve is converged in the grid.
- **The figures would not change.** The curves in Section 6's figure are indistinguishable, and its error panel would show the same ~1e-12 profile.
- **A lighter reference is still a real test.** It still changes both the span and the spacing of every solve. Because DVR errors fall exponentially with resolution and span, an under-resolved base solve would still show up as a visible difference.

**Suggested change.** Add a separate pair `CLASSICAL_REFERENCE_SPAN_FACTOR`/`CLASSICAL_REFERENCE_DX_FACTOR` = 1.5/1.5, a 3.7× speed-up with comfortable margin, or 1.25/1.25 for 6.4×. Keep Section 2's 2/2 for the spectrum reference, which costs under a minute. *Applied with 1.5/1.5 and a precise 2/2 option in section 8a.*

## 8. Follow-up: Section 6's lighter classical reference, and a quick initial scan

### 8a. Section 6's classical reference on 1.5 / 1.5, with a precise option

**`src/config.py`**
- New `CLASSICAL_REFERENCE_SPAN_FACTOR` = 1.5 and `CLASSICAL_REFERENCE_DX_FACTOR` = 1.5: the grid of every ξ²V solve in Section 6's classical reference. The comment records the section 7c measurement.
- New `CLASSICAL_REFERENCE_PRECISE` = `False`. `True` gives the classical reference Section 2's own factors (2 / 2): the stricter, slower option.
- New helper `reference_label(span_factor, dx_factor)`. `ref_label` is now built with it, with the same text as before.

**`src/Quantum_HO_Master.py`** (Section 6 only)
- Chooses the classical reference factors: `CLASSICAL_REFERENCE_*`, or Section 2's (`reference_result["span_factor"]`, `["dx_factor"]`, which respects interactive mode) when `CLASSICAL_REFERENCE_PRECISE` is on.
- Passes them, and a label naming them, to `run_cv_numerical_benchmark`.
- Header and Section 6 comments updated.

**`src/Cv_Numerical_Benchmark.py`**
- `run_cv_numerical_benchmark` and `print_cv_benchmark_summary` take an optional `classical_reference_label`. Its default is the quantum reference's label, which reproduces the old behavior.
- The classical figure's title uses that label, and the console names both reference grids.
- Both Section 6 figures reserve room for their two-line title (`tight_layout(rect=...)`). Before this, the title overlapped the legend and the top of the plot; this is visible in `figures/corrected_run/` here. Layout only; no number changes.
- Docstrings updated.

The quantum benchmark is unchanged: it still uses Section 2's 2 / 2 spectrum, which costs under a minute.

**Evidence.**
- **Same verification at a quarter of the cost:** the measurement in 7c, on the full 1000-temperature grid.
- **The wiring runs end to end.** `Quantum_HO_Master.py` was run unmodified with small settings (`NUM_STATES` 150, 30 temperatures, `BETA_MAX` 5, a three-variant Section 7; figures to a scratch folder):

  | Section 6 classical reference | outcome | Section 6 time | classical: max rel. difference to the base |
  |---|---|---|---|
  | default (1.5 / 1.5) | all seven sections, exit 0 | 315 s | 1.5e-13 |
  | `CLASSICAL_REFERENCE_PRECISE = True` (2 / 2) | stopped after Section 6 | 825 s | 2.0e-13 |

  - In both runs the console and the figure name the quantum reference (span×2, dx÷2) and the classical reference separately.
  - The quantum benchmark agreed to 3.7e-12 in both.
  - The title-layout fix came after these runs. It was checked by rendering both figures.

### 8b. Quick initial scan (`src/Quick_Scan.py`)

**What it is.** A first look at a candidate potential before a full run is committed to it. It draws only the quantum $C_v(T)$ (solid) and the classical limit (dashed), for one potential (`"single"`) or for the Section 7 sweep (`"sweep"`, each variant against its own classical limit), and prints a verdict per potential. Usage is in the README ("Quick Initial Scan"); the reasoning is in FINDINGS ("The Quick Scan").

**What it reuses (no new algorithm).**
- Classical limit: `Classical_Limit_Numerical.sweep_temperature_range` with a `ScaledSpectra` cache built at the chosen resolution (`thermal_coverage`, `dvr_tolerance`).
- Quantum $C_v$: `Quantum_Classical_Combined.compute_quantum_heat_capacity_curve`, on the ξ = 1 rung's spectrum (that of V itself, solved at the hottest temperature), so it costs no extra solve.
- Temperature grid: `Cv_AutoTune.resolve_beta_min`, exactly as the pipeline.
- Variants, formula and figure: `Cv_Coefficient_Sweep.generate_variant_params`, `format_potential_formula` and `plot_coefficient_sweep`. A failed variant is diagnosed with `_diagnose_variant_failure`, as in Section 7.
- Figure paths: `figures/output_paths.py`, category `quick_scan`.

**What is new in it.**
- `quick_settings` reads a preset and sizes the ladder so that it reaches `XI_MAX`.
- `compare_quantum_classical` classifies each temperature against the classical value's own error estimate ε: above (d > 2ε), below (d < −2ε) or unresolved. It also detects when the only unresolved band is the hot tail.
  - The margin `ERROR_MARGIN` = 2 was set after a first validation pass at margin 1. That pass found the true error up to 1.10ε, and found the coarse value always below the limit, which shifts d upward by about ε. With margin 1, a false "above" was therefore possible where the curves merge. That pass's run was stopped and repeated with the final code.
- `run_quick_scan` orchestrates; a small command-line interface (`--mode`, `--resolution`, `--beta-range`) overrides `config.py`.

**Other changed files.**
- **`src/config.py`:** `QUICK_SCAN_MODE`, `QUICK_SCAN_RESOLUTION`, `QUICK_SCAN_BETA_RANGE` and `QUICK_SCAN_PRESETS` (resolutions 1, 2 and 3; 3 is the full pipeline's settings).
- **`src/Cv_Coefficient_Sweep.py`:** `plot_coefficient_sweep` takes optional `title`, `category` and `name`, and returns the figure path. Section 7's call passes none of them, so its figure is unchanged and saved where it always was.

**Validation** ([`../quick_scan/validate_quick_scan.py`](../quick_scan/validate_quick_scan.py), [`../quick_scan/QUICK_SCAN_VALIDATION.txt`](../quick_scan/QUICK_SCAN_VALIDATION.txt)).

Compared with the full pipeline (the stored full-resolution Section 7 curves, [`../quick_scan/full_resolution_section7.npz`](../quick_scan/full_resolution_section7.npz), and each variant's 500-level quantum curve) and with the exact classical $C_v$ (`exact_classical_R` from this folder's audit script, used only as a check):

| resolution | run | time | quantum vs full pipeline | classical error vs exact | true error ÷ estimate | classical vs full-resolution ξ-scan | wrong resolved verdicts |
|---|---|---|---|---|---|---|---|
| 1 | six-variant sweep, 60 T | 5.0 min | ≤ 1.5e-7 | ≤ 5.2e-3 | 1.03–1.10 | ≤ 4.5e-3 | 0 |
| 2 | six-variant sweep, 150 T | 10.6 min | ≤ 1.8e-9 | ≤ 3.0e-3 | 1.01–1.04 | ≤ 2.4e-3 | 0 |
| 3 | b = −0.5, 1000 T | 9.9 min | 2.2e-12 | 8.4e-4 | 1.01 | 3.8e-12 | 0 |

- **Same conclusion as the full run.** Every variant is "below" at resolutions 1 and 2, the same as the ~44-minute full-resolution Section 7 run. The quick figures ([`../quick_scan/figures/`](../quick_scan/figures/)) show the same peaks.
- **Resolution 3 is the full pipeline.** Its classical curve is identical to the full run's to 4e-12.
- **Cheaper cold end.** The coarse tolerance reaches its plateau at a lower ξ: b = −0.9 needs ξ ≈ 985 at resolution 1, against 1,972 at full resolution.

**Open item: the WKB sizing overhead (not applied, needs your approval; resolved inside the quick scan in 8d).** Profiling resolution 1 on the base double well (17 solves, ~48 s) showed that about half of the time goes not to diagonalization but to `ScaledSpectra`'s WKB sizing. `_wkb_energy_for_count` bisects 60 times, and every step re-samples V on a 40,001-point window grown from scratch (`_allowed_region`). Sampling V once per cache, or bisecting to a relative 1e-6 instead of 60 halvings, would not change the method or its accuracy (the sizing only picks how many levels to solve for, with margins on top), and would make the quick scan roughly twice as fast (the full pipeline ~1 min faster). It is a change to the classical-limit engine, so it is left for your decision.

### 8c. Documentation

- **`README.md`:** the Section 6 description and runtimes (default 1.5 / 1.5, the precise option, full-run time); a new "Quick Initial Scan" section (what it is, what it reuses, the resolution presets with measured times, the verdict, usage); the repository tree, the `quick_scan` figure category, and a quick-scan step in "Running on a New Potential".
- **`FINDINGS.md`:** a physics subsection "The Quick Scan" (why a looser tolerance is cheap, why the hot tail stays unresolved); the quick-scan validation results; the new parameters; the Section 6 factor paragraph (now the default, with when to use the precise option); the quick-scan figure in the plot table; troubleshooting entries.
- **`HISTORY.md`:** the October 2026 chapter continues with the Section 6 change and the motivation, design, validation and profiling of the quick scan; four new version-table rows.
- **Not touched:** the `.tex` files.

### 8d. Follow-up: the sizing speed-up and a hot-end limit (quick scan only)

**Condition.** The sizing speed-up could be applied only if the full scan's classical-limit engine stayed untouched. `Classical_Limit_Numerical.py` is unchanged; everything below lives in the quick scan.

**`src/Quick_Scan.py`**
- **`RememberedPotential`.** It wraps V and returns a copy of the stored output whenever V is called again on an identical input array (the last 64 inputs). The quick scan hands it to the engine for every potential. The engine's WKB sizing re-samples V on the same 40,001-point windows dozens of times per solve, and those repeats are now free.
- **`_estimated_grid` and `hottest_feasible_T`.** The first estimates the DVR grid one ξ-scan solve needs at temperature T, from the engine's own sizing formula (1.15 × the WKB level count at ħ/ξ up to $E_0$ + 1.2 × coverage × $k_BT$, plus 10, times 4 points per level). The second bisects for the hottest temperature whose grid stays within `QUICK_SCAN_MAX_GRID`, checked at ξ = `xi_mult`^`min_stable_xi`, the earliest a plateau can end.
- **`scan_potential`.** Temperatures above that limit are skipped: both curves are NaN there, and a note gives the grid, level count, matrix size and hottest feasible T. If all temperatures are too hot, it raises an error with that explanation instead of attempting the solve.
- **`compare_quantum_classical`.** It takes `skipped`. Skipped temperatures are neither "missing" nor part of the hot-tail test.
- **Output.**
  - The summary has a "too hot" column.
  - No empty figure is saved when nothing was computed.
  - A `MemoryError` gets a one-line hint.
  - Differences print in scientific format: −0.0000 hid differences of 1e-5.
- **Docstring.** The module docstring has two new sections.

**`src/config.py`:** `QUICK_SCAN_MAX_GRID` = 6000. It is a soft limit, since the estimate can fall a few percent short.

**`audit/quick_scan/validate_quick_scan.py`:** pins the validated setup (b = −0.5, `BETA_MAX` 50, no zoom window, the b sweep) with `setattr` on `config`, so local experiments in `config.py` cannot change what it checks.

**Evidence.**
- **Bit-identical:** `cv_quantum`, `cv_classical`, `error_estimate` and the spectrum, base double well at resolution 1.
- **Speed:** 10.7 s instead of 26.0 s (Python 3.12 environment), 16.6 s instead of 31.6 s (Anaconda 3.7). Re-run of the validation: the six-variant sweep took 117 s instead of 299 s at resolution 1, and 277 s instead of 637 s at resolution 2. Resolution 3 on the base potential took 363 s instead of 595 s. The errors are unchanged, since the numbers are identical.
- **b = −10, T = 10⁴–10⁵ (the reported failure):** stops in under a second. At T = 10⁴ it would need ~83,000 grid points (~21,000 levels, a 56 GB matrix); the hottest feasible T is ≈ 1,170.
- **b = −10, T = 10–10⁴:** 41 of 60 temperatures computed and 19 skipped, in 264 s. The largest grid used was 6,257 points, hence "soft limit". Verdict: below at every computed temperature, the curves merging to within a few 1e-5 by T ≈ 1,000.

**Documentation.**
- **README:** the quick-scan section explains the speed-up and the "too hot" temperatures (with the b = −10 example), and the usage gains a "too hot" case.
- **FINDINGS:** the "Quick Scan" physics subsection covers the speed-up and why very hot temperatures are out of reach yet cost little. Results gain a b = −10 entry; parameters, `QUICK_SCAN_MAX_GRID`; troubleshooting, a "too hot" entry.
- **HISTORY:** the chapter continues with both changes; new version row `Quick_Scan` 1.1.
