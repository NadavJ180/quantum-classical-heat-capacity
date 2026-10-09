# History

How this pipeline got to its current state. For a meeting-by-meeting account (including the reasoning behind each change), see [`docs/summaries/Meetings_Summary.tex`](docs/summaries/Meetings_Summary.tex). This file is a condensed summary; full detail is in `git log`, and the classical-limit correction is documented in full in [`audit/classical_limit/`](audit/classical_limit/).

## Narrative

### Analytic systems, the DVR solver, and the original ξ-scan

The project began with the two textbook systems that have closed-form solutions, the **harmonic oscillator (HO)** and the **box potential**. They were used to validate the variance-based $C_v = k_B\beta^2\mathrm{Var}(E)$ formula and a "ξ-scaling" route to the classical limit.

That first ξ-scan was written for the box (commit `e9fa799`: `E_n = n**2 * E_g / xi**2`, the box spectrum with ħ → ħ/ξ). It was later generalized (commit `cb09603`) into a function that divides *any* computed spectrum by ξ² and scans ξ for a plateau. As described below, this generalization was only exact for the box, the HO and other single power laws. That went unnoticed until October 2026, because those were exactly the systems it was checked on.

A **Colbert–Miller sinc-DVR solver** was then built to obtain energy levels numerically, first tested against the HO. Early iterations tuned the grid by hand. This was replaced by an **automatic grid configurator** (turning-point root-finding plus a local Nyquist criterion).

Moving beyond the HO exposed several issues specific to anharmonic potentials:

- **Insufficient tail padding.** A fixed padding heuristic, calibrated for the HO, underestimated how far anharmonic wavefunction tails extend. It was replaced with iterative span widening that stops once eigenvalues stop changing.
- **Narrow ξ-plateau window.** The double well's levels scale as $E_n\propto n^{4/3}$, so the old scan's valid ξ-window was narrower than the HO's, and the ξ multiplier was reduced from 1.3 to 1.1. *This concern disappeared with the 2026 correction, whose scan has no finite-N window.*
- **Hard-wall potentials dropped.** The DVR core was restricted to smooth, everywhere-finite potentials.

The validation strategy also evolved. A **numerical reference generator** (an independently wider and finer DVR solve) replaced the HO's analytic formula as the ground truth for every run. The analytic HO comparison was kept as a one-time external certification in a separate module. The codebase was reorganized from a monolithic script into the modular `src/` layout.

### Auto-tuning, coefficient sweep, log-spaced temperatures

These additions came while exploring the quartic/cubic/quadratic double well:

- **Auto-tuning (`Cv_AutoTune.py`).** Sections 1 and 4 were merged into a closed loop that inspects the sweep's own diagnostics and grows `NUM_STATES` and the ξ-scan budget automatically. `BETA_MIN` gained an "auto" mode derived from $T_{\max}=10(E_1-E_0)/k_B$.
- **Coefficient sweep (`Cv_Coefficient_Sweep.py`, Section 7).** Overlays several variants' quantum $C_v(T)$ (one coefficient swept). Originally they were drawn against the **base** potential's classical-limit curve; this changed in 2026, see below.
  - Non-confining variants are caught and skipped instead of aborting the sweep. This also exposed and fixed a latent array-vs-scalar bug in `auto_configure_dvr`.
  - A second figure, `potential_comparison.png`, was added.
  - Reusing the base `NUM_STATES` for every variant was made safe by a per-variant hot-end coverage check (`solve_variant_with_hot_coverage`).
- **Log-spaced temperatures.** `beta_arr` switched from `np.linspace` to `np.geomspace`. With linear spacing, curves looked polygonal near the hot end on the log-T plots.

### October 2026: verifying the double-well "toy model" and correcting the classical limit

**The trigger.** The double well $V=\tfrac14x^4+bx^3-\tfrac12x^2$ looked like a candidate toy model. In the coefficient-sweep plot, its quantum $C_v$ showed a Schottky-like bump that rose above the "classical limit" and grew with $|b|$. Before investing in it, a verification plan was drawn up:
- code review;
- reproduction of Hasegawa's double-well results;
- analytic limits;
- quantitative Schottky signatures;
- a physical explanation.

**First check: the classical limit itself.** Each b's classical limit was computed with the pipeline and compared with the exact classical heat capacity of the same potential, $C_v^{cl}/k_B=\tfrac12+\beta^2\mathrm{Var}(V)$ (the one-dimensional phase-space integral; see [`OPTION_A_physics.md`](audit/classical_limit/OPTION_A_physics.md)).
- The pipeline's "classical limit" was a nearly flat ≈0.73 for every b.
- The exact classical $C_v$ of the base potential runs from 1.00 up to 1.45 and back down to 0.71, and depends strongly on b.

**The audit** (branch `correct-classical`, [`AUDIT.md`](audit/classical_limit/AUDIT.md)) found the cause:

1. **Which half was missing.** The paper the method comes from (Gelbwaser-Klimovsky et al., SI-III/IV, Eq. S7) scales **the potential and the temperature together**, V → ξ²V and T → ξ²T, which is ħ → ħ/ξ. The pipeline reused the spectrum of the *unscaled* V at every ξ, i.e. it applied only T → ξ²T, so each scan point was simply the quantum $C_v$ at the hotter temperature ξ²T. Its "plateau" was the system's high-temperature value, the same at every T. A search of every file and of the whole git history confirmed that the potential was never scaled anywhere.
2. **Why it went unnoticed.** Reusing the spectrum is exact precisely when every level gap scales by one common factor under V → ξ²V. That holds for the box (factor 1), the HO (factor ξ) and every |x|^k. It also holds only for those, which are exactly the systems whose classical $C_v$ does not depend on T.
3. **Why the HO checks passed.** For the HO the old scan was *identical* to the correct one with ξ read as ξ², to 7e-14. Moreover:
   - the HO's classical $C_v$ is the constant 1, so any method returning the high-T value passes;
   - the "machine-precision" classical-limit figure compared the base and reference grids running the *same* algorithm (a self-comparison);
   - the analytic HO benchmark only quantified the quantum curve.
4. **The decisive test.** The Pöschl–Teller well $V_0\tan^2x$ has an exact spectrum for every ξ²V and an analytic, T-dependent classical $C_v$. On it the old scan failed by 0.47, while the paper's prescription, evaluated with the pipeline's own `compute_cv`, converged onto the exact answer (4× per doubling of ξ, the ħ² law).
5. **Two smaller defects** in the same code:
   - the plateau picker could take its value from the collapsing finite-N tail;
   - when the ξ-scan failed, the quantum $C_v$ itself was stored as the "classical limit". In the HO validation run this was 90 of 500 points.

**The correction.** All code changes below are described file by file in [`CHANGES.md`](audit/classical_limit/CHANGES.md).

- **`Classical_Limit_Numerical.py`** was rewritten to follow the paper literally. At every ξ the DVR is **re-solved for ξ²V** and evaluated with the existing `compute_cv(E, β, ξ)`.
  - One solve per ξ serves every temperature. A cache (`ScaledSpectra`) sizes each solve to cover $E_0+$ `HOT_STATE_SAFETY`·$k_BT$, so finite-N collapse can no longer occur.
  - The plateau test uses the ħ² law to estimate each step's distance from the ξ → ∞ limit, and rejects "plateaus" below ½, the rigorous lower bound of the classical $C_v$. This stops the frozen-out regime at low T being mistaken for a plateau.
  - The reported value is the last point of the verified plateau, which fixes the picker. A failed scan gives NaN with a warning, which removes the fallback.
- **`DVR_Algorithm.auto_configure_dvr`** gained optional `energy_ceiling`, `turning_points` and `padding` arguments, with defaults unchanged. The scaled potentials need a different energy ceiling and a tighter starting window; the existing span-convergence loop and 3-pass check still validate every solve.
- **`config.py`:** new ξ settings:
  - `XI_START` = 1;
  - `XI_MULT` = 1.25;
  - `TOL_XI` = 2e-3, now a bound on the error of the classical value;
  - `MAX_XI_STEPS` = 35;
  - new hard cap `XI_MAX` = 2000.
- **Callers:** `Quantum_Classical_Combined.run`, the Section 6 benchmark (whose classical reference now re-solves every ξ²V on refined grids), `Cv_AutoTune` (failures now map to the knob that controls them) and the master script were adapted.
- **Docstrings:** `HO_Analytical`'s docstring, which described the classical limit as the T → ∞ value "by definition", was corrected.

**Verification of the correction:**
- **Against the exact classical $C_v$** (HO, x⁴ and the double well), the corrected scan converges at every temperature with errors ≤ 8.2e-4. Its own error estimate matches the true error to 6e-6.
- **The full pipeline** runs end to end. In Section 6, the base and refined-grid classical curves agree to 5e-12.
- **The audit**, now pinned to the pre-correction code, still reproduces its original results.

**Section 7 gets per-variant classical curves.** With a correct classical limit, comparing every variant's quantum $C_v$ with the *base* potential's classical curve became visibly misleading: b = −0.9 and −0.7 appeared to rise "above classical". Each variant is now drawn against its **own** classical limit, computed with the same scan and settings (quantum solid, classical dashed, same color).

**The answer for the toy model.** Run on all six variants (b = −0.9 … 0, 1000 temperatures, ~44 min; [`SECTION7_RESULTS.txt`](audit/classical_limit/SECTION7_RESULTS.txt)), every variant converged at every temperature.
- For **every** b, the quantum $C_v$ stays below its own classical $C_v$ at every temperature, by about 1e-3 at closest, at the hottest temperatures.
- The Schottky-like quantum bumps of b = −0.9 and −0.7 sit underneath taller classical peaks: those peaks come from a second region of configuration space opening up, and are not a quantum effect.

As a candidate for "quantum $C_v$ above classical", this potential family, in this range of b, does not show the effect.

**Section 6 cost and lighter grid factors (measured first, then applied; see below).** After the correction, Section 6's classical reference re-solves every ξ²V on a grid widened and refined by Section 2's factors (span×2, dx÷2). That took ~36 of the ~50 minutes of a full run. Lighter factors were measured on the full 1000-temperature grid ([`section6_grid_factor_check.py`](audit/classical_limit/section6_grid_factor_check.py), [`SECTION6_GRID_FACTORS.txt`](audit/classical_limit/SECTION6_GRID_FACTORS.txt)):

| span / dx factor | reference sweep | largest grid | max relative difference to the base curve |
|---|---|---|---|
| 2 / 2 (current) | 36.3 min | 12,785 pts | 5.5e-12 |
| 1.5 / 1.5 | 9.8 min | 7,192 pts | 5.4e-12 |
| 1.25 / 1.25 | 5.7 min | 4,995 pts | 6.9e-12 |

All three agree with the base at the round-off level. The classical curve is converged far beyond anything a factor can distinguish, so a lighter factor gives the same verification at a fraction of the cost.

### October 2026, continued: a lighter Section 6 reference and a quick initial scan

**Section 6's classical reference gets its own grid.** The measurement above was applied as a separate pair of factors, `CLASSICAL_REFERENCE_SPAN_FACTOR` / `CLASSICAL_REFERENCE_DX_FACTOR` = 1.5 / 1.5.
- **Cost.** Section 6 drops from ~36 to ~10 minutes, and a full run from ~1.4 h to ~1 h.
- **The precise option.** The stricter 2 / 2 reference is kept: `CLASSICAL_REFERENCE_PRECISE = True` gives the classical reference Section 2's own factors (or the ones entered at the prompt in interactive mode).
- **The quantum reference is unchanged.** Section 6's quantum benchmark still uses Section 2's 2 / 2 spectrum, which costs under a minute.
- **Labels.** Each Section 6 figure and console summary now names the grid of its own reference. Both figures' two-line titles had been overlapping the legend and the top of the plot; they now have room.

**Why a quick scan.** With a correct classical limit, the search for a potential whose quantum $C_v$ exceeds its own classical $C_v$ had a cost problem: each candidate took about an hour of the full pipeline. Most of that hour is verification (reference grids, Section 6, convergence figures), which matters only once a candidate is worth reporting. The requirements for a first-look tool were:
- only the quantum/classical $C_v(T)$ graph, for one potential or for a sweep of one coefficient;
- low resolution by default, so it runs quickly, with an option to raise the resolution when the result is not conclusive;
- no new algorithm: the existing pipeline functions and the fully numerical ξ method. The exact classical integral (option A) is not used, since that decision is still open.

**Design (`src/Quick_Scan.py`).** The quick scan is a thin driver over the pipeline's own functions:
- the classical limit is `sweep_temperature_range` (ξ²V re-solved at every ξ), run with a `ScaledSpectra` cache built at the chosen resolution;
- the quantum $C_v$ is `compute_quantum_heat_capacity_curve`;
- the variants and the figure are Section 7's `generate_variant_params` and `plot_coefficient_sweep`.

Only the settings change, through presets in `config.py` (`QUICK_SCAN_PRESETS`): fewer temperatures, a looser tolerance on the classical value, a coarser ξ ladder with a two-step plateau, fewer levels per solve and a looser DVR check. Resolution 3 is the full pipeline's own settings. The tolerance is the main lever, because the cold end needs the largest ξ and its cost grows as $\text{tol\_xi}^{-3/2}$ (FINDINGS, "The Quick Scan").

Two consequences of reusing the pipeline:
- **The quantum curve is free.** The ξ ladder starts at ξ = 1, and that rung, solved at the hottest temperature, is the spectrum of V itself.
- **Every classical value carries an error estimate.** The scan's own ε was shown, during the verification of the correction, to match the true error closely. The quick scan uses it to call each temperature "above" (d > 2ε), "below" (d < −2ε) or "unresolved", and says when the only unresolved band is the hot tail, where the two curves merge anyway.
  - The factor 2 came out of the validation. At resolution 1 the true error reached 1.10ε, and the coarse value always sits below the limit, which shifts d upward by about ε. A margin of 1 would therefore leave a false "above" possible exactly where the curves merge.

The only changes to existing code are default-preserving:
- `plot_coefficient_sweep` takes an optional title and save location, and returns the figure path, so the quick scan reuses Section 7's plot without overwriting its figure;
- `config.py` gained the `QUICK_SCAN_*` settings.

**Validation** ([`audit/quick_scan/`](audit/quick_scan/)).

The quick scan was run on the six-variant double-well sweep at resolutions 1 and 2, and on the base potential at resolution 3. It was compared with the full pipeline (the stored full-resolution Section 7 curves, and the pipeline's own 500-level quantum curves) and with the exact classical $C_v$:
- **Resolution 1** took 5 minutes for the sweep, against ~44 for Sections 4 + 7, and reached the same conclusion: every b below its own classical limit. Its classical curves were within 5.2e-3 of exact, and its quantum curves within 1.5e-7 of the full pipeline's.
- **Resolution 2** took 11 minutes, with classical errors ≤ 3.0e-3.
- **Resolution 3** reproduced the full pipeline's classical curve to 4e-12, as expected for the same computation.
- At every resolution the true classical error was 1.01–1.10 times the scan's own estimate, and no resolved verdict had the wrong sign.

A first validation pass also exposed two problems, both fixed before the pass was repeated:
- a margin of 1 in the verdict was too thin (see above);
- the validation script's interpolation assumed ascending temperatures.

**Profiling, and an optimization left for approval.** At resolution 1 about half of the run time is not diagonalization. It is the WKB sizing of each solve in `ScaledSpectra`: a bisection that re-samples the potential on 40,001-point grids dozens of times per solve. A cheaper sizing would not change the method or its accuracy (it only picks how many levels to solve for, with margins on top), and would make the quick scan roughly twice as fast (and the full pipeline ~1 min faster). It changes the classical-limit engine, however, so it is listed as an open item instead of being applied.

## Version History

Pre-reorganization filenames carried explicit version suffixes (e.g. `DVR_Algorithm_1_4.py`); current files no longer do. This table records the module-level changes that shaped the current design.

| File | Version | Key change |
|------|---------|-----------|
| `DVR_Algorithm` | 1.3 | Hard-wall support removed; smooth-only |
| `DVR_Algorithm` | 1.4 | Multiprocessing timer removed; replaced with inline per-pass timing |
| `DVR_Algorithm` | 1.5 | Adaptive span-expansion loop added to `auto_configure_dvr`; `E_ceiling` reverted from an oversized temporary hack back to `1.5 * num_levels` |
| `DVR_Algorithm` | 1.6 | Optional `energy_ceiling`, `turning_points`, `padding` in `auto_configure_dvr` (defaults unchanged), for the scaled potentials ξ²V |
| `Classical_Limit_Numerical` | 1.0 | Extracted from the combined file; fully general |
| `Classical_Limit_Numerical` | 2.0 | Classical-limit correction: ξ²V re-solved at every ξ (Eq. S7), ħ²-law error estimate, ½ lower bound, plateau picker fixed, no fallback |
| `Quantum_Classical_Combined` | 1.9 | System-agnostic Cv pipeline; xi/n engine extracted |
| `Quantum_Classical_Combined` | 2.0 | `run()` takes the potential; quantum Cv from the whole base spectrum; new ξ-diagnostic plot |
| `DVR_Reference_Generator` | 1.0 | New: numerical reference grid generation, replacing the analytic formula as ground truth |
| `DVR_Limit_Finder` | 1.0 | New: resolution and level-count limit searches |
| `DVR_Limit_Finder` | 1.1 | $\Delta x$ replaces point count as the resolution-search axis |
| `DVR_Limit_Finder` | 1.2 | Point annotations removed from the $\Delta x$ plot |
| `Cv_Numerical_Benchmark` | 1.0 | New: base vs. reference Cv comparison (quantum + classical) |
| `Cv_Numerical_Benchmark` | 1.1 | Classical reference re-solves every ξ²V on refined grids; data-scaled classical plot |
| `Cv_Numerical_Benchmark` | 1.2 | Separate label for the classical reference grid (`classical_reference_label`), shown on its figure and in the console; the two-line titles no longer overlap the top panels |
| `HO_Energy_Level_Error` | 1.1 | Docstring clarified — function is fully generic |
| `Quantum_HO_Master` | 1.5 | Analytical sections removed; fully numerical pipeline (numerical reference is now the sole ground truth for Sections 3, 5, 6) |
| `Quantum_HO_Master` | 1.6 | Section 6's classical reference on its own grid factors (`CLASSICAL_REFERENCE_*`, default 1.5/1.5), with the precise option `CLASSICAL_REFERENCE_PRECISE` (Section 2's 2/2) |
| `Cv_AutoTune` | 1.0 | New: `BETA_MIN` auto-fill from $T_{\max}=10\Delta E/k_B$, and escalation diagnostics for the `NUM_STATES`/`XI_START` closed loop |
| `Cv_AutoTune` | 1.1 | Escalation mapped to the corrected scan: `NUM_STATES` for the quantum curve only, ξ ladder on `max_steps`; `xi_cap`/`dvr_failed` reported, not escalated |
| `Cv_Coefficient_Sweep` | 1.0 | New: Section 7 coefficient-sweep comparison plot (quantum $C_v(T)$ per variant vs. the base run's reused classical limit) |
| `Cv_Coefficient_Sweep` | 1.1 | Each variant drawn against its own classical limit (quantum solid, classical dashed, same color) |
| `Cv_Coefficient_Sweep` | 1.2 | `plot_coefficient_sweep` takes an optional title and save location (`title`, `category`, `name`) and returns the figure path, so the quick scan can reuse it |
| `Quick_Scan` | 1.0 | New: quick, low-resolution quantum vs. classical $C_v(T)$ for one potential or a coefficient sweep, with a resolved/unresolved verdict; reuses the pipeline's functions |
