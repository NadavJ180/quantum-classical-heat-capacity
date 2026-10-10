# Quantum-Classical Heat Capacity

A modular, system-agnostic numerical pipeline for the heat capacity $C_v(T)$ of a quantum particle in a smooth, confining 1-D potential, and for its classical limit.

- **Quantum $C_v(T)$:** computed from an energy spectrum obtained with the Discrete Variable Representation (DVR).
- **Classical limit:** found with the ξ-scaling of Gelbwaser-Klimovsky et al.: scale the potential and the temperature together, V → ξ²V and T → ξ²T, re-solve the Schrödinger equation for ξ²V at every ξ, and follow $C_v$ to its plateau. This is ħ → ħ/ξ at fixed V and T.
- **Self-validating:** every result is checked against an independently generated, higher-resolution numerical reference rather than a closed-form formula, so the same checks work for potentials with no analytic solution.
- **Quick initial scan:** [`src/Quick_Scan.py`](src/Quick_Scan.py) draws only the quantum and classical $C_v(T)$ of one potential, or of a coefficient sweep, with the same ξ method at a lower resolution. It takes about 20 seconds per potential at the lowest resolution, so a candidate can be checked before a full run (about an hour) is committed to it.

**Status:**
- The harmonic oscillator (HO) is validated against its exact solution.
- The classical limit is verified against the exact classical heat capacity (HO, $x^4$, and the double well) to better than $10^{-3}$.
- The pipeline currently runs on the quartic/cubic/quadratic double well of `src/config.py`.

Other documents:
- [`FINDINGS.md`](FINDINGS.md): the physics behind each step, how to read the diagnostic plots, and current results.
- [`HISTORY.md`](HISTORY.md): how the pipeline got here, including the classical-limit correction.
- [`CHANGES_SUMMARY.md`](CHANGES_SUMMARY.md): a short list of every change since the classical-limit scaling error, each with where it is in the code and the physical principle behind it.
- [`audit/classical_limit/`](audit/classical_limit/): the audit, the change log and the verification scripts for that correction. [`CHANGES.md`](audit/classical_limit/CHANGES.md) there is the file-by-file change log of the whole branch, including the Section 6 grid and the quick scan.
- [`audit/quick_scan/`](audit/quick_scan/): the quick scan's validation against the full pipeline and the exact classical $C_v$, and the numerical check of the high-temperature merge (Wigner–Kirkwood).
- [`docs/summaries/IEEE_Summary.tex`](docs/summaries/IEEE_Summary.tex) (the project report) and [`docs/summaries/Meetings_Summary.tex`](docs/summaries/Meetings_Summary.tex) (meeting notes and derivations). The report still describes the classical-limit method as it was before the correction.

---

## Repository Structure

```
src/
├── config.py                     Single source of truth: potential, constants, control parameters
├── Quantum_HO_Master.py          Master driver — entry point, runs the full 7-section pipeline
├── Quantum_Classical_Combined.py Cv pipeline for one potential: quantum Cv(T) + classical limit + plots
├── Classical_Limit_Numerical.py  Classical-limit engine: xi-ladder, DVR re-solve of xi^2 V per xi,
│                                 plateau detection with error estimate, n-convergence check
├── Cv_Numerical_Benchmark.py     Cv comparison: base grid vs. numerical reference (quantum + classical)
├── Cv_AutoTune.py                BETA_MIN auto-fill + escalation diagnostics (see "Auto-Tuning")
├── Cv_Coefficient_Sweep.py       Quantum and own classical Cv(T) across a swept coefficient (see "Coefficient Sweep")
├── Quick_Scan.py                 Quick initial scan: quantum + classical Cv(T) only, low resolution,
│                                 one potential or a sweep (see "Quick Initial Scan")
├── DVR/
│   ├── DVR_Algorithm.py          Core DVR solver and automatic grid configuration
│   ├── DVR_Reference_Generator.py  Numerical reference grid generator
│   └── DVR_Limit_Finder.py       DVR accuracy limit searches (resolution and level count)
├── analytical/
│   ├── HO_Analytical.py          Closed-form HO energy levels and Cv(T)
│   └── HO_Benchmark.py           External benchmark: numerical pipeline vs. analytic HO
├── error/
│   └── error_energylevels.py     Energy-level comparison (base vs. reference)
└── figures/
    ├── output_paths.py           Figure save-path resolver (see "Figure Output")
    ├── plot_potential.py         Potential-shape figures (called from Section 1 automatically)
    └── pipeline_diagram.py       Workflow-diagram generator

audit/classical_limit/   Audit of the pre-correction classical limit, the correction's change log,
                         the physics of the exact classical heat capacity, and verification scripts
audit/quick_scan/        Validation of the quick scan against the full pipeline and the exact classical Cv
figures/                 Generated plots (see "Figure Output"). HO/ and SymmetricDoubleWell/ are
                         hand-saved figures referenced directly by IEEE_Summary.tex; runs no longer
                         write into them.
docs/summaries/          IEEE_Summary.tex (report), Meetings_Summary.tex (meeting notes + derivations)
```

## Pipeline Overview

`Quantum_HO_Master.py` runs seven sections, all driven by the single potential defined in `config.py`. Section numbers keep their original roles; 1 and 4 form one auto-tuned loop.

1. **& 4. DVR base solve + Cv pipeline (auto-tuned).**
   - Solves for the lowest `NUM_STATES` levels on an automatically configured grid and computes the quantum $C_v(T)$ from them.
   - Computes the classical limit at every temperature with the ξ-scan (see "Classical Limit").
   - After each round, `Cv_AutoTune.diagnose_escalation` checks the results and grows `NUM_STATES` and/or the ξ ladder if needed, up to `MAX_ESCALATION_ROUNDS` (see "Auto-Tuning").
   - Saves the potential-shape figures once the loop settles.
2. **Numerical reference solve.** The same spectrum on an independently wider and finer grid (`REFERENCE_SPAN_FACTOR`, `REFERENCE_DX_FACTOR`), used as ground truth by the later sections.
3. **Energy-level accuracy.** Base vs. reference eigenvalues, absolute and relative error.
5. **DVR limit analysis.** The solver's own resolution (Δx) and level-count (n) breakdown points, measured against the reference.
6. **Cv numerical benchmark.**
   - Quantum $C_v$ from Section 2's reference spectrum.
   - Classical limit re-computed with every ξ²V solve on a grid widened and refined by `CLASSICAL_REFERENCE_SPAN_FACTOR` / `CLASSICAL_REFERENCE_DX_FACTOR` (default 1.5 / 1.5).
   - **More precise option:** `CLASSICAL_REFERENCE_PRECISE = True` uses Section 2's own factors (`REFERENCE_SPAN_FACTOR` / `REFERENCE_DX_FACTOR`, 2 / 2) instead. It is a stricter grid test, at about 3.7× the cost. On the double well both agree with the base curve to ~5e-12 (FINDINGS, "Key Parameters").
   - Both curves are compared with Section 4. Each figure states the grid of its own reference.
7. **Coefficient sweep.** Several variants of the base potential that differ in one coefficient. Each variant's quantum $C_v(T)$ is plotted against its own classical limit (see "Coefficient Sweep").

Only the potential block of `config.py` needs editing to run on a new potential.

**Runtime (double well, default settings):**

| part | time |
|---|---|
| Sections 1 & 4 | ~7 min, of which ~6.5 min is the classical limit (~33 DVR solves of ξ²V) |
| Sections 2, 3, 5 | ~2.5 min |
| Section 6 | ~10 min with the default 1.5 / 1.5 classical reference; ~36 min with `CLASSICAL_REFERENCE_PRECISE = True` (2 / 2) |
| Section 7 | ~6–8 min per non-base variant (one classical sweep each); ~37 min for the default six-variant sweep |
| **Full run** | **~1 h** (~1.4 h with the precise Section 6 reference) |
| `Quick_Scan.py` (separate) | resolution 1: ~20 s per potential (2 min for the six-variant sweep); resolution 2: ~45 s per potential (4.4 min); resolution 3 (the full settings): ~6–7 min per potential |

## Classical Limit

The classical limit is $C_v$ as ħ → 0 with the potential and the temperature held fixed. Following Gelbwaser-Klimovsky et al. (Eq. S7: $E_n(\hbar,\xi^2V)=\xi^2E_n(\hbar/\xi,V)$), it is reached by scaling the potential and the temperature together:

1. **The ladder:** ξ = `XI_START`·`XI_MULT`^k, never above `XI_MAX`.
2. **Each rung:** the DVR is solved for **ξ²V**, and $C_v$ is evaluated at ξ²T (`compute_cv(E_n(ξ²V), β, ξ)`). One solve per ξ serves every temperature, and each solve keeps the levels needed up to $E_0 + $ `HOT_STATE_SAFETY`·$k_BT$, so the result is never truncation-limited.
3. **Plateau test:** quantum corrections shrink as 1/ξ², so each step's distance from the ξ → ∞ limit is estimated as $|\Delta C_v|/(\text{XI\_MULT}^2-1)$.
   - A step is stable when that estimate is below `TOL_XI` and $C_v \ge \tfrac12 - $ `TOL_XI` (the classical $C_v$ can never be below ½).
   - `MIN_STABLE_XI` consecutive stable steps form the plateau. The last point is reported, together with its error estimate.
4. **Failure:** if the ladder runs out, hits `XI_MAX`, or a DVR solve fails, that temperature's classical limit is NaN and a warning is printed.
5. **Companion check:** an n-scan on the converged ξ's spectrum confirms the value doesn't depend on truncation.

The ξ-convergence diagnostic plot shows this scan at the hardest temperature: a frozen-out region, a rise, then the verified plateau. See [`FINDINGS.md`](FINDINGS.md) for the physics and [`audit/classical_limit/OPTION_A_physics.md`](audit/classical_limit/OPTION_A_physics.md) for the exact classical formula $C_v^{cl}/k_B = \tfrac12 + \beta^2\mathrm{Var}(V)$ that this limit converges to.

## Auto-Tuning

The temperature range (`BETA_MAX`, and optionally `BETA_MIN`) is meant to be the only knob you change from run to run.

- **`BETA_MIN`:** leave it as `None` and it is derived from $T_{\max}=10\,(E_1-E_0)/k_B$, a practical "T → ∞" point. Set a float to choose the hot end yourself.
- **`NUM_STATES`** controls only the quantum curve. It grows (`NUM_STATES_GROWTH`, capped by `NUM_STATES_CAP`) when the base spectrum's top level would be thermally accessible at the hottest temperature, $E_{\max} < $ `HOT_STATE_SAFETY`·$k_BT_{\text{hot}}$.
- **The ξ ladder** (`XI_START`, `MAX_XI_STEPS`) grows (`XI_START_GROWTH`, `MAX_XI_STEPS_GROWTH`) when more than `ESCALATION_FRACTION_THRESHOLD` of the temperatures ran out of ladder before reaching a plateau.
  - A scan stopped by the hard cap `XI_MAX`, or by a failed DVR solve, is **not** escalated. It is reported with a warning.
  - Raising `XI_MAX` is a deliberate cost decision, because grid sizes grow with ξ.
- Escalation is bounded by `MAX_ESCALATION_ROUNDS`. If the diagnostics still fail at the last round, a `UserWarning` is printed and the last attempt is used.
- The cache of ξ²V solves is reused across rounds, so a round that only grows `NUM_STATES` repeats no classical solve.
- `AUTO_ESCALATE = False` runs a single pass with the config values.

## Coefficient Sweep

Section 7 varies one coefficient of the current potential (e.g. the double well's cubic term `b`) and produces two figures in their own `coefficient_sweep/` folder:

- **`cv_coefficient_sweep.png`:** for every variant, its quantum $C_v(T)$ (solid) and **its own classical limit** (dashed), in the same color. Each quantum curve is therefore compared only with the classical curve of the same potential; the classical $C_v$ itself depends strongly on the coefficient.
  - Every classical curve is computed with the same ξ-scan and settings as the base run. The base variant reuses the base curve.
  - The legend lists each variant (with its $E_1-E_0$ gap) plus the two line styles. The title shows the potential's formula with the swept coefficient symbolic.
- **`potential_comparison.png`** (or one file per variant): each variant's $V(x)$ with its low-lying spectrum, zoomed on the well's structure and colored to match the Cv plot.
  - With up to 6 variants: side-by-side panels.
  - With up to 15: one overlay of the $V(x)$ curves.
  - Beyond that: one figure per variant (`potential_<param>_<value>.png`).

Config knobs:
- `SCAN_PARAM`: the `POTENTIAL_PARAMS` key to vary.
- `SCAN_STEP`: the spacing between variants.
- `SCAN_COUNT`: the number of extra variants on each side of the base value (the base is always included).
- `SCAN_SYMMETRIC_VALUE`: a value that recovers the symmetric potential, added as an extra reference variant if not already in the sweep; `None` disables it.
- `POTENTIAL_FORMULA`: the formula template for the plot title, with `<<name>>` placeholders, one per additive term (see `format_potential_formula`).

Behaviors to know:
- **Mirror pairs:** for a potential whose odd-degree terms are the only symmetry breakers (like `b x³`), +b and −b are mirror images with identical spectra, so their curves coincide exactly.
- **Failed variants:** a variant whose potential is not confining (its DVR solve fails) is skipped with a diagnosis instead of aborting the sweep. A variant whose classical scan fails at some temperatures is kept, with gaps in its dashed curve and a † in the legend.
- **Per-variant NUM_STATES:** each variant's `NUM_STATES` is checked against the same hot-end coverage criterion as the base run and escalated per variant if needed. A variant that still falls short is flagged with `*` (see `Cv_Coefficient_Sweep.py`'s docstring).

## Quick Initial Scan

`src/Quick_Scan.py` gives a first look at a candidate potential before a full run is committed to it. It produces **one figure**: the quantum $C_v(T)$ (solid) and the classical-limit $C_v(T)$ (dashed, same color), either for the potential in `config.py` (`"single"`) or for every variant of the Section 7 sweep (`"sweep"`), each against its own classical limit. It also prints a verdict per potential. There are no convergence figures, reference grids, auto-tuning or potential panels.

**The same method, at a lower resolution.** Nothing in it is a new algorithm. It calls the pipeline's own functions:
- the classical limit is the ξ-scan with ξ²V re-solved at every ξ (`sweep_temperature_range`);
- the quantum $C_v$ is `compute_quantum_heat_capacity_curve`;
- the variants and the figure come from Section 7 (`generate_variant_params`, `plot_coefficient_sweep`).

Only the settings are coarser. The quantum spectrum costs nothing extra: the ξ ladder starts at ξ = 1, and that rung, solved at the hottest temperature, is the spectrum of V itself.

The engine re-evaluates V on the same grids many times while sizing each solve. The quick scan hands it a potential that remembers its recent evaluations (`RememberedPotential`), which makes it about 2.5 times as fast. The results are bit-identical, and the classical-limit engine itself is untouched.

**Where quantum and classical merge ($T_{\text{merge}}$).** At high T the quantum $C_v$ approaches the classical one as a power law: the difference is $-(\beta\hbar\omega)^2/12$ for a harmonic well, and ∝ $T^{-(1+2/k)}$ for a $|x|^k$ tail (FINDINGS, "Why Quantum and Classical $C_v$ Merge at High T"). Computing there is unnecessary, and increasingly expensive: every temperature needs all its thermally accessible levels, and a dense DVR solve costs ~grid³.

So when the window reaches above $T_0 = 10\max(E_1-E_0,\ \hbar\omega_{\min})$, the quick scan first finds this potential's $T_{\text{merge}}$:
- it probes $C_v^q-C_v^{cl}$ at $T_0, 2T_0, 4T_0,\dots$ with the same ξ-scan;
- it accepts the first pair that agrees within `tol_xi` (with the 2ε margin) and is falling at least as fast as 1/T.

Temperatures above $T_{\text{merge}}$ are not computed. The console prints that quantum = classical within `tol_xi` there, the evidence (the two probes and the slope), and how many levels they would have needed. A window lying wholly above $T_{\text{merge}}$ gets the verdict **merged**.

For example, b = −10 (deep well ħω ≈ 30, second region at $k_BT\sim10^4$) has $T_{\text{merge}}\approx301$: its probes give −6.6e-4 at T = 301 and −1.6e-4 at T = 601, falling as $T^{-2.0}$. The window T = 10⁴–10⁵, which would need 9,000–56,000 levels per solve, is answered "merged" in 72 s. Default windows end at $10(E_1-E_0) \le T_0$ and are never affected.

**Resolutions** (`QUICK_SCAN_PRESETS` in `config.py`). Measured on the double well:

| resolution | temperatures | `tol_xi` | ξ ladder | plateau | levels kept up to | DVR tolerance | time per potential | six-variant sweep | classical error (measured) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 60 | 1e-2 | ×1.5 | 2 steps | $E_0+12\,k_BT$ | 1e-3 | ~20 s | 2.0 min | ≤ 5.2e-3 |
| 2 | 150 | 5e-3 | ×1.35 | 2 steps | $E_0+15\,k_BT$ | 1e-4 | ~45 s | 4.4 min | ≤ 3.0e-3 |
| 3 | 1000 | 2e-3 | ×1.25 | 3 steps | $E_0+20\,k_BT$ | 1e-5 | ~6–7 min | ≈ Sections 4 + 7 (~44 min) | ≤ 8.4e-4 |

At every resolution, the quantum curve matched the full pipeline's to ≤ 1.5e-7, and the true classical error was 1.01–1.10 times the scan's own estimate. The details are in [`audit/quick_scan/QUICK_SCAN_VALIDATION.txt`](audit/quick_scan/QUICK_SCAN_VALIDATION.txt).

Resolution 3 uses the full pipeline's own settings. A higher resolution costs more mainly at the cold end, which needs the largest ξ.

**The verdict.** Every classical value carries the ξ-scan's own estimate ε of its distance to the ξ → ∞ limit. At each temperature, d = $C_v^{q}-C_v^{cl}$ is classified:
- **above** if d > 2ε;
- **below** if d < −2ε;
- **unresolved** if |d| ≤ 2ε.

The factor 2 (`ERROR_MARGIN` in `Quick_Scan.py`) is a safety margin. The estimate tracks the true error to within ~10%, but a coarse classical value sits *below* the limit by about ε, which shifts d upward by the same amount. With a margin of 1, a temperature where the curves have merged could be called "above".

The console prints one line per potential and a summary table. At high T the quantum $C_v$ approaches the classical one, so a band of unresolved temperatures at the hot end is expected at any resolution, and the verdict says when that band is the only unresolved region.

**Usage:**
1. Set the potential in `config.py` (and `SCAN_PARAM`/`SCAN_STEP`/`SCAN_COUNT`/`SCAN_SYMMETRIC_VALUE` for a sweep).
2. Set `QUICK_SCAN_MODE` (`"single"` or `"sweep"`) and start at `QUICK_SCAN_RESOLUTION = 1`. Then, from `src/`, run:

   ```bash
   python Quick_Scan.py
   ```

   Command-line flags override `config.py`, for example `python Quick_Scan.py --mode single --resolution 2 --beta-range 0.5 5`.
3. Read the verdict:
   - **below** everywhere: no quantum excess at this resolution. The printed closest approach tells you by how much.
   - **ABOVE**: a candidate. Confirm it at a higher resolution, then with the full pipeline.
   - **unresolved** outside the hot tail, or gaps in a classical curve (†): raise the resolution, or zoom into that window with `QUICK_SCAN_BETA_RANGE` (`--beta-range`) and raise it there. A window that leaves out the cold end is much cheaper.
   - **merged**, or temperatures skipped above $T_{\text{merge}}$: the quantum and classical $C_v$ agree within `tol_xi` there (shown numerically, see the note), so there is nothing to look for. A finer resolution (smaller `tol_xi`) moves $T_{\text{merge}}$ up, at a cost that grows with T.
4. Run the full pipeline (`Quantum_HO_Master.py`) on a candidate that holds up.

The figure is saved to `figures/<system>/<params>/quick_scan/quick_<mode>_res<N>.png`, with `_beta<min>_to_<max>` appended for a zoom window, so a quick scan never overwrites the pipeline's figures.

## Figure Output

Every figure the pipeline saves is written to

```
figures/<system>/<params>/<category>/<name>.png
```

- **`<system>`:** a slug of `SYSTEM_NAME`.
- **`<params>`:** a slug of `POTENTIAL_PARAMS`, e.g. `a-0p25_b-m0p5_c-m0p5_d-0`. Runs with different parameters never overwrite each other.
- **`<category>`:** one of
  - `energy_levels` (including the two potential-shape figures),
  - `convergence` (ξ- and n-convergence diagnostics),
  - `cv` (the Cv summary and all Cv benchmarks),
  - `dvr_limits`,
  - `coefficient_sweep` (Section 7),
  - `quick_scan` (`Quick_Scan.py`).

All path components use only `[A-Za-z0-9_-]` (`-` → `m`, `.` → `p`), so the tree can be used directly from LaTeX. Paths are set by [`src/figures/output_paths.py`](src/figures/output_paths.py): a driver calls `set_context(SYSTEM_NAME, POTENTIAL_PARAMS)` once, and every plotting function calls `save_figure(fig, category, name)`, which saves and closes the figure (there is no `plt.show()` anywhere). The workflow diagram (`pipeline_diagram.py`) is saved to `figures/fig_pipeline.png`.

## Running on a New Potential

1. In `src/config.py`, set `POTENTIAL_PARAMS` and `my_potential(x)`. The potential must be finite everywhere (no hard walls) and confining. `my_potential` should read its coefficients from `POTENTIAL_PARAMS`, which also names the figure folder.
2. Update `SYSTEM_NAME`, `T_UNITS_LABEL` and `POTENTIAL_FORMULA`.
3. Set `BETA_MAX` (the cold end) for the new energy scale, and `NUM_STATES` as a starting guess. Leave `BETA_MIN = None` and the ξ settings at their defaults; the auto-tune loop adjusts them.
4. Optionally set `SCAN_PARAM`/`SCAN_STEP`/`SCAN_COUNT` for Section 7.
5. Check the candidate first with the quick scan (`python Quick_Scan.py` from `src/`; see "Quick Initial Scan").
6. Run `python src/Quantum_HO_Master.py`.

**Tested environments:** the full pipeline under Anaconda Python 3.7 (numpy 1.18, matplotlib 3.1); the quick scan under both that and Python 3.12 (numpy 2.4, matplotlib 3.10.9). The figures need a matplotlib that still provides `matplotlib.cm.get_cmap` (used by `Cv_Coefficient_Sweep.py`).

## Acknowledgements

**Supervisor:** Dr. David Gelbwaser-Klimovsky
**Student:** Nadav Jean
