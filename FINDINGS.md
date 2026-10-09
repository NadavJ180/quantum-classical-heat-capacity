# Findings

Technical summary of the physics behind this pipeline, how to read its diagnostic output, and the results so far. For the full write-up see [`docs/summaries/IEEE_Summary.tex`](docs/summaries/IEEE_Summary.tex); note that its classical-limit section still describes the pre-correction method (see [`HISTORY.md`](HISTORY.md)).

## Physics Background

### Heat Capacity from the Partition Function

Given energy eigenvalues $\{E_n\}$ at inverse temperature $\beta=1/(k_BT)$, the canonical-ensemble heat capacity is

$$C_v(T) = k_B \beta^2 \left[\langle E^2 \rangle - \langle E \rangle^2\right] = k_B\beta^2\,\mathrm{Var}(E), \qquad \langle E^k \rangle = \frac{\sum_n E_n^k\, e^{-\beta E_n}}{Z},\ \ Z = \sum_n e^{-\beta E_n}.$$

This holds regardless of whether $\{E_n\}$ comes from a closed-form expression or a numerical diagonalization. Throughout, the system is assumed to be in thermal equilibrium with a bath at temperature $T$ (the canonical ensemble). This is the standard basis for a quasi-static engine-cycle picture, not a finite-time or non-equilibrium treatment.

### The DVR Method

The Discrete Variable Representation (DVR) solves the 1-D time-independent Schrödinger equation $\hat H\psi=E\psi$, $\hat H=\hat T+V(x)$, by representing the wavefunction by its values on an evenly spaced grid. This implementation uses the **Colbert–Miller sinc-DVR** (1992). Its kinetic-energy matrix on a grid of spacing $\Delta x$ is exact for band-limited functions sampled on that grid:

$$T_{ij} = \frac{\hbar^2}{2m\,\Delta x^2}\times\begin{cases}\dfrac{\pi^2}{3}, & i=j,\\[4pt]\dfrac{2\,(-1)^{i-j}}{(i-j)^2}, & i\neq j.\end{cases}$$

The potential energy is diagonal, $V_{ij}=V(x_i)\,\delta_{ij}$. $T_{ij}$ depends only on $|i-j|$ (a Toeplitz matrix), so it is built from a single row. The Hamiltonian is diagonalized with `scipy.linalg.eigvalsh` (`subset_by_index`), which extracts only the lowest $N_{\text{lev}}$ eigenvalues.

- **Scope:** $V(x)$ must be finite everywhere on the grid. Hard walls are not supported: a discontinuity has unbounded momentum content that no finite grid can resolve.
- **Grid limitations:** only roughly the lower half of a requested spectrum is trustworthy, and the dense diagonalization costs $\mathcal{O}(N^2)$ memory and $\mathcal{O}(N^3)$ time.

### The Classical Limit

**What it is.** The classical limit is $C_v$ as ħ → 0 with the potential, the mass and the temperature held fixed. For $H=p^2/2m+V(x)$ it equals

$$\frac{C_v^{cl}}{k_B} = \frac12 + \beta^2\Big(\langle V^2\rangle-\langle V\rangle^2\Big),\qquad \langle V^k\rangle=\frac{\int V^k e^{-\beta V}dx}{\int e^{-\beta V}dx}.$$

The ½ is the kinetic energy's equipartition share, and the second term is the fluctuation of the potential energy. The derivation and its meaning are in [`audit/classical_limit/OPTION_A_physics.md`](audit/classical_limit/OPTION_A_physics.md).

Two consequences:
- $C_v^{cl}\ge\tfrac12k_B$ always.
- $C_v^{cl}$ is independent of T only for single power laws:
  - HO: 1;
  - $|x|^k$: $\tfrac12+\tfrac1k$, e.g. ¾ for $x^4$;
  - box: ½.

  For any other potential it depends on T. The double well's classical $C_v$ runs 1 → 1.45 → 0.71 (see results below).

**How the pipeline reaches it: ξ-scaling.** Following Gelbwaser-Klimovsky et al. (Eq. S7),

$$E_n(\hbar,\ \xi^2V) = \xi^2\,E_n\!\left(\tfrac{\hbar}{\xi},\ V\right),$$

so scaling the potential **and** the temperature by ξ² is ħ → ħ/ξ at fixed V and T, and ξ → ∞ is the classical limit.

- **Each ladder step.** At every ξ on a geometric ladder, the DVR is re-solved for **ξ²V** and $C_v$ is evaluated at ξ²T (`compute_cv(E_n(ξ²V), β, ξ)`, whose weights are $e^{-\beta E/\xi^2}$).
  - The spectrum of ξ²V does not depend on T, so one solve per ξ serves every temperature (cached, processed hot → cold).
  - Each solve keeps the levels up to $E_0+$ `HOT_STATE_SAFETY`·$k_BT$, so the result is never truncation-limited.
- **Why ξ²V must be re-solved.** Re-using the spectrum of V at every ξ (the method before October 2026) applies only T → ξ²T, and converges to the high-temperature quantum $C_v$ at every T. That is the correct classical limit only when every level gap scales by one common factor under V → ξ²V, i.e. for single power laws. See [`HISTORY.md`](HISTORY.md) and [`audit/classical_limit/AUDIT.md`](audit/classical_limit/AUDIT.md).

**Convergence and its error estimate.** The leading quantum correction is of order ħ². The Wigner–Kirkwood expansion gives $F_q=F_{cl}+\frac{\hbar^2\beta}{24m}\langle V''\rangle_{cl}+O(\hbar^4)$; for the HO, $C_v/k_B=1-(\beta\hbar\omega)^2/12+\dots$. With ħ → ħ/ξ, the remaining distance to the limit falls as 1/ξ², so on a ladder with ratio $m$ (`XI_MULT`):

$$\varepsilon_k=\frac{|C_v(\xi_k)-C_v(\xi_{k-1})|}{m^2-1}$$

estimates how far step $k$ still is from the limit. The rules:
- A step is **stable** when $\varepsilon_k<$ `TOL_XI` **and** $C_v\ge\tfrac12-$ `TOL_XI`. The second condition rejects the frozen-out regime at low T, where $C_v\approx0$ for several steps.
- `MIN_STABLE_XI` consecutive stable steps form the plateau.
- The reported value is the plateau's last point, together with its $\varepsilon$.
- If the ladder runs out, hits `XI_MAX`, or a DVR solve fails, the value is NaN; it is never replaced by another quantity.
- A companion **n-scan** on the converged ξ's spectrum checks that the value does not depend on how many levels were kept.

In practice the reported value approaches the limit from below by almost exactly its own $\varepsilon$.

### Coefficient Sweep: Per-Variant Classical Limits and NUM_STATES Reuse

**Per-variant classical limits.** The classical $C_v$ depends on the potential, and strongly so for this family. Section 7 therefore compares each variant's quantum $C_v$ with **its own** classical limit, computed with the same scan and settings as the base run. Comparing with another potential's classical curve would say nothing about quantum effects.

**Reusing `NUM_STATES`.** Each variant's quantum $C_v$ is evaluated directly from its own spectrum of `NUM_STATES` levels. This is exact wherever the omitted levels carry negligible Boltzmann weight, $E_{\max}-E_0\gg k_BT_{\text{hot}}$ (`HOT_STATE_SAFETY`).

- **Why reuse can fail.** `NUM_STATES` is tuned for the base potential only. WKB quantization of a quartic-dominated well gives

  $$\oint p\,dx = \left(n+\tfrac12\right)\pi\hbar,\qquad p=\sqrt{2m(E-ax^4)}\;\Rightarrow\; E_n\sim a^{1/3}n^{4/3}.$$

  Sweeping a leading coefficient downward therefore shrinks $E_{\max}$ for the same level count.
- **The safeguard.** `solve_variant_with_hot_coverage` checks the criterion per variant, escalates that variant's `NUM_STATES` if needed, and flags (`*`) a variant that still falls short.

### The Quick Scan: Same Method, Coarser Settings

`Quick_Scan.py` runs exactly the ξ-scan above (and the same quantum $C_v$ formula) at coarser settings. Its speed comes from where the cost of the ξ-scan sits.

- **The cold end dominates.** The remaining distance to the limit falls as $a(T)/\xi^2$, so a tolerance `tol_xi` is met at $\xi\propto\sqrt{a(T)/\text{tol\_xi}}$. The ħ² coefficient $a(T)$ grows at low T (for the HO, $C_v/k_B=1-(\beta\hbar\omega)^2/12+\dots$), which is why $\xi_{\text{conv}}$ grows roughly like 1/T.
- **Cost per rung.** ξ²V has the level density of V at ħ/ξ, so the number of levels needed to cover $k_BT$, and with it the DVR grid, grows ∝ ξ. A dense diagonalization costs ∝ (grid)³.
- **Result.** The coldest rungs cost ∝ $\xi^3\propto\text{tol\_xi}^{-3/2}$. Loosening `tol_xi` from 2e-3 to 1e-2 lowers the largest ξ by about √5 and the coldest solves by roughly an order of magnitude. On the base double well, the largest ξ went from 1,262 to 657 and the largest grid from 3,197 to 989 points.
- **The other knobs:**
  - fewer temperatures only coarsen the curve, since the solves are per ξ and shared by all temperatures;
  - a coarser ladder (×1.5) and a two-step plateau need fewer rungs;
  - fewer levels per solve: a coverage of 12 $k_BT$ leaves the top level a Boltzmann weight of e⁻¹²;
  - the looser DVR tolerance only relaxes the pass/fail threshold of the 3-pass check.
- **The quantum curve is free.** The ladder starts at ξ = 1, and that rung, solved at the hottest temperature, is the spectrum of V itself, from the same DVR and 3-pass check.

**Reading "quantum above classical".** Each classical value carries the scan's estimate ε of its remaining distance to the limit, and the true error tracks it to within ~10% (results below). The difference d = $C_v^q-C_v^{cl}$ counts as resolved only where |d| > 2ε. The margin is needed because the approach is one-sided: a coarse classical value sits *below* the limit by about ε, so the measured d is shifted upward by about ε. Wherever the true d is close to zero, a margin of 1 could therefore turn it into a false "above". With 2, that would need the error to exceed twice its estimate.

At high T the two curves merge: the quantum correction falls off as ħ²β² (Wigner–Kirkwood), so d → 0. At the hottest temperatures |d| eventually drops below any finite ε, and an unresolved band at the hot end is expected at any resolution. Unresolved temperatures elsewhere mean the resolution is too low for that window.

## Findings So Far

- **HO.**
  - The base-grid energy levels agree with the numerical reference to machine precision.
  - The quantum $C_v(T)$ agrees with the exact Einstein formula to ~3e-13.
  - The classical limit agrees with the exact value 1 to within 8e-4 at every temperature.
- **DVR limits.** Characterized directly for the HO, where only roughly the lower half of a requested spectrum should be trusted. Section 5 measures the same limits for every run; for the base double well it found Δx ≤ 0.037 for 500 levels, and 582 trustworthy levels at the run's own spacing.
- **Classical-limit verification.** Against the exact classical $C_v$, the scan converges at every temperature for the HO, $x^4$ and the double well. The largest error is 8.2e-4 (below `TOL_XI` = 2e-3), and the error estimate matches the true error to 6e-6 ([`audit/classical_limit/VERIFY_CORRECTED.txt`](audit/classical_limit/VERIFY_CORRECTED.txt)).
- **Double well, base b = −0.5** ($V=\tfrac14x^4-\tfrac12x^3-\tfrac12x^2$):
  - The classical $C_v$ is ≈1.00 at low T (harmonic in the deep well), rises to 1.45 near T ≈ 0.55 as a second region of configuration space becomes accessible, and falls toward ¾ at high T (0.71 at the hottest grid T).
  - The quantum $C_v$ stays **below** the classical $C_v$ at every temperature, by −0.0014 at most.
  - The classical curve is converged in the DVR grid: base and reference grids agree to 5e-12.

- **Double well across b (Section 7, each variant against its own classical limit)** ([`audit/classical_limit/SECTION7_RESULTS.txt`](audit/classical_limit/SECTION7_RESULTS.txt), figure `figures/corrected_run/cv_coefficient_sweep_own_classical.png` there):

  | b | max(quantum − own classical) | quantum $C_v$ peak | own classical $C_v$ peak |
  |---|---|---|---|
  | −0.9 | −0.0021 | 1.28 at T ≈ 2.6 | 1.42 at T ≈ 2.3 |
  | −0.7 | −0.0013 | 1.13 at T ≈ 1.5 | 1.43 at T ≈ 1.2 |
  | −0.5 | −0.0010 | 0.86 at T ≈ 0.91 | 1.45 at T ≈ 0.54 |
  | −0.3 | −0.0010 | none (rises monotonically to 0.72) | 1.50 at T ≈ 0.23 |
  | −0.1 | −0.0009 | none | 1.64 at T ≈ 0.077 |
  | 0 (symmetric) | −0.0009 | none | 1.14 at T ≈ 0.051 |

  - **No quantum excess.** For every b the quantum $C_v$ lies **below its own classical $C_v$ at every temperature**; the two meet only at the hottest temperatures. The "quantum above classical" bumps of the earlier coefficient-sweep plot came entirely from comparing against the wrong classical curve.
  - **The bumps are classical.** The quantum peaks of b = −0.9 and −0.7 are softened versions of a peak the classical $C_v$ already has. It comes from configuration space: as T rises, a second region (the shoulder or shallow well) becomes accessible, and the potential energy fluctuates strongly.
    - The classical peak moves to lower T and grows as the two wells approach degeneracy (b → −0.1).
    - In the symmetric case (b = 0) both wells are equivalent, so that two-region contribution largely disappears.
  - **Ladder margin.** The stiffer the deep well, the larger the ξ the coldest temperatures need. b = −0.9 reached ξ ≈ 1972, the last rung of the default ladder, with all temperatures converged but no spare margin. More negative b, or a colder `BETA_MAX`, will need a longer ladder (`MAX_XI_STEPS`, `XI_MAX`). Section 7 variants do not auto-escalate.

- **Quick scan, validated on the same sweep** ([`audit/quick_scan/QUICK_SCAN_VALIDATION.txt`](audit/quick_scan/QUICK_SCAN_VALIDATION.txt)). It was compared with the full pipeline (the stored full-resolution Section 7 curves and the 500-level quantum curves) and with the exact classical $C_v$:

  | resolution | run | time | quantum vs full pipeline | classical error vs exact | true error ÷ estimate | wrong resolved verdicts |
  |---|---|---|---|---|---|---|
  | 1 | six-variant sweep | 5.0 min | ≤ 1.5e-7 | ≤ 5.2e-3 | 1.03–1.10 | 0 |
  | 2 | six-variant sweep | 10.6 min | ≤ 1.8e-9 | ≤ 3.0e-3 | 1.01–1.04 | 0 |
  | 3 | b = −0.5 only | 9.9 min | 2.2e-12 | 8.4e-4 | 1.01 | 0 |

  - **Resolution 1 reaches the full run's conclusion in 5 minutes instead of ~44:** every b below its own classical limit at every temperature. The closest approach is at the hottest temperature, by about 1e-3, still resolved there (|d| ≈ 4ε).
  - **The error estimate is honest at low resolution too.** The coarse classical value always sits below the limit, by about its own estimate.
  - **Resolution 3 is the full pipeline.** Its classical curve is identical to the full run's to 4e-12.

## Reading the Diagnostic Plots

| Plot | What to look for |
|---|---|
| **Energy-level comparison** (base vs. reference, full range + zoom near largest error) | Curves visually indistinguishable at full scale; the zoom panel shows the worst-case disagreement in context. |
| **Relative error vs. state index $n$** | A smooth, gently rising curve. Low-lying states are most accurate, the highest least. An abrupt spike at some $n^*$ flags where grid resolution first becomes insufficient. |
| **ξ-convergence diagnostic** (the hardest temperature) | $C_v(\xi)$ on a log-ξ axis. At low T: ≈0 for small ξ (frozen out; grey, not accepted), then a rise, then the **verified plateau** (green), whose last point is reported with its estimated distance to the limit. Yellow points passed the per-step test but did not complete a plateau. There is no collapse: each ξ's spectrum covers its temperature. |
| **n-convergence diagnostic** | Run on the converged ξ's spectrum. Near-zero at low $n$, a smooth rise, a flat tail. The converged point should sit well before the right edge. |
| **$C_v(T)$ summary** (quantum + classical limit + $\xi_{\text{conv}}(T)$/$n_{\text{conv}}(T)$) | The classical curve is the ħ → 0 heat capacity at each T and can exceed $k_B$. The quantum curve rises from 0 and joins it at high T. $\xi_{\text{conv}}$ grows roughly like 1/T: colder temperatures need a smaller effective ħ. |
| **DVR resolution/level-count limit plots** | A long machine-precision floor, then an abrupt cliff once Δx (or the requested n) crosses the solver's breakdown point. |
| **Cv benchmark plots** (base vs. reference) | Flat, featureless error well below the tolerances. Quantum: the reference spectrum. Classical: the same ξ-scan with every ξ²V solve on refined grids. This shows grid convergence; it is not an independent test of the method. |
| **Coefficient sweep — Cv** (`coefficient_sweep/cv_coefficient_sweep.png`) | One color per variant: quantum $C_v$ solid, **its own** classical limit dashed. A quantum excess means the solid line above the dashed line *of the same color*. A † marks a variant whose classical scan failed at some temperatures (gaps); a `*` marks marginal hot-end coverage of the quantum curve. Mirror-image coefficient values (±b) give identical curves. A non-confining variant is simply absent (see console). |
| **Coefficient sweep — potentials** (`potential_comparison.png` or `potential_<param>_<value>.png`) | Each variant's $V(x)$ with its low-lying spectrum, colored to match the Cv plot. Use it to connect changes in the Cv curves to changes in well depth, barrier and asymmetry. |
| **Quick scan** (`quick_scan/quick_<mode>_res<N>.png`) | Same layout as the coefficient-sweep Cv plot, at the quick scan's resolution (in the title). Read it together with the console verdict: a solid line above its dashed line counts only where the verdict calls it resolved. |

## Key Parameters (`src/config.py`)

| Parameter | Effect |
|-----------|--------|
| `NUM_STATES` | Levels of the base spectrum, used for the **quantum** curve only. A starting guess; grown automatically if the top level is thermally accessible at the hottest T. |
| `BETA_MIN` / `BETA_MAX` / `N_BETA` | Temperature window and point count. `BETA_MIN = None` derives the hot end from $T_{\max}=10(E_1-E_0)/k_B$. `beta_arr` is log-spaced; raise `N_BETA` if a curve looks polygonal. Colder `BETA_MAX` means a longer ξ ladder (ξ grows roughly like 1/T). |
| `XI_START`, `XI_MULT`, `MAX_XI_STEPS` | The ξ ladder `XI_START`·`XI_MULT`^k, k < `MAX_XI_STEPS`. Every rung is a DVR solve of ξ²V. The defaults are 1, 1.25 and 35, reaching ξ ≈ 1,970. |
| `XI_MAX` | Hard cap on ξ (default 2000). It is never raised automatically, because grid sizes grow with ξ. |
| `TOL_XI`, `MIN_STABLE_XI` | Plateau criterion. `TOL_XI` bounds the estimated distance of the reported value from the ξ → ∞ limit (default 2e-3; observed errors ≈ 8e-4). `MIN_STABLE_XI` = 3 consecutive stable steps. |
| `TOL_CV`, `MIN_STABLE_N` | n-scan stability tolerance and run length. |
| `HOT_STATE_SAFETY` | Thermal coverage: base spectrum $E_{\max}\ge$ safety·$k_BT_{\text{hot}}$, and every ξ²V solve keeps levels up to $E_0+$ safety·$k_BT$. Default 20, so the top level weighs ~e⁻²⁰. |
| `LIMIT_TOLERANCE` | Pass/fail threshold of the DVR limit searches (Section 5). |
| `REFERENCE_SPAN_FACTOR` / `REFERENCE_DX_FACTOR` | How much wider/finer the reference grid is (Section 2), used by Sections 3 and 5 and by Section 6's quantum benchmark. Default 2.0/2.0. |
| `CLASSICAL_REFERENCE_SPAN_FACTOR` / `CLASSICAL_REFERENCE_DX_FACTOR` | Grid widening/refinement of every ξ²V solve in Section 6's classical reference. Default 1.5/1.5. |
| `CLASSICAL_REFERENCE_PRECISE` | `True` gives Section 6's classical reference Section 2's own factors (2/2): a stricter grid test at ~3.7× the cost. Default `False`. |
| `AUTO_ESCALATE`, `MAX_ESCALATION_ROUNDS`, `NUM_STATES_GROWTH`, `NUM_STATES_CAP`, `XI_START_GROWTH`, `MAX_XI_STEPS_GROWTH`, `ESCALATION_FRACTION_THRESHOLD` | Control the auto-tune loop (see README). Meant to be touched rarely. |
| `SCAN_PARAM`, `SCAN_STEP`, `SCAN_COUNT`, `SCAN_SYMMETRIC_VALUE` | Section 7: which coefficient is swept, its spacing, the extra variants per side, and an optional symmetric reference value. The quick scan's `"sweep"` mode uses the same four. |
| `QUICK_SCAN_MODE`, `QUICK_SCAN_RESOLUTION`, `QUICK_SCAN_BETA_RANGE` | `Quick_Scan.py`: one potential (`"single"`) or the Section 7 sweep (`"sweep"`); which preset; an optional (β_min, β_max) zoom window (`None` = `BETA_MIN`/`BETA_MAX`). |
| `QUICK_SCAN_PRESETS` | Per resolution: `n_beta`, `tol_xi`, `xi_mult`, `min_stable_xi`, `thermal_coverage` and `dvr_tolerance`, with the same meanings as the pipeline's `N_BETA`, `TOL_XI`, `XI_MULT`, `MIN_STABLE_XI` and `HOT_STATE_SAFETY`, and the DVR's 3-pass tolerance. The ladder always runs from `XI_START` up to `XI_MAX`. Resolution 3 is the full pipeline's own settings. |

**Section 6 reference factors.** With 2/2, Section 6's classical reference sweep took ~36 min for the double well. Lighter factors agree with the base curve equally well, all at round-off. The default is therefore 1.5/1.5, with 2/2 kept as the `CLASSICAL_REFERENCE_PRECISE` option:

| factor (span / dx) | time | max relative difference |
|---|---|---|
| 2 / 2 | 36 min | 5.5e-12 |
| 1.5 / 1.5 | 9.8 min | 5.4e-12 |
| 1.25 / 1.25 | 5.7 min | 6.9e-12 |

See [`audit/classical_limit/SECTION6_GRID_FACTORS.txt`](audit/classical_limit/SECTION6_GRID_FACTORS.txt) and `figures/section6_grid_factors.png` there.

- **Why 1.5/1.5 is still a real test.** It changes both the span and the spacing of every ξ²V solve by 50%. DVR errors fall off exponentially with resolution and span, so a base solve that is too narrow or too coarse would still show up as a visible difference.
- **When to use the precise option.** When the Section 6 classical error panel is not at round-off level, or when a new potential's grid is in doubt: 2/2 perturbs the grid more strongly.
- **The quantum benchmark is unaffected.** It reuses Section 2's 2/2 spectrum, which costs under a minute.

## Troubleshooting

The auto-tune loop handles the two common cases automatically (README, "Auto-Tuning"). The manual remedies below still apply if escalation is disabled, if it runs out of rounds (a `UserWarning` is printed), or for Section 6 and Section 7, which reuse the base run's final settings.

- **Classical limit is NaN at some temperatures.** The console reports the stop reason per temperature:
  - `max_steps` (the ladder ran out; auto-escalated): raise `MAX_XI_STEPS`, or `XI_START`, which shifts the ladder up.
  - `xi_cap` (the ladder hit `XI_MAX`, typically at very cold T): raise `XI_MAX` deliberately, or use a smaller `BETA_MAX` (a warmer cold end). Cost and memory grow with ξ; the largest grid is printed after each sweep.
  - `dvr_failed`: a scaled solve failed its convergence check. The potential may be non-confining or pathological in that range.
- **Quantum $C_v$ falls off at the hottest temperatures, or the Cv benchmark error rises there.** The base spectrum is truncating thermally accessible levels. Increase `NUM_STATES` (auto-escalated when $E_{\max}<$ `HOT_STATE_SAFETY`·$k_BT_{\text{hot}}$).
- **Quantum $C_v$ cuts off abruptly at the edge of the window.** Extend the temperature range in that direction.
- **Cv benchmark error rises at low T.** The lowest eigenvalues are inaccurate; check the base grid's resolution or span.
- **Relative-error metrics blow up at very low T.** For the quantum curve, numerator and denominator both underflow; use the absolute error. The classical curve never approaches zero (≥ ½).
- **A curve looks like straight segments rather than smooth.** Too few temperatures in that stretch; increase `N_BETA`.
- **A long run seems stuck at the cold end.** The largest-ξ solves are the most expensive (several seconds each in Sections 1 & 4, and longer on Section 6's refined grids). The progress bar's rate slows there, which is expected.
- **Section 6's classical error panel is not at round-off level.** Re-run it with `CLASSICAL_REFERENCE_PRECISE = True` to see whether the difference persists on the stricter 2/2 grid. If it does, the base grid of the scaled solves is not converged.
- **Quick scan verdict "unresolved" outside the hot tail.** The difference between the curves there is within the classical value's error at this resolution. Raise `QUICK_SCAN_RESOLUTION`, or zoom with `QUICK_SCAN_BETA_RANGE` onto that window (cheaper without the cold end) and raise it there.
- **Quick scan classical curve has gaps (†).** The ladder hit `XI_MAX` (`xi_cap`) or a solve failed. Narrow the window to warmer temperatures, or raise `XI_MAX` deliberately.
