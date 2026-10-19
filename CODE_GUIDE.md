# How the Code Works

This file explains everything the code in `src/` does, why it does it that way, and how to run it. The code itself has only short comments; the explanations live here.

**The rule of this version:** you choose every numerical setting in `src/config.py`. The code never changes a setting by itself and never stops because a check failed. It computes, then **reports** where the results can be trusted. You decide what to change.

---

## 1. Quick start

From `src/`:

```bash
python Quick_Scan.py
```

```bash
python Quantum_HO_Master.py
```

- `Quick_Scan.py` is the **quick scan**: quantum and classical $C_v(T)$ of one potential (or a sweep of one parameter) at coarse settings, plus the trust report. Seconds to a few minutes.
- `Quantum_HO_Master.py` is the **full pipeline**: the same calculation at fine settings, plus every diagnostic (reference grid, level-by-level errors, DVR limits, convergence plots, coefficient sweep). It takes as long as the settings require; nothing is cut short.

A typical session:

1. Choose the potential and its parameters in `config.py` (section 1), and a grid that contains it (section 2).
2. Run the quick scan. Read the **suggestion** lines (what the grid should be) and the **trust report** (where each curve can be trusted).
3. Change `X_MIN`, `X_MAX`, `dx`, `num_states` or `xi_values` until the report says the temperatures you care about are trusted.
4. If the quick scan shows a candidate window where $C_v^{q} > C_v^{cl}$, run the full pipeline on it, with finer settings in `FULL_RUN`.

---

## 2. The physics, step by step

### 2.1 Heat capacity from a list of energy levels

For a particle with energy levels $E_n$ in contact with a heat bath at temperature $T$ (canonical ensemble, units $k_B = 1$):

$$
p_n = \frac{e^{-E_n/T}}{Z},\qquad Z = \sum_n e^{-E_n/T},\qquad U = \langle E\rangle = \sum_n p_n E_n .
$$

The heat capacity is $C_v = dU/dT$. Differentiating $U$ gives the fluctuation formula

$$
C_v = \frac{\langle E^2\rangle - \langle E\rangle^2}{T^2} = \frac{\mathrm{Var}(E)}{T^2}.
$$

`Classical_Limit_Numerical.compute_cv` does exactly this, with two small numerical safeguards:

- The Boltzmann factors are computed as $e^{-(E_n - E_{\min})/T}$. Shifting all energies by the same constant does not change any $p_n$, but it keeps the exponentials from overflowing at low $T$.
- The variance is computed as $\sum_n p_n (E_n - \langle E\rangle)^2$, not as $\langle E^2\rangle - \langle E\rangle^2$. The two are equal mathematically, but the second loses precision when the energies are large compared with their spread.

**Only the levels you give it are counted.** If a level that should be populated at temperature $T$ is missing, $C_v(T)$ is wrong. This is why the number of levels, `num_states`, matters (section 2.4).

### 2.2 Energy levels from the DVR

The energy levels come from the time-independent Schrödinger equation

$$
\left[-\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)\right]\psi = E\,\psi .
$$

The Discrete Variable Representation (DVR; Colbert & Miller, *J. Chem. Phys.* **96**, 1982 (1992)) turns this into a matrix eigenvalue problem on an evenly spaced grid $x_i = x_{\min} + i\,\Delta x$:

- **Potential energy:** diagonal, $H_{ii} \mathrel{+}= V(x_i)$. In the DVR the potential is just its value at each grid point.
- **Kinetic energy:** a full matrix that couples every pair of grid points. For a grid with open ends (`colbert_miller_dvr_1d`, Colbert & Miller eq. A7) it depends only on the distance $n = |i-j|$:

$$
T_{ii} = \frac{\hbar^2}{2m\Delta x^2}\,\frac{\pi^2}{3},\qquad
T_{ij} = \frac{\hbar^2}{2m\Delta x^2}\,\frac{2(-1)^{n}}{n^2}\quad(i\ne j).
$$

  Because it depends only on $|i-j|$, the code builds the first row and lets `scipy.linalg.toeplitz` copy it into the whole matrix.

- **With walls** (`box_dvr_1d`, Colbert & Miller eq. A6): the wavefunction must vanish at $x_{\min}$ and $x_{\max}$ (infinite walls). The grid points are strictly inside the box, $x_i = x_{\min} + iL/N$ for $i = 1,\dots,N-1$ with $L = x_{\max} - x_{\min}$, and

$$
T_{ii} = \frac{\hbar^2}{2m}\frac{\pi^2}{2L^2}\left[\frac{2N^2+1}{3} - \frac{1}{\sin^2(\pi i/N)}\right],\quad
T_{ij} = \frac{\hbar^2}{2m}\frac{\pi^2}{2L^2}(-1)^{i-j}\left[\frac{1}{\sin^2\!\frac{\pi(i-j)}{2N}} - \frac{1}{\sin^2\!\frac{\pi(i+j)}{2N}}\right].
$$

  The code fills this matrix one row at a time with a plain `for` loop.

- **The levels** are the eigenvalues of $H = T + V$. `scipy.linalg.eigvalsh` with `subset_by_index` returns only the lowest `num_levels` of them.

`solve_energy_levels(V, num_levels, x_min, x_max, dx, mass, hbar, walls)` is the one function the rest of the code calls. It turns `dx` into a number of grid points and calls one of the two solvers.

**When are the levels right?** Two physical conditions:

1. **Enough grid points per de Broglie wavelength.** A level with energy $E$ has local momentum $p(x) = \sqrt{2m(E - V(x))}$ and wavelength $\lambda = 2\pi\hbar/p$. The grid can only describe waves longer than about $2\Delta x$, so the highest levels (largest $p$) fail first when $\Delta x$ is too large. A narrow feature of $V$ (for example a dimple of width $s$) also needs $\Delta x \ll s$.
2. **(Open grid only) the grid ends must be where the particle never goes.** Without walls, the grid ends still act like walls. That is harmless only if $V$ is already so high there that the wavefunctions have died out. The highest levels reach furthest, so they fail first when the grid is too narrow.

With `WALLS = True` the walls are part of the physical system (a particle in a box), so condition 2 does not apply.

Checked on exact answers: the particle in a box and the harmonic oscillator come out to $10^{-13}$, and the two bound states of a Pöschl–Teller well ($-2$ and $-1/2$ for depth 3) to $10^{-14}$.

### 2.3 The classical limit: $\hbar \to \hbar/\xi$

The classical heat capacity is the $\hbar \to 0$ limit of the quantum one, with the same $V(x)$ and the same $T$. The code computes it the direct way: it solves the **same potential** with a **smaller Planck constant** $\hbar/\xi$, for every $\xi$ in the list `xi_values`, and computes $C_v(T)$ from those levels. As $\xi$ grows, the curves approach the classical $C_v(T)$.

**Why this is the same as Gelbwaser-Klimovsky's $\xi$-scaling.** Their prescription is $V \to \xi^2 V$ and $T \to \xi^2 T$. Divide the Schrödinger equation for $\xi^2 V$ by $\xi^2$:

$$
\left[-\frac{(\hbar/\xi)^2}{2m}\frac{d^2}{dx^2} + V\right]\psi = \frac{E}{\xi^2}\,\psi
\quad\Longrightarrow\quad E_n(\hbar,\ \xi^2V) = \xi^2\,E_n(\hbar/\xi,\ V)
$$

(their Eq. S7). In the Boltzmann factors, $E_n(\hbar,\xi^2V)/(\xi^2 T) = E_n(\hbar/\xi, V)/T$. So both give exactly the same $C_v$. The $\hbar/\xi$ form is simpler to read: one potential, one temperature, a smaller $\hbar$.

**Why each $\xi$ uses `dx/xi` and `num_states*xi` levels.**

- The de Broglie wavelength $2\pi\hbar/p$ shrinks by $\xi$, so the grid spacing shrinks by $\xi$ too. Every $\xi$ is then solved with the same number of points per wavelength as $\xi = 1$.
- The number of levels below an energy $E$ is, by Bohr–Sommerfeld, $N(E) \approx \frac{1}{\pi\hbar}\int p\,dx$, which grows as $1/\hbar$. Keeping `num_states*xi` levels covers the same energy window, so the same temperatures.

So if the trust checks pass at $\xi = 1$, the same grid quality holds at every $\xi$. The cost is that the matrix has $(N_{\text{points}}\,\xi)^2$ entries: memory grows as $\xi^2$ and time roughly as $\xi^3$. Measured on this machine: 8,000 grid points take about 45 s and 0.5 GB; 16,000 points would take about 6 min and 2 GB. The console prints the matrix size for every $\xi$ before it is solved.

**How far is the largest $\xi$ from the classical limit?** The leading quantum correction to the classical heat capacity is of order $\hbar^2$ (Wigner–Kirkwood expansion). So for large $\xi$,

$$
C_v(\xi) \approx C_v^{cl} + \frac{a}{\xi^2}.
$$

With the last two values $\xi_1 < \xi_2$ of the list, $C_v(\xi_1) - C_v(\xi_2) = a(1/\xi_1^2 - 1/\xi_2^2)$, and the distance still left at $\xi_2$ is

$$
\varepsilon = \frac{a}{\xi_2^2} = \frac{|C_v(\xi_2) - C_v(\xi_1)|}{(\xi_2/\xi_1)^2 - 1}.
$$

`classical_error_estimate` returns this $\varepsilon(T)$. The classical $C_v$ reported is simply $C_v$ at the largest $\xi$; nothing is extrapolated.

The $\hbar^2$ law only holds once $\hbar\omega/(\xi k_BT)$ is small. At low $T$ all the $\xi$ values can still be in the frozen quantum regime, where $C_v(\xi)$ is tiny and $\varepsilon$ is tiny too, although the value is far from classical. That is why the classical trust test is **relative**: $\varepsilon < \text{cv\_tol}\times C_v^{cl}$. A frozen curve, say $C_v(\xi_1) = 10^{-5}$ and $C_v(\xi_2) = 10^{-3}$, fails it ($\varepsilon/C_v \approx 0.3$), while a converged one ($C_v \approx 1$, $\varepsilon \sim 10^{-4}$) passes. Physically the classical $C_v$ of $p^2/2m + V(x)$ is never small: it is $\tfrac12 + \beta^2\mathrm{Var}(V) \ge \tfrac12$.

`plot_xi_convergence` (full pipeline) draws $C_v$ against $1/\xi^2$ at six temperatures. Near $1/\xi^2 = 0$ the $\hbar^2$ law makes the points fall on a straight line, and the dotted line through the last two points shows where it would end.

### 2.4 Trust checks: change one knob at a time

`Quantum_Classical_Combined.run` computes the quantum $C_v(T)$, then **recomputes it with one knob changed at a time** and records how much it moved at every temperature:

| check | what changes | what it detects |
|---|---|---|
| finer grid | `dx / REFERENCE_DX_FACTOR` | too few points per wavelength, or a narrow feature of V that is under-sampled |
| wider grid (open grid only) | width × `REFERENCE_SPAN_FACTOR`, same `dx` | grid ends inside the region the particle reaches |
| more levels | `num_states × REFERENCE_LEVELS_FACTOR`, same grid | populated levels missing from the sum |

Each check answers one physical question, so a failed check tells you which knob to change. The quantum $C_v$ is **trusted** at the temperatures where none of the changes moves it by more than `cv_tol`.

The classical curve uses the same grid quality at every $\xi$ (section 2.3), so it needs the quantum checks to pass too. On top of that, its $\varepsilon$ must be below `cv_tol` × $C_v^{cl}$.

Then, at every temperature, with $d = C_v^q - C_v^{cl}$ and uncertainty $u$ = (largest quantum change) + $\varepsilon$:

- **ABOVE:** both curves trusted and $d > \text{ERROR\_MARGIN}\times u$.
- **below:** both curves trusted and $d < -\text{ERROR\_MARGIN}\times u$.
- **undecided:** everything else. The report also lists the undecided temperatures where $C_v^q$ is above the current classical estimate. These are the places to look harder.

`ERROR_MARGIN = 2` because a classical value at finite $\xi$ is usually still a bit below its limit (for the harmonic oscillator $C_v(\xi)$ rises towards 1), which pushes $d$ up by about $\varepsilon$. Twice the uncertainty keeps that from being called an excess.

The checks never stop the run and never change a setting.

### 2.5 Grid suggestion (printed only)

`DVR_Algorithm.suggest_grid` estimates, from $V(x)$ alone, what the grid needs for temperatures up to `T_MAX`:

- **Highest energy that matters:** $E_{\text{top}} = V_{\min} + \text{THERMAL\_CUTOFF}\times k_BT_{\max}$. A level that high has Boltzmann weight below $e^{-20} \approx 2\times10^{-9}$ for the default cutoff of 20.
- **Region:** where $V(x) < E_{\text{top}}$, the classically allowed region at that energy. The grid must contain it, with room for the tails. The line is marked if $V$ is still below $E_{\text{top}}$ at a grid end.
- **Levels:** the Bohr–Sommerfeld count $\frac{1}{\pi\hbar}\int \sqrt{2m(E_{\text{top}} - V)}\,dx$.
- **dx:** four points per shortest wavelength, $\Delta x = \pi\hbar/(2p_{\max})$ with $p_{\max} = \sqrt{2m(E_{\text{top}} - V_{\min})}$.

It only looks at $V$ inside your grid and ignores narrow features. The trust checks are the real test.

---

## 3. Every setting in `config.py`

| setting | meaning | physical effect | cost |
|---|---|---|---|
| `MASS`, `HBAR` | particle mass, Planck constant | set the units (both 1 by default) | none |
| `POTENTIAL` | a function from `potentials.py` | the system | none |
| `POTENTIAL_PARAMS` | its parameters (a dict) | the system | none |
| `SYSTEM_NAME` | name for titles and the figure folder | none | none |
| `X_MIN`, `X_MAX` | grid ends | must contain everything the particle reaches at `T_MAX` (or are the box walls) | points ∝ width |
| `WALLS` | `True`: infinite walls at the grid ends | particle in a box + V; `False`: V confines the particle | none |
| `T_MIN`, `T_MAX` | temperature range ($k_BT$), log-spaced | which temperatures are computed | hotter needs more levels and a wider grid |
| `QUICK_SCAN["dx"]`, `FULL_RUN["dx"]` | grid spacing at $\xi = 1$ | points per wavelength | points ∝ 1/dx, time ∝ points³ |
| `["num_states"]` | levels kept at $\xi = 1$ | highest populated level | small (must be below the number of points) |
| `["xi_values"]` | the $\hbar/\xi$ values, ascending | how close the largest $\xi$ is to the classical limit | memory ∝ $\xi^2$, time ∝ $\xi^3$ |
| `["n_T"]` | number of temperatures | curve resolution | negligible |
| `["cv_tol"]` | trust tolerance on $C_v$ (units of $k_B$) | how strict the trust report is | none |
| `QUICK_SCAN_MODE` | `"single"` or `"sweep"` | which quick scan runs | sweep = one run per value |
| `REFERENCE_DX_FACTOR`, `REFERENCE_SPAN_FACTOR`, `REFERENCE_LEVELS_FACTOR` | how much each trust check changes its knob; the first two also build the reference grid | strictness of the checks | each check is one extra solve at $\xi = 1$ |
| `ERROR_MARGIN` | how many uncertainties $d$ must exceed | strictness of "ABOVE" / "below" | none |
| `THERMAL_CUTOFF` | for the grid suggestion only | none | none |
| `LIMIT_TOLERANCE` | energy tolerance for a "correct" level (full pipeline, Sections 3 and 5) | strictness | none |
| `LIMIT_DX_FACTORS`, `LIMIT_SPAN_FACTORS` | the grids tried in Section 5 | none | one solve each |
| `RUN_CLASSICAL_REFERENCE` | Section 6 also recomputes the classical curve on the reference grid | none | large: every $\xi$ solve gets `REFERENCE_SPAN_FACTOR × REFERENCE_DX_FACTOR` times more points |
| `SCAN_PARAM`, `SCAN_VALUES` | parameter to sweep and its values; `[]` skips Section 7 | none | one full calculation per value |

---

## 4. The files

```
src/
├── config.py                      every setting (section 3)
├── potentials.py                  library of V(x) functions
├── Quick_Scan.py                  quick scan (entry point)
├── Quantum_HO_Master.py           full pipeline (entry point)
├── Quantum_Classical_Combined.py  run(): one potential, quantum + classical + trust checks; report; Cv figures
├── Classical_Limit_Numerical.py   Cv from levels; classical limit (hbar/xi); error estimate; Cv vs number of levels
├── Cv_Coefficient_Sweep.py        run() for every value of one parameter; sweep figures
├── Cv_Numerical_Benchmark.py      Cv from the base grid vs from the reference grid
├── DVR/
│   ├── DVR_Algorithm.py           the two DVR solvers, solve_energy_levels, suggest_grid
│   ├── DVR_Reference_Generator.py the reference grid (wider and finer) and its spectrum
│   └── DVR_Limit_Finder.py        level errors for a list of dx values and grid widths
├── error/error_energylevels.py    level-by-level comparison of two spectra
├── analytical/
│   ├── HO_Analytical.py           exact harmonic-oscillator levels and Cv
│   └── HO_Benchmark.py            numerical vs exact (run by the full pipeline for potentials.harmonic)
└── figures/
    ├── output_paths.py            where figures are saved
    ├── plot_potential.py          V(x) with the levels drawn on it
    └── pipeline_diagram.py        workflow diagram (python figures/pipeline_diagram.py)
```

**Who calls what.**

```
Quick_Scan.py ──┬── DVR_Algorithm.suggest_grid
                ├── Quantum_Classical_Combined.run ──┬── DVR_Algorithm.solve_energy_levels
                │                                    ├── Classical_Limit_Numerical.heat_capacity_curve
                │                                    └── Classical_Limit_Numerical.classical_limit_curves
                ├── Quantum_Classical_Combined.print_report / plot_cv_summary
                └── Cv_Coefficient_Sweep.run_coefficient_sweep (sweep mode) ── run() for every value

Quantum_HO_Master.py ── the same, then
                        DVR_Reference_Generator, error_energylevels, DVR_Limit_Finder,
                        Cv_Numerical_Benchmark, Cv_Coefficient_Sweep, plot_potential
```

**Main functions.**

| function | file | what it does |
|---|---|---|
| `colbert_miller_dvr_1d` | DVR_Algorithm | levels on an open grid (section 2.2) |
| `box_dvr_1d` | DVR_Algorithm | levels with walls (section 2.2) |
| `number_of_grid_points` | DVR_Algorithm | grid points for a given dx |
| `solve_energy_levels` | DVR_Algorithm | the function everything uses: dx → grid → levels |
| `suggest_grid` | DVR_Algorithm | section 2.5 |
| `compute_cv`, `heat_capacity_curve` | Classical_Limit_Numerical | section 2.1, at one T / at every T |
| `classical_limit_curves` | Classical_Limit_Numerical | $C_v(T)$ for every $\xi$ (section 2.3) |
| `classical_error_estimate` | Classical_Limit_Numerical | $\varepsilon(T)$ (section 2.3) |
| `cv_vs_number_of_levels` | Classical_Limit_Numerical | $C_v$ when only the lowest n levels are kept |
| `run` | Quantum_Classical_Combined | everything for one potential (section 2.4); returns a dict |
| `print_report`, `print_grid_suggestion` | Quantum_Classical_Combined | console output |
| `plot_cv_summary`, `plot_xi_convergence`, `plot_n_convergence` | Quantum_Classical_Combined | figures |
| `run_coefficient_sweep` | Cv_Coefficient_Sweep | `run` for every value in `SCAN_VALUES` |
| `generate_reference_energies` | DVR_Reference_Generator | levels on the reference grid |
| `run_dvr_limit_analysis` | DVR_Limit_Finder | Section 5 |
| `run_cv_numerical_benchmark` | Cv_Numerical_Benchmark | Section 6 |

**What `run` returns** (a dict): `T_arr`; `energies` (the `num_states` levels) and `energies_all` (with the extra levels of the "more levels" check); `cv_quantum`; `curves` ({ξ: Cv(T)}) and `cv_classical` (the largest ξ); `checks` ({check name: |change in Cv(T)|}); `quantum_error` (largest change at each T); `classical_error` (ε); `quantum_ok` and `classical_ok` (True/False at each T); `difference` ($d$), `uncertainty` ($u$); `above`, `below`, `undecided` (True/False at each T); and the settings used.

---

## 5. The quick scan in detail

`python Quick_Scan.py` (or `--mode sweep`):

1. Prints the grid suggestion for `T_MAX` next to your settings.
2. Calls `run` with `QUICK_SCAN`:
   - one solve with `num_states × REFERENCE_LEVELS_FACTOR` levels; the first `num_states` give the quantum $C_v$, all of them give the "more levels" check;
   - one solve with the finer grid and, without walls, one with the wider grid;
   - one solve per $\xi$, printed as it runs with the number of points and the matrix size.
3. Prints the trust report (section 2.4) and saves one figure, `figures/<system>/<params>/quick_scan/quick_single.png`:
   - **top:** the quantum $C_v$ (blue), every $\xi$ curve (light to dark green), the classical estimate (dashed). Grey shading: quantum not trusted. Yellow: classical not converged. Red: quantum ABOVE classical;
   - **middle:** $d = C_v^q - C_v^{cl}$ with a band of ± `ERROR_MARGIN` × uncertainty;
   - **bottom:** every trust check against T, and $\varepsilon/C_v^{cl}$, with the `cv_tol` line. The curve that crosses the line tells you which knob to change.

In sweep mode it does step 2 for every value in `SCAN_VALUES` and saves `quick_sweep.png` (every variant's quantum and classical $C_v$, and $d$, solid where trusted) and `quick_sweep_potentials.png` (every $V(x)$ with its lowest levels).

Default settings, double well: 10–20 s for one potential (depending on the Python installation), about 1 min for the four-value sweep.

---

## 6. The full pipeline in detail

`python Quantum_HO_Master.py` runs, with `FULL_RUN`:

| section | what | figures (in `figures/<system>/<params>/`) |
|---|---|---|
| 1 & 4 | `run` (as in the quick scan) and the report; for `potentials.harmonic`, also the comparison with the exact formulas | `cv/cv_summary.png`, `convergence/xi_convergence.png`, `convergence/n_convergence.png`, `energy_levels/potential_*.png`, (`cv/analytic_cv_benchmark.png`) |
| 2 | levels on the reference grid: `REFERENCE_SPAN_FACTOR` wider, `REFERENCE_DX_FACTOR` finer (with walls only finer) | none |
| 3 | base levels vs reference levels | `energy_levels/energy_level_comparison.png` |
| 5 | the levels redone for every dx in `LIMIT_DX_FACTORS` × dx and every width in `LIMIT_SPAN_FACTORS`, each compared with the reference: how many levels are correct to `LIMIT_TOLERANCE` | `dvr_limits/level_errors_vs_dx.png`, `dvr_limits/level_errors_vs_span.png` |
| 6 | quantum $C_v$ from the base vs the reference levels (and the classical $C_v$ if `RUN_CLASSICAL_REFERENCE`) | `cv/quantum_cv_benchmark.png` (`cv/classical_cv_benchmark.png`) |
| 7 | the sweep of `SCAN_PARAM` over `SCAN_VALUES`, at the full settings | `coefficient_sweep/*.png` |

Section numbers are the old ones, so FINDINGS.md and the report still match. Sections 5 and 6 are plain loops: no bisection, no automatic search.

**The energy-level error is an energy, $|E_n - E_{n,\text{ref}}|$**, not a relative error. Boltzmann factors depend on $(E_n - E_0)/k_BT$, so an error in $E_n$ matters compared with $k_BT$, not with $E_n$ itself. $E_n$ can also be negative or zero (double well, dimple), which makes $|\Delta E_n|/E_n$ meaningless.

Default settings, double well: about 13 min in total. Sections 1 & 4 take 2.5 min (the largest matrix, at ξ = 32, has 11,201 points and needs 1 GB), Sections 2–6 a few seconds, and Section 7 about 10 min for the four values of `b`.

---

## 7. Adding a new potential

1. In `src/potentials.py`, add a function whose first argument is `x` and whose other arguments are the physical parameters:

   ```python
   def gaussian_well(x, depth, width):
       """V = -depth exp(-x^2 / (2 width^2))."""
       return -depth * np.exp(-x**2 / (2.0 * width**2))
   ```

2. In `src/config.py`:

   ```python
   POTENTIAL        = potentials.gaussian_well
   POTENTIAL_PARAMS = {"depth": 5.0, "width": 1.0}
   SYSTEM_NAME      = "1-D Gaussian well"
   ```

   This well flattens to 0 far away, so the particle is only bound below 0. Put it in a box: `WALLS = True`, with `X_MIN`, `X_MAX` the walls.

3. Run the quick scan and adjust the grid from the suggestion and the trust report.

`config.my_potential(x)` evaluates the chosen potential with `POTENTIAL_PARAMS`; `config.my_potential(x, params)` with other parameters (the sweep uses this).

---

## 8. Reading the report: what to change

| what the report says | physical cause | what to change |
|---|---|---|
| "finer grid" fails at high T | the highest populated levels have too few points per wavelength | smaller `dx` |
| "finer grid" fails at all T | a narrow feature of V is under-sampled | smaller `dx` (well below the feature width) |
| "wider grid" fails at high T | hot particles reach the grid ends | wider `X_MIN`, `X_MAX` (the suggestion gives the region) |
| "wider grid" fails everywhere | V does not confine the particle (it flattens or goes down far away) | `WALLS = True`, or a confining potential |
| "more levels" fails at high T | populated levels are missing | larger `num_states` |
| classical not converged at low T | $\hbar\omega/(\xi k_BT)$ is not small yet at the largest $\xi$ | add a larger $\xi$ (memory ∝ $\xi^2$), or accept that those temperatures are undecided |
| "undecided, of which Cv_q is above the classical estimate" | a possible excess, not confirmed | larger $\xi$ there, or a narrower T range to save time |
| `ValueError: ... levels need more than ... grid points` | more levels asked than the grid has points | larger grid or fewer levels |
| memory error | the largest $\xi$ matrix does not fit | smaller largest $\xi$, or larger `dx` (the console prints the size of each matrix) |

---

## 9. Limits of the method

- **Low temperatures need large $\xi$.** The classical limit at $T$ needs $\hbar\omega/(\xi k_BT) \ll 1$, where $\omega$ is the vibration frequency at the bottom of the well. Example: the Pöschl–Teller well of depth 1 in a box of width 100 (your analytic case, $\alpha L = 100$). Its excess window is at $k_BT \approx 0.1$, where $\hbar\omega/k_BT \approx 14$. With $\xi$ up to 16 the classical curve converges only above $T = 0.12$, so the window comes out **undecided, with $C_v^q$ above the classical estimate at T = 0.09–0.12**. The quantum side is right there: the DVR gives $C_v^q - C_v^{cl} = +0.343$ at $T = 0.105$ when $C_v^{cl}$ is taken from the exact classical integral, the same as the analytic $+0.34$. Converging the $\xi$ curve there would need $\xi \sim 100$, about 50,000 grid points: too big for a dense matrix.
- **Not implemented (your decision):**
  - **The exact classical $C_v$.** For $H = p^2/2m + V(x)$ the classical heat capacity is the 1-D integral $C_v^{cl} = \tfrac12 + \beta^2\mathrm{Var}(V)$, with weights $e^{-\beta V(x)}$ over the grid (CHANGES_SUMMARY.md, section F). It costs nothing and would settle the low-T region directly. It is not $\xi$-numerical, so it was left out.
  - **Richardson extrapolation** in $1/\xi^2$: report $C_v(\xi_2) - \varepsilon$ instead of $C_v(\xi_2)$. One line of code, but it only helps where the $\hbar^2$ law already holds.

---

## 10. What changed compared with the previous version

Removed, because they changed settings or stopped the run on their own:

- `auto_configure_dvr` (automatic grid) and `get_fully_converged_energy_levels` (3-pass check that raised an error). You set the grid; the trust checks report.
- The auto-tune loop (`Cv_AutoTune.py`, escalation of `NUM_STATES` and the $\xi$ ladder).
- The $\xi$-ladder plateau search (`converge_xi`, `XI_START`/`XI_MULT`/`MAX_XI_STEPS`/`XI_MAX`/`TOL_XI`/`MIN_STABLE_XI`) and its $C_v \ge \tfrac12$ guard. You give the list `xi_values`.
- `ScaledSpectra` (per-temperature WKB sizing and caching of the $\xi^2 V$ solves), `RememberedPotential`, `T_merge`, quick-scan resolution presets, the coefficient-sweep escalation and failure diagnosis, the formula templates (`POTENTIAL_FORMULA`), the interactive reference prompt, `DVR_DX_MAX`.
- The bisection searches of `DVR_Limit_Finder.py`, replaced by loops over lists you give.

Changed:

- $\xi^2V$ at $\xi^2T$ is now written as $\hbar/\xi$ at $T$ (the same numbers, section 2.3).
- Temperatures are given as `T_MIN`, `T_MAX` instead of `BETA_MIN`, `BETA_MAX`.
- The energy-level error is $|\Delta E_n|$ instead of $|\Delta E_n|/E_n$ (section 6).
- Potentials live in `potentials.py`, with explicit parameters.

Added: walls (`WALLS`, `box_dvr_1d`), the one-knob-at-a-time trust checks, the grid suggestion.

Kept, with the same names: the two entry points, the section numbers, `colbert_miller_dvr_1d` (same arguments), `compute_cv`, the figure folders.

Checked with this version:

| test | result |
|---|---|
| particle in a box, walls | levels exact to $3\times10^{-13}$ |
| harmonic oscillator, open grid and inside a wide box | levels exact to $10^{-13}$; with $\hbar/\xi$ ($\xi$ = 2, 4) to $3\times10^{-13}$ |
| harmonic oscillator, full pipeline | quantum $C_v$ within $4\times10^{-5}$ of the Einstein formula; classical within $10^{-3}$ of 1 where reported converged |
| Pöschl–Teller, depth 3, walls at ±30 | bound states $-2$, $-1/2$ to $10^{-14}$; box states independent of dx and on $kL + \delta(k) = j\pi$ |
| Pöschl–Teller, depth 1, box of width 100 | $C_v^q - C_v^{cl} = +0.343$ at $T = 0.105$ (analytic: $+0.34$) |
| double well, $b = -0.9 \dots 0$ | quantum below classical wherever both are trusted, as before |
| Python 3.7 (Anaconda) and 3.12 (`C:\Technion\myenv`) | both run |
