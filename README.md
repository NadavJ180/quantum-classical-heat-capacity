# Quantum-Classical Heat Capacity

A modular, system-agnostic numerical pipeline for computing the heat capacity $C_v(T)$ of a quantum particle in a smooth 1-D potential, and locating its classical limit, using the Discrete Variable Representation (DVR) method. The pipeline validates itself without ever depending on a closed-form solution: every result is checked against an independently generated, higher-resolution numerical reference rather than an analytic formula — the same check that will be needed once the project moves to potentials that have no closed-form answer at all.

**Status:** The harmonic oscillator (HO) is fully validated, reproducing the exact analytic result to machine precision. The pipeline is now being extended to a quartic double well.

For the full write-up — theory, validation strategy, results, and derivations — see [`docs/summaries/IEEE_Summary.tex`](docs/summaries/IEEE_Summary.tex) (the project's IEEE-style report). For a technical summary of the physics and how to read the diagnostic plots, see [`FINDINGS.md`](FINDINGS.md). For how the pipeline evolved over the semester, see [`HISTORY.md`](HISTORY.md). Meeting-by-meeting notes and worked derivations are in [`docs/summaries/Meetings_Summary.tex`](docs/summaries/Meetings_Summary.tex).

---

## Repository Structure

```
src/
├── config.py                     Single source of truth: potential, constants, control parameters
├── Quantum_HO_Master.py          Master driver — entry point, runs the full 7-section pipeline
├── Quantum_Classical_Combined.py General Cv pipeline (quantum Cv(T) + classical-limit scan)
├── Classical_Limit_Numerical.py  xi/n convergence engine (classical-limit search)
├── Cv_Numerical_Benchmark.py     Cv comparison: base grid vs. numerical reference
├── Cv_AutoTune.py                BETA_MIN auto-fill + NUM_STATES/XI_START escalation diagnostics (see "Auto-Tuning" below)
├── Cv_Coefficient_Sweep.py       Quantum Cv(T) comparison across a swept POTENTIAL_PARAMS coefficient (see "Coefficient Sweep" below)
├── DVR/
│   ├── DVR_Algorithm.py          Core DVR solver and automatic grid configuration
│   ├── DVR_Reference_Generator.py  Numerical reference grid generator
│   └── DVR_Limit_Finder.py       DVR accuracy limit searches (resolution and level-count)
├── analytical/
│   ├── HO_Analytical.py          Closed-form HO energy levels and Cv(T)
│   └── HO_Benchmark.py           External benchmark: numerical pipeline vs. analytic HO
├── error/
│   └── error_energylevels.py     Energy-level comparison (base vs. reference)
└── figures/
    ├── output_paths.py           Figure save-path resolver (see "Figure Output" below)
    ├── plot_potential.py         Potential-shape figure generator (called from Quantum_HO_Master.py's
    │                             Section 1 automatically; run standalone only if you want just this figure)
    └── pipeline_diagram.py       Workflow-diagram generator

figures/            Generated plots -- see "Figure Output" below for how a run's plots are
                    organized under here. HO/ and SymmetricDoubleWell/ [name predates the
                    b != 0 run in config.py, so that potential is currently asymmetric] are
                    older, hand-saved figures kept because docs/summaries/IEEE_Summary.tex
                    references them directly; new runs no longer write into either folder.
docs/summaries/      IEEE_Summary.tex (main report), Meetings_Summary.tex (meeting notes + derivations)
docs/SchottkyAnomaly/  Reference literature (Schottky-anomaly papers)
```

## Figure Output

Every figure the pipeline saves is written under:

```
figures/<system>/<params>/<category>/<name>.png
```

- `<system>` — a slug of `SYSTEM_NAME` (e.g. `1_d_asymmetric_double_well`).
- `<params>` — a slug of `POTENTIAL_PARAMS` (e.g. `a-0p25_b-m0p5_c-m0p5_d-0` for `{"a": 0.25, "b": -0.5, "c": -0.5, "d": 0.0}`), so that two runs of the *same* potential with *different* parameters never collide or overwrite each other. This is what makes a future bulk scan — e.g. sweeping the double well's `b` over a range of values — safe to run unattended: each parameter combination lands in its own folder automatically, with no manual bookkeeping.
- `<category>` — one of `energy_levels` (also where `plot_potential.py`'s two potential-shape figures land — `potential_full_spectrum.png`, zoomed just enough to show every computed level, and `potential_zoomed.png`, zoomed tightly on the well's own minima so its shape is actually visible even at the cost of most levels falling outside the frame), `convergence`, `cv` (the quantum/classical Cv summary, every base-vs-reference or vs-analytic Cv benchmark, and the coefficient-sweep comparison `cv_coefficient_sweep.png` — kept together rather than split by comparison source), or `dvr_limits` — matching the diagnostic categories in [`FINDINGS.md`](FINDINGS.md#reading-the-diagnostic-plots).

Every path component and figure filename is built only from `[A-Za-z0-9_-]` — no spaces, dots, commas, or `=` signs — so the whole `figures/` tree is safe to point `\graphicspath`/`\includegraphics` at directly, or copy wholesale into a LaTeX project's figures folder, with no renaming. A value's sign and decimal point survive as letters instead of being stripped (`-` → `m`, `.` → `p`), so `b=-0.5` and `b=0.5` still land in distinct folders (`b-m0p5` vs. `b-0p5`) rather than colliding.

This is implemented in [`src/figures/output_paths.py`](src/figures/output_paths.py). A driver script calls `set_context(SYSTEM_NAME, POTENTIAL_PARAMS)` once near the start of a run (`Quantum_HO_Master.py` and `plot_potential.py` both do this already); every plotting function then calls `save_figure(fig, category, name)` once its figure is finished. There is no `plt.show()` anywhere in the pipeline — `save_figure` saves the PNG and closes the figure immediately, so a run never blocks on a plot window and a long or bulk run never accumulates open figures. A re-run with identical parameters overwrites its own previous figures in place; a run with different parameters gets a fresh folder. A bulk-scan driver should call `set_context(...)` again at the top of each loop iteration, before that iteration's pipeline runs.

The hand-authored `fig_pipeline.png` workflow diagram (`src/figures/pipeline_diagram.py`) is not tied to any particular run and is saved directly to `figures/fig_pipeline.png`, unaffected by this scheme.

## Pipeline Overview

`Quantum_HO_Master.py` runs seven sequential sections, all driven by the single potential defined in `config.py`:

1. **& 4. DVR base solve + Cv pipeline (auto-tuned)** — solves for the lowest `NUM_STATES` energy levels on an automatically configured grid, then runs the full quantum $C_v(T)$ + $\xi$/$n$-convergence classical-limit sweep, in a closed loop: after each attempt, `Cv_AutoTune.diagnose_escalation` inspects the sweep's own convergence diagnostics for the two known truncation artifacts (see "Auto-Tuning" below) and grows `NUM_STATES` and/or `XI_START`/`MAX_XI_STEPS` before retrying, up to `MAX_ESCALATION_ROUNDS`. Also saves the two potential-shape figures (`V(x)` with the computed spectrum overlaid — see `plot_potential.py`) once the loop settles.
2. **Numerical reference solve** — the same spectrum on an independently wider/finer grid, generated once and shared by every later step as the ground truth.
3. **Energy-level accuracy** — base vs. reference eigenvalues, absolute and relative error.
5. **DVR limit analysis** — the DVR solver's own resolution ($\Delta x$) and level-count ($n$) breakdown points, measured against the reference.
6. **Cv numerical benchmark** — the full Cv pipeline re-run on the reference spectrum, closing the loop between eigenvalue accuracy and the final thermodynamic observable.
7. **Coefficient sweep** — the quantum $C_v(T)$ curves of several variants of the base potential (differing only in one named `POTENTIAL_PARAMS` coefficient) plotted against the base run's classical-limit curve — see "Coefficient Sweep" below.

(Sections are numbered to match their original six-section role; 1 and 4 are merged into one auto-tuned loop rather than run back-to-back as fixed, single-shot steps.)

Only Section 0 of `config.py` needs editing to run on a new potential — no other file changes.

## Auto-Tuning

The temperature range (`BETA_MAX`, and optionally `BETA_MIN`) is meant to be the only knob you touch run to run. Two numerical control parameters that the correctness of $C_v(T)$ actually depends on — `NUM_STATES` and the $\xi$-scan (`XI_START`/`MAX_XI_STEPS`) — are handled automatically instead:

- **`BETA_MIN`** — leave it as `None` (the default) and Section 1/4's loop derives it from `T_max = 10 \Delta E / k_B`, where $\Delta E = E_1 - E_0$ is the spectrum's fundamental gap: a practical "$T\to\infty$" checkpoint by which the classical-limit plateau should already be reached (not a hard limit — see [`Cv_AutoTune.py`](src/Cv_AutoTune.py)). Set it to a float instead to hand-pick the hot end, e.g. to zoom into a specific temperature window.
- **`NUM_STATES`/`XI_START`/`MAX_XI_STEPS`** — the values in `config.py` are only a round-0 starting guess. After each attempt, `diagnose_escalation` checks the sweep's own diagnostics on the matching half of the temperature range (hot half → `NUM_STATES`, per FINDINGS.md's "increase NUM_STATES" guidance for a truncated partition sum / numerical Schottky anomaly at high T; cold half → `XI_START`/`MAX_XI_STEPS`, per its "increase XI_START" guidance for a classical limit that drops because the $\xi$-scan ran out of budget at low T) and grows the relevant knob(s) before retrying.
- Escalation is bounded by `MAX_ESCALATION_ROUNDS` (default 4); if it's still failing at the last round, a `UserWarning` is printed and the pipeline proceeds with the last attempt rather than looping forever. Growth factors (`NUM_STATES_GROWTH`, `XI_START_GROWTH`, `MAX_XI_STEPS_GROWTH`, `NUM_STATES_CAP`) and the escalation trigger thresholds (`HOT_STATE_SAFETY`, `ESCALATION_FRACTION_THRESHOLD`) are all in `config.py`, meant to be touched rarely if ever.
- Set `AUTO_ESCALATE = False` to disable the loop entirely and run a single pass with the config's starting values, as before.

## Coefficient Sweep

To see how one coefficient of the current potential affects the shape of $C_v(T)$ (e.g. the double well's cubic term `b` and the size of its Schottky-anomaly-like bump), Section 7 sweeps that coefficient and plots every variant's quantum $C_v(T)$ against the base run's classical limit, in one figure (`cv/cv_coefficient_sweep.png`) — no xi/n-convergence search is redone per variant, and no other diagnostics are added to this plot.

- `SCAN_PARAM` — which key of `POTENTIAL_PARAMS` to vary (must be a real key of that dict).
- `SCAN_STEP` — spacing between consecutive variants.
- `SCAN_COUNT` — how many *extra* variants to add on each side of the value already in `config.py`, so the base potential is always included as one of the curves. Total curves plotted = `2*SCAN_COUNT + 1`.

Note: for a potential where only odd-degree terms break the $x\to-x$ symmetry (like the double well's `b x^3`), coefficient values equidistant from 0 in opposite signs (`+b`, `-b`) give mirror-image potentials with *identical* energy spectra and therefore identical $C_v(T)$ curves — if your swept range straddles such a pair, one curve will sit exactly underneath the other. This is real physics, not a bug.

## Generalising to a New Potential

1. Edit `POTENTIAL_PARAMS` and `my_potential(x)` in `src/config.py` to the new $V(x)$ (must be finite everywhere — no hard walls). `my_potential` should read its coefficients from `POTENTIAL_PARAMS` rather than hard-coding them twice, since that dict is also what names each run's figure folder (see "Figure Output" above).
2. Update `SYSTEM_NAME` and `T_UNITS_LABEL` for plot labelling.
3. Adjust `NUM_STATES` and `BETA_MAX` for the new system's energy scale — leave `BETA_MIN` as `None` unless you want to hand-pick the hot end (see "Auto-Tuning" above), and leave `XI_START` at its default; both are escalated automatically if the sweep shows a truncation artifact.
4. Optionally set `SCAN_PARAM`/`SCAN_STEP`/`SCAN_COUNT` to whichever coefficient you want the Section 7 comparison plot to vary.
5. Run `src/Quantum_HO_Master.py` — the grid auto-configurator, reference generator, and auto-tune loop all adapt automatically.

## Acknowledgements

**Supervisor:** Dr. David Gelbwaser-Klimovsky
**Student:** Nadav Jean
