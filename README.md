# Quantum-Classical Heat Capacity

A modular, system-agnostic numerical pipeline for computing the heat capacity $C_v(T)$ of a quantum particle in a smooth 1-D potential, and locating its classical limit, using the Discrete Variable Representation (DVR) method. The pipeline validates itself without ever depending on a closed-form solution: every result is checked against an independently generated, higher-resolution numerical reference rather than an analytic formula — the same check that will be needed once the project moves to potentials that have no closed-form answer at all.

**Status:** The harmonic oscillator (HO) is fully validated, reproducing the exact analytic result to machine precision. The pipeline is now being extended to a quartic double well.

For the full write-up — theory, validation strategy, results, and derivations — see [`docs/summaries/IEEE_Summary.tex`](docs/summaries/IEEE_Summary.tex) (the project's IEEE-style report). For a technical summary of the physics and how to read the diagnostic plots, see [`FINDINGS.md`](FINDINGS.md). For how the pipeline evolved over the semester, see [`HISTORY.md`](HISTORY.md). Meeting-by-meeting notes and worked derivations are in [`docs/summaries/Meetings_Summary.tex`](docs/summaries/Meetings_Summary.tex).

---

## Repository Structure

```
src/
├── config.py                     Single source of truth: potential, constants, control parameters
├── Quantum_HO_Master.py          Master driver — entry point, runs the full 6-section pipeline
├── Quantum_Classical_Combined.py General Cv pipeline (quantum Cv(T) + classical-limit scan)
├── Classical_Limit_Numerical.py  xi/n convergence engine (classical-limit search)
├── Cv_Numerical_Benchmark.py     Cv comparison: base grid vs. numerical reference
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
- `<category>` — one of `energy_levels` (also where `plot_potential.py`'s two potential-shape figures land — `potential_full_spectrum.png`, zoomed just enough to show every computed level, and `potential_zoomed.png`, zoomed tightly on the well's own minima so its shape is actually visible even at the cost of most levels falling outside the frame), `convergence`, `cv` (the quantum/classical Cv summary and every base-vs-reference or vs-analytic Cv benchmark, kept together rather than split by comparison source), or `dvr_limits` — matching the diagnostic categories in [`FINDINGS.md`](FINDINGS.md#reading-the-diagnostic-plots).

Every path component and figure filename is built only from `[A-Za-z0-9_-]` — no spaces, dots, commas, or `=` signs — so the whole `figures/` tree is safe to point `\graphicspath`/`\includegraphics` at directly, or copy wholesale into a LaTeX project's figures folder, with no renaming. A value's sign and decimal point survive as letters instead of being stripped (`-` → `m`, `.` → `p`), so `b=-0.5` and `b=0.5` still land in distinct folders (`b-m0p5` vs. `b-0p5`) rather than colliding.

This is implemented in [`src/figures/output_paths.py`](src/figures/output_paths.py). A driver script calls `set_context(SYSTEM_NAME, POTENTIAL_PARAMS)` once near the start of a run (`Quantum_HO_Master.py` and `plot_potential.py` both do this already); every plotting function then calls `save_figure(fig, category, name)` once its figure is finished. There is no `plt.show()` anywhere in the pipeline — `save_figure` saves the PNG and closes the figure immediately, so a run never blocks on a plot window and a long or bulk run never accumulates open figures. A re-run with identical parameters overwrites its own previous figures in place; a run with different parameters gets a fresh folder. A bulk-scan driver should call `set_context(...)` again at the top of each loop iteration, before that iteration's pipeline runs.

The hand-authored `fig_pipeline.png` workflow diagram (`src/figures/pipeline_diagram.py`) is not tied to any particular run and is saved directly to `figures/fig_pipeline.png`, unaffected by this scheme.

## Pipeline Overview

`Quantum_HO_Master.py` runs six sequential sections, all driven by the single potential defined in `config.py`:

1. **DVR base solve** — lowest `NUM_STATES` energy levels on an automatically configured grid, plus two potential-shape figures (`V(x)` with the computed spectrum overlaid: a full-spectrum overview and a version zoomed to actually show the well structure — see `plot_potential.py`).
2. **Numerical reference solve** — the same spectrum on an independently wider/finer grid, generated once and shared by every later step as the ground truth.
3. **Energy-level accuracy** — base vs. reference eigenvalues, absolute and relative error.
4. **Cv pipeline** — quantum $C_v(T)$ from the base spectrum, plus the numerical classical limit via the $\xi$/$n$-convergence scan.
5. **DVR limit analysis** — the DVR solver's own resolution ($\Delta x$) and level-count ($n$) breakdown points, measured against the reference.
6. **Cv numerical benchmark** — the full Cv pipeline re-run on the reference spectrum, closing the loop between eigenvalue accuracy and the final thermodynamic observable.

Only Section 0 of `config.py` needs editing to run on a new potential — no other file changes.

## Generalising to a New Potential

1. Edit `POTENTIAL_PARAMS` and `my_potential(x)` in `src/config.py` to the new $V(x)$ (must be finite everywhere — no hard walls). `my_potential` should read its coefficients from `POTENTIAL_PARAMS` rather than hard-coding them twice, since that dict is also what names each run's figure folder (see "Figure Output" above).
2. Update `SYSTEM_NAME` and `T_UNITS_LABEL` for plot labelling.
3. Adjust `NUM_STATES`, `BETA_MIN`/`BETA_MAX`, and `XI_START` for the new system's energy scale.
4. Run `src/Quantum_HO_Master.py` — the grid auto-configurator and reference generator adapt automatically.

## Acknowledgements

**Supervisor:** Dr. David Gelbwaser-Klimovsky
**Student:** Nadav Jean
