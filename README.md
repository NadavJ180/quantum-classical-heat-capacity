# Quantum-Classical Heat Capacity

Numerical heat capacity $C_v(T)$ of a quantum particle in a 1-D potential $V(x)$, compared with its classical limit. The goal is to find potentials and temperatures where the quantum $C_v$ is **above** the classical one.

- **Quantum $C_v(T)$:** from energy levels computed with the Discrete Variable Representation (DVR), on an open grid or between infinite walls.
- **Classical limit:** the same potential solved with a smaller Planck constant, ħ → ħ/ξ, for a list of ξ values you choose. This is the ξ-scaling of Gelbwaser-Klimovsky et al. (V → ξ²V, T → ξ²T), written at fixed V and T.
- **You control every setting** in [`src/config.py`](src/config.py). The code never changes a setting by itself; it reports where each curve can be trusted, by redoing the calculation with one setting changed at a time.

**How everything works, and why:** [`CODE_GUIDE.md`](CODE_GUIDE.md).

## Running

From `src/`:

```bash
python Quick_Scan.py
```

```bash
python Quantum_HO_Master.py
```

| | what | time (double well, default settings) |
|---|---|---|
| `Quick_Scan.py` | quantum and classical $C_v(T)$ at coarse settings, the trust report, one figure; `--mode sweep` does it for every value of one parameter | 10–20 s per potential |
| `Quantum_HO_Master.py` | the same at fine settings, plus the reference grid, level-by-level errors, DVR limits, convergence plots and the coefficient sweep | ~13 min (2.5 min without the sweep) |

**A new potential:** add a function to [`src/potentials.py`](src/potentials.py) and point `POTENTIAL` and `POTENTIAL_PARAMS` in `config.py` at it ([`CODE_GUIDE.md`](CODE_GUIDE.md), section 7).

**Tested environments:** Anaconda Python 3.7 (numpy 1.18, scipy 1.7, matplotlib 3.1) and Python 3.12 (numpy 2.4, scipy 1.17, matplotlib 3.10).

## Status

- Harmonic oscillator, particle in a box and Pöschl–Teller bound states: DVR levels exact to ~1e-13.
- Harmonic oscillator: quantum $C_v$ within 4e-5 of the Einstein formula; classical $C_v$ within 1e-3 of 1 wherever it is reported converged.
- Pöschl–Teller well in a box ($\alpha L = 100$): $C_v^q - C_v^{cl} = +0.343$ at $k_BT = 0.105$, the same as the analytic result (classical side from the exact integral; the ξ curves do not converge that cold, see [`CODE_GUIDE.md`](CODE_GUIDE.md), section 9).
- Quartic/cubic/quadratic double well: the quantum $C_v$ stays below the classical one wherever both are trusted ([`FINDINGS.md`](FINDINGS.md)).

## Repository Structure

```
src/
├── config.py                      every setting
├── potentials.py                  library of V(x) functions
├── Quick_Scan.py                  quick scan (entry point)
├── Quantum_HO_Master.py           full pipeline (entry point)
├── Quantum_Classical_Combined.py  one potential: quantum + classical Cv + trust checks; report; Cv figures
├── Classical_Limit_Numerical.py   Cv from energy levels; classical limit with hbar/xi; its error estimate
├── Cv_Coefficient_Sweep.py        the calculation for every value of one parameter
├── Cv_Numerical_Benchmark.py      Cv from the base grid vs the reference grid
├── DVR/
│   ├── DVR_Algorithm.py           DVR solvers (open grid, walls), grid suggestion
│   ├── DVR_Reference_Generator.py wider and finer reference grid
│   └── DVR_Limit_Finder.py        level errors for lists of grid spacings and widths
├── error/error_energylevels.py    level-by-level comparison of two spectra
├── analytical/                    exact harmonic-oscillator results and the comparison with them
└── figures/                       figure paths, potential plot, workflow diagram

audit/, homog/   checks of earlier versions (they use functions this version removed; see CODE_GUIDE.md)
figures/         generated plots: figures/<system>/<params>/<category>/<name>.png
docs/summaries/  IEEE_Summary.tex (report), Meetings_Summary.tex (meeting notes, derivations)
```

Other documents:
- [`CODE_GUIDE.md`](CODE_GUIDE.md): the code, the physics behind each step, every setting, how to read the output.
- [`FINDINGS.md`](FINDINGS.md): physics background and results so far.
- [`HISTORY.md`](HISTORY.md): how the code got here.
- [`CHANGES_SUMMARY.md`](CHANGES_SUMMARY.md): the classical-limit correction of October 2026, change by change.

## Figure Output

Every figure is saved to `figures/<system>/<params>/<category>/<name>.png`. `<system>` comes from `SYSTEM_NAME` and `<params>` from `POTENTIAL_PARAMS` (e.g. `a-0p25_b-m0p5_c-m0p5_d-0`), so runs with different parameters never overwrite each other. Categories: `quick_scan`, `cv`, `convergence`, `energy_levels`, `dvr_limits`, `coefficient_sweep`. Names use only `[A-Za-z0-9_-]` (`-` → `m`, `.` → `p`), so LaTeX can include them directly.

## Acknowledgements

**Supervisor:** Dr. David Gelbwaser-Klimovsky
**Student:** Nadav Jean
