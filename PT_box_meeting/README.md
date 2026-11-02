# PT well in a box — material for the 28.10.2026 meeting

A reflectionless Pöschl–Teller well (λ = 2) inside an infinite square well: an exactly solvable case of Cv_q > Cv_cl.

## Files

| File | What it is |
|---|---|
| `PT_Box_Summary.pdf` | The 2-page handout: the assumptions, then only equations, short derivation steps and page-specific references (with pointers "D§n" into the derivation). |
| `PT_Box_Derivation.pdf` | Every derivation behind the handout, step by step, with "Check" boxes, the mechanism heuristics, and page-specific references. For you to go over first. |
| `PT_Box_Explained.pdf` | Standalone, book-style version of the Cv_q and Cv_cl derivation for a reader with only Calculus 1–2 and linear algebra: a math and physics toolkit first (F1–F12, P1–P5), then every step, each marked "Source" (page-checked reference) or "Own step" (constructed by Claude; "no reference found" where none was found). §4.7 derives the levels a second way (associated Legendre functions) and §4.8 compares the two routes and recommends one. Strokes are not covered. |
| `check_explained.py` | Recomputes every number in `PT_Box_Explained.pdf` (checks C1–C19, cited in the text) and draws `figures/fig_explained_levels.pdf`. About 30 s; output in `results/check_explained.log`. |
| `cv_plot.py` | Plots Cv_q and Cv_cl vs T for several (alpha, L) at once. Edit the `CONSTANTS` block at the top. |
| `cv_interactive.py` | The same curves with sliders for alpha, L and gamma (the scale-invariant stroke), next to plots of the potential. |
| `check_formulas.py` | Prints every number quoted in the two PDFs (closed forms vs exact level sums, vs quadrature, stroke scaling). About 10 s. |
| `make_figures.py` | Rebuilds `figures/` for the PDFs. |
| `dvr_poly_well.py` | The same well in a smooth polynomial trap instead of hard walls, computed with the repo's DVR pipeline (`src/Quick_Scan.py`, unchanged). Writes `results/` and two figures: `figures/fig_dvr_poly_well.png` (the anomaly window) and `figures/fig_dvr_poly_well_wide.png` (all temperatures the saved levels cover). See below. |

## Running

```
C:\Technion\myenv\Scripts\python.exe cv_plot.py
C:\Technion\myenv\Scripts\python.exe cv_interactive.py
```

Run them from this folder: the last three scripts import the formulas from `cv_plot.py`. They also run on the Anaconda 3.7 `python`.

To rebuild a PDF after editing the `.tex`: `pdflatex PT_Box_Summary.tex` (twice for the derivation, to update cross-references).

## Using the plots

- **Units.** `HBAR = M = KB = 1`, and the T axis is plain k_B T. Only alpha*L and k_B T / (hbar^2 alpha^2 / 2m) matter; the console prints both temperatures for the maximum of Cv_q − Cv_cl.
- **To see the anomaly form:** in `cv_interactive.py`, start at L ≈ 30 (no window) and slide L up to 1000 or more. The shaded window first appears between αL = 70 and 80 (max Cv_q − Cv_cl = −0.04 k_B at 70, +0.10 k_B at 80) and keeps growing.
- **Exact values:** the box to the right of each slider shows its value. Click it, type a number (for example `1000` for L) and press Enter. A value outside the slider's range widens the slider; anything that is not a positive number is ignored.
- **Reset to defaults** (button at the bottom left): back to α = 1, L = 100, γ = 1, the original slider ranges, and all check boxes off.
- **The potential** (right-hand panels of `cv_interactive.py`): the whole box on top, and a zoom on the well below, where the thick bar is the bound state at −ħ²α²/2m. They follow all three sliders.
- **To see homogeneity:** move the γ slider. It applies the stroke α → α/γ **and** L → γL (so αL, the shape, is unchanged). The curves slide along T by 1/γ² without changing shape, and the potential gets wider and shallower. The grey curves and potential are the system before the stroke. Then tick **"plot vs γ²T"**: it plots the after-stroke curves against γ²T (and the potential as γ²V against x/γ). The coloured curves then lie exactly on the grey ones, which is the homogeneous scaling E_n → E_n/γ².
- **Empty box only** (`SHOW_BOX_ONLY` in `cv_plot.py`, check box in `cv_interactive.py`): the green curves are the box *without* the well. In the Cv panel they are its Cv_q (dotted) and Cv_cl = k_B/2 (dash-dot); in the difference panel, its Cv_q − Cv_cl. This is the wall effect, at most k_B/16, and it is what remains at high T. `cv_plot.py` also opens a second figure with the box alone against T/E_1.
- **Why αL and not L:** the well's only length is its width 1/α, so L can only matter through αL, the box length in well widths. Two systems with the same αL have identical Cv curves against k_B T/(ħ²α²/2m). See the summary's assumptions box and Sec. 1 of the derivation.

### Removing the box-only option

It is self-contained. Delete:
- in `cv_plot.py`: the block between `BOX-ONLY (delete from here ... )` and `... to here`, the three short blocks marked `BOX-ONLY (delete these ... lines ...)`, and the `SHOW_BOX_ONLY` constant;
- in `cv_interactive.py`: the block between `BOX-ONLY (delete this block ...)` and `end BOX-ONLY`, `BOX_ONLY` from the `LABELS` list (and from the line that defines it), and `cv_quantum_box_only` from the import line.

## DVR without walls (`dvr_poly_well.py`)

The hard walls are replaced by a polynomial:

    V(x) = (ħ²α²/2m)(2x/L)^10 − (ħ²α²/m)/cosh²(αx)

The polynomial is close to 0 for |x| < L/2 and reaches the binding energy ħ²α²/2m at x = ±L/2, so L still plays the role of the box length. Nothing here is analytic:
- the quantum Cv comes from the DVR spectrum of V;
- the classical Cv is the pipeline's ξ-scan (ħ → ħ/ξ, with ξ²V re-solved by the DVR at every ξ).

`src/` is not changed. The script sets the same names a potential block of `config.py` sets (`POTENTIAL_PARAMS`, `my_potential`, `SYSTEM_NAME`, ...) on the imported config module, then calls `Quick_Scan.run_quick_scan`. It also sets `DVR_DX_MAX = 0.25/α`, because the well (width 1/α) is narrower than the DVR's own grid spacing at ξ = 1.

```
C:\Technion\myenv\Scripts\python.exe dvr_poly_well.py               # quick scan (~70 min), then the figure
C:\Technion\myenv\Scripts\python.exe dvr_poly_well.py --plot-only   # figure only, from results/
```

- **Settings** (top of the script): α = 1, L = 100 (αL = 100), polynomial order 10, resolution 1, k_B T = 0.08–0.13.
- **The window is narrow because of cost.** The classical side needs ξ up to 86.5 there. On a trap 100 well-widths wide, that means DVR matrices of up to 13,595 points. The quick run took 70 min (12 DVR solves) and an estimated 5 GB of memory.
- **Two figures.** `fig_dvr_poly_well.png` is the quick-scan window (k_B T = 0.08–0.13). `fig_dvr_poly_well_wide.png` covers k_B T = 0.01–0.167 and is drawn from the same saved results, without a new DVR run:
  - Cv_q comes from the 59 saved DVR levels (up to E = 1.50). They reproduce the scan's Cv_q to 4 × 10⁻¹⁴. By the scan's own rule (levels up to E_0 + 12 k_B T) they cover k_B T ≤ 0.167, which sets the hot end.
  - Cv_cl is the exact classical integral, because the ξ-scan exists only in the window; its points are drawn on top there.
  - Going hotter than 0.167 needs more levels: one quantum-only DVR solve (about 500 levels for k_B T ≤ 3, a few seconds). The ξ-scan up there would take about 3.5 h.
- **The pipeline's own figure** is also saved, under `figures/1_d_pt_well_in_x10_trap/` at the repo root.
- **The dotted check line** is the exact classical Cv, 1/2 + β²Var(V) over e^{−βV}. It is not part of the pipeline and can be switched off with `SHOW_EXACT_CLASSICAL`.
- **Control:** the order-10 trap alone (no well) never exceeds its classical Cv of 1/2 + 1/10, so the excess comes from the well.

### Result of the quick run (28.10.2026)

- **Cv_q > Cv_cl for k_B T ≈ 0.093–0.121.** This is resolved at 33 of the 60 temperatures, and the curves cross at 0.0925 and 0.1216. The largest excess is +0.259 k_B at k_B T = 0.108 (Cv_q = 2.333, Cv_cl = 2.074).
- **The box gives nearly the same window.** The closed form at αL = 100 has Cv_q > Cv_cl between 0.088 and 0.120, with a maximum of +0.343 k_B at 0.105. The soft polynomial walls give the same anomaly with a somewhat smaller peak.
- **The bound state** from the DVR is E = −0.5000000000, the analytic −ħ²α²/2m.
- **Check of the classical side.** The ξ-scan lies below the exact classical integral at every temperature, by 2.8–4.5 × 10⁻³. That is the ħ²/ξ² correction left at the last ξ. Its own error estimate has the right size except at k_B T = 0.0935–0.0983, where it is 2–6 times too small (the console prints this). Measured against the exact classical value instead, the crossings move by 10⁻⁴ (to 0.0926 and 0.1215) and the peak is +0.255 instead of +0.259.
- **The pipeline's one "unresolved" temperature (k_B T = 0.122)** is where Cv_q − Cv_cl changes sign at the hot edge. The difference is zero there, so no resolution can resolve it.
- **Without zooming in** (the wide figure), the anomaly sits between two larger effects:
  - At cold T, Cv_q is frozen at 0, because the bound state lies 0.50 below the first trap level. The classical Cv is 1, from the harmonic bottom of the well.
  - The quantum curve rises first, which opens the window. It peaks at 3.52 at k_B T = 0.149.
  - The classical curve rises more steeply and peaks at 6.07 at k_B T = 0.176, just past the figure's end.

### Cost of the higher resolutions (same window)

Estimated from the engine's sizing rules; the estimate for resolution 1 (63 min, largest matrix 13,701 points) matched the run.

| Resolution | DVR solves | ξ up to | Largest matrix | Time | Peak memory |
|---|---|---|---|---|---|
| 1 (measured) | 12 | 86.5 | 13,595 | 70 min | ~5 GB |
| 2 | 17 | 122 | ~20,500 | ~5.6 h | ~11 GB |
| 3 (the full pipeline's settings) | 25 | 212 | ~44,600 | ~80 h | ~54 GB |

This machine has 16 GB of memory.

## Validity

The closed forms drop terms of order e^{−αL}, which are negligible for αL ≳ 10. For αL ≥ 60 they agree with sums over the exact levels to about 2e-13 k_B. The classical integral stays accurate down to k_B T ≈ depth / 20000.
