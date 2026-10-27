# PT well in a box — material for the 28.10.2026 meeting

A reflectionless Pöschl–Teller well (λ = 2) inside an infinite square well: an exactly solvable case of Cv_q > Cv_cl.

## Files

| File | What it is |
|---|---|
| `PT_Box_Summary.pdf` | The 2-page handout: potential, levels, closed-form Cv_q and Cv_cl, the heuristic for why Cv_q > Cv_cl, and homogeneous strokes. |
| `PT_Box_Derivation.pdf` | Every derivation behind the handout, step by step, with "Check" boxes and a list of references to confirm (Sec. 8). For you to go over first. |
| `cv_plot.py` | Plots Cv_q and Cv_cl vs T for several (alpha, L) at once. Edit the `CONSTANTS` block at the top. |
| `cv_interactive.py` | The same curves with sliders for alpha, L and gamma (the scale-invariant stroke). |
| `check_formulas.py` | Prints every number quoted in the two PDFs (closed forms vs exact level sums, vs quadrature, stroke scaling). About 10 s. |
| `make_figures.py` | Rebuilds `figures/` for the PDFs. |

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
- **To see homogeneity:** move the γ slider. The curves slide along T (by 1/γ²) without changing shape. The grey curves are the system before the stroke.
- **Empty box only** (`SHOW_BOX_ONLY` in `cv_plot.py`, check box in `cv_interactive.py`): the dotted curve is Cv_q − Cv_cl of the box *without* the well. It is the wall effect, at most k_B/16, and it is what remains at high T. `cv_plot.py` also opens a second figure with the box alone against T/E_1.

### Removing the box-only option

It is self-contained. Delete:
- in `cv_plot.py`: the block between `BOX-ONLY (delete from here ... )` and `... to here`, the three short blocks marked `BOX-ONLY (delete these ... lines ...)`, and the `SHOW_BOX_ONLY` constant;
- in `cv_interactive.py`: the 3 lines marked `BOX-ONLY`, the `"empty box only"` check box entry, and `cv_quantum_box_only` from the import line.

## Validity

The closed forms drop terms of order e^{−αL}, which are negligible for αL ≳ 10. For αL ≥ 60 they agree with sums over the exact levels to about 2e-13 k_B. The classical integral stays accurate down to k_B T ≈ depth / 20000.
