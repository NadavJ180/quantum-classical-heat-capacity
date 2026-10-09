# Quick-scan validation

Does [`src/Quick_Scan.py`](../../src/Quick_Scan.py) give the same answer as the full pipeline, and is its error estimate honest? Run on the double well of `config.py` (b = −0.9 … 0, as Section 7).

| file | what it is |
|---|---|
| [`validate_quick_scan.py`](validate_quick_scan.py) | Runs the quick scan at resolutions 1 and 2 (sweep) and 3 (b = −0.5), and compares each run with the full pipeline and with the exact classical $C_v$ (~20 min). Pins its own potential and window, so local edits of `config.py` don't affect it. |
| [`QUICK_SCAN_VALIDATION.txt`](QUICK_SCAN_VALIDATION.txt) | The results: time, DVR solves, quantum and classical errors, true error ÷ estimate, verdicts. |
| [`full_resolution_section7.npz`](full_resolution_section7.npz) | Per-variant quantum and classical curves of the full-resolution Section 7 run (config.py settings, 1000 temperatures, 2026-10-09; summarized in [`../classical_limit/SECTION7_RESULTS.txt`](../classical_limit/SECTION7_RESULTS.txt)). Stored so the comparison doesn't need a ~44-minute rerun. |
| `figures/.../quick_scan/` | The three quick-scan figures of that run. |

**Result** (with `RememberedPotential`; the first run, before it, took 5.0, 10.6 and 9.9 min with bit-identical numbers).

| resolution | run | time | quantum vs full pipeline | classical error vs exact | true error ÷ estimate | wrong resolved verdicts |
|---|---|---|---|---|---|---|
| 1 | six-variant sweep | 2.0 min | ≤ 1.5e-7 | ≤ 5.2e-3 | 1.03–1.10 | 0 |
| 2 | six-variant sweep | 4.6 min | ≤ 1.8e-9 | ≤ 3.0e-3 | 1.01–1.04 | 0 |
| 3 | b = −0.5 | 6.1 min | 2.2e-12 | 8.4e-4 | 1.01 | 0 |

Every b came out "below" (no quantum excess), the same conclusion as the full run. The exact classical $C_v$ is used here only as an independent check; the quick scan itself is fully numerical (the ξ-scan).
