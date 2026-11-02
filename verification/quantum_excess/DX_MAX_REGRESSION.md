# DVR_DX_MAX regression check

Written by `dx_max_regression_check.py`. Baseline: the code before the change (`dx_max_regression_baseline.npz`). Bit-identical = every element equal (NaN = NaN).

| run | quantity | bit-identical to baseline |
|---|---|---|
| DVR_DX_MAX = None | energies | yes |
| DVR_DX_MAX = None | grid | yes |
| DVR_DX_MAX = None | quick_cv_classical | yes |
| DVR_DX_MAX = None | quick_cv_quantum | yes |
| DVR_DX_MAX = None | quick_error_estimate | yes |
| DVR_DX_MAX = None | sweep_cv_classical | yes |
| DVR_DX_MAX = None | sweep_error_estimate | yes |
| DVR_DX_MAX = None | sweep_ref15_cv_classical | yes |
| DVR_DX_MAX = None | sweep_ref15_error_estimate | yes |
| DVR_DX_MAX = None | sweep_ref15_xi_conv | yes |
| DVR_DX_MAX = None | sweep_xi_conv | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | energies | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | grid | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | quick_cv_classical | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | quick_cv_quantum | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | quick_error_estimate | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | sweep_cv_classical | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | sweep_error_estimate | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | sweep_ref15_cv_classical | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | sweep_ref15_error_estimate | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | sweep_ref15_xi_conv | yes |
| DVR_DX_MAX = 10 (looser than the momentum rule) | sweep_xi_conv | yes |

All identical: yes. Run time 192 s.
