"""
config.py
Every knob of the solver, in one place. Nothing is changed automatically:
the scans only report where the results can be trusted (CODE_GUIDE.md).

    python Quick_Scan.py          quick look, a few minutes at most
    python Quantum_HO_Master.py   full diagnostics, can take hours

Units: hbar = m = k_B = 1 (change MASS / HBAR if you need other units).
"""

import potentials

# =====================================================================
# 1. The physical system
# =====================================================================
MASS, HBAR = 1.0, 1.0

POTENTIAL        = potentials.double_well
POTENTIAL_PARAMS = {"a": 0.25, "b": -0.5, "c": -0.5, "d": 0.0}
SYSTEM_NAME      = "1-D asymmetric double well"

# Other systems (copy the lines over the ones above, and adapt the grid below):
#   POTENTIAL = potentials.harmonic;       POTENTIAL_PARAMS = {"k": 1.0}
#   POTENTIAL = potentials.flat;           POTENTIAL_PARAMS = {}                         (WALLS = True)
#   POTENTIAL = potentials.poschl_teller;  POTENTIAL_PARAMS = {"depth": 1.0, "width": 1.0} (WALLS = True)
#   POTENTIAL = potentials.dimple_in_trap; POTENTIAL_PARAMS = {"D": 37.5, "s": 0.1, "L": 10.0}


def my_potential(x, p=None):
    """V(x) of the system above; p replaces POTENTIAL_PARAMS (used by the coefficient sweep)."""
    if p is None:
        p = POTENTIAL_PARAMS
    return POTENTIAL(x, **p)


# =====================================================================
# 2. Grid: where the particle can be
# =====================================================================
# WALLS = False: V(x) confines the particle. V must be large at X_MIN and X_MAX
#                (the grid ends act as walls, so they must be where the particle never goes).
# WALLS = True:  infinite walls at X_MIN and X_MAX (particle in a box, plus V).
X_MIN, X_MAX = -6.5, 7.5
WALLS        = False

# =====================================================================
# 3. Temperatures (k_B T, log-spaced)
# =====================================================================
T_MIN, T_MAX  = 0.05, 10.0
T_UNITS_LABEL = r"$k_B T$"

# =====================================================================
# 4. Quick scan (Quick_Scan.py) and full pipeline (Quantum_HO_Master.py)
# =====================================================================
# dx         : grid spacing for xi = 1
# num_states : number of energy levels kept for xi = 1
# xi_values  : hbar -> hbar/xi for the classical limit (ascending; the largest is the
#              classical estimate). Each xi uses dx/xi and num_states*xi levels, so the
#              matrix has (number of grid points * xi)^2 entries: memory grows as xi^2.
# n_T        : number of temperatures between T_MIN and T_MAX
# cv_tol     : trust tolerance on Cv (section 5 below)
QUICK_SCAN = {"dx": 0.05,  "num_states": 80,  "xi_values": [1, 2, 4, 8, 16], "n_T": 60,  "cv_tol": 1e-2}
FULL_RUN   = {"dx": 0.04,  "num_states": 120, "xi_values": [1, 2, 4, 8, 16, 32], "n_T": 400, "cv_tol": 1e-3}

QUICK_SCAN_MODE = "single"      # "single": the system above; "sweep": every value in SCAN_VALUES

# =====================================================================
# 5. Trust checks (reported, never acted on)
# =====================================================================
# The quantum Cv is recomputed with ONE knob changed at a time (dx / REFERENCE_DX_FACTOR,
# span * REFERENCE_SPAN_FACTOR, num_states * REFERENCE_LEVELS_FACTOR). It is trusted where
# none of these changes moves Cv by more than cv_tol. The classical Cv is trusted where its
# estimated distance to the hbar -> 0 limit is below cv_tol * Cv.
REFERENCE_DX_FACTOR     = 2.0
REFERENCE_SPAN_FACTOR   = 2.0
REFERENCE_LEVELS_FACTOR = 2.0
# Cv_q - Cv_cl counts as resolved only if it is larger than ERROR_MARGIN times its uncertainty.
ERROR_MARGIN            = 2.0
# Grid suggestions: levels more than THERMAL_CUTOFF * k_B T above V_min are treated as empty.
THERMAL_CUTOFF          = 20.0

# Full pipeline, Section 5 (DVR limit analysis): the grid is redone with each of these
# factors (dx * factor, span * factor) and every level is compared with the reference grid.
# A level is "correct" if it is within LIMIT_TOLERANCE (energy units) of the reference.
LIMIT_TOLERANCE    = 1e-6
LIMIT_DX_FACTORS   = [0.75, 1.0, 1.5, 2.0, 3.0]
LIMIT_SPAN_FACTORS = [0.7, 0.85, 1.0, 1.25]

# Full pipeline, Section 6: also recompute the classical limit on the reference grid.
# Expensive: every xi solve gets REFERENCE_SPAN_FACTOR * REFERENCE_DX_FACTOR times more points.
RUN_CLASSICAL_REFERENCE = False

# =====================================================================
# 6. Coefficient sweep (Quick_Scan.py in "sweep" mode, full pipeline Section 7)
# =====================================================================
# One POTENTIAL_PARAMS entry takes every value in SCAN_VALUES; everything else stays fixed.
# Set SCAN_VALUES = [] to skip Section 7 of the full pipeline.
SCAN_PARAM  = "b"
SCAN_VALUES = [-0.9, -0.5, -0.2, 0.0]
