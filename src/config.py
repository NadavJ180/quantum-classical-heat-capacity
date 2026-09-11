"""
config.py
=====================================================================
Single source of truth for the current run's physical setup: the
potential, its constants, and the numerical control parameters.

This is exactly "Section 0" of Quantum_HO_Master.py, pulled out into
its own module so that OTHER scripts (the figures/ tools, tests, a
notebook, ...) can import the potential / bounds / parameters WITHOUT
importing Quantum_HO_Master.py itself -- which would otherwise run
the full six-section pipeline (DVR solve, reference generation, Cv
sweep, ...) as a side effect of the import.

To switch systems: edit `my_potential`, `SYSTEM_NAME`, and
`T_UNITS_LABEL` below (swap which block is commented out), exactly as
before. Nothing else in the repository needs to change.

POTENTIAL_PARAMS is the single source of truth for this run's named
potential parameters: `my_potential` reads its coefficients from it
(so there is exactly one place to change a value -- no separate label
that can drift out of sync with the actual formula), and it is also
what figures/output_paths.py uses to name each run's figure folder
(figures/<system>/<params>/<category>/...). This is what makes a
future bulk scan -- e.g. over a range of the double well's `b` --
land every parameter combination in its own, unambiguous folder: each
iteration just sets POTENTIAL_PARAMS (or calls
figures.output_paths.set_context directly) before re-running the
pipeline.
=====================================================================
"""

# --- Physical constants (dimensionless: ℏ = m = ω = k_B = 1) ---
MASS, HBAR, OMEGA = 1.0, 1.0, 1.0

# --- Potential function (swap this for any smooth V(x)) ---
r'''
POTENTIAL_PARAMS = {"mass": MASS, "omega": OMEGA}

def my_potential(x, p=POTENTIAL_PARAMS):
    """1-D harmonic oscillator: V(x) = ½ m ω² x²."""
    return 0.5 * p["mass"] * (p["omega"]**2) * x**2

SYSTEM_NAME   = "1-D Harmonic Oscillator"
T_UNITS_LABEL = r"$k_B T \,/\, \hbar\omega$"
'''

# a, b, c, d: coefficients of V(x) = a*x^4 + b*x^3 + c*x^2 + d*x.
# b != 0 breaks the x -> -x symmetry (b = 0 recovers the symmetric well).
POTENTIAL_PARAMS = {"a": 0.25, "b": -0.5, "c": -0.5, "d": 0.0}

def my_potential(x, p=POTENTIAL_PARAMS):
        """1-D asymmetric double well: V(x) = a x^4 + b x^3 + c x^2 + d x."""
        return p["a"] * (x ** 4) + p["b"] * (x ** 3) + p["c"] * (x ** 2) + p["d"] * x

SYSTEM_NAME   = "1-D asymmetric double well"
T_UNITS_LABEL = r"$k_B T \,/\, \hbar\omega$"

# --- DVR base grid: number of energy levels ---
# The n-convergence diagnostic (Section 4) will confirm the exact
# number needed; 500 gives comfortable headroom for this system.
NUM_STATES = 500

# --- Temperature sweep ---
BETA_MIN, BETA_MAX, N_BETA = 0.01, 50.0, 1000

# --- xi / n convergence parameters ---
# XI_START = 3.0: first probe is already at effective T/9, allowing
# the classical-limit plateau to be found at much colder temperatures
# than XI_START = 1.0 would permit.
XI_START      = 3.0
TOL_XI        = 5e-3
MIN_STABLE_XI = 3
XI_MULT       = 1.1
MAX_XI_STEPS  = 80
TOL_CV        = 1e-4
MIN_STABLE_N  = 3

# --- DVR limit analysis tolerance (Section 5) ---
LIMIT_TOLERANCE = 1e-6

# --- Numerical reference scaling (Sections 2, 3, 5, 6) ---
# REFERENCE_SPAN_FACTOR: multiply base span by this (2.0 = double L)
# REFERENCE_DX_FACTOR:   divide base dx   by this (2.0 = halve Δx)
# Set INTERACTIVE = True to be prompted at runtime instead.
INTERACTIVE_REFERENCE_SCALING = False
REFERENCE_SPAN_FACTOR         = 2.0
REFERENCE_DX_FACTOR           = 2.0

# Human-readable label built from the scaling factors (used in plots).
ref_label = (f"numerical reference  "
             f"(span\u00d7{REFERENCE_SPAN_FACTOR:.2g}, "
             f"dx\u00f7{REFERENCE_DX_FACTOR:.2g})")