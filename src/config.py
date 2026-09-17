"""
config.py
=====================================================================
Single source of truth for the current run's physical setup: the
potential, its constants, and the numerical control parameters.

This is exactly "Section 0" of Quantum_HO_Master.py, pulled out into
its own module so that OTHER scripts (the figures/ tools, tests, a
notebook, ...) can import the potential / bounds / parameters WITHOUT
importing Quantum_HO_Master.py itself -- which would otherwise run
the full seven-section pipeline (DVR solve, reference generation, Cv
sweep, coefficient sweep, ...) as a side effect of the import.

To switch systems: edit `my_potential`, `SYSTEM_NAME`, and
`T_UNITS_LABEL` below (swap which block is commented out). The only
other values that normally need adjusting for a new system are
`NUM_STATES` and `BETA_MAX` (a starting guess and the cold end of the
sweep, respectively) -- see the "Auto-tune escalation" and
"Temperature sweep" sections below for why the rest doesn't.

POTENTIAL_PARAMS is the single source of truth for this run's named
potential parameters: `my_potential` reads its coefficients from it
(so there is exactly one place to change a value -- no separate label
that can drift out of sync with the actual formula), and it is also
what figures/output_paths.py uses to name each run's figure folder
(figures/<system>/<params>/<category>/...). SCAN_PARAM/SCAN_STEP/
SCAN_COUNT below drive one built-in use of this: Cv_Coefficient_Sweep.py
sweeps a single named coefficient and plots every variant's quantum
Cv(T) in one shared figure, without giving each variant its own
figure folder. A bulk scan that instead needs each parameter
combination's own full set of diagnostics (all seven sections, in
its own folder) is a different, heavier pattern: loop over
POTENTIAL_PARAMS values, calling figures.output_paths.set_context(...)
again at the top of each iteration before re-running the pipeline.
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

# Mathtext shown on Cv_Coefficient_Sweep.py's plot annotation. Plain
# text (no <<param>> placeholders) is fine when the formula doesn't
# decompose into one additive term per coefficient, as here (mass and
# omega multiply inside a single term rather than each owning their
# own term) -- format_potential_formula() only substitutes placeholder
# tokens that are actually present, so a template with none is left
# as-is. Placeholders use <<name>> rather than Python's {name} so they
# never collide with LaTeX's own braces (e.g. \frac{1}{2} below).
POTENTIAL_FORMULA = r"$V(x) = \frac{1}{2} m \omega^2 x^2$"
'''

# a, b, c, d: coefficients of V(x) = a*x^4 + b*x^3 + c*x^2 + d*x.
# b != 0 breaks the x -> -x symmetry (b = 0 recovers the symmetric well).
POTENTIAL_PARAMS = {"a": 0.25, "b": -0.5, "c": -0.5, "d": 0.0}

def my_potential(x, p=POTENTIAL_PARAMS):
        return p["a"] * (x ** 4) + p["b"] * (x ** 3) + p["c"] * (x ** 2) + p["d"] * x

SYSTEM_NAME   = "1-D asymmetric double well"
T_UNITS_LABEL = r"$k_B T \,/\, \hbar\omega$"

# Mathtext template for the formula shown on Cv_Coefficient_Sweep.py's
# comparison plot (Section 7): a template with one <<name>> placeholder
# per POTENTIAL_PARAMS key, one term per coefficient with NO operators
# between placeholders -- each substituted value supplies its own
# leading sign/spacing (see Cv_Coefficient_Sweep.format_potential_formula),
# so terms can be reordered/added freely without touching that function.
# <<name>> (not Python's {name}) so placeholders never collide with
# LaTeX's own braces. Must be kept in sync with my_potential by hand,
# same as SYSTEM_NAME/T_UNITS_LABEL above.
POTENTIAL_FORMULA = r"$V(x) = <<a>>x^4<<b>>x^3<<c>>x^2<<d>>x$"

# --- DVR base grid: number of energy levels ---
# Starting guess only -- Quantum_HO_Master.py's auto-tune loop (see
# Cv_AutoTune.py) will grow this if the hot end of the sweep shows
# signs of thermal-occupation truncation (the "numerical Schottky
# anomaly"). 500 gives comfortable headroom for this system already.
NUM_STATES = 500

# --- Temperature sweep ---
# BETA_MIN = None means "auto": Quantum_HO_Master.py derives it from
# T_max = 10 * (E1 - E0) / k_B, the practical "T -> infinity" checkpoint
# by which the classical limit should already be reached. Set it to a
# float instead to hand-pick the hot end (e.g. to zoom into a specific
# temperature window) -- BETA_MAX always stays a plain, user-set float,
# since that is the knob actually meant to be changed run to run.
BETA_MIN, BETA_MAX, N_BETA = None, 50.0, 1000

# --- Auto-tune escalation (Cv_AutoTune.py) ---
# Starting guesses for NUM_STATES/XI_START/MAX_XI_STEPS above are
# escalated automatically if the sweep shows a truncation artifact at
# the hot end (not enough levels) or a ran-out-of-budget xi-scan at the
# cold end (not enough xi steps) -- see FINDINGS.md's Troubleshooting
# table for the physical rationale. These knobs control that loop and
# are meant to be touched rarely, if ever.
AUTO_ESCALATE               = True
MAX_ESCALATION_ROUNDS       = 4
NUM_STATES_GROWTH           = 1.6
NUM_STATES_CAP              = 4000
XI_START_GROWTH             = 2.0
MAX_XI_STEPS_GROWTH         = 1.5
HOT_STATE_SAFETY            = 20.0   # target E_max / (k_B * T_hot)
ESCALATION_FRACTION_THRESHOLD = 0.05  # fraction of a sweep half allowed to fail before escalating

# --- Coefficient sweep (Cv_Coefficient_Sweep.py) ---
# Which POTENTIAL_PARAMS key to vary, the spacing between consecutive
# values, and how many EXTRA potentials to add on each side of the
# value already set above (so the base value/potential is always one
# of the plotted curves). Total curves plotted = 2*SCAN_COUNT + 1.
SCAN_PARAM = "b"
SCAN_STEP  = 0.2
SCAN_COUNT = 2

# The value of SCAN_PARAM that recovers the SYMMETRIC potential (b = 0
# for this well -- see my_potential's comment above). None of the
# variants produced by SCAN_STEP/SCAN_COUNT around the current base
# value (b = -0.5) land on b = 0, so without this the coefficient
# sweep would never actually include a genuinely symmetric reference
# curve. run_coefficient_sweep adds it as an extra variant (flagged
# "(symmetric)" in the Cv plot legend) whenever it isn't already
# covered by the regular sweep. Set to None to disable.
SCAN_SYMMETRIC_VALUE = 0.0

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