"""
output_paths.py  (src/figures/)
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Resolves where a generated figure should be saved on disk, and keeps
the repo-root figures/ folder organized as runs accumulate:

    figures/<system>/<params>/<category>/<name>.png

<system>  -- slug of SYSTEM_NAME (e.g. "1_d_asymmetric_double_well").
<params>  -- LaTeX-safe slug of the potential's parameter values (e.g.
             "a-0p25_b-m0p5_c-m0p5_d-0" for {"a":0.25,"b":-0.5,
             "c":-0.5,"d":0.0}), so that two runs of the SAME potential
             with DIFFERENT parameters (the double well's `b`, say)
             land in separate folders instead of overwriting each
             other. This is what makes a future bulk scan over a range
             of parameter values (e.g. b = -1.0, -0.9, ..., 1.0) safe
             to run unattended: every (system, params) combination
             gets its own directory tree, keyed only by the values
             actually used, with no manual bookkeeping required.
<category>-- which kind of diagnostic the figure is (energy_levels,
             convergence, cv, cv_benchmark, dvr_limits, potential),
             matching the categories already named in FINDINGS.md's
             "Reading the Diagnostic Plots" table.

Every path component (system, params, category, and every figure
filename) is built only from [A-Za-z0-9_-] -- no spaces, dots, commas,
equals signs, or other characters that LaTeX's \includegraphics /
\graphicspath can be fussy about -- so the whole figures/ tree can be
pointed at directly from a .tex file (or copied wholesale into a
LaTeX figures/ folder) without renaming anything by hand. See _encode
below for exactly how a value is turned into a filename-safe token
("-" for a minus sign, "p" for a decimal point -- so 0.25 -> "0p25"
and -0.5 -> "m0p5").

USAGE
---------------------------------------------------------------------
A driver script calls `set_context(SYSTEM_NAME, POTENTIAL_PARAMS)`
once, near the start of a run (Quantum_HO_Master.py and
figures/plot_potential.py both do this from config.py's values).
Every plotting function then calls `save_figure(fig, category, name)`
once its figure is finished; it resolves the current run's directory
tree from that context, creates it if needed, writes the PNG, and (by
default) closes the figure -- no plt.show() anywhere in this pipeline,
so a run never blocks on a plot window and a long/bulk run never
accumulates open figures.

If no context has been set (e.g. a plotting function is called
directly, outside any driver script -- common when iterating in a
notebook), figures fall back to figures/_unspecified_run/<category>/
rather than raising.

NOTE ON BULK SCANS: this module holds its context as a plain module-
level dict, not per-process/thread-isolated state. A bulk-scan driver
that loops over parameter values in a single process should call
`set_context(...)` again at the start of each iteration, before that
iteration's pipeline runs -- the same pattern already used once per
script by Quantum_HO_Master.py and plot_potential.py, just repeated
in a loop.
=====================================================================
"""

import os
import re

import matplotlib.pyplot as plt

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FIGURES_ROOT = os.path.join(REPO_ROOT, "figures")

# Current run's (system_name, params) context, set via `set_context`.
_context = {"system_name": None, "params": None}


def set_context(system_name, params=None):
    """
    Declare which system/parameters the figures saved from now on
    (until the next call) belong to.

    Parameters
    ----------
    system_name : str
        E.g. config.SYSTEM_NAME.
    params : dict or None, optional
        E.g. config.POTENTIAL_PARAMS -- the named values that define
        this specific run of the potential (default None = no
        parameters to distinguish, all runs of this system share one
        "default_params" folder).
    """
    _context["system_name"] = system_name
    _context["params"] = dict(params) if params else {}


def _slugify(text):
    """Aggressive slug for the system-name folder: lowercase, alphanumerics
    and underscores only. Matches the convention already used for
    fig_potential_<slug>.png filenames elsewhere in this project."""
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", str(text).strip()).strip("_").lower()
    return slug or "unnamed_system"


# Whatever survives _encode() below is stripped down to this set as a final
# safety net, so a stray character in a future parameter name/value can
# never leak an unsafe character into a path LaTeX will be pointed at.
_NOT_LATEX_SAFE = re.compile(r"[^A-Za-z0-9_-]+")


def _encode(text):
    """
    Turn arbitrary text into a token built only from [A-Za-z0-9_-] --
    safe to use anywhere in a LaTeX \\includegraphics / \\graphicspath
    argument with no escaping. "-" (minus sign) becomes "m" and "."
    (decimal point) becomes "p" -- in that order, and as a blanket
    substitution rather than only at the start of the string, so this
    handles a minus sign anywhere it appears (e.g. in "1e-05") -- so a
    negative value never collides with its positive counterpart the
    way plain digit-stripping would (-0.5 -> "m0p5", 0.5 -> "0p5":
    still distinct).
    """
    encoded = str(text).replace("-", "m").replace(".", "p")
    return _NOT_LATEX_SAFE.sub("_", encoded)


def _fmt_value(v):
    """Format a single parameter value before encoding: floats use %g
    (drops trailing zeros: 0.2500 -> "0.25"), everything else str()."""
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def _params_dirname(params):
    if not params:
        return "default_params"
    parts = [f"{_encode(k)}-{_encode(_fmt_value(v))}" for k, v in params.items()]
    return "_".join(parts)


def figure_dir(category):
    """
    Return the directory this run's figures of `category` should be
    saved into, creating it (and its parents) if it doesn't exist yet.
    """
    system_name = _context["system_name"] or "_unspecified_run"
    d = os.path.join(
        FIGURES_ROOT,
        _slugify(system_name),
        _params_dirname(_context["params"]),
        _encode(category),
    )
    os.makedirs(d, exist_ok=True)
    return d


def save_figure(fig, category, name, dpi=220, close=True):
    """
    Save `fig` as <name>.png under this run's <category> folder
    (figures/<system>/<params>/<category>/<name>.png), overwriting any
    previous figure of the same name for this exact (system, params,
    category) combination -- a re-run with identical parameters is
    expected to refresh its own figures in place, while a run with
    different parameters lands in a different folder entirely.

    By default this also closes `fig` (via plt.close) after writing
    it, so a long pipeline run or bulk parameter scan never leaves an
    ever-growing pile of open figure objects around, doesn't need a
    display, and never blocks waiting for a plot window to be closed.
    Pass close=False if the caller still needs `fig` afterward.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save (the caller's `fig, ax = plt.subplots(...)`).
    category : str
        One of "energy_levels", "convergence", "cv", "cv_benchmark",
        "dvr_limits", "potential" (or any other short, consistent tag).
    name : str
        Filename without extension, e.g. "cv_summary".
    dpi : int, optional
        Output resolution (default 220, matching the rest of the
        project's figure-generating scripts).
    close : bool, optional
        Close `fig` after saving (default True).

    Returns
    -------
    path : str
        The full path the figure was written to.
    """
    path = os.path.join(figure_dir(category), f"{_encode(name)}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(fig)
    return path
