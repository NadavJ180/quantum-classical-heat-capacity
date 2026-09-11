"""
output_paths.py  (src/figures/)
=====================================================================
WHAT THIS FILE DOES
---------------------------------------------------------------------
Resolves where a generated figure should be saved on disk, and keeps
the repo-root figures/ folder organized as runs accumulate:

    figures/<system>/<params>/<category>/<name>.png

<system>  -- slug of SYSTEM_NAME (e.g. "1_d_asymmetric_double_well").
<params>  -- slug of the potential's parameter values (e.g.
             "a=0.25,b=-0.5,c=-0.5,d=0"), so that two runs of the SAME
             potential with DIFFERENT parameters (the double well's
             `b`, say) land in separate folders instead of overwriting
             each other. This is what makes a future bulk scan over a
             range of parameter values (e.g. b = -1.0, -0.9, ..., 1.0)
             safe to run unattended: every (system, params) combination
             gets its own directory tree, keyed only by the values
             actually used, with no manual bookkeeping required.
<category>-- which kind of diagnostic the figure is (energy_levels,
             convergence, cv, cv_benchmark, dvr_limits, potential),
             matching the categories already named in FINDINGS.md's
             "Reading the Diagnostic Plots" table.

USAGE
---------------------------------------------------------------------
A driver script calls `set_context(SYSTEM_NAME, POTENTIAL_PARAMS)`
once, near the start of a run (Quantum_HO_Master.py and
figures/plot_potential.py both do this from config.py's values).
Every plotting function then calls `save_figure(fig, category, name)`
right before `plt.show()`; it resolves the current run's directory
tree from that context, creates it if needed, and writes the PNG.

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


# Only characters that are actually unsafe/awkward in a path component
# are replaced here -- unlike _slugify, "=", ",", "-", and "." are kept
# AS IS, because stripping "-" would make "b=-0.5" and "b=0.5" collide
# into the same folder name, which is exactly wrong for a parameter scan.
_UNSAFE_PATH_CHARS = re.compile(r'[<>:"/\\|?*\s]+')


def _fmt_value(v):
    """Format a single parameter value for the folder name: floats use
    %g (drops trailing zeros: 0.2500 -> "0.25"), everything else str()."""
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def _params_dirname(params):
    if not params:
        return "default_params"
    parts = [f"{k}={_fmt_value(v)}" for k, v in params.items()]
    return _UNSAFE_PATH_CHARS.sub("_", ",".join(parts))


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
        category,
    )
    os.makedirs(d, exist_ok=True)
    return d


def save_figure(fig, category, name, dpi=220):
    """
    Save `fig` as <name>.png under this run's <category> folder
    (figures/<system>/<params>/<category>/<name>.png), overwriting any
    previous figure of the same name for this exact (system, params,
    category) combination -- a re-run with identical parameters is
    expected to refresh its own figures in place, while a run with
    different parameters lands in a different folder entirely.

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

    Returns
    -------
    path : str
        The full path the figure was written to.
    """
    path = os.path.join(figure_dir(category), f"{name}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
