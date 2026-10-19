"""
output_paths.py  (src/figures/)
Where figures are saved:  figures/<system>/<params>/<category>/<name>.png
<params> is built from the potential's parameters (e.g. "a-0p25_b-m0p5"), so runs
with different parameters never overwrite each other. Names use only
[A-Za-z0-9_-] ("-" -> "m", "." -> "p") so LaTeX can include them directly.
"""

import os
import re

import matplotlib.pyplot as plt

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FIGURES_ROOT = os.path.join(REPO_ROOT, "figures")

_context = {"system_name": None, "params": None}


def set_context(system_name, params=None):
    """Declare which system / parameters the next figures belong to."""
    _context["system_name"] = system_name
    _context["params"] = dict(params) if params else {}


def _slugify(text):
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", str(text).strip()).strip("_").lower()
    return slug or "unnamed_system"


def _encode(text):
    encoded = str(text).replace("-", "m").replace(".", "p")
    return re.sub(r"[^A-Za-z0-9_-]+", "_", encoded)


def _params_dirname(params):
    if not params:
        return "default_params"
    parts = []
    for key, value in params.items():
        text = f"{value:g}" if isinstance(value, float) else str(value)
        parts.append(f"{_encode(key)}-{_encode(text)}")
    return "_".join(parts)


def figure_dir(category):
    """Folder for this run's figures of one category (created if needed)."""
    folder = os.path.join(FIGURES_ROOT, _slugify(_context["system_name"] or "_unspecified_run"),
                          _params_dirname(_context["params"]), _encode(category))
    os.makedirs(folder, exist_ok=True)
    return folder


def save_figure(fig, category, name, dpi=220, close=True):
    """Save fig as <name>.png in this run's <category> folder and return the path."""
    path = os.path.join(figure_dir(category), f"{_encode(name)}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(fig)
    return path
