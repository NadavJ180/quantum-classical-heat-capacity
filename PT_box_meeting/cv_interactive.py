"""
cv_interactive.py -- Cv_q and Cv_cl vs T with sliders, next to the potential.

    python cv_interactive.py

Sliders:
    alpha  - inverse width of the well (its depth is hbar^2 alpha^2 / m)
    L      - box length
    gamma  - scale-invariant stroke: alpha -> alpha/gamma AND L -> gamma*L,
             so alpha*L (the shape) is kept and every level goes E_n -> E_n/gamma^2.
             The grey curves are the system before the stroke (gamma = 1).
    The box to the right of each slider shows its value. Click it, type an exact
    value and press Enter. A value outside the slider's range widens the slider.
Reset to defaults button: alpha = 1, L = 100, gamma = 1, all check boxes off.
Check boxes:
    log T axis
    plot vs gamma^2 T - undoes the stroke on the axes (T -> gamma^2 T, and x -> x/gamma,
                        V -> gamma^2 V in the potential plots). The coloured curves then
                        lie exactly on the grey ones: the scaling is homogeneous.
    empty box only    - the box without the well (green), in both Cv panels

The formulas are the ones in cv_plot.py (same folder).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, TextBox, Button

from cv_plot import HBAR, M, KB, cv_quantum, cv_classical, cv_quantum_box_only

C_Q, C_CL, C_BOX, C_GHOST = "#2a78d6", "#eb6834", "#1baf7a", "0.6"
LOG_T, RESCALE, BOX_ONLY = "log T axis", r"plot vs $\gamma^2 T$", "empty box only"
LABELS = [LOG_T, RESCALE, BOX_ONLY]

fig = plt.figure(figsize=(12.5, 8))
ax = fig.add_axes([0.06, 0.53, 0.50, 0.41])
ax_d = fig.add_axes([0.06, 0.30, 0.50, 0.19], sharex=ax)
ax_box = fig.add_axes([0.65, 0.72, 0.33, 0.22])
ax_well = fig.add_axes([0.65, 0.42, 0.33, 0.22])

s_alpha = Slider(fig.add_axes([0.12, 0.17, 0.38, 0.025]), r"$\alpha$", 0.5, 2.0, valinit=1.0)
s_logL = Slider(fig.add_axes([0.12, 0.12, 0.38, 0.025]), r"$L$", 1.0, 3.5, valinit=2.0)     # slider position = log10(L)
s_gamma = Slider(fig.add_axes([0.12, 0.07, 0.38, 0.025]), r"$\gamma$ (stroke)", 0.5, 2.0, valinit=1.0)
box_alpha = TextBox(fig.add_axes([0.52, 0.165, 0.06, 0.035]), "")
box_L = TextBox(fig.add_axes([0.52, 0.115, 0.06, 0.035]), "")
box_gamma = TextBox(fig.add_axes([0.52, 0.065, 0.06, 0.035]), "")
fig.text(0.52, 0.21, "type + Enter", fontsize=8)
reset_button = Button(fig.add_axes([0.12, 0.012, 0.14, 0.04]), "Reset to defaults")
checks = CheckButtons(fig.add_axes([0.83, 0.05, 0.16, 0.22]), LABELS, [False] * len(LABELS))
info = fig.text(0.60, 0.31, "", fontsize=8, va="top", family="monospace")

SLIDERS = [s_alpha, s_logL, s_gamma]
RANGES = [(s.valmin, s.valmax) for s in SLIDERS]      # restored by the reset button
for s in SLIDERS:
    s.valtext.set_visible(False)                      # the text boxes show the values instead
state = {"busy": False}                               # True while the code itself moves the widgets


def windows(T, d):
    """Temperature intervals where d > 0."""
    out = []
    start = None
    for i in range(len(T)):
        if d[i] > 0 and start is None:
            start = T[i]
        if d[i] <= 0 and start is not None:
            out.append((start, T[i - 1]))
            start = None
    if start is not None:
        out.append((start, T[-1]))
    return out


def potential(x, alpha, L):
    """V(x) inside the box, nan outside (the walls are drawn separately)."""
    V = -HBAR**2 * alpha**2 / M / np.cosh(np.minimum(np.abs(alpha * x), 300))**2   # cosh(300)^2 is already ~1e260
    V[np.abs(x) > L / 2] = np.nan
    return V


def plot_potential(axp, u, alpha, L, sx, sV, style, walls_top):
    """Draw V and its walls at display coordinates u = sx * x, with V shown as sV * V."""
    axp.plot(u, sV * potential(u / sx, alpha, L), **style)
    for side in (-1, 1):
        axp.plot([side * sx * L / 2] * 2, [0, walls_top], **style)          # hard walls, V = +infinity
    Eb = -HBAR**2 * alpha**2 / (2 * M)                   # the bound state
    xb = np.arccosh(np.sqrt(2.0)) / alpha                # where V = E_b
    bar = dict(style)
    bar["lw"] = style["lw"] + 1.5
    axp.plot([-sx * xb, sx * xb], [sV * Eb, sV * Eb], **bar)


def show_values():
    """Write the slider values into the text boxes (state['busy'] stops this from counting as typing)."""
    state["busy"] = True
    box_alpha.set_val("%.6g" % s_alpha.val)
    box_L.set_val("%.6g" % 10**s_logL.val)
    box_gamma.set_val("%.6g" % s_gamma.val)
    state["busy"] = False


def typed(slider, text, log=False):
    """Move a slider to a typed value. A value outside the slider's range widens the slider."""
    if state["busy"]:
        return
    try:
        value = float(text)
    except ValueError:
        value = -1.0
    if value <= 0:                                    # not a positive number: put the old value back
        show_values()
        return
    if log:
        value = np.log10(value)
    if value < slider.valmin or value > slider.valmax:
        slider.valmin = min(slider.valmin, value)
        slider.valmax = max(slider.valmax, value)
        slider.ax.set_xlim(slider.valmin, slider.valmax)
    slider.set_val(value)                             # this calls update()


def reset(_=None):
    """Back to the default values, slider ranges and check boxes."""
    state["busy"] = True
    for s, (lo, hi) in zip(SLIDERS, RANGES):
        s.valmin, s.valmax = lo, hi
        s.ax.set_xlim(lo, hi)
        s.reset()
    for i, on in enumerate(checks.get_status()):
        if on:
            checks.set_active(i)                      # switches it off
    state["busy"] = False
    update()


def update(_=None):
    if state["busy"]:
        return
    alpha, L, gamma = s_alpha.val, 10**s_logL.val, s_gamma.val
    opt = dict(zip(LABELS, checks.get_status()))
    stroke = abs(gamma - 1) > 1e-3

    depth = HBAR**2 * alpha**2 / M                    # well depth before the stroke
    T = depth / KB * np.logspace(-3.5, 1, 800)
    alpha_s, L_s = alpha / gamma, gamma * L           # after the stroke

    cq = cv_quantum(T, alpha_s, L_s)
    ccl = cv_classical(T, alpha_s, L_s)
    d = cq - ccl
    xT = gamma**2 * T if opt[RESCALE] else T          # where the after-stroke curves are drawn

    ax_d.set_xscale("linear")
    ax.clear()
    ax_d.clear()
    if stroke:
        ax.plot(T, cv_quantum(T, alpha, L) / KB, color=C_GHOST, lw=4, alpha=0.4, label=r"before stroke ($\gamma$ = 1)")
        ax.plot(T, cv_classical(T, alpha, L) / KB, "--", color=C_GHOST, lw=4, alpha=0.4)
    ax.plot(xT, cq / KB, color=C_Q, label="quantum")
    ax.plot(xT, ccl / KB, "--", color=C_CL, label="classical")
    ax.fill_between(xT, ccl / KB, cq / KB, where=d > 0, color=C_Q, alpha=0.2, label=r"$C_v^q>C_v^{cl}$")

    ax_d.axhline(0, color="0.3", lw=0.8)
    ax_d.plot(xT, d / KB, color="0.2")
    ax_d.fill_between(xT, 0, d / KB, where=d > 0, color=C_Q, alpha=0.2)

    # ---- BOX-ONLY (delete this block to remove) ----
    if opt[BOX_ONLY]:
        cq_box = cv_quantum_box_only(T, L_s)
        ax.plot(xT, cq_box / KB, ":", color=C_BOX, lw=1.8, label="empty box: quantum")
        ax.plot(xT, 0.5 + 0 * T, "-.", color=C_BOX, lw=1, label="empty box: classical")
        ax_d.plot(xT, cq_box / KB - 0.5, ":", color=C_BOX, lw=1.8, label="empty box only")
        ax_d.legend(loc="lower right", fontsize=8)
    # ---- end BOX-ONLY ----

    ax.set_title(r"Pöschl–Teller well ($\lambda=2$) in a box,  $\alpha L$ = %.0f" % (alpha * L))
    ax.set_ylabel(r"$C_v/k_B$")
    ax.legend(loc="upper right", fontsize=8)
    ax_d.set_ylabel(r"$(C_v^q-C_v^{cl})/k_B$")
    xlabel = r"$\gamma^2 k_BT$" if opt[RESCALE] else r"$k_BT$"
    ax_d.set_xlabel(xlabel + r"   ($\hbar$ = %g, m = %g, $k_B$ = %g)" % (HBAR, M, KB))
    if opt[LOG_T]:
        ax_d.set_xscale("log")
        ax_d.set_xlim(T[0], T[-1])
    else:
        ax_d.set_xlim(0, 0.6 * depth / KB)
    for a_ in (ax, ax_d):
        a_.grid(alpha=0.3)
    plt.setp(ax.get_xticklabels(), visible=False)

    # ---- the potential: whole box (top) and a zoom on the well (bottom) ----
    sx, sV = (1 / gamma, gamma**2) if opt[RESCALE] else (1.0, 1.0)
    D = max(depth, sV * HBAR**2 * alpha_s**2 / M)
    top = 0.35 * D
    X = 0.55 * max(L, sx * L_s)
    w = max(1 / alpha, sx / alpha_s)                    # well width on the plot
    u = np.sort(np.concatenate([np.linspace(-X, X, 3000), np.linspace(-8 * w, 8 * w, 1500)]))
    for axp, title in [(ax_box, "potential: whole box"), (ax_well, "potential: zoom on the well")]:
        axp.clear()
        if stroke:
            plot_potential(axp, u, alpha, L, 1.0, 1.0, {"color": C_GHOST, "lw": 4, "alpha": 0.4}, top)
        plot_potential(axp, u, alpha_s, L_s, sx, sV, {"color": "k", "lw": 1.3}, top)
        axp.set_ylim(-1.15 * D, top)
        axp.set_title(title, fontsize=9)
        axp.grid(alpha=0.3)
        axp.tick_params(labelsize=8)
    ax_box.set_xlim(-X, X)
    ax_well.set_xlim(-6 * w, 6 * w)
    ax_well.set_xlabel(r"$x/\gamma$" if opt[RESCALE] else r"$x$", fontsize=9)
    ylabel = r"$\gamma^2 V$" if opt[RESCALE] else r"$V$"
    ax_box.set_ylabel(ylabel, fontsize=9)
    ax_well.set_ylabel(ylabel + "  (bar: bound state)", fontsize=9)

    # ---- numbers ----
    j = np.argmax(d)
    eps = HBAR**2 * alpha_s**2 / (2 * M)              # hbar^2 alpha^2 / 2m after the stroke
    text = "alpha L = %.0f\n" % (alpha_s * L_s)
    if stroke:
        text += "stroke: alpha -> alpha/%.2f, L -> %.2f L\n" % (gamma, gamma)
        text += "  every E_n -> E_n/gamma^2 = %.3f E_n\n" % (1 / gamma**2)
    text += "\nmax (Cv_q - Cv_cl) = %+.3f k_B\n" % (d[j] / KB)
    text += "  at k_B T = %.4g\n" % T[j]
    text += "     = %.3f hbar^2 alpha^2/2m\n" % (KB * T[j] / eps)
    text += "Cv_q > Cv_cl for k_B T in:\n"
    for (t1, t2) in windows(T, d)[:3]:
        edge = "  (walls)" if t2 == T[-1] else ""
        text += "  [%.3g, %.3g]%s\n" % (t1, t2, edge)
    info.set_text(text)
    show_values()
    fig.canvas.draw_idle()


s_alpha.on_changed(update)
s_logL.on_changed(update)
s_gamma.on_changed(update)
checks.on_clicked(update)
box_alpha.on_submit(lambda text: typed(s_alpha, text))
box_L.on_submit(lambda text: typed(s_logL, text, log=True))
box_gamma.on_submit(lambda text: typed(s_gamma, text))
reset_button.on_clicked(reset)
update()

if __name__ == "__main__":
    plt.show()
