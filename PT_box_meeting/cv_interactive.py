"""
cv_interactive.py -- Cv_q and Cv_cl vs T with sliders.

    python cv_interactive.py

Sliders:
    alpha  - inverse width of the well (its depth is hbar^2 alpha^2 / m)
    L      - box length
    gamma  - scale-invariant stroke: alpha -> alpha/gamma, L -> gamma*L.
             Every level goes E_n -> E_n / gamma^2, so the curves only slide
             along T (by 1/gamma^2) and keep their shape.  The grey curves are
             the system before the stroke (gamma = 1).
Check boxes:
    log T axis
    empty box only - dotted: Cv_q - Cv_cl of the box without the well (wall effect)

The formulas are the ones in cv_plot.py (same folder).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons

from cv_plot import HBAR, M, KB, cv_quantum, cv_classical, cv_quantum_box_only

C_Q, C_CL, C_BOX, C_GHOST = "#2a78d6", "#eb6834", "#1baf7a", "0.7"

fig = plt.figure(figsize=(10, 7.5))
ax = fig.add_axes([0.08, 0.50, 0.64, 0.44])
ax_d = fig.add_axes([0.08, 0.27, 0.64, 0.19], sharex=ax)

s_alpha = Slider(fig.add_axes([0.16, 0.15, 0.50, 0.03]), r"$\alpha$", 0.5, 2.0, valinit=1.0)
s_logL = Slider(fig.add_axes([0.16, 0.10, 0.50, 0.03]), r"$L$", 1.0, 3.5, valinit=2.0)
s_gamma = Slider(fig.add_axes([0.16, 0.05, 0.50, 0.03]), r"$\gamma$ (stroke)", 0.5, 2.0, valinit=1.0)
checks = CheckButtons(fig.add_axes([0.76, 0.05, 0.22, 0.13]), ["log T axis", "empty box only"], [False, False])
info = fig.text(0.76, 0.93, "", fontsize=9, va="top", family="monospace")


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


def update(_=None):
    alpha, L, gamma = s_alpha.val, 10**s_logL.val, s_gamma.val
    s_logL.valtext.set_text("%.0f" % L)
    log_T = checks.get_status()[0]

    depth = HBAR**2 * alpha**2 / M                    # well depth before the stroke
    T = depth / KB * np.logspace(-3.5, 1, 800)
    alpha_s, L_s = alpha / gamma, gamma * L           # after the stroke

    cq = cv_quantum(T, alpha_s, L_s)
    ccl = cv_classical(T, alpha_s, L_s)
    d = cq - ccl

    ax_d.set_xscale("linear")
    ax.clear()
    ax_d.clear()
    if abs(gamma - 1) > 1e-3:
        ax.plot(T, cv_quantum(T, alpha, L) / KB, color=C_GHOST, lw=1, label=r"before stroke ($\gamma$ = 1)")
        ax.plot(T, cv_classical(T, alpha, L) / KB, "--", color=C_GHOST, lw=1)
    ax.plot(T, cq / KB, color=C_Q, label="quantum")
    ax.plot(T, ccl / KB, "--", color=C_CL, label="classical")
    ax.fill_between(T, ccl / KB, cq / KB, where=d > 0, color=C_Q, alpha=0.2, label=r"$C_v^q>C_v^{cl}$")

    ax_d.axhline(0, color="0.3", lw=0.8)
    ax_d.plot(T, d / KB, color="0.2")
    ax_d.fill_between(T, 0, d / KB, where=d > 0, color=C_Q, alpha=0.2)

    # ---- BOX-ONLY (delete these 3 lines to remove) ----
    if checks.get_status()[1]:
        ax_d.plot(T, cv_quantum_box_only(T, L_s) / KB - 0.5, ":", color=C_BOX, lw=1.5, label="empty box only")
        ax_d.legend(loc="lower right", fontsize=8)

    ax.set_title(r"Pöschl–Teller well ($\lambda=2$) in a box,  $\alpha L$ = %.0f" % (alpha * L))
    ax.set_ylabel(r"$C_v/k_B$")
    ax.legend(loc="upper right", fontsize=8)
    ax_d.set_ylabel(r"$(C_v^q-C_v^{cl})/k_B$")
    ax_d.set_xlabel(r"$k_BT$   ($\hbar$ = %g, m = %g, $k_B$ = %g)" % (HBAR, M, KB))
    if log_T:
        ax_d.set_xscale("log")
        ax_d.set_xlim(T[0], T[-1])
    else:
        ax_d.set_xscale("linear")
        ax_d.set_xlim(0, 0.6 * depth / KB)
    for a_ in (ax, ax_d):
        a_.grid(alpha=0.3)
    plt.setp(ax.get_xticklabels(), visible=False)

    j = np.argmax(d)
    eps = HBAR**2 * alpha_s**2 / (2 * M)              # hbar^2 alpha^2 / 2m after the stroke
    text = "alpha L = %.0f\n" % (alpha_s * L_s)
    if abs(gamma - 1) > 1e-3:
        text += "stroke: alpha L unchanged,\nevery T scaled by 1/gamma^2 = %.3f\n" % (1 / gamma**2)
    text += "\n"
    text += "max (Cv_q - Cv_cl)\n  = %+.3f k_B\n" % (d[j] / KB)
    text += "  at k_B T = %.4g\n" % T[j]
    text += "  = %.3f hbar^2 alpha^2/2m\n\n" % (KB * T[j] / eps)
    text += "Cv_q > Cv_cl for k_B T in:\n"
    for (t1, t2) in windows(T, d)[:3]:
        text += "  [%.3g, %.3g]\n" % (t1, t2)
    text += "(a window reaching the right\n edge is the wall effect,\n see 'empty box only')"
    info.set_text(text)
    fig.canvas.draw_idle()


s_alpha.on_changed(update)
s_logL.on_changed(update)
s_gamma.on_changed(update)
checks.on_clicked(update)
update()

if __name__ == "__main__":
    plt.show()
