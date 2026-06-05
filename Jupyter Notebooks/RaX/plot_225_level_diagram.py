"""Render the 225RaF A(0)-X(0) hyperfine cycling + branching level diagram.

Reads frozen numbers from data_225_level_diagram.py (no Structure env needed).
Energy NOT to scale (schematic); intra-manifold F1 ordering honored.
F1-grouped (branching summed over F). Paper/talk sizing.

    python plot_225_level_diagram.py
Outputs figures/225RaF_XA_hf_cycling_v1.{pdf,svg,png}
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

import data_225_level_diagram as D

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, 'figures')

# ----------------------------------------------------------------- STYLE knobs
STYLE = dict(
    figsize=(9.6, 6.6),
    bar_w=1.05,                 # half-width of a level bar
    bar_lw=2.6,
    c_cycle='#1f4e8c',          # cycling levels (X F1=1, A F1=0)
    c_dark='#9aa0a6',           # dark / spectator levels
    c_level='#202124',          # ordinary levels
    c_pump='#c0392b',           # 752.8 nm pump
    c_decay='#37474f',          # spontaneous decay wiggles
    c_leak='#8e44ad',           # G=0 / N=3 leaks
    label_fs=12,
    f_fs=11,
    state_fs=13,
)

# ----------------------------------------------------------------- geometry (data coords)
# schematic two-band layout (energy NOT to scale; no break line -- the 752.8 nm
# pump label conveys the electronic gap). F1 ordering within each band is honored.
Y_X_F1_1 = 0.00            # X N=1 G=1 F1=1 (cycle ground)
Y_X_F1_0 = 0.70            # X N=1 G=1 F1=0  (~317 MHz above F1=1)
Y_X_F1_2 = 1.10            # X N=1 G=1 F1=2  (~nearly degenerate with F1=0; raised for clarity)
Y_X_G0   = 2.35            # X N=1 G=0 (~16.8 GHz above G=1; schematic leak target)
Y_A_F1_1 = 4.30            # A 2Pi1/2 J'=1/2 F1=1 (context)
Y_A_F1_0 = 5.05            # A 2Pi1/2 J'=1/2 F1=0 (cycling upper; ~2.46 GHz above F1=1)

X_CYC    = 1.85            # cycle column (X F1=1, A F1=0 directly above)
X_F0     = 3.45            # X F1=0
X_F2     = 5.50            # X F1=2
X_G0     = 7.25            # X G=0 (leak)
X_A_CTX  = 0.55            # A F1=1 context bar (left, short)


def _fmt_f(v):
    if v >= 0.1:
        return f'{v:.3f}'
    if v >= 1e-9:
        m, e = f'{v:.1e}'.split('e')
        return rf'{m}\times10^{{{int(e)}}}'
    return r'<10^{-9}'


def level_bar(ax, xc, y, color, lw=None, label=None, lab_dx=0.0, lab_dy=0.0,
              ha='left', va='center', label_color=None, fs=None, w=None):
    w = STYLE['bar_w'] if w is None else w
    ax.plot([xc - w, xc + w], [y, y], color=color, lw=lw or STYLE['bar_lw'],
            solid_capstyle='round', zorder=3)
    if label:
        ax.text(xc + w + lab_dx, y + lab_dy, label, color=label_color or color,
                ha=ha, va=va, fontsize=fs or STYLE['label_fs'], zorder=4)


def wavy(ax, x0, y0, x1, y1, color, lw, n=9, amp=0.05, alpha=1.0, ls='-'):
    t = np.linspace(0, 1, 400)
    x = x0 + (x1 - x0) * t
    y = y0 + (y1 - y0) * t
    dx, dy = x1 - x0, y1 - y0
    L = np.hypot(dx, dy)
    px, py = -dy / L, dx / L           # unit perpendicular
    env = np.sin(np.pi * t)            # taper to 0 at both ends
    w = amp * env * np.sin(2 * np.pi * n * t)
    ax.plot(x + px * w, y + py * w, color=color, lw=lw, alpha=alpha, ls=ls, zorder=2)
    # arrowhead at the end
    ax.add_patch(FancyArrowPatch((x[-3], y[-3]), (x1, y1), arrowstyle='-|>',
                 mutation_scale=12, color=color, alpha=alpha, lw=0, zorder=2))


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=STYLE['figsize'], constrained_layout=True)

    # ---- A 2Pi1/2 J'=1/2 levels (top band) ----
    level_bar(ax, X_A_CTX, Y_A_F1_1, STYLE['c_dark'], w=0.5,
              label=r'$F_1{=}1$', lab_dx=0.12, label_color=STYLE['c_dark'], fs=11)
    level_bar(ax, X_CYC, Y_A_F1_0, STYLE['c_cycle'],
              label=r"$|J'{=}1/2,\,F_1{=}0\rangle$", lab_dx=0.18)

    # ---- X 2Sigma N=1 G=1 levels (bottom band) ----
    level_bar(ax, X_CYC, Y_X_F1_1, STYLE['c_cycle'],
              label=r'$|G{=}1,\,F_1{=}1\rangle$', lab_dx=0.18)
    level_bar(ax, X_F0, Y_X_F1_0, STYLE['c_level'], w=0.8,
              label=r'$F_1{=}0$', lab_dx=0.15)
    level_bar(ax, X_F2, Y_X_F1_2, STYLE['c_level'], w=0.8,
              label=r'$F_1{=}2$', lab_dx=0.15)
    ax.text(X_F2 + 0.95, Y_X_F1_2 - 0.34, r'($F{=}5/2$ dark)', color=STYLE['c_dark'],
            ha='left', va='top', fontsize=9)

    # ---- X G=0 leak manifold (compact, upper right) ----
    level_bar(ax, X_G0, Y_X_G0, STYLE['c_leak'], w=0.75,
              label=r'$|G{=}0\rangle$' '\n' r'$+16.8$ GHz', lab_dx=0.15, fs=10)

    # ---- pump (752.8 nm), left of the cycle column ----
    ax.add_patch(FancyArrowPatch((X_CYC - 0.32, Y_X_F1_1), (X_CYC - 0.32, Y_A_F1_0),
                 arrowstyle='-|>', mutation_scale=16, color=STYLE['c_pump'],
                 lw=2.4, zorder=2))
    ax.text(X_CYC - 0.62, (Y_X_F1_1 + Y_A_F1_0) / 2,
            '752.8 nm\ncooling', color=STYLE['c_pump'], ha='center', va='center',
            fontsize=11, rotation=90)

    # ---- spontaneous decays from A(F1=0); label sits clear of the stems ----
    # (start_x, target, branching, color, lw, ls, label_pos)
    decays = [
        (X_CYC + 0.20, (X_CYC + 0.32, Y_X_F1_1), D.BRANCHING_BY_F1[(1, 1, 1)],
         STYLE['c_decay'], 2.8, '-', (X_CYC + 0.78, 2.55)),
        (X_CYC + 0.45, (X_F0, Y_X_F1_0), D.BRANCHING_BY_F1[(1, 1, 0)],
         STYLE['c_decay'], 1.7, '-', (X_F0 - 0.55, Y_X_F1_0 + 0.62)),
        (X_CYC + 0.70, (X_F2, Y_X_F1_2), D.BRANCHING_BY_F1[(1, 1, 2)],
         STYLE['c_decay'], 1.7, '-', (X_F2 + 0.05, Y_X_F1_2 + 0.78)),
        (X_CYC + 0.95, (X_G0, Y_X_G0), D.BRANCHING_BY_F1[(1, 0, 1)],
         STYLE['c_leak'], 1.2, ':', (0.50 * (X_CYC + 0.95) + 0.50 * X_G0, 3.55)),
    ]
    for sx, (tx, ty), f, col, lw, ls, (lx, ly) in decays:
        wavy(ax, sx, Y_A_F1_0 - 0.06, tx, ty, col, lw, ls=ls,
             alpha=1.0 if f > 1e-3 else 0.95)
        ax.text(lx, ly, rf'$f={_fmt_f(f)}$', color=col,
                ha='center', va='center', fontsize=STYLE['f_fs'], zorder=5,
                bbox=dict(boxstyle='round,pad=0.14', fc='white', ec='none', alpha=0.85))

    # ---- N=3 leak as a small annotation ----
    f_n3 = D.BRANCHING_BY_F1[(3, 1, 2)]
    ax.text(X_G0, Y_X_G0 - 0.95, rf'$\to N{{=}}3:\ f={_fmt_f(f_n3)}$',
            ha='center', va='top', fontsize=9, color=STYLE['c_dark'])

    # ---- state titles ----
    ax.text(-0.3, Y_A_F1_0 + 0.55, r"$A\,^2\Pi_{1/2}\ (v{=}0,\ J'{=}1/2,\ +)$",
            fontsize=STYLE['state_fs'], fontweight='bold', color=STYLE['c_cycle'])
    ax.text(-0.3, -0.55, r'$X\,^2\Sigma^+\ (v{=}0,\ N{=}1)$',
            fontsize=STYLE['state_fs'], fontweight='bold', color=STYLE['c_level'])

    # ---- footnote ----
    ax.text(-0.3, -1.25,
            r"$^{225}$RaF  $A(0)\!\to\!X(0)$ hyperfine branching, grouped by $F_1$ "
            r"(summed over $F$); fractions conditional on $v''{=}0$. Energy not to scale.",
            fontsize=9.5, color='0.3', ha='left', va='top')

    ax.set_xlim(-0.7, X_G0 + 2.4)
    ax.set_ylim(-1.7, Y_A_F1_0 + 1.05)
    ax.axis('off')

    base = os.path.join(FIGDIR, '225RaF_XA_hf_cycling_v1')
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(base + '.' + ext, dpi=300, bbox_inches='tight',
                    facecolor='white')
    print('WROTE', base + '.{pdf,svg,png}')
    plt.close(fig)


if __name__ == '__main__':
    main()
