"""Figures for the WRN weekly reports.

Three panels, one per report, drawn in the WRN palette (rule.txt):
Science Blue #1C3177, Deep Indigo #4B0082, Mist Lavender #E8EAF1,
Midnight Ink #101426 on Optic White.

The kernel figures read the frozen v3 kernel; nothing here fits anything.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from theory_a1_generalization import Kernel, DATA_WINDOW
from theory_a2_multiplicity import critical_values

PRIMARY = '#1C3177'
ACCENT = '#4B0082'
PALE = '#E8EAF1'
DARK = '#101426'

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', 'WRN', 'image')

plt.rcParams.update({
    'font.size': 9,
    'axes.edgecolor': DARK,
    'axes.labelcolor': DARK,
    'text.color': DARK,
    'xtick.color': DARK,
    'ytick.color': DARK,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'legend.frameon': True,
    'legend.framealpha': 1.0,
    'legend.edgecolor': DARK,
})


def fig_zero_contrast(path):
    """Apparent width of the dark line as a function of detection threshold."""
    # C crosses zero linearly across the contour; x is in units of the
    # crossing scale L = C_0 / |dC/dx|, so the curve is threshold-free.
    x = np.linspace(-1.0, 1.0, 801)
    contrast = np.abs(x)

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    ax.fill_between(x, 0.0, contrast, color=PALE, zorder=1)
    ax.plot(x, contrast, color=PRIMARY, lw=1.8, zorder=3,
            label=r'$|C|/C_0$ across the contour')
    ax.axvline(0.0, color=ACCENT, lw=1.2, ls='--', zorder=2,
               label=r'$C=0$ contour')

    for threshold, style in ((0.5, '-'), (0.2, ':')):
        ax.plot([-threshold, threshold], [threshold, threshold],
                color=DARK, lw=1.0, ls=style, zorder=4)
        ax.annotate('', xy=(-threshold, threshold), xytext=(threshold, threshold),
                    arrowprops=dict(arrowstyle='<->', color=DARK, lw=0.9))
        ax.text(0.0, threshold + 0.025, f'${2 * threshold:.1f}\\,L$',
                ha='center', va='bottom', fontsize=8, color=DARK, zorder=5,
                bbox=dict(facecolor='white', edgecolor='none', pad=1.0))

    ax.set_xlabel(r'Distance across the contour  $x / L$,   $L=C_0/|\mathrm{d}C/\mathrm{d}x|$')
    ax.set_ylabel(r'Contrast  $|C| / C_0$')
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.30), ncol=2,
              fontsize=8, borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_ladder(path):
    """Multiplicity ladder N(I/I_c): the step boundaries to be propagated."""
    kernel = Kernel()
    criticals = [c for c in critical_values(kernel) if c.power_ratio <= 6.0]
    criticals.sort(key=lambda c: c.power_ratio)

    edges = [1.0] + [c.power_ratio for c in criticals] + [6.0]
    counts, current = [], 2
    for critical in criticals:
        counts.append(current)
        current += critical.delta_multiplicity
    counts.append(current)

    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    for index, count in enumerate(counts):
        lo, hi = edges[index], edges[index + 1]
        ax.plot([lo, hi], [count, count], color=PRIMARY, lw=2.2, solid_capstyle='butt',
                label='multiplicity $N$ (frozen kernel, 120 GPa)' if index == 0 else None)
    for index in range(len(counts) - 1):
        ax.plot([edges[index + 1]] * 2, [counts[index], counts[index + 1]],
                color=PRIMARY, lw=0.8, ls=':')

    for critical in criticals:
        ax.axvline(critical.power_ratio, color=ACCENT, lw=0.8, alpha=0.55)
    ax.axvline(criticals[0].power_ratio, color=ACCENT, lw=0.8, alpha=0.55,
               label='step boundary $I_k/I_c = A_{\\max}/A_k$')

    ax.set_xscale('log')
    ax.set_xlim(1.0, 6.0)
    ax.set_ylim(0.5, 7.0)
    ax.set_xticks([1, 1.5, 2, 3, 4, 6])
    ax.set_xticklabels(['1', '1.5', '2', '3', '4', '6'])
    ax.set_yticks([1, 2, 3, 4, 5, 6])
    ax.set_xlabel(r'Excitation power ratio  $I / I_c$')
    ax.set_ylabel(r'Multiplicity  $N$')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=1,
              fontsize=8, borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_level_set(path, level=0.60):
    """Kernel with one level line: the members the (M) test must compare."""
    kernel = Kernel()
    lo, hi = DATA_WINDOW
    grid = np.arange(lo, hi + 1e-9, 0.01)
    values = kernel.a(grid)
    members = kernel.level_set(level)

    fig, ax = plt.subplots(figsize=(4.4, 3.3))
    ax.fill_between(grid, 0.0, values, color=PALE, zorder=1)
    ax.plot(grid, values, color=PRIMARY, lw=1.6, zorder=3,
            label=r'kernel $a(\lambda)$, 120 GPa')
    ax.axhline(level, color=DARK, lw=1.0, ls='--', zorder=2,
               label=f'level $a={level:.2f}$  ($I/I_c={1 / level:.2f}$)')
    ax.plot(members, [level] * len(members), 'o', ms=5.5, color=ACCENT, zorder=4,
            label='level-set members')

    ax.annotate(f'{members[0]:.1f} nm', xy=(members[0], level),
                xytext=(members[0] + 2, level + 0.30), fontsize=8, color=ACCENT,
                arrowprops=dict(arrowstyle='-', color=ACCENT, lw=0.8))
    ax.annotate(f'{members[-2]:.2f} / {members[-1]:.2f} nm',
                xy=(members[-1], level), xytext=(members[-1] - 2, level + 0.30),
                ha='right', fontsize=8, color=ACCENT,
                arrowprops=dict(arrowstyle='-', color=ACCENT, lw=0.8))

    ax.set_xlim(lo, hi)
    ax.set_ylim(0.0, 1.10)
    ax.set_xlabel(r'Excitation wavelength  $\lambda$ (nm)')
    ax.set_ylabel(r'Normalised absorbed-photon rate  $a(\lambda)$')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.30), ncol=2,
              fontsize=8, borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    fig_zero_contrast(os.path.join(OUT, 'wrn_zero_contrast_width.png'))
    fig_ladder(os.path.join(OUT, 'wrn_ladder_steps.png'))
    fig_level_set(os.path.join(OUT, 'wrn_level_set_members.png'))

    kernel = Kernel()
    print('level-set members at a=0.60:',
          [round(x, 2) for x in kernel.level_set(0.60)])
    print('step boundaries <= 6:',
          [round(c.power_ratio, 4) for c in critical_values(kernel)
           if c.power_ratio <= 6.0])


if __name__ == '__main__':
    main()
