"""Figures for Addendum A2.

(a) how the level set sweeps down through the kernel as power rises;
(b) the multiplicity ladder N(I) with the predicted transition powers;
(c) the gauge plane 2c + s + 2w = E, and the E > 1 splitting criterion;
(d) gauge-equivalent models: identical eta, different C, R and dnu.

Writes `a2_multiplicity_120GPa.png`.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from theory_a1_generalization import DATA_WINDOW, Kernel
from theory_a2_multiplicity import (
    GaugeResponse, critical_values, gauge_family, match_transitions,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'a2_multiplicity_120GPa.png')


def panel_levels(ax, kernel):
    grid = np.arange(*DATA_WINDOW, 0.02)
    ax.plot(grid, kernel.a(grid), color='k', lw=1.4)
    colours = plt.cm.viridis(np.linspace(0.05, 0.85, 4))
    for colour, ratio in zip(colours, (1.20, 1.45, 2.50, 4.00)):
        level = 1.0 / ratio
        roots = kernel.level_set(level, DATA_WINDOW)
        ax.axhline(level, color=colour, lw=1.0, ls='--')
        ax.plot(roots, [level] * len(roots), 'o', color=colour, ms=5)
        ax.annotate(f'$I/I_c$={ratio:.2f}  N={len(roots)}',
                    (DATA_WINDOW[0] + 1.5, level + 0.015), fontsize=7.5,
                    color=colour)
    for critical in critical_values(kernel):
        if critical.kind in ('max', 'min') and critical.level > 0.2:
            ax.plot([critical.wavelength], [critical.level], 's', ms=4,
                    color='tab:red')
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel(r'$A/A_{\max}$')
    ax.set_title('(a) the optimum set is a level set of $A$', fontsize=9.5)
    ax.set_xlim(*DATA_WINDOW)
    ax.set_ylim(0.0, 1.08)


def panel_ladder(ax, kernel, matched):
    ladder = matched['ladder']
    ax.step(ladder['ratios'], ladder['counts'], where='post', color='k', lw=1.5)
    for row in matched['rows']:
        ax.axvline(row['predicted_power_ratio'], color='tab:red', ls=':', lw=1.0)
        ax.annotate(f"{row['predicted_from_nm']:.0f} nm",
                    (row['predicted_power_ratio'], 6.5), rotation=90,
                    fontsize=7, color='tab:red', ha='right', va='top')
    ax.set_xscale('log')
    ax.set_xlabel(r'$I/I_c$')
    ax.set_ylabel('number of equally optimal $\\lambda$')
    ax.set_title('(b) multiplicity ladder; steps predicted, not fitted',
                 fontsize=9.5)
    ax.set_ylim(0, 7)
    ax.set_yticks(range(0, 8))
    ax.annotate('N=6 spans only $\\times$1.003\n(not resolvable)',
                (1.55, 1.2), fontsize=7.5, color='0.3')


def panel_gauge_plane(ax):
    """Everything collapses onto one axis: the splitting exponent E."""
    exponent = np.linspace(1.02, 5.0, 400)
    ax.plot(exponent, 1.0 / (exponent - 1.0), color='k', lw=1.6)
    ax.axvline(1.0, color='tab:red', ls='--', lw=1.2)
    ax.axvspan(0.5, 1.0, color='tab:red', alpha=0.10)
    ax.annotate('$E\\leq1$: no split\nat any power', (0.55, 3.2), fontsize=7.5,
                color='tab:red')
    named = (('saturation alone', 1.0, 'tab:purple'),
             ('$\\Delta\\nu$ broadening alone', 1.0, 'tab:purple'),
             ('contrast collapse', 2.0, 'tab:red'),
             ('sat + broadening', 2.0, 'tab:green'),
             ('all three', 4.0, 'tab:blue'))
    offsets = {1.0: (0.06, 4.3), 2.0: (2.12, 1.05), 4.0: (4.05, 0.36)}
    seen = set()
    for label, value, colour in named:
        if value > 1.0:
            ax.plot([value], [1.0 / (value - 1.0)], '*', ms=13, color=colour)
    ax.annotate('contrast collapse\nsat + broadening\n(same $E$, same $\\rho^*$)',
                offsets[2.0], fontsize=7.5, color='0.2')
    ax.annotate('all three', offsets[4.0], fontsize=7.5, color='tab:blue')
    ax.set_xlabel('$E = 2c + s + 2w$')
    ax.set_ylabel(r'$\rho^*=\Gamma_p^*/\Gamma = 1/(E-1)$')
    ax.set_title('(c) one scalar decides the split and its power',
                 fontsize=9.5)
    ax.set_xlim(0.5, 5.0)
    ax.set_ylim(0.0, 5.0)


def panel_gauge_equivalence(ax):
    gamma = np.logspace(-2, 2, 300)
    first = gauge_family()['contrast collapse alone']
    second = gauge_family()['saturation + broadening']
    ax.loglog(gamma, first.phi(gamma) / first.phi(gamma).max(), color='k',
              lw=3.0, label=r'$\Phi$ (both models)')
    ax.loglog(gamma, second.phi(gamma) / second.phi(gamma).max(), color='w',
              lw=1.0, ls='--')
    for response, marker, name in ((first, '-', 'contrast collapse'),
                                   (second, '--', 'sat + broadening')):
        ax.loglog(gamma, response.contrast(gamma), marker, color='tab:red',
                  lw=1.2)
        ax.loglog(gamma, response.linewidth(gamma), marker, color='tab:blue',
                  lw=1.2)
        ax.loglog(gamma, response.rate(gamma) / response.rate(gamma).max(),
                  marker, color='tab:green', lw=1.2)
    ax.plot([], [], '-', color='tab:red', label='$C$')
    ax.plot([], [], '-', color='tab:blue', label=r'$\Delta\nu$')
    ax.plot([], [], '-', color='tab:green', label='$R$')
    ax.plot([], [], '-', color='0.5', label='solid / dashed = the two models')
    ax.set_xlabel(r'$\Gamma_p$')
    ax.set_ylabel('normalised response')
    ax.set_title(r'(d) same $\Phi$, different $C$, $\Delta\nu$, $R$',
                 fontsize=9.5)
    ax.legend(fontsize=6.5, loc='lower left')
    ax.set_ylim(1e-3, 3.0)


def main():
    kernel = Kernel()
    matched = match_transitions(kernel)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.8))
    panel_levels(axes[0, 0], kernel)
    panel_ladder(axes[0, 1], kernel, matched)
    panel_gauge_plane(axes[1, 0])
    panel_gauge_equivalence(axes[1, 1])
    fig.suptitle('Addendum A2: multiplicity ladder and gauge degeneracy, '
                 '120 GPa', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT, dpi=170)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
