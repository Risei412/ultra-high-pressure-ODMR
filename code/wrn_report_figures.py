"""Figures for the WRN weekly reports.

One panel per report, drawn in the WRN palette (rule.txt):
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


def fig_optimum_penalty(path, band_penalty=1.05, ceiling=2.0):
    """Optical-limit penalty eta/eta_min = sqrt(A_max/A): the answer under (M)."""
    kernel = Kernel()
    lo, hi = DATA_WINDOW
    grid = np.arange(lo, hi + 1e-9, 0.01)
    penalty = 1.0 / np.sqrt(kernel.a(grid))
    blue, red = kernel.tolerance_band(band_penalty)

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    ax.plot(grid, penalty, color=PRIMARY, lw=1.6, zorder=3,
            label=r'$\eta/\eta_{\min}$, optical limit, 120 GPa')
    ax.axvspan(blue, red, color=PALE, zorder=1,
               label=f'{100 * (band_penalty - 1):.0f}% tolerance band')
    ax.axvline(kernel.lam_abs, color=ACCENT, lw=1.2, ls='--', zorder=2,
               label=f'optimum {kernel.lam_abs:.2f} nm')

    ax.set_xlim(lo, hi)
    ax.set_ylim(1.0, ceiling)
    ax.set_xlabel(r'Excitation wavelength  $\lambda$ (nm)')
    ax.set_ylabel(r'Sensitivity penalty  $\eta(\lambda)/\eta_{\min}$')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=2,
              fontsize=8, borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_split_visibility(path, sigma_eta=0.01, targets=(1.8646, 4.2989)):
    """Smallest eta bump between adjacent optima against power (L2, L5)."""
    from ladder_detection import SIGMA_MULTIPLE, plateau_structure

    ratios = np.unique(np.concatenate([
        np.geomspace(1.01, 6.0, 120),
        np.array(targets) * 0.999, np.array(targets) * 1.001,
    ]))
    bumps = []
    for ratio in ratios:
        structure = plateau_structure(float(ratio))
        bumps.append(min(structure['bumps']) if structure['bumps'] else np.nan)
    bumps = 100.0 * np.array(bumps)
    # A bump of exactly zero means the separating maximum is unresolvable on the
    # 0.01 nm grid; it has no place on a log axis.
    bumps[bumps <= 0.0] = np.nan

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    ax.plot(ratios, bumps, color=PRIMARY, lw=1.6, zorder=3,
            label='smallest $\\eta$ bump between optima, 120 GPa')
    ax.axhline(100.0 * SIGMA_MULTIPLE * sigma_eta, color=DARK, lw=1.0, ls='--',
               zorder=2, label=f'{SIGMA_MULTIPLE:.0f}$\\sigma$ at '
                               f'{100 * sigma_eta:.0f}% precision on $\\eta$')
    for index, target in enumerate(targets):
        ax.axvline(target, color=ACCENT, lw=1.0, ls=':', zorder=1,
                   label='target transitions' if index == 0 else None)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1.0, 6.0)
    ax.set_ylim(1e-3, 1e3)
    ax.set_xticks([1, 1.5, 2, 3, 4, 6])
    ax.set_xticklabels(['1', '1.5', '2', '3', '4', '6'])
    ax.set_xlabel(r'Excitation power ratio  $I / I_c$')
    ax.set_ylabel(r'Smallest $\eta$ bump (%)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=2,
              fontsize=8, borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    fig_zero_contrast(os.path.join(OUT, 'wrn_zero_contrast_width.png'))
    fig_ladder(os.path.join(OUT, 'wrn_ladder_steps.png'))
    fig_level_set(os.path.join(OUT, 'wrn_level_set_members.png'))
    fig_optimum_penalty(os.path.join(OUT, 'wrn_optimum_penalty.png'))
    fig_split_visibility(os.path.join(OUT, 'wrn_split_visibility.png'))

    kernel = Kernel()
    print('level-set members at a=0.60:',
          [round(x, 2) for x in kernel.level_set(0.60)])
    print('step boundaries <= 6:',
          [round(c.power_ratio, 4) for c in critical_values(kernel)
           if c.power_ratio <= 6.0])


if __name__ == '__main__':
    main()
