"""Theorem G figure: one sensitivity surface, five different physical models.

`eta = 1/(G sqrt(R))` depends on the ODMR response only through the single
combination Phi = C^2 R / dnu^2.  Writing R ~ Gamma_p (1+rho)^-s,
C ~ (1+rho)^-c and dnu ~ (1+rho)^w, Theorem G
(`theory_a2_multiplicity.py`) says Phi = Gamma_p (1+rho)^-E with
E = 2c + s + 2w, so every response on a plane of constant E produces the
*identical* sensitivity surface, the identical Gamma_p*, and the identical
multiplicity ladder.

Panel (a) plots Phi for the five representatives of the plane 2c + s + 2w = 2
that `gauge_family()` defines.  They coincide to machine precision -- the
curves are drawn with different dashes and markers only so the reader can see
that five of them are present.  Panel (b) plots what the same five models
predict for the three individually measurable factors R, C and dnu at
Gamma_p = 1.  Those differ by up to a factor of two.

The two panels together are the identifiability statement: a wavelength scan
of eta, at any number of powers, cannot say which mechanism is acting; only a
separate measurement of the rate and of one of (C, dnu) can.

Numbers belong in the caption, not in the raster.  The residual, the factor
of two and the pairwise identifiability gaps are printed by `main()` for the
caption to quote; nothing is written inside the axes.  The five members share
one figure-level legend, because they are the same five in both panels.

What this figure must NOT be used to claim
------------------------------------------
* **Panel (a) is not evidence that the models agree with data.**  It is an
  algebraic identity between models, verified numerically.  No measurement
  enters this figure at all -- neither the Ho kernel nor any ODMR data.
* **It does not show that mechanism attribution is impossible in general.**
  It shows that *eta alone* cannot do it.  Panel (b) is precisely the recipe
  for doing it with more spectra.
* The vertical scale of (a) is normalised to the shared maximum; the residual
  `main()` prints is the largest relative deviation over the whole decade
  range, and is a floating-point residual, not a physical spread.
* The five members are representatives, not an exhaustive list: the plane
  2c + s + 2w = 2 is two-dimensional and every point on it behaves the same.

Run for the figure and the numbers quoted in the caption:

    python fig_g_gauge_degeneracy.py
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import fig_style
from theory_a2_multiplicity import (gauge_degeneracy, gauge_family,
                                    identifiability)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT_PNG = HERE / 'outputs' / 'figures' / 'fig_g_gauge_degeneracy.png'
OUT_PDF = ROOT / 'paper' / 'figures' / 'fig_g_gauge_degeneracy.pdf'

# Okabe-Ito, colourblind safe; each colour also carries a dash pattern, a
# marker and a hatch, so nothing depends on colour alone.
BLACK = '#000000'
VERMILLION = '#D55E00'

MEMBER_STYLE = {
    'contrast collapse alone': ('#000000', '-', 'o', ''),
    'saturation + broadening': ('#0072B2', (0, (6, 2)), 's', '///'),
    'linewidth^1 (strong)': ('#009E73', (0, (1, 1.4)), '^', '\\\\\\'),
    'mixed a': ('#E69F00', (0, (5, 1.5, 1, 1.5)), 'D', 'xxx'),
    'mixed b': ('#CC79A7', (0, (3, 1, 1, 1, 1, 1)), 'v', '...'),
}

# The keys above are the model names the tests and Table II use.  What the
# legend draws is kept separate so that 'linewidth^1' does not reach the
# page as a literal caret.
DISPLAY_LABEL = {'linewidth^1 (strong)': r'linewidth$^{1}$ (strong)'}


def _label(name):
    return DISPLAY_LABEL.get(name, name)


# Gamma_p at which the individually measurable factors are compared.
PROBE = 1.0
# Same grid gauge_degeneracy() uses to report the residual.
GAMMA_GRID = np.logspace(-3, 3, 600)


def phi_curves(gamma=GAMMA_GRID):
    """Normalised Phi for every member of the gauge plane."""
    members = gauge_family()
    reference = members['contrast collapse alone'].phi(gamma)
    scale = float(np.max(reference))
    return {name: response.phi(gamma) / scale
            for name, response in members.items()}


def measurable_factors(gamma_p=PROBE):
    """R, C and dnu at `gamma_p` -- what the gauge freedom does not hide."""
    rows = {row['name']: row for row in gauge_degeneracy()}
    return {name: {'R': rows[name]['rate_at_gamma_1'],
                   'C': rows[name]['contrast_at_gamma_1'],
                   'dnu': rows[name]['linewidth_at_gamma_1']}
            for name in gauge_family()}


def _legend_handles():
    """One set of proxies for both panels: same five models, same colours."""
    return [Line2D([], [], color=colour, ls=dashes, lw=1.3, marker=marker,
                   ms=3.6, mfc='white', mew=0.8, label=_label(name))
            for name, (colour, dashes, marker, _) in MEMBER_STYLE.items()]


def _phi_panel(ax, gamma, curves):
    for index, (name, (colour, dashes, marker, _)) in enumerate(
            MEMBER_STYLE.items()):
        ax.plot(gamma, curves[name], color=colour, ls=dashes, lw=1.3,
                marker=marker, ms=3.2, markevery=(28 + 22 * index, 140),
                mfc='white', mew=0.8, label=name, zorder=3 + index)

    ax.axvline(1.0, color='0.55', lw=0.7, ls=':', zorder=1)

    ax.set_xscale('log')
    ax.set_xlim(gamma[0], gamma[-1])
    ax.set_ylim(0.0, 1.12)
    ax.set_xlabel(r'Pump rate $\Gamma_p$ (units of $\Gamma$)')
    ax.set_ylabel(r'$\Phi=C^{2}R/\Delta\nu^{2}$  (normalized)')
    ax.set_title('(a)', loc='left')
    ax.grid(alpha=0.22, lw=0.5)


def _factor_panel(ax, factors):
    keys = ('R', 'C', 'dnu')
    labels = ('$R$\n(PL rate)', '$C$\n(contrast)',
              r'$\Delta\nu$' '\n(linewidth)')
    names = list(MEMBER_STYLE)
    width = 0.16
    centres = np.arange(len(keys), dtype=float)

    for index, name in enumerate(names):
        colour, _, _, hatch = MEMBER_STYLE[name]
        offset = (index - (len(names) - 1) / 2.0) * width
        ax.bar(centres + offset, [factors[name][k] for k in keys],
               width=width * 0.92, color=colour, alpha=0.85, hatch=hatch,
               edgecolor='0.15', lw=0.5, label=name, zorder=3)

    ax.set_xticks(centres, labels)
    ax.set_xlim(-0.55, len(keys) - 0.45)
    ax.set_ylim(0.0, 2.2)
    ax.set_ylabel(r'value at $\Gamma_p=%g$ (arb. units)' % PROBE)
    ax.set_title('(b)', loc='left')
    ax.grid(axis='y', alpha=0.22, lw=0.5)


def main():
    curves = phi_curves()
    rows = gauge_degeneracy()
    residual = max(row['max_relative_phi_difference'] for row in rows)
    factors = measurable_factors()
    gaps = identifiability(PROBE)

    fig_style.use()
    fig, axes = plt.subplots(1, 2, figsize=(fig_style.FULL_WIDTH_IN, 3.35))
    fig.subplots_adjust(left=0.078, right=0.985, bottom=0.225, top=0.935,
                        wspace=0.26)

    _phi_panel(axes[0], GAMMA_GRID, curves)
    _factor_panel(axes[1], factors)

    fig.legend(handles=_legend_handles(), frameon=False, ncol=5,
               loc='lower center', bbox_to_anchor=(0.5, 0.005),
               handlelength=2.6, columnspacing=1.3, borderaxespad=0.0)

    print(fig_style.save(fig, OUT_PNG, OUT_PDF) + '\n')
    print('panel (a): gauge plane 2c + s + 2w = 2')
    print('    model                       c     s     w      E      '
          'max |dPhi/Phi|   Gamma_p*')
    for row in rows:
        print('    %-26s %.2f  %.2f  %.2f   %.3f    %.2e      %.4f'
              % (row['name'], row['c'], row['s'], row['w'], row['exponent'],
                 row['max_relative_phi_difference'], row['gamma_star']))
    print('  largest residual across the family: %.3e' % residual)
    peak = {name: float(np.max(values)) for name, values in curves.items()}
    print('  normalised Phi maxima: '
          + ', '.join('%.6f' % value for value in peak.values()))

    print('\npanel (b): measurable factors at Gamma_p = %g' % PROBE)
    print('    model                        R        C       dnu')
    for name, row in factors.items():
        print('    %-26s %7.4f  %7.4f  %7.4f'
              % (name, row['R'], row['C'], row['dnu']))
    for key in ('R', 'C', 'dnu'):
        values = [row[key] for row in factors.values()]
        print('  %-4s spread across the family: x%.4f (%.4f to %.4f)'
              % (key, max(values) / min(values), min(values), max(values)))

    print('\nidentifiability gaps at Gamma_p = %g' % PROBE)
    for key, row in gaps.items():
        verdict = ('degenerate' if row['gap'] < 1e-12
                   else 'separates (%.0f %% apart)' % (row['gap'] * 100))
        print('    %-12s %-22s %s' % (key, row['quantity'], verdict))


if __name__ == '__main__':
    main()
