"""Theorem X figure: the branch ratio is a density, so P* is a function of W.

Theorem X (`theory_a3_branch_exchange.py`) says the global sensitivity optimum
switches branch -- zero-phonon line to phonon sideband -- at the pressure where
the two branches' weighted peaks cross, and that the crossing is unique
whenever the log ratio is monotone.  Erratum E3 (`figure_validation.py`) found
that the ZPL peak heights in the extracted CSV are a clipping artefact of the
source figure, so the crossing must be rebuilt from the published Huang-Rhys
and Debye-Waller panels instead.

Doing so changes the *type* of the answer.  The sideband contributes a
lineshape density (per eV) and the ZPL contributes a dimensionless weight, so
their ratio r(P) = lambda_SB sigma_SB / (lambda_ZPL DWF_abs) carries units of
1/eV.  The measurable branch ratio is r(P) W, with W the excitation bandwidth.

Panel (a) shows r(P): monotone increasing over 0-120 GPa, which is Theorem X's
antecedent, so the crossing is still unique.  The right-hand axis is the same
quantity read as the critical bandwidth 1/r, the laser linewidth at which the
two branches are equally good at that pressure.  Panel (b) shows the
consequence: P* slides across the whole accessible range as W changes, and
leaves it entirely outside a narrow band of bandwidths.

What this figure must NOT be used to claim
------------------------------------------
* **Never use the panel (e) ZPL peak heights.**  All seven spikes are drawn to
  within one unit of the axis top; the CSV heights 6.364 ... 1.054 are a
  digitisation artefact.  Every ZPL quantity here comes from the published
  Debye-Waller factor in panel (c) instead.  The withdrawn numbers -- the
  "x6.04 ZPL collapse" and the single exchange pressure P* = 87.9 GPa -- must
  not be quoted from this figure or any other.
* **P* = 87.9 GPa is not recoverable from panel (b).**  There is no single
  exchange pressure.  Quoting one requires stating the excitation bandwidth it
  belongs to.
* **The shaded regions are not predictions of "no exchange".**  They mark
  bandwidths whose crossing falls outside the pressure range Ho et al.
  published; the extrapolation beyond 120 GPa is not made here.
* r(P) is a ratio of an extracted sideband peak to a published DWF.  Its
  absolute scale inherits the reconstruction; its monotonicity, which is what
  Theorem X needs, does not.

Run for the figure and the numbers quoted in the caption:

    python fig_x_branch_exchange.py
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import fig_style
from figure_validation import (branch_ratio_density,
                               drivers_from_published_panels,
                               exchange_pressure_at_bandwidth)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT_PNG = HERE / 'outputs' / 'figures' / 'fig_x_branch_exchange.png'
OUT_PDF = ROOT / 'paper' / 'figures' / 'fig_x_branch_exchange.pdf'

# Okabe-Ito, colourblind safe.
BLACK = '#000000'
BLUE = '#0072B2'
GREEN = '#009E73'
VERMILLION = '#D55E00'
GREY = '0.55'

# The accessible pressure range of the published series.
P_MIN, P_MAX = 0.0, 120.0
# Bandwidths called out with markers in panel (b), in meV.  They have to
# lie inside the accessible window, which the corrected panel (c) extraction
# moved from 1.49-9.68 meV to 2.39-10.30 meV.
MARKED_BANDWIDTHS_MEV = (2.5, 3.0, 4.0, 5.0, 7.0, 9.0)


def exchange_curve(n=361, lo_mev=2.0, hi_mev=11.0):
    """P*(W) over a bandwidth range that brackets the accessible window."""
    widths = np.linspace(lo_mev, hi_mev, n)
    stars = np.array([exchange_pressure_at_bandwidth(w * 1e-3)
                      for w in widths])
    return widths, stars


def accessible_bandwidths(data=None):
    """The bandwidth window inside which the crossing lands in 0-120 GPa.

    r is monotone increasing, so the crossing r(P) = 1/W sits inside the
    published range exactly when 1/W lies between r(120) and r(0); the
    endpoints are the critical bandwidths at the two ends of the series.
    """
    data = data or branch_ratio_density()
    widths = data['critical_bandwidth_meV']
    return float(np.min(widths)), float(np.max(widths))


def _ratio_panel(ax, data):
    pressure, ratio = data['pressure'], data['r']
    ax.semilogy(pressure, ratio, '-o', color=BLACK, lw=1.3, ms=4.0,
                mfc='white', mew=1.1, zorder=4,
                label=r'$r(P)=\lambda_{\rm SB}\sigma_{\rm SB}/'
                      r'(\lambda_{\rm ZPL}\,{\rm DWF_{abs}})$')
    ax.set_ylim(70.0, 1000.0)

    for index, (offset, align, vertical) in ((0, ((7, -10), 'left', 'top')),
                                             (-1, ((-6, 7), 'right', 'bottom'))):
        ax.annotate('%.1f  (%.2f meV)' % (ratio[index],
                                          data['critical_bandwidth_meV'][index]),
                    xy=(pressure[index], ratio[index]), xytext=offset,
                    textcoords='offset points', fontsize=fig_style.ANNOT,
                    color=BLACK, ha=align, va=vertical)

    # The growth factor, the monotonicity and what they buy Theorem X are
    # stated in the caption.  A single curve needs no legend either: the
    # ordinate label and Eq. (r) already define what is drawn.
    ax.set_xlim(-4.0, 124.0)
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel(r'branch ratio density $r$  (eV$^{-1}$)')
    ax.set_title('(a)', loc='left')
    ax.grid(which='both', alpha=0.20, lw=0.5)

    # The same axis read the other way round: 1/r is the laser linewidth at
    # which the two branches are equally good.  On a log scale that is an
    # exact reflection, so no second curve is drawn.
    twin = ax.twinx()
    twin.set_yscale('log')
    bottom, top = ax.get_ylim()
    twin.set_ylim(1e3 / bottom, 1e3 / top)
    twin.set_ylabel(r'critical bandwidth $1/r$  (meV)', color=BLUE)
    twin.tick_params(axis='y', colors=BLUE)
    twin.spines['right'].set_color(BLUE)


def _exchange_panel(ax, widths, stars, marked, window):
    ax.axhspan(P_MAX, 200.0, color=GREY, alpha=0.20, lw=0, zorder=1)
    ax.axhspan(-200.0, P_MIN, color=GREY, alpha=0.20, lw=0, zorder=1)
    ax.axhline(P_MAX, color=GREY, lw=0.7, ls='-', zorder=2)
    ax.axhline(P_MIN, color=GREY, lw=0.7, ls='-', zorder=2)

    good = np.isfinite(stars)
    ax.plot(widths[good], stars[good], '-', color=BLACK, lw=1.3, zorder=4,
            label=r'$P^{*}$ from $r(P^{*})\,W_{\mathrm{eff}}=1$')
    ax.plot([row[0] for row in marked], [row[1] for row in marked], 's',
            color=VERMILLION, ms=4.2, mfc='white', mew=1.3, zorder=5,
            label='quoted bandwidths')
    for width_mev, star in marked:
        ax.annotate('%.0f' % star, xy=(width_mev, star), xytext=(5, 4),
                    textcoords='offset points', fontsize=fig_style.ANNOT,
                    color=VERMILLION)

    # The bandwidth window and what the shading means are in the caption; the
    # dashed verticals get a legend entry instead of a sentence in the axes.
    lo, hi = window
    for edge in (lo, hi):
        ax.axvline(edge, color=GREEN, lw=0.9, ls=(0, (4, 1.5)), zorder=3,
                   label='_' if edge != lo else 'bandwidth window')

    ax.set_xlim(widths[0], widths[-1])
    ax.set_ylim(-22.0, 142.0)
    ax.set_yticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_xlabel(r'Effective excitation bandwidth $W_{\mathrm{eff}}$ (meV)')
    ax.set_ylabel('Exchange pressure $P^{*}$ (GPa)')
    ax.set_title('(b)', loc='left')
    ax.grid(alpha=0.22, lw=0.5)
    # The window rule at the upper bandwidth used to run straight through
    # the first legend row, and the legend sat unframed on the grey band.
    ax.legend(loc='upper right', handlelength=2.0, borderaxespad=0.6,
              labelspacing=0.35, frameon=True, framealpha=1.0,
              facecolor='white', edgecolor='0.8', fancybox=False)


def main():
    data = branch_ratio_density()
    drivers = drivers_from_published_panels()
    widths, stars = exchange_curve()
    marked = [(w, exchange_pressure_at_bandwidth(w * 1e-3))
              for w in MARKED_BANDWIDTHS_MEV]
    window = accessible_bandwidths(data)

    fig_style.use()
    fig, axes = plt.subplots(1, 2, figsize=(fig_style.FULL_WIDTH_IN, 3.35))
    # Panel (a) carries a right-hand axis, so the gutter has to hold that
    # label and panel (b)'s ordinate label without the two colliding.
    fig.subplots_adjust(left=0.100, right=0.945, bottom=0.155, top=0.930,
                        wspace=0.52)

    _ratio_panel(axes[0], data)
    _exchange_panel(axes[1], widths, stars, marked, window)

    print(fig_style.save(fig, OUT_PNG, OUT_PDF) + '\n')
    print('E3-corrected drivers (published panels (b), (c), never panel (e))')
    print('  S_abs   %.3f -> %.3f (+%.0f %%), monotone %s'
          % (drivers['S_abs_start'], drivers['S_abs_end'],
             drivers['S_abs_relative_growth'] * 100, drivers['S_abs_monotone']))
    print('  DWF_abs %.5f -> %.5f (falls x%.2f), monotone %s'
          % (drivers['DWF_start'], drivers['DWF_end'],
             drivers['DWF_fall_factor'], drivers['DWF_monotone']))

    print('\npanel (a): branch ratio density')
    print('    P [GPa]    r [1/eV]    1/r [meV]')
    for pressure, value, width in zip(data['pressure'], data['r'],
                                      data['critical_bandwidth_meV']):
        print('    %5.0f      %8.1f     %7.2f' % (pressure, value, width))
    print('  monotone increasing: %s, total growth x%.2f'
          % (data['monotone'], data['growth_factor']))

    print('\npanel (b): exchange pressure against excitation bandwidth')
    print('    W [meV]    P* [GPa]')
    for width_mev, star in marked:
        print('    %6.1f     %8.1f' % (width_mev, star))
    print('  crossing stays inside %.1f-%.1f GPa only for W in '
          '%.3f-%.3f meV' % (P_MIN, P_MAX, window[0], window[1]))
    for width_mev in (1.0, 10.0):
        star = exchange_pressure_at_bandwidth(width_mev * 1e-3)
        print('    W = %4.1f meV: %s' % (width_mev,
                                         'outside 0-120 GPa'
                                         if not np.isfinite(star)
                                         else '%.1f GPa' % star))


if __name__ == '__main__':
    main()
