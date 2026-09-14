"""Figures for Addendum A3.

(a) the two branches tracked across pressure, in wavelength;
(b) the log-ratio crossing zero -- the exchange itself;
(c) the driver: the Franck-Condon displacement and the asymmetric decay;
(d) the resulting optimal wavelength, discontinuous at P*.

Writes `a3_branch_exchange.png`.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

from theory_a3_branch_exchange import (
    coupling_growth, exchange_pressure, identify_branches, monotonicity,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'a3_branch_exchange.png')
ZPL_COLOUR, SB_COLOUR = 'tab:red', 'tab:blue'


def panel_branches(ax, branches, star):
    pressures = np.array([b.pressure for b in branches['zpl']])
    zpl = np.array([b.wavelength_nm for b in branches['zpl']])
    sideband = np.array([b.wavelength_nm for b in branches['sideband']])
    fine = np.linspace(pressures[0], pressures[-1], 400)
    for values, colour, label in ((zpl, ZPL_COLOUR, 'ZPL branch'),
                                  (sideband, SB_COLOUR, 'sideband branch')):
        ax.plot(fine, CubicSpline(pressures, values)(fine), color=colour, lw=1.4)
        ax.plot(pressures, values, 'o', color=colour, ms=5, label=label)
    ax.axvline(star['pressure'], color='k', ls='--', lw=1.1)
    ax.annotate(f"$P^*$ = {star['pressure']:.0f} GPa", (star['pressure'] - 3, 455),
                rotation=90, fontsize=8, ha='right', va='bottom')
    ax.annotate('637 nm\n(NV$^-$ ZPL)', (2, 637), fontsize=7.5, color=ZPL_COLOUR,
                va='center')
    ax.set_xlabel('pressure [GPa]')
    ax.set_ylabel('branch peak wavelength [nm]')
    ax.set_title('(a) two branches, tracked in the raw samples', fontsize=9.5)
    ax.legend(fontsize=7.5, loc='upper right')


def panel_crossing(ax, branches, star):
    mono = monotonicity(branches)
    pressures, log_ratio = mono['pressures'], mono['log_ratio']
    fine = np.linspace(pressures[0], pressures[-1], 400)
    ax.plot(fine, CubicSpline(pressures, log_ratio)(fine), color='k', lw=1.4)
    ax.plot(pressures, log_ratio, 'o', color='k', ms=5)
    ax.axhline(0.0, color='0.5', lw=0.9)
    ax.axvline(star['pressure'], color='k', ls='--', lw=1.1)
    ax.fill_between(fine, -2, 0, color=ZPL_COLOUR, alpha=0.10)
    ax.fill_between(fine, 0, 2, color=SB_COLOUR, alpha=0.10)
    ax.annotate('ZPL branch wins', (4, -0.75), fontsize=8, color=ZPL_COLOUR)
    ax.annotate('sideband wins', (91, 0.44), fontsize=8, color=SB_COLOUR)
    ax.annotate(f"$P^*$ = {star['pressure']:.1f} GPa",
                (star['pressure'] - 3, -1.05), rotation=90, fontsize=8, ha='right')
    ax.set_xlabel('pressure [GPa]')
    ax.set_ylabel(r'$\ln\,(A_{\rm SB}/A_{\rm ZPL})$')
    ax.set_title('(b) monotone, one sign change: a unique exchange',
                 fontsize=9.5)
    ax.set_ylim(-1.2, 0.6)


def panel_driver(ax, branches):
    growth = coupling_growth(branches)
    pressures = growth['pressures']
    ax.plot(pressures, growth['gap_ev'], 'o-', color='tab:green', lw=1.4, ms=5)
    ax.set_xlabel('pressure [GPa]')
    ax.set_ylabel(r'$S\,\hbar\omega$  [eV]', color='tab:green')
    ax.tick_params(axis='y', labelcolor='tab:green')
    ax.set_title('(c) the driver: coupling grows, ZPL collapses', fontsize=9.5)
    ax.set_ylim(0.20, 0.44)

    twin = ax.twinx()
    for name, colour in (('zpl', ZPL_COLOUR), ('sideband', SB_COLOUR)):
        sigma = np.array([b.sigma for b in branches[name]])
        twin.semilogy(pressures, sigma / sigma[0], 's--', color=colour, ms=4,
                      lw=1.1)
    twin.set_ylabel(r'$\sigma / \sigma(0)$ (dashed)', fontsize=8.5)
    twin.tick_params(labelsize=8)
    twin.annotate(r'ZPL: $\times$1/6.0', (52, 0.30), fontsize=7.5,
                  color=ZPL_COLOUR)
    twin.annotate(r'sideband: $\times$1/1.3', (14, 0.70), fontsize=7.5,
                  color=SB_COLOUR)


def panel_optimum(ax, branches, star):
    pressures = np.array([b.pressure for b in branches['zpl']])
    zpl = CubicSpline(pressures, [b.wavelength_nm for b in branches['zpl']])
    sideband = CubicSpline(pressures,
                           [b.wavelength_nm for b in branches['sideband']])
    left = np.linspace(pressures[0], star['pressure'], 200)
    right = np.linspace(star['pressure'], pressures[-1], 200)
    ax.plot(left, zpl(left), color=ZPL_COLOUR, lw=2.4)
    ax.plot(right, sideband(right), color=SB_COLOUR, lw=2.4)
    # the branches that are present but no longer optimal
    ax.plot(left, sideband(left), color=SB_COLOUR, lw=1.0, ls=':', alpha=0.7)
    ax.plot(right, zpl(right), color=ZPL_COLOUR, lw=1.0, ls=':', alpha=0.7)
    ax.plot([star['pressure']] * 2, [star['sideband_nm'], star['zpl_nm']],
            'k--', lw=1.2)
    ax.plot([star['pressure']] * 2, [star['sideband_nm'], star['zpl_nm']],
            'ko', ms=6)
    ax.annotate(f"jump {star['jump_nm']:.0f} nm",
                (star['pressure'] + 2, 0.5 * (star['zpl_nm'] + star['sideband_nm'])),
                fontsize=8)
    ax.annotate('degenerate pair\nat zero power',
                (star['pressure'] - 4, 585), fontsize=7.5, ha='right')
    ax.set_xlabel('pressure [GPa]')
    ax.set_ylabel(r'$\lambda$ of the global optimum [nm]')
    ax.set_title('(d) the optimum is discontinuous in pressure', fontsize=9.5)


def main():
    branches = identify_branches()
    star = exchange_pressure(branches)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.8))
    panel_branches(axes[0, 0], branches, star)
    panel_crossing(axes[0, 1], branches, star)
    panel_driver(axes[1, 0], branches)
    panel_optimum(axes[1, 1], branches, star)
    fig.suptitle('Addendum A3: pressure exchanges which branch carries the '
                 'sensitivity optimum', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT, dpi=170)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
