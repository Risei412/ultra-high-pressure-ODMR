"""Verification figures for the esa article: is the Ho kernel fit to be used?

The article's numerical tables all descend from one external input, the
120 GPa absorption curve reconstructed from Ho et al. Fig. 1(e).  Nothing in
the article's own logic can tell whether that reconstruction is sound, so the
check has to come from the source paper itself.  Two independent checks exist
and this module draws one figure for each.

The cross-figure check.  Fig. 1(e) is reconstructed without ever looking at
Fig. 5(b), then compared with the independently digitised Fig. 5(b) curves at
the two wavelengths that figure plots.  Agreement there is evidence about the
reconstruction, not about the physics: it is a reproduction of published
curves, not an independent calculation.  See ``repro_yield.py`` for the audit
and ``docs/theory/sanity_check_pl_yield.md`` for the gates.

The internal checks K1 to K4 from ``kernel_sanity_checks.py``.  Fig. 5(b)
probes only 532 and 457 nm and says nothing about the zero-phonon line, which
is exactly where the extraction is known to be clipped (erratum E3).  K1 to K4
close that gap using material internal to the same paper, and they mark the
two places where the kernel must NOT be used: above the blue edge of the
plotted window, and for anything ZPL-derived above ~40 GPa.

Why two figures rather than one.  Both are full-width `figure*` floats, and a
`figure*` is included at \\textwidth, about 7.06 in.  The six panels used to be
drawn on one 15.6 in canvas, so the page rescaled every letter in them by
0.46: 10 pt labels arrived at 4.6 pt and the legends at 3.9 pt.  Each figure
is now drawn at the width it is printed at, and the numbers each panel used to
letter into itself live in the captions.

Run for both figures and the numbers quoted in the captions:

    python esa_figv_kernel_sanity.py
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import fig_style
from ho_spectrum_model import HoPublishedSpectrumModel
from kernel_sanity_checks import (ZPL_TRUSTED_MAX_GPA, jahn_teller_coupling,
                                  normalisation, zpl_area, zpl_width)
from repro_yield import (EXPT_MAX_FRACTIONAL_RMS, EXPT_MAX_PEAK_ERROR_GPA,
                         HO_MAX_FRACTIONAL_RMS, HO_MAX_PEAK_ERROR_GPA,
                         compare_experiment, compare_ho_theory, load,
                         predict_absorption)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIGURES = HERE / 'outputs' / 'figures'
PAPER_FIGURES = ROOT / 'paper' / 'figures'

OUT_PNG_CROSS = FIGURES / 'esa_figv_kernel_cross.png'
OUT_PDF_CROSS = PAPER_FIGURES / 'esa_figv_kernel_cross.pdf'
OUT_PNG_INTERNAL = FIGURES / 'esa_figv_kernel_internal.png'
OUT_PDF_INTERNAL = PAPER_FIGURES / 'esa_figv_kernel_internal.pdf'

COLOURS = {532: '#238b45', 457: '#2166ac'}
RED = '#b2182b'
GREEN = '#238b45'


def _scaled_reconstruction(model, wavelength, pressure, reference):
    """The one multiplicative scale the audit allows; no shape freedom."""
    prediction = predict_absorption(model, wavelength, pressure)
    scale = np.sum(prediction * reference) / np.sum(prediction * prediction)
    return scale * prediction


# --------------------------------------------------- the cross-figure check

def _cross_figure_panel(ax, model, data, wavelength, tag, peak_reference):
    pressure, reference = data['theory%d_ho' % wavelength]
    measured = data['expt%d' % wavelength][0]
    window = (pressure >= measured.min()) & (pressure <= measured.max())
    pressure, reference = pressure[window], reference[window]
    reconstruction = _scaled_reconstruction(
        model, wavelength, pressure, reference)

    ax.plot(pressure, reference, color='0.15', lw=1.8)
    ax.plot(pressure, reconstruction, color=COLOURS[wavelength], lw=1.6,
            ls='--')
    ax.axvline(peak_reference, color='0.15', alpha=0.25, lw=0.9)
    # The RMS residual and the two peak pressures are in the caption; the
    # title says only which wavelength this is.
    ax.set_title('(%s) %d nm' % (tag, wavelength), loc='left')
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('Absorption (scaled a.u.)')
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.22, lw=0.5)


def _gate_panel(ax, theory, experiment):
    labels = ['Ho\n532', 'Ho\n457', 'Expt\n532', 'Expt\n457']
    rms = [theory['532']['fractional_rms'], theory['457']['fractional_rms'],
           experiment['532']['fractional_rms'],
           experiment['457']['fractional_rms']]
    peak = [theory['532']['peak_error_GPa'], theory['457']['peak_error_GPa'],
            experiment['532']['peak_error_GPa'],
            experiment['457']['peak_error_GPa']]
    ceilings = [HO_MAX_FRACTIONAL_RMS] * 2 + [EXPT_MAX_FRACTIONAL_RMS] * 2
    peak_ceilings = [HO_MAX_PEAK_ERROR_GPA] * 2 + [EXPT_MAX_PEAK_ERROR_GPA] * 2

    x = np.arange(len(labels))
    bars = ax.bar(x, [value * 100 for value in rms],
                  color=['#238b45', '#2166ac', '#74c476', '#6baed6'],
                  zorder=3)
    ax.scatter(x, [value * 100 for value in ceilings], marker='_', s=330,
               linewidths=2.2, color=RED, zorder=4)
    for bar, value in zip(bars, rms):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 100 + 0.4,
                '%.1f%%' % (value * 100), ha='center', va='bottom',
                fontsize=fig_style.ANNOT, zorder=5)
    ax.set_xticks(x, ['%s\n%.0f/%.0f' % (label, error, ceiling)
                      for label, error, ceiling
                      in zip(labels, peak, peak_ceilings)])
    ax.set_xlabel('peak error / ceiling (GPa)')
    ax.set_ylabel('Fractional RMS error (%)')
    ax.set_title('(c)', loc='left')
    ax.set_ylim(0, 19)
    ax.grid(axis='y', alpha=0.22, lw=0.5)


def cross_figure(model, data, theory, experiment):
    """Panels (a)-(c): the reconstruction against the curves it never saw."""
    fig, axes = plt.subplots(1, 3, figsize=(fig_style.FULL_WIDTH_IN, 2.95))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.300, top=0.910,
                        wspace=0.38)

    _cross_figure_panel(axes[0], model, data, 532, 'a',
                        theory['532']['peak_reference_GPa'])
    _cross_figure_panel(axes[1], model, data, 457, 'b',
                        theory['457']['peak_reference_GPa'])
    _gate_panel(axes[2], theory, experiment)

    handles = [
        Line2D([], [], color='0.15', lw=1.8, label='Ho $et$ $al.$ Fig. 5(b), calculated'),
        Line2D([], [], color=COLOURS[532], lw=1.6, ls='--',
               label='reconstructed, 532 nm'),
        Line2D([], [], color=COLOURS[457], lw=1.6, ls='--',
               label='reconstructed, 457 nm'),
        Line2D([], [], color=RED, lw=0, marker='_', ms=9, mew=2.2,
               label='PASS ceiling'),
    ]
    fig.legend(handles=handles, frameon=False, ncol=4, loc='lower center',
               bbox_to_anchor=(0.5, 0.005), handlelength=2.2,
               columnspacing=1.6, borderaxespad=0.0)
    return fig


# ------------------------------------------------------ the internal checks

def _k1_panel(ax, norm):
    pressures = np.array(sorted(norm))
    area = np.array([norm[p]['area'] for p in pressures])
    top = np.array([norm[p]['top_fraction'] for p in pressures])

    # No legend: each curve is drawn in the colour of the axis that scales
    # it, which is what the two ordinate labels say.
    ax.plot(pressures, area, 'o-', color='0.15', lw=1.6, ms=3.6)
    ax.axhline(1.0, color='0.15', ls=':', lw=0.9)
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('K1 integrated area (norm.)')
    ax.set_ylim(0.88, 1.03)
    ax.set_title('(a) K1', loc='left')
    ax.grid(alpha=0.22, lw=0.5)

    twin = ax.twinx()
    twin.plot(pressures, top * 100, 's--', color=RED, lw=1.5, ms=3.4)
    twin.set_ylabel('$a(E>3$ eV$)$ / peak (%)', color=RED)
    twin.tick_params(axis='y', colors=RED)
    twin.spines['right'].set_color(RED)
    twin.set_ylim(0, 100)


def _k2_panel(ax, width):
    pressures = np.array(sorted(width))
    values = np.array([width[p] for p in pressures]) * 1000.0

    ax.plot(pressures, values, 'o-', color='0.15', lw=1.6, ms=3.6,
            label='apparent FWHM')
    ax.axhline(values.mean(), color=RED, ls='--', lw=1.3,
               label='mean %.2f meV' % values.mean())
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('Apparent ZPL FWHM (meV)')
    ax.set_ylim(0, 4)
    ax.set_title('(b) K2', loc='left')
    ax.grid(alpha=0.22, lw=0.5)
    ax.legend(frameon=False, loc='lower left', handlelength=2.0,
              borderaxespad=0.5, labelspacing=0.35)


def _k3_k4_panel(ax, area, jt):
    pressures = np.array(sorted(area))
    ratio = np.array([area[p]['ratio'] for p in pressures])
    total = np.array([jt[p]['S_total'] for p in pressures])

    # Shading: green, the pressures below which the kernel ZPL may be used;
    # grey, the +/-5 % band around parity.  The caption says so.
    ax.axvspan(0, ZPL_TRUSTED_MAX_GPA, color='#c7e9c0', alpha=0.55, lw=0)
    ax.plot(pressures, ratio, 'o-', color='0.15', lw=1.6, ms=3.6)
    ax.axhline(1.0, color='0.15', ls=':', lw=0.9)
    ax.axhspan(0.95, 1.05, color='0.15', alpha=0.10, lw=0)
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('K3 ZPL area / DWF')
    ax.set_ylim(0.9, 1.55)
    ax.set_title('(c) K3, K4', loc='left')
    ax.grid(alpha=0.22, lw=0.5)

    twin = ax.twinx()
    twin.plot(pressures, total, 's--', color=RED, lw=1.5, ms=3.4)
    twin.set_ylabel(r'K4 $S_{\mathrm{total}}$', color=RED)
    twin.tick_params(axis='y', colors=RED)
    twin.spines['right'].set_color(RED)
    twin.set_ylim(3.0, 7.0)


def internal_figure(norm, width, area, jt):
    """Panels (a)-(c): K1 to K4, the checks internal to the source paper."""
    fig, axes = plt.subplots(1, 3, figsize=(fig_style.FULL_WIDTH_IN, 2.70))
    # Two of the three panels carry a right-hand axis, so the gutters have
    # to hold two ordinate labels each without them colliding.
    fig.subplots_adjust(left=0.086, right=0.914, bottom=0.200, top=0.895,
                        wspace=0.86)

    _k1_panel(axes[0], norm)
    _k2_panel(axes[1], width)
    _k3_k4_panel(axes[2], area, jt)
    return fig


def main():
    model = HoPublishedSpectrumModel()
    data = load()
    theory = compare_ho_theory(model, data)
    experiment = compare_experiment(model, data, collection=False)
    norm, width = normalisation(model), zpl_width(model)
    area, jt = zpl_area(model), jahn_teller_coupling()

    fig_style.use()
    cross = cross_figure(model, data, theory, experiment)
    internal = internal_figure(norm, width, area, jt)

    print(fig_style.save(cross, OUT_PNG_CROSS, OUT_PDF_CROSS))
    print(fig_style.save(internal, OUT_PNG_INTERNAL, OUT_PDF_INTERNAL) + '\n')

    print('caption numbers, cross-figure check')
    for lam in ('532', '457'):
        print('  Ho %s nm : RMS %.1f%%, peak %.0f/%.0f GPa'
              % (lam, theory[lam]['fractional_rms'] * 100,
                 theory[lam]['peak_model_GPa'],
                 theory[lam]['peak_reference_GPa']))
    for lam in ('532', '457'):
        print('  PL %s nm : RMS %.1f%%, peak error %.0f GPa'
              % (lam, experiment[lam]['fractional_rms'] * 100,
                 experiment[lam]['peak_error_GPa']))

    print('\ncaption numbers, internal checks')
    print('  K1 area 0 GPa %.4f, 120 GPa %.4f, top fraction 120 GPa %.1f%%'
          % (norm[0.0]['area'], norm[120.0]['area'],
             norm[120.0]['top_fraction'] * 100))
    values = np.array([width[p] for p in sorted(width)]) * 1000.0
    print('  K2 FWHM mean %.2f meV, spread %.1f%%'
          % (values.mean(),
             (values.max() - values.min()) / values.mean() * 100))
    print('  K3 ratio 40 GPa %.2f, 120 GPa %.2f'
          % (area[40.0]['ratio'], area[120.0]['ratio']))
    print('  K4 S_total %.2f -> %.2f, S_JT %.2f -> %.2f'
          % (jt[0.0]['S_total'], jt[120.0]['S_total'],
             jt[0.0]['S_JT'], jt[120.0]['S_JT']))


if __name__ == '__main__':
    main()
