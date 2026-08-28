"""Figure for the cross-figure reproduction of Ho et al. (2026) Fig. 5(b).

The freeze claims that the kernel reconstructed from Fig. 1(e) reproduces the
independently published Fig. 5(b) absorption curves.  `repro_yield.py` checks
that numerically but draws nothing, so the claim had no figure.  This is it.

Two panels, one per laser line, because that is the whole claim: Ho's
calculated absorption against the reconstruction, over the full published
0-120 GPa axis.  The v1 phenomenological model is drawn on the same axes to
show what a failing curve looks like -- at 532 nm it puts the maximum at
38 GPa instead of 17 and anticorrelates (r = -0.51) inside the audit window,
so the agreement below is not the kind any smooth model would achieve.

To be precise about what v1 is: it is anchored to Ho's *reported scalars*
(dE_ZPL(120 GPa) = 0.400 eV, S_abs 3.08 -> 4.61) but no parameter of it was
ever fitted to the Fig. 5(b) curves it is compared against here.  Nor does
fitting rescue it -- calibrating dE120 against the measured PL moves it to
0.555 eV and still leaves 34.8% pooled RMS against Ho's calculated
absorption, with r = -0.26 at 532 nm.  The disagreement is structural, not a
missing fit.

Each panel carries a residual strip beneath it.  The curves in the main panel
overlap too closely to read a 1.5% disagreement off them, and the residual is
not flat noise: at 532 nm it drifts monotonically from +2% to -3% across the
window, which is a systematic tilt worth seeing rather than a number to
assert in a title.

The measured PL of the same published figure is deliberately NOT drawn here.
It is a different observable -- a detection passband, the emission quantum
yield and charge-state conversion all enter it -- and `repro_yield.py` keeps
it in a separate comparison for that reason.  Its numbers still print below.

Two pressure ranges appear in every panel and must not be confused.

* Ho publishes both calculated curves over the **whole 0-120 GPa axis**, and
  that is the axis drawn here.
* `repro_yield.compare_ho_theory` scores each line only inside the range where
  Ho actually has measured points -- 4.7-51 GPa at 532 nm, 51-113.8 GPa at
  457 nm -- so that model and data are compared on the same footing.  That
  **audit window is shaded**, and the quoted RMS belongs to it.  The full-range
  RMS is annotated separately so neither number is mistaken for the other.

The 532 nm reference curve decays to zero and its digitised trace crosses
slightly below zero above 101 GPa, so fractional residuals there are
meaningless; the residual strip is drawn only where the reference exceeds
`RESIDUAL_FLOOR` of its own maximum.

Every curve in Fig. 5(b) is separately normalised in the source, so only the
shape of each is meaningful.  Each model curve therefore carries the single
least-squares multiplicative scale the audit allows, and nothing else.

Writes `repro_ho_fig5b.png`.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ho_spectrum_model import HoPublishedSpectrumModel
from nv_model import NVModel
from repro_yield import (
    HO_MAX_FRACTIONAL_RMS, LINES, _scaled_residual_from_prediction, load,
    predict_absorption,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'repro_ho_fig5b.png')

COLOUR = {532.0: 'tab:green', 457.0: 'tab:blue'}
REF_COLOUR = '0.35'
WINDOW_COLOUR = '0.90'
# Fraction of its own maximum below which the digitised reference is too close
# to zero for a fractional residual to mean anything.
RESIDUAL_FLOOR = 0.05


def audit_window(data, lam):
    """The pressure span the audit scores: where Ho has measured points."""
    measured = data['expt%d' % lam][0]
    return float(measured.min()), float(measured.max())


def scaled(prediction, reference):
    """Prediction after the one allowed multiplicative scale."""
    residual, scale = _scaled_residual_from_prediction(prediction, reference)
    return scale * np.asarray(prediction, float), residual


def main():
    data = load()
    published = HoPublishedSpectrumModel()
    frozen = NVModel(T=90.0)

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 5.8), sharex='col',
                             gridspec_kw=dict(height_ratios=(3.0, 1.0),
                                              hspace=0.08))
    summary, pooled_window = {}, []

    for column, lam in enumerate(LINES):
        top, strip = axes[0][column], axes[1][column]
        pressure, reference = data['theory%d_ho' % lam]
        lo, hi = audit_window(data, lam)
        inside = (pressure >= lo) & (pressure <= hi)
        scorable = reference > RESIDUAL_FLOOR * reference.max()
        norm = reference.max()

        prediction = predict_absorption(published, lam, pressure)
        v1_raw = predict_absorption(frozen, lam, pressure)

        # the audit's own scale, fitted inside the audit window only
        recon_window, residual_window = scaled(prediction[inside],
                                               reference[inside])
        recon = recon_window.max() / prediction[inside].max() * prediction
        v1 = (scaled(v1_raw[inside], reference[inside])[0].max()
              / v1_raw[inside].max()) * v1_raw
        residual = (recon[scorable] - reference[scorable]) / reference[scorable]

        top.axvspan(lo, hi, color=WINDOW_COLOUR, zorder=0,
                    label='audit window (measured range)')
        top.plot(pressure, reference / norm, color=REF_COLOUR, lw=4.0,
                 alpha=0.35, label='Ho Fig. 5(b), published')
        top.plot(pressure, recon / norm, color=COLOUR[lam], lw=1.8,
                 label='Fig. 1(e) reconstruction')
        top.plot(pressure, v1 / norm, color='tab:red', lw=1.2, ls='--',
                 alpha=0.8, label='v1 phenomenological model')

        rms_window = float(np.sqrt(np.mean(residual_window ** 2)))
        rms_full = float(np.sqrt(np.mean(residual ** 2)))
        corr = float(np.corrcoef(recon, reference)[0, 1])
        corr_window = float(np.corrcoef(recon[inside], reference[inside])[0, 1])
        v1_corr = float(np.corrcoef(v1[inside], reference[inside])[0, 1])
        summary[lam] = (lo, hi, rms_window, rms_full, corr, corr_window,
                        v1_corr, float(pressure[scorable].min()),
                        float(pressure[scorable].max()),
                        float(pressure[np.argmax(recon)]),
                        float(pressure[np.argmax(reference)]),
                        float(pressure[np.argmax(v1)]))
        pooled_window.append(residual_window)

        top.set_title(f'({"ab"[column]}) {lam:.0f} nm  --  audit window '
                      f'{rms_window*100:.1f}% RMS, r = {corr_window:.2f}',
                      fontsize=11)
        top.text(0.03 if lam == 532.0 else 0.60, 0.06,
                 f'full 0-120 GPa axis: {rms_full*100:.1f}% RMS\n'
                 f'(scored over {pressure[scorable].min():.0f}-'
                 f'{pressure[scorable].max():.0f} GPa)',
                 transform=top.transAxes, fontsize=8.5, va='bottom',
                 bbox=dict(fc='white', ec='0.7', alpha=0.85, pad=3.0))
        top.set_ylabel(r'$\sigma_{\rm abs}$, scaled to the published curve',
                       fontsize=9.5)
        top.legend(fontsize=8.5,
                   loc='upper right' if lam == 532.0 else 'upper left')
        top.grid(alpha=0.25)

        strip.axvspan(lo, hi, color=WINDOW_COLOUR, zorder=0)
        strip.axhspan(-100.0 * HO_MAX_FRACTIONAL_RMS,
                      100.0 * HO_MAX_FRACTIONAL_RMS,
                      color='0.75', alpha=0.35, zorder=0)
        strip.axhline(0.0, color='k', lw=0.8)
        strip.plot(pressure[scorable], residual * 100.0, color=COLOUR[lam],
                   lw=1.0, alpha=0.55)
        band = inside & scorable
        strip.plot(pressure[band],
                   ((recon[band] - reference[band]) / reference[band]) * 100.0,
                   color=COLOUR[lam], lw=2.2)
        strip.set_xlim(0.0, 120.0)
        strip.set_ylim(-16.0, 16.0)
        strip.set_yticks((-10, 0, 10))
        strip.set_xlabel('pressure [GPa]')
        strip.set_ylabel('residual [%]', fontsize=9.5)
        strip.grid(alpha=0.25)

    pooled = float(np.sqrt(np.mean(np.concatenate(pooled_window) ** 2)))
    axes[1][0].text(0.03, 0.86, 'grey band: 10% audit tolerance;  '
                    'thick: audit window,  thin: rest of the axis',
                    transform=axes[1][0].transAxes, fontsize=8.0, va='top')

    fig.suptitle('Cross-figure reproduction: the kernel reconstructed from Ho '
                 f'Fig. 1(e), tested against Ho Fig. 5(b)  --  pooled '
                 f'{pooled*100:.1f}% RMS', fontsize=12.5)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.87, bottom=0.10,
                        wspace=0.19, hspace=0.08)
    fig.savefig(OUT, dpi=180)

    print(f'wrote {OUT}')
    print('  published axis: 0-120 GPa for both calculated curves')
    for lam in LINES:
        (lo, hi, rms_w, rms_f, corr, corr_w, v1_corr, s_lo, s_hi,
         peak_model, peak_ref, peak_v1) = summary[lam]
        print(f'  {lam:.0f} nm  audit window {lo:.1f}-{hi:.1f} GPa: '
              f'RMS {rms_w*100:.1f}%, r = {corr_w:.3f}, '
              f'peak {peak_model:.0f}/{peak_ref:.0f} GPa')
        print(f'  {lam:.0f} nm  full axis, scored {s_lo:.0f}-{s_hi:.0f} GPa '
              f'(reference above {RESIDUAL_FLOOR:.0%} of max): '
              f'RMS {rms_f*100:.1f}%, r = {corr:.3f}')
        print(f'  {lam:.0f} nm  v1 model, same window: r = {v1_corr:.2f}, '
              f'peak {peak_v1:.0f} GPa')
    print(f'  pooled audit-window RMS: {pooled*100:.1f}%')

    print('\n  measured PL of the same figure is a different observable and is '
          'not drawn:')
    for lam in LINES:
        pressure, values = data['expt%d' % lam]
        _, residual = scaled(predict_absorption(published, lam, pressure),
                             values)
        print(f'    {lam:.0f} nm: RMS '
              f'{float(np.sqrt(np.mean(residual ** 2)))*100:.1f}% '
              '(see repro_yield.compare_experiment)')


if __name__ == '__main__':
    main()
