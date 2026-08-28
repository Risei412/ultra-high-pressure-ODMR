"""Figure for the cross-figure reproduction of Ho et al. (2026) Fig. 5(b).

The freeze claims that the kernel reconstructed from Fig. 1(e) reproduces the
independently published Fig. 5(b) absorption curves.  `repro_yield.py` checks
that numerically but draws nothing, so the claim had no figure.  This is it.

(a),(b) the reconstruction against Ho's calculated absorption at the two laser
        lines, over the full published 0-120 GPa axis, with the v1
        phenomenological model on the same axes to show what a failing curve
        looks like;
(c)     the fractional residual of (a),(b);
(d)     the measured PL of the same figure, kept in its own panel because it is
        a different observable: a detection passband may enter it, so it must
        not be pooled with the absorption comparison.

Two pressure ranges appear in every panel and must not be confused.

* Ho publishes both calculated curves over the **whole 0-120 GPa axis**, and
  that is the axis drawn here.
* `repro_yield.compare_ho_theory` scores each line only inside the range where
  Ho actually has measured points -- 4.7-51 GPa at 532 nm, 51-113.8 GPa at
  457 nm -- because the audit exists to compare model and data on the same
  footing.  That **audit window is shaded**, and the quoted RMS belongs to it.
  The full-range RMS is annotated separately so neither number is mistaken for
  the other.

The 532 nm reference curve decays to zero and its digitised trace crosses
slightly below zero above 101 GPa, so fractional residuals there are
meaningless; panel (c) shows the residual only where the reference exceeds
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


def published_curve(data, lam):
    """Ho's calculated absorption over the full published pressure axis."""
    return data['theory%d_ho' % lam]


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

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    residuals, summary = {}, {}

    for ax, lam in zip(axes[0], LINES):
        pressure, reference = published_curve(data, lam)
        lo, hi = audit_window(data, lam)
        inside = (pressure >= lo) & (pressure <= hi)
        scorable = reference > RESIDUAL_FLOOR * reference.max()
        norm = reference.max()

        prediction = predict_absorption(published, lam, pressure)
        v1_raw = predict_absorption(frozen, lam, pressure)

        # the audit's own scale, fitted inside the audit window only
        recon, residual_window = scaled(prediction[inside], reference[inside])
        recon_full = recon.max() / prediction[inside].max() * prediction
        v1_full = (scaled(v1_raw[inside], reference[inside])[0].max()
                   / v1_raw[inside].max()) * v1_raw

        # full-range residual, on the part of the axis where it is meaningful
        residual_full = ((recon_full[scorable] - reference[scorable])
                         / reference[scorable])
        residuals[lam] = (pressure, inside, scorable, recon_full, reference)

        ax.axvspan(lo, hi, color=WINDOW_COLOUR, zorder=0,
                   label='audit window (measured range)')
        ax.plot(pressure, reference / norm, color=REF_COLOUR, lw=4.0,
                alpha=0.35, label='Ho Fig. 5(b), published')
        ax.plot(pressure, recon_full / norm, color=COLOUR[lam], lw=1.8,
                label='Fig. 1(e) reconstruction')
        ax.plot(pressure, v1_full / norm, color='tab:red', lw=1.2, ls='--',
                alpha=0.8, label='v1 phenomenological model')

        rms_window = float(np.sqrt(np.mean(residual_window ** 2)))
        rms_full = float(np.sqrt(np.mean(residual_full ** 2)))
        corr = float(np.corrcoef(recon_full, reference)[0, 1])
        summary[lam] = (lo, hi, rms_window, rms_full, corr,
                        float(pressure[scorable].min()),
                        float(pressure[scorable].max()))

        ax.set_title(f'({"ab"[LINES.index(lam)]}) {lam:.0f} nm  --  '
                     f'audit window {rms_window*100:.1f}% RMS, r = {corr:.2f}',
                     fontsize=11)
        ax.text(0.02, 0.03,
                f'full 0-120 GPa axis: {rms_full*100:.1f}% RMS\n'
                f'(scored over {pressure[scorable].min():.0f}-'
                f'{pressure[scorable].max():.0f} GPa)',
                transform=ax.transAxes, fontsize=8.5, va='bottom',
                bbox=dict(fc='white', ec='0.7', alpha=0.85, pad=3.0))
        ax.set_xlabel('pressure [GPa]')
        ax.set_ylabel(r'$\sigma_{\rm abs}$, scaled to the published curve')
        ax.set_xlim(0.0, 120.0)
        ax.legend(fontsize=8.5,
                  loc='upper right' if lam == 532.0 else 'upper left')
        ax.grid(alpha=0.25)

    ax = axes[1][0]
    pooled_window = []
    for lam in LINES:
        pressure, inside, scorable, recon_full, reference = residuals[lam]
        shown = scorable
        full = (recon_full[shown] - reference[shown]) / reference[shown]
        ax.plot(pressure[shown], full * 100.0, color=COLOUR[lam], lw=1.0,
                alpha=0.55)
        band = inside & scorable
        ax.plot(pressure[band],
                ((recon_full[band] - reference[band]) / reference[band]) * 100.0,
                color=COLOUR[lam], lw=2.2, label=f'{lam:.0f} nm')
        pooled_window.append(
            (recon_full[band] - reference[band]) / reference[band])
    pooled = float(np.sqrt(np.mean(np.concatenate(pooled_window) ** 2)))
    ax.axhspan(-100.0 * HO_MAX_FRACTIONAL_RMS, 100.0 * HO_MAX_FRACTIONAL_RMS,
               color='0.80', alpha=0.45, zorder=0,
               label=f'audit tolerance ({HO_MAX_FRACTIONAL_RMS*100:.0f}% RMS)')
    ax.axhline(0.0, color='k', lw=0.8)
    ax.set_title('(c) fractional residual  --  audit windows (thick) '
                 f'pooled {pooled*100:.1f}% RMS', fontsize=11)
    ax.set_xlabel('pressure [GPa]')
    ax.set_ylabel('(reconstruction - published) / published  [%]')
    ax.set_xlim(0.0, 120.0)
    ax.set_ylim(-16.0, 16.0)
    ax.legend(fontsize=8.5, loc='upper right')
    ax.grid(alpha=0.25)

    ax = axes[1][1]
    expt_rms, floor = {}, 1.0
    for lam in LINES:
        pressure, values = data['expt%d' % lam]
        norm = values.max()
        prediction, residual = scaled(
            predict_absorption(published, lam, pressure), values)
        expt_rms[lam] = float(np.sqrt(np.mean(residual ** 2)))
        ax.plot(pressure, values / norm, 'o', color=COLOUR[lam], ms=4.0,
                mfc='none', label=f'{lam:.0f} nm, measured PL')
        ax.plot(pressure, prediction / norm, color=COLOUR[lam], lw=1.8,
                label=f'{lam:.0f} nm, reconstruction')
        floor = min(floor, float(min(prediction.min(), values.min())) / norm)
    ax.set_title('(d) different observable: measured PL  --  '
                 f'{expt_rms[532.0]*100:.0f}% / {expt_rms[457.0]*100:.0f}% RMS',
                 fontsize=11)
    ax.set_xlabel('pressure [GPa]')
    ax.set_ylabel('normalised to each series maximum')
    ax.set_xlim(0.0, 120.0)
    ax.set_ylim(floor - 0.05, 1.06)
    ax.legend(fontsize=8.5, loc='lower left')
    ax.grid(alpha=0.25)

    fig.suptitle('Cross-figure reproduction: the kernel reconstructed from '
                 'Ho Fig. 1(e), tested against Ho Fig. 5(b)', fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT, dpi=180)
    print(f'wrote {OUT}')
    print('  published axis: 0-120 GPa for both calculated curves')
    for lam in LINES:
        lo, hi, rms_w, rms_f, corr, s_lo, s_hi = summary[lam]
        print(f'  {lam:.0f} nm  audit window {lo:.1f}-{hi:.1f} GPa: '
              f'RMS {rms_w*100:.1f}%, r = {corr:.3f}')
        print(f'  {lam:.0f} nm  full axis, scored {s_lo:.0f}-{s_hi:.0f} GPa '
              f'(reference above {RESIDUAL_FLOOR:.0%} of max): '
              f'RMS {rms_f*100:.1f}%')
    print(f'  pooled audit-window RMS: {pooled*100:.1f}%')
    for lam in LINES:
        print(f'  {lam:.0f} nm measured PL: RMS {expt_rms[lam]*100:.1f}% '
              '(separate observable, not pooled)')


if __name__ == '__main__':
    main()
