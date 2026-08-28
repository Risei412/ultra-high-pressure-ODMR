"""Figure for the cross-figure reproduction of Ho et al. (2026) Fig. 5(b).

The freeze claims that the kernel reconstructed from Fig. 1(e) reproduces the
independently published Fig. 5(b) absorption curves.  `repro_yield.py` checks
that numerically but draws nothing, so the claim had no figure.  This is it.

(a),(b) the reconstruction against Ho's calculated absorption at the two laser
        lines, with the v1 phenomenological model on the same axes to show what
        a failing curve looks like;
(c)     the fractional residual of (a),(b) -- the 1.0% pooled agreement;
(d)     the measured PL of the same figure, kept in its own panel because it is
        a different observable: a detection passband may enter it, so it must
        not be pooled with the absorption comparison.

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


def audit_window(data, lam):
    """The pressure span the audit uses: Ho's curve, clipped to the data."""
    pressure, values = data['theory%d_ho' % lam]
    measured = data['expt%d' % lam][0]
    mask = (pressure >= measured.min()) & (pressure <= measured.max())
    return pressure[mask], values[mask]


def scaled(prediction, reference):
    """Prediction after the one allowed multiplicative scale."""
    residual, scale = _scaled_residual_from_prediction(prediction, reference)
    return scale * np.asarray(prediction, float), residual


def main():
    data = load()
    published = HoPublishedSpectrumModel()
    frozen = NVModel(T=90.0)

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.0))
    residuals = {}

    for ax, lam in zip(axes[0], LINES):
        pressure, reference = audit_window(data, lam)
        norm = reference.max()

        recon, residual = scaled(
            predict_absorption(published, lam, pressure), reference)
        v1, _ = scaled(predict_absorption(frozen, lam, pressure), reference)
        residuals[lam] = (pressure, residual)

        ax.plot(pressure, reference / norm, color=REF_COLOUR, lw=4.0,
                alpha=0.35, label='Ho Fig. 5(b), published')
        ax.plot(pressure, recon / norm, color=COLOUR[lam], lw=1.8,
                label='Fig. 1(e) reconstruction')
        ax.plot(pressure, v1 / norm, color='tab:red', lw=1.2, ls='--',
                alpha=0.8, label='v1 phenomenological model')

        rms = float(np.sqrt(np.mean(residual ** 2)))
        corr = float(np.corrcoef(recon, reference)[0, 1])
        ax.set_title(f'({"ab"[LINES.index(lam)]}) {lam:.0f} nm  '
                     f'--  RMS {rms*100:.1f}%, r = {corr:.2f}', fontsize=11)
        ax.set_xlabel('pressure [GPa]')
        ax.set_ylabel(r'$\sigma_{\rm abs}$, scaled to the published curve')
        ax.legend(fontsize=8.5, loc='best')
        ax.grid(alpha=0.25)

    ax = axes[1][0]
    for lam in LINES:
        pressure, residual = residuals[lam]
        ax.plot(pressure, residual * 100.0, color=COLOUR[lam], lw=1.6,
                label=f'{lam:.0f} nm')
    pooled = float(np.sqrt(np.mean(np.concatenate(
        [residuals[lam][1] for lam in LINES]) ** 2)))
    ax.axhspan(-100.0 * HO_MAX_FRACTIONAL_RMS, 100.0 * HO_MAX_FRACTIONAL_RMS,
               color='0.88', zorder=0,
               label=f'audit tolerance ({HO_MAX_FRACTIONAL_RMS*100:.0f}% RMS)')
    ax.axhline(0.0, color='k', lw=0.8)
    ax.set_title(f'(c) fractional residual  --  pooled RMS {pooled*100:.1f}%',
                 fontsize=11)
    ax.set_xlabel('pressure [GPa]')
    ax.set_ylabel('(reconstruction - published) / published  [%]')
    ax.set_ylim(-13.0, 13.0)
    ax.legend(fontsize=8.5, loc='best')
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
    ax.set_ylim(floor - 0.05, 1.06)
    ax.legend(fontsize=8.5, loc='best')
    ax.grid(alpha=0.25)

    fig.suptitle('Cross-figure reproduction: the kernel reconstructed from '
                 'Ho Fig. 1(e), tested against Ho Fig. 5(b)', fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT, dpi=180)
    print(f'wrote {OUT}')
    for lam in LINES:
        pressure, residual = residuals[lam]
        print(f'  {lam:.0f} nm absorption: window '
              f'{pressure.min():.1f}-{pressure.max():.1f} GPa, '
              f'RMS {np.sqrt(np.mean(residual**2))*100:.1f}%')
    print(f'  pooled absorption RMS: {pooled*100:.1f}%')
    for lam in LINES:
        print(f'  {lam:.0f} nm measured PL: RMS {expt_rms[lam]*100:.1f}% '
              '(separate observable, not pooled)')


if __name__ == '__main__':
    main()
