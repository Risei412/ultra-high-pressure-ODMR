"""Figure for the cross-figure reproduction of Ho et al. (2026) Fig. 5(b).

The freeze claims that the kernel reconstructed from Fig. 1(e) reproduces the
independently published Fig. 5(b) absorption curves.  `repro_yield.py` checks
that numerically but draws nothing, so the claim had no figure.  This is it.

Panels (a),(b) are one per laser line: Ho's calculated absorption against the
reconstruction, over the full published 0-120 GPa axis.  Panel (c) is there
because (a),(b) alone would now mislead -- see below.

What this figure does NOT show is our own physics predicting Ho's result.
The reconstruction is Ho's own Fig. 1(e) spectra, digitised and interpolated
(`ho_spectrum_model.py` -- no DFT, no Jahn-Teller), so the agreement is a
consistency check between two figures of the same paper: it establishes that
the extraction feeding the whole v3 chain is faithful, and nothing about
whether we could have predicted the curve.  The one curve here that IS our
own physics is v1, drawn on the same axes to show what a failing curve looks
like -- at 532 nm it puts the maximum at 38 GPa instead of 17 and
anticorrelates (r = -0.51) inside the audit window.  That contrast is the
point of drawing it: it shows the agreement is not the kind any smooth curve
would achieve, while making clear which of the two curves is ours.

To be precise about what v1 is: it is anchored to Ho's *reported scalars*
(dE_ZPL(120 GPa) = 0.400 eV, S_abs 3.08 -> 4.61) but no parameter of it was
ever fitted to the Fig. 5(b) curves it is compared against here.

The anticorrelation is not, however, a structural failure of the model.
`v1_diagnosis.py` shows that two of v1's input constants can be derived from
Ho's own Fig. 1(b),(e) -- the effective phonon energy is 101.1 meV, not the
65 meV of Kehayias et al., and the kernel's ZPL shift is 0.464 eV, not the
0.400 eV bound v1 is anchored to -- and that correcting both, with nothing
fitted to Fig. 5(b), takes v1 to pooled 1.0%, r = 1.00 and peaks at 17/17 and
89/88 GPa.  It then passes both gates, scoring exactly what the reconstruction
scores on this test.

So the red curve is not evidence that a single-mode model cannot reproduce
Ho.  It is v1 as frozen, which is what the freeze compares against, and what
it shows is the size of a two-constant error.  The reason v3 still takes the
kernel from outside is elsewhere: scanned in wavelength at 120 GPa the
corrected model has one local maximum where the kernel has four, and every
Addendum A2 result lives on those four.

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
from nv_model import NVModel, nm2eV
from v1_diagnosis import phonon_energy, zpl_shift
from repro_yield import (
    EXPT_MAX_FRACTIONAL_RMS, EXPT_MAX_PEAK_ERROR_GPA, HO_MAX_FRACTIONAL_RMS,
    HO_MAX_PEAK_ERROR_GPA, LINES, _scaled_residual_from_prediction,
    compare_experiment, compare_ho_theory, load, predict_absorption,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'repro_ho_fig5b.png')

COLOUR = {532.0: 'tab:green', 457.0: 'tab:blue'}
REF_COLOUR = '0.35'
WINDOW_COLOUR = '0.90'
# Fraction of its own maximum below which the digitised reference is too close
# to zero for a fractional residual to mean anything.
RESIDUAL_FLOOR = 0.05

# Panel (c): the wavelength scan that says why the kernel is taken from
# outside, at the pressure the whole v3 chain is quoted at.
STRUCTURE_PRESSURE_GPA = 120.0
STRUCTURE_WINDOW_NM = (402.0, 517.0)
V1_COLOUR, V1_FIXED_COLOUR = 'tab:red', 'tab:purple'


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

    _, _, hw = phonon_energy()
    dE120 = zpl_shift()
    repaired = NVModel(T=90.0, hw=hw, dE120=dE120)

    fig = plt.figure(figsize=(16.4, 5.8))
    grid = fig.add_gridspec(2, 3, height_ratios=(3.0, 1.0), hspace=0.08,
                            width_ratios=(1.0, 1.0, 1.05))
    axes = [[fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])],
            [None, None]]
    axes[1][0] = fig.add_subplot(grid[1, 0], sharex=axes[0][0])
    axes[1][1] = fig.add_subplot(grid[1, 1], sharex=axes[0][1])
    structure_ax = fig.add_subplot(grid[:, 2])
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
        fixed_raw = predict_absorption(repaired, lam, pressure)

        # the audit's own scale, fitted inside the audit window only
        recon_window, residual_window = scaled(prediction[inside],
                                               reference[inside])
        recon = recon_window.max() / prediction[inside].max() * prediction
        v1 = (scaled(v1_raw[inside], reference[inside])[0].max()
              / v1_raw[inside].max()) * v1_raw
        fixed_window, fixed_residual = scaled(fixed_raw[inside],
                                              reference[inside])
        fixed = fixed_window.max() / fixed_raw[inside].max() * fixed_raw
        residual = (recon[scorable] - reference[scorable]) / reference[scorable]

        top.axvspan(lo, hi, color=WINDOW_COLOUR, zorder=0,
                    label='audit window (measured range)')
        top.plot(pressure, reference / norm, color=REF_COLOUR, lw=4.0,
                 alpha=0.35, label='Ho Fig. 5(b), published')
        top.plot(pressure, recon / norm, color=COLOUR[lam], lw=1.8,
                 label='Fig. 1(e) reconstruction')
        top.plot(pressure, v1 / norm, color=V1_COLOUR, lw=1.0, ls='--',
                 alpha=0.55,
                 label='OURS, superseded: v1 as frozen (E4)')
        top.plot(pressure, fixed / norm, color=V1_FIXED_COLOUR, lw=2.0,
                 ls=(0, (1.2, 1.2)),
                 label='OURS, current: v1 + Ho-derived constants')

        # the peak gate: both markers coincide, so only one line is visible
        peak_ref = float(pressure[inside][np.argmax(reference[inside])])
        peak_recon = float(pressure[inside][np.argmax(recon[inside])])
        top.axvline(peak_ref, color=REF_COLOUR, lw=1.0, alpha=0.6)
        top.axvline(peak_recon, color=COLOUR[lam], lw=1.0, ls=':', alpha=0.9)
        top.annotate(f'peak {peak_recon:.0f}/{peak_ref:.0f} GPa\n'
                     f'(gate: within {HO_MAX_PEAK_ERROR_GPA:.0f} GPa)',
                     xy=(peak_ref, 0.30), xycoords=('data', 'axes fraction'),
                     ha='center', va='center', fontsize=8.0, color='0.3',
                     bbox=dict(fc='white', ec='0.75', alpha=0.85, pad=2.5))

        rms_window = float(np.sqrt(np.mean(residual_window ** 2)))
        rms_fixed = float(np.sqrt(np.mean(fixed_residual ** 2)))
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
        top.text(0.03 if lam == 532.0 else 0.56, 0.06,
                 f'full 0-120 GPa axis: {rms_full*100:.1f}% RMS\n'
                 f'(scored over {pressure[scorable].min():.0f}-'
                 f'{pressure[scorable].max():.0f} GPa)\n'
                 f'v1 repaired, audit window: {rms_fixed*100:.1f}% RMS',
                 transform=top.transAxes, fontsize=8.5, va='bottom',
                 bbox=dict(fc='white', ec='0.7', alpha=0.85, pad=3.0))
        top.set_ylabel(r'$\sigma_{\rm abs}$, scaled to the published curve',
                       fontsize=9.5)
        top.legend(fontsize=7.8, framealpha=0.92,
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
        strip.plot(pressure[inside], fixed_residual * 100.0,
                   color=V1_FIXED_COLOUR, lw=1.4, ls=(0, (1, 1.4)))
        strip.set_xlim(0.0, 120.0)
        strip.set_ylim(-16.0, 16.0)
        strip.set_yticks((-10, 0, 10))
        strip.set_xlabel('pressure [GPa]')
        strip.set_ylabel('(curve - Ho) / Ho  [%]', fontsize=8.5)
        strip.grid(alpha=0.25)

    pooled = float(np.sqrt(np.mean(np.concatenate(pooled_window) ** 2)))
    for column in range(2):
        axes[1][column].text(
            0.015, 0.93,
            'residual of the panel above, magnified\n'
            + ('solid = reconstruction, dotted = ours\n'
               'band = 10% gate, thick = audit window'
               if column == 0 else 'both curves are inside the 10% gate'),
            transform=axes[1][column].transAxes, fontsize=7.0, va='top')

    # ---- (c) why the kernel is still taken from outside --------------------
    lam_nm = np.arange(STRUCTURE_WINDOW_NM[0], STRUCTURE_WINDOW_NM[1] + 0.025,
                       0.05)
    energy = nm2eV(lam_nm)
    kernel = np.asarray(published.sigma_abs(energy, STRUCTURE_PRESSURE_GPA),
                        float)
    ours = np.array([repaired.sigma_abs(value, STRUCTURE_PRESSURE_GPA)
                     for value in energy])
    kernel, ours = kernel / kernel.max(), ours / ours.max()
    interior = lambda y: np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1

    structure_ax.plot(lam_nm, kernel, color=REF_COLOUR, lw=2.4,
                      label='Ho Fig. 1(e) kernel')
    structure_ax.plot(lam_nm, ours, color=V1_FIXED_COLOUR, lw=1.6,
                      ls=(0, (1.2, 1.2)),
                      label='OURS, current: same constants')
    peaks = interior(kernel)
    structure_ax.plot(lam_nm[peaks], kernel[peaks], 'v', color=REF_COLOUR,
                      ms=8, zorder=5)
    for index in peaks:
        structure_ax.axvline(lam_nm[index], color=REF_COLOUR, lw=0.7,
                             alpha=0.35)
    structure_ax.plot(lam_nm[interior(ours)], ours[interior(ours)], 'v',
                      color=V1_FIXED_COLOUR, ms=8, zorder=5)
    structure_ax.set_title(f'(c) the same repaired model in WAVELENGTH at '
                           f'{STRUCTURE_PRESSURE_GPA:.0f} GPa\n'
                           f'{len(peaks)} local maxima against '
                           f'{len(interior(ours))}', fontsize=11)
    structure_ax.set_xlabel('excitation wavelength [nm]')
    structure_ax.set_ylabel('absorption, each normalised to its own maximum')
    structure_ax.set_xlim(*STRUCTURE_WINDOW_NM)
    structure_ax.set_ylim(0.0, 1.52)
    structure_ax.legend(fontsize=8.5, loc='center left')
    structure_ax.grid(alpha=0.25)
    structure_ax.text(
        0.03, 0.985,
        'passing (a),(b) does not make the model a substitute for the\n'
        'kernel: Addendum A2 -- the level sets, the ladder\n'
        r'$N=2\to4\to6\to4\to3\to5\to3$, the transition powers --'
        '\nlives on these maxima, and one mode cannot produce four.\n'
        'The 514 nm member is the ZPL, whose HEIGHT erratum E3 withdraws;\n'
        'three interior maxima survive that, and one is still not three.',
        transform=structure_ax.transAxes, fontsize=8.0, va='top',
        bbox=dict(fc='white', ec='0.75', alpha=0.9, pad=3.0))

    fig.suptitle('Cross-figure reproduction: the kernel reconstructed from Ho '
                 f'Fig. 1(e), tested against Ho Fig. 5(b)  --  pooled '
                 f'{pooled*100:.1f}% RMS', fontsize=12.5)
    fig.text(0.5, 0.925, 'the reconstruction is Ho\u2019s own published '
             'spectra interpolated, not an independent calculation; the two '
             'dashed/dotted curves are ours \u2014 read the purple one',
             ha='center', fontsize=9.0, color='0.35')
    fig.subplots_adjust(left=0.048, right=0.988, top=0.845, bottom=0.10,
                        wspace=0.21, hspace=0.08)
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

    # The gates the freeze actually quotes PASS/FAIL against.  They are numbers,
    # not curves, so they are printed rather than given a panel of their own.
    ho = compare_ho_theory(published, data)
    expt = compare_experiment(published, data, collection=True)
    print('\n  reproduction gates (ceiling in brackets)')
    print(f'  {"comparison":<28}{"shape RMS":>18}{"peak error":>20}')
    for label, audit, rms_max, peak_max in (
            ('vs Ho calculated abs.', ho, HO_MAX_FRACTIONAL_RMS,
             HO_MAX_PEAK_ERROR_GPA),
            ('vs measured PL', expt, EXPT_MAX_FRACTIONAL_RMS,
             EXPT_MAX_PEAK_ERROR_GPA)):
        for lam in LINES:
            metrics = audit[str(int(lam))]
            rms, peak = metrics['fractional_rms'], metrics['peak_error_GPa']
            print(f'  {label + f" {lam:.0f} nm":<28}'
                  f'{rms*100:8.1f}% [{rms_max*100:4.0f}%] '
                  f'{"PASS" if rms <= rms_max else "FAIL"}'
                  f'{peak:9.0f} [{peak_max:3.0f}] GPa '
                  f'{"PASS" if peak <= peak_max else "FAIL"}')
    print('  the measured-PL rows are a different observable -- a detection '
          'passband,\n  the emission quantum yield and charge-state '
          'conversion all enter them --\n  which is why they are gated '
          'separately and are not drawn above.')


if __name__ == '__main__':
    main()
