"""The headline result: where to put the laser at 120 GPa.

Panel (a) is the optical-limit answer -- the penalty a given excitation
wavelength pays against the best one, so the y axis is the quantity the
decision is actually made on rather than the absorption it derives from.  In
the optical limit every non-optical factor is wavelength independent, so

    eta(lambda) / eta_min = sqrt( A_max / A(lambda) ),

and the frozen tolerance band is the level set of that at 1.05.

Panel (b) is what stops (a) from being the whole answer.  Above the splitting
power the minimiser of eta is not a point but a LEVEL SET of A, so "the
optimal wavelength" becomes a set whose size is a step function of power
(Theorem M).  A figure that showed only (a) would licence a single-wavelength
experiment, which is exactly the prediction this theory contradicts.

Three things on (a) are worth saying out loud because they are counterintuitive
and each is load bearing:

* the band is NOT symmetric -- 14.2 nm to the blue, 17.3 nm to the red -- so
  a symmetric quoted tolerance is wrong on both sides;
* the zero-phonon line at 514.5 nm costs almost exactly what 473 nm costs,
  despite lying 74 nm away, because the kernel is multimodal;
* the 402 nm edge is where Ho's figure stops, not where the band stops (K1),
  so nothing here licences a claim about the deeper blue.

Writes `optimum_wavelength_120GPa.png`.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from theory_a1_generalization import (
    DATA_WINDOW, MAIN_BAND, Kernel, MediatedResponse, eta_at,
    sensitivity_optima,
)
from v1_diagnosis import phonon_energy, wavelength_structure, zpl_shift

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'optimum_wavelength_120GPa.png')

# Lines a high-pressure ODMR bench can actually produce, and the ZPL, which is
# not a laser line but is where the kernel's fourth maximum sits.
CANDIDATES = ((457.0, 'tab:green'), (473.0, 'tab:orange'), (488.0, 'tab:red'))
BAND_PENALTY = 1.05
PENALTY_CEILING = 2.0


def penalty(kernel, lam):
    """Optical-limit sensitivity penalty against the best wavelength."""
    return 1.0 / np.sqrt(kernel.a(lam))


def panel_answer(ax, kernel):
    grid = np.arange(DATA_WINDOW[0], DATA_WINDOW[1] + 0.005, 0.05)
    values = penalty(kernel, grid)
    blue, red = kernel.tolerance_band(BAND_PENALTY)

    ax.axvspan(blue, red, color='#cfe3f5', zorder=0,
               label=f'within {100*(BAND_PENALTY-1):.0f}% of the optimum: '
                     f'[{blue:.2f}, {red:.2f}] nm')
    ax.plot(grid, values, color='k', lw=2.0, zorder=3)
    ax.axhline(1.0, color='0.6', lw=0.8)
    ax.axhline(BAND_PENALTY, color='0.6', lw=0.8, ls=':')

    ax.axvline(kernel.lam_abs, color='tab:blue', lw=1.6)
    ax.annotate(f'{kernel.lam_abs:.2f} nm', (kernel.lam_abs, 1.005),
                fontsize=12, fontweight='bold', color='tab:blue',
                ha='center', va='bottom')

    for lam, colour in CANDIDATES:
        cost = float(penalty(kernel, lam))
        ax.plot([lam], [cost], 'o', color=colour, ms=7, zorder=4)
        ax.annotate(f'{lam:.0f} nm\n+{100*(cost-1):.1f}%', (lam, cost),
                    textcoords='offset points', xytext=(0, 9), fontsize=8.5,
                    color=colour, ha='center', va='bottom')

    zpl = max(kernel.local_maxima(), key=lambda pair: pair[0])[0]
    zpl_cost = float(penalty(kernel, zpl))
    ax.plot([zpl], [zpl_cost], 'D', color='tab:purple', ms=6, zorder=4)
    ax.annotate(f'ZPL {zpl:.1f} nm\n+{100*(zpl_cost-1):.1f}%', (zpl, zpl_cost),
                textcoords='offset points', xytext=(-6, 10), fontsize=8.5,
                color='tab:purple', ha='right', va='bottom')

    # robustness: the repaired v1 envelope lands within 1.5 nm
    structure = wavelength_structure(phonon_energy()[2], zpl_shift())
    ours = structure['lambda_opt_model']
    ax.plot([ours], [float(penalty(kernel, ours))], '|', color='tab:purple',
            ms=14, mew=2.0, zorder=4)
    ax.annotate(f'our own repaired model (E4) puts the optimum at '
                f'{ours:.1f} nm,\n{abs(ours - kernel.lam_abs):.1f} nm away: '
                'the answer is not an artefact of the reconstruction',
                (0.02, 0.965), xycoords='axes fraction', fontsize=8.0,
                color='tab:purple', ha='left', va='top',
                bbox=dict(fc='white', ec='0.8', alpha=0.9, pad=3.0))

    ax.annotate('402 nm is where Ho’s figure stops,\n'
                'not where the band stops (K1)',
                (0.02, 0.80), xycoords='axes fraction', fontsize=7.5,
                color='0.35', ha='left', va='top')

    ax.set_xlim(*DATA_WINDOW)
    ax.set_ylim(0.98, PENALTY_CEILING)
    ax.set_xlabel('excitation wavelength [nm]')
    ax.set_ylabel(r'sensitivity penalty  $\eta(\lambda)\,/\,\eta_{\min}$'
                  '\n(1.00 = best possible; lower is better)')
    ax.set_title('(a) the optical-limit answer  --  '
                 f'{kernel.lam_abs:.2f} nm, with an asymmetric tolerance '
                 f'({kernel.lam_abs-blue:.1f} nm blue, '
                 f'{red-kernel.lam_abs:.1f} nm red)', fontsize=11)
    ax.legend(fontsize=8.5, loc='lower right')
    ax.grid(alpha=0.25)


def panel_level_set(ax, kernel):
    response = MediatedResponse(gamma_contrast=1.0)
    star = response.gamma_star()
    grid = np.arange(DATA_WINDOW[0], DATA_WINDOW[1] + 0.005, 0.02)
    counts = []
    for ratio, colour in ((0.90, 'tab:blue'), (0.60, 'tab:green'),
                          (0.50, 'tab:red')):
        gamma_max = star / ratio
        eta = eta_at(kernel, response, grid, gamma_max)
        result = sensitivity_optima(kernel, response, gamma_max, DATA_WINDOW)
        members = result['optima']
        counts.append((1.0 / ratio, len(members), result['truncated_blue']))
        ax.plot(grid, eta / eta.min(), color=colour, lw=1.4,
                label=fr'$I/I_c = {1.0/ratio:.2f}$   $N = {len(members)}$'
                      + ('  (one lost off 402 nm)'
                         if result['truncated_blue'] else ''))
        for member in members:
            ax.plot([member], [1.0], 'v', color=colour, ms=7, zorder=5)

    ax.axvline(kernel.lam_abs, color='tab:blue', lw=1.0, alpha=0.5)
    ax.set_xlim(*DATA_WINDOW)
    ax.set_ylim(0.997, 1.075)
    ax.set_xlabel('excitation wavelength [nm]')
    ax.set_ylabel(r'$\eta/\eta_{\min}$ at finite power')
    ax.set_title('(b) above $I_c$ the optimum is a SET, not a point',
                 fontsize=11)
    ax.legend(fontsize=8.5, loc='upper left')
    ax.grid(alpha=0.25)
    return counts


def main():
    kernel = Kernel()
    fig, axes = plt.subplots(1, 2, figsize=(15.2, 5.6),
                             gridspec_kw=dict(width_ratios=(1.65, 1.0)))
    panel_answer(axes[0], kernel)
    counts = panel_level_set(axes[1], kernel)
    counts = sorted(counts)

    fig.suptitle('Optimal excitation wavelength for NV ODMR sensitivity at '
                 '120 GPa', fontsize=13.5)
    fig.text(0.5, 0.925,
             'optical kernel reconstructed from Ho et al. (2026) Fig. 1(e); '
             'the penalty axis assumes only that non-optical factors are '
             'wavelength independent at low power',
             ha='center', fontsize=9.0, color='0.35')
    fig.text(
        0.5, 0.035,
        '(b) every marked wavelength is EXACTLY as good as the others: they '
        r'are the level set $\{\lambda: A(\lambda)=I_c/I\}$, whose size is '
        'a step function of power (Theorem M).\n'
        f'Raising the power by {counts[2][0]/counts[0][0]:.1f}x takes '
        f'N = {counts[0][1]} to {counts[1][1]} to {counts[2][1]} and moves the '
        'optimum onto the zero-phonon line, discontinuously.  So (a) alone '
        'would licence a single-wavelength experiment, which is the '
        'prediction this theory contradicts.',
        ha='center', va='bottom', fontsize=8.5)
    fig.subplots_adjust(left=0.058, right=0.99, top=0.845, bottom=0.19,
                        wspace=0.20)
    fig.savefig(OUT, dpi=180)

    blue, red = kernel.tolerance_band(BAND_PENALTY)
    print(f'wrote {OUT}')
    print(f'  optimum          {kernel.lam_abs:.2f} nm')
    print(f'  {100*(BAND_PENALTY-1):.0f}% band       [{blue:.2f}, {red:.2f}] nm'
          f'  (blue {kernel.lam_abs-blue:.2f}, red {red-kernel.lam_abs:.2f})')
    for lam, _ in CANDIDATES:
        print(f'  {lam:.0f} nm           +{100*(penalty(kernel, lam)-1):.1f}%')
    zpl = max(kernel.local_maxima(), key=lambda pair: pair[0])[0]
    print(f'  ZPL {zpl:.2f} nm     +{100*(penalty(kernel, zpl)-1):.1f}%')
    for ratio, count, truncated in counts:
        print(f'  I/I_c = {ratio:.2f}     {count} equally optimal wavelengths'
              + ('  (one lost off the 402 nm window edge)' if truncated else ''))


if __name__ == '__main__':
    main()
