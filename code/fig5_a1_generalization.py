"""Figures for the numerical execution of Addendum A1.

(a) the reconstructed kernel A(lambda) at 120 GPa, its multimodal structure,
    the frozen 5 % band, and where the planned laser lines actually fall;
(b) the exactly degenerate doublet opening as power rises, against the
    single-peaked PL spectrum;
(c) P2 tested against exact argmax, including the jump to the ZPL;
(d) the P4 mechanism decomposition.

Writes `a1_generalization_120GPa.png`.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from theory_a1_generalization import (
    DATA_WINDOW, MAIN_BAND, Kernel, MediatedResponse, eta_at,
    section2_split_formula, sensitivity_optima,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'a1_generalization_120GPa.png')
PLANNED = {457.0: 'tab:green', 473.0: 'tab:orange', 488.0: 'tab:red'}


def panel_kernel(ax, kernel):
    grid = np.arange(*DATA_WINDOW, 0.02)
    ax.plot(grid, kernel.a(grid), color='k', lw=1.4, label=r'$A(\lambda)$')
    lo, hi = kernel.tolerance_band(1.05)
    ax.axvspan(lo, hi, color='tab:blue', alpha=0.15,
               label=f'5 % band [{lo:.0f}, {hi:.0f}] nm')
    for lam, value in kernel.local_maxima():
        ax.plot([lam], [value], 'o', ms=5, color='tab:blue')
        ax.annotate(f'{lam:.1f}', (lam, value), textcoords='offset points',
                    xytext=(0, 7), ha='center', fontsize=7.5)
    for lam, colour in PLANNED.items():
        ax.axvline(lam, color=colour, ls='--', lw=1.1)
        ax.annotate(f'{lam:.0f} nm', (lam, 0.04), rotation=90, fontsize=7.5,
                    color=colour, ha='right', va='bottom')
    ax.axvline(475.51, color='tab:purple', ls=':', lw=1.4)
    ax.annotate('v1 as frozen\n475.5 nm\n(superseded, E4:\nrepaired v1 gives\n439.1 nm)',
                (478.0, 1.06), fontsize=7.0,
                color='tab:purple', ha='left', va='top')
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel(r'$A/A_{\max}$')
    ax.set_title('(a) reconstructed kernel is multimodal', fontsize=9.5)
    ax.legend(fontsize=7.5, loc='lower left')
    ax.set_xlim(*DATA_WINDOW)
    ax.set_ylim(0.0, 1.12)


def panel_doublet(ax, kernel):
    response = MediatedResponse(gamma_contrast=1.0)
    star = response.gamma_star()
    grid = np.arange(*MAIN_BAND, 0.05)
    for ratio, colour in ((0.98, 'tab:blue'), (0.9, 'tab:green'),
                          (0.7, 'tab:red')):
        gamma_max = star / ratio
        eta = eta_at(kernel, response, grid, gamma_max)
        ax.plot(grid, eta / eta.min(), color=colour, lw=1.3,
                label=fr'$I/I_c={1.0 / ratio:.2f}$')
        for member in sensitivity_optima(kernel, response, gamma_max,
                                         MAIN_BAND)['optima']:
            ax.plot([member], [1.0], 'v', color=colour, ms=6)
    ax.axvline(kernel.lam_abs, color='k', lw=0.7, alpha=0.4)
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel(r'$\eta/\eta_{\min}$')
    ax.set_title('(b) sensitivity splits while PL stays single-peaked',
                 fontsize=9.5)
    ax.set_ylim(0.995, 1.05)
    ax.legend(fontsize=7.5, loc='upper left')

    # The PL spectrum lives on a far larger scale; give it its own axis so the
    # doublet in eta stays legible.
    twin = ax.twinx()
    rate = response.rate(star / 0.7 * kernel.a(grid))
    twin.plot(grid, rate / rate.max(), color='k', ls='--', lw=1.1)
    twin.set_ylabel(r'$R/R_{\max}$ (dashed)', fontsize=8.5)
    twin.set_ylim(0.0, 1.05)
    twin.tick_params(labelsize=8)


def panel_p2(ax):
    result = section2_split_formula()
    rows = sorted(result['rows'], key=lambda r: r['l_G_per_nm'])
    l_g = np.array([r['l_G_per_nm'] for r in rows]) * 100.0
    predicted = np.array([r['predicted_shift_nm'] for r in rows])
    exact = np.array([r['local_shift_nm'] for r in rows])
    jumped = np.array([r['jumped_to_zpl'] for r in rows])
    ax.plot(l_g, predicted, 'k--', lw=1.2, label=r'P2: $2\ell_G/\kappa_R$')
    ax.plot(l_g[~jumped], exact[~jumped], 'o-', color='tab:blue', ms=4,
            lw=1.2, label='exact argmax (main band)')
    ax.plot(l_g[jumped], exact[jumped], 'o', color='tab:red', ms=5,
            label='global optimum at ZPL')
    ax.axvline(result['zpl_jump_l_G'] * 100.0, color='tab:red', ls=':', lw=1.2)
    ax.annotate(f"ZPL jump\n{result['zpl_jump_l_G'] * 100:.2f} %/nm",
                (result['zpl_jump_l_G'] * 100.0, -35), fontsize=7.5,
                color='tab:red', ha='left')
    ax.axhspan(-15.6, 15.6, color='tab:blue', alpha=0.12, label='5 % band')
    ax.set_xlabel(r'$\ell_G = \mathrm{d}\ln(C/\Delta\nu)/\mathrm{d}\lambda$  [%/nm]')
    ax.set_ylabel(r'$\lambda_\eta-\lambda_{\rm PL}$  [nm]')
    ax.set_title('(c) P2 against exact argmax', fontsize=9.5)
    ax.legend(fontsize=7.5, loc='upper left')


def panel_mechanisms(ax):
    gamma = np.logspace(-2, 2, 400)
    # The curves coincide in pairs to machine precision, so they are drawn with
    # offset dashes: that degeneracy is the result, not a plotting artefact.
    cases = (
        ('rate saturation alone', MediatedResponse(gamma_sat=1.0),
         'tab:blue', '-', 3.2),
        ('power broadening alone', MediatedResponse(gamma_width=1.0),
         'tab:green', (0, (4, 4)), 1.5),
        ('saturation + broadening',
         MediatedResponse(gamma_sat=1.0, gamma_width=1.0), 'tab:orange', '-', 3.2),
        ('contrast collapse alone',
         MediatedResponse(gamma_contrast=1.0), 'tab:red', (0, (4, 4)), 1.5),
    )
    for label, response, colour, style, width in cases:
        ax.semilogx(gamma, response.dlnphi_dlngamma(gamma), color=colour,
                    lw=width, ls=style, label=label)
        star = response.gamma_star()
        if np.isfinite(star):
            ax.plot([star], [0.0], 'o', color=colour, ms=6)
    ax.axhline(0.0, color='k', lw=0.8)
    ax.annotate('the two pairs coincide exactly:\n$\\eta$ cannot tell them apart',
                (0.012, 0.45), fontsize=7.5, color='0.25')
    ax.set_xlabel(r'$\Gamma_p$ [half-scale units]')
    ax.set_ylabel(r'$\mathrm{d}\ln\Phi/\mathrm{d}\ln\Gamma_p$')
    ax.set_title('(d) splitting mechanisms are degenerate in $\\eta$', fontsize=9.5)
    ax.legend(fontsize=7.5, loc='lower left')
    ax.set_ylim(-1.2, 1.2)


def main():
    kernel = Kernel()
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6))
    panel_kernel(axes[0, 0], kernel)
    panel_doublet(axes[0, 1], kernel)
    panel_p2(axes[1, 0])
    panel_mechanisms(axes[1, 1])
    fig.suptitle('Addendum A1 executed against the frozen Ho kernel, 120 GPa',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT, dpi=170)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
