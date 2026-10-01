"""Plot the current model's alpha-invariant 120 GPa optical baseline.

This is a handoff figure, not an orientation-resolved optical prediction.
Panel (a) shows that the axisymmetric stress tensor changes over
alpha = sigma_xx/sigma_zz = sigma_yy/sigma_zz. Panel (b) shows why the current
answer does not: every alpha uses the same Ho hydrostatic 120 GPa absorption
kernel because no calibrated anisotropic-stress optical kernel is available.

Run from ``code/``:

    python fig_alpha_invariance.py
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import fig_style
from nv_model_111 import NV111AxialModel


HERE = Path(__file__).resolve().parent
# This standalone handoff places generated figures beside ``codes/`` rather
# than assuming the parent repository's ``docs/experiment`` hierarchy.
OUT_DIR = HERE.parent / 'figures'
OUT_PNG = OUT_DIR / 'alpha_invariance_120gpa.png'
OUT_PDF = OUT_DIR / 'alpha_invariance_120gpa.pdf'

SIGMA_Z_GPA = 120.0
OPTICAL_PRESSURE_GPA = 120.0
ALPHA_MIN = 0.5
ALPHA_MAX = 0.7

BLACK = '#000000'
BLUE = '#0072B2'
GREEN = '#009E73'
VERMILLION = '#D55E00'
PURPLE = '#CC79A7'


def alpha_sweep(n=81):
    """Return the stress invariants and unchanged hydrostatic optical optimum."""
    alpha = np.linspace(ALPHA_MIN, ALPHA_MAX, n)
    sigma_xy = SIGMA_Z_GPA * alpha
    mean = SIGMA_Z_GPA * (1.0 + 2.0 * alpha) / 3.0
    axial_deviator = SIGMA_Z_GPA * (1.0 - alpha)
    numerical = np.array([
        NV111AxialModel(alpha=value).lambda_opt(OPTICAL_PRESSURE_GPA)
        for value in alpha
    ])
    reported = np.array([
        NV111AxialModel(alpha=value).reported_lambda_opt(
            OPTICAL_PRESSURE_GPA)
        for value in alpha
    ])
    return {
        'alpha': alpha,
        'sigma_xy_GPa': sigma_xy,
        'mean_GPa': mean,
        'axial_deviator_GPa': axial_deviator,
        'numerical_optimum_nm': numerical,
        'reported_optimum_nm': reported,
    }


def make_figure(data=None):
    data = data or alpha_sweep()
    alpha = data['alpha']

    fig_style.use()
    fig, axes = plt.subplots(1, 2, figsize=(fig_style.FULL_WIDTH_IN, 3.45))
    fig.subplots_adjust(left=0.095, right=0.975, bottom=0.285, top=0.91,
                        wspace=0.31)

    ax = axes[0]
    ax.plot(alpha, np.full_like(alpha, SIGMA_Z_GPA), color=BLACK, lw=1.3,
            label=r'$\sigma_{zz}$')
    ax.plot(alpha, data['sigma_xy_GPa'], color=BLUE, lw=1.3,
            label=r'$\sigma_{xx}=\sigma_{yy}$')
    ax.plot(alpha, data['mean_GPa'], color=GREEN, lw=1.3,
            label=r'$P_{\rm mean}={\rm tr}(\sigma)/3$')
    ax.plot(alpha, data['axial_deviator_GPa'], color=PURPLE, lw=1.3,
            label=r'$q=\sigma_{zz}-(\sigma_{xx}+\sigma_{yy})/2$')
    ax.set_xlim(ALPHA_MIN, ALPHA_MAX)
    ax.set_ylim(30.0, 125.0)
    ax.set_xlabel(r'Stress ratio $\alpha$')
    ax.set_ylabel('Stress (GPa)')
    ax.set_title('(a) Stress tensor changes', loc='left')
    ax.grid(alpha=0.22, lw=0.5)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2,
              frameon=True, framealpha=1.0, facecolor='white',
              edgecolor='0.8', fancybox=False, columnspacing=1.2)

    ax = axes[1]
    numerical = data['numerical_optimum_nm']
    ax.plot(alpha, numerical, color=VERMILLION, lw=1.4, marker='o',
            markevery=10, ms=3.5, mfc='white', mew=1.0,
            label='reconstructed-curve maximum')
    ax.axhline(441.0, color=BLACK, lw=0.8, ls=(0, (4, 2)),
               label=r'reported value $\approx441$ nm')
    ax.set_xlim(ALPHA_MIN, ALPHA_MAX)
    ax.set_ylim(438.5, 443.0)
    ax.set_xlabel(r'Stress ratio $\alpha$')
    ax.set_ylabel(r'Optimal excitation wavelength $\lambda_{\rm opt}$ (nm)')
    ax.set_title('(b) Hydrostatic optical baseline is unchanged', loc='left')
    ax.grid(alpha=0.22, lw=0.5)
    ax.legend(loc='upper right', frameon=True, framealpha=1.0,
              facecolor='white', edgecolor='0.8', fancybox=False)
    return fig


def main():
    data = alpha_sweep()
    fig = make_figure(data)
    print(fig_style.save(fig, OUT_PNG, OUT_PDF))
    print('\nalpha range: %.2f-%.2f' % (ALPHA_MIN, ALPHA_MAX))
    print('sigma_zz: %.1f GPa' % SIGMA_Z_GPA)
    print('hydrostatic optical-kernel pressure: %.1f GPa'
          % OPTICAL_PRESSURE_GPA)
    print('reconstructed-curve numerical maximum: %.2f-%.2f nm'
          % (np.min(data['numerical_optimum_nm']),
             np.max(data['numerical_optimum_nm'])))
    print('reported optimum: %.0f nm' % data['reported_optimum_nm'][0])
    print('max alpha-dependent numerical change: %.3g nm'
          % np.ptp(data['numerical_optimum_nm']))
    print('status: hydrostatic baseline, not an orientation-resolved prediction')


if __name__ == '__main__':
    main()
