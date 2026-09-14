"""Concept figure for the esa article on the branching theorem (post 16).

Reproduces the three stages of Addendum A1's ordering proposition P5:

    (a) low power        Phi and the PL hill peak at the same wavelength
    (b) intermediate     Phi has already split in two; the PL hill is still
                         single-peaked, and its summit is now the locally
                         *worst* choice between the two optima
    (c) high power       the PL itself saturates into a plateau

This is deliberately a **concept** figure, so the absorption band is the
Gaussian approximation the article itself uses,

    a(lambda) = exp(-kappa (lambda - lambda_abs)^2 / 2),
    kappa = 8e-4 nm^-2,  lambda_abs = 440.65 nm,

not the reconstructed Ho kernel.  That keeps it consistent with the article's
own Table 2 and Table 5 -- panel (b) is drawn at a*/a_max = 0.9, whose
predicted separation is 32 nm -- and keeps the multimodal structure and the
clipped zero-phonon line of the real reconstruction out of a diagram whose job
is to show a mechanism.  For the real kernel see `fig5_a1_generalization.py`
(panel b) and `fig6_a2_multiplicity.py`.

Writes `esa_fig0_split_concept.png`.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

matplotlib.rcParams['font.family'] = 'IPAGothic'
matplotlib.rcParams['axes.unicode_minus'] = False

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'esa_fig0_split_concept.png')

LAMBDA_ABS = 440.65      # nm, the frozen optical-limit optimum
KAPPA = 8.0e-4           # nm^-2, the article's Table 2 low-power curvature

# Response half-scales, in units of the pump rate.  Rate saturation is set well
# below the contrast and linewidth scales so that panel (c) shows the PL
# plateau the caption calls for while the doublet stays on the plot.
GAMMA_SAT, GAMMA_CONTRAST, GAMMA_WIDTH = 0.15, 1.2, 1.2

PL_COLOUR, PHI_COLOUR = '0.55', 'tab:red'


def band(lam):
    """Gaussian absorption proxy, normalised to unity at the band maximum."""
    return np.exp(-KAPPA * (np.asarray(lam, float) - LAMBDA_ABS) ** 2 / 2.0)


def rate(gamma_p):
    return gamma_p / (1.0 + gamma_p / GAMMA_SAT)


def contrast(gamma_p):
    return 1.0 / (1.0 + gamma_p / GAMMA_CONTRAST)


def linewidth(gamma_p):
    return np.sqrt(1.0 + gamma_p / GAMMA_WIDTH)


def phi(gamma_p):
    """Phi = (C/dnu)^2 R = 1/eta^2."""
    gamma_p = np.asarray(gamma_p, float)
    return (contrast(gamma_p) / linewidth(gamma_p)) ** 2 * rate(gamma_p)


def dlnphi(gamma_p):
    rho_s = gamma_p / GAMMA_SAT
    rho_c = gamma_p / GAMMA_CONTRAST
    rho_w = gamma_p / GAMMA_WIDTH
    return (1.0 - rho_s / (1.0 + rho_s) - 2.0 * rho_c / (1.0 + rho_c)
            - rho_w / (1.0 + rho_w))


def gamma_star():
    """Pump rate at the interior maximum of Phi (solved in log space)."""
    root = brentq(lambda x: dlnphi(np.exp(x)), np.log(1e-6), np.log(1e6))
    return float(np.exp(root))


def doublet(level):
    """Wavelengths where a(lambda) = level, from the Gaussian band."""
    half = np.sqrt(2.0 * np.log(1.0 / level) / KAPPA)
    return LAMBDA_ABS - half, LAMBDA_ABS + half


def panel(ax, lam, gamma_max, title, level=None, annotate=None):
    gamma_p = gamma_max * band(lam)
    pl = rate(gamma_p)
    sensitivity = phi(gamma_p)

    # unsaturated reference: what R would look like with no rate saturation,
    # so that the flattening in panel (c) is visible despite per-panel scaling
    reference = band(lam)
    ax.plot(lam, reference, color=PL_COLOUR, lw=1.0, ls='--', alpha=0.75,
            label=r'未飽和の $R$（参照）')
    ax.plot(lam, pl / pl.max(), color=PL_COLOUR, lw=2.6,
            label=r'発光 $R(\lambda)$')
    ax.plot(lam, sensitivity / sensitivity.max(), color=PHI_COLOUR, lw=2.2,
            label=r'$\Phi=(C/\Delta\nu)^2R$')
    ax.axvline(LAMBDA_ABS, color='k', lw=0.8, ls=':', alpha=0.55)

    if level is None:
        ax.plot([LAMBDA_ABS], [1.0], 'v', color=PHI_COLOUR, ms=9, zorder=5)
    else:
        low, high = doublet(level)
        for member in (low, high):
            ax.plot([member], [1.0], 'v', color=PHI_COLOUR, ms=9, zorder=5)
        ax.annotate('', xy=(low, 1.07), xytext=(high, 1.07),
                    arrowprops=dict(arrowstyle='<->', color=PHI_COLOUR, lw=1.2))
        ax.text(LAMBDA_ABS, 1.10, f'{high - low:.0f} nm', ha='center',
                fontsize=9, color=PHI_COLOUR)
        # the PL summit is now the worst choice between the two optima
        centre = float(phi(gamma_max * band(LAMBDA_ABS)) / sensitivity.max())
        ax.plot([LAMBDA_ABS], [centre], 'o', color=PHI_COLOUR, ms=6,
                mfc='white', mew=1.6, zorder=5)

    if annotate:
        ax.annotate(annotate[0], xy=annotate[1], xytext=annotate[2],
                    fontsize=9, ha='center',
                    arrowprops=dict(arrowstyle='->', color='0.3', lw=1.0))

    ax.set_title(title, fontsize=11)
    ax.set_xlabel('励起波長 [nm]')
    ax.set_xlim(lam[0], lam[-1])
    ax.set_ylim(0.0, 1.22)


def main():
    lam = np.linspace(380.0, 501.0, 1400)
    star = gamma_star()

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5), sharey=True)

    panel(axes[0], lam, 0.15 * star, '(a) 低パワー — 一致',
          annotate=('発光の山と感度の山が\n同じ位置',
                    (LAMBDA_ABS, 1.0), (LAMBDA_ABS, 0.42)))
    panel(axes[1], lam, star / 0.9, '(b) 中間パワー — 感度だけが割れる',
          level=0.9,
          annotate=('発光極大は\n2 つの最適点のあいだの\n局所的な最悪点',
                    (LAMBDA_ABS, 0.965), (LAMBDA_ABS, 0.36)))
    panel(axes[2], lam, star / 0.53, '(c) 高パワー — 発光も頭打ち',
          level=0.53,
          annotate=('実線が破線より平ら\n= 発光が飽和',
                    (LAMBDA_ABS - 22, 0.93), (LAMBDA_ABS - 40, 0.45)))

    axes[0].set_ylabel('各パネルの最大値で規格化')
    axes[0].legend(fontsize=9, loc='lower left')
    fig.suptitle('発光最適点と感度最適点が分かれる過程'
                 '（概念図・ガウス近似、$\\kappa_A=8\\times10^{-4}\\,$nm$^{-2}$）',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT, dpi=180)
    print(f'wrote {OUT}')
    print(f'  Gamma_p* = {star:.4f}')
    for name, level in (('(b)', 0.9), ('(c)', 0.53)):
        low, high = doublet(level)
        print(f'  {name} a*/a_max = {level}: '
              f'{low:.1f} / {high:.1f} nm, separation {high - low:.1f} nm')


if __name__ == '__main__':
    main()
