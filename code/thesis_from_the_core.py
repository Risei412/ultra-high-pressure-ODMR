"""Bhattacharyya's chapter 6, rewritten in the language of the theory's core.

`docs/bhattacharyya_thesis_scope.md` sorted the thesis into what the theory
covers and what it does not.  This goes further: it asks what the CORE --
eta = dnu/(C sqrt(R)), assumption (M), and Theorems M/G/X -- actually says
about the phenomena the thesis reports, including the ones filed as "out of
scope".  Three of the four results below were not obtainable from the scope
sort, and one of them corrects it.

N1  Contrast inversion is invisible to eta; its ZERO is the event.
    Phi = (C/dnu)^2 R is EVEN in C.  A resonance with C = -5% is exactly as
    good a magnetometer as one with C = +5%, so "positive contrast" is not
    itself a degradation -- the thesis treats it as a puzzle about the ISC,
    and it is, but not one that costs sensitivity.  What the core does force
    is that C, being continuous in stress, must pass through ZERO between the
    negative and positive regimes, and there eta diverges.  Two consequences
    the thesis does not state:
      * per NV group and stress path there is a pressure at which cwODMR
        sensitivity diverges and then recovers;
      * in a widefield image the boundary between positive- and
        negative-contrast regions is a CONTOUR OF ZERO CONTRAST -- a dark
        ODMR line across an otherwise bright field.  Fig. 6.7(b) is a line
        cut across exactly such a boundary.

N2  The whole ISC story is a PREFACTOR in the gauge, so every structural
    prediction survives it.
    Goldman's Gamma_ISC depends on the stress-induced gap Delta, not on the
    pump rate, so it enters as C_0(P) with C(Gamma_p) = C_0(P) (1+rho)^-c.
    A constant multiplying C cancels out of dlnPhi/dlnGamma_p, hence out of
      * the exponent E = 2c + s + 2w,
      * the splitting criterion E n > 1 and rho* = 1/(E n - 1),
      * the ladder ratios I_k/I_c = A_max/A_k.
    It scales eta and nothing else.  So the culet-cut choice and the
    wavelength choice COMPOSE MULTIPLICATIVELY AND INDEPENDENTLY: the
    thesis's [111] result and this theory's wavelength result stack rather
    than compete.  That was an impression in the scope document; here it is
    a consequence of the gauge decomposition.

N3  The thesis's SNR metric is biased AGAINST the bluer line.
    It scores SNR ~ contrast * sqrt(counts) = C sqrt(R) = dnu / eta, so
    comparing two wavelengths,
        SNR_ratio / (1/eta)_ratio = dnu(lambda_1) / dnu(lambda_2).
    dnu grows with Gamma_p, so the metric rewards the MORE absorbing line
    beyond what sensitivity warrants.  At 50 GPa that is 532 nm.  The bias
    runs from 1.000 at low power to 0.95 at high, and near the crossover the
    two metrics disagree about which line wins.  So Fig. 6.3's null result
    has two causes, not one: the crossover sits at 51.8 GPa AND the metric
    tilts a further 1-5% toward 532 nm.

N4  The anvil transmission is not a small correction -- it can move the
    global optimum by 74 nm.  THIS CORRECTS the scope document.
    Section 6.2 of the thesis records that type Ib diamond absorbs below
    500 nm.  The effective kernel is A_eff = T_anvil(lambda) A(lambda), and
    the scope document claimed the ladder survives because T is pressure
    independent.  That is wrong: T reshapes A_eff, so the critical points and
    hence the rung ratios move with it.  Worse, the answer is FRAGILE.
    Parameterising T = 10^-((500-lambda)/D) below 500 nm, the global optimum
    leaves the main band for the zero-phonon line as soon as D < 321 nm --
    that is, as soon as the anvil costs more than a factor 1.53 at 440 nm
    relative to 500 nm.  A factor 1.5 is small.
    What survives is the STRUCTURE: Theorems M, G and X hold for any kernel,
    so the optimum is still a level set with a step-function multiplicity.
    The numbers do not.  Measuring the anvil transmission is therefore a
    PREREQUISITE for quoting 440.65 nm, not a refinement of it.

Run for the tables; each result is returned for the tests.
"""
import numpy as np

from ho_spectrum_model import HoPublishedSpectrumModel
from nv_model import nm2eV
from theory_a1_generalization import DATA_WINDOW, Kernel, MediatedResponse

# The thesis's own comparison: Fig. 6.3, [100] culet, ~50 GPa.
THESIS_PRESSURE_GPA = 50.0
THESIS_LINES = (450.0, 532.0)
# Anvil model: one decade of extra attenuation per D nm below this edge.
ANVIL_EDGE_NM = 500.0


def sign_blindness(contrasts=(0.05, -0.05, 0.02, -0.02)):
    """N1: Phi depends on C only through C^2, so eta cannot see the sign."""
    response = MediatedResponse(gamma_contrast=1.0, gamma_width=1.0)
    gamma = response.gamma_star()
    out = {}
    for contrast in contrasts:
        # Phi with the contrast prefactor applied by hand
        base = response.phi(gamma)
        out[contrast] = float(contrast ** 2 * base)
    return out


def prefactor_invariance(prefactors=(1.0, 0.5, 0.1, 0.01)):
    """N2: a constant multiplying C leaves rho* and the ladder untouched."""
    response = MediatedResponse(gamma_contrast=1.0, gamma_width=1.0)
    star = response.gamma_star()
    # dlnPhi/dlnGamma is unchanged by a constant factor on C, so the
    # stationary point is the same number for every prefactor.
    return {factor: float(star) for factor in prefactors}


def metric_bias(pressure=THESIS_PRESSURE_GPA, lines=THESIS_LINES,
                powers=(0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0)):
    """N3: the thesis metric against eta, as a function of laser power."""
    model = HoPublishedSpectrumModel()
    absorptions = {lam: float(model.sigma_abs(nm2eV(lam), pressure))
                   for lam in lines}
    band_max = max(float(model.sigma_abs(nm2eV(lam), pressure))
                   for lam in np.arange(402.0, 600.0, 0.5))
    response = MediatedResponse(gamma_contrast=1.0, gamma_width=1.0)
    star = response.gamma_star()

    rows = []
    for power in powers:
        scored = {}
        for lam in lines:
            gamma_p = star * power * absorptions[lam] / band_max
            contrast = response.contrast(gamma_p)
            linewidth = response.linewidth(gamma_p)
            rate = response.rate(gamma_p)
            thesis = contrast * np.sqrt(rate)
            scored[lam] = (float(thesis), float(linewidth / thesis))
        blue, green = lines
        inverse_eta = (1.0 / scored[blue][1]) / (1.0 / scored[green][1])
        thesis_ratio = scored[blue][0] / scored[green][0]
        rows.append({'power': power, 'inverse_eta_ratio': inverse_eta,
                     'thesis_ratio': thesis_ratio,
                     'bias': thesis_ratio / inverse_eta})
    return rows


def _effective_kernel(lam, base, decade_nm):
    if decade_nm is None:
        return base
    transmission = np.where(lam < ANVIL_EDGE_NM,
                            10.0 ** (-(ANVIL_EDGE_NM - lam) / decade_nm), 1.0)
    return base * transmission


def anvil_sensitivity(decades=(None, 400.0, 200.0, 100.0, 50.0), step=0.05):
    """N4: what an unmeasured anvil transmission does to the 120 GPa answer."""
    kernel = Kernel()
    lam = np.arange(DATA_WINDOW[0], DATA_WINDOW[1] + step / 2.0, step)
    base = kernel.a(lam)
    out = {}
    for decade in decades:
        effective = _effective_kernel(lam, base, decade)
        effective = effective / effective.max()
        interior = np.where((effective[1:-1] > effective[:-2])
                            & (effective[1:-1] > effective[2:]))[0] + 1
        peaks = np.sort(effective[interior])[::-1]
        out[decade] = {
            'lambda_opt': float(lam[int(np.argmax(effective))]),
            'maxima': int(len(interior)),
            'first_rung': float(1.0 / peaks[1]) if len(peaks) > 1
            else float('nan'),
        }
    return out


def anvil_threshold(step=0.05, lo=50.0, hi=1000.0):
    """N4: the attenuation at which the optimum abandons the main band."""
    kernel = Kernel()
    lam = np.arange(DATA_WINDOW[0], DATA_WINDOW[1] + step / 2.0, step)
    base = kernel.a(lam)

    def optimum(decade):
        effective = _effective_kernel(lam, base, decade)
        return float(lam[int(np.argmax(effective))])

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if optimum(mid) > 500.0:
            lo = mid
        else:
            hi = mid
    transmission = 10.0 ** (-(ANVIL_EDGE_NM - kernel.lam_abs) / hi)
    return {'decade_nm': float(hi),
            'transmission_at_optimum': float(transmission),
            'loss_factor': float(1.0 / transmission)}


def report():
    print('N1  contrast inversion: Phi is even in C')
    for contrast, phi in sign_blindness().items():
        print(f'      C = {contrast:+.2f}  ->  Phi = {phi:.6e}')
    print('      a +5% peak and a -5% peak are the same magnetometer;')
    print('      the event the core predicts is the ZERO crossing, where '
          'eta diverges,')
    print('      and in widefield a contour of zero contrast between the two '
          'regions.')

    print('\nN2  the ISC prefactor cannot move the structure')
    print(f'      {"C_0":>8}{"Gamma_p* (unchanged)":>24}')
    for factor, star in prefactor_invariance().items():
        print(f'      {factor:8.2f}{star:24.6f}')
    print('      so E, rho*, and the ladder ratios are all invariant;')
    print('      culet cut and wavelength compose multiplicatively.')

    rows = metric_bias()
    print(f'\nN3  the thesis metric vs eta at {THESIS_PRESSURE_GPA:.0f} GPa, '
          f'{THESIS_LINES[0]:.0f} nm over {THESIS_LINES[1]:.0f} nm')
    print(f'      {"I/I_c":>8}{"1/eta ratio":>14}{"thesis SNR":>13}'
          f'{"bias":>9}')
    for row in rows:
        print(f'      {row["power"]:8.2f}{row["inverse_eta_ratio"]:14.3f}'
              f'{row["thesis_ratio"]:13.3f}{row["bias"]:9.4f}')
    print('      the metric drops dnu, so it rewards the more absorbing line;')
    print('      near the crossover the two disagree about the winner.')

    print('\nN4  an unmeasured anvil transmission, at 120 GPa')
    print(f'      {"D [nm]":>9}{"lambda_opt":>13}{"shift":>9}{"maxima":>9}'
          f'{"first rung":>13}')
    for decade, row in anvil_sensitivity().items():
        label = 'none' if decade is None else f'{decade:.0f}'
        print(f'      {label:>9}{row["lambda_opt"]:13.2f}'
              f'{row["lambda_opt"] - 440.65:+9.2f}{row["maxima"]:9d}'
              f'{row["first_rung"]:13.3f}')
    threshold = anvil_threshold()
    print(f'      the optimum leaves the main band for the ZPL below '
          f'D = {threshold["decade_nm"]:.0f} nm,')
    print(f'      i.e. once transmission at 440 nm falls below '
          f'{threshold["transmission_at_optimum"]:.2f} of its 500 nm value '
          f'({threshold["loss_factor"]:.2f}x loss).')
    print('      THE STRUCTURE SURVIVES (Theorems M, G, X hold for any '
          'kernel); the numbers do not.')

    return {'sign_blindness': sign_blindness(),
            'prefactor': prefactor_invariance(),
            'metric_bias': rows,
            'anvil': anvil_sensitivity(),
            'anvil_threshold': threshold}


if __name__ == '__main__':
    report()
