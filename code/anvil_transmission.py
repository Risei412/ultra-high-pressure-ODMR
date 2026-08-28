"""Bound the diamond-anvil transmission from data already in the repository.

`thesis_from_the_core.py` (N4) found the 120 GPa answer to be fragile: with
T(lambda) = 10^-((500-lambda)/D) below 500 nm, the global optimum abandons
440.65 nm for the zero-phonon line as soon as D < 321 nm -- as soon as
transmission at 440 nm falls below 0.65 of its value at 500 nm.  A(lambda) in
this theory is the NV cross section only, so nothing in it constrained T, and
the headline number was left provisional pending a measurement.

The measurement turns out to exist.  Bhattacharyya's Fig. 6.3(b) plots
cwODMR SNR at 532 nm and 450 nm against pressure through a type Ia anvil, and
the figure is VECTOR: marker centres and error bars are exact path
coordinates (`data/bhattacharyya_fig63b_snr.csv`).  Six pressures, two lines.

The extraction reads:

    P [GPa]     2.9    27.0    42.0    50.0    60.0    70.0
    450/532    0.245   0.689   0.774   0.898   1.129   1.110

which crosses unity between 50 and 60 GPa -- interpolating, at 54.4 GPa.
The kernel, knowing nothing of this figure, puts the crossover at 54.37 GPa.

Better than the crossover alone: an anvil transmission that is
pressure-independent multiplies the observed ratio by a CONSTANT
sqrt(T_450/T_532) at every pressure, so the ratio of observed to predicted is
a one-parameter fit with five degrees of freedom to spend on testing it.  It
fits, chi2/dof = 1.08, and returns

    T(450)/T(532) = 0.94,  68% CI [0.80, 1.11]

i.e. consistent with no attenuation at all.  Mapped onto N4's parameterisation
that is D ~ 2100 nm against a threshold of 321 nm: the optimum stays at
440.65 nm, and the threshold is 2.9 sigma away.

So 440.65 nm survives -- at 2.9 sigma, on a type Ia anvil, and conditional on
the kernel being right, since the fit uses it.  At -1 sigma the optimum moves
to 450.85 nm, still inside the tolerance band; at -2 sigma it jumps 74 nm to
the ZPL.  That is not comfortable enough to close the question: it downgrades
a direct measurement of the anvil spectrum from PREREQUISITE to STRONGLY
RECOMMENDED, and it says nothing about type Ib anvils, which carry more
nitrogen and are what the NV-bearing side of a DAC usually is.

Run for the tables.
"""
import os

import numpy as np

from ho_spectrum_model import HoPublishedSpectrumModel
from nv_model import nm2eV
from theory_a1_generalization import DATA_WINDOW, Kernel

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'data', 'bhattacharyya_fig63b_snr.csv')

# Calibration read off the tick stubs of Fig. 6.3(b), in PDF points.
X_TICK_PT, X_TICK_GPA = (359.17, 404.93, 450.69, 496.46), (0.0, 20.0, 40.0, 60.0)
Y_DECADE_PT = 41.46           # points per decade on the log SNR axis
Y_AT_100_PT = 201.97          # y of SNR = 10^2
# A bar this wide in log10 is assumed where the figure omits one.
MISSING_BAR_LOG10 = 0.02
ANVIL_EDGE_NM = 500.0
# N4's threshold, recomputed here rather than hard-coded.
THRESHOLD_LINES = (450.0, 532.0)


def _pressure(x_pt):
    slope = (X_TICK_GPA[-1] - X_TICK_GPA[0]) / (X_TICK_PT[-1] - X_TICK_PT[0])
    return (x_pt - X_TICK_PT[0]) * slope + X_TICK_GPA[0]


def _snr(y_pt):
    return 10.0 ** (2.0 + (Y_AT_100_PT - y_pt) / Y_DECADE_PT)


def load(path=DATA):
    """Digitised Fig. 6.3(b), converted from path coordinates to physics."""
    out = {}
    with open(path) as stream:
        for line in stream:
            if line.startswith('#') or line.startswith('line_nm'):
                continue
            lam, x_pt, y_pt, top, bottom = line.rstrip('\n').split(',')
            x_pt, y_pt = float(x_pt), float(y_pt)
            if top and bottom:
                sigma = (float(bottom) - float(top)) / 2.0 / Y_DECADE_PT
            else:
                sigma = MISSING_BAR_LOG10
            out.setdefault(float(lam), []).append(
                {'pressure_GPa': _pressure(x_pt), 'snr': _snr(y_pt),
                 'sigma_log10': sigma})
    return out


def observed_ratio(data=None, blue=450.0, green=532.0):
    """SNR(blue)/SNR(green) at each pressure, with a log10 uncertainty."""
    data = load() if data is None else data
    rows = []
    for point_blue, point_green in zip(data[blue], data[green]):
        assert abs(point_blue['pressure_GPa']
                   - point_green['pressure_GPa']) < 1e-6
        rows.append({
            'pressure_GPa': point_blue['pressure_GPa'],
            'ratio': point_blue['snr'] / point_green['snr'],
            'sigma_log10': float(np.hypot(point_blue['sigma_log10'],
                                          point_green['sigma_log10'])),
        })
    return rows


def predicted_ratio(pressure, blue=450.0, green=532.0, model=None):
    """Optical-limit SNR ratio from the NV cross section alone.

    At FIXED OPTICAL POWER -- the convention the thesis states for Fig. 6.3
    ("similar laser and microwave powers") -- the absorbed rate is
    A = lambda sigma, not sigma; see the freeze, Erratum E5.
    """
    model = HoPublishedSpectrumModel() if model is None else model
    return float(np.sqrt(blue * model.sigma_abs(nm2eV(blue), pressure)
                         / (green * model.sigma_abs(nm2eV(green), pressure))))


def kernel_crossover(blue=450.0, green=532.0, low=1.0, high=119.0):
    """Pressure at which the kernel alone puts `blue` ahead of `green`."""
    from scipy.optimize import brentq
    return float(brentq(lambda p: predicted_ratio(p, blue, green) - 1.0,
                        low, high))


def crossover(rows):
    """Pressure at which the observed ratio passes unity, by interpolation."""
    for lower, upper in zip(rows, rows[1:]):
        if lower['ratio'] <= 1.0 < upper['ratio']:
            span = upper['ratio'] - lower['ratio']
            fraction = (1.0 - lower['ratio']) / span
            return float(lower['pressure_GPa'] + fraction
                         * (upper['pressure_GPa'] - lower['pressure_GPa']))
    return float('nan')


def fit_transmission(rows=None, skip_ambient=True):
    """One-parameter fit of a pressure-independent T(blue)/T(green)."""
    rows = observed_ratio() if rows is None else rows
    used = rows[1:] if skip_ambient else rows
    residual = np.array([np.log10(row['ratio']
                                  / predicted_ratio(row['pressure_GPa']))
                         for row in used])
    weight = np.array([1.0 / row['sigma_log10'] ** 2 for row in used])
    mean = float((weight * residual).sum() / weight.sum())
    error = float(1.0 / np.sqrt(weight.sum()))
    chi2 = float((weight * (residual - mean) ** 2).sum())
    return {'log10_scale': mean, 'log10_error': error,
            'transmission_ratio': 10.0 ** (2.0 * mean),
            'ci_low': 10.0 ** (2.0 * (mean - error)),
            'ci_high': 10.0 ** (2.0 * (mean + error)),
            'chi2_per_dof': chi2 / (len(used) - 1), 'n_points': len(used)}


def _decade_from_ratio(transmission_ratio, blue=450.0):
    """D in T = 10^-((500-lambda)/D) implied by a measured T(blue)/T(green)."""
    return float((ANVIL_EDGE_NM - blue) / (-np.log10(transmission_ratio)))


def optimum_under(decade_nm, step=0.05):
    """Optimal wavelength once the anvil is folded into the kernel."""
    kernel = Kernel()
    lam = np.arange(DATA_WINDOW[0], DATA_WINDOW[1] + step / 2.0, step)
    effective = kernel.a(lam) * np.where(
        lam < ANVIL_EDGE_NM,
        10.0 ** (-(ANVIL_EDGE_NM - lam) / decade_nm), 1.0)
    return float(lam[int(np.argmax(effective))])


def verdict(fit=None, threshold_decade_nm=321.0):
    """Does the measured anvil keep the optimum inside the main band?"""
    fit = fit_transmission() if fit is None else fit
    blue, green = THRESHOLD_LINES
    critical = 10.0 ** (-(ANVIL_EDGE_NM - blue) / threshold_decade_nm)
    sigma_away = (np.log10(critical) / 2.0 - fit['log10_scale']) \
        / -fit['log10_error']
    return {'critical_ratio': float(critical),
            'sigma_above_threshold': float(sigma_away),
            'optimum_central': optimum_under(
                _decade_from_ratio(fit['transmission_ratio'])),
            'optimum_minus_1sigma': optimum_under(
                _decade_from_ratio(fit['ci_low'])),
            'optimum_minus_2sigma': optimum_under(_decade_from_ratio(
                10.0 ** (2.0 * (fit['log10_scale']
                                - 2.0 * fit['log10_error']))))}


def report():
    rows = observed_ratio()
    print('Fig. 6.3(b), digitised from the vector content')
    print(f'  {"P [GPa]":>9}{"observed":>11}{"+-1sig":>9}{"predicted":>12}'
          f'{"obs/pred":>11}')
    for row in rows:
        predicted = predicted_ratio(row['pressure_GPa'])
        print(f'  {row["pressure_GPa"]:9.1f}{row["ratio"]:11.3f}'
              f'{10 ** row["sigma_log10"] - 1:8.0%}{predicted:12.3f}'
              f'{row["ratio"] / predicted:11.3f}')
    print(f'  observed crossover  {crossover(rows):.1f} GPa')
    print(f'  predicted crossover '
          f'{kernel_crossover():.1f} GPa   (kernel alone, '
          'thesis_crosscheck T1)')

    fit = fit_transmission(rows)
    print(f'\nOne-parameter fit of a pressure-independent anvil factor '
          f'({fit["n_points"]} points, ambient dropped)')
    print(f'  T(450)/T(532) = {fit["transmission_ratio"]:.2f}   '
          f'68% CI [{fit["ci_low"]:.2f}, {fit["ci_high"]:.2f}]')
    print(f'  chi2/dof = {fit["chi2_per_dof"]:.2f}  '
          f'-> a constant factor describes the data')

    check = verdict(fit)
    print('\nAgainst N4\'s threshold')
    print(f'  the optimum leaves the main band below T(450)/T(532) = '
          f'{check["critical_ratio"]:.2f}')
    print(f'  measured value is {check["sigma_above_threshold"]:.1f} sigma '
          f'above it')
    print(f'  lambda_opt at the fit           : '
          f'{check["optimum_central"]:.2f} nm')
    print(f'  lambda_opt at -1 sigma          : '
          f'{check["optimum_minus_1sigma"]:.2f} nm')
    print(f'  lambda_opt at -2 sigma          : '
          f'{check["optimum_minus_2sigma"]:.2f} nm')
    print(f'  440.65 nm survives, but on a type Ia anvil and at '
          f'{check["sigma_above_threshold"]:.1f} sigma;')
    print('  a direct transmission measurement is still strongly recommended,')
    print('  and this says nothing about the type Ib side of a DAC.')
    return {'rows': rows, 'fit': fit, 'verdict': check,
            'crossover': crossover(rows)}


if __name__ == '__main__':
    report()
