"""Audit Ho et al. (2026) Fig. 5(b) without mixing observables.

``expt*`` contains detected, spectrally integrated PL. A passband factor may
enter this observable, so absorption-only and detected-yield proxies are kept
separate. ``theory*_ho`` contains Ho et al.'s calculated absorption cross
section and must be compared with ``sigma_abs`` only.
"""
import os

import numpy as np

from nv_model import NVModel, nm2eV
from ho_spectrum_model import HoPublishedSpectrumModel

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'data', 'pl_yield_vs_pressure.csv')
LINES = (532.0, 457.0)

# Generous audit tolerances for smooth vector curves and experimental sampling.
HO_MAX_FRACTIONAL_RMS = 0.10
HO_MAX_PEAK_ERROR_GPA = 5.0
EXPT_MAX_FRACTIONAL_RMS = 0.16
# Experimental sampling has a 16 GPa gap around the blue dome; 10 GPa is a
# resolution-aware peak tolerance. The smooth Ho-theory curve keeps 5 GPa.
EXPT_MAX_PEAK_ERROR_GPA = 10.0


def load(path=DATA):
    """Return all experimental and Ho-theory series from the digitised CSV."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            key, pressure, value = line.split(',')
            out.setdefault(key, []).append((float(pressure), float(value)))
    return {
        key: (np.array([row[0] for row in values]),
              np.array([row[1] for row in values]))
        for key, values in out.items()
    }


def predict_absorption(model, lam, pressure):
    """Absorption cross section at a fixed excitation wavelength."""
    pressure = np.atleast_1d(np.asarray(pressure, float))
    return np.array([model.sigma_abs(nm2eV(lam), p) for p in pressure])


def predict_detected_yield(model, lam, pressure):
    """Detected-PL proxy for the explicitly assumed fixed passband."""
    pressure = np.atleast_1d(np.asarray(pressure, float))
    return predict_absorption(model, lam, pressure) * np.asarray(
        model.eta_col(pressure), float)


def predict(model, lam, pressure, collection=True):
    """Legacy wrapper; new analyses should call an explicit proxy."""
    if collection:
        return predict_detected_yield(model, lam, pressure)
    return predict_absorption(model, lam, pressure)


def _scaled_residual_from_prediction(prediction, reference):
    """Fractional residual after the one allowed multiplicative scale."""
    prediction = np.asarray(prediction, float)
    reference = np.asarray(reference, float)
    if prediction.shape != reference.shape:
        raise ValueError('prediction and reference must have the same shape')
    if np.any(reference <= 0.0):
        raise ValueError('reference values must be positive in the audit window')
    denom = float(np.sum(prediction * prediction))
    if denom <= 0.0:
        raise ValueError('prediction has zero norm')
    scale = float(np.sum(prediction * reference) / denom)
    residual = (scale * prediction - reference) / reference
    return residual, scale


def _observed_peak(pressure, values, window=5):
    """Peak pressure of a noisy series, estimated with a running mean."""
    order = np.argsort(pressure)
    pressure, values = pressure[order], values[order]
    kernel = np.ones(window) / window
    smoothed = np.convolve(values, kernel, 'valid')
    return float(pressure[window // 2 + int(np.argmax(smoothed))])


def _curve_metrics(pressure, prediction, reference, noisy=False):
    """Scale-free shape and peak metrics for one curve."""
    pressure = np.asarray(pressure, float)
    prediction = np.asarray(prediction, float)
    reference = np.asarray(reference, float)
    residual, scale = _scaled_residual_from_prediction(prediction, reference)
    peak_ref = (_observed_peak(pressure, reference) if noisy else
                float(pressure[int(np.argmax(reference))]))
    peak_model = float(pressure[int(np.argmax(prediction))])
    return {
        'fractional_rms': float(np.sqrt(np.mean(residual ** 2))),
        'scale': scale,
        'correlation': float(np.corrcoef(prediction, reference)[0, 1]),
        'peak_model_GPa': peak_model,
        'peak_reference_GPa': peak_ref,
        'peak_error_GPa': abs(peak_model - peak_ref),
        '_residual': residual,
    }


def compare_experiment(model, data, collection):
    """Compare with measured PL; ``collection`` must be stated explicitly."""
    out, pooled = {}, []
    for lam in LINES:
        pressure, values = data['expt%d' % lam]
        prediction = predict(model, lam, pressure, collection=collection)
        metrics = _curve_metrics(pressure, prediction, values, noisy=True)
        pooled.append(metrics.pop('_residual'))
        out['%d' % lam] = metrics
    out['pooled_fractional_rms'] = float(
        np.sqrt(np.mean(np.concatenate(pooled) ** 2)))
    return out


def compare_ho_theory(model, data):
    """Compare ``sigma_abs`` with Ho's calculated absorption curves only."""
    out, pooled = {}, []
    for lam in LINES:
        pressure, values = data['theory%d_ho' % lam]
        measured_pressure = data['expt%d' % lam][0]
        mask = ((pressure >= measured_pressure.min())
                & (pressure <= measured_pressure.max()))
        pressure, values = pressure[mask], values[mask]
        prediction = predict_absorption(model, lam, pressure)
        metrics = _curve_metrics(pressure, prediction, values, noisy=False)
        pooled.append(metrics.pop('_residual'))
        out['%d' % lam] = metrics
    out['pooled_fractional_rms'] = float(
        np.sqrt(np.mean(np.concatenate(pooled) ** 2)))
    return out


def compare(model, data, collection=True):
    """Legacy compact result for callers that consume the old API."""
    audit = compare_experiment(model, data, collection=collection)
    return {
        '532': audit['532']['fractional_rms'],
        '457': audit['457']['fractional_rms'],
        'pooled': audit['pooled_fractional_rms'],
    }


def fit_dE120(data, T=300.0, grid=np.arange(0.35, 0.70, 0.005),
              collection=True, target='experiment', **model_kw):
    """Return a conditional one-parameter calibration, not a measurement."""
    if target not in ('experiment', 'ho_theory'):
        raise ValueError("target must be 'experiment' or 'ho_theory'")
    best = None
    for dE in grid:
        model = NVModel(T=T, dE120=float(dE), **model_kw)
        if target == 'experiment':
            score = compare_experiment(
                model, data, collection=collection)['pooled_fractional_rms']
        else:
            score = compare_ho_theory(model, data)['pooled_fractional_rms']
        if best is None or score < best[2]:
            best = (float(dE), model, score)
    return best


def peak_pressure(model, lam, pressure=np.linspace(0.0, 200.0, 401),
                  collection=True):
    """Peak of the explicitly selected observable at a fixed laser line."""
    prediction = predict(model, lam, pressure, collection=collection)
    return float(pressure[int(np.argmax(prediction))])


def reproduction_gate(model, data):
    """Return independent Ho-theory and experiment-shape gate results."""
    ho = compare_ho_theory(model, data)
    expt = compare_experiment(model, data, collection=True)
    ho_pass = all(
        ho[str(int(lam))]['fractional_rms'] <= HO_MAX_FRACTIONAL_RMS
        and ho[str(int(lam))]['peak_error_GPa'] <= HO_MAX_PEAK_ERROR_GPA
        for lam in LINES)
    expt_pass = all(
        expt[str(int(lam))]['fractional_rms'] <= EXPT_MAX_FRACTIONAL_RMS
        and expt[str(int(lam))]['peak_error_GPa'] <= EXPT_MAX_PEAK_ERROR_GPA
        for lam in LINES)
    return {'ho_theory': ho_pass, 'experiment_shape': expt_pass,
            'overall': ho_pass and expt_pass, 'ho': ho, 'experiment': expt}


def _print_metrics(label, audit):
    print(label)
    print(f'  {"":10}{"frac RMS":>12}{"corr":>9}{"peak model/ref":>20}')
    for lam in LINES:
        metrics = audit[str(int(lam))]
        print(f'  {int(lam):>4} nm  {metrics["fractional_rms"]*100:10.1f}%'
              f'{metrics["correlation"]:9.2f}'
              f'{metrics["peak_model_GPa"]:8.0f}/'
              f'{metrics["peak_reference_GPa"]:<8.0f}')
    print(f'  pooled fractional RMS: {audit["pooled_fractional_rms"]*100:.1f}%')


def report(T=300.0):
    data = load()
    frozen = NVModel(T=T)
    dE, refit, _ = fit_dE120(data, T=T, collection=True,
                             target='experiment')

    print(f'Ho Fig. 5(b) observable audit (T = {T:.0f} K)')
    _print_metrics('\nFrozen vs experiment: detected-yield proxy',
                   compare_experiment(frozen, data, collection=True))
    _print_metrics('\nFrozen vs experiment: absorption-only proxy',
                   compare_experiment(frozen, data, collection=False))
    _print_metrics('\nFrozen vs Ho calculated absorption',
                   compare_ho_theory(frozen, data))
    _print_metrics(f'\nExperiment-calibrated (dE120={dE:.3f} eV) '
                   'vs experiment: detected-yield proxy',
                   compare_experiment(refit, data, collection=True))
    _print_metrics('\nSame calibrated model vs Ho calculated absorption',
                   compare_ho_theory(refit, data))

    published = HoPublishedSpectrumModel()
    _print_metrics('\nFig. 1(e) reconstruction vs Ho Fig. 5(b) absorption',
                   compare_ho_theory(published, data))
    _print_metrics('\nFig. 1(e) reconstruction vs experimental PL',
                   compare_experiment(published, data, collection=False))

    gate = reproduction_gate(refit, data)
    print('\nReproduction gates for the experiment-calibrated model')
    print(f'  Ho theory       : {"PASS" if gate["ho_theory"] else "FAIL"}')
    print(f'  experiment shape: {"PASS" if gate["experiment_shape"] else "FAIL"}')
    print(f'  overall         : {"PASS" if gate["overall"] else "FAIL"}')
    print(f'  lambda_opt(120 GPa): {refit.lambda_opt(120):.1f} nm')
    print('  dE120 above is a conditional fit parameter, not a direct ZPL '
          'measurement.')

    published_gate = reproduction_gate(published, data)
    print('\nPublished-curve cross-figure reproduction')
    print(f'  Ho theory       : '
          f'{"PASS" if published_gate["ho_theory"] else "FAIL"}')
    print(f'  experiment shape: '
          f'{"PASS" if published_gate["experiment_shape"] else "FAIL"}')
    print(f'  overall         : '
          f'{"PASS" if published_gate["overall"] else "FAIL"}')
    print(f'  lambda_opt(120 GPa): {published.lambda_opt(120):.1f} nm')
    print('  This reproduces published curves; it is not an independent DFT/JT '
          'calculation.')


if __name__ == '__main__':
    report(300.0)
    print()
    report(90.0)

