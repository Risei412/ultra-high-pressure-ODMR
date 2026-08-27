"""
repro_yield.py
--------------
The one comparison that tests the wavelength claim against measurement.

lambda_opt is made of sigma_abs(lambda, P) and almost nothing else -- the
charge-state fraction is flat across the blue window and the ODMR contrast is
wavelength independent by construction.  So a figure that compares MEASURED
ODMR contrast with the model tests the absolute sensitivity scale, not the
wavelength recommendation.  The quantity that does test it is the one in
Ho et al. 2026 Fig. 5(b): integrated photoluminescence yield versus pressure at
a FIXED excitation wavelength, taken at constant laser power and integration
time.  As the ZPL blue shifts, the absorption band sweeps past the fixed laser
line, so the yield rises, peaks, and falls; the pressure at which it peaks is a
direct read-out of where the absorption maximum is, and it is invariant under
the arbitrary vertical scale of the measurement.

Two lines were measured: 532 nm from 4.7 to 51 GPa and 457 nm from 51 to
114 GPa.  Between them they track the absorption maximum across 2.33-2.71 eV.

What this script does:

    - loads data/pl_yield_vs_pressure.csv
    - compares it against the FROZEN model (dE120 = 400 meV)
    - re-fits the single anchor the comparison actually constrains, the ZPL
      shift at 120 GPa, and reports what that does to lambda_opt

Run:  python repro_yield.py
"""
import os

import numpy as np

from nv_model import NVModel, nm2eV

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'data', 'pl_yield_vs_pressure.csv')
LINES = (532.0, 457.0)


def load(path=DATA):
    """{'expt532': (P, y), 'expt457': ..., 'theory532_ho': ...}"""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            k, p, y = line.split(',')
            out.setdefault(k, []).append((float(p), float(y)))
    return {k: (np.array([r[0] for r in v]), np.array([r[1] for r in v]))
            for k, v in out.items()}


def predict(model, lam, P, collection=True):
    """Detected PL yield, in the model's own arbitrary units.

    What the camera integrates is not sigma_abs but the emission that survives
    a fixed long-pass filter, and the emission band blue shifts out of that
    filter as the ZPL does.  eta_col carries that, is wavelength independent,
    and therefore cannot move lambda_opt -- but leaving it out here would
    compare the model to a different quantity than the one measured.
    """
    P = np.atleast_1d(np.asarray(P, float))
    s = np.array([model.sigma_abs(nm2eV(lam), p) for p in P])
    return s * np.array([model.eta_col(p) for p in P]) if collection else s


def _scaled_residual(model, lam, P, y, collection=True):
    """Fractional residuals after the one free scale each series is allowed."""
    m = predict(model, lam, P, collection)
    k = float(np.sum(m * y) / np.sum(m * m))       # least squares in y
    return (k * m - y) / y, k


def compare(model, data, collection=True):
    """{'532': rms fractional residual, ...} plus the pooled value."""
    out, pooled = {}, []
    for lam in LINES:
        P, y = data['expt%d' % lam]
        r, _ = _scaled_residual(model, lam, P, y, collection)
        out['%d' % lam] = float(np.sqrt(np.mean(r ** 2)))
        pooled.append(r)
    out['pooled'] = float(np.sqrt(np.mean(np.concatenate(pooled) ** 2)))
    return out


def fit_dE120(data, T=300.0, grid=np.arange(0.35, 0.70, 0.005), **kw):
    """Re-anchor the ZPL shift on these data.  Returns (dE120, model, rms).

    This is the sanctioned exception to the freeze: dE120 is a measured input,
    not a model form, and this is a measurement of it.
    """
    best = None
    for dE in grid:
        m = NVModel(T=T, dE120=float(dE), **kw)
        c = compare(m, data)['pooled']
        if best is None or c < best[2]:
            best = (float(dE), m, c)
    return best


def peak_pressure(model, lam, P=np.linspace(0.0, 200.0, 401)):
    """Pressure at which the yield at this fixed line is maximal."""
    return float(P[int(np.argmax(predict(model, lam, P)))])


def _observed_peak(P, y, window=5):
    """Peak pressure of a noisy series, from a short running mean."""
    order = np.argsort(P)
    P, y = P[order], y[order]
    k = np.ones(window) / window
    sm = np.convolve(y, k, 'valid')
    return float(P[window // 2 + int(np.argmax(sm))])


def report(T=300.0):
    data = load()
    frozen = NVModel(T=T)
    dE, refit, _ = fit_dE120(data, T)

    print(f'PL yield vs pressure at a fixed line   (T = {T:.0f} K)')
    print(f'  {"":22}{"532 nm":>10}{"457 nm":>10}{"pooled":>10}')
    for name, m in (('frozen  dE120=0.400', frozen),
                    (f'refit   dE120={dE:.3f}', refit)):
        c = compare(m, data)
        print(f'  {name:22}{c["532"]*100:9.1f}%{c["457"]*100:9.1f}%'
              f'{c["pooled"]*100:9.1f}%')

    print('\n  pressure of maximum yield (GPa)')
    print(f'  {"":22}{"532 nm":>10}{"457 nm":>10}')
    for name, m in (('frozen', frozen), ('refit', refit)):
        print(f'  {name:22}{peak_pressure(m, 532.0):10.0f}'
              f'{peak_pressure(m, 457.0):10.0f}')
    obs = [_observed_peak(*data['expt%d' % lam]) for lam in LINES]
    print(f'  {"measured":22}{obs[0]:10.0f}{obs[1]:10.0f}')

    print(f'\n  lambda_opt(120 GPa)   frozen {frozen.lambda_opt(120):.1f} nm'
          f'   ->   refit {refit.lambda_opt(120):.1f} nm')
    print('  The 457 nm branch is what moves it: the frozen anchor puts the\n'
          '  absorption maximum at 469 nm at 120 GPa, and these data put it\n'
          '  bluer.  Nothing here touches hbar-omega, which stays the dominant\n'
          '  model-form uncertainty.')


if __name__ == '__main__':
    report(300.0)
    print()
    report(90.0)

