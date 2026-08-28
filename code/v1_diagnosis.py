"""Why v1 fails against Ho et al. (2026), and how much of it is one constant.

`repro_yield.py` records that the v1 phenomenological model misses Ho's
calculated absorption badly -- pooled 39.5% RMS, and at 532 nm it is
*anti*correlated with the published curve (r = -0.51) with its maximum at
38 GPa where Ho's is at 17.  That was read as a structural failure of the
single-mode Franck--Condon picture under compression.  It is not.

Two of v1's input constants can be *derived* from Ho's own published panels,
independently of the Fig. 5(b) curves the model is scored against, and both
turn out to be wrong:

D1  Effective phonon energy.  In a single-mode Franck--Condon band the
    absorption maximum sits p* phonons above the zero-phonon line, where p*
    solves the Pekarian stationarity condition

        psi(p* + 1) = ln S            (psi = digamma)

    so that hw(P) = (E_sideband(P) - E_ZPL(P)) / p*(S_abs(P)), with
    E_sideband and E_ZPL read off Fig. 1(e) and S_abs off Fig. 1(b).  Ho's
    spectra imply hw = 101.1 meV, roughly pressure independent.  v1 uses
    65 meV, the ambient-pressure value of Kehayias et al.  The deficit is
    what puts the 532 nm maximum at 38 GPa instead of 17: the fixed-energy
    absorption peaks when E_ZPL + p* hw sweeps through the laser line, so an
    undersized hw delays the crossing.

    NOTE, because it cost a wrong answer once: the continuum shortcut
    p* = S is NOT good enough here.  exp(-S) S^p / Gamma(p+1) peaks half a
    phonon below S -- S - p* = 0.51, near-constant over 3.0 < S < 4.6 -- and
    at ~100 meV per phonon that half quantum is a 50 meV error in the band
    maximum, which leaves a visible 7 GPa residual in the peak pressure.
    Using p* = S gives hw = 87.9 meV and a model that still fails both
    gates; using the correct p* gives 101.1 meV and one that passes.

D2  ZPL shift at 120 GPa.  The reconstructed kernel puts the zero-phonon line
    at 1.946 eV at ambient and 2.410 eV at 120 GPa, a shift of 0.464 eV.  v1
    is anchored to 0.400 eV, which is what Ho's text states as a bound
    (">400 meV"), not as the value.

D3  Applying both -- and nothing else, with no parameter fitted to Fig. 5(b)
    -- moves v1 from pooled 39.5% to 1.0%, from r = -0.51 to r = +1.00 at
    532 nm, and puts both peak pressures on Ho's (17/17 and 89/88 GPa).  It
    passes both reproduction gates.  That is the same score the reconstructed
    kernel achieves.  The failure was two wrong constants, not a wrong model.

D4  This does not reinstate v1 as the optical kernel, and the reason is that
    Fig. 5(b) is a two-wavelength slice.  Scanned in WAVELENGTH at 120 GPa the
    corrected model has ONE local maximum where the kernel has FOUR, and the
    fractional disagreement is 52% even though the correlation is +0.993: the
    gross envelope is right and the structure is absent.  Every Addendum A2
    result -- the level sets, the ladder N = 2->4->6->4->3->5->3, the
    transition powers -- lives on those four maxima, and a single-mode
    Pekarian cannot produce them at any parameter value.

    One thing does survive the model change: the optical-limit optimum comes
    out at 439.10 nm against the kernel's 440.60 nm, a 1.5 nm difference,
    where uncorrected v1 gave 475.5 nm.  The frozen 440.65 nm is therefore
    robust to the choice of envelope once the constants are right.

`nv_model.py` is NOT edited.  Its defaults are frozen and `tests/test_freeze.py`
pins them; more to the point, D4 is the reason v3 takes the kernel from
outside.  The point of D1-D3 is diagnostic: the abandoned model failed for a
locatable reason, and saying so is what makes "the envelope must come from
outside" an argued choice -- argued on the structure of D4, not on a model
that could not be made to work at all.

Run for the table; `report()` returns it for the tests.
"""
import os

import numpy as np

from scipy.optimize import brentq
from scipy.special import digamma

from nv_model import NVModel, nm2eV
from ho_spectrum_model import HoPublishedSpectrumModel
from repro_yield import (
    HO_MAX_FRACTIONAL_RMS, HO_MAX_PEAK_ERROR_GPA, LINES, compare_ho_theory,
    load,
)
from theory_a3_branch_exchange import identify_branches

PANELS_BC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'data', 'ho_fig1_panels_bc.csv')

# v1's frozen inputs, for reference: nv_model.HW and NVModel's dE120 default.
V1_HW_EV = 0.065          # Kehayias et al., ambient pressure
V1_DE120_EV = 0.400       # Ho et al., quoted as a bound


def _panels_bc(path=PANELS_BC):
    """Return (pressure, S_abs, DWF_abs) digitised from Fig. 1(b),(c)."""
    rows = []
    with open(path) as stream:
        for line in stream:
            if not line[:1].isdigit():
                continue
            pressure, s_abs, dwf = line.split(',')
            rows.append((float(pressure), float(s_abs), float(dwf)))
    rows.sort()
    return (np.array([row[0] for row in rows]),
            np.array([row[1] for row in rows]),
            np.array([row[2] for row in rows]))


def pekarian_peak(s_abs):
    """Phonon order p* at which exp(-S) S^p / Gamma(p+1) is stationary."""
    return brentq(lambda p: digamma(p + 1.0) - np.log(s_abs), 0.1, 40.0)


def phonon_energy(continuum=False):
    """D1: the effective mode energy Ho's own panels imply, in eV.

    Returns (pressure, hw_per_pressure, hw_mean).  The ambient point is
    excluded from the mean: at 0 GPa the sideband maximum is the least well
    resolved of the seven, and it is the only one that departs from the rest.

    `continuum=True` reproduces the wrong shortcut p* = S, kept so that the
    tests can pin the size of the error it causes.
    """
    branches = identify_branches()
    pressure = np.array([point.pressure for point in branches['zpl']])
    zpl = np.array([point.energy_ev for point in branches['zpl']])
    sideband = np.array([point.energy_ev for point in branches['sideband']])
    panel_pressure, s_abs, _ = _panels_bc()
    if not np.array_equal(pressure, panel_pressure):
        raise ValueError('Fig. 1(e) and Fig. 1(b) are on different pressures')
    order = (s_abs if continuum
             else np.array([pekarian_peak(value) for value in s_abs]))
    hw = (sideband - zpl) / order
    return pressure, hw, float(hw[1:].mean())


def wavelength_structure(hw, dE120, pressure=120.0, T=90.0,
                         window=(402.0, 517.0), step=0.05):
    """D4: local maxima and optimum in WAVELENGTH, model against kernel."""
    lam = np.arange(window[0], window[1] + step / 2.0, step)
    energy = nm2eV(lam)
    model = NVModel(T=T, hw=hw, dE120=dE120)
    ours = np.array([model.sigma_abs(value, pressure) for value in energy])
    theirs = np.asarray(HoPublishedSpectrumModel().sigma_abs(energy, pressure),
                        float)
    ours, theirs = ours / ours.max(), theirs / theirs.max()
    interior = lambda y: int(np.sum((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])))
    return {
        'lambda_opt_model': float(lam[np.argmax(ours)]),
        'lambda_opt_kernel': float(lam[np.argmax(theirs)]),
        'maxima_model': interior(ours),
        'maxima_kernel': interior(theirs),
        'fractional_rms': float(np.sqrt(np.mean(
            ((ours - theirs) / np.maximum(theirs, 1e-3)) ** 2))),
        'correlation': float(np.corrcoef(ours, theirs)[0, 1]),
    }


def zpl_shift():
    """D2: the 0 -> 120 GPa zero-phonon shift carried by the kernel, in eV."""
    branches = identify_branches()
    zpl = [point.energy_ev for point in branches['zpl']]
    return float(zpl[-1] - zpl[0])


def audit(hw, dE120, data=None, T=90.0):
    """Score one v1 parameterisation against Ho's calculated absorption."""
    data = load() if data is None else data
    return compare_ho_theory(NVModel(T=T, hw=hw, dE120=dE120), data)


def report(T=90.0):
    """D3, D4: the table.  Returns it so the tests can assert on it."""
    data = load()
    pressure, hw_curve, hw_mean = phonon_energy()
    shift = zpl_shift()

    print('D1  effective phonon energy implied by Ho, '
          '(E_sideband - E_ZPL) / p*(S_abs)')
    for value, energy in zip(pressure, hw_curve):
        print(f'      {value:5.0f} GPa: {energy*1000:6.1f} meV')
    print(f'      mean over 20-120 GPa: {hw_mean*1000:.1f} meV     '
          f'(v1 uses {V1_HW_EV*1000:.1f} meV)')

    print(f'\nD2  kernel ZPL shift 0 -> 120 GPa: {shift:.3f} eV     '
          f'(v1 anchored to {V1_DE120_EV:.3f} eV)')

    _, _, hw_continuum = phonon_energy(continuum=True)
    cases = (
        ('v1 as frozen', V1_HW_EV, V1_DE120_EV),
        ('hw corrected only', hw_mean, V1_DE120_EV),
        ('dE120 corrected only', V1_HW_EV, shift),
        ('continuum shortcut p*=S', hw_continuum, shift),
        ('both, from Ho panels', hw_mean, shift),
    )
    results = {}
    print('\nD3  scored against Ho Fig. 5(b) calculated absorption')
    print(f'    {"":24}{"pooled":>9}'
          f'{"532 r":>9}{"532 peak":>11}{"457 r":>9}{"457 peak":>11}')
    for label, hw, dE120 in cases:
        metrics = audit(hw, dE120, data=data, T=T)
        results[label] = metrics
        print(f'    {label:24}{metrics["pooled_fractional_rms"]*100:8.1f}%'
              f'{metrics["532"]["correlation"]:+9.2f}'
              f'{metrics["532"]["peak_model_GPa"]:8.0f}/17'
              f'{metrics["457"]["correlation"]:+9.2f}'
              f'{metrics["457"]["peak_model_GPa"]:8.0f}/88')

    best = results['both, from Ho panels']
    passes = all(
        best[str(int(lam))]['fractional_rms'] <= HO_MAX_FRACTIONAL_RMS
        and best[str(int(lam))]['peak_error_GPa'] <= HO_MAX_PEAK_ERROR_GPA
        for lam in LINES)
    print(f'\n    gates ({HO_MAX_FRACTIONAL_RMS*100:.0f}% RMS, '
          f'{HO_MAX_PEAK_ERROR_GPA:.0f} GPa peak): '
          f'{"PASS" if passes else "FAIL"}')
    for lam in LINES:
        metrics = best[str(int(lam))]
        print(f'      {lam:.0f} nm: {metrics["fractional_rms"]*100:.1f}% RMS, '
              f'peak error {metrics["peak_error_GPa"]:.0f} GPa')

    structure = wavelength_structure(hw_mean, shift)
    print('\nD4  the same corrected model, scanned in wavelength at 120 GPa')
    print(f'      local maxima:  model {structure["maxima_model"]}   '
          f'kernel {structure["maxima_kernel"]}')
    print(f'      lambda_opt:    model {structure["lambda_opt_model"]:.2f} nm  '
          f'kernel {structure["lambda_opt_kernel"]:.2f} nm')
    print(f'      fractional RMS {structure["fractional_rms"]*100:.0f}%, '
          f'correlation {structure["correlation"]:+.3f}')
    print('    the envelope is right and the structure is absent; A2 lives on '
          'the maxima.')
    return {'hw_mean': hw_mean, 'hw_continuum': hw_continuum,
            'zpl_shift': shift, 'results': results, 'gates_pass': passes,
            'structure': structure}


if __name__ == '__main__':
    report()
