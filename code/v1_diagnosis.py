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
    absorption maximum sits S phonons above the zero-phonon line, so

        hw(P) = (E_sideband(P) - E_ZPL(P)) / S_abs(P)

    with E_sideband and E_ZPL read off Fig. 1(e) and S_abs off Fig. 1(b).
    Ho's spectra imply hw = 87.9 meV, roughly pressure independent.  v1 uses
    65 meV, the ambient-pressure value of Kehayias et al.  The 35% deficit
    is what puts the 532 nm maximum at 38 GPa instead of 17: the fixed-energy
    absorption peaks when E_ZPL + S*hw sweeps through the laser line, so an
    undersized hw delays the crossing.

D2  ZPL shift at 120 GPa.  The reconstructed kernel puts the zero-phonon line
    at 1.946 eV at ambient and 2.410 eV at 120 GPa, a shift of 0.464 eV.  v1
    is anchored to 0.400 eV, which is what Ho's text states as a bound
    (">400 meV"), not as the value.

D3  Applying both -- and nothing else, with no parameter fitted to Fig. 5(b)
    -- moves v1 from pooled 39.5% to 10.6% and from r = -0.51 to r = +0.89 at
    532 nm.  The anticorrelation was a wrong constant, not a wrong model.

D4  What is left is genuine.  The corrected v1 still misses both peak
    pressures (24 vs 17 GPa, 101 vs 88 GPa) and so still fails the 5 GPa peak
    gate, and 10.6% still just misses the 10% shape gate.  The single-mode
    picture reproduces the trend but not the position.

This does NOT license editing `nv_model.py`.  Its defaults are frozen and
`tests/test_freeze.py` pins them; more to the point, v3 does not use v1's
absorption envelope for anything -- the optical kernel is Ho's own spectra.
The point of the correction is diagnostic: it says the abandoned model failed
for a locatable reason rather than an unlocatable one, which is what makes
"the envelope must come from outside" an argued choice instead of a retreat.

Run for the table; `report()` returns it for the tests.
"""
import numpy as np

from nv_model import NVModel
from repro_yield import (
    HO_MAX_FRACTIONAL_RMS, HO_MAX_PEAK_ERROR_GPA, LINES, compare_ho_theory,
    load,
)
from theory_a3_branch_exchange import identify_branches

PANELS_BC = 'data/ho_fig1_panels_bc.csv'

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


def phonon_energy():
    """D1: the effective mode energy Ho's own panels imply, in eV.

    Returns (pressure, hw_per_pressure, hw_mean).  The ambient point is
    excluded from the mean: at 0 GPa the sideband maximum is the least well
    resolved of the seven, and it is the only one that departs from the rest.
    """
    branches = identify_branches()
    pressure = np.array([point.pressure for point in branches['zpl']])
    zpl = np.array([point.energy_ev for point in branches['zpl']])
    sideband = np.array([point.energy_ev for point in branches['sideband']])
    panel_pressure, s_abs, _ = _panels_bc()
    if not np.array_equal(pressure, panel_pressure):
        raise ValueError('Fig. 1(e) and Fig. 1(b) are on different pressures')
    hw = (sideband - zpl) / s_abs
    return pressure, hw, float(hw[1:].mean())


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
          '(E_sideband - E_ZPL) / S_abs')
    for value, energy in zip(pressure, hw_curve):
        print(f'      {value:5.0f} GPa: {energy*1000:6.1f} meV')
    print(f'      mean over 20-120 GPa: {hw_mean*1000:.1f} meV     '
          f'(v1 uses {V1_HW_EV*1000:.1f} meV)')

    print(f'\nD2  kernel ZPL shift 0 -> 120 GPa: {shift:.3f} eV     '
          f'(v1 anchored to {V1_DE120_EV:.3f} eV)')

    cases = (
        ('v1 as frozen', V1_HW_EV, V1_DE120_EV),
        ('hw corrected only', hw_mean, V1_DE120_EV),
        ('dE120 corrected only', V1_HW_EV, shift),
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
    print(f'\nD4  the corrected model still fails both gates '
          f'({HO_MAX_FRACTIONAL_RMS*100:.0f}% RMS, '
          f'{HO_MAX_PEAK_ERROR_GPA:.0f} GPa peak):')
    for lam in LINES:
        metrics = best[str(int(lam))]
        print(f'      {lam:.0f} nm: {metrics["fractional_rms"]*100:.1f}% RMS, '
              f'peak error {metrics["peak_error_GPa"]:.0f} GPa')
    print('    the trend is recovered; the peak position is not.')
    return {'hw_mean': hw_mean, 'zpl_shift': shift, 'results': results}


if __name__ == '__main__':
    report()
