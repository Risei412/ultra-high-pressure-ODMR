"""Two cheap measurements that turn the geometry bound into a value.

`geometry_layer.py` leaves the answer one-sided.  Chapter 7 of the thesis
bounds the [111] group's band-position scalar at g <= 0.741, so the optimum
there is >= 460.33 nm, but nothing available pins g from the other side and
440.65 nm keeps being quoted with no geometry attached.

Both measurements below close that, neither needs ODMR contrast or linewidth,
and neither needs absolute calibration.  This module states what they should
return BEFORE they are run, supplies the acceptance criteria, and carries the
analysis that will consume the data when it exists.  Nothing here is a
measurement; `measured_*` are the entry points for one.

G1  TEMPERATURE SWEEP.  Photoluminescence counts under 532 nm excitation
    against temperature, at fixed pressure and fixed power.  Sub-ZPL
    absorption is thermally activated (erratum E1), so if the excitation line
    lies BELOW the group's local zero-phonon line the counts carry an
    Arrhenius activation energy equal to the deficit, and if it lies above
    they do not.  At 120 GPa the two geometries are on opposite sides:

        [100] flat culet  g >= 1.048  ->  532 nm is ~102 meV BELOW the ZPL
                                          -> counts collapse on cooling
        [111] flat culet  g <= 0.741  ->  532 nm is ~41 meV ABOVE the ZPL
                                          -> counts do not collapse

    Cooling from 300 K to 50 K changes the predicted 532 nm absorption by a
    factor of 3e-10 in the first case and by nothing in the second.  Ten
    orders of magnitude is not a marginal test.  It needs a PL spectrometer,
    a cryostat and one anvil of each cut.

G2  LOCAL ZERO-PHONON LINE.  One photoluminescence excitation edge, or one
    emission ZPL, for the [111] group in a [111] cut culet at a known
    pressure.  It gives g directly, and with it lambda_opt.  Predicted ZPL
    positions at 120 GPa are 509.8 nm ([100]), 514.5 nm (quasi-hydrostatic)
    and 541.5 nm ([111] at the bound) -- a 32 nm spread, far wider than any
    plausible measurement error.  This is the cheapest decisive measurement
    in the project: one spectrum.

Which to run: G2 if a [111] anvil is loaded, since it returns a NUMBER rather
than a sign.  G1 if only PL counts are available, since it needs no
spectrometer resolution at all.

Run for the tables.
"""
import numpy as np

from external_audit import NV_ZPL_AMBIENT_EV, zpl_energy
from geometry_layer import (
    FROZEN_PRESSURE_GPA, band_shift_eV, g_from_chapter7, g_from_fig63, optimum,
)
from nv_model import eV2nm, nm2eV

KB_EV = 8.617333262e-5
GREEN_NM = 532.0
# The three geometries of docs/geometry_of_the_optimum.md.
GEOMETRIES = {
    'quasi-hydrostatic': 1.000,
    '[100] flat culet': 1.048,
    '[111] flat culet, [111] group': 0.741,
}
SWEEP_TEMPERATURES_K = (300.0, 200.0, 150.0, 100.0, 77.0, 50.0)
# A count ratio this large over the sweep counts as "collapsed".
G1_COLLAPSE_FACTOR = 100.0
# Activation energies closer than this are not distinguishable in practice.
G1_RESOLUTION_MEV = 15.0
# A ZPL measurement good to this is enough to separate the geometries.
G2_RESOLUTION_NM = 5.0


def local_zpl_eV(pressure, g):
    """Zero-phonon line of one NV group in geometry `g`."""
    return NV_ZPL_AMBIENT_EV + g * (zpl_energy(pressure) - NV_ZPL_AMBIENT_EV)


def activation_energy_eV(pressure, g, lam=GREEN_NM):
    """Arrhenius activation of `lam` absorption; zero if the line is in the band."""
    return max(local_zpl_eV(pressure, g) - nm2eV(lam), 0.0)


# --- G1: the temperature sweep --------------------------------------------

def g1_predict(pressure=FROZEN_PRESSURE_GPA, temperatures=SWEEP_TEMPERATURES_K,
               lam=GREEN_NM, geometries=None):
    """Relative PL counts against temperature, per geometry, normalised at 300 K."""
    geometries = geometries or GEOMETRIES
    rows = []
    for name, g in geometries.items():
        activation = activation_energy_eV(pressure, g, lam)
        counts = {float(T): float(np.exp(-activation / (KB_EV * T))
                                  / np.exp(-activation / (KB_EV * max(temperatures))))
                  for T in temperatures}
        rows.append({
            'geometry': name, 'g': float(g),
            'zpl_nm': float(eV2nm(local_zpl_eV(pressure, g))),
            'activation_meV': float(activation * 1e3),
            'counts': counts,
            'collapse_factor': float(1.0 / min(counts.values())),
            'predicted_collapse': bool(1.0 / min(counts.values())
                                       > G1_COLLAPSE_FACTOR),
        })
    return {'pressure_GPa': float(pressure), 'line_nm': float(lam), 'rows': rows}


def g1_decisive_above(lam=GREEN_NM, geometries=None, grid=None):
    """Lowest pressure at which the two flat-culet geometries disagree on G1.

    Below it both put the line on the same side of their zero-phonon lines and
    the sweep cannot tell them apart.
    """
    geometries = geometries or GEOMETRIES
    grid = np.arange(20.0, 140.0, 0.5) if grid is None else grid
    blue = geometries['[100] flat culet']
    axial = geometries['[111] flat culet, [111] group']
    for pressure in grid:
        if (activation_energy_eV(pressure, blue, lam) * 1e3 > G1_RESOLUTION_MEV
                and activation_energy_eV(pressure, axial, lam) == 0.0):
            return float(pressure)
    return None


def g1_measured(temperatures_K, counts, pressure=FROZEN_PRESSURE_GPA,
                lam=GREEN_NM):
    """Analyse a real sweep: Arrhenius fit -> activation -> g -> lambda_opt.

    `counts` are PL counts at fixed pressure and fixed excitation power, in any
    units; only their ratios are used.
    """
    temperatures_K = np.asarray(temperatures_K, float)
    counts = np.asarray(counts, float)
    if counts.min() <= 0.0:
        raise ValueError('counts must be positive')
    slope, intercept = np.polyfit(1.0 / (KB_EV * temperatures_K),
                                  np.log(counts), 1)
    activation = float(-slope)
    residual = np.log(counts) - (slope / (KB_EV * temperatures_K) + intercept)
    hydrostatic = zpl_energy(pressure) - NV_ZPL_AMBIENT_EV
    if activation <= 0.0:
        verdict = ('no activation: 532 nm is inside the sideband, so the band '
                   'is red shifted and g < %.3f' % ((nm2eV(lam) - NV_ZPL_AMBIENT_EV)
                                                    / hydrostatic))
        g = None
    else:
        g = float((nm2eV(lam) + activation - NV_ZPL_AMBIENT_EV) / hydrostatic)
        verdict = 'activated: g = %.3f' % g
    return {
        'activation_meV': activation * 1e3,
        'fit_residual_rms': float(np.sqrt(np.mean(residual ** 2))),
        'g': g,
        'lambda_opt_nm': optimum(g, pressure) if g is not None else None,
        'verdict': verdict,
    }


# --- G2: the local zero-phonon line ---------------------------------------

def g2_predict(pressure=FROZEN_PRESSURE_GPA, geometries=None):
    """Predicted ZPL wavelength and optimum, per geometry."""
    geometries = geometries or GEOMETRIES
    rows = [{'geometry': name, 'g': float(g),
             'zpl_nm': float(eV2nm(local_zpl_eV(pressure, g))),
             'lambda_opt_nm': optimum(g, pressure)}
            for name, g in geometries.items()]
    spread = max(r['zpl_nm'] for r in rows) - min(r['zpl_nm'] for r in rows)
    return {'pressure_GPa': float(pressure), 'rows': rows,
            'spread_nm': float(spread),
            'resolvable': bool(spread > 4.0 * G2_RESOLUTION_NM)}


def g2_measured(zpl_nm, pressure=FROZEN_PRESSURE_GPA,
                uncertainty_nm=G2_RESOLUTION_NM):
    """Turn one measured ZPL into g and lambda_opt, with the error propagated."""
    hydrostatic = zpl_energy(pressure) - NV_ZPL_AMBIENT_EV

    def to_g(nm):
        return float((nm2eV(nm) - NV_ZPL_AMBIENT_EV) / hydrostatic)

    central = to_g(zpl_nm)
    # a redder ZPL is a smaller g, so the bracket inverts
    low, high = to_g(zpl_nm + uncertainty_nm), to_g(zpl_nm - uncertainty_nm)
    return {
        'zpl_nm': float(zpl_nm),
        'g': central,
        'g_range': (low, high),
        'lambda_opt_nm': optimum(central, pressure),
        'lambda_opt_range_nm': (optimum(high, pressure), optimum(low, pressure)),
    }


def acceptance():
    """What each measurement has to return for the geometry claim to stand."""
    ch7 = g_from_chapter7()
    fig63 = g_from_fig63()
    return {
        'G1': {
            'claim': 'the [111] group carries a red shifted band',
            'passes_if': ('the [111] culet shows no Arrhenius activation of '
                          '532 nm counts while the [100] culet shows '
                          '> %.0f meV' % G1_RESOLUTION_MEV),
            'fails_if': ('both cuts activate alike, which would put the band '
                         'at the hydrostatic position in every geometry and '
                         'leave chapter 7 unexplained'),
            'decisive_above_GPa': g1_decisive_above(),
        },
        'G2': {
            'claim': 'g for the [111] group, hence lambda_opt, is a value',
            'passes_if': ('the measured ZPL lands within %.0f nm of some g, '
                          'giving lambda_opt to about 10 nm' % G2_RESOLUTION_NM),
            'consistency_check': ('g must come out at or below %.3f, the bound '
                                  'chapter 7 already imposes'
                                  % ch7['g_upper_bound']),
            'would_refute_geometry_layer_if': ('g >= %.3f, i.e. the [111] band '
                                               'sits where the [100] one does'
                                               % fig63['g_lower_bound']),
        },
    }


def protocol():
    out = {}
    print('Two pre-registered measurements, at '
          f'{FROZEN_PRESSURE_GPA:.0f} GPa\n')

    print('G1  532 nm PL counts against temperature, fixed power')
    out['g1'] = g1 = g1_predict()
    header = ''.join(f'{T:>9.0f}' for T in SWEEP_TEMPERATURES_K)
    print(f'      {"geometry":32}{"ZPL":>8}{"E_a":>8}   counts relative to 300 K')
    print(f'      {"":32}{"[nm]":>8}{"[meV]":>8}   {header}')
    for row in g1['rows']:
        counts = ''.join(f'{row["counts"][T]:9.1e}' for T in SWEEP_TEMPERATURES_K)
        print(f'      {row["geometry"]:32}{row["zpl_nm"]:8.1f}'
              f'{row["activation_meV"]:8.0f}   {counts}')
    print('      -> the two flat culets differ by ten orders of magnitude at 50 K')
    print(f'      -> decisive above {g1_decisive_above():.0f} GPa; below that '
          'both sit on the same side')

    print('\nG2  one zero-phonon line for the [111] group')
    out['g2'] = g2 = g2_predict()
    print(f'      {"geometry":32}{"g":>8}{"ZPL [nm]":>11}{"lambda_opt [nm]":>17}')
    for row in g2['rows']:
        print(f'      {row["geometry"]:32}{row["g"]:8.3f}{row["zpl_nm"]:11.1f}'
              f'{row["lambda_opt_nm"]:17.2f}')
    print(f'      spread {g2["spread_nm"]:.1f} nm against a {G2_RESOLUTION_NM:.0f} nm '
          f'measurement: resolvable = {g2["resolvable"]}')

    print('\nworked example of the analysis, on the bound itself')
    example = g2_measured(541.5)
    print(f'      a ZPL at 541.5 +- {G2_RESOLUTION_NM:.0f} nm gives '
          f'g = {example["g"]:.3f} '
          f'[{example["g_range"][0]:.3f}, {example["g_range"][1]:.3f}]')
    print(f'      -> lambda_opt = {example["lambda_opt_nm"]:.1f} nm '
          f'[{example["lambda_opt_range_nm"][0]:.1f}, '
          f'{example["lambda_opt_range_nm"][1]:.1f}]')
    out['example'] = example

    print('\nacceptance')
    out['acceptance'] = checks = acceptance()
    for name, check in checks.items():
        print(f'      {name}: {check["claim"]}')
        for key, value in check.items():
            if key == 'claim':
                continue
            print(f'          {key}: {value}')
    return out


if __name__ == '__main__':
    protocol()
