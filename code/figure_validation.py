"""Validation of the extracted kernel against the source figure, and the
correction it forces on Addendum A3.

The source Fig. 1 was checked pixel by pixel against
`data/ho_fig1e_absorption.csv`.  Two results, one reassuring and one not.

**The sideband branch is exact.**  Tracing each coloured curve in panel (e) and
taking its maximum on the absorption side reproduces the extracted CSV to
better than 1 % in height at all seven pressures, and to 0.015 eV in position
at six of them.  Everything A3 says about the sideband stands.

**The zero-phonon-line peak heights are not physical: they are clipped.**  Every
ZPL spike in panel (e) rises from its own baseline to within a few pixels of the
top of the axis --- measured tops are 14.93, 14.95, 14.67, 14.03, 14.09, 14.62
and 14.94 on an axis that ends at about 15.  They are drawn as near-delta lines
and the plot simply cuts them off.  The heights in the CSV (6.364 down to 1.054)
are therefore a digitisation artefact of a clipped feature, not a measurement,
and the "ZPL cross-section falls by x6.04" claim in A3 has no support.

The fix is to stop reading the ZPL off panel (e) at all.  Panels (b) and (c)
publish exactly what is needed --- the absorption Huang-Rhys factor S_abs and
the absorption Debye-Waller factor DWF_abs, both as theory curves with markers,
both calibrated from their own axis ticks and neither clipped.  Those are in
`data/ho_fig1_panels_bc.csv`.

What the published panels give us:

* S_abs rises 3.023 -> 4.554 over 0-120 GPa, +51 %, dS/dP = 12.7 milli/GPa,
  monotone.  **Theorem X's first driver is confirmed directly.**
* DWF_abs falls 0.0205 -> 0.00226, a factor **9.09**, monotone.  The ZPL weight
  really does collapse --- more steeply than the clipped spikes suggested.

But the comparison between a narrow line and a broad band is **bandwidth
dependent**, which the clipped-peak version silently hid.  The sideband enters
through a lineshape density (per eV); the ZPL enters through a dimensionless
weight that must be divided by whatever bandwidth actually samples it --- the
laser linewidth, or the ZPL's own width, whichever is larger.  Writing

    r(P) = lambda_SB sigma_SB / (lambda_ZPL DWF_abs)      [1/eV]

the branch ratio at excitation bandwidth W is A_SB/A_ZPL = r(P) W, so the
exchange happens where r(P) = 1/W.  r is monotone increasing (x6.5 over
0-120 GPa), so **the crossing is still unique** --- Theorem X survives intact.
What does not survive is a single number for P*: it moves with W.
"""
import os

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from theory_a3_branch_exchange import identify_branches

PANELS_BC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'data', 'ho_fig1_panels_bc.csv')

# Measured tops of the panel (e) ZPL spikes, in figure units above each curve's
# own baseline, from the pixel trace.  The axis ends near 15.
CLIPPED_ZPL_TOPS = {0: 14.929, 20: 12.947, 40: 10.667, 60: 8.026,
                    80: 6.092, 100: 4.623, 120: 2.938}
ZPL_BASELINES = {0: 0.0, 20: 2.0, 40: 4.0, 60: 6.0, 80: 8.0, 100: 10.0, 120: 12.0}
AXIS_TOP = 15.0

# Sideband peaks read off the figure by tracing, for comparison with the CSV.
FIGURE_SIDEBAND = {0: (2.1703, 2.352), 20: (2.2913, 2.191), 40: (2.4520, 2.095),
                   60: (2.5408, 1.999), 80: (2.6252, 1.902),
                   100: (2.7293, 1.830), 120: (2.8083, 1.774)}


def load_panels_bc(path=PANELS_BC):
    """S_abs and DWF_abs against pressure, from the published panels (b), (c)."""
    rows = []
    with open(path) as stream:
        for line in stream:
            text = line.strip()
            if not text or text.startswith('#') or text.startswith('pressure_'):
                continue
            rows.append([float(v) for v in text.split(',')])
    data = np.array(sorted(rows), float)
    return {'pressure': data[:, 0], 'S_abs': data[:, 1], 'DWF_abs': data[:, 2]}


def zpl_spikes_are_clipped(tolerance=1.1):
    """Every panel (e) ZPL spike reaches the top of the axis."""
    tops = {p: ZPL_BASELINES[p] + h for p, h in CLIPPED_ZPL_TOPS.items()}
    return {
        'absolute_tops': tops,
        'axis_top': AXIS_TOP,
        'all_within_tolerance_of_axis_top': all(
            AXIS_TOP - value < tolerance for value in tops.values()),
        'spread': float(max(tops.values()) - min(tops.values())),
    }


def sideband_agreement():
    """Compare the traced figure sideband peaks with the extracted CSV."""
    branches = identify_branches()
    rows = []
    for point in branches['sideband']:
        figure_e, figure_h = FIGURE_SIDEBAND[int(point.pressure)]
        rows.append({
            'pressure': point.pressure,
            'csv_eV': point.energy_ev, 'figure_eV': figure_e,
            'delta_eV': point.energy_ev - figure_e,
            'csv_sigma': point.sigma, 'figure_sigma': figure_h,
            'height_ratio': point.sigma / figure_h,
        })
    return rows


def drivers_from_published_panels():
    """Theorem X's antecedents, read off panels (b) and (c) rather than (e)."""
    panels = load_panels_bc()
    pressure, s_abs, dwf = panels['pressure'], panels['S_abs'], panels['DWF_abs']
    return {
        'S_abs_start': float(s_abs[0]), 'S_abs_end': float(s_abs[-1]),
        'S_abs_monotone': bool(np.all(np.diff(s_abs) > 0)),
        'S_abs_relative_growth': float(s_abs[-1] / s_abs[0] - 1.0),
        'dS_dP_milli_per_GPa': float(np.polyfit(pressure, s_abs, 1)[0] * 1e3),
        'DWF_start': float(dwf[0]), 'DWF_end': float(dwf[-1]),
        'DWF_monotone': bool(np.all(np.diff(dwf) < 0)),
        'DWF_fall_factor': float(dwf[0] / dwf[-1]),
        # The single-effective-mode estimate, for comparison only.
        'exp_minus_S_fall_factor': float(np.exp(s_abs[-1] - s_abs[0])),
    }


def branch_ratio_density():
    """r(P) = lambda_SB sigma_SB / (lambda_ZPL DWF_abs), in units of 1/eV.

    The sideband contributes a lineshape density; the ZPL contributes a
    dimensionless weight.  Their ratio therefore carries units of 1/eV, and the
    dimensionless branch ratio is r(P) times the excitation bandwidth.
    """
    panels = load_panels_bc()
    branches = identify_branches()
    pressure = panels['pressure']
    ratio = []
    for index, point in enumerate(branches['sideband']):
        zpl = branches['zpl'][index]
        sideband_merit = point.wavelength_nm * point.sigma
        zpl_merit = zpl.wavelength_nm * panels['DWF_abs'][index]
        ratio.append(sideband_merit / zpl_merit)
    ratio = np.array(ratio)
    return {'pressure': pressure, 'r': ratio,
            'monotone': bool(np.all(np.diff(ratio) > 0)),
            'growth_factor': float(ratio[-1] / ratio[0]),
            'critical_bandwidth_meV': 1e3 / ratio}


def exchange_pressure_at_bandwidth(bandwidth_ev):
    """Pressure at which the branches trade places, for a given bandwidth.

    Returns nan when the crossing falls outside the published pressure range.
    """
    data = branch_ratio_density()
    curve = CubicSpline(data['pressure'], np.log(data['r'] * bandwidth_ev))
    lo, hi = data['pressure'][0], data['pressure'][-1]
    if curve(lo) * curve(hi) > 0.0:
        return float('nan')
    return float(brentq(curve, lo, hi))


def main():
    bar = '=' * 78
    print(bar)
    print('Validation of the extracted kernel against the source figure')
    print(bar)

    print('\n[V1] Sideband branch: CSV against a direct trace of panel (e)')
    print('    P     CSV [eV]  fig [eV]   dE [eV]   CSV sigma  fig    ratio')
    for row in sideband_agreement():
        print(f"  {row['pressure']:4.0f}    {row['csv_eV']:.4f}   "
              f"{row['figure_eV']:.4f}   {row['delta_eV']:+.4f}    "
              f"{row['csv_sigma']:.3f}     {row['figure_sigma']:.3f}  "
              f"{row['height_ratio']:.3f}")
    print('    -> agreement better than 1 % in height at every pressure.')

    clip = zpl_spikes_are_clipped()
    print('\n[V2] ZPL spikes in panel (e) are clipped by the axis')
    print(f"    axis top ~ {clip['axis_top']:.1f}; measured spike tops:")
    print('     ' + '  '.join(f'{p}:{v:.2f}' for p, v in
                              sorted(clip['absolute_tops'].items())))
    print(f"    spread across all seven = {clip['spread']:.2f} figure units")
    print(f"    all within 1.1 of the axis top: "
          f"{clip['all_within_tolerance_of_axis_top']}")
    print('    -> the CSV ZPL heights (6.364 ... 1.054) are artefacts of a')
    print('       clipped feature.  A3\'s "x6.04 ZPL collapse" is withdrawn.')

    drivers = drivers_from_published_panels()
    print('\n[V3] Theorem X drivers, taken from panels (b) and (c) instead')
    print(f"    S_abs   {drivers['S_abs_start']:.3f} -> {drivers['S_abs_end']:.3f}  "
          f"(+{drivers['S_abs_relative_growth'] * 100:.0f} %), monotone: "
          f"{drivers['S_abs_monotone']}, dS/dP = "
          f"{drivers['dS_dP_milli_per_GPa']:.2f} milli/GPa")
    print(f"    DWF_abs {drivers['DWF_start']:.5f} -> {drivers['DWF_end']:.5f}  "
          f"(falls x{drivers['DWF_fall_factor']:.2f}), monotone: "
          f"{drivers['DWF_monotone']}")
    print(f"    (single-mode exp(-S) would give only x"
          f"{drivers['exp_minus_S_fall_factor']:.2f})")

    data = branch_ratio_density()
    print('\n[V4] Corrected branch ratio, now a density r(P) [1/eV]')
    print('    P      r [1/eV]   critical bandwidth 1/r [meV]')
    for pressure, value, width in zip(data['pressure'], data['r'],
                                      data['critical_bandwidth_meV']):
        print(f'  {pressure:4.0f}    {value:8.1f}       {width:6.2f}')
    print(f"    monotone increasing: {data['monotone']}, "
          f"total growth x{data['growth_factor']:.2f}")
    print('    -> Theorem X\'s antecedent holds; the crossing is still unique.')

    print('\n[V5] The crossing pressure now depends on excitation bandwidth')
    print('    bandwidth [meV]    P* [GPa]')
    for width_mev in (1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 9.0, 10.0):
        star = exchange_pressure_at_bandwidth(width_mev * 1e-3)
        text = 'outside 0-120' if not np.isfinite(star) else f'{star:.1f}'
        print(f'      {width_mev:5.1f}          {text}')
    print('    -> a single number for P* is not defensible; 87.9 GPa is withdrawn.')


if __name__ == '__main__':
    main()
