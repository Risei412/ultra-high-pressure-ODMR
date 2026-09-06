"""Extract Ho et al. Fig. 1(a) from the published vector figure.

    pdftocairo -svg -f 2 -l 2 2606.02399v1.pdf page2.svg
    python extract_ho_fig1a_zpl.py page2.svg --output data/ho_fig1a_zpl.csv

Panel (a) is the ZPL shift relative to ambient pressure.  Its legend reads
``expt.`` / ``thr. alpha = 1`` / ``thr. alpha = 0.95``, and the running text
says:

    "We calculated the ZPL position for perfect hydrostatic stress (orange
     solid line) and the lower branch of the splitted ZPL for alpha = 0.95
     (orange dashed line)."

Why this file exists
--------------------
Addendum A4 assumed the anisotropy parameter alpha could not enter the frozen
optical layer at all.  It can: **Ho et al. computed the anisotropic case
themselves.**  The solid/dashed pair in panel (a) is a two-point calibration of
the deviatoric optical coupling at the same level of theory that produced the
absorption kernel, and it is the only such calibration that does not have to be
imported from a different experiment.

Both curves are single polylines of 15 vertices drawn in the same orange, the
solid one for alpha = 1 and the dashed one for alpha = 0.95, so the vertices
ARE the computed values and no curve fitting is involved.  The vertices land on
0, ..., 120 GPa in 15 equal steps once the axis calibration below is applied.

Axis calibration comes from the major tick marks of panel (a) in the same
coordinate system as the data: x ticks at 0/25/50/75/100 GPa and y ticks at
0.0/0.2/0.4 eV.

`AUDIT` records the one published anchor available for this panel: erratum E4
of the freeze independently derived the corrected hydrostatic ZPL shift at
120 GPa as 0.464 eV, from Fig. 1(b),(e) plus the Pekarian stationarity
condition -- a completely different route through the paper.  `main()` checks
the extracted alpha = 1 endpoint against it, so a future re-extraction cannot
silently drift.

The extraction needs only poppler's `pdftocairo`; it is not imported by
anything and is not needed to use the committed CSV.
"""
import argparse
import csv
import re
import sys

# The orange used for both theory curves, as pdftocairo writes it.
ORANGE = '94.116211%,34.901428%,11.764526%'

# Panel (a) major tick marks, in the local coordinate system shared by the
# curves and the axes (transform matrix(0.88196,0,0,-0.88196,54,243.379)).
X_TICK_0_GPA = 33.6918
X_TICK_100_GPA = 102.8914
Y_TICK_0_EV = 111.5990
Y_TICK_0P4_EV = 172.4410

N_VERTICES = 15
P_MAX_GPA = 120.0

# Freeze erratum E4, derived independently of this panel.
AUDIT = {'hydrostatic_dZPL_120GPa_eV': 0.464, 'tolerance': 0.01}


def _paths(svg_text):
    for path in re.findall(r'<path[^>]*/>', svg_text):
        if ORANGE not in path:
            continue
        data = re.search(r'\sd="([^"]*)"', path)
        if not data:
            continue
        points = [(float(x), float(y)) for x, y in
                  re.findall(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', data.group(1))]
        if len(points) != N_VERTICES:
            continue
        yield ('alpha_0.95' if 'dasharray' in path else 'alpha_1.0'), points


def extract(svg_path):
    """Return {'pressure_GPa': [...], 'alpha_1.0': [...], 'alpha_0.95': [...]}."""
    with open(svg_path, encoding='utf-8') as handle:
        curves = dict(_paths(handle.read()))
    missing = {'alpha_1.0', 'alpha_0.95'} - set(curves)
    if missing:
        raise SystemExit('panel (a) curves not found: %s' % ', '.join(sorted(missing)))

    dx_per_gpa = (X_TICK_100_GPA - X_TICK_0_GPA) / 100.0
    dy_per_ev = (Y_TICK_0P4_EV - Y_TICK_0_EV) / 0.4

    out = {}
    for name, points in curves.items():
        out[name] = [(y - Y_TICK_0_EV) / dy_per_ev for _, y in points]
    out['pressure_GPa'] = [(x - X_TICK_0_GPA) / dx_per_gpa
                           for x, _ in curves['alpha_1.0']]
    return out


def deviatoric_coupling(data):
    """Return the alpha = 0.95 offset at 120 GPa and its raw effective k.

    At sigma_zz = 120 GPa and alpha = 0.95 the stress invariants are
    P_mean = 116 GPa and q = 6 GPa.  ``k_effective`` attributes the WHOLE
    offset from the hydrostatic curve to an effective hydrostatic pressure.
    For a (100) anvil that is an upper bound in magnitude, because the dashed
    curve is the lower branch of the split ZPL and therefore carries half the
    3E orbital splitting on top of the symmetric (A1) deviatoric shift.  For
    [111]-axisymmetric stress on a [111]-aligned NV the splitting vanishes by
    symmetry and only the A1 part applies.
    """
    import numpy as np

    pressure = np.array(data['pressure_GPa'])
    hydrostatic = np.array(data['alpha_1.0'])
    alpha = 0.95
    mean = P_MAX_GPA * (1.0 + 2.0 * alpha) / 3.0
    dev = P_MAX_GPA * (1.0 - alpha)

    shift = data['alpha_0.95'][-1]
    at_mean = float(np.interp(mean, pressure, hydrostatic))
    effective = float(np.interp(shift, hydrostatic, pressure))
    return {
        'P_mean_GPa': mean,
        'q_GPa': dev,
        'dZPL_hydrostatic_at_P_mean_eV': at_mean,
        'dZPL_alpha095_eV': shift,
        'offset_meV': 1e3 * (shift - at_mean),
        'P_effective_GPa': effective,
        'k_effective': (effective - mean) / dev,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('svg', help='page 2 of the paper, as SVG')
    parser.add_argument('--output', default=None)
    args = parser.parse_args(argv)

    data = extract(args.svg)

    hydrostatic_120 = data['alpha_1.0'][-1]
    error = abs(hydrostatic_120 - AUDIT['hydrostatic_dZPL_120GPa_eV'])
    print('hydrostatic dZPL(120 GPa) = %.4f eV  (E4 audit %.3f eV, error %.4f)'
          % (hydrostatic_120, AUDIT['hydrostatic_dZPL_120GPa_eV'], error))
    if error > AUDIT['tolerance']:
        raise SystemExit('extraction disagrees with the E4 audit anchor')

    coupling = deviatoric_coupling(data)
    for key in ('P_mean_GPa', 'q_GPa', 'dZPL_hydrostatic_at_P_mean_eV',
                'dZPL_alpha095_eV', 'offset_meV', 'P_effective_GPa',
                'k_effective'):
        print('%-32s %+.4f' % (key, coupling[key]))

    if args.output:
        with open(args.output, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle)
            writer.writerow(['pressure_GPa', 'dZPL_alpha_1.0_eV',
                             'dZPL_alpha_0.95_eV'])
            for row in zip(data['pressure_GPa'], data['alpha_1.0'],
                           data['alpha_0.95']):
                writer.writerow(['%.4f' % value for value in row])
        print('wrote %s' % args.output)
    return 0


if __name__ == '__main__':
    sys.exit(main())
