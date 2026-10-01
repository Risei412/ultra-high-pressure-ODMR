"""Extract Ho et al. Fig. 1(b) and 1(c) from the published vector figure.

    pdftocairo -svg -f 2 -l 2 2606.02399v1.pdf page2.svg
    python extract_ho_fig1_panels_bc.py page2.svg --output data/ho_fig1_panels_bc.csv

Panel (b) gives the absorption Huang-Rhys factor S_abs of the JT-inactive
modes; panel (c) gives the absorption Debye-Waller factor DWF_abs.  Both are
drawn as a polyline through seven markers at 0, 20, ... 120 GPa, so the seven
polyline vertices ARE the computed values and no curve fitting is involved.

Why this file exists
--------------------
The committed CSV was previously produced by pixel analysis of a raster
rendering, and it carried a constant baseline offset: S_abs low by 0.060 and
DWF_abs low by 0.00138 at every pressure.  On panel (c) that is invisible at
ambient pressure, where DWF_abs is 0.022, and it is a 60 % error at 120 GPa,
where DWF_abs is 0.0036 -- 7 % of a linear full scale.  It put the ZPL weight
collapse at a factor 9.07 where the source's own running text says
"more than a fivefold reduction", and it propagated into r(P), into the
critical bandwidths, and into the whole worked example of Theorem X.

Reading the vector paths removes the baseline question entirely: the axis
calibration comes from the major tick marks in the same coordinate system as
the data, and both panels then reproduce the endpoints the source states in
its text (S_abs 3.08 -> 4.61; DWF_abs 2.2 % -> 0.36 %) to better than 0.5 %.
`AUDIT` below records those four published anchors and `main()` checks against
them, so a future re-extraction cannot silently drift again.

The extraction needs only poppler's `pdftocairo`; it is not imported by
anything and is not needed to use the committed CSV.
"""
import argparse
import csv
import re
import sys

# The crimson used for both absorption curves, as pdftocairo writes it.
CRIMSON = '65.097046%,2.352905%,15.686035%'
BLACK = 'stroke:rgb(0%,0%,0%)'

PRESSURES = (0, 20, 40, 60, 80, 100, 120)

# Where each panel sits on page 2, in PDF points, and the axis values its two
# outermost major ticks carry.  The boxes only have to separate the panels.
# Boxes are cut around each panel's LEFT spine only: panel (b) carries a
# second ordinate on its right spine, whose major ticks would otherwise be
# mixed in with the S_abs ones and corrupt the scale.
PANELS = {
    'b': {'box': (60, 150, 100, 230), 'major': (5.0, 3.0), 'column': 'S_abs'},
    'c': {'box': (195, 85, 215, 150), 'major': (0.04, 0.00), 'column': 'DWF_abs'},
}

# Values Ho et al. state in their running text, as (panel, index, value).
# Agreement with these is the check that the calibration is right.
AUDIT = (('b', 0, 3.08, 0.01), ('b', -1, 4.61, 0.01),
         ('c', 0, 0.022, 0.03), ('c', -1, 0.0036, 0.03))

PATH_RE = re.compile(
    r'<path style="([^"]*)"\s+d="([^"]*)"\s+transform="matrix\(([^)]*)\)"')


def _transform(matrix, x, y):
    a, b, c, d, e, f = (float(v) for v in matrix.split(','))
    return a * x + c * y + e, b * x + d * y + f


def _paths(svg):
    for style, data, matrix in PATH_RE.findall(svg):
        points = [_transform(matrix, float(x), float(y))
                  for x, y in re.findall(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', data)]
        yield style, points


def _inside(point, box):
    left, top, right, bottom = box
    return left <= point[0] <= right and top <= point[1] <= bottom


def curve(svg, box):
    """The seven polyline vertices of the solid crimson curve inside `box`.

    Panel (b) carries two crimson curves: S_abs, drawn solid, and the squared
    vibronic coupling K^2, drawn dashed against a second ordinate we do not
    use.  The caption's own distinction -- "red solid line" against "red
    dashed line" -- is the one applied here.
    """
    found = [points for style, points in _paths(svg)
             if CRIMSON in style and 'stroke:none' not in style
             and 'dasharray' not in style
             and len(points) == len(PRESSURES) and _inside(points[0], box)]
    if len(found) != 1:
        raise AssertionError('expected one solid 7-point crimson polyline in '
                             '%s, found %d' % (box, len(found)))
    return found[0]


def major_ticks(svg, box, tick_max_pt=10.0):
    """Ordinates of the major ticks in `box`, top first.

    Ticks are two-point horizontal segments springing from a spine, and the
    major ones are drawn about twice the length of the minor ones, so a
    threshold on the run of lengths separates them.  The horizontal spines of
    the panel are also two-point horizontal segments, an order of magnitude
    longer, so they are discarded first; without that they dominate the
    threshold and no tick survives it.  Both spines of a panel carry ticks at
    the same ordinates, so the result is deduplicated.
    """
    lengths = []
    for style, points in _paths(svg):
        if BLACK not in style or len(points) != 2:
            continue
        (x0, y0), (x1, y1) = points
        if abs(y0 - y1) > 0.2 or not _inside(points[0], box):
            continue
        length = abs(x1 - x0)
        if 0.0 < length <= tick_max_pt:
            lengths.append((y0, length))
    if not lengths:
        raise AssertionError('no ticks found in %s' % (box,))
    longest = max(length for _, length in lengths)
    ticks = sorted({round(y, 3) for y, length in lengths
                    if length > 0.75 * longest})
    if len(ticks) < 2:
        raise AssertionError('fewer than two major ticks in %s: %s'
                             % (box, ticks))
    return ticks


def panel_values(svg, panel):
    """Read one panel into data units."""
    spec = PANELS[panel]
    ticks = major_ticks(svg, spec['box'])
    top_value, bottom_value = spec['major']
    top, bottom = ticks[0], ticks[-1]
    scale = (top_value - bottom_value) / (bottom - top)
    return [bottom_value + (bottom - y) * scale
            for _, y in curve(svg, spec['box'])]


def extract(svg_text):
    values = {panel: panel_values(svg_text, panel) for panel in PANELS}
    return [dict(pressure_GPa=p,
                 S_abs=values['b'][i],
                 DWF_abs=values['c'][i])
            for i, p in enumerate(PRESSURES)]


def audit(rows):
    """Check the four values Ho et al. state in their text."""
    report = []
    for panel, index, published, tolerance in AUDIT:
        column = PANELS[panel]['column']
        got = rows[index][column]
        error = abs(got / published - 1.0)
        report.append((column, PRESSURES[index], got, published, error,
                       error <= tolerance))
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('svg', help='page 2 of the source, as SVG')
    parser.add_argument('--output', default=None,
                        help='CSV to write; prints to stdout if omitted')
    args = parser.parse_args(argv)

    with open(args.svg, encoding='utf-8', errors='replace') as stream:
        rows = extract(stream.read())

    ok = True
    print('published anchors')
    for column, pressure, got, published, error, passed in audit(rows):
        ok &= passed
        print('  %-8s %3d GPa  %.5f  vs %.5f  (%.1f %%)  %s'
              % (column, pressure, got, published, error * 100,
                 'ok' if passed else 'FAIL'))
    if not ok:
        print('calibration disagrees with the published text', file=sys.stderr)
        return 1

    print('\n pressure   S_abs     DWF_abs')
    for row in rows:
        print('  %5d    %.4f    %.5f'
              % (row['pressure_GPa'], row['S_abs'], row['DWF_abs']))
    print('  S_abs  %.4f -> %.4f, dS/dP = %.5f /GPa'
          % (rows[0]['S_abs'], rows[-1]['S_abs'],
             (rows[-1]['S_abs'] - rows[0]['S_abs']) / 120.0))
    print('  DWF    %.5f -> %.5f, falls x%.2f'
          % (rows[0]['DWF_abs'], rows[-1]['DWF_abs'],
             rows[0]['DWF_abs'] / rows[-1]['DWF_abs']))

    if args.output:
        with open(args.output, 'w', newline='', encoding='utf-8') as stream:
            stream.write(
                '# Huang-Rhys factor and Debye-Waller factor for ABSORPTION,\n'
                '# read from the vector paths of the published Fig. 1 panels\n'
                '# (b) and (c) by extract_ho_fig1_panels_bc.py.  The seven\n'
                '# polyline vertices are the computed values; the axis scale\n'
                '# comes from the major ticks in the same coordinate system.\n'
                '# Reproduces the endpoints the source states in its text,\n'
                '# S_abs 3.08 -> 4.61 and DWF_abs 2.2 % -> 0.36 %, to better\n'
                '# than 0.5 %.  Panel (e) zero-phonon peak heights are clipped\n'
                '# and are never used; these two panels are not affected.\n')
            writer = csv.DictWriter(
                stream, fieldnames=['pressure_GPa', 'S_abs', 'DWF_abs'])
            writer.writeheader()
            for row in rows:
                writer.writerow({'pressure_GPa': row['pressure_GPa'],
                                 'S_abs': '%.4f' % row['S_abs'],
                                 'DWF_abs': '%.5f' % row['DWF_abs']})
        print('\nwrote %s' % args.output)
    return 0


if __name__ == '__main__':
    sys.exit(main())
