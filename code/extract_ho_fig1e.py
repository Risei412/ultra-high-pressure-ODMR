"""Extract Ho et al. Fig. 1(e) absorption curves from the arXiv vector PDF.

The source archive is available at https://arxiv.org/src/2606.02399v1 .
Extract figures/theory1ev_V4.pdf and run:

    python extract_ho_fig1e.py theory1ev_V4.pdf \
        --output data/ho_fig1e_absorption.csv

The SHA256 guard makes the path indices below fail loudly if the source figure
is revised. PyMuPDF is needed only to regenerate the committed CSV, not to use
the reference model.
"""
import argparse
import csv
import hashlib
import sys

import numpy as np

EXPECTED_SHA256 = (
    'ba18ac4c15f68f024d5edad3915d1d4ff40d0ea95b9f9adc50516613b92e4382'
)

# Matplotlib vector paths for the absorption half of the 0--120 GPa spectra.
PATH_INDEX = {0: 486, 20: 482, 40: 478, 60: 474,
              80: 470, 100: 466, 120: 462}

# Published axes of panel (e), in PDF points and electronvolts.
X_LEFT, X_RIGHT = 288.38, 563.50
E_LEFT, E_RIGHT = 1.4, 3.2


def _line_points(drawing):
    points = []
    for item in drawing['items']:
        if item[0] != 'l':
            continue
        if not points:
            points.append((item[1].x, item[1].y))
        points.append((item[2].x, item[2].y))
    if len(points) < 20:
        raise RuntimeError('expected a sampled Matplotlib line path')
    points = np.asarray(points, float)
    if np.any(np.diff(points[:, 0]) <= 0.0):
        raise RuntimeError('energy-axis path is not strictly monotone')
    return points


def extract(pdf_path):
    try:
        import fitz
    except ImportError as exc:
        raise RuntimeError('PyMuPDF is required to regenerate the CSV') from exc

    with open(pdf_path, 'rb') as stream:
        digest = hashlib.sha256(stream.read()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise RuntimeError(
            f'unrecognised source figure SHA256 {digest}; '
            'audit path indices before extraction')

    page = fitz.open(pdf_path)[0]
    drawings = page.get_drawings()
    raw, baselines = {}, []
    for pressure, index in PATH_INDEX.items():
        points = _line_points(drawings[index])
        raw[pressure] = points
        baselines.append((pressure, points[0, 1]))

    # Curves are vertically offset by P/10 on the published y axis. Infer the
    # PDF-point scale from all seven baselines rather than one pair.
    slope, intercept = np.polyfit(*np.asarray(baselines).T, 1)
    y_points_per_unit = -10.0 * slope
    if not (10.3 < y_points_per_unit < 10.5):
        raise RuntimeError('unexpected y-axis calibration')

    rows = []
    for pressure in sorted(raw):
        x, y = raw[pressure].T
        energy = E_LEFT + (x - X_LEFT) * (
            (E_RIGHT - E_LEFT) / (X_RIGHT - X_LEFT))
        baseline = slope * pressure + intercept
        absorption = np.clip(
            (baseline - y) / y_points_per_unit, 0.0, None)
        for e, value in zip(energy, absorption):
            rows.append((pressure, e, value))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pdf', help='path to figures/theory1ev_V4.pdf')
    parser.add_argument('--output', help='CSV path; stdout when omitted')
    args = parser.parse_args()

    stream = open(args.output, 'w', newline='') if args.output else sys.stdout
    try:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['pressure_GPa', 'energy_eV', 'absorption_au'])
        for pressure, energy, value in extract(args.pdf):
            writer.writerow([f'{pressure:d}', f'{energy:.8f}', f'{value:.8f}'])
    finally:
        if args.output:
            stream.close()


if __name__ == '__main__':
    main()

