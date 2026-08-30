"""What the multiplicity ladder looks like to an experiment, not to the theory.

Theorem M says the optimum is a level set whose size jumps with power, and A2
puts the jumps at I/I_c = 1.441, 1.514, 1.519, 1.865, 3.608, 4.299 with
N = 2, 4, 6, 4, 3, 5, 3.  `ladder_detection.py` then asks what survives two
filters -- members closer than the wavelength resolution merge, and a rung is
only counted if the eta bump BETWEEN its members clears the measurement error
-- and the answer is a much shorter ladder.  This figure is that comparison.

Panel (a) puts N_math and N_obs on the same power axis and shades the rungs
that are readable at 1% precision on eta, with the recommended calibration-free
grid (ratio 1.227, ten points) drawn underneath so the two can be read against
each other.

Panel (b) is why panel (a) loses rungs: the smallest separating bump on each
rung, against the 3-sigma detection floor at several precisions.  It spans four
orders of magnitude, so width is not what binds -- the x1.003 rung and the
x1.192 rung both fail on bump depth, not on being stepped over.

Writes `../report/image/ladder_detection.png`.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ladder_detection import (
    PLATEAU_EDGES, SIGMA_MULTIPLE, ladder, power_grid, split_threshold,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, 'report', 'image', 'ladder_detection.png')

SIGMA_ETA = 0.01
I_MAX = 6.0
# Precisions worth drawing as floors in panel (b): the recommendation, the
# next step up, and what "we saw the optimum divide" would need.
FLOORS = ((0.05, 'tab:red'), (0.01, 'tab:blue'), (0.001, 'tab:purple'))


def panel_ladder(ax, rows, grid):
    for row in rows:
        if row['readable']:
            ax.axvspan(row['low'], row['high'], color='#d9ecd9', zorder=0)

    edges = np.array(PLATEAU_EDGES)
    for key, colour, style, label in (
            ('n_math', '0.55', ':', r'$N_{\rm math}$  (Theorem M level set)'),
            ('n_observable', 'k', '-',
             r'$N_{\rm obs}$  (members within 0.5 nm merged)')):
        values = [row[key] for row in rows]
        ax.step(edges, np.append(values, values[-1]),
                where='post', color=colour, ls=style, lw=2.0, label=label)

    for ratio in grid['spacing'] ** np.arange(grid['n_points']):
        if ratio <= I_MAX:
            ax.plot([ratio], [0.7], '|', color='tab:blue', ms=11, mew=1.8,
                    zorder=4)
    ax.annotate(f'recommended grid: {grid["n_points"]} log-spaced powers, '
                f'ratio {grid["spacing"]:.3f}\n'
                'calibration free -- only the ratio between points is needed',
                (1.02, 0.95), fontsize=8.0, color='tab:blue', ha='left',
                va='bottom')

    split = split_threshold(SIGMA_ETA)
    ax.axvline(split, color='tab:orange', lw=1.4, ls='--')
    ax.annotate(f'the split itself only clears the noise\nfrom '
                f'$I/I_c = {split:.2f}$ at this precision', (split, 4.35),
                textcoords='offset points', xytext=(6, 0), fontsize=8.0,
                color='tab:orange', ha='left', va='center')

    ax.set_xscale('log')
    ax.set_xlim(1.0, I_MAX)
    ax.set_ylim(0.5, 6.6)
    ax.set_xticks([1.0, 1.441, 1.865, 2.5, 3.608, 4.299, 6.0])
    ax.set_xticklabels(['1.0', '1.44', '1.87', '2.5', '3.61', '4.30', '6.0'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(range(1, 7))
    ax.set_xlabel(r'excitation power  $I/I_c$  (ratio only; no absolute '
                  'calibration)')
    ax.set_ylabel('size of the optimum set')
    ax.set_title('(a) computed vs. observable ladder: '
                 f'{grid["readable_plateaus"]} readable rungs at '
                 f'{SIGMA_ETA:.0%} precision', fontsize=11)
    ax.legend(fontsize=8.5, loc='upper left', framealpha=0.95)
    ax.grid(alpha=0.25)


def panel_bumps(ax, rows):
    positions = np.arange(len(rows))
    bumps = np.array([100.0 * row['smallest_bump'] for row in rows])
    colours = ['tab:green' if row['readable'] else '0.75' for row in rows]

    ax.bar(positions, np.maximum(bumps, 1e-4), color=colours, width=0.62,
           zorder=3)
    for pos, row, bump in zip(positions, rows, bumps):
        ax.annotate(f'{bump:.2f}%' if bump >= 0.005 else '0.00%',
                    (pos, max(bump, 1e-4)), textcoords='offset points',
                    xytext=(0, 3), fontsize=8.0, ha='center', va='bottom')

    for sigma, colour in FLOORS:
        ax.axhline(100.0 * SIGMA_MULTIPLE * sigma, color=colour, lw=1.2,
                   ls='--' if sigma != SIGMA_ETA else '-',
                   label=fr'{SIGMA_MULTIPLE:.0f}$\sigma$ floor at '
                         fr'$\sigma_\eta$ = {sigma:.1%}')
    ax.annotate('bars above the solid line are the rungs panel (a) shades;\n'
                'the two lost rungs fail on bump DEPTH, not on width',
                (0.98, 0.97), xycoords='axes fraction', fontsize=8.0,
                color='0.25', ha='right', va='top')

    ax.set_yscale('log')
    ax.set_ylim(3e-5, 1e4)
    ax.set_xticks(positions)
    ax.set_xticklabels([f'{row["low"]:.3f}$-${row["high"]:.3f}\n'
                        fr'($\times${row["width"]:.3f})' for row in rows],
                       fontsize=8.0)
    ax.set_xlabel(r'rung, in $I/I_c$   (in brackets: rung width as a ratio)')
    ax.set_ylabel(r'smallest $\eta$ bump separating two members [%]')
    ax.set_title('(b) what has to be resolved to count the members',
                 fontsize=11)
    ax.legend(fontsize=8.0, loc='upper left')
    ax.grid(alpha=0.25, axis='y')


def main():
    rows = ladder(SIGMA_ETA)
    grid = power_grid(SIGMA_ETA, i_max=I_MAX)
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.2))
    panel_ladder(axes[0], rows, grid)
    panel_bumps(axes[1], rows)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(f'wrote {os.path.normpath(OUT)}')


if __name__ == '__main__':
    main()
