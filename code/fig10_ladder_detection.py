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

The figure carries no explanatory text: everything a caption or the report's
interpretation can say has been moved out of it, and the type is sized for the
0.92\linewidth it is printed at rather than for a screen.

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

# The figure is printed at 0.92\linewidth on A4 with 22 mm margins, so it is
# reduced to about 60% of the size drawn here; the type is set accordingly.
plt.rcParams.update({
    'font.size': 13.0,
    'axes.titlesize': 13.0,
    'axes.labelsize': 12.5,
    'xtick.labelsize': 11.5,
    'ytick.labelsize': 11.5,
    'legend.fontsize': 11.0,
    'lines.linewidth': 2.0,
})


def panel_ladder(ax, rows, grid):
    for row in rows:
        if row['readable']:
            ax.axvspan(row['low'], row['high'], color='#d9ecd9', zorder=0)

    edges = np.array(PLATEAU_EDGES)
    for key, colour, style, label in (
            ('n_math', '0.55', ':', r'$N_{\rm math}$  (Theorem M level set)'),
            ('n_observable', 'k', '-',
             r'$N_{\rm obs}$  (0.5 nm resolution)')):
        values = [row[key] for row in rows]
        ax.step(edges, np.append(values, values[-1]),
                where='post', color=colour, ls=style, lw=2.6, label=label)

    marks = grid['spacing'] ** np.arange(grid['n_points'])
    ax.plot(marks[marks <= I_MAX], np.full(sum(marks <= I_MAX), 0.7), '|',
            color='tab:blue', ms=12, mew=2.0, ls='none', zorder=4,
            label=f'grid: {grid["n_points"]} points, '
                  f'ratio {grid["spacing"]:.3f}')

    split = split_threshold(SIGMA_ETA)
    ax.axvline(split, color='tab:orange', lw=2.0, ls='--',
               label=f'split resolvable, $I/I_c \\geq {split:.2f}$')

    ax.set_xscale('log')
    ax.set_xlim(1.0, I_MAX)
    ax.set_ylim(0.5, 9.2)
    ax.set_xticks([1.0, 1.441, 1.865, 2.5, 4.299, 6.0])
    ax.set_xticklabels(['1.0', '1.44', '1.87', '2.5', '4.30', '6.0'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(range(1, 7))
    ax.set_xlabel(r'excitation power  $I/I_c$')
    ax.set_ylabel('size of the optimum set')
    ax.set_title(f'(a) {grid["readable_plateaus"]} readable rungs at '
                 f'$\\sigma_\\eta$ = {SIGMA_ETA:.0%} (shaded)')
    ax.legend(loc='upper left', framealpha=0.95)
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
                    xytext=(0, 3), fontsize=10.5, ha='center', va='bottom')

    for sigma, colour in FLOORS:
        floor = 100.0 * SIGMA_MULTIPLE * sigma
        ax.axhline(floor, color=colour, lw=2.0,
                   ls='--' if sigma != SIGMA_ETA else '-')
        ax.annotate(fr'$\sigma_\eta$ = {sigma:.1%}', (6.85, floor),
                    textcoords='offset points', xytext=(0, 3), fontsize=10.5,
                    color=colour, ha='left', va='bottom')
    ax.set_yscale('log')
    ax.set_ylim(3e-5, 3e4)
    ax.set_xticks(positions)
    ax.set_xticklabels([f'{row["low"]:.2f}' for row in rows], fontsize=11.5)
    ax.set_xlim(-0.6, 8.6)
    ax.set_xlabel(r'rung, from $I/I_c$')
    ax.set_ylabel(r'smallest separating $\eta$ bump [%]')
    ax.set_title('(b) what must be resolved to count members')
    ax.grid(alpha=0.25, axis='y')


def main():
    rows = ladder(SIGMA_ETA)
    grid = power_grid(SIGMA_ETA, i_max=I_MAX)
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0),
                             constrained_layout=True)
    panel_ladder(axes[0], rows, grid)
    panel_bumps(axes[1], rows)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches='tight')
    print(f'wrote {os.path.normpath(OUT)}')


if __name__ == '__main__':
    main()
