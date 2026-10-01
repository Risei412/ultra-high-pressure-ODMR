"""How many power points, and how precise, to see the multiplicity ladder.

PLAN.md names this the calculation to do next, because it fixes the power grid
of the proposed experiment.  Theorem M predicts the optimum set is a level set
whose size is a step function of power; A2 computes the steps at 120 GPa as
N = 2 -> 4 -> 6 -> 4 -> 3 -> 5 -> 3 with transitions at
I/I_c = 1.441, 1.514, 1.519, 1.865, 3.608, 4.299.  Whether an experiment can
SEE that is a different question, and it has two parts.

L1  The mathematical N is not the observable N.
    Two members of the level set separated by less than the wavelength
    resolution cannot be addressed separately, and the pair straddling the
    zero-phonon line is separated by 0.1 nm.  Merging members closer than
    `RESOLUTION_NM` collapses the ladder to N_obs = 2, 3, 5, 3, 2, 5, 3.

L2  Counting minima means resolving the maximum BETWEEN them.
    Two adjacent members are, by construction, exactly degenerate; what
    distinguishes "two optima" from "one broad optimum" is the bump in eta
    between them.  That bump is the observable, and it ranges from 0.00% to
    172% across the ladder.  A plateau is readable only if every one of its
    separating bumps clears `SIGMA_MULTIPLE` times the fractional measurement
    error on eta.

Together these are much more restrictive than the ×1.003 width of the N = 6
rung, which is what the freeze had recorded as the only unreachable step.  At
1% precision on eta the readable ladder is 2 -> 3 -> 3 -> 2, with the
N = 6 plateau and the five-member plateau near 500 nm both lost -- the latter
to bump depth, not to width.

L3  The power grid.
    A log-spaced grid of ratio r lands at least once inside a plateau of ratio
    w only if r <= w; with the grid phase unknown relative to the plateau
    edges, the miss probability is max(0, 1 - ln w / ln r).  The binding
    plateau is the narrowest READABLE one, not the narrowest one.

L4  What this buys.
    Because the transition powers are ratios I_k/I_c = A_max/A_k, the grid is
    calibration free: it does not need the absolute laser power, only that the
    points be log-spaced by a known factor.  What it does need is enough span,
    and I_max/I_c >= 4.3 to reach the last readable rung.

L5  The headline claim is the hardest thing to measure.  Worth saying plainly.
    The split itself opens continuously from zero at I_c, so the bump is
    0.005% at I/I_c = 1.02 and only reaches 1.7% at 1.44.  At 1% precision on
    eta the split first becomes visible around I/I_c = 1.6.  The robust
    observables are not the split but the LATER transitions, at 1.865 and
    4.299, where the zero-phonon line enters and the bumps exceed 100%.
    An experiment aiming at "we saw the optimum divide" needs sub-0.1%
    photometry; one aiming at "we saw the multiplicity change at the predicted
    calibration-free power ratios" needs 1% and ten points.

Run for the tables.
"""
import numpy as np

from theory_a1_generalization import (
    DATA_WINDOW, Kernel, MediatedResponse, eta_at, sensitivity_optima,
)

# Two wavelengths closer than this cannot be addressed separately by a tunable
# source, so the experiment sees one optimum where the theory counts two.
RESOLUTION_NM = 0.5
# A bump must clear this many standard errors to count as resolved.
SIGMA_MULTIPLE = 3.0
# Transition powers from Addendum A2, I_k/I_c.
TRANSITIONS = (1.4414, 1.5145, 1.5191, 1.8646, 3.6078, 4.2989)
PLATEAU_EDGES = (1.0,) + TRANSITIONS + (6.0,)
GRID_STEP_NM = 0.01


def _response():
    return MediatedResponse(gamma_contrast=1.0)


def plateau_structure(power_ratio, resolution_nm=RESOLUTION_NM,
                      kernel=None, response=None):
    """Observable optima and the eta bumps that separate them, at one power."""
    kernel = Kernel() if kernel is None else kernel
    response = _response() if response is None else response
    grid = np.arange(DATA_WINDOW[0], DATA_WINDOW[1], GRID_STEP_NM)
    gamma_max = response.gamma_star() * power_ratio
    eta = eta_at(kernel, response, grid, gamma_max)
    eta = eta / eta.min()

    members = list(sensitivity_optima(kernel, response, gamma_max,
                                      DATA_WINDOW)['optima'])
    merged = []
    for member in members:
        if merged and member - merged[-1][-1] < resolution_nm:
            merged[-1].append(member)
        else:
            merged.append([member])
    observable = [float(np.mean(group)) for group in merged]

    bumps = []
    for left, right in zip(observable, observable[1:]):
        window = (grid > left) & (grid < right)
        bumps.append(float(eta[window].max() - 1.0) if window.sum() > 2
                     else 0.0)
    return {'n_math': len(members), 'n_observable': len(observable),
            'members_nm': observable, 'bumps': bumps}


def ladder(sigma_eta, resolution_nm=RESOLUTION_NM,
           sigma_multiple=SIGMA_MULTIPLE):
    """The ladder as measured, at a given fractional precision on eta."""
    kernel, response = Kernel(), _response()
    rows = []
    for low, high in zip(PLATEAU_EDGES, PLATEAU_EDGES[1:]):
        structure = plateau_structure(float(np.sqrt(low * high)),
                                      resolution_nm, kernel, response)
        smallest = min(structure['bumps']) if structure['bumps'] else np.inf
        rows.append({
            'low': low, 'high': high, 'width': high / low,
            'n_math': structure['n_math'],
            'n_observable': structure['n_observable'],
            'members_nm': structure['members_nm'],
            'smallest_bump': smallest,
            'readable': bool(smallest >= sigma_multiple * sigma_eta),
        })
    return rows


def power_grid(sigma_eta, i_max=6.0, **kwargs):
    """Grid spacing, point count and miss probability for a readable ladder."""
    rows = ladder(sigma_eta, **kwargs)
    readable = [row for row in rows if row['readable']]
    if not readable:
        return {'readable_plateaus': 0, 'spacing': float('nan'),
                'n_points': 0, 'misses': []}
    binding = min(readable, key=lambda row: row['width'])
    spacing = binding['width']
    n_points = int(np.ceil(np.log(i_max) / np.log(spacing))) + 1
    misses = [{'low': row['low'], 'high': row['high'],
               'miss_probability': float(max(0.0, 1.0 - np.log(row['width'])
                                             / np.log(spacing)))}
              for row in rows]
    return {'readable_plateaus': len(readable), 'binding': binding,
            'spacing': float(spacing), 'n_points': n_points,
            'misses': misses}


def split_threshold(sigma_eta, sigma_multiple=SIGMA_MULTIPLE,
                    grid=np.arange(1.01, 1.90, 0.01)):
    """L5: lowest I/I_c at which the split clears the noise."""
    for ratio in grid:
        structure = plateau_structure(float(ratio))
        bumps = structure['bumps']
        if bumps and min(bumps) >= sigma_multiple * sigma_eta:
            return float(ratio)
    return float('nan')


def report(precisions=(0.05, 0.02, 0.01, 0.005, 0.001)):
    print('L1/L2  the ladder at 120 GPa, as computed and as observable')
    print(f'  {"plateau I/I_c":>21}{"width":>8}{"N_math":>8}{"N_obs":>7}'
          f'{"smallest bump":>16}')
    for row in ladder(0.01):
        print(f'  {row["low"]:9.3f}-{row["high"]:<9.3f}{row["width"]:8.3f}'
              f'{row["n_math"]:8d}{row["n_observable"]:7d}'
              f'{row["smallest_bump"]*100:15.2f}%')
    print(f'  (members merged below {RESOLUTION_NM} nm; the zero-phonon pair '
          'is 0.1 nm apart)')

    print(f'\nL3  what each precision on eta buys ({SIGMA_MULTIPLE:.0f} sigma '
          'to call a bump resolved)')
    print(f'  {"sigma_eta":>10}{"readable":>10}{"binding width":>15}'
          f'{"log points to I/I_c = 6":>25}')
    out = {}
    for sigma in precisions:
        grid = power_grid(sigma)
        out[sigma] = grid
        width = (f'{grid["binding"]["width"]:.3f}'
                 if grid['readable_plateaus'] else '--')
        print(f'  {sigma:9.1%}{grid["readable_plateaus"]:10d}{width:>15}'
              f'{grid["n_points"]:25d}')

    chosen = out[0.01]
    print('\nL4  the recommendation at 1% precision on eta')
    print(f'  {chosen["n_points"]} log-spaced powers, ratio '
          f'{chosen["spacing"]:.3f}, spanning I/I_c = 1 to 6')
    print('  readable rungs and the risk of stepping over each:')
    for row, miss in zip(ladder(0.01), chosen['misses']):
        mark = 'readable' if row['readable'] else 'not readable'
        print(f'    {row["low"]:6.3f}-{row["high"]:<6.3f} N_obs='
              f'{row["n_observable"]}  {mark:<12} '
              f'miss {miss["miss_probability"]:5.1%}')
    print('  the grid is calibration free: only the RATIO between points has '
          'to be known.')

    print('\nL5  the split itself is the hardest thing to see')
    print(f'  {"sigma_eta":>10}{"split visible from I/I_c":>26}')
    thresholds = {}
    for sigma in precisions:
        thresholds[sigma] = split_threshold(sigma)
        print(f'  {sigma:9.1%}{thresholds[sigma]:26.2f}')
    print('  the bump opens continuously from zero at I_c: 0.005% at 1.02, '
          '1.7% at 1.44.')
    print('  "we saw the optimum divide" needs sub-0.1% photometry;')
    print('  "the multiplicity changed at the predicted ratio" needs 1% and '
          'ten points.')
    return {'grids': out, 'split_thresholds': thresholds}


if __name__ == '__main__':
    report()
