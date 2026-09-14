"""Tests for the ladder detection-power calculation."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ladder_detection import (  # noqa: E402
    PLATEAU_EDGES, ladder, plateau_structure, power_grid, report,
    split_threshold,
)


@pytest.fixture(scope='module')
def summary():
    return report()


@pytest.fixture(scope='module')
def rows():
    return ladder(0.01)


def test_merging_shrinks_the_ladder(rows):
    """The zero-phonon pair is 0.1 nm apart and cannot be addressed twice."""
    assert any(row['n_observable'] < row['n_math'] for row in rows)
    assert [row['n_math'] for row in rows] == [2, 4, 6, 4, 3, 5, 3]
    assert [row['n_observable'] for row in rows] == [2, 3, 5, 3, 2, 5, 3]


def test_the_narrow_plateau_is_also_flat(rows):
    """The N = 6 rung fails on bump depth as well as on width."""
    narrow = rows[2]
    assert narrow['width'] == pytest.approx(1.003, abs=0.002)
    assert narrow['smallest_bump'] < 1e-4
    assert not narrow['readable']


def test_a_plateau_can_be_wide_and_still_unreadable(rows):
    """3.608-4.299 is 19% wide and lost to a 0.1% bump, not to width."""
    row = rows[5]
    assert row['width'] > 1.15
    assert row['smallest_bump'] < 0.002
    assert not row['readable']


def test_one_percent_precision_leaves_three_readable_plateaus(rows):
    assert sum(row['readable'] for row in rows) == 3


def test_more_precision_never_loses_a_plateau():
    counts = [sum(row['readable'] for row in ladder(sigma))
              for sigma in (0.05, 0.02, 0.01, 0.005, 0.001)]
    assert counts == sorted(counts)


def test_the_recommended_grid(summary):
    grid = summary['grids'][0.01]
    assert grid['spacing'] == pytest.approx(1.227, abs=0.01)
    assert grid['n_points'] == 10


def test_tighter_precision_costs_many_more_points(summary):
    """Buying the 1.441-1.514 rung takes the count from 10 to 38."""
    assert summary['grids'][0.005]['n_points'] > 3 * \
        summary['grids'][0.01]['n_points']


def test_readable_plateaus_are_never_missed(summary):
    grid = summary['grids'][0.01]
    for row, miss in zip(ladder(0.01), grid['misses']):
        if row['readable']:
            assert miss['miss_probability'] == pytest.approx(0.0, abs=1e-9)


def test_the_split_itself_needs_more_than_the_ladder(summary):
    """L5: the headline claim is the hardest measurement in the experiment."""
    thresholds = summary['split_thresholds']
    assert thresholds[0.01] > 1.5
    assert thresholds[0.001] < thresholds[0.01]


def test_the_bump_opens_continuously_from_the_critical_power():
    depths = [min(plateau_structure(x)['bumps'] or [0.0])
              for x in (1.02, 1.05, 1.11, 1.20, 1.30)]
    assert depths == sorted(depths)
    assert depths[0] < 1e-4


def test_plateau_edges_match_the_frozen_transitions():
    assert PLATEAU_EDGES[1:-1] == (1.4414, 1.5145, 1.5191, 1.8646, 3.6078,
                                   4.2989)
