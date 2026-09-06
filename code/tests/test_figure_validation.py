"""Tests for the source-figure validation and the correction it forces on A3.

These lock in both halves of the check: that the sideband extraction is exact,
and that the zero-phonon-line peak heights are a clipping artefact.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from figure_validation import (  # noqa: E402
    branch_ratio_density, drivers_from_published_panels,
    exchange_pressure_at_bandwidth, load_panels_bc, sideband_agreement,
    zpl_spikes_are_clipped,
)


# ------------------------------------------------------- the sideband is sound

def test_sideband_heights_match_the_figure_to_one_percent():
    """V1: the extracted CSV reproduces a direct trace of panel (e)."""
    rows = sideband_agreement()
    assert len(rows) == 7
    for row in rows:
        assert abs(row['height_ratio'] - 1.0) < 0.01, row


def test_sideband_positions_match_except_one_shoulder():
    """V1: six of seven agree to 0.02 eV; 20 GPa picks a nearby shoulder."""
    rows = {row['pressure']: row for row in sideband_agreement()}
    for pressure in (0, 40, 60, 80, 100, 120):
        assert abs(rows[pressure]['delta_eV']) < 0.02, rows[pressure]
    assert abs(rows[20]['delta_eV']) < 0.08


# --------------------------------------------------------- the ZPL is clipped

def test_every_zpl_spike_reaches_the_axis_top():
    """V2: the spikes are drawn to the axis limit, so their heights are not data."""
    clip = zpl_spikes_are_clipped()
    assert clip['all_within_tolerance_of_axis_top'] is True
    assert clip['spread'] < 1.0
    for value in clip['absolute_tops'].values():
        assert value > clip['axis_top'] - 1.1


def test_clipping_is_independent_of_pressure():
    """V2: a real signal would not put all seven tops within 1 unit."""
    tops = np.array(sorted(zpl_spikes_are_clipped()['absolute_tops'].values()))
    assert float(np.std(tops)) < 0.5


# ----------------------------------------- drivers, from the published panels

def test_huang_rhys_factor_rises_monotonically():
    """V3: Theorem X's first driver, read straight off panel (b)."""
    drivers = drivers_from_published_panels()
    assert drivers['S_abs_monotone'] is True
    # The two values Ho et al. state in their running text.
    assert drivers['S_abs_start'] == pytest.approx(3.08, abs=0.03)
    assert drivers['S_abs_end'] == pytest.approx(4.61, abs=0.03)
    assert drivers['S_abs_relative_growth'] == pytest.approx(0.50, abs=0.02)
    assert drivers['dS_dP_milli_per_GPa'] == pytest.approx(12.7, abs=0.3)


def test_debye_waller_factor_collapses():
    """V3: the ZPL weight collapses by x6.0, from panel (c).

    The endpoints are the two the source states in its own running text,
    2.2 % and 0.36 %, and it describes the collapse as "more than a fivefold
    reduction".  An earlier raster extraction of this panel carried a constant
    baseline offset and reported x9.07 instead.
    """
    drivers = drivers_from_published_panels()
    assert drivers['DWF_monotone'] is True
    assert drivers['DWF_start'] == pytest.approx(0.0220, rel=0.02)
    assert drivers['DWF_end'] == pytest.approx(0.0036, rel=0.02)
    assert drivers['DWF_fall_factor'] == pytest.approx(6.01, rel=0.02)


def test_dwf_falls_faster_than_the_single_mode_estimate():
    """V3: exp(-S) gives x4.63 against the published x6.01, so the collapse
    is multi-mode.  The margin is 1.30, not the 2.0 the earlier extraction
    reported; the conclusion survives, with less room.
    """
    drivers = drivers_from_published_panels()
    assert drivers['DWF_fall_factor'] > 1.25 * drivers['exp_minus_S_fall_factor']
    assert (drivers['DWF_fall_factor'] / drivers['exp_minus_S_fall_factor']
            == pytest.approx(1.30, abs=0.03))


def test_panels_bc_cover_the_published_pressures():
    panels = load_panels_bc()
    assert list(panels['pressure']) == [0, 20, 40, 60, 80, 100, 120]


# ------------------------------------------- Theorem X survives, P* does not

def test_branch_ratio_density_is_monotone_so_the_crossing_stays_unique():
    """V4: the antecedent of Theorem X holds on published quantities."""
    data = branch_ratio_density()
    assert data['monotone'] is True
    assert np.all(np.diff(data['r']) > 0)
    assert data['growth_factor'] == pytest.approx(4.31, rel=0.02)


def test_critical_bandwidth_spans_a_few_meV():
    """V4: the ZPL only competes for a laser narrower than ~2.4-10.3 meV."""
    widths = branch_ratio_density()['critical_bandwidth_meV']
    assert widths[0] == pytest.approx(10.30, rel=0.02)
    assert widths[-1] == pytest.approx(2.39, rel=0.02)
    assert np.all(np.diff(widths) < 0)


def test_exchange_pressure_moves_with_bandwidth():
    """V5: P* is a function of excitation bandwidth, not a single number."""
    assert exchange_pressure_at_bandwidth(3.0e-3) == pytest.approx(101.2, abs=1.0)
    assert exchange_pressure_at_bandwidth(5.0e-3) == pytest.approx(56.8, abs=1.0)
    assert exchange_pressure_at_bandwidth(7.0e-3) == pytest.approx(31.0, abs=1.0)
    # Monotone: a broader laser favours the sideband, so the crossing moves down.
    widths = [2.0e-3, 3.0e-3, 5.0e-3, 7.0e-3]
    stars = [exchange_pressure_at_bandwidth(w) for w in widths]
    assert stars == sorted(stars, reverse=True)


def test_crossing_leaves_the_published_range_for_extreme_bandwidths():
    """V5: outside roughly 2.4-10.3 meV there is no crossing within 0-120 GPa."""
    assert not np.isfinite(exchange_pressure_at_bandwidth(2.0e-3))
    assert not np.isfinite(exchange_pressure_at_bandwidth(11.0e-3))
