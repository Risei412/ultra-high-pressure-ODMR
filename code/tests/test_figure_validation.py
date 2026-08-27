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
    assert drivers['S_abs_start'] == pytest.approx(3.02, abs=0.05)
    assert drivers['S_abs_end'] == pytest.approx(4.55, abs=0.05)
    assert drivers['S_abs_relative_growth'] == pytest.approx(0.51, abs=0.02)
    assert drivers['dS_dP_milli_per_GPa'] == pytest.approx(12.7, abs=0.3)


def test_debye_waller_factor_collapses():
    """V3: the ZPL weight really does collapse -- by x9.07, from panel (c)."""
    drivers = drivers_from_published_panels()
    assert drivers['DWF_monotone'] is True
    assert drivers['DWF_start'] == pytest.approx(0.0205, rel=0.02)
    assert drivers['DWF_end'] == pytest.approx(0.00226, rel=0.02)
    assert drivers['DWF_fall_factor'] == pytest.approx(9.07, rel=0.02)


def test_dwf_falls_faster_than_the_single_mode_estimate():
    """V3: exp(-S) would give only x4.6, so DWF is genuinely multi-mode."""
    drivers = drivers_from_published_panels()
    assert drivers['DWF_fall_factor'] > 1.8 * drivers['exp_minus_S_fall_factor']


def test_panels_bc_cover_the_published_pressures():
    panels = load_panels_bc()
    assert list(panels['pressure']) == [0, 20, 40, 60, 80, 100, 120]


# ------------------------------------------- Theorem X survives, P* does not

def test_branch_ratio_density_is_monotone_so_the_crossing_stays_unique():
    """V4: the antecedent of Theorem X holds on published quantities."""
    data = branch_ratio_density()
    assert data['monotone'] is True
    assert np.all(np.diff(data['r']) > 0)
    assert data['growth_factor'] == pytest.approx(6.51, rel=0.02)


def test_critical_bandwidth_spans_a_few_meV():
    """V4: the ZPL only competes for a laser narrower than ~1.5-10 meV."""
    widths = branch_ratio_density()['critical_bandwidth_meV']
    assert widths[0] == pytest.approx(9.68, rel=0.02)
    assert widths[-1] == pytest.approx(1.49, rel=0.02)
    assert np.all(np.diff(widths) < 0)


def test_exchange_pressure_moves_with_bandwidth():
    """V5: P* is a function of excitation bandwidth, not a single number."""
    assert exchange_pressure_at_bandwidth(3.0e-3) == pytest.approx(77.7, abs=1.0)
    assert exchange_pressure_at_bandwidth(5.0e-3) == pytest.approx(45.2, abs=1.0)
    assert exchange_pressure_at_bandwidth(7.0e-3) == pytest.approx(24.1, abs=1.0)
    # Monotone: a broader laser favours the sideband, so the crossing moves down.
    widths = [2.0e-3, 3.0e-3, 5.0e-3, 7.0e-3]
    stars = [exchange_pressure_at_bandwidth(w) for w in widths]
    assert stars == sorted(stars, reverse=True)


def test_crossing_leaves_the_published_range_for_extreme_bandwidths():
    """V5: outside roughly 1.5-9.7 meV there is no crossing within 0-120 GPa."""
    assert not np.isfinite(exchange_pressure_at_bandwidth(1.0e-3))
    assert not np.isfinite(exchange_pressure_at_bandwidth(10.0e-3))
