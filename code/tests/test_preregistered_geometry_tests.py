"""Tests for the two pre-registered geometry measurements (G1, G2)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from preregistered_geometry_tests import (  # noqa: E402
    GEOMETRIES, GREEN_NM, KB_EV, SWEEP_TEMPERATURES_K, acceptance,
    activation_energy_eV, g1_decisive_above, g1_measured, g1_predict,
    g2_measured, g2_predict, local_zpl_eV, protocol,
)
from geometry_layer import FROZEN_PRESSURE_GPA  # noqa: E402
from nv_model import nm2eV  # noqa: E402


@pytest.fixture(scope='module')
def summary():
    return protocol()


# --- the physics the protocol rests on ------------------------------------

def test_a_red_shifted_band_puts_green_inside_the_sideband():
    axial = GEOMETRIES['[111] flat culet, [111] group']
    assert local_zpl_eV(FROZEN_PRESSURE_GPA, axial) < nm2eV(GREEN_NM)
    assert activation_energy_eV(FROZEN_PRESSURE_GPA, axial) == 0.0


def test_a_blue_shifted_band_leaves_it_outside():
    flat = GEOMETRIES['[100] flat culet']
    assert activation_energy_eV(FROZEN_PRESSURE_GPA, flat) * 1e3 > 90.0


# --- G1 -------------------------------------------------------------------

def test_g1_separates_the_two_flat_culets_by_ten_orders(summary):
    rows = {r['geometry']: r for r in summary['g1']['rows']}
    flat = rows['[100] flat culet']['collapse_factor']
    axial = rows['[111] flat culet, [111] group']['collapse_factor']
    assert flat > 1e8
    assert axial == pytest.approx(1.0)


def test_g1_predicts_collapse_only_where_the_line_is_sub_zpl(summary):
    for row in summary['g1']['rows']:
        assert row['predicted_collapse'] == (row['activation_meV'] > 0.0)


def test_g1_is_decisive_only_above_the_zpl_crossing():
    """Same 92 GPa the sub-ZPL audit (external_audit X2) found independently."""
    assert g1_decisive_above() == pytest.approx(92.0, abs=2.0)


def test_g1_recovers_a_planted_activation():
    """Round trip: synthesise a sweep at a known g, fit it back."""
    truth = 1.048
    activation = activation_energy_eV(FROZEN_PRESSURE_GPA, truth)
    temperatures = np.array(SWEEP_TEMPERATURES_K)
    counts = 1e6 * np.exp(-activation / (KB_EV * temperatures))
    result = g1_measured(temperatures, counts)
    assert result['activation_meV'] == pytest.approx(activation * 1e3, rel=1e-6)
    assert result['g'] == pytest.approx(truth, abs=0.005)
    assert result['fit_residual_rms'] < 1e-9


def test_g1_reports_no_activation_as_a_bound_not_a_value():
    temperatures = np.array(SWEEP_TEMPERATURES_K)
    result = g1_measured(temperatures, np.full(len(temperatures), 5.0e5))
    assert result['g'] is None
    assert 'no activation' in result['verdict']


def test_g1_rejects_impossible_counts():
    with pytest.raises(ValueError):
        g1_measured([300.0, 100.0], [1.0, 0.0])


# --- G2 -------------------------------------------------------------------

def test_g2_separates_the_geometries_by_far_more_than_it_needs(summary):
    assert summary['g2']['spread_nm'] > 25.0
    assert summary['g2']['resolvable'] is True


def test_g2_round_trips_each_geometry(summary):
    for row in summary['g2']['rows']:
        recovered = g2_measured(row['zpl_nm'])
        assert recovered['g'] == pytest.approx(row['g'], abs=0.002)
        assert recovered['lambda_opt_nm'] == pytest.approx(
            row['lambda_opt_nm'], abs=0.05)


def test_g2_error_bars_bracket_the_central_value(summary):
    example = summary['example']
    low, high = example['g_range']
    assert low < example['g'] < high
    lo_nm, hi_nm = example['lambda_opt_range_nm']
    assert lo_nm < example['lambda_opt_nm'] < hi_nm


def test_a_five_nm_zpl_measurement_gives_the_optimum_to_ten(summary):
    lo, hi = summary['example']['lambda_opt_range_nm']
    assert hi - lo < 10.0


# --- acceptance criteria are stated before the data exists ----------------

def test_both_tests_state_what_would_refute_them(summary):
    checks = summary['acceptance']
    assert 'fails_if' in checks['G1']
    assert 'would_refute_geometry_layer_if' in checks['G2']


def test_g2s_consistency_check_is_the_chapter_seven_bound(summary):
    assert '0.741' in summary['acceptance']['G2']['consistency_check']
