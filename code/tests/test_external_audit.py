"""Regression tests for the Ho / thesis explanatory audit (X1-X5)."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from external_audit import (  # noqa: E402
    geometry_degeneracy, green_at_megabar, power_convention, report,
    residual_trend, sub_zpl_coverage, zpl_energy,
)


@pytest.fixture(scope='module')
def summary():
    return report()


# --- X1 -------------------------------------------------------------------

def test_the_power_convention_moves_the_crossover_by_two_and_a_half_GPa(summary):
    x1 = summary['power_convention']
    assert x1['crossover_fixed_flux_GPa'] == pytest.approx(51.82, abs=0.1)
    assert x1['crossover_fixed_power_GPa'] == pytest.approx(54.37, abs=0.1)


def test_the_stated_convention_is_the_one_that_hits(summary):
    """The thesis compared at fixed laser power, and that is the row that lands."""
    x1 = summary['power_convention']
    observed = x1['crossover_observed_GPa']
    assert abs(x1['crossover_fixed_power_GPa'] - observed) < 0.2
    assert abs(x1['crossover_fixed_flux_GPa'] - observed) > 2.0


def test_the_tie_at_fifty_GPa_survives_the_correction(summary):
    """Fig. 6.3's null result is still the prediction, at 0.87 rather than 0.94."""
    x1 = summary['power_convention']
    assert 0.8 < x1['snr_ratio_50GPa_fixed_power'] < 1.0


# --- X2 -------------------------------------------------------------------

def test_e1_reaches_down_to_one_hundred_GPa(summary):
    x2 = summary['sub_zpl']
    assert x2['last_sampled_GPa'] == pytest.approx(80.0)
    flags = {row['pressure_GPa']: row['bracketed_by_real_samples']
             for row in x2['rows']}
    assert flags[80.0] is True
    assert flags[100.0] is False
    assert flags[120.0] is False


def test_the_empty_interval_opens_where_the_zpl_passes_the_line(summary):
    x2 = summary['sub_zpl']
    assert 80.0 < x2['zpl_crossing_GPa'] < 100.0


def test_a_line_inside_the_band_is_bracketed_at_every_pressure():
    """457 nm never goes sub-ZPL in the published range, so E1 never applies."""
    coverage = sub_zpl_coverage(457.0)
    assert all(row['bracketed_by_real_samples'] for row in coverage['rows'])


# --- X3 -------------------------------------------------------------------

def test_the_kernel_zpl_extrapolates_above_the_published_range():
    assert zpl_energy(120.0) == pytest.approx(2.410, abs=0.002)
    assert zpl_energy(140.0) > zpl_energy(130.0) > zpl_energy(120.0)


def test_green_is_deeply_sub_zpl_at_chapter_sevens_pressures(summary):
    rows = summary['green_at_megabar']['rows']
    assert all(row['deficit_meV'] > 100.0 for row in rows)


def test_cryogenic_hot_band_absorption_cannot_account_for_chapter_seven(summary):
    """Whatever made it work, it was not thermal activation over that gap."""
    for row in summary['green_at_megabar']['rows']:
        assert row['hot_band'][90.0] < 1e-5
        assert row['hot_band'][50.0] < 1e-9


def test_the_bound_puts_the_optimum_red_of_the_frozen_answer(summary):
    x3 = summary['green_at_megabar']
    assert x3['required_red_shift_at_120GPa_meV'] > 90.0
    assert x3['lambda_opt_lower_bound_nm'] > 455.0
    assert x3['lambda_opt_lower_bound_nm'] > x3['frozen_lambda_opt_nm']


# --- X4 -------------------------------------------------------------------

def test_the_flat_culet_correction_alone_would_break_the_retrodiction(summary):
    x4 = summary['geometry']
    assert x4['g_micropillar'] == pytest.approx(1.0)
    assert x4['g_flat_culet'] == pytest.approx(0.564, abs=0.002)
    assert x4['crossover_at_flat_culet_GPa'] > 85.0


def test_band_position_and_anvil_factor_are_degenerate(summary):
    """Lower g fits better, but only by asking the anvil to favour blue."""
    rows = {round(row['g'], 3): row for row in summary['geometry']['rows']}
    assert rows[0.564]['chi2'] < rows[1.0]['chi2']
    assert rows[0.564]['transmission_ratio'] > 2.5
    assert rows[1.0]['transmission_ratio'] == pytest.approx(1.12, abs=0.05)


def test_requiring_a_physical_anvil_pins_the_band_position(summary):
    """No g below 1 keeps T(450)/T(532) at or under unity."""
    for row in summary['geometry']['rows']:
        if row['g'] < 0.99:
            assert row['transmission_ratio'] > 1.0


# --- X5 -------------------------------------------------------------------

def test_the_residual_declines_with_pressure(summary):
    x5 = summary['residual_trend']
    assert x5['slope_per_GPa'] < 0.0
    assert 1.5 < x5['slope_sigma'] < 3.0


def test_but_the_linear_form_is_not_the_right_one(summary):
    """It extrapolates to an anvil that transmits blue better than green."""
    assert summary['residual_trend']['ambient_extrapolation'] > 1.0


def test_the_improvement_is_suggestive_not_decisive(summary):
    assert 2.0 < summary['residual_trend']['delta_chi2'] < 9.0
