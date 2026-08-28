"""Tests for the anvil-transmission bound taken from Fig. 6.3(b)."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from anvil_transmission import (  # noqa: E402
    crossover, fit_transmission, load, observed_ratio, optimum_under,
    predicted_ratio, report, verdict,
)


@pytest.fixture(scope='module')
def summary():
    return report()


def test_the_calibration_lands_on_round_pressures():
    """The check that the tick-stub calibration is right."""
    pressures = [row['pressure_GPa'] for row in load()[532.0]]
    for value, expected in zip(pressures, (2.9, 27.0, 42.0, 50.0, 60.0, 70.0)):
        assert value == pytest.approx(expected, abs=0.15)


def test_both_lines_are_digitised_at_the_same_pressures():
    data = load()
    for blue, green in zip(data[450.0], data[532.0]):
        assert blue['pressure_GPa'] == pytest.approx(green['pressure_GPa'])


def test_the_observed_ratio_rises_and_crosses_unity(summary):
    ratios = [row['ratio'] for row in summary['rows']]
    assert ratios[0] < 0.3
    assert ratios[-2] > 1.0
    assert 50.0 < summary['crossover'] < 60.0


def test_the_observed_crossover_matches_the_kernel(summary):
    """The kernel put it at 51.8 GPa without seeing this figure."""
    assert abs(summary['crossover'] - 51.8) < 4.0


def test_a_constant_anvil_factor_fits(summary):
    assert summary['fit']['chi2_per_dof'] < 2.0
    assert summary['fit']['n_points'] == 5


def test_the_fitted_transmission_is_consistent_with_no_attenuation(summary):
    fit = summary['fit']
    assert fit['transmission_ratio'] == pytest.approx(0.94, abs=0.03)
    assert fit['ci_low'] < 1.0 < fit['ci_high']


def test_the_optimum_survives_at_the_fit(summary):
    assert summary['verdict']['optimum_central'] == pytest.approx(440.65,
                                                                  abs=0.1)


def test_but_only_at_about_two_sigma(summary):
    """Not comfortable enough to close the question."""
    check = summary['verdict']
    assert 1.5 < check['sigma_above_threshold'] < 2.5
    assert check['optimum_minus_1sigma'] < 460.0
    assert check['optimum_minus_2sigma'] > 500.0


def test_no_attenuation_reproduces_the_frozen_answer():
    assert optimum_under(1e9) == pytest.approx(440.65, abs=0.1)


def test_predicted_ratio_is_monotone_in_pressure():
    values = [predicted_ratio(float(p)) for p in (10, 30, 50, 70, 90)]
    assert values == sorted(values)
