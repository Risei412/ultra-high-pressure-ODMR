"""Regression tests for the Bhattacharyya-thesis retrodictions."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ho_spectrum_model import HoPublishedSpectrumModel  # noqa: E402
from thesis_crosscheck import (  # noqa: E402
    REFERENCE_NM, crossover, ratio, report, snr_ratio,
)


@pytest.fixture(scope='module')
def summary():
    return report()


@pytest.fixture(scope='module')
def model():
    return HoPublishedSpectrumModel()


def test_450_and_532_are_tied_where_the_thesis_measured(summary):
    """Fig. 6.3's null result is the prediction, not an accommodation."""
    assert summary['fig63_snr_ratio'] == pytest.approx(0.94, abs=0.03)


def test_the_450_crossover_lands_on_their_measurement_pressure(summary):
    assert summary['fig63_crossover'] == pytest.approx(51.8, abs=1.0)


def test_450_would_have_won_higher_up(model):
    """Half a hit: they compressed to ~70 GPa and still saw no advantage."""
    assert snr_ratio(model, 450.0, 70.0) > 1.8


def test_405_spans_four_orders_between_their_two_pressures(summary):
    low = summary['fig62_ratio_17']
    high = summary['fig62_ratio_100']
    assert low < 0.02
    assert high > 50.0
    assert high / low > 1e3


def test_405_crossover_is_above_both_thesis_pressures_but_one(summary):
    assert summary['fig62_crossover'] == pytest.approx(72.0, abs=1.5)


def test_532_survives_the_zpl_crossing(summary):
    """Doherty projected 532 nm to be precluded at ~66 GPa; Dai ran to 140."""
    assert summary['a532_fraction_66'] > 0.2
    assert summary['a532_fraction_120'] > 0.0


def test_532_absorption_is_never_negative_in_the_used_range(model):
    for pressure in range(0, 121, 10):
        assert ratio(model, REFERENCE_NM, float(pressure)) == 1.0


def test_crossover_is_monotone_in_wavelength(model):
    """Bluer lines take over later, because they start further off the band."""
    pressures = [crossover(model, lam) for lam in (473.0, 457.0, 450.0, 405.0)]
    assert pressures == sorted(pressures)
