"""Regression tests for the 120--400 GPa micro-NV decision figures."""

import sys
from pathlib import Path

import pytest


CODE = Path(__file__).resolve().parents[1]
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import fig_microNV_wavelength_decision as low  # noqa: E402
import fig_microNV_wavelength_decision_300_400gpa as high  # noqa: E402


def test_120_200_conditional_ranges():
    assert low.alpha_range_120() == pytest.approx((466.08, 488.10), abs=0.06)
    assert low.alpha_range_200() == pytest.approx((430.01, 514.30), abs=0.06)


@pytest.mark.parametrize(
    "pressure, expected",
    [(300.0, (380.37, 509.96)), (400.0, (346.12, 508.39))],
)
def test_300_400_model_envelopes(pressure, expected):
    assert high.wavelength_envelope(pressure) == pytest.approx(expected, abs=0.06)


def test_488_is_inside_every_reported_envelope():
    ranges = [low.alpha_range_120(), low.alpha_range_200()]
    ranges.extend(high.wavelength_envelope(P) for P in high.PRESSURES)
    assert all(lo <= 488.0 <= hi for lo, hi in ranges)

