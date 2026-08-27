"""Tests for the Ho-integrated conditional ODMR sensitivity layer."""
import numpy as np
import pytest

from ho_odmr_sensitivity import (HoIntegratedODMRModel, ODMRResponse,
                                 optical_limit_summary)


def test_optical_limit_reproduces_published_curve_optimum():
    summary = optical_limit_summary(120.0)
    assert summary['optimum_nm'] == pytest.approx(440.65, abs=0.1)
    assert summary['penalty_457'] == pytest.approx(1.04494, rel=2e-4)


def test_green_penalty_and_blue_green_advantage_are_consistent():
    summary = optical_limit_summary(120.0)
    assert summary['penalty_532'] == pytest.approx(12.527, rel=2e-3)
    assert summary['advantage_457_over_532'] == pytest.approx(
        11.9893, rel=2e-4)


def test_contrast_threshold_is_inverse_optical_advantage():
    summary = optical_limit_summary(120.0)
    assert summary['minimum_C457_over_C532_for_457_to_win'] == pytest.approx(
        1.0 / summary['advantage_457_over_532'], rel=1e-10)


def test_nonoptical_response_can_reverse_the_ranking():
    def contrast(lam, pressure):
        lam = np.asarray(lam, float)
        return np.where(np.isclose(lam, 457.0), 0.05, 1.0)

    model = HoIntegratedODMRModel(response=ODMRResponse(contrast=contrast))
    assert model.pair_advantage(457.0, 532.0, 120.0) < 1.0


def test_saturation_flattens_but_does_not_move_optimum_when_other_factors_flat():
    model = HoIntegratedODMRModel()
    assert model.optimum(120.0, relative_power=1e-6) == pytest.approx(
        model.optimum(120.0, relative_power=10.0), abs=0.1)
    assert model.penalty(457.0, 120.0, relative_power=10.0) < 1.01


def test_invalid_response_factor_is_rejected():
    model = HoIntegratedODMRModel(
        response=ODMRResponse(charge_yield=lambda lam, p: 0.0))
    with pytest.raises(ValueError, match='charge_yield'):
        model.sensitivity(457.0, 120.0)
