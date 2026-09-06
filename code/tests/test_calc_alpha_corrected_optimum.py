"""Regression tests for the conditional E5 anisotropic correction."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import calc_alpha_corrected_optimum as correction  # noqa: E402
import nv_model  # noqa: E402


def test_alpha_06_at_120_gpa():
    assert correction.corrected_optimum_nm(0.60) == pytest.approx(488.10, abs=0.02)


def test_hydrostatic_anchor_is_unchanged():
    assert correction.corrected_optimum_nm(1.0) == pytest.approx(440.65, abs=1e-9)


def test_calculation_does_not_mutate_legacy_model_constant():
    previous = nv_model.ZPL_DEV_K
    correction.corrected_optimum_nm(0.60)
    assert nv_model.ZPL_DEV_K == previous


# ---------------------------------------------------------------------------
# The additive anchor is a 120 GPa quantity (added 2026-09-02)
# ---------------------------------------------------------------------------

def test_anchor_pressure_matches_the_kernel_limit():
    """HO_HYDROSTATIC_OPTIMUM_NM was read off the reconstructed kernel, which
    refuses above 120 GPa, so the anchor pressure is not a free choice."""
    from ho_spectrum_model import HoPublishedSpectrumModel
    assert correction.HO_ANCHOR_PRESSURE_GPA == pytest.approx(120.0)
    kernel = HoPublishedSpectrumModel()
    assert kernel.lambda_opt(correction.HO_ANCHOR_PRESSURE_GPA) == pytest.approx(
        correction.HO_HYDROSTATIC_OPTIMUM_NM, abs=0.05)
    with pytest.raises(ValueError, match='limited to 0--120 GPa'):
        kernel.sigma_abs(2.5, 200.0)


@pytest.mark.parametrize('P', [60.0, 130.0, 200.0, 400.0])
def test_corrected_optimum_refuses_away_from_the_anchor(P):
    """Mixing a 120 GPa anchor with a differential taken elsewhere produced an
    optimum that moved RED with pressure.  It must now refuse instead."""
    with pytest.raises(ValueError, match='only defined there'):
        correction.corrected_optimum_nm(0.60, P)


def test_corrected_optimum_unchanged_at_the_anchor():
    assert correction.corrected_optimum_nm(0.60) == pytest.approx(488.10, abs=0.02)
    assert correction.corrected_optimum_nm(0.60, 120.0) == pytest.approx(
        488.10, abs=0.02)
    assert correction.corrected_optimum_nm(1.0) == pytest.approx(440.65, abs=1e-9)


def test_shift_is_the_decomposition_of_the_corrected_optimum():
    """corrected_optimum = anchor + shift, exactly, at the anchor pressure."""
    for alpha in (0.56, 0.60, 0.70, 0.95):
        assert (correction.HO_HYDROSTATIC_OPTIMUM_NM
                + correction.deviatoric_shift_nm(alpha)) == pytest.approx(
                    correction.corrected_optimum_nm(alpha), abs=1e-9)


def test_shift_is_defined_at_any_pressure_and_grows_with_it():
    """The shift alone is well defined away from the anchor, and it grows --
    which is exactly why freezing the reference while advancing the shift was
    wrong."""
    at120 = correction.deviatoric_shift_nm(0.50, 120.0)
    at200 = correction.deviatoric_shift_nm(0.50, 200.0)
    assert at120 == pytest.approx(60.8, abs=0.5)
    assert at200 == pytest.approx(71.4, abs=0.5)
    assert at200 > at120


def test_shift_vanishes_at_hydrostatic():
    for P in (120.0, 200.0):
        assert correction.deviatoric_shift_nm(1.0, P) == pytest.approx(0.0, abs=1e-9)


def test_shift_keeps_the_cancellation_floor():
    floor = correction.cancellation_alpha()
    with pytest.raises(ValueError, match='changes sign'):
        correction.deviatoric_shift_nm(floor * 0.9)
