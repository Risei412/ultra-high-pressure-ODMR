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
