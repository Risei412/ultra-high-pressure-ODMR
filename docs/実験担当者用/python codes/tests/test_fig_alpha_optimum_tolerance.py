"""Regression tests for the conditional alpha-tolerance figure."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import fig_alpha_optimum_tolerance as figure  # noqa: E402
import fig_style  # noqa: E402


def test_optima_and_common_five_percent_band():
    data = figure.alpha_band(n=10)
    assert data["optimum_nm"][0] == pytest.approx(501.50, abs=0.02)
    assert data["optimum_nm"][-1] == pytest.approx(446.10, abs=0.02)
    low, high = figure.common_band(data)
    assert low == pytest.approx(487.28, abs=0.08)
    assert high == pytest.approx(492.61, abs=0.08)
    assert low < 488.0 < high

    axial = figure.alpha_band(
        n=10, deviatoric_ratio=figure.HILBERER_RATIO_EQUAL_AXIAL_STRESS)
    axial_low, axial_high = figure.common_band(axial)
    assert axial_low == pytest.approx(470.4, abs=0.1)
    assert axial_high == pytest.approx(483.3, abs=0.1)


def test_no_line_is_five_percent_safe_across_both_normalisations():
    wavelengths = figure.np.linspace(460.0, 505.0, 451)
    penalty = figure.worst_case_penalty(wavelengths)
    assert penalty.min() > 1.05


def test_figure_renders(tmp_path, monkeypatch):
    monkeypatch.setattr(figure, "OUT_PNG", tmp_path / "figure.png")
    monkeypatch.setattr(figure, "OUT_PDF", tmp_path / "figure.pdf")
    figure.main()
    assert (tmp_path / "figure.png").stat().st_size > 20_000
    assert (tmp_path / "figure.pdf").stat().st_size > 5_000
    width = fig_style.verify_width(tmp_path / "figure.pdf",
                                   fig_style.FULL_WIDTH_IN)
    assert width == pytest.approx(fig_style.FULL_WIDTH_IN * 72.0, abs=0.75)
