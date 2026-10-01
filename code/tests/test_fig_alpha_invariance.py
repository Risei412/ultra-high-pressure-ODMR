"""Regression tests for the experimental handoff alpha-invariance figure."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import fig_alpha_invariance as figure  # noqa: E402
import fig_style  # noqa: E402


def test_stress_changes_but_current_optimum_does_not():
    data = figure.alpha_sweep(n=5)
    assert data['sigma_xy_GPa'] == pytest.approx([60., 66., 72., 78., 84.])
    assert data['mean_GPa'] == pytest.approx([80., 84., 88., 92., 96.])
    assert data['axial_deviator_GPa'] == pytest.approx([60., 54., 48., 42., 36.])
    assert np.ptp(data['numerical_optimum_nm']) == pytest.approx(0.0, abs=1e-12)
    assert data['reported_optimum_nm'] == pytest.approx([441.] * 5)


def test_figure_renders_both_handoff_formats(tmp_path, monkeypatch):
    monkeypatch.setattr(figure, 'OUT_PNG', tmp_path / 'figure.png')
    monkeypatch.setattr(figure, 'OUT_PDF', tmp_path / 'figure.pdf')
    figure.main()
    assert (tmp_path / 'figure.png').stat().st_size > 20_000
    assert (tmp_path / 'figure.pdf').stat().st_size > 5_000
    width = fig_style.verify_width(tmp_path / 'figure.pdf',
                                   fig_style.FULL_WIDTH_IN)
    assert width == pytest.approx(fig_style.FULL_WIDTH_IN * 72.0, abs=0.75)
