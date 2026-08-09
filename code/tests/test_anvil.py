"""
Tests for the culet-orientation model (anvil_orientation.py).

Run from the `code/` directory:   python -m pytest tests/ -q
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import anvil_orientation as ao                                  # noqa: E402

P, T = 120.0, 5.0


def test_hydrostatic_stress_gives_no_splitting_and_one_line():
    """t = 0: every family sees pure hydrostatic stress, so E = 0 and D is common."""
    for culet in ('[111]', '[100]', '[110]'):
        fam = ao.families(P, 0.0, culet)
        assert all(f['E'] == pytest.approx(0.0, abs=1e-9) for f in fam)
        assert np.allclose([f['D'] for f in fam], fam[0]['D'])


def test_111_axial_family_has_exactly_zero_E():
    """The symmetry argument for a [111] culet: the axial NV keeps C3v."""
    fam = ao.families(P, T, '[111]')
    axial = [f for f in fam if f['theta_deg'] < 1e-6]
    assert len(axial) == 1
    assert axial[0]['E'] == pytest.approx(0.0, abs=1e-9)
    assert axial[0]['dEdt'] == pytest.approx(0.0, abs=1e-6)


def test_111_off_axis_families_are_degenerate_with_E_from_c():
    fam = [f for f in ao.families(P, T, '[111]') if f['theta_deg'] > 1.0]
    assert len(fam) == 3
    assert np.allclose([f['E'] for f in fam], fam[0]['E'])
    assert fam[0]['E'] == pytest.approx(4.0 / 3.0 * abs(ao.C_) * T, rel=1e-9)


def test_100_families_are_all_equivalent_with_E_from_b():
    fam = ao.families(P, T, '[100]')
    assert np.allclose([f['theta_deg'] for f in fam], 54.7356, atol=1e-3)
    assert np.allclose([f['E'] for f in fam], fam[0]['E'])
    assert fam[0]['E'] == pytest.approx(2.0 * abs(ao.B_) * T, rel=1e-9)


def test_100_decouples_hydrostatic_and_deviatoric_stress():
    """No shear in the crystal frame => D reads mean stress only, E reads t only."""
    for f in ao.families(P, T, '[100]'):
        assert f['dDdt'] == pytest.approx(0.0, abs=1e-6)
        assert f['dEdt'] == pytest.approx(2.0 * abs(ao.B_), rel=1e-6)


def test_the_E_zero_line_is_the_most_stress_sensitive_one():
    """The central surprise: E = 0 does not mean stress-insensitive."""
    lines = ao.odmr_lines(P, T, '[111]', dt=0.0)
    axial = max(lines, key=lambda L: abs(L['dnudt']))
    assert axial['dnudt'] == pytest.approx(-2.0 * ao.A2, rel=1e-6)     # 7.40 MHz/GPa
    for L in ao.odmr_lines(P, T, '[100]', dt=0.0):
        assert abs(L['dnudt']) < abs(axial['dnudt'])


def test_total_photon_rate_is_the_same_for_both_culets():
    """R does not drive the comparison: contrast and linewidth do."""
    w111 = sum(f['weight'] for f in ao.families(P, T, '[111]'))
    w100 = sum(f['weight'] for f in ao.families(P, T, '[100]'))
    assert w111 == pytest.approx(w100, rel=1e-12)
    assert w111 == pytest.approx(8.0 / 3.0, rel=1e-12)


def test_axial_line_never_beats_100():
    """Using the E = 0 line is worse than a [100] culet at every dt."""
    for dt in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0):
        lines = ao.odmr_lines(P, T, '[111]', dt)
        axl = max(lines, key=lambda L: L['dnudt'])
        e_ax = axl['width'] / (axl['contrast'] * np.sqrt(axl['total_pl']))
        assert e_ax > ao.best_line(P, T, '[100]', dt)[0]


def test_crossover_scales_linearly_with_the_intrinsic_linewidth():
    from scipy.optimize import brentq
    xs = []
    for dnu0 in (2.0, 5.0, 10.0, 20.0):
        xs.append(brentq(lambda d: ao.eta_ratio(P, T, d, dnu0=dnu0) - 1.0, 0.01, 60.0) / dnu0)
    assert np.allclose(xs, xs[0], rtol=1e-3)
    assert xs[0] == pytest.approx(0.42, abs=0.02)


def test_eta_ratio_is_independent_of_the_mean_deviatoric_stress():
    """Only the SPREAD dt broadens the lines; t itself just moves them."""
    for t in (1.0, 5.0, 20.0):
        assert ao.eta_ratio(P, t, 1.0) == pytest.approx(ao.eta_ratio(P, 5.0, 1.0), rel=1e-9)
