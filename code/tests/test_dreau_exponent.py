"""Tests anchoring A2's splitting exponent in Dreau et al., PRB 84, 195204.

These lock in the claim that A2's second antecedent -- Phi has a unique
interior maximum -- is not assumed but measured, and that the ladder's first
rung sits at an ordinary operating power.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dreau_exponent import (  # noqa: E402
    OMEGA_R, contrast, first_rung_power, headroom, kappa, linewidth,
    optimum_pump, phi, rate, robustness_over_rabi, splitting_exponent,
)


def test_published_forms_give_E_equals_three():
    """D1: c = 1, s_exp = 0, w = 1/2 read off Dreau's Eqs. (10)-(12)."""
    split = splitting_exponent()
    assert split['E'] == pytest.approx(3.0)
    assert split['n'] == pytest.approx(2.0)
    assert split['E_times_n'] == pytest.approx(6.0)


def test_splitting_is_satisfied_with_margin():
    """D1: the criterion is E n > 1; the published model gives 6."""
    split = splitting_exponent()
    assert split['splits'] is True
    assert split['E_times_n'] > 5.0  # not marginal


def test_rho_star_closed_form():
    """D1/D2: rho* = 1/(E n - 1) = 1/5, confirmed numerically."""
    split = splitting_exponent()
    assert split['rho_star'] == pytest.approx(0.2)
    star = optimum_pump()
    assert star['rho_star_measured'] == pytest.approx(0.2, rel=1e-4)


def test_x_star_matches_sqrt_kappa_over_five():
    """D2: the interior maximum is at x* = sqrt(kappa/5)."""
    star = optimum_pump()
    assert star['x_star'] == pytest.approx(star['x_star_closed_form'], rel=1e-4)
    assert star['x_star'] == pytest.approx(np.sqrt(kappa() / 5.0), rel=1e-4)


def test_phi_really_peaks_there():
    """D2: Phi is larger at s* than on either side -- a genuine interior max."""
    s_star = optimum_pump()['s_star']
    assert phi(s_star) > phi(s_star / 3.0)
    assert phi(s_star) > phi(s_star * 3.0)


def test_rho_star_is_independent_of_the_microwave_setting():
    """D4: every Rabi frequency Dreau reports gives the same rho* = 1/5.

    The split condition is a property of the functional forms, not of the
    microwave tuning, so it cannot be tuned away.
    """
    rows = robustness_over_rabi()
    assert len(rows) == 5
    for row in rows:
        assert row['rho_star'] == pytest.approx(0.2, rel=1e-4)
    # s* itself does move with Omega_R, monotonically.
    star_values = [row['s_star'] for row in rows]
    assert star_values == sorted(star_values)


def test_first_rung_sits_at_an_ordinary_operating_power():
    """D3: the 120 GPa ladder's first rung is at ~10 % of saturation power."""
    rung = first_rung_power(0.6938)
    assert rung['I_over_Ic'] == pytest.approx(1.441, abs=0.002)
    assert rung['s_at_first_rung'] == pytest.approx(0.104, abs=0.005)
    # Well below saturation: this is not an exotic high-power regime.
    assert rung['s_at_first_rung'] < 0.2


def test_many_rungs_are_reachable_within_a_normal_power_sweep():
    """D3: s up to 10 P_sat covers I/I_c far beyond the highest rung (4.30)."""
    reach = headroom()['max_I_over_Ic']
    assert reach > 100.0
    assert reach > 4.30


def test_dreau_forms_reproduce_their_qualitative_statements():
    """Sanity: contrast falls and linewidth grows with optical power."""
    low, high = 0.02, 2.0
    assert contrast(low) > contrast(high)
    assert linewidth(low) < linewidth(high)
    assert rate(low) < rate(high)
