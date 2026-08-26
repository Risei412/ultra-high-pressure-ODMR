"""
The two numbers a run needs, and the property that makes them trustworthy:
step 1 must reproduce the frozen anchor when handed the frozen ZPL, and step 2
must be expressible as a ratio, so that no absolute intensity is required.
"""
import numpy as np
import pytest

from experiment_plan import (dE120_from_zpl, lambda_opt_from_zpl,
                             line_penalties, zpl_range_for_line,
                             power_working_point, analyse_power_sweep,
                             u_nonlinear, u_optimal, ZPL_AMBIENT)
from nv_model import NVModel, eV2nm
from nv_model_power import NVModelPower


ANCHOR_ZPL = float(eV2nm(ZPL_AMBIENT + 0.400))      # 529 nm, the frozen anchor


# ---------------------------------------------------------------- step 1
def test_the_frozen_anchor_round_trips():
    assert dE120_from_zpl(ANCHOR_ZPL) == pytest.approx(0.400, abs=1e-3)
    assert lambda_opt_from_zpl(ANCHOR_ZPL) == pytest.approx(475.5, abs=0.3)


def test_inversion_works_away_from_120_GPa():
    """The ZPL can be read at whatever pressure the cell happens to be at."""
    m = NVModel(dE120=0.372)
    for P in (60.0, 90.0, 120.0):
        zpl = float(eV2nm(m.ZPL(P)))
        assert dE120_from_zpl(zpl, P) == pytest.approx(0.372, abs=2e-3)


def test_a_redder_zpl_pushes_the_optimum_red():
    """A (111) culet red-shifts the axial group; the recommendation must follow."""
    assert (lambda_opt_from_zpl(545.) > lambda_opt_from_zpl(529.)
            > lambda_opt_from_zpl(512.))


def test_473_is_the_best_commercial_line_at_the_anchor():
    best, penalty = line_penalties(ANCHOR_ZPL)[0]
    assert best == 473.0
    assert penalty < 1.01


def test_the_go_no_go_window_contains_the_anchor():
    lo, hi = zpl_range_for_line(473.0)
    assert lo < ANCHOR_ZPL < hi
    assert lo == pytest.approx(512, abs=3)
    assert hi == pytest.approx(541, abs=3)


def test_488_covers_the_redder_half():
    lo, hi = zpl_range_for_line(488.0)
    assert lo > zpl_range_for_line(473.0)[0]      # 488 is the redder choice
    assert hi > zpl_range_for_line(473.0)[1]


# ---------------------------------------------------------------- step 2
def test_the_recipe_is_a_ratio():
    """Neither end of the recipe may need an absolute intensity."""
    wp = power_working_point(473.0)
    assert wp['ratio'] == pytest.approx(wp['u_star'] / wp['u_knee'])
    assert 1.5 < wp['ratio'] < 2.5


def test_the_working_point_is_where_the_model_puts_it():
    assert u_nonlinear(473.0) == pytest.approx(0.060, abs=0.01)
    assert u_optimal(473.0) == pytest.approx(0.114, abs=0.02)


def test_the_knee_is_below_the_optimum():
    """If it were the other way round the recipe would read backwards."""
    wp = power_working_point(473.0)
    assert wp['u_knee'] < wp['u_star']


def test_a_synthetic_sweep_recovers_its_own_knee():
    mp = NVModelPower()
    u = np.logspace(-2.2, 0.4, 30)
    R = np.array([float(np.asarray(mp.eta_lambda_u(473., 120., uu)[2]))
                  for uu in u])
    got = analyse_power_sweep(u, R, 473., 120.)
    assert got['knee_power'] == pytest.approx(got['u_knee'], rel=0.15)
    assert got['operating_power'] == pytest.approx(got['u_star'], rel=0.15)


def test_a_sweep_that_never_saturates_says_so():
    mp = NVModelPower()
    u = np.logspace(-4, -3, 8)
    R = np.array([float(np.asarray(mp.eta_lambda_u(473., 120., uu)[2]))
                  for uu in u])
    assert np.isnan(analyse_power_sweep(u, R, 473., 120.)['knee_power'])
