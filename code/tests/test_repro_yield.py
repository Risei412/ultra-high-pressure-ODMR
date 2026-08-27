"""The measured-yield comparison: does it load, and does it say what it says."""
import numpy as np
import pytest

from nv_model import NVModel
from repro_yield import (load, predict, compare, fit_dE120, peak_pressure,
                         _observed_peak, LINES)


@pytest.fixture(scope='module')
def data():
    return load()


def test_both_series_present_and_disjoint_in_pressure(data):
    p532 = data['expt532'][0]
    p457 = data['expt457'][0]
    assert len(p532) == 32 and len(p457) == 29
    # the green run stops where the blue one starts, which is why one figure
    # covers 0-120 GPa at all
    assert p532.max() == pytest.approx(p457.min(), abs=0.1)


def test_yields_are_positive_and_on_the_plotted_scale(data):
    for lam in LINES:
        y = data['expt%d' % lam][1]
        assert np.all(y > 0) and y.max() < 6.5


def test_each_series_has_an_interior_maximum(data):
    """The whole content of the figure: the band sweeps PAST each line."""
    for lam, lo, hi in ((532.0, 15.0, 40.0), (457.0, 70.0, 95.0)):
        P, y = data['expt%d' % lam]
        assert lo < _observed_peak(P, y) < hi


def test_scale_is_free_so_only_shape_is_tested(data):
    """Doubling a series must not change its residual."""
    m = NVModel(T=300.0)
    P, y = data['expt532']
    a = compare(m, data)['532']
    doubled = dict(data, expt532=(P, 2.0 * y))
    assert compare(m, doubled)['532'] == pytest.approx(a, rel=1e-9)


def test_frozen_model_fits_green_but_not_blue(data):
    c = compare(NVModel(T=300.0), data)
    assert c['532'] < 0.20
    assert c['457'] > 0.30


def test_frozen_blue_branch_never_turns_over_but_the_data_do(data):
    frozen = NVModel(T=300.0)
    assert peak_pressure(frozen, 457.0) > 115.0
    assert _observed_peak(*data['expt457']) < 95.0


def test_refit_moves_one_anchor_and_fixes_both_branches(data):
    dE, m, rms = fit_dE120(data, T=300.0)
    assert 0.52 < dE < 0.58
    assert rms < 0.16
    c = compare(m, data)
    assert c['532'] < 0.16 and c['457'] < 0.16
    assert peak_pressure(m, 457.0) == pytest.approx(
        _observed_peak(*data['expt457']), abs=5.0)


def test_refit_pulls_the_optimum_blue(data):
    _, m, _ = fit_dE120(data, T=300.0)
    frozen = NVModel(T=300.0)
    assert frozen.lambda_opt(120) == pytest.approx(475.5, abs=0.5)
    assert 440.0 < m.lambda_opt(120) < 460.0


def test_conclusion_does_not_depend_on_the_temperature_assumed(data):
    """The measurement was at 90 K; the frozen numbers are at 300 K."""
    warm = fit_dE120(data, T=300.0)[0]
    cold = fit_dE120(data, T=90.0)[0]
    assert abs(warm - cold) < 0.02


def test_dropping_the_collection_factor_makes_the_fit_worse(data):
    """eta_col is not a fudge: leaving it out degrades both branches."""
    m = NVModel(T=300.0)
    with_col = compare(m, data, collection=True)['pooled']
    without = compare(m, data, collection=False)['pooled']
    assert without > with_col


def test_predict_is_shape_preserving():
    m = NVModel(T=300.0)
    assert predict(m, 532.0, [0.0, 60.0, 120.0]).shape == (3,)
