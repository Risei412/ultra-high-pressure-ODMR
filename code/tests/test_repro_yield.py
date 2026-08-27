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


# --- is it a wrong constant or a wrong theory? ------------------------------

SINGLE_KNOB = (('dE120', np.arange(0.40, 0.70, 0.01)),
               ('hw', np.arange(0.065, 0.14, 0.002)),
               ('S_slope', np.arange(1.5, 6.5, 0.10)),
               ('alpha', np.arange(0.6, 1.5, 0.02)))


def _best_single(data, name, values, T=300.0):
    best = None
    for v in values:
        m = NVModel(T=T, **{name: float(v)})
        c = compare(m, data)['pooled']
        if best is None or c < best[0]:
            best = (c, m)
    return best


@pytest.mark.parametrize('name,values', SINGLE_KNOB)
def test_every_admissible_knob_lands_the_optimum_near_450(data, name, values):
    """The answer must not depend on which input is blamed.

    If the structure of the model were wrong, no single input would reconcile
    both branches, and the ones that came closest would disagree about where
    the optimum is.  They agree.
    """
    rms, m = _best_single(data, name, values)
    assert rms < 0.16
    assert 440.0 < m.lambda_opt(120) < 460.0


def test_fixing_blue_does_not_cost_green(data):
    """No trade-off: the signature of a mis-set constant, not a wrong form."""
    frozen = compare(NVModel(T=300.0), data)
    _, refit, _ = fit_dE120(data, T=300.0)
    after = compare(refit, data)
    assert after['457'] < frozen['457'] / 2.0
    assert after['532'] <= frozen['532']


def test_the_green_advantage_is_not_pinned_by_these_data(data):
    """Honest limit: eta(532)/eta_opt varies by ~5x across the knobs that fit."""
    ratios = []
    for name, values in SINGLE_KNOB:
        m = _best_single(data, name, values)[1]
        opt = m.lambda_opt(120)
        e = lambda lam: float(np.asarray(m.eta_lambda(lam, 120)[0]))
        ratios.append(e(532.0) / e(opt))
    assert max(ratios) / min(ratios) > 3.0


def test_the_ZPL_slope_alone_cannot_rescue_the_blue_branch(data):
    """One knob that does NOT work, so the test above is not vacuous."""
    rms, _ = _best_single(data, 'slope0', np.arange(4.0e-3, 1.0e-2, 2e-4))
    assert rms > 0.20


def test_457_is_safe_under_every_model_that_fits(data):
    """The engineering decision does not depend on resolving the attribution."""
    models = [NVModel(T=300.0)]
    models += [_best_single(data, n, v)[1] for n, v in SINGLE_KNOB]
    for m in models:
        opt = m.lambda_opt(120)
        e = lambda lam: float(np.asarray(m.eta_lambda(lam, 120)[0]))
        assert e(457.0) / e(opt) < 1.15
        assert e(532.0) / e(opt) > 2.0        # and blue always beats green

