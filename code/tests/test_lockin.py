"""
Tests for the modulated-detection simulation (lockin_sim.py).

Run from the `code/` directory:   python -m pytest tests/ -q
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import lockin_sim as ls                                          # noqa: E402


def test_noise_psd_normalisation_is_white_when_knee_is_zero():
    rng = np.random.default_rng(0)
    fs, n, S0 = 20_000.0, 1 << 16, 3e-8
    v = np.var([ls.technical_noise(n, fs, S0, 0.0, rng) for _ in range(8)])
    assert v == pytest.approx(S0 * fs / 2.0, rel=0.10)


def test_one_over_f_variance_matches_the_band_integral():
    """var = int S0 (1 + f_knee/f) df over [fs/n, fs/2] = S0[fs/2 + f_knee ln(n/2)]."""
    rng = np.random.default_rng(0)
    fs, n, S0, f_knee = 20_000.0, 1 << 16, 3e-8, 1000.0
    expected = S0 * (fs / 2.0 + f_knee * np.log(n / 2.0))
    v = np.mean([np.var(ls.technical_noise(n, fs, S0, f_knee, rng)) for _ in range(8)])
    assert v == pytest.approx(expected, rel=0.10)
    assert v > 1.8 * S0 * fs / 2.0                    # clearly above the white floor


def test_each_dwell_holds_a_whole_number_of_modulation_cycles():
    for f_mod in (100.0, 500.0, 2000.0, 5000.0):
        m, p = ls.dwell_samples(ls.FS, ls.T_TOTAL, ls.N_PTS, f_mod)
        assert p % 2 == 0                      # exact 50% duty cycle
        assert m % p == 0                      # exact DC rejection on demodulation
        assert m >= p


def test_noise_free_lineshapes_sit_on_the_true_centre():
    rng = np.random.default_rng(0)
    nu, dc = ls.simulate_sweep(0.0, rng, 'dc', tech_gain=0.0, shot=False)
    _, oo = ls.simulate_sweep(0.0, rng, 'onoff', tech_gain=0.0, shot=False)
    _, fm = ls.simulate_sweep(0.0, rng, 'fm', tech_gain=0.0, shot=False)
    assert nu[np.argmin(dc)] == pytest.approx(ls.D_TRUE, abs=1e-6)     # dip
    assert nu[np.argmax(oo)] == pytest.approx(ls.D_TRUE, abs=1e-6)     # peak
    i = int(np.argmin(np.abs(nu - ls.D_TRUE)))
    assert abs(fm[i]) < 1e-6 * np.max(np.abs(fm))                      # zero crossing
    assert fm[i - 3] * fm[i + 3] < 0                                   # dispersive


def test_fm_harmonic_is_odd_about_the_centre():
    d = np.linspace(-20.0, 20.0, 41)
    y = ls.fm_harmonic(ls.D_TRUE + d, ls.D_TRUE, ls.LINEWIDTH, 0.5)
    assert np.allclose(y, -y[::-1], atol=1e-12)


def test_all_estimators_are_unbiased_without_technical_noise():
    for scheme in ('dc', 'onoff', 'fm'):
        rng = np.random.default_rng(11)
        est = np.array([ls.estimate_D(0.0, rng, scheme, tech_gain=0.0) for _ in range(60)])
        est = est[np.isfinite(est)]
        bias = est.mean() - ls.D_TRUE
        assert abs(bias) < 3.0 * est.std() / np.sqrt(est.size) + 1e-3


def test_dc_degrades_with_the_knee_while_modulation_stays_flat():
    lo = ls.sigma_D(1.0, n_mc=60, seed=2)
    hi = ls.sigma_D(3000.0, n_mc=60, seed=2)
    assert hi['dc'] > 10.0 * lo['dc']                    # DC collapses
    assert hi['onoff'] < 2.0 * lo['onoff']               # modulation barely moves
    assert hi['fm'] < 2.0 * lo['fm']


def test_modulation_costs_duty_cycle_when_there_is_no_technical_noise():
    """The reason a break-even knee exists at all."""
    s = ls.sigma_D(0.0, n_mc=120, seed=4)
    assert s['dc'] < s['onoff']
    assert s['dc'] < s['fm']


def test_higher_modulation_frequency_helps_up_to_saturation():
    slow = ls.sigma_D(500.0, n_mc=100, seed=6, schemes=('onoff',), f_mod=105.0)['onoff']
    fast = ls.sigma_D(500.0, n_mc=100, seed=6, schemes=('onoff',), f_mod=2500.0)['onoff']
    assert fast < slow
