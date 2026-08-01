"""
Regression / sanity tests for the background-aware model (Part E).

Run from the `code/` directory:   python -m pytest tests/ -q
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import background as bg                                        # noqa: E402
from nv_model import NVModel                                   # noqa: E402
from nv_bg import (NVModelBG, optimum_wavelength,               # noqa: E402
                   tolerance_band, crossover_pressure)

LAM = np.linspace(402, 560, 1580)


# ------------------------------------------------------------ background ---
def test_shapes_normalised_at_532():
    for k in bg.CHANNELS:
        assert bg.g_channel(k, bg.LAM_REF) == pytest.approx(1.0)
    assert float(bg.g_total(bg.LAM_REF)) == pytest.approx(1.0)


def test_ruby_has_transmission_gap_in_the_blue():
    """Cr3+ absorbs at 410 and 555 nm, so 460-500 nm must be a minimum."""
    lam = np.linspace(402, 560, 800)
    gap = lam[int(np.argmin(bg.g_channel('ruby', lam)))]
    assert 460.0 < gap < 500.0
    assert bg.g_channel('ruby', gap) < 0.5 * bg.g_channel('ruby', 532.0)


def test_n3_penalises_the_near_uv():
    assert bg.g_channel('n3', 405.0) > bg.g_channel('n3', 457.0) > bg.g_channel('n3', 532.0)


def test_diamond_raman_never_reaches_the_detection_window():
    for lam, r1, r2, in1, in2 in bg.raman_report():
        assert not in1, f'1st-order Raman of {lam} nm at {r1:.0f} nm entered the window'
        assert not in2, f'2nd-order Raman of {lam} nm at {r2:.0f} nm entered the window'


# ------------------------------------------------------- consistency -------
def test_rho0_zero_reproduces_the_baseline_model():
    base, m = NVModel(), NVModelBG()
    P = np.linspace(0, 140, 71)
    for beams in ([(532.0, 1.0)], [(457.0, 1.0)], [(532.0, 0.5), (457.0, 0.5)]):
        e0 = base.eta(beams, P)[0]
        eb = m.eta_bg(beams, P, 0.0)[0]
        assert np.allclose(eb, e0, rtol=1e-12, atol=0.0)


def test_background_always_costs_sensitivity():
    """Monotonic increase in rho0 -- catches the sqrt(1+rho) sign error."""
    m = NVModelBG()
    rhos = np.logspace(-2, 3, 40)
    for lam in (457.0, 475.0, 532.0):
        e = np.array([float(m.eta_bg([(lam, 1.0)], 120.0, r)[0]) for r in rhos])
        assert np.all(np.diff(e) > 0)
        assert float(m.eta_bg([(lam, 1.0)], 120.0, 0.0)[0]) < e[0]


def test_penalty_is_exactly_sqrt_one_plus_rho():
    m = NVModelBG()
    for rho0 in (0.3, 3.0, 30.0):
        eta, _, R, _, B, rho = m.eta_bg([(475.0, 1.0)], 120.0, rho0)
        eta0 = m.eta([(475.0, 1.0)], 120.0)[0]
        assert float(rho) == pytest.approx(float(B / R))
        assert float(eta) == pytest.approx(float(eta0 * np.sqrt(1.0 + rho)))


# ---------------------------------------------------------- regressions ----
def test_baseline_regression_values():
    m = NVModelBG()
    assert optimum_wavelength(m, 120.0, 0.0, LAM) == pytest.approx(475.0, abs=1.0)
    assert optimum_wavelength(m, 100.0, 0.0, LAM) == pytest.approx(487.0, abs=1.5)
    assert crossover_pressure(m, 0.0) == pytest.approx(86.0, abs=1.5)


def test_optimum_moves_red_and_stays_inside_the_tolerance_band():
    """H2: background pushes lambda_opt to the red, but only marginally."""
    m = NVModelBG()
    lo0 = optimum_wavelength(m, 120.0, 0.0, LAM)
    lo_hi = optimum_wavelength(m, 120.0, 1e3, LAM)
    assert lo_hi > lo0                       # red shift
    assert lo_hi - lo0 < 10.0                # but small
    band0 = tolerance_band(m, 120.0, 0.0, LAM)
    assert band0[0] <= lo_hi <= band0[1]     # still within the zero-background 5% band


def test_blue_never_loses_at_120_GPa_however_bad_the_background():
    """The green/blue decision saturates: no rho0 flips it at the target pressure."""
    m = NVModelBG()
    for rho0 in (0.0, 1.0, 1e3, 1e6):
        lo = optimum_wavelength(m, 120.0, rho0, LAM)
        ratio = float(m.eta_bg([(lo, 1.0)], 120.0, rho0)[0]
                      / m.eta_bg([(532.0, 1.0)], 120.0, rho0)[0])
        assert ratio < 1.0
