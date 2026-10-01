"""
test_freeze_a4_handoff.py
-------------------------
Freeze tests for Addendum A4 of ``docs/theory/theory_freeze_v3_ho_integrated.md``
-- the [111] axisymmetric-stress layer and the experimental handoff.

Four things are locked here:

  1. THE PRACTICAL BAND.  The 5% and 2% optical-limit tolerance bands and the
     per-line penalty table that the experiment is designed against (A4.1).

  2. THE FIXED LINE.  445 nm leaves 4.6% of the 5% invariance budget; 457 nm
     leaves 0.48%.  That margin, not the raw penalty, is why the criterion is
     evaluated on 445 nm (A4.2).  If a future edit inverts this, the frozen
     protocol is void and this test says so.

  3. THE A2 IDENTITY.  A tolerance factor f on sensitivity is the level set
     A = A_max/f^2, i.e. the A2 degenerate doublet at I/I_c = f^2 (A4.3).

  4. THE ALPHA IDENTITY.  d(lambda_opt)/d(alpha) = 0 is a property of the
     model, not a measurement: alpha never reaches the optical kernel (A4.0).
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ho_odmr_sensitivity import HoIntegratedODMRModel          # noqa: E402
from nv_model_111 import NV111AxialModel                       # noqa: E402
from theory_a1_generalization import Kernel                    # noqa: E402


TOLERANCE = 1.05          # the pre-registered invariance criterion
FIXED_LINE = 445.0        # A4.2
SECONDARY_LINE = 457.0


@pytest.fixture(scope='module')
def model():
    return HoIntegratedODMRModel()


@pytest.fixture(scope='module')
def kernel():
    return Kernel()


# ---------------------------------------------------------------------------
# 1. the practical band (A4.1)
# ---------------------------------------------------------------------------

def test_reported_optimum_is_441_nm(model):
    """440.65 nm is the interpolant maximum; 441 nm is what may be reported."""
    assert model.optimum(120.0) == pytest.approx(440.65, abs=0.01)
    assert NV111AxialModel(alpha=0.56).reported_lambda_opt(120.0) == 441.0


@pytest.mark.parametrize('factor,lo,hi', [(1.05, 426.43, 457.90),
                                          (1.02, 432.04, 452.34)])
def test_tolerance_bands_are_frozen(kernel, factor, lo, hi):
    band_lo, band_hi = kernel.tolerance_band(factor)
    assert band_lo == pytest.approx(lo, abs=0.02)
    assert band_hi == pytest.approx(hi, abs=0.02)


# The A4.1 table.  Penalties are 1/sqrt(relative absorbed rate).
FROZEN_PENALTIES = {
    405.0: 1.30601,
    441.0: 1.00030,
    445.0: 1.00379,
    457.0: 1.04494,
    473.0: 1.20538,
    488.0: 1.53772,
    505.0: 3.02145,
    532.0: 12.52813,
}


@pytest.mark.parametrize('line', sorted(FROZEN_PENALTIES))
def test_handoff_penalty_table_is_frozen(model, line):
    optimum_rate = model.absorbed_photon_proxy(model.optimum(120.0), 120.0)
    rate = model.absorbed_photon_proxy(line, 120.0)
    penalty = float(np.sqrt(optimum_rate / rate))
    assert penalty == pytest.approx(FROZEN_PENALTIES[line], rel=2e-4)


def test_every_line_inside_the_band_is_inside_the_band(model, kernel):
    """Self-consistency: the table and the band edges are the same object."""
    lo, hi = kernel.tolerance_band(TOLERANCE)
    for line, penalty in FROZEN_PENALTIES.items():
        assert (lo <= line <= hi) == (penalty <= TOLERANCE), line


def test_the_zpl_line_is_a_candidate_the_earlier_lists_omitted(model):
    """E2.1: 514.46 nm ties 473 nm on penalty and sits on real extracted data."""
    optimum_rate = model.absorbed_photon_proxy(model.optimum(120.0), 120.0)
    zpl = float(np.sqrt(optimum_rate
                        / model.absorbed_photon_proxy(514.46, 120.0)))
    assert zpl == pytest.approx(1.2006, rel=5e-3)
    assert abs(zpl - FROZEN_PENALTIES[473.0]) < 0.01


# ---------------------------------------------------------------------------
# 2. why the fixed line is 445 nm and not 457 nm (A4.2)
# ---------------------------------------------------------------------------

def _margin(model, line):
    """Budget left over after the hydrostatic baseline has been paid for."""
    optimum_rate = model.absorbed_photon_proxy(model.optimum(120.0), 120.0)
    penalty = float(np.sqrt(optimum_rate
                            / model.absorbed_photon_proxy(line, 120.0)))
    return TOLERANCE / penalty


def test_445_leaves_room_for_the_invariance_test(model):
    assert _margin(model, FIXED_LINE) == pytest.approx(1.0460, rel=1e-3)
    assert _margin(model, FIXED_LINE) - 1.0 > 0.04


def test_457_does_not_leave_room_for_the_invariance_test(model):
    """
    A4.2: at 457 nm the hydrostatic baseline alone spends x1.04494 of a x1.05
    budget.  The criterion would then fail on any real alpha dependence
    regardless of whether the invariance holds, so 457 nm cannot be the test
    line.  This is the whole reason the frozen protocol moved to 445 nm.
    """
    assert _margin(model, SECONDARY_LINE) == pytest.approx(1.0048, rel=1e-3)
    assert _margin(model, SECONDARY_LINE) - 1.0 < 0.01
    assert _margin(model, FIXED_LINE) > 5.0 * _margin(model, SECONDARY_LINE) - 4.0


def test_the_fixed_line_costs_less_than_half_a_percent(model):
    assert FROZEN_PENALTIES[FIXED_LINE] < 1.005


# ---------------------------------------------------------------------------
# 3. the band is the A2 level set (A4.3)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('factor', [1.02, 1.05])
def test_tolerance_band_equals_the_a2_doublet_at_power_ratio_f_squared(
        kernel, factor):
    """
    eta <= f*eta_min  <=>  A >= A_max/f^2, so the band edges ARE the P3(c)
    degenerate doublet at I/I_c = f^2.  An identity, and it must hold to the
    node resolution of the kernel.
    """
    lo, hi = kernel.tolerance_band(factor)
    doublet = NV111AxialModel(alpha=0.56).finite_power_optima(
        120.0, factor ** 2)
    assert len(doublet) == 2
    assert doublet[0] == pytest.approx(lo, abs=0.02)
    assert doublet[1] == pytest.approx(hi, abs=0.02)


def test_the_doublet_straddles_the_optimum(kernel):
    lo, hi = kernel.tolerance_band(TOLERANCE)
    assert lo < kernel.lam_abs < hi


# ---------------------------------------------------------------------------
# 4. alpha invariance is an identity, not a result (A4.0)
# ---------------------------------------------------------------------------

def test_alpha_never_reaches_the_optical_kernel():
    """
    Reporting rule 8.  The zero shift over alpha = 0.5-0.7 is the only value
    the code can return, because the same hydrostatic 120 GPa kernel is used
    at every alpha.  Locked so it is never re-reported as a finding.
    """
    optima = [NV111AxialModel(alpha=a).lambda_opt(120.0)
              for a in np.linspace(0.5, 0.7, 21)]
    assert np.ptp(optima) == 0.0
    assert all(NV111AxialModel(alpha=a).is_orientation_resolved_prediction
               is False for a in (0.5, 0.6, 0.7))


def test_stress_bookkeeping_over_the_handoff_alpha_range():
    """A4.0: P_mean = 80-96 GPa and q = 60-36 GPa at sigma_zz = 120 GPa."""
    for alpha, mean, dev in ((0.5, 80.0, 60.0), (0.7, 96.0, 36.0)):
        invariants = NV111AxialModel(alpha=alpha).stress_invariants(120.0)
        assert invariants['sigma_mean'] == pytest.approx(mean)
        assert invariants['sigma_dev_z_minus_xy'] == pytest.approx(dev)
        # [111] axisymmetric loading leaves no symmetry-breaking E component
        assert invariants['transverse_E1'] == pytest.approx(0.0, abs=1e-9)
        assert invariants['transverse_E2'] == pytest.approx(0.0, abs=1e-9)


def test_naive_mean_stress_substitution_is_not_a_small_shift(model):
    """
    A4.4: substituting P_mean into the Ho pressure coordinate -- which the
    frozen model refuses to do -- moves the optimum by 18 nm at alpha = 0.7 and
    branch-exchanges at alpha = 0.5.  Locked because it is the justification
    for keeping 473/514/532 nm in the scan.
    """
    assert model.optimum(96.0) == pytest.approx(458.70, abs=0.1)
    assert model.optimum(80.0) > 500.0          # branch exchange, Theorem X
    assert model.optimum(80.0) == pytest.approx(540.70, abs=0.1)


if __name__ == '__main__':
    raise SystemExit(pytest.main([os.path.abspath(__file__), '-q']))
