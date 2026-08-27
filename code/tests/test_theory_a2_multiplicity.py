"""Regression tests for Addendum A2 (multiplicity ladder, gauge degeneracy).

A2 is the kernel-independent layer, so most of these assertions are about
exact structure rather than about numbers read off the Ho reconstruction.  The
120 GPa ladder is included because it is the worked example the experiment is
designed against.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from theory_a1_generalization import DATA_WINDOW, Kernel  # noqa: E402
from theory_a2_multiplicity import (  # noqa: E402
    GaugeResponse, a1_mechanism_table, critical_values, gauge_degeneracy,
    gauge_family, identifiability, ladder_is_gauge_invariant, match_transitions,
    multiplicity, plateau_widths,
)


@pytest.fixture(scope='module')
def kernel():
    return Kernel()


@pytest.fixture(scope='module')
def matched(kernel):
    return match_transitions(kernel)


# ------------------------------------------------------------- Theorem M

def test_every_observed_transition_is_predicted_by_a_critical_value(matched):
    """M2: the ladder is not fitted -- each step falls on a critical value."""
    assert len(matched['rows']) == 6
    for row in matched['rows']:
        assert row['matched'], row
        assert row['relative_error'] < 1e-3


def test_transition_powers_are_the_reciprocal_critical_levels(kernel, matched):
    """M1: I_k/I_c = A_max/A_k exactly, so the ladder needs no calibration."""
    for row in matched['rows']:
        assert row['observed_power_ratio'] == pytest.approx(
            row['predicted_power_ratio'], rel=1e-3)


def test_multiplicity_changes_by_the_morse_rule(matched):
    """M1: +2 at an interior maximum, -2 at a minimum, -1 at a window edge."""
    expected = {'max': +2, 'min': -2, 'edge-blue': -1, 'edge-red': -1}
    for row in matched['rows']:
        assert row['delta'] == expected[row['kind']], row


def test_the_120GPa_ladder(matched):
    """M2: the specific sequence the experiment is designed to see."""
    steps = [(round(row['observed_power_ratio'], 3), row['before'], row['after'])
             for row in matched['rows']]
    assert steps == [
        (1.441, 2, 4),   # ZPL at 514.46 nm enters
        (1.514, 4, 6),   # 475.55 nm local maximum enters
        (1.519, 6, 4),   # 474.47 nm local minimum annihilates a pair
        (1.865, 4, 3),   # blue member leaves the data window at 402 nm
        (3.608, 3, 5),   # 500.19 nm local maximum enters
        (4.299, 5, 3),   # 497.86 nm local minimum annihilates a pair
    ]


def test_the_six_fold_plateau_is_too_narrow_to_test(kernel):
    """M3: N = 6 survives only a 0.3 % power range; it must not be claimed."""
    rows = {row['multiplicity']: row for row in plateau_widths(kernel)}
    assert rows[6]['width_factor'] == pytest.approx(1.003, abs=0.002)
    # The plateaux the experiment can actually resolve.
    widths = [row['width_factor'] for row in plateau_widths(kernel)
              if np.isfinite(row['width_factor'])]
    assert max(widths) == pytest.approx(1.935, abs=0.02)


def test_zpl_is_the_first_rung_so_omitting_it_hides_the_effect(kernel):
    """M1: the lowest transition comes from the ZPL, not from the main band."""
    first = critical_values(kernel)[0]
    assert first.wavelength == pytest.approx(514.46, abs=0.05)
    assert first.power_ratio == pytest.approx(1.441, abs=0.005)


def test_multiplicity_is_two_just_above_the_critical_power(kernel):
    """M: below the first rung the classic A1 doublet is all there is."""
    assert multiplicity(kernel, 1.0 / 1.05, DATA_WINDOW) == 2
    assert multiplicity(kernel, 1.0 / 1.40, DATA_WINDOW) == 2


# ------------------------------------------------------------- Theorem G

def test_splitting_is_decided_by_the_single_exponent():
    """G0: A1's four-row mechanism table is the sign of E - 1."""
    table = {row['mechanism']: row for row in a1_mechanism_table()}
    assert table['rate saturation alone']['exponent'] == pytest.approx(1.0)
    assert table['rate saturation alone']['splits'] is False
    assert table['power broadening alone']['exponent'] == pytest.approx(1.0)
    assert table['power broadening alone']['splits'] is False
    assert table['saturation + broadening']['exponent'] == pytest.approx(2.0)
    assert table['contrast collapse alone']['exponent'] == pytest.approx(2.0)
    for name in ('saturation + broadening', 'contrast collapse alone'):
        assert table[name]['splits'] is True
        assert table[name]['rho_star'] == pytest.approx(1.0)


def test_all_three_mechanisms_reproduce_the_a1_numeric():
    """G0: E = 4 gives rho* = 1/3, matching A1's Gamma_p* = 0.3333."""
    row = {r['mechanism']: r for r in a1_mechanism_table()}['all three']
    assert row['exponent'] == pytest.approx(4.0)
    assert row['rho_star'] == pytest.approx(1.0 / 3.0, rel=1e-9)


def test_gauge_plane_members_share_phi_exactly():
    """G1: same E means the same sensitivity surface, to machine precision."""
    for row in gauge_degeneracy():
        assert row['exponent'] == pytest.approx(2.0)
        assert row['max_relative_phi_difference'] < 1e-12
        assert row['gamma_star'] == pytest.approx(1.0, rel=1e-9)


def test_gauge_equivalent_models_differ_in_the_measured_factors():
    """G2: what eta hides is visible in R, C and dnu separately."""
    rows = {row['name']: row for row in gauge_degeneracy()}
    first = rows['contrast collapse alone']
    second = rows['saturation + broadening']
    assert first['rate_at_gamma_1'] == pytest.approx(1.0)
    assert second['rate_at_gamma_1'] == pytest.approx(0.5)
    assert first['contrast_at_gamma_1'] == pytest.approx(0.5)
    assert second['contrast_at_gamma_1'] == pytest.approx(1.0)
    assert second['linewidth_at_gamma_1'] == pytest.approx(np.sqrt(2.0))


def test_identifiability_ladder():
    """G3: eta fixes only Phi; the rate and one more spectrum fix the rest."""
    result = identifiability()
    assert result['eta alone']['gap'] < 1e-12
    for key in ('+ PL rate', '=> gives G', '+ contrast'):
        assert result[key]['gap'] > 0.1


def test_ladder_is_invariant_across_the_gauge_plane(kernel):
    """G4: the multiplicity prediction survives not knowing the mechanism."""
    rows = ladder_is_gauge_invariant(kernel)
    assert len(rows) == len(gauge_family())
    assert all(row['identical_to_reference'] for row in rows)
    assert rows[0]['counts'] == [2, 4, 4, 3, 5, 3]


def test_exponent_below_one_never_splits():
    """G: E <= 1 has no interior maximum, at any power."""
    for exponent, response in ((0.5, GaugeResponse(c=0.0, s=0.5, w=0.0)),
                               (1.0, GaugeResponse(c=0.0, s=1.0, w=0.0))):
        assert response.exponent == pytest.approx(exponent)
        assert response.splits is False
        assert not np.isfinite(response.gamma_star())


def test_a2_reproduces_the_published_dreau_anchor():
    """G: A2's generalised rule at (c,s,w,n) = (1, 0, 1/2, 2) gives rho* = 1/5.

    Those exponents are read straight off Dreau et al.'s fitted CW-ODMR forms,
    so A2's splitting antecedent is measured rather than assumed.
    """
    response = GaugeResponse(c=1.0, s=0.0, w=0.5, n=2.0)
    assert response.exponent == pytest.approx(3.0)
    assert response.splits is True
    assert response.rho_star == pytest.approx(0.2, rel=1e-9)
    # And it really is the maximiser of Phi.
    star = response.gamma_star()
    assert response.phi(star) > response.phi(star / 3.0)
    assert response.phi(star) > response.phi(star * 3.0)


def test_pump_nonlinearity_defaults_to_linear():
    """G: n = 1 recovers A1's rho* = 1/(E-1), so A2 stays backward compatible."""
    linear = GaugeResponse(c=1.0, s=0.0, w=0.0)
    assert linear.n == pytest.approx(1.0)
    assert linear.rho_star == pytest.approx(1.0)
    assert linear.gamma_star() == pytest.approx(1.0)


def test_rho_star_formula():
    """G: rho* = 1/(E-1) for every member with E > 1."""
    for c, s, w in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.5), (0.5, 1.0, 0.25),
                    (1.0, 1.0, 0.5), (2.0, 0.0, 0.0)):
        response = GaugeResponse(c=c, s=s, w=w)
        if response.splits:
            assert response.gamma_star() == pytest.approx(
                1.0 / (response.exponent - 1.0), rel=1e-9)
            # Confirm it really is the maximiser of Phi.
            grid = response.gamma_star() * np.array([0.5, 1.0, 2.0])
            values = response.phi(grid)
            assert values[1] > values[0] and values[1] > values[2]
