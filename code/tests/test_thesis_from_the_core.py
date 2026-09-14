"""Tests for the core-language reading of Bhattacharyya's chapter 6."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from theory_a1_generalization import Kernel  # noqa: E402
from thesis_from_the_core import (  # noqa: E402
    anvil_sensitivity, anvil_threshold, metric_bias, prefactor_invariance,
    report, sign_blindness,
)


@pytest.fixture(scope='module')
def summary():
    return report()


# ---- N1 -------------------------------------------------------------------

def test_eta_cannot_see_the_sign_of_the_contrast(summary):
    """A -5% resonance is exactly the magnetometer a +5% one is."""
    values = summary['sign_blindness']
    assert values[0.05] == pytest.approx(values[-0.05], rel=1e-12)
    assert values[0.02] == pytest.approx(values[-0.02], rel=1e-12)


def test_smaller_contrast_is_worse_regardless_of_sign(summary):
    values = summary['sign_blindness']
    assert values[0.02] < values[0.05]


# ---- N2 -------------------------------------------------------------------

def test_isc_prefactor_leaves_the_stationary_point_alone(summary):
    stars = list(summary['prefactor'].values())
    assert max(stars) == pytest.approx(min(stars), rel=1e-12)


# ---- N3 -------------------------------------------------------------------

def test_thesis_metric_is_biased_against_the_bluer_line(summary):
    for row in summary['metric_bias']:
        assert row['bias'] <= 1.0


def test_the_bias_grows_with_power(summary):
    biases = [row['bias'] for row in summary['metric_bias']]
    assert biases == sorted(biases, reverse=True)
    assert biases[-1] < 0.96


def test_the_two_metrics_disagree_about_the_winner_near_the_crossover(summary):
    """At I/I_c = 2 eta says blue wins and their SNR says green does."""
    row = next(r for r in summary['metric_bias'] if r['power'] == 2.0)
    assert row['inverse_eta_ratio'] > 1.0
    assert row['thesis_ratio'] < 1.0


# ---- N4 -------------------------------------------------------------------

def test_anvil_transmission_moves_the_optimum(summary):
    """Corrects the scope document: the ladder does NOT survive unchanged."""
    rows = summary['anvil']
    assert rows[None]['lambda_opt'] == pytest.approx(Kernel().lam_abs, abs=0.1)
    assert rows[200.0]['lambda_opt'] > 500.0
    assert rows[None]['first_rung'] != pytest.approx(
        rows[200.0]['first_rung'], rel=0.05)


def test_the_answer_is_fragile_to_a_modest_attenuation(summary):
    threshold = summary['anvil_threshold']
    assert threshold['decade_nm'] == pytest.approx(321.0, abs=5.0)
    assert threshold['loss_factor'] < 2.0


def test_a_74_nm_jump_is_available(summary):
    rows = summary['anvil']
    jump = rows[100.0]['lambda_opt'] - rows[None]['lambda_opt']
    assert jump > 70.0
