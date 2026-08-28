"""Regression tests for the N1 test against Fig. 6.7(b).  The result is null."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fig67_contrast_map import (  # noqa: E402
    BASELINE_WINDOWS, DEFAULT_WINDOW, PANEL_B_SHAPE, THESIS_PDF,
    baseline_sensitivity, bridge, contrast_map, load_panel_b, report,
    sign_weight_profile, tracks, verdict,
)

pytestmark = pytest.mark.skipif(not os.path.exists(THESIS_PDF),
                                reason='thesis PDF not present')


@pytest.fixture(scope='module')
def panel():
    return load_panel_b()


@pytest.fixture(scope='module')
def summary():
    return report()


# --- the figure inverts exactly -------------------------------------------

def test_the_palette_inverts_without_guessing(panel):
    """An indexed colour space means the colour map is IN the PDF."""
    indices, palette, exact = panel
    assert exact == 1.0
    assert len(palette) == 254
    assert indices.shape == PANEL_B_SHAPE


def test_the_palette_is_monotone_and_reversed(panel):
    """Index 0 is the bright end, so v = (253 - index)/253 rises with counts."""
    _, palette, _ = panel
    luminance = palette.astype(float) @ np.array([0.2126, 0.7152, 0.0722])
    assert luminance[0] > luminance[-1]
    assert np.corrcoef(np.arange(len(palette)), luminance)[0, 1] < -0.9


# --- the baseline check that catches a false positive ---------------------

def test_the_brightest_positive_feature_is_a_baseline_artefact(panel):
    """+0.49 at a 401-px window, +0.02 at 1201: the median dragged between dips."""
    sensitivity = baseline_sensitivity(panel[0])
    bright = sensitivity['bright top feature']
    assert bright[401] > 0.4
    assert bright[1201] < 0.05
    assert bright['stable'] is False


def test_the_other_positive_lines_survive_every_window(panel):
    sensitivity = baseline_sensitivity(panel[0])
    for name in ('positive line, mid', 'positive line, right'):
        values = [sensitivity[name][w] for w in BASELINE_WINDOWS]
        assert min(values) > 0.1
        assert max(values) - min(values) < 0.1 * max(values) + 0.02
        assert sensitivity[name]['stable'] is True


def test_the_deep_negative_line_is_not_a_baseline_effect(panel):
    sensitivity = baseline_sensitivity(panel[0])
    assert all(sensitivity['deep negative line'][w] < -0.6
               for w in BASELINE_WINDOWS)


# --- what the map contains ------------------------------------------------

def test_both_signs_are_present(summary):
    signs = summary['verdict']['signs']
    assert signs['negative'] >= 5
    assert signs['positive'] >= 3


def test_no_resonance_changes_sign_along_the_cut(summary):
    """The one observable N1 names."""
    assert summary['verdict']['any_track_changes_sign'] is False
    assert summary['verdict']['signs']['mixed'] == 0


def test_the_best_candidate_is_a_frequency_crossing(summary):
    """Two positive segments at one frequency, with a deep dip in between."""
    result = summary['bridge']
    assert result['min'] < -0.5
    assert result['is_zero_crossing'] is False


def test_the_negatives_fade_at_the_ends(summary):
    profile = summary['profile']
    assert profile['negative_ends'] < 0.5 * profile['negative_middle']


def test_but_no_positive_resonance_is_identified_there(summary):
    """A zero contour would kill the negatives and leave the positives.

    Instead every identified positive track lives in the interior: none
    reaches the first 250 or the last 80 rows, where the negatives are gone.
    """
    positives = [t for t in summary['found'] if t['sign'] == 'positive']
    assert positives
    assert min(t['y'][0] for t in positives) > 250
    assert max(t['y'][-1] for t in positives) < 1450


def test_the_positive_weight_at_the_ends_is_near_the_noise_floor(summary):
    check = summary['verdict']
    assert check['positive_weight_noise_floor'] is not None
    assert summary['profile']['positive_ends'] < 3.0 * check[
        'positive_weight_noise_floor']


# --- the verdict ----------------------------------------------------------

def test_n1_is_not_confirmed(summary):
    assert summary['verdict']['n1_confirmed'] is False


def test_the_summed_weight_would_have_given_a_false_positive(summary):
    """Why the verdict is not scored on total positive weight."""
    profile = summary['profile']
    assert profile['positive_ends'] > 1.5 * profile['positive_middle']
    assert summary['verdict']['n1_confirmed'] is False


def test_the_contrast_map_is_reproducible(panel):
    D = contrast_map(panel[0], DEFAULT_WINDOW)
    assert D.shape == PANEL_B_SHAPE
    assert D.min() < -0.6
    assert D.max() > 0.1
    assert abs(float(np.median(D))) < 0.02
