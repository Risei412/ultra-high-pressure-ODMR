"""Regression tests for the internal sanity checks on the Ho kernel."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kernel_sanity_checks import (  # noqa: E402
    ZPL_TRUSTED_MAX_GPA, jahn_teller_coupling, normalisation, report,
    zpl_area, zpl_width,
)


@pytest.fixture(scope='module')
def checks():
    return report()


# ---- K1 -------------------------------------------------------------------

def test_spectra_are_unit_normalised_where_they_are_complete(checks):
    for pressure in (0.0, 20.0, 40.0, 60.0):
        assert checks['normalisation'][pressure]['area'] == pytest.approx(
            1.0, abs=0.005)


def test_high_pressure_area_is_lost_to_the_blue_truncation(checks):
    norm = checks['normalisation']
    assert norm[120.0]['area'] < 0.93
    # the deficit tracks the weight still present at the top of the window
    fractions = [norm[p]['top_fraction'] for p in sorted(norm)]
    assert fractions == sorted(fractions)
    assert norm[0.0]['top_fraction'] < 0.01
    assert norm[120.0]['top_fraction'] > 0.7


def test_the_blue_edge_is_a_figure_limit_not_a_band_edge(checks):
    """At 120 GPa the band is still most of its peak where the figure stops."""
    assert checks['normalisation'][120.0]['top_fraction'] > 0.5


# ---- K2 -------------------------------------------------------------------

def test_zpl_width_does_not_move_over_120_GPa(checks):
    values = np.array(list(checks['zpl_width'].values()))
    assert (values.max() - values.min()) / values.mean() < 0.05


def test_zpl_width_is_at_the_plot_resolution(checks):
    values = np.array(list(checks['zpl_width'].values()))
    assert values.mean() == pytest.approx(0.00256, abs=1e-4)


def test_zpl_broadening_cannot_be_read_from_the_kernel(checks):
    """Theorem X's second driver is unavailable from Fig. 1(e), like E3's."""
    width = checks['zpl_width']
    assert width[120.0] / width[0.0] == pytest.approx(1.0, abs=0.05)


# ---- K3 -------------------------------------------------------------------

def test_kernel_zpl_area_validates_against_published_dwf_at_low_pressure(
        checks):
    for pressure, row in checks['zpl_area'].items():
        if pressure <= ZPL_TRUSTED_MAX_GPA:
            assert row['ratio'] == pytest.approx(0.90, abs=0.03)


def test_kernel_zpl_sits_uniformly_low_against_the_published_dwf(checks):
    """K3 against the corrected panel (c): a constant offset, not a trend.

    Read from the old raster extraction of panel (c), this ratio climbed from
    0.96 to 1.43 and was reported as the reconstruction over-weighting the
    zero-phonon line at high pressure.  Against the vector extraction it is
    0.88-0.91 at every pressure with no monotone trend, which is a calibration
    constant rather than a pressure-dependent artefact.  Everything in the
    manuscript that rests on a x1.43 over-weight rests on the old number.
    """
    ratios = [checks['zpl_area'][p]['ratio'] for p in sorted(checks['zpl_area'])]
    assert max(ratios) - min(ratios) < 0.05
    assert all(0.85 < ratio < 0.95 for ratio in ratios)
    assert checks['zpl_area'][120.0]['ratio'] == pytest.approx(0.89, abs=0.03)


def test_kernel_zpl_and_published_dwf_fall_together(checks):
    """K3 against the corrected panel (c): the two collapses agree to 1.3 %.

    The kernel's own zero-phonon area falls by x6.08 over 0-120 GPa and the
    published Debye-Waller factor by x6.01.  Read from the old raster
    extraction the published fall was x9.07 against the kernel's x6.3, a 44 %
    discrepancy, and that is what the manuscript's "the reconstruction
    over-weights the ZPL" passages were built on.
    """
    area = checks['zpl_area']
    kernel = area[0.0]['kernel_dwf'] / area[120.0]['kernel_dwf']
    published = area[0.0]['published_dwf'] / area[120.0]['published_dwf']
    assert published == pytest.approx(6.01, abs=0.1)
    assert kernel == pytest.approx(6.08, abs=0.1)
    assert abs(kernel / published - 1.0) < 0.03


# ---- K4 -------------------------------------------------------------------

def test_dwf_is_not_exp_minus_s_abs(checks):
    for row in checks['jahn_teller'].values():
        assert row['S_JT'] > 0.0


def test_jahn_teller_remainder_is_smooth_and_monotone(checks):
    jt = checks['jahn_teller']
    values = [jt[p]['S_JT'] for p in sorted(jt)]
    assert values == sorted(values)
    # smooth: no second difference larger than the first differences
    first = np.diff(values)
    assert np.all(np.abs(np.diff(first)) < first.min())


def test_total_coupling_endpoints(checks):
    jt = checks['jahn_teller']
    assert jt[0.0]['S_total'] == pytest.approx(3.83, abs=0.01)
    assert jt[120.0]['S_total'] == pytest.approx(5.62, abs=0.01)
