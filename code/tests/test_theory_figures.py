"""Regression tests for the three theory figures (Theorem M, G and X).

A figure is a claim.  These tests pin the numbers the three publication
figures draw and print, so that a change in the kernel, in A2's ladder, or in
the E3-corrected branch ratio cannot leave a caption asserting something the
figure no longer shows.

They also pin the two negative claims the figures are built to make:
the N = 6 plateau is too narrow to test, and there is no single exchange
pressure P*.
"""
import os
import sys

import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import fig_g_gauge_degeneracy as fig_g  # noqa: E402
import fig_style  # noqa: E402
import fig_m_level_set_ladder as fig_m  # noqa: E402
import fig_x_branch_exchange as fig_x  # noqa: E402
from figure_validation import branch_ratio_density  # noqa: E402
from theory_a2_multiplicity import (gauge_degeneracy,  # noqa: E402
                                    match_transitions, plateau_widths)


# The `kernel` fixture is session scoped in tests/conftest.py so that this
# file and test_theory_a2_multiplicity.py share one instance; see the note
# there. Do not redefine it here.


# ------------------------------------------------- figure M, panel (a)

def test_the_four_interior_maxima_the_figure_marks(kernel):
    """(a): the kernel is multimodal, and these are the four peaks drawn."""
    maxima = kernel.local_maxima()
    assert len(maxima) == 4
    expected = [(440.64, 1.0000), (475.55, 0.6603),
                (500.19, 0.2772), (514.46, 0.6938)]
    for (lam, value), (lam_ref, value_ref) in zip(maxima, expected):
        assert lam == pytest.approx(lam_ref, abs=0.02)
        assert value == pytest.approx(value_ref, abs=5e-4)


def test_the_drawn_level_line_really_has_four_members(kernel):
    """(a): the point of the panel is that a level meets the curve >2 times."""
    level, members = fig_m.level_line(kernel)
    assert level == pytest.approx(0.6603, abs=2e-3)
    assert len(members) == 4
    assert members == sorted(members)
    # One on the blue flank, one on the red flank of the main band, and the
    # near-degenerate pair straddling the zero-phonon line.
    assert members[0] == pytest.approx(409.5, abs=0.5)
    assert members[1] == pytest.approx(474.4, abs=0.5)
    assert members[3] - members[2] < 0.2
    assert all(m > 505.0 for m in members[2:])


def test_the_level_line_sits_inside_the_four_fold_plateau(kernel):
    """(a) and (b) are the same statement: 1/a* lands on the N = 4 rung."""
    level, members = fig_m.level_line(kernel)
    plateau = [row for row in plateau_widths(kernel)
               if row['from_power_ratio'] < 1.0 / level < row['to_power_ratio']]
    assert len(plateau) == 1
    assert plateau[0]['multiplicity'] == len(members) == 4


# ------------------------------------------------- figure M, panel (b)

def test_the_six_rung_positions(kernel):
    """(b): the six vertical lines, in power ratio and generating wavelength."""
    rows = match_transitions(kernel)['rows']
    expected = [(1.4414, 514.46, 'max', 2, 4),
                (1.5145, 475.55, 'max', 4, 6),
                (1.5191, 474.47, 'min', 6, 4),
                (1.8646, 402.00, 'edge-blue', 4, 3),
                (3.6078, 500.19, 'max', 3, 5),
                (4.2989, 497.86, 'min', 5, 3)]
    assert len(rows) == len(expected)
    for row, (ratio, lam, kind, before, after) in zip(rows, expected):
        assert row['observed_power_ratio'] == pytest.approx(ratio, abs=5e-4)
        assert row['predicted_power_ratio'] == pytest.approx(ratio, abs=5e-4)
        assert row['predicted_from_nm'] == pytest.approx(lam, abs=0.02)
        assert row['kind'] == kind
        assert (row['before'], row['after']) == (before, after)


def test_the_annotated_six_fold_plateau_is_not_resolvable(kernel):
    """(b): the annotation says x1.003, and it must stay true."""
    rows = plateau_widths(kernel)
    narrow = min(rows, key=lambda row: row['width_factor'])
    assert narrow['multiplicity'] == 6
    assert narrow['width_factor'] == pytest.approx(1.003, abs=0.002)
    # Every other resolvable plateau is at least an order of magnitude wider
    # in log-power, which is why only this one carries a warning.
    others = [row['width_factor'] - 1.0 for row in rows
              if row['multiplicity'] != 6 and np.isfinite(row['width_factor'])]
    assert min(others) > 10.0 * (narrow['width_factor'] - 1.0)


# ------------------------------------------------------------- figure G

def test_the_gauge_family_is_exactly_degenerate_in_phi():
    """(a): the five curves coincide to floating-point noise, not to the eye."""
    rows = gauge_degeneracy()
    assert len(rows) == 5
    for row in rows:
        assert row['exponent'] == pytest.approx(2.0)
        assert row['max_relative_phi_difference'] < 1e-12
        assert row['gamma_star'] == pytest.approx(1.0, rel=1e-9)
    residual = max(row['max_relative_phi_difference'] for row in rows)
    assert residual < 1e-12
    # It is a rounding residual, not a real difference: a few ulp.
    assert residual < 1e-14


def test_the_plotted_phi_curves_are_one_curve():
    """(a): tested on the drawn arrays, not only on the reported summary."""
    curves = fig_g.phi_curves()
    assert len(curves) == 5
    reference = curves['contrast collapse alone']
    assert float(np.max(reference)) == pytest.approx(1.0, rel=1e-12)
    for values in curves.values():
        assert np.max(np.abs(values - reference) / reference) < 1e-12


def test_the_measurable_factors_differ_by_up_to_a_factor_two():
    """(b): the identifiability point -- what eta hides is not small."""
    factors = fig_g.measurable_factors()
    assert factors['contrast collapse alone'] == pytest.approx(
        {'R': 1.0, 'C': 0.5, 'dnu': 1.0}, rel=1e-9)
    assert factors['saturation + broadening'] == pytest.approx(
        {'R': 0.5, 'C': 1.0, 'dnu': np.sqrt(2.0)}, rel=1e-9)
    assert factors['linewidth^1 (strong)'] == pytest.approx(
        {'R': 1.0, 'C': 1.0, 'dnu': 2.0}, rel=1e-9)
    for key in ('R', 'C', 'dnu'):
        values = [row[key] for row in factors.values()]
        assert max(values) / min(values) == pytest.approx(2.0, rel=1e-9)


# ------------------------------------------------------------- figure X

def test_r_of_p_is_monotone_with_the_growth_the_caption_quotes():
    """(a): E3-corrected branch ratio density, from panels (b) and (c)."""
    data = branch_ratio_density()
    expected = [97.1, 121.9, 162.1, 208.0, 267.1, 329.1, 418.4]
    assert list(data['pressure']) == [0, 20, 40, 60, 80, 100, 120]
    for value, reference in zip(data['r'], expected):
        assert value == pytest.approx(reference, rel=2e-3)
    assert data['monotone'] is True
    assert np.all(np.diff(data['r']) > 0)
    assert 4.2 < data['growth_factor'] < 4.4


def test_the_critical_bandwidth_axis_endpoints():
    """(a): the right-hand axis runs 10.30 meV down to 2.39 meV."""
    widths = branch_ratio_density()['critical_bandwidth_meV']
    assert widths[0] == pytest.approx(10.30, rel=2e-3)
    assert widths[-1] == pytest.approx(2.39, rel=3e-3)
    assert np.all(np.diff(widths) < 0)


def test_the_six_marked_exchange_pressures():
    """(b): the (W, P*) pairs the panel marks and labels."""
    expected = {2.5: 116.8, 3.0: 101.2, 4.0: 74.5,
                5.0: 56.8, 7.0: 31.0, 9.0: 13.0}
    assert set(fig_x.MARKED_BANDWIDTHS_MEV) == set(expected)
    for width_mev in fig_x.MARKED_BANDWIDTHS_MEV:
        star = fig_x.exchange_pressure_at_bandwidth(width_mev * 1e-3)
        assert star == pytest.approx(expected[width_mev], abs=0.3)
    # A broader laser favours the sideband, so P* falls with W.
    stars = [fig_x.exchange_pressure_at_bandwidth(w * 1e-3)
             for w in fig_x.MARKED_BANDWIDTHS_MEV]
    assert stars == sorted(stars, reverse=True)


def test_the_shaded_region_is_where_no_crossing_exists():
    """(b): outside ~2.39-10.30 meV the crossing leaves the published range."""
    lo, hi = fig_x.accessible_bandwidths()
    assert lo == pytest.approx(2.39, rel=3e-3)
    assert hi == pytest.approx(10.30, rel=3e-3)
    widths, stars = fig_x.exchange_curve()
    inside = np.isfinite(stars)
    assert np.all(stars[inside] >= fig_x.P_MIN - 1e-6)
    assert np.all(stars[inside] <= fig_x.P_MAX + 1e-6)
    assert np.all(~inside[widths < lo - 0.02])
    assert np.all(~inside[widths > hi + 0.02])


def test_the_figure_never_touches_the_clipped_zpl_peak_heights():
    """(a) and (b) both derive the ZPL weight from the published DWF.

    Checked on the figure module's own source: it may reach the sideband
    through ``figure_validation``, which is audited, but it must never read a
    ZPL height out of panel (e) itself.
    """
    with open(fig_x.__file__, encoding='utf-8') as stream:
        source = stream.read()
    assert 'CLIPPED_ZPL_TOPS' not in source
    assert 'identify_branches' not in source
    assert 'from theory_a3_branch_exchange' not in source
    assert 'import theory_a3_branch_exchange' not in source
    assert 'DWF' in source


# ------------------------------------------------------- the files render

@pytest.mark.parametrize('module', [fig_m, fig_g, fig_x])
def test_each_figure_renders_both_formats(module, tmp_path, monkeypatch,
                                          kernel):
    monkeypatch.setattr(module, 'OUT_PNG', tmp_path / 'figure.png')
    monkeypatch.setattr(module, 'OUT_PDF', tmp_path / 'figure.pdf')
    if hasattr(module, 'Kernel'):
        # The ladder scan costs about 95 s; hand the figure module the
        # session fixture's kernel so its cached result is shared with every
        # other test in the suite (tests/conftest.py).
        monkeypatch.setattr(module, 'Kernel', lambda: kernel)
    module.main()
    assert (tmp_path / 'figure.png').stat().st_size > 20_000
    assert (tmp_path / 'figure.pdf').stat().st_size > 5_000
    # Each of these is a full-width `figure*` included at 	extwidth, so the
    # file has to BE 	extwidth: any other width is rescaled on the page, and
    # rescaling a figure rescales its lettering.  fig_style.save() enforces
    # this on write; this assertion is what fails if that is ever removed.
    width_pt = fig_style.verify_width(tmp_path / 'figure.pdf',
                                      fig_style.FULL_WIDTH_IN)
    assert width_pt == pytest.approx(fig_style.FULL_WIDTH_IN * 72.0, abs=0.75)
