"""The esa verification figure must keep quoting the audited numbers.

The caption in section 2.3 of the esa article states the reproduction and
sanity-check results as fixed numbers.  These tests pin the ones the caption
quotes, so a change in the extraction or the audit cannot silently leave the
caption saying something the figure no longer shows.
"""
import matplotlib

matplotlib.use('Agg')

import pytest  # noqa: E402

import esa_figv_kernel_sanity as figure  # noqa: E402
from ho_spectrum_model import HoPublishedSpectrumModel  # noqa: E402
from repro_yield import (EXPT_MAX_FRACTIONAL_RMS,  # noqa: E402
                         EXPT_MAX_PEAK_ERROR_GPA, HO_MAX_FRACTIONAL_RMS,
                         HO_MAX_PEAK_ERROR_GPA, compare_experiment,
                         compare_ho_theory, load)


@pytest.fixture(scope='module')
def audit():
    model = HoPublishedSpectrumModel()
    data = load()
    return (compare_ho_theory(model, data),
            compare_experiment(model, data, collection=False))


def test_cross_figure_numbers_in_the_caption(audit):
    theory, _ = audit
    assert theory['532']['fractional_rms'] == pytest.approx(0.015, abs=0.002)
    assert theory['457']['fractional_rms'] == pytest.approx(0.004, abs=0.002)
    assert theory['532']['peak_error_GPa'] == 0.0
    assert theory['457']['peak_error_GPa'] == 0.0


def test_experimental_shape_numbers_in_the_caption(audit):
    _, experiment = audit
    assert experiment['532']['fractional_rms'] == pytest.approx(0.100,
                                                                abs=0.005)
    assert experiment['457']['fractional_rms'] == pytest.approx(0.102,
                                                                abs=0.005)
    assert experiment['532']['peak_error_GPa'] == pytest.approx(8.0, abs=1.0)
    assert experiment['457']['peak_error_GPa'] == pytest.approx(7.0, abs=1.0)


def test_every_bar_stays_under_its_ceiling(audit):
    theory, experiment = audit
    for lam in ('532', '457'):
        assert theory[lam]['fractional_rms'] <= HO_MAX_FRACTIONAL_RMS
        assert theory[lam]['peak_error_GPa'] <= HO_MAX_PEAK_ERROR_GPA
        assert experiment[lam]['fractional_rms'] <= EXPT_MAX_FRACTIONAL_RMS
        assert experiment[lam]['peak_error_GPa'] <= EXPT_MAX_PEAK_ERROR_GPA


def test_the_one_allowed_freedom_is_a_scale(audit):
    """The reconstruction may be rescaled; its shape carries no fit freedom."""
    model = HoPublishedSpectrumModel()
    data = load()
    pressure, reference = data['theory457_ho']
    scaled = figure._scaled_reconstruction(model, 457, pressure, reference)
    doubled = figure._scaled_reconstruction(model, 457, pressure,
                                            2.0 * reference)
    assert doubled == pytest.approx(2.0 * scaled, rel=1e-12)


def test_both_figures_render_at_the_width_they_are_printed_at(tmp_path,
                                                              monkeypatch):
    """Two full-width figures, and neither is rescaled by the page.

    Both are `figure*` floats included at `\\textwidth`, so a file that is not
    `\\textwidth` wide is scaled on the page and its lettering is scaled with
    it.  That is the regression that once put 4 pt labels in print, and the
    earlier version of this test could not have caught it: it measured
    `Figure.get_size_inches()`, the canvas *before* saving, whereas the width
    that matters is the one written into the file.  A `bbox_inches='tight'`
    save could therefore pass this test and still deliver 7.29 in.
    """
    import fig_style

    names = ('OUT_PNG_CROSS', 'OUT_PDF_CROSS',
             'OUT_PNG_INTERNAL', 'OUT_PDF_INTERNAL')
    written = {}
    for name in names:
        suffix = '.png' if 'PNG' in name else '.pdf'
        written[name] = tmp_path / (name.lower() + suffix)
        monkeypatch.setattr(figure, name, written[name])
    figure.main()

    for name in names:
        assert written[name].stat().st_size > 0
    expected_pt = fig_style.FULL_WIDTH_IN * 72.0
    for name in ('OUT_PDF_CROSS', 'OUT_PDF_INTERNAL'):
        width_pt = fig_style.verify_width(written[name],
                                          fig_style.FULL_WIDTH_IN)
        assert width_pt == pytest.approx(expected_pt, abs=0.75)
