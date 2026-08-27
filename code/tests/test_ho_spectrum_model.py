"""Cross-figure reproduction tests for the published Ho spectra."""
import numpy as np
import pytest

from ho_spectrum_model import HoPublishedSpectrumModel, HBARC
from repro_yield import (LINES, compare_experiment, compare_ho_theory, load,
                         reproduction_gate)


@pytest.fixture(scope='module')
def model():
    return HoPublishedSpectrumModel()


@pytest.fixture(scope='module')
def data():
    return load()


def test_reference_pressure_grid_is_complete(model):
    assert np.array_equal(model.pressures, np.arange(0.0, 121.0, 20.0))


def test_interpolation_round_trips_extracted_path_points(model):
    for pressure in model.pressures:
        energy, reference = model.spectra[pressure]
        calculated = model.sigma_abs(energy, pressure)
        assert calculated == pytest.approx(reference, abs=1e-10)


def test_scalar_and_vector_pressure_calls_agree(model):
    energy = HBARC / 457.0
    pressure = np.array([0.0, 37.0, 88.0, 120.0])
    vector = model.sigma_abs(energy, pressure)
    scalar = np.array([model.sigma_abs(energy, p) for p in pressure])
    assert vector == pytest.approx(scalar, rel=1e-12)


def test_fig1e_reconstructs_fig5b_theory(model, data):
    """Independent cross-figure validation of extraction and interpolation."""
    comparison = compare_ho_theory(model, data)
    assert comparison['532']['fractional_rms'] < 0.02
    assert comparison['457']['fractional_rms'] < 0.02
    assert comparison['532']['peak_error_GPa'] <= 1.0
    assert comparison['457']['peak_error_GPa'] <= 1.0
    assert reproduction_gate(model, data)['ho_theory']


def test_published_model_matches_experimental_dome_shape(model, data):
    comparison = compare_experiment(model, data, collection=False)
    assert comparison['532']['fractional_rms'] < 0.16
    assert comparison['457']['fractional_rms'] < 0.16
    assert comparison['532']['peak_error_GPa'] < 10.0
    assert comparison['457']['peak_error_GPa'] < 10.0
    assert reproduction_gate(model, data)['experiment_shape']


def test_published_curve_optimum_at_120_gpa_is_near_441_nm(model):
    assert model.lambda_opt(120.0) == pytest.approx(440.6, abs=0.5)
    # The published spectrum supports 457 nm as a practical fixed line, while
    # distinguishing it from the reconstructed optimum.
    optimum_lam = model.lambda_opt(120.0)
    optimum_rate = optimum_lam * model.sigma_abs(
        HBARC / optimum_lam, 120.0)
    blue_rate = 457.0 * model.sigma_abs(HBARC / 457.0, 120.0)
    assert blue_rate / optimum_rate > 0.90


def test_reference_model_never_applies_collection_factor(model):
    assert np.all(model.eta_col(np.array([0.0, 60.0, 120.0])) == 1.0)


def test_lines_used_for_cross_figure_validation_are_unchanged():
    assert LINES == (532.0, 457.0)

