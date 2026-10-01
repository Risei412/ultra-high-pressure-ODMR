"""Tests for the literature-anchored [111] low-power model."""

import numpy as np
import pytest

from nv_model_111 import (
    NV111AxialModel,
    axisymmetric_111_stress,
    stress_in_nv111_frame,
)


def test_axisymmetric_tensor_has_requested_principal_stresses():
    values = np.linalg.eigvalsh(axisymmetric_111_stress(120.0, 0.56))
    assert values == pytest.approx([67.2, 67.2, 120.0])


def test_aligned_nv_frame_has_no_transverse_or_shear_stress():
    sigma = stress_in_nv111_frame(120.0, 0.56)
    assert sigma == pytest.approx(np.diag([67.2, 67.2, 120.0]), abs=1e-12)


def test_default_reproduces_published_kernel_low_power_optimum():
    model = NV111AxialModel(alpha=0.56)
    assert model.lambda_opt(120.0) == pytest.approx(441.0, abs=0.5)
    assert model.reported_lambda_opt(120.0) == 441.0
    assert model.is_orientation_resolved_prediction is False


def test_alpha_does_not_masquerade_as_an_optical_stress_prediction():
    """The only available optical kernel is hydrostatic and alpha-independent."""
    values = [NV111AxialModel(alpha=alpha).reported_lambda_opt(120.0)
              for alpha in (0.0, 0.56, 0.95, 1.0)]
    assert values == [441.0] * 4


def test_wavelength_independent_contrast_does_not_move_optimum():
    model = NV111AxialModel(alpha=0.56)
    wavelengths = np.arange(400.0, 600.0001, 0.05)
    optima = []
    for contrast in (0.03, 0.10, 0.30):
        eta = model.low_power_sensitivity(
            wavelengths, 120.0, contrast_scale=contrast)
        optima.append(float(wavelengths[np.argmin(eta)]))
    assert optima == pytest.approx([441.0, 441.0, 441.0], abs=0.5)


def test_invalid_alpha_is_rejected():
    with pytest.raises(ValueError):
        NV111AxialModel(alpha=1.1)


def test_finite_power_optimum_splits_above_ic():
    model = NV111AxialModel(alpha=0.56)
    assert model.finite_power_optima(120.0, 1.0) == pytest.approx((440.64,), abs=0.01)
    assert model.finite_power_optima(120.0, 1.01) == pytest.approx(
        (437.49, 446.36), abs=0.02)
    assert model.finite_power_optima(120.0, 1.10) == pytest.approx(
        (426.61, 457.69), abs=0.02)


def test_multimodal_kernel_adds_zpl_pair_at_higher_power():
    model = NV111AxialModel(alpha=0.56)
    assert model.finite_power_optima(120.0, 1.50) == pytest.approx(
        (409.83, 474.06, 514.43, 514.48), abs=0.03)


def test_invalid_power_ratio_is_rejected():
    with pytest.raises(ValueError):
        NV111AxialModel().finite_power_optima(120.0, 0.0)
