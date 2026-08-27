"""Ho-published absorption kernel extended to an ODMR sensitivity model.

The optical kernel is reconstructed from Ho et al. Fig. 1(e).  Everything
after absorption is deliberately explicit and replaceable: charge-state yield,
ODMR contrast, linewidth, and saturation.  Consequently the module supports
two different statements without mixing them:

1. ``optical_limit``: a reproducible, low-power result with wavelength-
   independent non-optical factors;
2. ``sensitivity``: a conditional ODMR prediction once those factors are
   supplied or experimentally calibrated.

This is not a microscopic reimplementation of Ho et al.'s dynamical
Jahn--Teller calculation.
"""
from dataclasses import dataclass

import numpy as np

from ho_spectrum_model import HBARC, HoPublishedSpectrumModel


def _one(wavelength_nm, pressure_gpa):
    """Default wavelength-independent response factor."""
    wavelength_nm, pressure_gpa = np.broadcast_arrays(
        np.asarray(wavelength_nm, float), np.asarray(pressure_gpa, float))
    return np.ones_like(wavelength_nm)


@dataclass(frozen=True)
class ODMRResponse:
    """Non-optical response functions that must be measured or justified.

    All callables accept ``(wavelength_nm, pressure_gpa)`` and may return a
    scalar or broadcastable array.  ``charge_yield`` includes NV- population,
    radiative yield, and wavelength-dependent detection effects.  The
    saturation scale is dimensionless because incident power is reported
    relative to the absorbed-power scale at the optical optimum.
    """

    charge_yield: object = _one
    contrast: object = _one
    linewidth: object = _one
    saturation_scale: object = _one


class HoIntegratedODMRModel:
    """Convert Ho's absorption spectrum into conditional ODMR sensitivity."""

    def __init__(self, optical_model=None, response=None):
        self.optical = optical_model or HoPublishedSpectrumModel()
        self.response = response or ODMRResponse()

    def absorbed_photon_proxy(self, wavelength_nm, pressure_gpa):
        """Absorbed photons at fixed incident optical power, in relative units."""
        wavelength_nm = np.asarray(wavelength_nm, float)
        return wavelength_nm * self.optical.sigma_abs(
            HBARC / wavelength_nm, pressure_gpa)

    def _response_value(self, name, wavelength_nm, pressure_gpa):
        value = getattr(self.response, name)(wavelength_nm, pressure_gpa)
        value = np.asarray(value, float)
        if np.any(value <= 0.0):
            raise ValueError(f'{name} must be strictly positive')
        return value

    def detected_rate(self, wavelength_nm, pressure_gpa,
                      relative_power=1e-6):
        """Conditional detected rate including saturation and charge yield.

        ``relative_power`` is referenced so that the optical optimum at the
        selected pressure has an unsaturated excitation parameter equal to
        ``relative_power`` when ``saturation_scale == 1``.
        """
        wavelength_nm = np.asarray(wavelength_nm, float)
        if relative_power <= 0.0:
            raise ValueError('relative_power must be positive')
        optical_optimum = self.optical.lambda_opt(float(pressure_gpa))
        reference = self.absorbed_photon_proxy(optical_optimum, pressure_gpa)
        excitation = (relative_power
                      * self.absorbed_photon_proxy(wavelength_nm, pressure_gpa)
                      / reference)
        saturation = self._response_value(
            'saturation_scale', wavelength_nm, pressure_gpa)
        charge_yield = self._response_value(
            'charge_yield', wavelength_nm, pressure_gpa)
        return charge_yield * excitation / (1.0 + excitation / saturation)

    def sensitivity(self, wavelength_nm, pressure_gpa,
                    relative_power=1e-6):
        """Shot-noise CW-ODMR sensitivity in a common arbitrary scale."""
        rate = self.detected_rate(wavelength_nm, pressure_gpa, relative_power)
        contrast = self._response_value('contrast', wavelength_nm, pressure_gpa)
        linewidth = self._response_value(
            'linewidth', wavelength_nm, pressure_gpa)
        return linewidth / (contrast * np.sqrt(rate))

    def optimum(self, pressure_gpa, relative_power=1e-6,
                lam_min=400.0, lam_max=600.0, step=0.05):
        wavelengths = np.arange(lam_min, lam_max + step / 2.0, step)
        eta = self.sensitivity(wavelengths, pressure_gpa, relative_power)
        return float(wavelengths[int(np.argmin(eta))])

    def penalty(self, wavelength_nm, pressure_gpa, relative_power=1e-6,
                lam_min=400.0, lam_max=600.0, step=0.05):
        """Sensitivity divided by the best sensitivity on the scan grid."""
        wavelengths = np.arange(lam_min, lam_max + step / 2.0, step)
        eta = self.sensitivity(wavelengths, pressure_gpa, relative_power)
        target = self.sensitivity(
            float(wavelength_nm), pressure_gpa, relative_power)
        return float(target / np.min(eta))

    def pair_advantage(self, better_nm, reference_nm, pressure_gpa,
                       relative_power=1e-6):
        """Return eta(reference)/eta(better); values above one favour better."""
        eta_better = self.sensitivity(
            better_nm, pressure_gpa, relative_power)
        eta_reference = self.sensitivity(
            reference_nm, pressure_gpa, relative_power)
        return float(eta_reference / eta_better)

    def contrast_ratio_threshold(self, candidate_nm, reference_nm,
                                 pressure_gpa):
        """Low-power C_candidate/C_reference threshold for candidate to win.

        This diagnostic assumes equal linewidth and charge yield.  A measured
        contrast ratio above the returned number makes ``candidate_nm`` more
        sensitive than ``reference_nm`` in the optical-kernel limit.
        """
        candidate_rate = self.absorbed_photon_proxy(candidate_nm, pressure_gpa)
        reference_rate = self.absorbed_photon_proxy(reference_nm, pressure_gpa)
        return float(np.sqrt(reference_rate / candidate_rate))


def optical_limit_summary(pressure_gpa=120.0):
    """Return the frozen, assumption-labelled optical-limit numbers."""
    model = HoIntegratedODMRModel()
    optimum = model.optimum(pressure_gpa)
    threshold = model.contrast_ratio_threshold(457.0, 532.0, pressure_gpa)
    return {
        'pressure_GPa': float(pressure_gpa),
        'optimum_nm': optimum,
        'penalty_457': model.penalty(457.0, pressure_gpa),
        'penalty_532': model.penalty(532.0, pressure_gpa),
        # Exact unsaturated optical limit.  Do not evaluate this through a
        # merely small finite-power saturation approximation.
        'advantage_457_over_532': 1.0 / threshold,
        'minimum_C457_over_C532_for_457_to_win':
            threshold,
    }


if __name__ == '__main__':
    for key, value in optical_limit_summary().items():
        print(f'{key}: {value:.8g}')
