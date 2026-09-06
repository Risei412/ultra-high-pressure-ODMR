"""Literature-anchored low-power model for a [111]-aligned NV in a DAC.

Validated and conditional pieces are kept separate.

Validated optical input
-----------------------
``HoPublishedSpectrumModel`` reconstructs the hydrostatic absorption spectra of
Ho et al., "Optical Stability and Photophysics of NV Centers in Diamond up to
120 GPa" (2026), arXiv:2606.02399. At fixed incident optical power and low
power, the wavelength objective is ``lambda * sigma_abs``.

Literature support for [111] stress
-----------------------------------
Huang et al., "Elucidating the Inter-system Crossing of the Nitrogen-Vacancy
Center up to Megabar Pressures" (2026), arXiv:2511.20750, treat general stress
and find that uniaxial [111] loading preserves the NV symmetry and gives the
largest contrast. Wang et al., Nature Communications 15, 8843 (2024),
doi:10.1038/s41467-024-52272-y, observe high ODMR contrast in a (111)-cut DAC.
Those results constrain ISC/contrast, but do not publish a 120-GPa,
orientation-resolved absorption spectrum from which a new excitation optimum
can be reconstructed. No stress-induced optical shift is invented here: the
calculation reports the reproducible Ho optical baseline and labels the missing
[111] optical coupling explicitly.
"""

from __future__ import annotations

import numpy as np

from ho_spectrum_model import (HBARC, REPORTED_WAVELENGTH_RESOLUTION_NM,
                               HoPublishedSpectrumModel)
from theory_a1_generalization import DATA_WINDOW, Kernel


N111 = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
DEFAULT_PRESSURE_GPA = 120.0
DEFAULT_ALPHA = 0.56
DEFAULT_POWER_RATIOS = (
    1.00, 1.01, 1.05, 1.10, 1.20, 1.40, 1.45,
    1.50, 1.52, 1.60, 2.00, 4.00, 5.00,
)

__all__ = [
    "DEFAULT_ALPHA",
    "DEFAULT_POWER_RATIOS",
    "DEFAULT_PRESSURE_GPA",
    "N111",
    "NV111AxialModel",
    "REPORTED_WAVELENGTH_RESOLUTION_NM",
    "axisymmetric_111_stress",
    "nv111_basis",
    "stress_in_nv111_frame",
    "stress_invariants",
]


def nv111_basis() -> np.ndarray:
    """Return columns (x_NV, y_NV, z_NV) in cubic crystal coordinates."""
    z_nv = N111
    x_nv = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    y_nv = np.cross(z_nv, x_nv)
    return np.column_stack((x_nv, y_nv, z_nv))


def axisymmetric_111_stress(P: float, alpha: float) -> np.ndarray:
    """Return stress in cubic [100]/[010]/[001] coordinates, in GPa.

    The principal stresses in the aligned NV/DAC frame are
    ``(alpha*P, alpha*P, P)``. Thus alpha=1 is hydrostatic and alpha=0 is
    purely uniaxial [111], matching the convention used by Huang et al.
    Compression is positive.
    """
    if P < 0.0:
        raise ValueError("P must be non-negative")
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must be between 0 (uniaxial) and 1 (hydrostatic)")
    return (alpha * P * np.eye(3)
            + (1.0 - alpha) * P * np.outer(N111, N111))


def stress_in_nv111_frame(P: float, alpha: float) -> np.ndarray:
    """Rotate the cubic-coordinate stress tensor into the aligned NV frame."""
    basis = nv111_basis()
    return basis.T @ axisymmetric_111_stress(P, alpha) @ basis


def stress_invariants(P: float, alpha: float) -> dict[str, float]:
    """Return axial, mean, and transverse symmetry components of stress."""
    sigma = stress_in_nv111_frame(P, alpha)
    mean = float(np.trace(sigma) / 3.0)
    return {
        "sigma_x": float(sigma[0, 0]),
        "sigma_y": float(sigma[1, 1]),
        "sigma_z": float(sigma[2, 2]),
        "sigma_mean": mean,
        "sigma_dev_z_minus_xy": float(
            sigma[2, 2] - 0.5 * (sigma[0, 0] + sigma[1, 1])),
        "transverse_E1": float(0.5 * (sigma[0, 0] - sigma[1, 1])),
        "transverse_E2": float(sigma[0, 1]),
    }


class NV111AxialModel:
    """Low-power optical calculation with an explicit [111] stress tensor.

    Wavelength-independent contrast changes sensitivity magnitude but cannot
    change the low-power optimum, so contrast is not fitted here.
    """

    def __init__(
        self,
        *,
        alpha: float = DEFAULT_ALPHA,
    ):
        if not 0.0 <= alpha <= 1.0:
            raise ValueError("alpha must be between 0 and 1")
        self.alpha = float(alpha)
        self.optical = HoPublishedSpectrumModel()
        self._kernel_cache = {}

    @property
    def is_orientation_resolved_prediction(self) -> bool:
        """False: the literature lacks a calibrated [111] absorption kernel."""
        return False

    def stress_tensor(self, P: float, frame: str = "crystal") -> np.ndarray:
        if frame == "crystal":
            return axisymmetric_111_stress(P, self.alpha)
        if frame in {"nv", "dac"}:
            return stress_in_nv111_frame(P, self.alpha)
        raise ValueError("frame must be 'crystal', 'nv', or 'dac'")

    def stress_invariants(self, P: float) -> dict[str, float]:
        return stress_invariants(P, self.alpha)

    def sigma_abs(self, wavelength_nm, hydrostatic_pressure_gpa):
        """Reproduced hydrostatic Ho kernel, not an anisotropic-stress model.

        ``hydrostatic_pressure_gpa`` is the pressure coordinate of Ho's
        published spectra.  It must not be identified with one principal
        stress of the tensor returned by :meth:`stress_tensor`: the literature
        does not provide the optical coupling needed to make that conversion.
        """
        wavelength_nm = np.asarray(wavelength_nm, float)
        return self.optical.sigma_abs(
            HBARC / wavelength_nm, hydrostatic_pressure_gpa)

    def absorbed_photon_proxy(self, wavelength_nm, hydrostatic_pressure_gpa):
        """Absorbed rate at fixed incident optical power, up to a constant."""
        wavelength_nm = np.asarray(wavelength_nm, float)
        return wavelength_nm * self.sigma_abs(
            wavelength_nm, hydrostatic_pressure_gpa)

    def low_power_sensitivity(self, wavelength_nm, P, contrast_scale=1.0):
        """Shot-noise scale for wavelength-independent non-optical factors."""
        if contrast_scale <= 0.0:
            raise ValueError("contrast_scale must be positive")
        rate = np.clip(self.absorbed_photon_proxy(wavelength_nm, P), 1e-300, None)
        return 1.0 / (contrast_scale * np.sqrt(rate))

    def lambda_opt(self, hydrostatic_pressure_gpa, lam_min=400.0,
                   lam_max=600.0, step=0.05) -> float:
        """Numerical optimum of the hydrostatic baseline.

        The fine scan locates the maximum of the interpolant; it does not add
        wavelength information to the digitised source curve.  Report this
        value to no finer than ``REPORTED_WAVELENGTH_RESOLUTION_NM``.
        """
        wavelengths = np.arange(lam_min, lam_max + step / 2.0, step)
        objective = self.absorbed_photon_proxy(
            wavelengths, hydrostatic_pressure_gpa)
        return float(wavelengths[int(np.argmax(objective))])

    def reported_lambda_opt(self, hydrostatic_pressure_gpa) -> float:
        """Hydrostatic baseline rounded to the supported reporting precision."""
        return self.optical.reported_lambda_opt(
            hydrostatic_pressure_gpa, lam_min=400.0, lam_max=600.0)

    def finite_power_optima(self, P, power_ratio, window=DATA_WINDOW):
        """Return the optimal wavelength set at normalized power ``I/Ic``.

        This applies the existing mediated-response level-set result. ``Ic`` is
        the power at which the maximum-absorption wavelength reaches the
        response-optimal pumping rate. Its absolute value is not calibrated for
        a [111] DAC. At and below ``Ic`` the optimum is one point; above it the
        optimum is every solution of ``A(lambda)/Amax = Ic/I``.
        """
        if power_ratio <= 0.0:
            raise ValueError("power_ratio must be positive")
        key = (float(P), tuple(float(value) for value in window))
        if key not in self._kernel_cache:
            self._kernel_cache[key] = Kernel(pressure=P, window=window)
        kernel = self._kernel_cache[key]
        if power_ratio <= 1.0:
            return (kernel.lam_abs,)
        return tuple(kernel.level_set(1.0 / power_ratio, window, step=0.01))
