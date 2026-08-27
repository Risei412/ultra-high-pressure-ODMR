"""Published-curve surrogate for Ho et al. (2026) absorption spectra.

This is not a DFT or dynamical-Jahn--Teller implementation. It reconstructs the
seven absorption spectra published in Fig. 1(e), interpolates them in energy,
and uses the cubic pressure interpolation described for Fig. 5(b). Fig. 5(b)
therefore provides an independent cross-figure validation of the extraction.
"""
import os

import numpy as np
from scipy.interpolate import CubicSpline

HBARC = 1239.84  # eV nm
DEFAULT_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'data', 'ho_fig1e_absorption.csv')


class HoPublishedSpectrumModel:
    """Interpolate the absorption spectra published at 20 GPa intervals."""

    def __init__(self, data_path=DEFAULT_DATA):
        rows = []
        with open(data_path) as stream:
            for line in stream:
                if (not line.strip() or line.startswith('#')
                        or line.startswith('pressure_')):
                    continue
                pressure, energy, absorption = line.split(',')
                rows.append((float(pressure), float(energy), float(absorption)))

        grouped = {}
        for pressure, energy, absorption in rows:
            grouped.setdefault(pressure, []).append((energy, absorption))
        self.pressures = np.array(sorted(grouped), float)
        if not np.array_equal(self.pressures, np.arange(0.0, 121.0, 20.0)):
            raise ValueError('expected spectra at 0, 20, ..., 120 GPa')

        self.spectra = {}
        for pressure, values in grouped.items():
            values = np.asarray(values, float)
            order = np.argsort(values[:, 0])
            energy, absorption = values[order].T
            if np.any(np.diff(energy) <= 0.0):
                raise ValueError(f'non-monotone energy grid at {pressure:g} GPa')
            self.spectra[pressure] = (energy, absorption)

    def _at_all_reference_pressures(self, energy):
        """Return shape (7, *energy.shape) sampled on every reference curve."""
        energy = np.asarray(energy, float)
        flat = energy.ravel()
        values = np.vstack([
            np.interp(flat, *self.spectra[pressure], left=0.0, right=0.0)
            for pressure in self.pressures
        ])
        return values.reshape((len(self.pressures),) + energy.shape)

    def sigma_abs(self, energy, pressure):
        """Absorption in published arbitrary units, with NumPy broadcasting."""
        energy = np.asarray(energy, float)
        pressure = np.asarray(pressure, float)
        if np.any((pressure < 0.0) | (pressure > 120.0)):
            raise ValueError('published-curve interpolation is limited to 0--120 GPa')

        energy, pressure = np.broadcast_arrays(energy, pressure)
        if pressure.ndim == 0 or np.all(pressure == pressure.flat[0]):
            values = self._at_all_reference_pressures(energy)
            result = CubicSpline(self.pressures, values, axis=0)(
                float(pressure.flat[0]))
            return np.clip(result, 0.0, None)

        if energy.ndim == 0 or np.all(energy == energy.flat[0]):
            values = self._at_all_reference_pressures(float(energy.flat[0]))
            result = CubicSpline(self.pressures, values)(pressure)
            return np.clip(result, 0.0, None)

        # General pairwise broadcast case. It is uncommon and intentionally
        # explicit; the two vectorised branches above cover all project uses.
        out = np.empty(energy.size, float)
        for index, (e_value, p_value) in enumerate(
                zip(energy.ravel(), pressure.ravel())):
            values = self._at_all_reference_pressures(float(e_value))
            out[index] = CubicSpline(self.pressures, values)(p_value)
        return np.clip(out.reshape(energy.shape), 0.0, None)

    @staticmethod
    def eta_col(pressure):
        """No detector factor: the reconstructed quantity is absorption."""
        return np.ones_like(np.asarray(pressure, float))

    def lambda_opt(self, pressure, lam_min=400.0, lam_max=850.0,
                   step=0.05, fixed_optical_power=True):
        """Maximise absorbed photon rate over the published spectral window."""
        wavelengths = np.arange(lam_min, lam_max + step / 2.0, step)
        absorption = self.sigma_abs(HBARC / wavelengths, pressure)
        objective = wavelengths * absorption if fixed_optical_power else absorption
        return float(wavelengths[int(np.argmax(objective))])


if __name__ == '__main__':
    model = HoPublishedSpectrumModel()
    for pressure in (0, 20, 40, 60, 80, 100, 120):
        print(f'{pressure:3d} GPa  lambda_opt = '
              f'{model.lambda_opt(pressure):.1f} nm')

