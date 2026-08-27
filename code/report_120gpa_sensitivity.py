"""Print the frozen 120 GPa optical-limit sensitivity table as CSV."""
import argparse

import numpy as np

from ho_odmr_sensitivity import HoIntegratedODMRModel


LINES = (405.0, 440.65, 445.0, 457.0, 473.0, 488.0, 505.0, 532.0)


def rows(pressure=120.0):
    model = HoIntegratedODMRModel()
    optimum = model.optimum(pressure)
    optimum_rate = model.absorbed_photon_proxy(optimum, pressure)
    for wavelength in LINES:
        rate = model.absorbed_photon_proxy(wavelength, pressure)
        yield {
            'pressure_GPa': pressure,
            'wavelength_nm': wavelength,
            'relative_absorbed_photon_rate': float(rate / optimum_rate),
            'optical_limit_sensitivity_penalty': float(
                np.sqrt(optimum_rate / rate)),
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pressure', type=float, default=120.0)
    args = parser.parse_args()
    print('pressure_GPa,wavelength_nm,relative_absorbed_photon_rate,'
          'optical_limit_sensitivity_penalty')
    for row in rows(args.pressure):
        print(f'{row["pressure_GPa"]:.1f},{row["wavelength_nm"]:.2f},'
              f'{row["relative_absorbed_photon_rate"]:.8f},'
              f'{row["optical_limit_sensitivity_penalty"]:.8f}')


if __name__ == '__main__':
    main()
