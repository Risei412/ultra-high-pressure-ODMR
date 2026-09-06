"""Report the literature-anchored [111] low-power calculation.

Run with: ``python report_111_anisotropic.py --pressure 120 --alpha 0.56``
"""

from __future__ import annotations

import argparse
import numpy as np

from nv_model_111 import (
    DEFAULT_ALPHA,
    DEFAULT_POWER_RATIOS,
    DEFAULT_PRESSURE_GPA,
    NV111AxialModel,
)


def _fmt_matrix(matrix: np.ndarray) -> str:
    return "\n".join(
        "  [" + "  ".join(f"{value:8.3f}" for value in row) + "]"
        for row in matrix)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pressure", type=float, default=DEFAULT_PRESSURE_GPA,
                        help="pressure coordinate of the hydrostatic optical "
                             "kernel, in GPa")
    parser.add_argument("--axial-stress", type=float, default=None,
                        help="optional sigma_z for the illustrative [111] "
                             "stress tensor; no optical shift is inferred")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA,
                        help="sigma_x/sigma_z = sigma_y/sigma_z")
    parser.add_argument(
        "--power-ratios", type=float, nargs="+",
        default=DEFAULT_POWER_RATIOS,
        help="normalized incident powers I/Ic for the finite-power optimum set",
    )
    args = parser.parse_args()

    model = NV111AxialModel(alpha=args.alpha)
    optical_pressure = args.pressure
    axial_stress = (optical_pressure if args.axial_stress is None
                    else args.axial_stress)
    inv = model.stress_invariants(axial_stress)
    numerical_optimum_nm = model.lambda_opt(optical_pressure)
    reported_optimum_nm = model.reported_lambda_opt(optical_pressure)

    print("[111]-aligned, axisymmetric DAC: literature-anchored low-power report")
    print(f"hydrostatic optical-kernel pressure = {optical_pressure:.3f} GPa")
    print(f"illustrative sigma_z = {axial_stress:.3f} GPa; "
          f"alpha = {args.alpha:.3f}")
    print("stress tensor in cubic [100]/[010]/[001] coordinates (GPa):")
    print(_fmt_matrix(model.stress_tensor(axial_stress, "crystal")))
    print("stress tensor in the aligned NV/DAC frame (GPa):")
    print(_fmt_matrix(model.stress_tensor(axial_stress, "nv")))
    print(f"mean stress = {inv['sigma_mean']:.3f} GPa")
    print(f"axial deviatoric excess = {inv['sigma_dev_z_minus_xy']:.3f} GPa")
    print(f"transverse E stress = ({inv['transverse_E1']:.3e}, "
          f"{inv['transverse_E2']:.3e}) GPa")
    print()
    print(f"low-power optimum excitation = approximately "
          f"{reported_optimum_nm:.0f} nm")
    print(f"  interpolant maximum = {numerical_optimum_nm:.2f} nm; the extra "
          "digits are not source-data precision")
    print("status: reproduced hydrostatic Ho optical baseline; alpha and the")
    print("        illustrative stress tensor do not enter this wavelength.")
    print("        No conversion from a principal/mean stress to the optical")
    print("        pressure coordinate is justified by the available data.")
    print("        No unmeasured [111] stress-induced absorption shift has")
    print("        been added.")
    print("        [111] contrast can change sensitivity magnitude but not this")
    print("        optimum if it is wavelength independent.")
    print()
    print("finite-power optimum set (conditional on the mediated-response model):")
    print("  I/Ic     N    optimal wavelengths [nm]")
    for ratio in args.power_ratios:
        optima = model.finite_power_optima(optical_pressure, ratio)
        values = ", ".join(f"{value:.2f}" for value in optima)
        print(f"  {ratio:5.2f}    {len(optima):d}    {values}")
    print("Ic is not an absolute power prediction: it must be calibrated by a")
    print("single-wavelength power sweep in the actual [111] DAC.")


if __name__ == "__main__":
    main()
