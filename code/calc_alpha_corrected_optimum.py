"""Conditional anisotropic-stress correction to the Ho 120 GPa optimum.

Ho et al. Fig. 1(e) supplies the hydrostatic (alpha = 1) absorption kernel.
The alpha = 0.95 value describes the micropillar validation geometry; it is
not the alpha assigned to that calculated absorption kernel.  This module
therefore anchors 440.65 nm at alpha = 1 and adds only the alpha-dependent
difference predicted by the phenomenological NV model.

This is a competing hypothesis to the static-kernel result (440.65 nm for
every alpha), not a replacement for a measured anisotropic absorption spectrum.

The deviatoric ratio is a RANGE, not a number
---------------------------------------------
Hilberer et al. report the ZPL slope in two geometries, -434 and
-769 meV/(cm^3 mol^-1), and attribute the difference to deviatoric stress.
Converting that ratio into k_q/k_h needs to know what was held equal between
the two geometries, and the published Letter does not settle it:

  equal compression (equal P_mean)   -> -0.734
  equal axial stress (equal sigma_zz) -> -0.391

The paper plots the ZPL against compressed diamond molar volume, which argues
for the first; but that volume is DERIVED from the calibrated diamond Raman
gauge at the anvil tip through a hydrostatic equation of state, and the Raman
gauge is itself referenced to the axial stress, which argues for the second.
See erratum E5 section 3.  Only the SIGN is settled.

``HILBERER_DEVIATORIC_RATIO`` keeps the larger-magnitude end as the default so
a single-number call stays on the conservative (more red-shifted) side; use
:func:`corrected_optimum_range_nm` to carry both.
"""

from contextlib import contextmanager
from functools import lru_cache

import numpy as np

import nv_model


HO_HYDROSTATIC_OPTIMUM_NM = 440.65
# The pressure at which HO_HYDROSTATIC_OPTIMUM_NM was obtained.  It is not a
# free parameter: the reconstructed kernel refuses to be evaluated anywhere
# else, so the additive anchor is only meaningful here.
HO_ANCHOR_PRESSURE_GPA = 120.0

# Both readings of the Hilberer 434/769 ratio; see the module docstring.
HILBERER_RATIO_EQUAL_COMPRESSION = -0.734
HILBERER_RATIO_EQUAL_AXIAL_STRESS = -0.39211
HILBERER_RATIO_RANGE = (HILBERER_RATIO_EQUAL_COMPRESSION,
                        HILBERER_RATIO_EQUAL_AXIAL_STRESS)

HILBERER_DEVIATORIC_RATIO = HILBERER_RATIO_EQUAL_COMPRESSION


@contextmanager
def _deviatoric_ratio(value):
    """Temporarily select the reviewed stress coefficient without global leak."""
    previous = nv_model.ZPL_DEV_K
    nv_model.ZPL_DEV_K = value
    try:
        yield
    finally:
        nv_model.ZPL_DEV_K = previous


def cancellation_alpha(deviatoric_ratio=HILBERER_DEVIATORIC_RATIO):
    """Return the alpha at which the modeled net ZPL shift vanishes.

    The shift factor is f(alpha) = (1 + 2 alpha)/3 + k (1 - alpha), the
    hydrostatic blueshift plus the deviatoric redshift.  Below the root of
    f the two-point Hilberer calibration predicts a 120 GPa ZPL BELOW the
    ambient one, which is unphysical: it is a linear extrapolation carried far
    outside the alpha = 0.56 / 0.95 interval it was fitted on.

    Roots: 0.286 for the equal-compression reading, 0.056 for equal axial
    stress.  Huang et al. recommend alpha -> 0 for contrast, which is on the
    wrong side of both -- this model cannot evaluate that proposal.
    """
    k = deviatoric_ratio
    return -(3.0 * k + 1.0) / (2.0 - 3.0 * k)


@lru_cache(maxsize=1024)
def deviatoric_shift_nm(alpha, sigma_zz_gpa=HO_ANCHOR_PRESSURE_GPA,
                        temperature_k=300.0,
                        deviatoric_ratio=HILBERER_DEVIATORIC_RATIO):
    """Return only the alpha-dependent shift, in nm, at any pressure.

        d_lambda(alpha, P) = lambda_opt(alpha, P) - lambda_opt(1, P)

    This is the half of the correction that is well defined away from the
    120 GPa anchor, because it is a difference taken at one pressure.  Adding
    it to a hydrostatic reference is the caller's job, and above 120 GPa the
    caller has to say which stand-in reference it used -- see
    ``extrapolation_bounds.py``.

    The shift GROWS with pressure: at alpha = 0.5 it is +60.8 nm at 120 GPa
    and +71.4 nm at 200 GPa (equal-compression reading).  That growth is real;
    what is not is adding it to a reference frozen at 120 GPa.
    """
    if not 0.0 < alpha <= 1.0:
        raise ValueError("alpha must satisfy 0 < alpha <= 1")
    floor = cancellation_alpha(deviatoric_ratio)
    if alpha <= floor:
        raise ValueError(
            "alpha = %.3f is at or below the %.3f where this model's net ZPL "
            "shift changes sign; the two-point calibration cannot be "
            "extrapolated there" % (alpha, floor))
    with _deviatoric_ratio(deviatoric_ratio):
        hydrostatic = nv_model.NVModel(T=temperature_k, alpha=1.0).lambda_opt(
            sigma_zz_gpa
        )
        anisotropic = nv_model.NVModel(T=temperature_k, alpha=alpha).lambda_opt(
            sigma_zz_gpa
        )
    return anisotropic - hydrostatic


def corrected_optimum_nm(alpha, sigma_zz_gpa=HO_ANCHOR_PRESSURE_GPA,
                         temperature_k=300.0,
                         deviatoric_ratio=HILBERER_DEVIATORIC_RATIO):
    """Return Ho hydrostatic optimum plus the modeled alpha-dependent shift.

    Raises for alpha at or below :func:`cancellation_alpha`, where the linear
    stress model has extrapolated past a sign change and stops meaning
    anything.  Hilberer calibrated it on alpha = 0.56 and 0.95; Huang et al.
    report NV-DAC experiments at alpha ~ 0.57, 0.73 and 1.  Trust it there.

    Raises for any pressure other than 120 GPa.  The additive constant
    ``HO_HYDROSTATIC_OPTIMUM_NM`` is Ho's hydrostatic optimum AT 120 GPa, and
    the reconstructed kernel that produced it cannot be evaluated anywhere
    else (``ho_spectrum_model.sigma_abs`` raises above 120 GPa).  Evaluating
    the differential at some other pressure and adding it to that constant
    mixes anchors: it returns a value that moves RED with pressure (501.5 nm
    at 120 GPa becoming 512.1 nm at 200 GPa for alpha = 0.5), because the
    hydrostatic reference is held at its 120 GPa value while only the
    alpha-dependent part advances.  Compression blue shifts the optimum, so
    that is unphysical.

    To work at another pressure, take the SHIFT alone --
    :func:`deviatoric_shift_nm` -- and add it to a hydrostatic reference
    obtained at that pressure.  Above 120 GPa no such reference exists in the
    published record; ``extrapolation_bounds.py`` gives the band that has to
    stand in for one, and the choice of stand-in must be stated.
    """
    if not 0.0 < alpha <= 1.0:
        raise ValueError("alpha must satisfy 0 < alpha <= 1")
    if not np.isclose(sigma_zz_gpa, HO_ANCHOR_PRESSURE_GPA):
        raise ValueError(
            "the additive anchor %.2f nm is Ho's hydrostatic optimum at "
            "%.0f GPa, so this function is only defined there; got "
            "sigma_zz = %g GPa. Use deviatoric_shift_nm(alpha, sigma_zz_gpa) "
            "and add your own hydrostatic reference for that pressure."
            % (HO_HYDROSTATIC_OPTIMUM_NM, HO_ANCHOR_PRESSURE_GPA,
               sigma_zz_gpa))
    floor = cancellation_alpha(deviatoric_ratio)
    if alpha <= floor:
        raise ValueError(
            "alpha = %.3f is at or below the %.3f where this model's net ZPL "
            "shift changes sign; the two-point calibration cannot be "
            "extrapolated there" % (alpha, floor))
    with _deviatoric_ratio(deviatoric_ratio):
        hydrostatic = nv_model.NVModel(T=temperature_k, alpha=1.0).lambda_opt(
            sigma_zz_gpa
        )
        anisotropic = nv_model.NVModel(T=temperature_k, alpha=alpha).lambda_opt(
            sigma_zz_gpa
        )
    return HO_HYDROSTATIC_OPTIMUM_NM + anisotropic - hydrostatic


def corrected_optimum_range_nm(alpha, sigma_zz_gpa=120.0,
                               temperature_k=300.0):
    """Return (blue_end, red_end) over both readings of the Hilberer ratio.

    The width of this interval is the unresolved normalisation of E5 section 3,
    not a statistical uncertainty.  It collapses to a point the moment the
    experiment measures the ZPL splitting at known alpha (E5 section 9).
    """
    values = [corrected_optimum_nm(alpha, sigma_zz_gpa, temperature_k, ratio)
              for ratio in HILBERER_RATIO_RANGE]
    return (min(values), max(values))


if __name__ == "__main__":
    print("sigma_zz = 120 GPa")
    print("hydrostatic-kernel-only alternative = 440.65 nm (every alpha)")
    print()
    print("alpha   equal-sigma_zz   equal-compression        range")
    for alpha in (0.50, 0.56, 0.60, 0.70, 0.95):
        low, high = corrected_optimum_range_nm(alpha)
        axial = corrected_optimum_nm(
            alpha, deviatoric_ratio=HILBERER_RATIO_EQUAL_AXIAL_STRESS)
        compression = corrected_optimum_nm(
            alpha, deviatoric_ratio=HILBERER_RATIO_EQUAL_COMPRESSION)
        print(f" {alpha:.2f}    {axial:9.2f} nm    {compression:9.2f} nm"
              f"    {low:6.1f}-{high:6.1f}")
