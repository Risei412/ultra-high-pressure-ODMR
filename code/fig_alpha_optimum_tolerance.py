"""Plot conditional 120 GPa optimum and sensitivity-tolerance band vs alpha.

The optimum follows the E5 differential correction.  The band and wavelength
penalty additionally assume that the measured hydrostatic sensitivity curve is
translated rigidly in wavelength, without changing its shape.  This makes the
figure useful for laser selection while keeping the unmeasured assumption
explicit.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from calc_alpha_corrected_optimum import (
    HILBERER_RATIO_EQUAL_AXIAL_STRESS,
    HILBERER_RATIO_EQUAL_COMPRESSION,
    corrected_optimum_nm,
)
import fig_style
from ho_spectrum_model import HBARC, HoPublishedSpectrumModel


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if ROOT.name == "python codes":
    OUT_DIR = ROOT / "figures"
else:
    OUT_DIR = ROOT / "docs" / "experiment" / "figures"
OUT_PNG = OUT_DIR / "alpha_optimum_tolerance_120gpa.png"
OUT_PDF = OUT_DIR / "alpha_optimum_tolerance_120gpa.pdf"

PRESSURE_GPA = 120.0
ALPHA_MIN = 0.50
ALPHA_MAX = 0.95
PRACTICAL_MIN = 0.50
PRACTICAL_MAX = 0.70
TOLERANCE = 1.05
LASER_LINES_NM = (445.0, 457.0, 473.0, 488.0, 505.0, 532.0)

BLUE = "#0072B2"
GREEN = "#009E73"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
BLACK = "#000000"


def hydrostatic_penalty(wavelength_nm, model=None):
    """Sensitivity ratio under the optical-limit hydrostatic kernel."""
    model = model or HoPublishedSpectrumModel()
    wavelength_nm = np.asarray(wavelength_nm, float)
    grid = np.arange(400.0, 850.0001, 0.02)
    objective = grid * model.sigma_abs(HBARC / grid, PRESSURE_GPA)
    maximum = float(np.max(objective))
    sampled = wavelength_nm * model.sigma_abs(
        HBARC / wavelength_nm, PRESSURE_GPA
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(maximum / sampled)


def hydrostatic_tolerance_band(model=None, tolerance=TOLERANCE):
    """Return the contiguous tolerance interval containing the global maximum."""
    model = model or HoPublishedSpectrumModel()
    grid = np.arange(400.0, 600.0001, 0.01)
    penalty = hydrostatic_penalty(grid, model)
    optimum_index = int(np.nanargmin(penalty))
    accepted = penalty <= tolerance
    left = optimum_index
    right = optimum_index
    while left > 0 and accepted[left - 1]:
        left -= 1
    while right + 1 < len(grid) and accepted[right + 1]:
        right += 1
    return float(grid[left]), float(grid[right])


def alpha_band(n=19, deviatoric_ratio=HILBERER_RATIO_EQUAL_COMPRESSION):
    """Conditional optimum and rigid-shifted 5% band over alpha."""
    alpha = np.linspace(ALPHA_MIN, ALPHA_MAX, n)
    optimum = np.array([
        corrected_optimum_nm(value, deviatoric_ratio=deviatoric_ratio)
        for value in alpha
    ])
    hydro_optimum = corrected_optimum_nm(1.0)
    hydro_low, hydro_high = hydrostatic_tolerance_band()
    shift = optimum - hydro_optimum
    return {
        "alpha": alpha,
        "optimum_nm": optimum,
        "lower_5pct_nm": hydro_low + shift,
        "upper_5pct_nm": hydro_high + shift,
        "hydrostatic_band_nm": (hydro_low, hydro_high),
    }


def worst_case_penalty(wavelength_nm, alpha_min=PRACTICAL_MIN,
                       alpha_max=PRACTICAL_MAX, n_alpha=11,
                       deviatoric_ratios=(HILBERER_RATIO_EQUAL_COMPRESSION,
                                          HILBERER_RATIO_EQUAL_AXIAL_STRESS)):
    """Worst penalty over alpha under the rigid wavelength-shift assumption."""
    wavelengths = np.asarray(wavelength_nm, float)
    alpha = np.linspace(alpha_min, alpha_max, n_alpha)
    model = HoPublishedSpectrumModel()
    penalties = []
    for ratio in deviatoric_ratios:
        shifts = np.array([
            corrected_optimum_nm(value, deviatoric_ratio=ratio)
            for value in alpha
        ]) - 440.65
        penalties.extend(
            hydrostatic_penalty(wavelengths - shift, model) for shift in shifts
        )
    penalties = np.vstack(penalties)
    return np.nanmax(penalties, axis=0)


def common_band(data, alpha_min=PRACTICAL_MIN, alpha_max=PRACTICAL_MAX):
    """Intersection of the shifted tolerance bands in an alpha interval."""
    mask = (data["alpha"] >= alpha_min) & (data["alpha"] <= alpha_max)
    lower = float(np.max(data["lower_5pct_nm"][mask]))
    upper = float(np.min(data["upper_5pct_nm"][mask]))
    return lower, upper


def make_figure(data=None, axial_data=None):
    data = data or alpha_band()
    axial_data = axial_data or alpha_band(
        deviatoric_ratio=HILBERER_RATIO_EQUAL_AXIAL_STRESS)
    common_low, common_high = common_band(data)
    axial_low, axial_high = common_band(axial_data)

    fig_style.use()
    fig, axes = plt.subplots(1, 2, figsize=(fig_style.FULL_WIDTH_IN, 4.05))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.34, top=0.91,
                        wspace=0.30)

    ax = axes[0]
    ax.axvspan(PRACTICAL_MIN, PRACTICAL_MAX, color="0.94", zorder=0,
               label=r"common DAC range ($0.5$--$0.7$)")
    ax.fill_between(data["alpha"], data["lower_5pct_nm"],
                    data["upper_5pct_nm"], color=BLUE, alpha=0.22,
                    label=r"5% band: equal $P_{\rm mean}$")
    ax.plot(data["alpha"], data["optimum_nm"], color=BLUE, lw=1.6,
            label=r"$\lambda_{\rm opt}$: equal $P_{\rm mean}$")
    ax.fill_between(axial_data["alpha"], axial_data["lower_5pct_nm"],
                    axial_data["upper_5pct_nm"], color=PURPLE, alpha=0.16,
                    label=r"5% band: equal $\sigma_{zz}$")
    ax.plot(axial_data["alpha"], axial_data["optimum_nm"], color=PURPLE,
            lw=1.25, label=r"$\lambda_{\rm opt}$: equal $\sigma_{zz}$")
    for wavelength in LASER_LINES_NM[:-1]:
        color = VERMILLION if wavelength == 488.0 else "0.60"
        width = 1.15 if wavelength == 488.0 else 0.65
        ax.axhline(wavelength, color=color, lw=width, ls=(0, (4, 2)))
        ax.text(ALPHA_MAX + 0.006, wavelength, f"{wavelength:.0f}",
                va="center", ha="left", fontsize=6.5, color=color,
                clip_on=False)
    ax.set_xlim(ALPHA_MIN, ALPHA_MAX)
    ax.set_ylim(425.0, 522.0)
    ax.set_xlabel(r"Stress ratio $\alpha$")
    ax.set_ylabel("Excitation wavelength (nm)")
    ax.set_title("(a) Optimum and conditional tolerance", loc="left")
    ax.grid(alpha=0.20, lw=0.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2,
              frameon=True, framealpha=1.0,
              facecolor="white", edgecolor="0.8", fancybox=False)

    ax = axes[1]
    wavelengths = np.linspace(450.0, 525.0, 751)
    penalty = worst_case_penalty(wavelengths)
    ax.plot(wavelengths, penalty, color=GREEN, lw=1.6,
            label=r"worst case over $\alpha=0.5$--$0.7$")
    ax.axhline(TOLERANCE, color=BLACK, lw=0.9, ls=(0, (4, 2)),
               label="5% criterion")
    ax.axvspan(common_low, common_high, color=BLUE, alpha=0.16,
               label=f"equal $P_{{mean}}$: {common_low:.1f}--{common_high:.1f} nm")
    ax.axvspan(axial_low, axial_high, color=PURPLE, alpha=0.14,
               label=fr"equal $\sigma_{{zz}}$: {axial_low:.1f}--{axial_high:.1f} nm")
    for wavelength in LASER_LINES_NM:
        if 450.0 <= wavelength <= 525.0:
            value = float(worst_case_penalty([wavelength])[0])
            color = VERMILLION if wavelength == 488.0 else PURPLE
            ax.plot(wavelength, value, "o", ms=4.0, color=color)
            ax.annotate(f"{wavelength:.0f}", (wavelength, value),
                        xytext=(0, 5), textcoords="offset points", ha="center",
                        fontsize=6.5, color=color)
    ax.set_xlim(450.0, 525.0)
    ax.set_ylim(0.995, 1.42)
    ax.set_xlabel("Fixed laser wavelength (nm)")
    ax.set_ylabel(r"Worst-case $\eta/\eta_{\rm opt}$")
    ax.set_title("(b) Fixed-laser decision", loc="left")
    ax.grid(alpha=0.20, lw=0.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25),
              frameon=True, framealpha=1.0,
              facecolor="white", edgecolor="0.8", fancybox=False)
    return fig


def main():
    data = alpha_band()
    axial_data = alpha_band(deviatoric_ratio=HILBERER_RATIO_EQUAL_AXIAL_STRESS)
    common_low, common_high = common_band(data)
    axial_low, axial_high = common_band(axial_data)
    fig = make_figure(data, axial_data)
    print(fig_style.save(fig, OUT_PNG, OUT_PDF))
    print(f"conditional common 5% band at alpha=0.5--0.7: "
          f"{common_low:.2f}--{common_high:.2f} nm")
    print(f"alternative-normalisation common 5% band: "
          f"{axial_low:.2f}--{axial_high:.2f} nm")
    for wavelength in LASER_LINES_NM:
        value = float(worst_case_penalty([wavelength])[0])
        print(f"{wavelength:6.1f} nm: worst-case penalty x{value:.4f}")
    print("worst case includes both Hilberer normalisations")
    print("assumption: hydrostatic sensitivity curve is rigidly shifted")


if __name__ == "__main__":
    main()
