"""Decision plot for micro-NV excitation at 120 and 200 GPa.

The 120 GPa interval is the conditional alpha=0.6--0.7 correction to the
440.65 nm hydrostatic Ho kernel.  The 200 GPa interval combines that correction
with the consistent extrapolation variants in extrapolation_bounds.py.  It is
an experimental design envelope, not a confidence interval.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fig_style
from calc_alpha_corrected_optimum import (
    HILBERER_RATIO_EQUAL_AXIAL_STRESS,
    HILBERER_RATIO_EQUAL_COMPRESSION,
    deviatoric_shift_nm,
)
from extrapolation_bounds import build, lambda_ion, VARIANTS


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if ROOT.name == "python codes":
    OUT_DIR = ROOT / "figures"
else:
    OUT_DIR = ROOT / "docs" / "experiment" / "figures"
OUT_PNG = OUT_DIR / "microNV_wavelength_decision_120_200gpa.png"
OUT_PDF = OUT_DIR / "microNV_wavelength_decision_120_200gpa.pdf"

ALPHAS = np.linspace(0.60, 0.70, 21)
RATIOS = (HILBERER_RATIO_EQUAL_COMPRESSION, HILBERER_RATIO_EQUAL_AXIAL_STRESS)
LASERS = (441.0, 457.0, 488.0, 532.0)

BLUE = "#0072B2"
GREEN = "#009E73"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"


def alpha_range_120():
    values = []
    for ratio in RATIOS:
        values.extend(
            440.65 + deviatoric_shift_nm(float(a), 120.0,
                                          deviatoric_ratio=ratio)
            for a in ALPHAS
        )
    return min(values), max(values)


def alpha_range_200():
    values = []
    for _, kwargs, consistent in VARIANTS:
        if not consistent:
            continue
        hydro = build(**kwargs)
        reference = float(hydro.lambda_opt(200.0, lo=280.0, hi=700.0))
        for ratio in RATIOS:
            values.extend(
                reference + deviatoric_shift_nm(float(a), 200.0,
                                                deviatoric_ratio=ratio)
                for a in ALPHAS
            )
    return min(values), max(values)


def make_figure():
    ranges = {120: alpha_range_120(), 200: alpha_range_200()}
    # The conservative wall is the longest wavelength among the two allowed
    # extrapolations: saturating IP gives 405.2 nm at both pressures.
    walls = {120: 405.2, 200: 405.2}
    linear_wall_200 = float(lambda_ion(build(zpl="linear", anchors="linear"), 200.0))

    fig_style.use()
    fig, axes = plt.subplots(2, 1, figsize=(fig_style.FULL_WIDTH_IN, 3.9),
                             sharex=True)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.16, top=0.78,
                        hspace=0.72)

    for ax, pressure in zip(axes, (120, 200)):
        lo, hi = ranges[pressure]
        ax.axvspan(390, walls[pressure], color=VERMILLION, alpha=0.10,
                   label="photoionization-risk region" if pressure == 120 else None)
        ax.axvline(walls[pressure], color=VERMILLION, lw=1.0,
                   label=r"ionization wall ($\lambda_{\rm ion}$)" if pressure == 120 else None)
        ax.plot([lo, hi], [0, 0], color=BLUE, lw=10, solid_capstyle="butt",
                label=r"conditional $\alpha=0.6$--$0.7$ range" if pressure == 120 else None)
        ax.plot([lo, hi], [0, 0], color="white", lw=5, solid_capstyle="butt")
        ax.plot([lo, hi], [0, 0], color=BLUE, lw=2.0, solid_capstyle="butt")
        ax.text(0.5 * (lo + hi), 0.14, f"{lo:.0f}--{hi:.0f} nm",
                ha="center", va="bottom", color=BLUE,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.90,
                          pad=0.8))
        for laser in LASERS:
            color = GREEN if laser == 488.0 else "0.45"
            lw = 1.8 if laser == 488.0 else 0.9
            ax.axvline(laser, color=color, lw=lw, ls=(0, (4, 2)))
            ax.text(laser, -0.14, f"{laser:.0f}", ha="center", va="top",
                    color=color, fontsize=7)
        ax.text(0.01, 0.78, f"{pressure} GPa", transform=ax.transAxes,
                fontsize=9, fontweight="bold")
        if pressure == 200:
            ax.text(0.99, 0.78, f"linear wall: {linear_wall_200:.0f} nm",
                    transform=ax.transAxes, ha="right", color=VERMILLION,
                    fontsize=7)
        ax.set_ylim(-0.25, 0.55)
        ax.set_yticks([])
        ax.grid(axis="x", alpha=0.20, lw=0.5)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=3,
                    frameon=True, framealpha=1.0, facecolor="white",
                    edgecolor="0.8", fancybox=False)
    axes[1].set_xlabel("Excitation wavelength (nm)")
    axes[1].set_xlim(390, 545)
    fig.text(0.12, 0.995,
             r"Micro-NV in DAC sample chamber: wavelength decision envelope",
             ha="left", va="top", fontsize=9)
    return fig, ranges, linear_wall_200


def main():
    fig, ranges, linear_wall_200 = make_figure()
    print(fig_style.save(fig, OUT_PNG, OUT_PDF))
    for pressure, values in ranges.items():
        print(f"{pressure} GPa: conditional alpha range = "
              f"{values[0]:.1f}--{values[1]:.1f} nm")
    print(f"200 GPa linear-continuation wall = {linear_wall_200:.1f} nm")
    print("recommended common laser for the 200 GPa target: 488 nm")


if __name__ == "__main__":
    main()
