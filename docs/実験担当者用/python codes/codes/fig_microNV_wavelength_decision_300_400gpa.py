"""Excitation-wavelength design envelope for micro-NV at 300--400 GPa.

This figure deliberately does not plot sensitivity.  The available optical
kernel is hydrostatic and anchored only to 120 GPa; above that pressure the
calculation can provide an extrapolation envelope and an ionization constraint,
but not a validated anisotropic-stress sensitivity optimum.
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
from extrapolation_bounds import VARIANTS, build, lambda_ion


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if ROOT.name == "python codes":
    OUT_DIR = ROOT / "figures"
else:
    OUT_DIR = ROOT / "docs" / "experiment" / "figures"
OUT_PNG = OUT_DIR / "microNV_wavelength_decision_300_400gpa.png"
OUT_PDF = OUT_DIR / "microNV_wavelength_decision_300_400gpa.pdf"

PRESSURES = (300.0, 400.0)
ALPHAS = np.linspace(0.60, 0.70, 21)
RATIOS = (HILBERER_RATIO_EQUAL_COMPRESSION, HILBERER_RATIO_EQUAL_AXIAL_STRESS)
LASERS = (441.0, 457.0, 488.0, 532.0)

BLUE = "#0072B2"
GREEN = "#009E73"
VERMILLION = "#D55E00"


def wavelength_envelope(pressure):
    values = []
    for _, kwargs, consistent in VARIANTS:
        if not consistent:
            continue
        reference = float(build(**kwargs).lambda_opt(pressure, lo=280.0, hi=700.0))
        for ratio in RATIOS:
            values.extend(
                reference + deviatoric_shift_nm(float(alpha), pressure,
                                                deviatoric_ratio=ratio)
                for alpha in ALPHAS
            )
    return min(values), max(values)


def make_figure():
    envelopes = {P: wavelength_envelope(P) for P in PRESSURES}
    walls = {
        P: float(lambda_ion(build(), P)) for P in PRESSURES
    }
    linear_walls = {
        P: float(lambda_ion(build(zpl="linear", anchors="linear"), P))
        for P in PRESSURES
    }

    fig_style.use()
    fig, axes = plt.subplots(2, 1, figsize=(fig_style.FULL_WIDTH_IN, 3.65),
                             sharex=True)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.17, top=0.80,
                        hspace=0.70)

    for ax, pressure in zip(axes, PRESSURES):
        lo, hi = envelopes[pressure]
        wall = walls[pressure]
        linear_wall = linear_walls[pressure]
        ax.axvspan(300.0, wall, color=VERMILLION, alpha=0.10,
                   label="ionization-risk region" if pressure == 300.0 else None)
        ax.axvline(wall, color=VERMILLION, lw=1.1,
                   label=r"conservative wall ($\lambda_{\rm ion}$)" if pressure == 300.0 else None)
        ax.axvline(linear_wall, color=VERMILLION, lw=0.9, ls=(0, (2, 2)),
                   label="linear-continuation wall" if pressure == 300.0 else None)
        ax.plot([lo, hi], [0, 0], color=BLUE, lw=10, solid_capstyle="butt",
                label=r"model envelope ($\alpha=0.6$--$0.7$)" if pressure == 300.0 else None)
        safe_lo = max(lo, wall)
        ax.plot([safe_lo, hi], [0, 0], color=GREEN, lw=4,
                solid_capstyle="butt",
                label="conservative usable part" if pressure == 300.0 else None)
        ax.text(0.5 * (lo + hi), 0.15, f"{lo:.0f}--{hi:.0f} nm",
                ha="center", va="bottom", color=BLUE,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.90,
                          pad=0.8))
        for laser in LASERS:
            if laser == 488.0:
                color, lw = GREEN, 1.9
            else:
                color, lw = "0.45", 0.85
            ax.axvline(laser, color=color, lw=lw, ls=(0, (4, 2)))
            ax.text(laser, -0.14, f"{laser:.0f}", ha="center", va="top",
                    color=color, fontsize=7)
        ax.text(wall + 2.0, -0.14, f"{wall:.0f} wall", ha="left", va="top",
                color=VERMILLION, fontsize=7)
        ax.text(0.01, 0.62, f"{pressure:.0f} GPa", transform=ax.transAxes,
                fontsize=9, fontweight="bold")
        ax.text(linear_wall + 2.0, -0.14,
                f"{linear_wall:.0f} linear wall", ha="left", va="top",
                color=VERMILLION, fontsize=7)
        ax.set_ylim(-0.24, 0.52)
        ax.set_yticks([])
        ax.grid(axis="x", alpha=0.20, lw=0.5)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, 1.17), ncol=3,
                    frameon=True, framealpha=1.0, facecolor="white",
                    edgecolor="0.8", fancybox=False)
    axes[1].set_xlabel("Excitation wavelength (nm)")
    axes[1].set_xlim(300.0, 545.0)
    fig.text(0.12, 0.995,
             r"Micro-NV in DAC sample chamber: 300--400 GPa design envelope",
             ha="left", va="top", fontsize=9)
    return fig, envelopes, walls, linear_walls


def main():
    fig, envelopes, walls, linear_walls = make_figure()
    print(fig_style.save(fig, OUT_PNG, OUT_PDF))
    for pressure in PRESSURES:
        lo, hi = envelopes[pressure]
        print(f"{pressure:.0f} GPa: model envelope = {lo:.1f}--{hi:.1f} nm; "
              f"conservative usable part = {max(lo, walls[pressure]):.1f}--{hi:.1f} nm; "
              f"walls = {walls[pressure]:.1f} / {linear_walls[pressure]:.1f} nm")
    print("recommended common operating laser: 488 nm")
    print("sensitivity: not calculated; anisotropic 300--400 GPa kernel is unavailable")


if __name__ == "__main__":
    main()

