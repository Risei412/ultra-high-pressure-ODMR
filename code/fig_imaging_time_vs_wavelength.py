"""Fig. 5 (grant): relative wide-field image acquisition time vs excitation wavelength.

In wide-field ODMR magnetic imaging the integration time needed to reach a
fixed SNR scales as the square of the magnetic sensitivity,

    T / T_min = ( eta(lambda) / eta_min )^2 ,

and in the optical (shot-noise) limit, where contrast and linewidth carry no
wavelength dependence, eta ~ 1/sqrt(R), so the same quantity is just the
reciprocal relative absorbed-photon rate R_opt / R(lambda).

NOTHING NEW IS COMPUTED HERE.  The curves are the frozen 120 GPa hydrostatic
(alpha = 1) wavelength envelope, re-plotted on a time axis.

UNITS WARNING (this is the trap this figure exists to avoid).  The number
"157" quoted for 532 nm in PLAN.md and in
theory_limits_and_required_measurements.md is ALREADY the squared quantity:

    relative absorbed rate at 532 nm = 0.006371
    eta   penalty  = sqrt(1/0.006371) =  12.53   <- SENSITIVITY
    eta^2 penalty  =      1/0.006371  = 157.0    <- ACQUISITION TIME

So 157 is the acquisition-time penalty, not the sensitivity penalty, and
squaring it again (157^2 = 2.5e4) double-counts.  See main() for the assertion
that pins this.

Run:  python fig_imaging_time_vs_wavelength.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fig_style
from ho_odmr_sensitivity import HoIntegratedODMRModel
from nv_model import NVModel

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT_DIR = ROOT / "docs" / "experiment" / "figures"
OUT_PNG = OUT_DIR / "imaging_time_vs_wavelength_120gpa.png"
OUT_PDF = OUT_DIR / "imaging_time_vs_wavelength_120gpa.pdf"

PRESSURE = 120.0
LAM = np.arange(400.0, 560.0 + 0.05, 0.05)

# Model-spread band for the optimum (B1 diagnostic, 2026-09-05).
BAND_LO, BAND_HI = 441.0, 465.0
KERNEL_OPT = 440.65        # reconstructed Ho kernel
SINGLE_OPT = 475.5         # single-mode envelope
ION_WALL = 405.2           # hc / IP(3A2) at 120 GPa

# The narrow notch near 514.5 nm is the ZPL absorption line of the kernel, not
# numerical noise: the kernel's maximum sits ON the ZPL branch below ~80 GPa and
# crosses to the phonon sideband by ~100 GPa, so at 120 GPa the ZPL survives as a
# narrow resonance on the far red flank.  Its penalty (x1.21 in eta) matches the
# x1.2006 quoted in theory_limits_and_required_measurements.md.  It is annotated
# because a reviewer would otherwise read it either as noise or as "there is a
# good wavelength next to 532 nm after all" -- it is too narrow to operate on.
ZPL_NOTCH = 514.45

LASERS = ((445.0, "445"), (487.0, "487"), (532.0, "532"))

BLUE = "#0072B2"
GREEN = "#009E73"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"


def kernel_time():
    """Relative acquisition time from the reconstructed Ho kernel."""
    model = HoIntegratedODMRModel()
    rate = np.array([model.absorbed_photon_proxy(float(w), PRESSURE)
                     for w in LAM])
    return rate.max() / rate


def single_mode_time():
    """Relative acquisition time from the single-mode envelope, alpha = 1."""
    eta = np.asarray(NVModel(alpha=1.0).eta_lambda(LAM, PRESSURE)[0])
    return (eta / eta.min()) ** 2


def make_figure(t_kernel, t_single):
    fig_style.use()
    fig, ax = plt.subplots(figsize=(fig_style.FULL_WIDTH_IN, 3.35))
    fig.subplots_adjust(left=0.088, right=0.988, bottom=0.135, top=0.795)

    ax.axvspan(400.0, ION_WALL, color=VERMILLION, alpha=0.10)
    ax.axvline(ION_WALL, color=VERMILLION, lw=1.0, ls="-")
    ax.text(ION_WALL + 2.0, 2.6e2, "ionization wall 405.2 nm",
            color=VERMILLION, fontsize=fig_style.ANNOT, ha="left", va="top")

    ax.axvspan(BAND_LO, BAND_HI, color=BLUE, alpha=0.15, lw=0,
               label=f"model spread of $\\lambda_{{\\rm opt}}$ ({BAND_LO:.0f}--{BAND_HI:.0f} nm)")

    ax.plot(LAM, t_kernel, color=BLUE, lw=1.7,
            label="reconstructed Ho kernel")
    ax.plot(LAM, t_single, color=PURPLE, lw=1.7, ls=(0, (5, 2)),
            label="single-mode envelope")

    for value, color, label in ((KERNEL_OPT, BLUE, None),
                                (SINGLE_OPT, PURPLE, None)):
        ax.plot([value], [1.0], marker="v", color=color, ms=5.5,
                clip_on=False, zorder=6, label=label)
    ax.annotate(f"{KERNEL_OPT:.1f}", xy=(KERNEL_OPT, 1.0), xytext=(-9, -13),
                textcoords="offset points", color=BLUE,
                fontsize=fig_style.ANNOT, ha="right", va="center")
    ax.annotate(f"{SINGLE_OPT:.1f}", xy=(SINGLE_OPT, 1.0), xytext=(9, -13),
                textcoords="offset points", color=PURPLE,
                fontsize=fig_style.ANNOT, ha="left", va="center")

    # (x offset, y offset) in points, hand-placed to clear the curves.
    offsets = {445.0: (-6, 26), 487.0: (-4, 12), 532.0: (0, 13)}
    for wavelength, tag in LASERS:
        index = int(np.abs(LAM - wavelength).argmin())
        y = t_kernel[index]
        ax.axvline(wavelength, color="0.45", lw=0.8, ls=(0, (2, 2)), zorder=1)
        ax.plot([wavelength], [y], marker="o", ms=4.5, color=GREEN,
                markeredgecolor="white", markeredgewidth=0.6, zorder=7)
        text = f"{tag} nm\n$\\times${y:.0f}" if y >= 10 else f"{tag} nm\n$\\times${y:.2f}"
        ax.annotate(text, xy=(wavelength, y), xytext=offsets[wavelength],
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=fig_style.ANNOT, color=GREEN, zorder=8,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85,
                              pad=0.6))

    notch_index = int(np.abs(LAM - ZPL_NOTCH).argmin())
    notch_y = float(t_kernel[notch_index])
    ax.annotate(
        f"ZPL line {ZPL_NOTCH:.1f} nm ($\\times${notch_y:.2f})\n"
        "real resonance, ~1 nm wide:\nnot a usable operating point",
        xy=(ZPL_NOTCH, notch_y), xytext=(0, -14),
        textcoords="offset points", ha="center", va="top",
        fontsize=fig_style.ANNOT, color="0.30", zorder=8,
        arrowprops=dict(arrowstyle="-", color="0.45", lw=0.8,
                        shrinkA=0.0, shrinkB=2.0),
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.6))

    ax.set_yscale("log")
    ax.set_xlim(400.0, 560.0)
    ax.set_ylim(0.22, 5e2)   # headroom below the curves for the ZPL-notch note
    ax.set_xlabel("Excitation wavelength (nm)")
    ax.set_ylabel(r"Relative acquisition time  $T/T_{\min}=(\eta/\eta_{\min})^2$")
    ax.grid(alpha=0.20, lw=0.5, which="both")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.175), ncol=3,
              frameon=True, framealpha=1.0, facecolor="white",
              edgecolor="0.8", fancybox=False, handlelength=1.8,
              columnspacing=1.2)
    fig.text(0.088, 0.985,
             r"Wide-field NV magnetic imaging at 120 GPa (hydrostatic, $\alpha=1$)",
             ha="left", va="top", fontsize=9.0)
    return fig


def main():
    t_kernel, t_single = kernel_time(), single_mode_time()

    def at(array, wavelength):
        return float(array[int(np.abs(LAM - wavelength).argmin())])

    # Pin the units: 157 is the TIME penalty, 12.53 is the SENSITIVITY penalty.
    time_532 = at(t_kernel, 532.0)
    assert abs(time_532 - 157.0) < 1.0, time_532
    assert abs(np.sqrt(time_532) - 12.53) < 0.05, np.sqrt(time_532)

    print(fig_style.save(fig := make_figure(t_kernel, t_single), OUT_PNG, OUT_PDF))
    plt.close(fig)
    print()
    print(f"120 GPa, hydrostatic (alpha = 1).  Grid {LAM[0]:.0f}-{LAM[-1]:.0f} nm.")
    print(f"{'lambda':>8}  {'T/Tmin kernel':>14}  {'T/Tmin single':>14}  "
          f"{'eta pen kernel':>15}")
    for wavelength in (405.2, 441.0, 445.0, 457.0, 465.0, 475.5, 487.0, 532.0):
        tk, ts = at(t_kernel, wavelength), at(t_single, wavelength)
        print(f"{wavelength:8.1f}  {tk:14.2f}  {ts:14.2f}  {np.sqrt(tk):15.2f}")
    print()
    print(f"532 nm: acquisition time x{time_532:.1f}, sensitivity "
          f"x{np.sqrt(time_532):.2f}.")
    print("The '157' in PLAN.md is the TIME penalty; the sensitivity penalty "
          "is 12.5.")


if __name__ == "__main__":
    main()
