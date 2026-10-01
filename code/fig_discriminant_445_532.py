"""Grant Fig. 2: the measurable discriminant  D = eta(445 nm) / eta(532 nm).

SIGN CONVENTION.  eta is a SENSITIVITY: smaller is better.  Therefore
    D < 1  =>  445 nm is the more sensitive line (blue favoured)
    D > 1  =>  532 nm is the more sensitive line (green favoured)
A "penalty" quoted anywhere in this project is eta(lambda)/eta(best) >= 1, so a
LARGER penalty is WORSE.  Keep both conventions straight when reading the plot.

THE TEST.  At 120 GPa the two frozen hypotheses predict:

    hydrostatic-kernel-fixed (H2)   lambda_opt = 441 nm for every alpha
                                    ->  D = 0.080   (alpha-independent)
    deviatoric-corrected    (H2')   lambda_opt = 506-497 nm at alpha 0.5-0.7
                                    ->  D = 2.63-1.41

The decision boundary is NOT 1.0.  H2' dips below 1 once alpha exceeds ~0.8
(and at 200 GPa it does so already at alpha ~ 0.72), so "is D above or below 1"
does not separate the hypotheses.  The boundary used here is the geometric mean
of the two predictions at the least favourable point of the expected range
(alpha = 0.7, 120 GPa):

    D_threshold = sqrt(0.0801 * 1.414) = 0.337  ->  quoted as 0.34

At 120 GPa this separates the hypotheses over the WHOLE alpha range 0.5-1.0
(H2' never falls below 0.409), and inside the expected range alpha = 0.5-0.7 it
carries a margin of at least x4.2 on both sides.

Being a RATIO of two lines measured on the same sample in the same run, D
cancels absolute collection efficiency, NV density and incident power
calibration -- the quantities hardest to control in a DAC.

VALIDITY.  The reconstructed Ho kernel raises above 120 GPa, so H2 can only be
evaluated at 120 GPa: THE DECISIVE TEST IS AT 120 GPa.  H2' is also shown at
200 GPa, but only to extend the lambda_opt(alpha, P) map -- it is NOT a second
test, because H2 has no computable counterpart there.  Above 120 GPa NVModel
runs clipped (Sabs, DWF and IP held at their 120 GPa values, only the ZPL shift
moves), so that curve is a design envelope, not a prediction.

Run:  python fig_discriminant_445_532.py
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
OUT_DIR = HERE.parent / "docs" / "experiment" / "figures"
OUT_PNG = OUT_DIR / "discriminant_445_532.png"
OUT_PDF = OUT_DIR / "discriminant_445_532.pdf"

BLUE_LINE = 445.0
GREEN_LINE = 532.0
ALPHAS = np.linspace(0.50, 1.00, 51)
EXPECTED_LO, EXPECTED_HI = 0.50, 0.70   # what a DAC sample chamber reaches
PRESSURES = (120.0, 200.0)

BLUE = "#0072B2"
GREEN = "#009E73"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
GREY = "0.35"


def h2_ratio():
    """Hydrostatic-kernel-fixed: alpha-independent, 120 GPa only."""
    model = HoIntegratedODMRModel()
    blue = model.absorbed_photon_proxy(BLUE_LINE, 120.0)
    green = model.absorbed_photon_proxy(GREEN_LINE, 120.0)
    return float(np.sqrt(green / blue))


def h2prime_ratio(pressure):
    """Deviatoric-corrected single-mode envelope, as a function of alpha."""
    out = []
    for alpha in ALPHAS:
        model = NVModel(alpha=float(alpha))
        blue = model.eta_lambda(BLUE_LINE, pressure)[0]
        green = model.eta_lambda(GREEN_LINE, pressure)[0]
        out.append(float(blue / green))
    return np.asarray(out)


def make_figure(h2, curves, threshold):
    fig_style.use()
    fig, ax = plt.subplots(figsize=(fig_style.FULL_WIDTH_IN, 3.7))
    fig.subplots_adjust(left=0.098, right=0.988, bottom=0.128, top=0.775)

    ax.axvspan(EXPECTED_LO, EXPECTED_HI, color=GREY, alpha=0.11, lw=0,
               label="expected DAC sample-chamber range")

    # Reference line: equal sensitivity.  NOT the decision boundary.
    ax.axhline(1.0, color=GREY, lw=0.8, ls=(0, (1, 2)))
    ax.text(0.998, 1.05, r"$\eta$ equal: 445 nm and 532 nm as sensitive",
            ha="right", va="bottom", fontsize=fig_style.ANNOT, color=GREY)

    # The actual decision boundary.
    ax.axhline(threshold, color=GREEN, lw=1.6, ls=(0, (6, 2)))
    ax.text(0.998, threshold * 1.10,
            f"decision boundary  $D={threshold:.2f}$",
            ha="right", va="bottom", fontsize=fig_style.ANNOT, color=GREEN)

    ax.axhline(h2, color=BLUE, lw=1.9,
               label="H2  hydrostatic kernel fixed (120 GPa)")
    ax.text(0.998, h2 * 0.83, f"$D={h2:.3f}$  (445 nm $\\times${1.0 / h2:.1f} more sensitive)",
            ha="right", va="top", fontsize=fig_style.ANNOT, color=BLUE)

    styles = {120.0: ("-", VERMILLION, 1.9, "H2$'$  deviatoric-corrected, 120 GPa"),
              200.0: ((0, (5, 2)), PURPLE, 1.4,
                      "H2$'$  200 GPa (map only: H2 not evaluable)")}
    for pressure, values in curves.items():
        style, colour, width, label = styles[pressure]
        ax.plot(ALPHAS, values, color=colour, lw=width, ls=style, label=label)

    ax.set_yscale("log")
    ax.set_xlim(0.50, 1.00)
    ax.set_ylim(0.05, 6.0)
    ax.set_xlabel(r"Stress anisotropy  $\alpha=\sigma_{xx}/\sigma_{zz}$"
                  "   (1 = hydrostatic)")
    ax.set_ylabel(r"Discriminant  $D=\eta(445\,\mathrm{nm})/\eta(532\,\mathrm{nm})$"
                  "\n" r"($\eta$ = sensitivity, smaller is better)")
    ax.grid(alpha=0.18, lw=0.5, which="both")

    ax.annotate("445 nm more sensitive", xy=(0.512, 0.098),
                fontsize=fig_style.ANNOT, color=BLUE, ha="left", va="center")
    ax.annotate("532 nm more sensitive", xy=(0.512, 4.3),
                fontsize=fig_style.ANNOT, color=VERMILLION, ha="left",
                va="center")

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.215), ncol=2,
              frameon=True, framealpha=1.0, facecolor="white",
              edgecolor="0.8", fancybox=False, handlelength=2.0,
              columnspacing=1.2, fontsize=fig_style.ANNOT)
    fig.text(0.098, 0.988,
             "The decisive test at 120 GPa: one ratio, separated by "
             r"$\geq\times4.2$ on both sides of $D=0.34$",
             ha="left", va="top", fontsize=9.0)
    return fig


def main():
    h2 = h2_ratio()
    curves = {p: h2prime_ratio(p) for p in PRESSURES}

    window = (ALPHAS >= EXPECTED_LO) & (ALPHAS <= EXPECTED_HI)
    worst = float(curves[120.0][window].min())          # alpha = 0.7, 120 GPa
    threshold = float(np.sqrt(h2 * worst))

    assert abs(h2 - 0.0801) < 5e-4, h2
    # The boundary must sit strictly between the two predictions, with equal
    # log-margin at the least favourable point of the expected range.
    assert h2 < threshold < worst, (h2, threshold, worst)
    assert abs((worst / threshold) - (threshold / h2)) < 1e-6
    assert worst / threshold > 4.0, worst / threshold
    # And it must still separate them for every alpha at 120 GPa, not only in
    # the expected window -- otherwise a better-than-expected cell breaks it.
    assert curves[120.0].min() > threshold, curves[120.0].min()

    print(fig_style.save(fig := make_figure(h2, curves, threshold),
                         OUT_PNG, OUT_PDF))
    plt.close(fig)
    print()
    print("eta is a SENSITIVITY: smaller is better.")
    print("  D < 1 -> 445 nm more sensitive;  D > 1 -> 532 nm more sensitive.")
    print()
    print(f"H2  (kernel fixed, 120 GPa, alpha-independent): D = {h2:.4f}")
    print(f"Decision boundary D = sqrt({h2:.4f} * {worst:.3f}) = {threshold:.3f}"
          f"  -> quoted as {threshold:.2f}")
    print()
    print(f"{'alpha':>7}  {'H2prime 120':>12}  {'vs H2':>9}  {'vs bound':>9}"
          f"  {'H2prime 200':>12}")
    for alpha in (0.50, 0.60, 0.70, 0.80, 0.90, 1.00):
        index = int(np.abs(ALPHAS - alpha).argmin())
        a120, a200 = curves[120.0][index], curves[200.0][index]
        print(f"{alpha:7.2f}  {a120:12.3f}  {a120 / h2:8.1f}x  "
              f"{a120 / threshold:8.2f}x  {a200:12.3f}")
    print()
    print("The decisive test is at 120 GPa (H2 cannot be evaluated above it).")
    print("Inside the expected range alpha = 0.5-0.7 the margin is at least")
    print(f"x{worst / threshold:.1f} on both sides of D = {threshold:.2f}.")
    print("D = 1 is NOT the boundary: H2' falls below 1 for alpha > ~0.8.")


if __name__ == "__main__":
    main()
