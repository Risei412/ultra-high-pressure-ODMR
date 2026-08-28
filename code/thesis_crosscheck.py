"""Retrodictions against Bhattacharyya's thesis (UC Berkeley, 2022).

`docs/ref/Principle_and_Applications_of_.pdf`, chapter 6, contains two
wavelength comparisons that were made before this theory existed and that
nothing in the v3 chain was tuned against.  They are the only excitation-
wavelength data in the repository that is independent of Ho et al., so they
are the sharpest external test available.

T1  Fig. 6.3.  532 nm against 450 nm in a [100] cut culet at ~50 GPa, scored
    as SNR ~ contrast * sqrt(counts).  The finding is a NULL result: "no
    distinct advantage of working with higher energy excitation at high
    pressures".  In the optical limit our kernel gives an SNR ratio of
    sqrt(A_450/A_532), and at 50 GPa that is 0.94 -- the two lines are within
    6% of each other because 50 GPa is essentially the crossover, which the
    kernel puts at 51.8 GPa.  The null result is predicted, not merely
    accommodated: it is what you get for measuring at the one pressure where
    the answer is a tie.

T2  Fig. 6.2.  405 nm against 532 nm in hydrostatically compressed
    microdiamonds: no clear resonance under 405 nm at ~17 GPa, a clear one at
    ~100 GPa.  The kernel gives an absorption ratio A_405/A_532 of 0.007 at
    17 GPa and 54.6 at 100 GPa, with the crossover at 74.4 GPa.  Four orders
    of magnitude between the two pressures is why one works and the other
    does not.

T3  Section 6.1, note.  Doherty et al. projected the ZPL to cross 532 nm at
    ~66 GPa, "potentially precluding standard excitation schemes".  Dai et al.
    then ran 532 nm cwODMR to 140 GPa.  The kernel resolves this without
    appealing to the nonlinear ZPL shift the thesis cites: absorption at a
    fixed line is a SIDEBAND process, so it does not switch off when the ZPL
    crosses the line.  A_532 falls to 0.6% of the band maximum by 100 GPa but
    stays finite and positive.  (This is the same physics as erratum E1 and
    revision C-1, which removed a hard cut below the ZPL from v1.)

Erratum E5 (2026-08-28) changed the figure of merit used here.  This module
originally compared bare cross sections, sigma(lambda), which is the fixed
PHOTON FLUX convention.  The thesis states that its two lines were compared at
"similar laser and microwave powers", i.e. at fixed OPTICAL POWER, where the
absorbed rate carries the extra photon-energy factor and the figure of merit is
A = lambda sigma -- the same A the frozen record uses everywhere else.  Under
the correct convention the Fig. 6.3 crossover moves from 51.82 to 54.37 GPa,
against 54.4 GPa read off Fig. 6.3(b) itself.  See the freeze, Erratum E5.

Caveat carried by all three: the thesis scores SNR as contrast * sqrt(counts),
which drops the linewidth.  In A2's decomposition E = 2c + s + 2w the
linewidth carries 2w = 1 of E = 3, a third of the exponent, so their metric
and eta = dnu / (C sqrt(R)) are not the same quantity at finite power.  The
comparisons below are made in the optical limit, where they coincide.

Writes nothing; run it for the table.
"""
import numpy as np
from scipy.optimize import brentq

from ho_spectrum_model import HoPublishedSpectrumModel
from nv_model import nm2eV

REFERENCE_NM = 532.0
# The pressures and lines the thesis actually reports.
THESIS_FIG63_GPA = 50.0
THESIS_FIG63_LINE = 450.0
THESIS_FIG62_GPA = (17.0, 100.0)
THESIS_FIG62_LINE = 405.0
# Highest reference curve of Fig. 1(e) whose real samples still bracket 532 nm.
# Above it erratum E1's empty sub-ZPL interval swallows the line (the ZPL
# passes 532 nm at 92.1 GPa), so A_532 there is interpolation, not extraction.
E1_LAST_SAMPLED_GPA = 80.0


def absorption(model, lam, pressure):
    """Bare cross section: the fixed-photon-flux figure of merit."""
    return float(model.sigma_abs(nm2eV(lam), float(pressure)))


def merit(model, lam, pressure):
    """A = lambda sigma: the fixed-OPTICAL-POWER figure of merit (E5).

    The thesis compared its lines at "similar laser and microwave powers", so
    this, not `absorption`, is the quantity to compare with its SNR.
    """
    return lam * absorption(model, lam, pressure)


def ratio(model, lam, pressure, reference=REFERENCE_NM):
    """A(lam)/A(reference) at one pressure, at fixed optical power (E5)."""
    return merit(model, lam, pressure) / merit(model, reference, pressure)


def snr_ratio(model, lam, pressure, reference=REFERENCE_NM):
    """Optical-limit SNR advantage of `lam` over `reference`."""
    return float(np.sqrt(ratio(model, lam, pressure, reference)))


def crossover(model, lam, reference=REFERENCE_NM, bracket=(1.0, 119.0)):
    """Pressure at which `lam` overtakes `reference`, in GPa."""
    return float(brentq(lambda p: ratio(model, lam, p, reference) - 1.0,
                        *bracket))


def report():
    model = HoPublishedSpectrumModel()
    out = {}

    print('T1  Fig. 6.3 -- 450 nm vs 532 nm, [100] culet, ~50 GPa')
    out['fig63_snr_ratio'] = snr_ratio(model, THESIS_FIG63_LINE,
                                       THESIS_FIG63_GPA)
    out['fig63_crossover'] = crossover(model, THESIS_FIG63_LINE)
    print(f'      thesis: no distinct advantage for 450 nm')
    print(f'      kernel: SNR ratio {out["fig63_snr_ratio"]:.2f}x, '
          f'crossover at {out["fig63_crossover"]:.1f} GPa')
    print('      the measurement sits on the crossover, so a tie is the '
          'prediction')
    print(f'      {"P [GPa]":>9}{"450/532 SNR":>14}')
    for pressure in (30, 40, 50, 55, 60, 65, 70):
        print(f'      {pressure:9d}{snr_ratio(model, 450.0, pressure):13.2f}x')
    print('      they compressed to ~70 GPa and report no change above '
          '~50 GPa;\n      the kernel predicts 1.8x by 70 GPa, so this is '
          'only half a hit --\n      see the anvil-transmission gap in the '
          'scope document.')

    print('\nT2  Fig. 6.2 -- 405 nm vs 532 nm, hydrostatic microdiamonds')
    for pressure in THESIS_FIG62_GPA:
        out[f'fig62_ratio_{pressure:.0f}'] = ratio(model, THESIS_FIG62_LINE,
                                                   pressure)
        print(f'      {pressure:5.0f} GPa: A_405/A_532 = '
              f'{out[f"fig62_ratio_{pressure:.0f}"]:.3f}')
    out['fig62_crossover'] = crossover(model, THESIS_FIG62_LINE)
    print(f'      crossover at {out["fig62_crossover"]:.1f} GPa')
    print('      thesis: no resonance at 17 GPa, clear resonance at 100 GPa')

    print('\nT3  Section 6.1 -- does 532 nm die when the ZPL crosses it?')
    print(f'      {"P [GPa]":>9}{"A_532(P)/A_532(0)":>20}   provenance')
    ambient = absorption(model, REFERENCE_NM, 0.0)
    for pressure in (0, 66, 80, 100, 120):
        fraction = absorption(model, REFERENCE_NM, pressure) / ambient
        out[f'a532_fraction_{pressure}'] = fraction
        flag = ('extracted' if pressure <= E1_LAST_SAMPLED_GPA
                else 'INTERPOLATED (E1)')
        print(f'      {pressure:9d}{fraction*100:19.2f}%   {flag}')
    out['e1_last_sampled_GPa'] = E1_LAST_SAMPLED_GPA
    print('      absorption at a fixed line is a sideband process, so it does '
          'not\n      switch off at the ZPL crossing Doherty projected for '
          '~66 GPa.\n      Dai et al. ran 532 nm to 140 GPa; the kernel keeps '
          'A_532 positive.')
    print(f'      But erratum E1 applies above {E1_LAST_SAMPLED_GPA:.0f} GPa: '
          'the ZPL passes 532 nm at 92.1 GPa,\n      and beyond the last '
          'reference curve that brackets the line with real\n      samples '
          '(80 GPa) the value is interpolation across an empty interval.\n'
          '      The 100 and 120 GPa rows must not be quoted as data.')
    return out


if __name__ == '__main__':
    report()
