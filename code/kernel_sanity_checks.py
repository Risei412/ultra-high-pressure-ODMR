"""Internal sanity checks on the reconstructed Ho kernel.

`repro_yield.py` validates the extraction against Ho's Fig. 5(b), but that
figure probes only two photon energies, 2.331 eV (532 nm) and 2.713 eV
(457 nm), and both sit 2.2 to 8.1 effective phonons above the zero-phonon
line throughout their audit windows.  The cross-figure check therefore says
nothing about the ZPL region -- which is exactly where erratum E3 found the
published spectra to be clipped.  These four checks close that gap using only
material internal to the same paper.

K1  Area normalisation and the blue truncation.
    Ho's spectra are drawn to unit area; the digitised curves integrate to
    1.005, 1.004, 1.002 and 0.997 at 0-60 GPa, confirming it.  Above that the
    integral falls to 0.914 at 120 GPa, and the reason is visible in the data:
    the absorption at the top of the plotted window is still 72.7% of its
    peak at 120 GPa, against 0.0% at ambient.  The band does not return to
    zero inside the figure.  So the blue edge of the analysis window is a
    property of the figure, not of the defect, and area-based quantities are
    underestimated at high pressure by the missing tail.

K2  The ZPL width is resolution limited, not physical.
    The apparent full width at half maximum of the zero-phonon spike is
    2.53-2.59 meV at every one of the seven pressures -- constant to 1% over
    120 GPa.  A real ZPL broadens under compression; that is the second
    driver of Theorem X.  A width that does not move is the plotting
    resolution.  Taken with E3 (all seven apexes clipped at the axis top),
    the only trustworthy ZPL information in Fig. 1(e) is its POSITION.

K3  ZPL area against the published Debye-Waller factor.
    Integrating the background-subtracted spike gives the kernel's own DWF.
    It matches the published panel (c) to within 4% up to 40 GPa -- an
    independent validation of the extraction in the region Fig. 5(b) cannot
    reach -- and then diverges monotonically, over-weighting the ZPL by 1.43x
    at 120 GPa.  Hence: trust the kernel's ZPL below ~40 GPa, use the
    published DWF above it, which is what E3 already prescribed.

K4  Panels (b) and (c) against each other.
    DWF is not exp(-S_abs); the ratio grows smoothly from 2.37 to 4.66.  That
    is expected -- Ho defines S_abs over the Jahn-Teller-INACTIVE modes only
    -- and the residual is the JT-active coupling, S_JT = 0.86 -> 1.54, so
    S_total = 3.89 -> 6.09.  Two panels digitised separately yielding a
    smooth monotone derived quantity is itself evidence that neither
    extraction is broken.

Verdict: the extraction is sound where it is used.  The withdrawn Addendum A3
numbers cannot be recovered from Fig. 1(e) -- K2 shows the width is as
unusable as E3 showed the height to be -- so the published S_abs and DWF
remain the only source for the drivers of Theorem X.

Run for the tables; each check returns its numbers for the tests.
"""
import numpy as np

from ho_spectrum_model import HoPublishedSpectrumModel
from theory_a3_branch_exchange import identify_branches

PANELS_BC = 'data/ho_fig1_panels_bc.csv'

# Half-window for the ZPL spike, ~8x its apparent width, and the fraction of
# the plotted energy range counted as "the top of the window" in K1.
SPIKE_HALF_WINDOW_EV = 0.020
TOP_OF_WINDOW_EV = 3.0
# Pressure below which the kernel's own ZPL is usable (K3).
ZPL_TRUSTED_MAX_GPA = 40.0


def _panels_bc(path=PANELS_BC):
    rows = []
    with open(path) as stream:
        for line in stream:
            if not line[:1].isdigit():
                continue
            pressure, s_abs, dwf = line.split(',')
            rows.append((float(pressure), float(s_abs), float(dwf)))
    rows.sort()
    return (np.array([row[0] for row in rows]),
            np.array([row[1] for row in rows]),
            np.array([row[2] for row in rows]))


def _zpl_apex(model, pressure, energy_ev):
    energy, absorption = model.spectra[pressure]
    index = int(np.argmin(np.abs(energy - energy_ev)))
    return energy, absorption, index


def normalisation(model=None):
    """K1: integrated area and the weight left at the top of the window."""
    model = HoPublishedSpectrumModel() if model is None else model
    out = {}
    for pressure in model.pressures:
        energy, absorption = model.spectra[pressure]
        top = absorption[energy >= TOP_OF_WINDOW_EV]
        out[pressure] = {
            'area': float(np.trapezoid(absorption, energy)),
            'top_fraction': float(top.max() / absorption.max()),
        }
    return out


def zpl_width(model=None):
    """K2: apparent FWHM of the zero-phonon spike, in eV."""
    model = HoPublishedSpectrumModel() if model is None else model
    branches = identify_branches()
    out = {}
    for point in branches['zpl']:
        energy, absorption, index = _zpl_apex(model, point.pressure,
                                              point.energy_ev)
        half = absorption[index] / 2.0
        left = index
        while left > 0 and absorption[left] > half:
            left -= 1
        blue_edge = np.interp(half, [absorption[left], absorption[left + 1]],
                              [energy[left], energy[left + 1]])
        right = index
        while right < len(absorption) - 1 and absorption[right] > half:
            right += 1
        red_edge = np.interp(half, [absorption[right], absorption[right - 1]],
                             [energy[right], energy[right - 1]])
        out[point.pressure] = float(red_edge - blue_edge)
    return out


def zpl_area(model=None):
    """K3: background-subtracted spike area against the published DWF."""
    model = HoPublishedSpectrumModel() if model is None else model
    branches = identify_branches()
    pressures, _, published = _panels_bc()
    published = dict(zip(pressures, published))
    out = {}
    for point in branches['zpl']:
        energy, absorption, _ = _zpl_apex(model, point.pressure,
                                          point.energy_ev)
        window = ((energy >= point.energy_ev - SPIKE_HALF_WINDOW_EV)
                  & (energy <= point.energy_ev + SPIKE_HALF_WINDOW_EV))
        e_win, a_win = energy[window], absorption[window]
        background = np.interp(e_win, [e_win[0], e_win[-1]],
                               [a_win[0], a_win[-1]])
        area = float(np.trapezoid(np.clip(a_win - background, 0.0, None),
                                  e_win))
        out[point.pressure] = {
            'kernel_dwf': area,
            'published_dwf': published[point.pressure],
            'ratio': area / published[point.pressure],
        }
    return out


def jahn_teller_coupling():
    """K4: the coupling left over when DWF is compared with exp(-S_abs)."""
    pressures, s_abs, dwf = _panels_bc()
    s_jt = np.log(np.exp(-s_abs) / dwf)
    return {float(p): {'S_abs': float(a), 'DWF': float(d), 'S_JT': float(j),
                       'S_total': float(a + j)}
            for p, a, d, j in zip(pressures, s_abs, dwf, s_jt)}


def report():
    """Print all four checks; returns them so the tests can assert."""
    model = HoPublishedSpectrumModel()
    norm, width = normalisation(model), zpl_width(model)
    area, jt = zpl_area(model), jahn_teller_coupling()

    print('K1  area normalisation and the blue truncation')
    print(f'    {"P [GPa]":>9}{"integral":>11}{"a(>3.0 eV)/peak":>19}')
    for pressure in sorted(norm):
        row = norm[pressure]
        print(f'    {pressure:9.0f}{row["area"]:11.4f}'
              f'{row["top_fraction"]*100:18.1f}%')

    print('\nK2  ZPL apparent width -- constant, therefore resolution limited')
    print(f'    {"P [GPa]":>9}{"FWHM [meV]":>13}')
    for pressure in sorted(width):
        print(f'    {pressure:9.0f}{width[pressure]*1000:13.2f}')
    values = np.array([width[p] for p in sorted(width)])
    print(f'    spread {(values.max()-values.min())/values.mean()*100:.1f}% '
          'over 120 GPa')

    print('\nK3  kernel ZPL area against the published Debye-Waller factor')
    print(f'    {"P [GPa]":>9}{"kernel":>11}{"published":>12}{"ratio":>9}')
    for pressure in sorted(area):
        row = area[pressure]
        print(f'    {pressure:9.0f}{row["kernel_dwf"]:11.5f}'
              f'{row["published_dwf"]:12.5f}{row["ratio"]:9.2f}')
    print(f'    usable to {ZPL_TRUSTED_MAX_GPA:.0f} GPa; above it use the '
          'published DWF (erratum E3)')

    print('\nK4  panels (b) and (c): the Jahn-Teller-active remainder')
    print(f'    {"P [GPa]":>9}{"S_abs":>9}{"DWF":>10}{"S_JT":>8}{"S_total":>10}')
    for pressure in sorted(jt):
        row = jt[pressure]
        print(f'    {pressure:9.0f}{row["S_abs"]:9.3f}{row["DWF"]:10.5f}'
              f'{row["S_JT"]:8.3f}{row["S_total"]:10.3f}')

    return {'normalisation': norm, 'zpl_width': width, 'zpl_area': area,
            'jahn_teller': jt}


if __name__ == '__main__':
    report()
