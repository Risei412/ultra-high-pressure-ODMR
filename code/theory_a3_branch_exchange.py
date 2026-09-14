"""Addendum A3: pressure-driven branch exchange of the sensitivity optimum.

A2 showed that at one pressure the sensitivity-optimal set is a level set of
the absorption kernel, and that raising the power walks that level down through
the kernel's critical values.  A3 asks the orthogonal question: what does
*pressure* do?

The answer is that pressure does something A2's ladder cannot.  An absorption
spectrum built on a Franck-Condon progression carries two structurally distinct
branches -- the zero-phonon line and the phonon sideband -- whose peak heights
scale differently with the Huang-Rhys factor S.  Pressure raises S, and the
branches therefore move relative to one another.  Where their weighted peaks
cross, the global optimum switches branch discontinuously.

Two degeneracy mechanisms, then, with different causes and different signatures:

    A2 ladder          driven by power     needs I > I_c     within one branch
    A3 branch exchange driven by pressure  needs P = P*      between branches

The A3 degeneracy exists at *zero* power, which is what makes it distinguishable
in practice: at P* the low-power optimum is already twofold, before any of A2's
power-induced structure appears.

Theorem X (branch exchange)
    Let branches i, j have peak values A_i(P), A_j(P) of the fixed-power figure
    of merit A = lambda sigma_abs.  The global optimum switches branch at any
    P* where D(P) = ln(A_j/A_i) changes sign, and the optimal wavelength jumps
    by lambda_i(P*) - lambda_j(P*) with no intermediate values.  For a
    Franck-Condon pair,

        A_ZPL  ~  e^{-S}  / Gamma_ZPL      (weight collapses as S grows)
        A_SB   ~ (1-e^{-S}) / Gamma_SB     (weight saturates at 1)

    so dD/dP = dS/dP + d/dP ln(Gamma_ZPL/Gamma_SB) + ... > 0 whenever pressure
    strengthens electron-phonon coupling and broadens the ZPL faster than the
    sideband.  D is then monotone and the crossing is unique: a generic,
    once-only exchange rather than an accident of one material.

The worked example is Ho et al.'s published 120 GPa series, in which both
branches are present as separate *raw extracted* maxima at all seven pressures,
so the exchange is read off the data rather than out of an interpolation.

.. warning::

   **Erratum E3 supersedes every ZPL-derived number in this module.**  Checking
   the extraction against the source figure showed that the panel (e)
   zero-phonon-line spikes are clipped by the axis: all seven rise to within
   1 unit of the axis top.  The ZPL peak heights in the CSV, and hence
   ``exchange_pressure``'s P* = 87.9 GPa and the "x6.04 ZPL collapse", are
   artefacts of that clipping and are withdrawn.

   The sideband half of this module is exact --- it reproduces a direct trace
   of the figure to better than 1 %.  Theorem X also survives: taking the ZPL
   weight from the published Debye-Waller factor instead, the branch ratio is
   still monotone, so the crossing is still unique.  What does not survive is a
   single value for P*, because comparing a narrow line with a broad band is
   bandwidth dependent.  See ``figure_validation.py`` for the corrected
   treatment and ``docs/theory_a3_branch_exchange.md`` for the write-up.
"""
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from ho_spectrum_model import DEFAULT_DATA, HBARC

# Branches are separated by the Franck-Condon displacement S*hbar_omega, which
# is ~0.23-0.40 eV in this system; 0.20 eV cleanly separates them without
# catching sideband substructure.
BRANCH_GAP_EV = 0.20
# Above this the ZPL has already left the sampled window at low pressure.
ZPL_SEARCH_MAX_EV = 2.50


@dataclass(frozen=True)
class BranchPoint:
    """One branch's peak at one pressure, from the raw extracted samples."""

    pressure: float
    energy_ev: float
    sigma: float

    @property
    def wavelength_nm(self):
        return HBARC / self.energy_ev

    @property
    def figure_of_merit(self):
        """A = lambda * sigma, the fixed-incident-power absorbed-photon proxy."""
        return self.wavelength_nm * self.sigma


def load_raw_spectra(path=DEFAULT_DATA):
    """The extracted Fig. 1(e) samples, grouped by pressure.  No interpolation."""
    grouped = {}
    with open(path) as stream:
        for line in stream:
            text = line.strip()
            if not text or text.startswith('#') or text.startswith('pressure_'):
                continue
            pressure, energy, absorption = (float(v) for v in text.split(','))
            grouped.setdefault(pressure, []).append((energy, absorption))
    return {p: np.array(sorted(v), float) for p, v in grouped.items()}


def raw_local_maxima(values):
    """Local maxima of the extracted samples themselves."""
    energy, sigma = values[:, 0], values[:, 1]
    return [(float(energy[i]), float(sigma[i]))
            for i in range(1, len(energy) - 1)
            if sigma[i] >= sigma[i - 1] and sigma[i] > sigma[i + 1]]


def identify_branches(path=DEFAULT_DATA, gap_ev=BRANCH_GAP_EV):
    """Track the ZPL and sideband branches across pressure.

    The ZPL is the strongest raw maximum below `ZPL_SEARCH_MAX_EV`; the
    sideband is the strongest raw maximum at least `gap_ev` above it.  Both are
    genuine extracted maxima at every pressure, so neither branch is an artefact
    of the pressure interpolation.
    """
    spectra = load_raw_spectra(path)
    zpl, sideband = [], []
    for pressure in sorted(spectra):
        maxima = raw_local_maxima(spectra[pressure])
        low = [m for m in maxima if m[0] < ZPL_SEARCH_MAX_EV]
        if not low:
            raise ValueError(f'no ZPL candidate at {pressure:g} GPa')
        zpl_point = max(low, key=lambda m: m[1])
        high = [m for m in maxima if m[0] > zpl_point[0] + gap_ev]
        if not high:
            raise ValueError(f'no sideband candidate at {pressure:g} GPa')
        sideband_point = max(high, key=lambda m: m[1])
        zpl.append(BranchPoint(pressure, *zpl_point))
        sideband.append(BranchPoint(pressure, *sideband_point))
    return {'zpl': zpl, 'sideband': sideband}


def branch_table(branches=None):
    branches = branches or identify_branches()
    rows = []
    for z, s in zip(branches['zpl'], branches['sideband']):
        rows.append({
            'pressure': z.pressure,
            'zpl_nm': z.wavelength_nm, 'zpl_ev': z.energy_ev, 'zpl_sigma': z.sigma,
            'sb_nm': s.wavelength_nm, 'sb_ev': s.energy_ev, 'sb_sigma': s.sigma,
            'gap_ev': s.energy_ev - z.energy_ev,
            'sigma_ratio': s.sigma / z.sigma,
            'A_ratio': s.figure_of_merit / z.figure_of_merit,
        })
    return rows


def _spline(pressures, values):
    return CubicSpline(np.asarray(pressures, float), np.asarray(values, float))


def exchange_pressure(branches=None, weighted=True):
    """Solve ln(A_SB/A_ZPL) = 0 for the branch-exchange pressure P*.

    SUPERSEDED by erratum E3: A_ZPL here comes from the clipped panel (e)
    spikes, so the returned P* (87.9 GPa at fixed optical power) is not
    defensible.  Retained because it is a faithful computation of what the CSV
    contains, and because the tests pin it as a regression.  Use
    ``figure_validation.exchange_pressure_at_bandwidth`` instead.

    `weighted=True` uses A = lambda sigma (fixed incident optical power);
    `weighted=False` uses sigma alone (fixed incident photon flux).
    """
    branches = branches or identify_branches()
    pressures = [b.pressure for b in branches['zpl']]
    if weighted:
        ratio = [s.figure_of_merit / z.figure_of_merit
                 for z, s in zip(branches['zpl'], branches['sideband'])]
    else:
        ratio = [s.sigma / z.sigma
                 for z, s in zip(branches['zpl'], branches['sideband'])]
    curve = _spline(pressures, np.log(ratio))
    if curve(pressures[0]) * curve(pressures[-1]) > 0.0:
        return float('nan')
    star = float(brentq(curve, pressures[0], pressures[-1]))
    zpl_nm = float(_spline(pressures, [b.wavelength_nm for b in branches['zpl']])(star))
    sb_nm = float(_spline(pressures, [b.wavelength_nm for b in branches['sideband']])(star))
    return {'pressure': star, 'zpl_nm': zpl_nm, 'sideband_nm': sb_nm,
            'jump_nm': zpl_nm - sb_nm, 'weighted': weighted}


def monotonicity(branches=None):
    """Theorem X's antecedent: is ln(A_SB/A_ZPL) monotone in pressure?"""
    branches = branches or identify_branches()
    pressures = np.array([b.pressure for b in branches['zpl']])
    ratio = np.array([s.figure_of_merit / z.figure_of_merit
                      for z, s in zip(branches['zpl'], branches['sideband'])])
    log_ratio = np.log(ratio)
    steps = np.diff(log_ratio)
    return {
        'pressures': pressures, 'log_ratio': log_ratio,
        'increments': steps,
        'monotone_increasing': bool(np.all(steps > 0.0)),
        'total_change': float(log_ratio[-1] - log_ratio[0]),
        'n_sign_changes': int(np.sum(np.diff(np.sign(log_ratio)) != 0)),
    }


def coupling_growth(branches=None):
    """The physical driver: the Franck-Condon displacement S*hbar_omega."""
    branches = branches or identify_branches()
    pressures = np.array([b.pressure for b in branches['zpl']])
    gap = np.array([s.energy_ev - z.energy_ev
                    for z, s in zip(branches['zpl'], branches['sideband'])])
    slope = np.polyfit(pressures, gap, 1)[0]
    return {'pressures': pressures, 'gap_ev': gap,
            'growth_meV_per_GPa': float(slope * 1e3),
            'relative_growth': float(gap[-1] / gap[0] - 1.0),
            'monotone': bool(np.all(np.diff(gap) > 0.0))}


def branch_shift_rates(branches=None):
    """Each branch has its own pressure coefficient -- why they can cross."""
    branches = branches or identify_branches()
    out = {}
    for name in ('zpl', 'sideband'):
        points = branches[name]
        pressures = np.array([b.pressure for b in points])
        energies = np.array([b.energy_ev for b in points])
        out[name] = {
            'mean_meV_per_GPa': float(np.polyfit(pressures, energies, 1)[0] * 1e3),
            'local_meV_per_GPa': (np.diff(energies) / np.diff(pressures) * 1e3),
            'sigma_falls_by': float(points[0].sigma / points[-1].sigma),
        }
    return out


def fixed_wavelength_crossover(lam_a=457.0, lam_b=532.0, path=DEFAULT_DATA):
    """The 457/532 comparison, to show it is NOT the branch exchange.

    A crossover between two fixed probe wavelengths happens whenever a single
    peak sweeps past their midpoint; it needs no branch structure at all.
    """
    from ho_spectrum_model import HoPublishedSpectrumModel
    model = HoPublishedSpectrumModel(path)

    def merit(lam, pressure):
        return lam * float(model.sigma_abs(HBARC / lam, pressure))

    def gap(pressure):
        return np.log(merit(lam_a, pressure) / merit(lam_b, pressure))

    if gap(20.0) * gap(120.0) > 0.0:
        return float('nan')
    return float(brentq(gap, 20.0, 120.0))


def zero_power_degeneracy(branches=None):
    """At P* the low-power optimum is already twofold, with no power applied.

    The *structure* stands -- an exchange produces a degenerate pair at zero
    power -- but the pressure and separation reported here inherit E3's clipped
    P*, so treat them as illustrative only.

    This is what separates A3 from A2: the ladder needs I > I_c, whereas the
    branch exchange produces an exactly degenerate pair in the low-power limit.
    """
    star = exchange_pressure(branches)
    if not isinstance(star, dict):
        return None
    return {'pressure': star['pressure'],
            'members_nm': (star['sideband_nm'], star['zpl_nm']),
            'separation_nm': star['jump_nm'],
            'requires_power': False}


def general_condition():
    """Theorem X's antecedent, stated so it can be checked in any system.

    d/dP ln(A_SB/A_ZPL) > 0 holds whenever, with growing pressure,
      (i)  the Huang-Rhys factor S increases (weight moves ZPL -> sideband), and
      (ii) the ZPL broadens at least as fast as the sideband.
    Both are generic for a colour centre under compression, so the exchange is
    expected rather than exceptional; what varies between materials is only
    where P* falls relative to the accessible pressure range.
    """
    return {
        'driver_1': 'dS/dP > 0  (Franck-Condon weight transfer ZPL -> sideband)',
        'driver_2': 'd/dP ln(Gamma_ZPL/Gamma_SB) >= 0  (ZPL broadens faster)',
        'consequence': 'ln(A_SB/A_ZPL) monotone increasing -> unique crossing',
        'observable': 'optimal wavelength jumps discontinuously at P*',
    }


def main():
    bar = '=' * 78
    print(bar)
    print('Addendum A3 - pressure-driven branch exchange of the optimum')
    print(bar)

    branches = identify_branches()

    print('\n[X1] Branches tracked in the RAW extracted samples')
    print('    P    ZPL [nm]  sigma     SB [nm]  sigma     gap [eV]  A_SB/A_ZPL')
    for row in branch_table(branches):
        print(f"  {row['pressure']:4.0f}   {row['zpl_nm']:8.1f} {row['zpl_sigma']:7.3f}   "
              f"{row['sb_nm']:8.1f} {row['sb_sigma']:7.3f}   {row['gap_ev']:7.4f}   "
              f"{row['A_ratio']:8.3f}")

    mono = monotonicity(branches)
    print('\n[X2] Theorem X antecedent: is the log-ratio monotone?')
    print(f"  monotone increasing : {mono['monotone_increasing']}")
    print(f"  total change        : {mono['total_change']:+.3f} in ln, "
          f"i.e. x{np.exp(mono['total_change']):.2f}")
    print(f"  sign changes of ln  : {mono['n_sign_changes']} (1 = a unique exchange)")

    star = exchange_pressure(branches)
    unweighted = exchange_pressure(branches, weighted=False)
    print('\n[X3] Branch exchange pressure')
    print(f"  fixed optical power (A = lambda sigma): P* = {star['pressure']:.1f} GPa")
    print(f"    ZPL branch at {star['zpl_nm']:.1f} nm, sideband at "
          f"{star['sideband_nm']:.1f} nm")
    print(f"    -> the optimum jumps {star['jump_nm']:.1f} nm discontinuously")
    print(f"  fixed photon flux   (sigma alone)     : P* = "
          f"{unweighted['pressure']:.1f} GPa")

    growth = coupling_growth(branches)
    print('\n[X4] Physical driver: Franck-Condon displacement S*hbar_omega')
    print(f"  gap grows {growth['gap_ev'][0]:.4f} -> {growth['gap_ev'][-1]:.4f} eV "
          f"({growth['relative_growth'] * 100:+.0f} %), monotone: {growth['monotone']}")
    print(f"  growth rate = {growth['growth_meV_per_GPa']:.2f} meV/GPa")

    rates = branch_shift_rates(branches)
    print('\n[X5] The two branches are not the same object')
    for name in ('zpl', 'sideband'):
        row = rates[name]
        print(f"  {name:9s}: mean shift {row['mean_meV_per_GPa']:.2f} meV/GPa, "
              f"sigma falls by x{row['sigma_falls_by']:.2f} over 0-120 GPa")

    print('\n[X6] A3 is not the 457/532 crossover')
    print(f"  457/532 fixed-wavelength crossover : "
          f"{fixed_wavelength_crossover():.1f} GPa")
    print(f"  branch exchange                    : {star['pressure']:.1f} GPa")
    print('  -> different phenomena, tens of GPa apart in the same kernel')

    degenerate = zero_power_degeneracy(branches)
    print('\n[X7] The A3 degeneracy needs no power at all')
    print(f"  at P* = {degenerate['pressure']:.1f} GPa the low-power optimum is "
          f"twofold: {degenerate['members_nm'][0]:.1f} and "
          f"{degenerate['members_nm'][1]:.1f} nm")
    print(f"  separation {degenerate['separation_nm']:.1f} nm, "
          f"requires_power = {degenerate['requires_power']}")

    print('\n[X8] General condition for any high-pressure optical sensor')
    for key, text in general_condition().items():
        print(f"  {key:12s}: {text}")


if __name__ == '__main__':
    main()
