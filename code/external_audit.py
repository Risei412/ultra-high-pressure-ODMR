"""What the theory explains in the two external sources, and what it does not.

Two documents outside this project carry results the theory can be scored
against.  Neither was written with it in mind.

  HO   K. O. Ho et al., arXiv:2606.02399 (2026).  Supplies the optical kernel.
       The question is not whether we reproduce it -- we take it -- but which
       of Ho's OWN results follow from our layer and which are inputs we
       cannot derive.  `docs/ho_results_audit.md`.

  THE  P. Bhattacharyya, PhD thesis, UC Berkeley 2022
       (`docs/ref/Principle_and_Applications_of_.pdf`).  Independent of Ho and
       used nowhere in the calibration.  `docs/bhattacharyya_thesis_scope.md`
       sorted its chapter 6; this module carries the whole-thesis pass and the
       corrections it forced.

Five results, each returned as a dict for the tests.

X1  POWER CONVENTION (erratum E5).  `thesis_crosscheck` and
    `anvil_transmission` compared bare cross sections, sigma(lambda) -- fixed
    photon flux.  The thesis states fixed optical power ("similar laser and
    microwave powers"), where the absorbed rate is A = lambda sigma, which is
    also the A the frozen record uses everywhere else.  Correcting it moves
    the predicted Fig. 6.3 crossover from 51.82 to 54.37 GPa against an
    observed 54.4 GPa, and the anvil margin from 1.9 to 2.9 sigma.  The
    repository's single external hit was better than it was reported to be.

X2  E1's EMPTY INTERVAL REACHES 100 GPa, not just 120.  Erratum E1 named
    120 GPa as the pressure at which 532 nm falls into the gap between the
    axis-edge baseline sample and the ZPL onset.  Checking every reference
    curve, the last one whose real samples still bracket 532 nm is 80 GPa;
    at 100 GPa the nearest sample below is already the 885.7 nm axis edge.
    The kernel's ZPL passes 532 nm at 92.1 GPa, which is why.  So
    `thesis_crosscheck` T3's 100 GPa row was as much an interpolation
    artefact as its 120 GPa row.

X3  CHAPTER 7 IS OUTSIDE THE HYDROSTATIC KERNEL, and says by how much.  The
    hydride work ran cwODMR at 130-140 GPa on a [111] cut anvil at cryogenic
    temperature, through the confocal setup of Fig. 1.5, i.e. at 532 nm.
    Extrapolating the kernel's ZPL, 532 nm then sits 107-135 meV BELOW it, and
    E1's own thermal-activation factor for sub-ZPL absorption is 1e-7 at 90 K
    and 1e-13 at 50 K.  It plainly worked, so the local optical gap of the
    [111] NV group must be red shifted by at least that much -- which is the
    SAME red shift the thesis invokes in section 6.4 to explain why [111] cut
    culets retain contrast.  Contrast retention and green viability are one
    effect, and our layer adds a third consequence: the absorption band, hence
    the optimum, moves red with it.  The bound gives lambda_opt >= 458 nm at
    120 GPa in that geometry, not 440.65 nm.

X4  THE GEOMETRY AND THE ANVIL ARE DEGENERATE in Fig. 6.3(b).  `nv_model`'s
    C-4 constant says a standard flat culet (alpha = 0.56) carries only 0.564
    of the quasi-hydrostatic ZPL shift.  Applied to the kernel that moves the
    450/532 crossover from 54.4 to 94 GPa and destroys the retrodiction --
    UNLESS a constant anvil factor absorbs it, and it can: g = 0.564 fits
    better than g = 1, at T(450)/T(532) = 3.0.  A type Ia anvil passing blue
    three times better than green is not physical, so requiring T <= 1 is what
    pins g ~ 1.  The margin quoted for 440.65 nm rests on that requirement,
    not on the fit alone.

X5  A SECOND ATTENUATOR, and this one is pressure dependent.  Section 2.4 and
    Fig. 2.7(b) of the thesis record that CsI metallises near 50-65 GPa and
    absorbs the short wavelengths -- "the sample chamber appears red".  The
    excitation path of the Fig. 6.3 loading crosses the chamber.  Adding a
    linear pressure term to the constant-anvil fit improves chi2 by 4.0
    (2.0 sigma) in exactly that direction.  Not significant, and the
    extrapolated ambient value is unphysical, so the linear form is wrong --
    but a pressure-dependent attenuator does not merely rescale the ratio, it
    MOVES the crossover, and the recommended re-run at 65-70 GPa therefore has
    to name a non-metallising medium or it tests the medium, not the NV.

Run for the tables.
"""
import numpy as np
from scipy.optimize import brentq

import anvil_transmission as anvil
import thesis_crosscheck as thesis
from ho_spectrum_model import HoPublishedSpectrumModel
from nv_model import _alpha_factor, eV2nm, nm2eV
from theory_a3_branch_exchange import branch_table, identify_branches

KB_EV = 8.617333262e-5          # Boltzmann constant, eV/K
NV_ZPL_AMBIENT_EV = 1.945       # the NV- zero-phonon line at ambient

# Chapter 7: pressures, and the temperatures its Meissner data were taken at.
CH7_PRESSURES_GPA = (130.0, 140.0)
CH7_TEMPERATURES_K = (300.0, 90.0, 50.0)
GREEN_NM = 532.0
BLUE_NM = 450.0
# Hilberer's two geometries, as carried by nv_model's C-4 block.
ALPHA_FLAT_CULET = 0.56
ALPHA_MICROPILLAR = 0.95


def _zpl_track():
    """Kernel ZPL against pressure, from the raw extracted samples."""
    rows = branch_table(identify_branches())
    pressure = np.array([row['pressure'] for row in rows])
    energy = np.array([row['zpl_ev'] for row in rows])
    return pressure, energy


def zpl_energy(pressure):
    """Kernel ZPL at `pressure`, linearly extrapolated above 120 GPa."""
    grid, energy = _zpl_track()
    slope = (energy[-1] - energy[-2]) / (grid[-1] - grid[-2])
    pressure = float(pressure)
    if pressure <= grid[-1]:
        return float(np.interp(pressure, grid, energy))
    return float(energy[-1] + slope * (pressure - grid[-1]))


# --- X1  the power convention --------------------------------------------

def power_convention(model=None):
    """Fixed photon flux against fixed optical power, on the thesis's figure."""
    model = HoPublishedSpectrumModel() if model is None else model
    sigma = lambda lam, p: float(model.sigma_abs(nm2eV(lam), float(p)))
    merit = lambda lam, p: lam * sigma(lam, p)

    def cross(f):
        return float(brentq(lambda p: f(BLUE_NM, p) / f(GREEN_NM, p) - 1.0,
                            1.0, 119.0))

    observed = anvil.crossover(anvil.observed_ratio())
    return {
        'crossover_fixed_flux_GPa': cross(sigma),
        'crossover_fixed_power_GPa': cross(merit),
        'crossover_observed_GPa': observed,
        'snr_ratio_50GPa_fixed_flux': np.sqrt(sigma(BLUE_NM, 50.0)
                                              / sigma(GREEN_NM, 50.0)),
        'snr_ratio_50GPa_fixed_power': np.sqrt(merit(BLUE_NM, 50.0)
                                               / merit(GREEN_NM, 50.0)),
    }


# --- X2  how far E1's empty interval reaches ------------------------------

def sub_zpl_coverage(lam=GREEN_NM, model=None):
    """Which reference curves still bracket `lam` with real extracted samples.

    Erratum E1 is not a property of 120 GPa; it is a property of the line
    sitting below the ZPL, and it switches on where the ZPL passes the line.
    """
    model = HoPublishedSpectrumModel() if model is None else model
    energy = nm2eV(lam)
    rows = []
    for pressure in model.pressures:
        grid, _ = model.spectra[pressure]
        below = grid[grid < energy]
        above = grid[grid >= energy]
        nearest_below = float(below.max()) if len(below) else None
        nearest_above = float(above.min()) if len(above) else None
        # A sample within 0.05 eV on the low side means the line is bracketed
        # by real data rather than by the figure's axis-edge baseline.
        bracketed = (nearest_below is not None
                     and energy - nearest_below < 0.05)
        rows.append({
            'pressure_GPa': float(pressure),
            'nearest_below_nm': (None if nearest_below is None
                                 else float(eV2nm(nearest_below))),
            'nearest_above_nm': (None if nearest_above is None
                                 else float(eV2nm(nearest_above))),
            'bracketed_by_real_samples': bool(bracketed),
        })
    sampled = [r['pressure_GPa'] for r in rows if r['bracketed_by_real_samples']]
    grid, zpl = _zpl_track()
    return {
        'line_nm': float(lam),
        'rows': rows,
        'last_sampled_GPa': max(sampled) if sampled else None,
        'zpl_crossing_GPa': float(np.interp(energy, zpl, grid)),
    }


# --- X3  chapter 7 against the hydrostatic kernel -------------------------

def green_at_megabar(pressures=CH7_PRESSURES_GPA,
                     temperatures=CH7_TEMPERATURES_K, lam=GREEN_NM):
    """How far below the hydrostatic ZPL chapter 7's excitation line sits.

    Returns the thermal-activation factors E1 prescribes for sub-ZPL
    absorption, and the red shift of the local optical gap that would be
    needed to put the line back into the sideband.
    """
    energy = nm2eV(lam)
    rows = []
    for pressure in pressures:
        deficit = zpl_energy(pressure) - energy
        rows.append({
            'pressure_GPa': float(pressure),
            'zpl_eV': zpl_energy(pressure),
            'zpl_nm': float(eV2nm(zpl_energy(pressure))),
            'deficit_meV': float(deficit * 1e3),
            'hot_band': {float(T): float(np.exp(-deficit / (KB_EV * T)))
                         for T in temperatures},
            'required_red_shift_meV': float(max(deficit, 0.0) * 1e3),
        })
    # The deviatoric stress driving the shift scales with pressure, so carry
    # the requirement back to 120 GPa before turning it into a wavelength.
    worst = max(rows, key=lambda r: r['required_red_shift_meV'])
    shift_120 = (worst['required_red_shift_meV'] * 1e-3
                 * 120.0 / worst['pressure_GPa'])
    band_peak_120 = nm2eV(440.65)
    return {
        'rows': rows,
        'required_red_shift_at_120GPa_meV': float(shift_120 * 1e3),
        'lambda_opt_lower_bound_nm': float(eV2nm(band_peak_120 - shift_120)),
        'frozen_lambda_opt_nm': 440.65,
    }


# --- X4  band position and anvil factor are degenerate --------------------

def _shifted_merit(model, lam, pressure, g):
    """A = lambda sigma with the band scaled to `g` of the hydrostatic shift."""
    shift = (g - 1.0) * (zpl_energy(pressure) - NV_ZPL_AMBIENT_EV)
    return lam * float(model.sigma_abs(nm2eV(lam) - shift, float(pressure)))


def geometry_degeneracy(scan=(1.0, 0.9, 0.8, 0.7), model=None):
    """Fit a constant anvil factor at several band positions `g`.

    `g` = 1 is Ho's quasi-hydrostatic geometry; `g` = 0.564 is what C-4 gives
    for a standard flat culet.  Both fit Fig. 6.3(b) -- the data cannot tell
    them apart without a bound on the anvil.
    """
    model = HoPublishedSpectrumModel() if model is None else model
    rows = anvil.observed_ratio()[1:]           # ambient dropped, as the fit does
    pressure = np.array([r['pressure_GPa'] for r in rows])
    observed = np.log10([r['ratio'] for r in rows])
    weight = 1.0 / np.array([r['sigma_log10'] for r in rows]) ** 2

    def fit(g):
        predicted = np.array([
            0.5 * np.log10(_shifted_merit(model, BLUE_NM, p, g)
                           / _shifted_merit(model, GREEN_NM, p, g))
            for p in pressure])
        residual = observed - predicted
        offset = float(np.sum(weight * residual) / np.sum(weight))
        return (float(np.sum(weight * (residual - offset) ** 2)),
                float(10.0 ** (2.0 * offset)))

    values = list(scan) + [_alpha_factor(ALPHA_FLAT_CULET)]
    table = []
    for g in values:
        chi2, transmission = fit(g)
        table.append({'g': float(g), 'chi2': chi2,
                      'transmission_ratio': transmission})
    return {
        'rows': table,
        'g_flat_culet': float(_alpha_factor(ALPHA_FLAT_CULET)),
        'g_micropillar': float(_alpha_factor(ALPHA_MICROPILLAR)),
        'crossover_at_flat_culet_GPa': float(brentq(
            lambda p: (_shifted_merit(model, BLUE_NM, p,
                                      _alpha_factor(ALPHA_FLAT_CULET))
                       / _shifted_merit(model, GREEN_NM, p,
                                        _alpha_factor(ALPHA_FLAT_CULET))) - 1.0,
            1.0, 119.0)),
    }


# --- X5  is the leftover attenuation pressure dependent? ------------------

def residual_trend():
    """Test a pressure-dependent term against the constant anvil factor."""
    rows = anvil.observed_ratio()[1:]
    pressure = np.array([r['pressure_GPa'] for r in rows])
    residual = np.log10([r['ratio'] for r in rows]) - np.log10(
        [anvil.predicted_ratio(r['pressure_GPa']) for r in rows])
    weight = 1.0 / np.array([r['sigma_log10'] for r in rows]) ** 2

    offset = np.sum(weight * residual) / np.sum(weight)
    chi2_constant = float(np.sum(weight * (residual - offset) ** 2))

    design = np.vstack([np.ones_like(pressure), pressure]).T
    normal = design.T @ np.diag(weight) @ design
    covariance = np.linalg.inv(normal)
    beta = covariance @ (design.T @ (weight * residual))
    chi2_linear = float(np.sum(weight * (residual - design @ beta) ** 2))
    error = np.sqrt(np.diag(covariance))
    return {
        'pressures_GPa': pressure.tolist(),
        'residual_log10': residual.tolist(),
        'chi2_constant': chi2_constant,
        'chi2_linear': chi2_linear,
        'delta_chi2': chi2_constant - chi2_linear,
        'slope_per_GPa': float(beta[1]),
        'slope_sigma': float(abs(beta[1] / error[1])),
        'ambient_extrapolation': float(10.0 ** (2.0 * beta[0])),
    }


def report():
    out = {}

    print('X1  power convention (erratum E5)')
    out['power_convention'] = x1 = power_convention()
    print(f'      fixed photon flux : crossover '
          f'{x1["crossover_fixed_flux_GPa"]:.2f} GPa, '
          f'SNR(50 GPa) {x1["snr_ratio_50GPa_fixed_flux"]:.3f}')
    print(f'      fixed optical power: crossover '
          f'{x1["crossover_fixed_power_GPa"]:.2f} GPa, '
          f'SNR(50 GPa) {x1["snr_ratio_50GPa_fixed_power"]:.3f}')
    print(f'      observed (Fig. 6.3(b)) : '
          f'{x1["crossover_observed_GPa"]:.2f} GPa')
    print('      the thesis states fixed power, so the second row is the '
          'comparable one')

    print('\nX2  how far erratum E1 reaches at 532 nm')
    out['sub_zpl'] = x2 = sub_zpl_coverage()
    print(f'      {"P [GPa]":>9}{"nearest sample below":>23}   provenance')
    for row in x2['rows']:
        below = ('none' if row['nearest_below_nm'] is None
                 else f'{row["nearest_below_nm"]:.1f} nm')
        flag = ('extracted' if row['bracketed_by_real_samples']
                else 'axis-edge baseline -> E1')
        print(f'      {row["pressure_GPa"]:9.0f}{below:>23}   {flag}')
    print(f'      ZPL passes 532 nm at {x2["zpl_crossing_GPa"]:.1f} GPa; '
          f'last usable curve {x2["last_sampled_GPa"]:.0f} GPa')

    print('\nX3  chapter 7: 532 nm at 130-140 GPa, cryogenic, [111] culet')
    out['green_at_megabar'] = x3 = green_at_megabar()
    for row in x3['rows']:
        factors = '  '.join(f'{T:.0f} K: {v:.1e}'
                            for T, v in sorted(row['hot_band'].items(),
                                               reverse=True))
        print(f'      {row["pressure_GPa"]:5.0f} GPa  ZPL {row["zpl_nm"]:.1f} nm'
              f'  532 nm is {row["deficit_meV"]:.0f} meV below it   {factors}')
    print(f'      so the [111] group\'s gap must be red shifted by >= '
          f'{x3["required_red_shift_at_120GPa_meV"]:.0f} meV at 120 GPa,')
    print(f'      which puts the optimum at >= '
          f'{x3["lambda_opt_lower_bound_nm"]:.1f} nm, not '
          f'{x3["frozen_lambda_opt_nm"]:.2f} nm')

    print('\nX4  band position against anvil factor, on Fig. 6.3(b)')
    out['geometry'] = x4 = geometry_degeneracy()
    print(f'      {"g":>7}{"chi2 (dof 3)":>14}{"T(450)/T(532)":>16}')
    for row in x4['rows']:
        print(f'      {row["g"]:7.3f}{row["chi2"]:14.2f}'
              f'{row["transmission_ratio"]:16.2f}')
    print(f'      C-4 puts a flat culet at g = {x4["g_flat_culet"]:.3f}, which '
          f'alone moves the crossover to\n'
          f'      {x4["crossover_at_flat_culet_GPa"]:.0f} GPa; it survives '
          f'only by demanding an anvil that passes blue better\n'
          f'      than green.  T <= 1 is what pins g ~ 1.')

    print('\nX5  is the leftover attenuation pressure dependent?')
    out['residual_trend'] = x5 = residual_trend()
    print(f'      constant  chi2 = {x5["chi2_constant"]:.2f} (dof 4)')
    print(f'      + linear  chi2 = {x5["chi2_linear"]:.2f} (dof 3), '
          f'slope {x5["slope_sigma"]:.1f} sigma, '
          f'delta chi2 = {x5["delta_chi2"]:.2f}')
    print(f'      extrapolated to ambient the linear form gives '
          f'T(450)/T(532) = {x5["ambient_extrapolation"]:.1f} > 1, so the '
          f'form is wrong;\n      what the data say is that the residual '
          f'DECLINES with pressure, which is what a\n      metallising '
          f'pressure medium does (thesis section 2.4, CsI at ~50-65 GPa).')
    return out


if __name__ == '__main__':
    report()
