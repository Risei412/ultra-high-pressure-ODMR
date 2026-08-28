"""Which geometry does 440.65 nm belong to?

The frozen optical-limit optimum was computed from Ho et al.'s kernel, and Ho's
NV centres sit in a FIB-machined micropillar with stress anisotropy
alpha ~ 0.95 -- quasi-hydrostatic by construction.  Megabar ODMR is not done
there.  It is done on a [111] cut flat culet, where the whole point of the
geometry (Bhattacharyya, section 6.4) is that the 3E orbital of the [111] NV
group RED shifts, closing the intersystem-crossing gap Delta and retaining
contrast.  A red shifted 3E is a red shifted absorption band, and a red shifted
band has a redder optimum.  `external_audit.py` X3 showed the effect is not
small: chapter 7 ran 532 nm at 130-140 GPa in that geometry, which the
hydrostatic kernel says is 107-135 meV BELOW the zero-phonon line.

So the optimum is a function of geometry, and the frozen number is one value of
it.  This module makes that dependence explicit.

## The parameter

One scalar per (culet cut, NV group):

    g = (band shift of this group at chamber pressure P)
        / (quasi-hydrostatic band shift at the same P)

g = 1 is Ho's geometry, by definition.  g < 1 is a group whose 3E red shifts
relative to hydrostatic; g > 1 is one that blue shifts further.  The band is
carried rigidly, so the shift in energy is

    delta(P, g) = (g - 1) * [E_ZPL^hydro(P) - 1.945 eV].

## Why this is NOT an extension of nv_model's C-4

C-4 already carries a stress-anisotropy factor, from Hilberer et al.'s two
measured ZPL shift rates (-769 micropillar against -434 standard anvil, per
unit compressed volume), giving g = 0.564 for a standard flat culet.  The
obvious move was to extend it per NV group.  **That would have been wrong.**

C-4 is a single scalar applied to all four groups.  The two independent
anchors below need OPPOSITE signs for two different groups in two different
cuts, which no single scalar can supply -- and applying C-4's 0.564 to
Bhattacharyya's [100] culet destroys the Fig. 6.3(b) retrodiction unless the
anvil passes blue three times better than green (`external_audit` X4).  C-4 is
a group-averaged, hydrostatic-equivalent quantity.  It is left frozen and
unextended; `test_freeze.py` still pins it.

## Three estimates of g, and what each is worth

  A  QUASI-HYDROSTATIC (micropillar; microdiamonds loose in the chamber).
     g = 1 by definition.  This is the geometry of the frozen 440.65 nm and of
     Ho's own Fig. 5(b) validation.

  B  [100] FLAT CULET, all four groups.  Measured here, from Fig. 6.3(b) plus
     the requirement that a type Ia anvil cannot transmit blue better than
     green.  Comes out at g >~ 1.05, i.e. hydrostatic or slightly bluer --
     which is the sign Fig. 6.11(a) predicts, since [100] stress blue shifts
     all four groups.

  C  [111] FLAT CULET, [111] group.  Bounded here, from chapter 7 having
     worked at all: 532 nm must lie in the sideband at 140 GPa, so
     g <= 0.741.  This is an UPPER bound on g and hence a LOWER bound on the
     optimum; nothing in the available data bounds it from the other side.

B and C disagree in sign, in the direction Davies & Hamer and the thesis's
Fig. 6.11 require.  That is the evidence that g is group-resolved and C-4 is
not.

## The answer

At 120 GPa, carrying the band rigidly:

    geometry                          g          lambda_opt
    quasi-hydrostatic (Ho)          1.000        440.65 nm   <- the frozen value
    [100] flat culet               >=1.048       <=437.20 nm
    [111] flat culet, [111] group  <=0.741       >=460.33 nm

A second shift model -- mapping to the effective hydrostatic pressure that
produces the same ZPL shift, so that the Huang-Rhys factor moves too -- puts
the [111] bound at 473.07 nm instead.  The rigid shift is the conservative one,
so the bound is quoted from it.

The [100] row moves the wrong way by 3.4 nm, which is inside the 5% tolerance
band and not worth quoting on its own; it is here because its SIGN is the
check.  A single scalar cannot give 437 nm and 460 nm at the same pressure, so
the two rows together are what rules out a group-averaged geometry factor.

**457 nm, the line Ho switched to above 55 GPa, is inside the [111] band.**
The frozen statement that 457 nm costs 4.5% is a quasi-hydrostatic statement.

Run for the tables.
"""
import numpy as np
from scipy.optimize import brentq

import anvil_transmission as anvil
from external_audit import NV_ZPL_AMBIENT_EV, zpl_energy
from ho_spectrum_model import HoPublishedSpectrumModel
from nv_model import _alpha_factor, eV2nm, nm2eV
from theory_a1_generalization import DATA_WINDOW, Kernel
from theory_a2_multiplicity import critical_values

FROZEN_PRESSURE_GPA = 120.0
FROZEN_OPTIMUM_NM = 440.65
GREEN_NM = 532.0
BLUE_NM = 450.0
# Chapter 7's pressures.  The binding one is the highest: 532 nm has to be in
# the sideband there too.
CH7_PRESSURES_GPA = (130.0, 140.0)
# Hilberer's standard flat culet, as carried by nv_model's C-4 block.  Kept for
# comparison only -- see the module docstring on why it is not extended.
ALPHA_FLAT_CULET = 0.56
# Above g ~ 1.25 a blue-shifted band carries 532 nm out of the extracted
# window at the lower pressures of Fig. 6.3(b) and the fit is undefined.
G_SCAN_LIMITS = (1.0, 1.20)


def band_shift_eV(pressure, g):
    """Rigid energy shift of the absorption band for geometry `g`.

    Negative for a group whose 3E red shifts relative to hydrostatic.
    """
    return (float(g) - 1.0) * (zpl_energy(pressure) - NV_ZPL_AMBIENT_EV)


class GeometryKernel(Kernel):
    """`Kernel` with the whole band carried rigidly by `band_shift_eV`.

    The data window travels with the band: the figure supplies absorption over
    a fixed ENERGY interval, so shifting the band shifts the wavelengths at
    which that interval can be sampled.
    """

    def __init__(self, g=1.0, pressure=FROZEN_PRESSURE_GPA, window=DATA_WINDOW):
        self.g = float(g)
        self.shift = band_shift_eV(pressure, g)
        lo, hi = (float(x) for x in window)
        # window in energy, carried by the shift, then back to wavelength
        e_hi, e_lo = nm2eV(lo) + self.shift, nm2eV(hi) + self.shift
        super().__init__(pressure=pressure,
                         window=(float(eV2nm(e_hi)), float(eV2nm(e_lo))))

    def _raw(self, lam):
        lam = np.asarray(lam, float)
        return lam * self.optical.sigma_abs(nm2eV(lam) - self.shift,
                                            self.pressure)


def optimum(g=1.0, pressure=FROZEN_PRESSURE_GPA, step=0.01):
    """Optical-limit optimum at fixed incident optical power, for geometry `g`."""
    kernel = GeometryKernel(g, pressure)
    grid = np.arange(kernel.lam_min, kernel.lam_max + step / 2.0, step)
    return float(grid[int(np.argmax(kernel.a(grid)))])


def optimum_by_compression(g=1.0, pressure=FROZEN_PRESSURE_GPA, step=0.01):
    """The same question under the other shift model, as a robustness check.

    Instead of carrying the band rigidly, map to the hydrostatic pressure that
    produces the same ZPL shift.  The Huang-Rhys factor then moves too, so the
    Franck-Condon displacement shrinks along with the ZPL.  This gives a REDDER
    optimum than the rigid shift, so quoting bounds from the rigid shift is
    conservative.
    """
    target = NV_ZPL_AMBIENT_EV + float(g) * (zpl_energy(pressure)
                                             - NV_ZPL_AMBIENT_EV)
    effective = brentq(lambda p: zpl_energy(p) - target, 0.0, 120.0)
    kernel = Kernel(pressure=effective)
    grid = np.arange(kernel.lam_min, kernel.lam_max + step / 2.0, step)
    return {'effective_pressure_GPa': float(effective),
            'lambda_opt_nm': float(grid[int(np.argmax(kernel.a(grid)))])}


# --- estimate C: the [111] group, bounded by chapter 7 --------------------

def g_from_chapter7(pressures=CH7_PRESSURES_GPA, lam=GREEN_NM):
    """Upper bound on g for the [111] group, from chapter 7 having worked.

    Chapter 7 ran cwODMR at these pressures with `lam` excitation at cryogenic
    temperature.  Sub-ZPL absorption is thermally activated (erratum E1) and is
    1e-11 at 50 K over the hydrostatic deficit, so the line must actually lie
    at or above the group's local ZPL.  That is a bound on the red shift.
    """
    energy = nm2eV(lam)
    rows = []
    for pressure in pressures:
        hydrostatic = zpl_energy(pressure) - NV_ZPL_AMBIENT_EV
        deficit = zpl_energy(pressure) - energy
        rows.append({
            'pressure_GPa': float(pressure),
            'hydrostatic_shift_eV': float(hydrostatic),
            'required_red_shift_meV': float(deficit * 1e3),
            'g_upper_bound': float(1.0 - deficit / hydrostatic),
        })
    bound = min(row['g_upper_bound'] for row in rows)
    return {
        'rows': rows,
        'g_upper_bound': bound,
        'lambda_opt_lower_bound_nm': optimum(bound),
        'lambda_opt_lower_bound_compression_nm':
            optimum_by_compression(bound)['lambda_opt_nm'],
    }


# --- estimate B: the [100] culet, measured against Fig. 6.3(b) ------------

def _fit_transmission_at(g, model=None):
    """Anvil factor Fig. 6.3(b) demands if the band sits at `g`."""
    model = HoPublishedSpectrumModel() if model is None else model

    def merit(lam, pressure):
        shift = band_shift_eV(pressure, g)
        value = float(model.sigma_abs(nm2eV(lam) - shift, pressure))
        if value <= 0.0:
            raise ValueError(
                f'the band at g = {g:.3f} carries {lam:.0f} nm out of the '
                f'extracted window at {pressure:.0f} GPa; keep the scan inside '
                f'{G_SCAN_LIMITS}')
        return lam * value

    rows = anvil.observed_ratio()[1:]           # ambient dropped, as the fit does
    pressure = np.array([row['pressure_GPa'] for row in rows])
    observed = np.log10([row['ratio'] for row in rows])
    weight = 1.0 / np.array([row['sigma_log10'] for row in rows]) ** 2
    predicted = np.array([0.5 * np.log10(merit(BLUE_NM, p) / merit(GREEN_NM, p))
                          for p in pressure])
    residual = observed - predicted
    offset = float(np.sum(weight * residual) / np.sum(weight))
    error = float(1.0 / np.sqrt(np.sum(weight)))
    return {
        'g': float(g),
        'chi2': float(np.sum(weight * (residual - offset) ** 2)),
        'transmission_ratio': float(10.0 ** (2.0 * offset)),
        'transmission_ratio_plus_1sigma': float(10.0 ** (2.0 * (offset + error))),
    }


def g_from_fig63(scan=(1.20, 1.10, 1.05, 1.00, 0.90, 0.80, 0.74)):
    """Lower bound on g for the [100] culet, from T(450)/T(532) <= 1.

    Fig. 6.3(b) constrains only the PRODUCT of the band position and the anvil
    factor (`external_audit` X4).  A type Ia anvil absorbs below 500 nm, so it
    cannot transmit blue better than green; that requirement is what turns the
    product into a bound on g alone.
    """
    table = [_fit_transmission_at(g) for g in scan]
    physical = [row for row in table if row['transmission_ratio'] <= 1.0]
    boundary = brentq(
        lambda g: _fit_transmission_at(g)['transmission_ratio'] - 1.0,
        *G_SCAN_LIMITS, xtol=1e-3)
    return {
        'rows': table,
        'g_lower_bound': float(boundary),
        'g_lower_bound_1sigma': float(brentq(
            lambda g: _fit_transmission_at(g)['transmission_ratio_plus_1sigma']
            - 1.0, *G_SCAN_LIMITS, xtol=1e-3)),
        'n_physical_in_scan': len(physical),
    }


# --- estimate A', for comparison only: Hilberer's C-4 --------------------

def g_from_hilberer():
    """C-4's flat-culet factor, and what Fig. 6.3(b) says about applying it."""
    g = float(_alpha_factor(ALPHA_FLAT_CULET))
    return {
        'g': g,
        'lambda_opt_nm': optimum(g),
        'transmission_demanded': _fit_transmission_at(g)['transmission_ratio'],
        'note': 'group-averaged; demands an unphysical anvil on the [100] data',
    }


def geometry_table(pressure=FROZEN_PRESSURE_GPA):
    """The three geometries, and the optimum in each."""
    ch7 = g_from_chapter7()
    fig63 = g_from_fig63()
    return {
        'pressure_GPa': float(pressure),
        'quasi_hydrostatic': {
            'g': 1.0,
            'lambda_opt_nm': optimum(1.0, pressure),
            'basis': 'definition; Ho micropillar, alpha ~ 0.95',
        },
        'flat_culet_100': {
            'g_lower_bound': fig63['g_lower_bound'],
            'lambda_opt_upper_bound_nm': optimum(fig63['g_lower_bound'],
                                                 pressure),
            'basis': 'Fig. 6.3(b) with T(450)/T(532) <= 1',
        },
        'flat_culet_111_axial': {
            'g_upper_bound': ch7['g_upper_bound'],
            'lambda_opt_lower_bound_nm': optimum(ch7['g_upper_bound'],
                                                 pressure),
            'basis': 'chapter 7 ran 532 nm at 140 GPa, cryogenic',
        },
        'frozen_optimum_nm': FROZEN_OPTIMUM_NM,
    }


def structure_under_geometry(g, pressure=FROZEN_PRESSURE_GPA):
    """What the geometry does to A2's ladder, not just to the optimum."""
    kernel = GeometryKernel(g, pressure)
    window = (kernel.lam_min, kernel.lam_max)
    maxima = kernel.local_maxima()
    critical = critical_values(kernel, window)
    return {
        'g': float(g),
        'n_interior_maxima': len(maxima),
        'maxima_nm': [lam for lam, _ in maxima],
        'n_critical_values': len(critical),
        'lowest_rung_I_over_Ic': (float(critical[0].power_ratio)
                                  if critical else None),
        'lowest_rung_nm': float(critical[0].wavelength) if critical else None,
    }


def report():
    out = {}

    print('Which geometry does 440.65 nm belong to?')
    print(f'  at {FROZEN_PRESSURE_GPA:.0f} GPa, band carried rigidly\n')

    print('estimate C  [111] flat culet, [111] group -- bounded by chapter 7')
    out['chapter7'] = ch7 = g_from_chapter7()
    for row in ch7['rows']:
        print(f'      {row["pressure_GPa"]:5.0f} GPa: 532 nm needs a red shift '
              f'of {row["required_red_shift_meV"]:6.1f} meV  ->  '
              f'g <= {row["g_upper_bound"]:.3f}')
    print(f'      binding bound  g <= {ch7["g_upper_bound"]:.3f}')
    print(f'      lambda_opt >= {ch7["lambda_opt_lower_bound_nm"]:.2f} nm '
          f'(rigid shift; {ch7["lambda_opt_lower_bound_compression_nm"]:.2f} nm '
          f'under the compression model, so the rigid one is conservative)')

    print('\nestimate B  [100] flat culet -- measured against Fig. 6.3(b)')
    out['fig63'] = fig63 = g_from_fig63()
    print(f'      {"g":>7}{"chi2":>9}{"T(450)/T(532)":>16}   physical?')
    for row in fig63['rows']:
        mark = 'yes' if row['transmission_ratio'] <= 1.0 else 'no'
        print(f'      {row["g"]:7.3f}{row["chi2"]:9.2f}'
              f'{row["transmission_ratio"]:16.3f}   {mark}')
    print(f'      T = 1 at g = {fig63["g_lower_bound"]:.3f}, so g >= that; '
          f'at +1 sigma, g >= {fig63["g_lower_bound_1sigma"]:.3f}')
    print('      [100] stress blue shifts all four groups (Fig. 6.11(a)), '
          'which is this sign')

    print('\ncomparison  C-4\'s flat-culet factor, applied per group')
    out['hilberer'] = hil = g_from_hilberer()
    print(f'      g = {hil["g"]:.3f} -> lambda_opt {hil["lambda_opt_nm"]:.2f} nm,'
          f' but demands T(450)/T(532) = {hil["transmission_demanded"]:.2f}')
    print('      so C-4 is not extended per group; nv_model stays frozen')

    print('\nthe table')
    out['table'] = table = geometry_table()
    print(f'      {"geometry":34}{"g":>12}{"lambda_opt [nm]":>18}')
    print(f'      {"quasi-hydrostatic (Ho, frozen)":34}'
          f'{table["quasi_hydrostatic"]["g"]:12.3f}'
          f'{table["quasi_hydrostatic"]["lambda_opt_nm"]:18.2f}')
    print(f'      {"[100] flat culet":34}'
          f'{">= " + format(table["flat_culet_100"]["g_lower_bound"], ".3f"):>12}'
          f'{"<= " + format(table["flat_culet_100"]["lambda_opt_upper_bound_nm"], ".2f"):>18}')
    print(f'      {"[111] flat culet, [111] group":34}'
          f'{"<= " + format(table["flat_culet_111_axial"]["g_upper_bound"], ".3f"):>12}'
          f'{">= " + format(table["flat_culet_111_axial"]["lambda_opt_lower_bound_nm"], ".2f"):>18}')

    print('\nwhat survives structurally')
    out['structure'] = {}
    for g in (1.0, ch7['g_upper_bound']):
        out['structure'][round(g, 3)] = row = structure_under_geometry(g)
        print(f'      g = {g:.3f}: {row["n_interior_maxima"]} interior maxima at '
              + ', '.join(f'{lam:.1f}' for lam in row['maxima_nm']) + ' nm')
        print(f'                 {row["n_critical_values"]} critical values, '
              f'lowest rung I/Ic = {row["lowest_rung_I_over_Ic"]:.4f} '
              f'at {row["lowest_rung_nm"]:.1f} nm')
    print('      Theorems M, G and X hold for any kernel, so the level-set and '
          'step-function\n      structure is untouched; the NUMBERS on it are '
          'geometry dependent.')
    return out


if __name__ == '__main__':
    report()
