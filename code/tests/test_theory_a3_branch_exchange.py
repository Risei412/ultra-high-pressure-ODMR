"""Regression tests for Addendum A3 (pressure-driven branch exchange).

The branches are identified in the *raw extracted samples*, so these tests
also serve as a guard that the effect is not an artefact of the pressure
interpolation the A1 audit found unreliable.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from theory_a3_branch_exchange import (  # noqa: E402
    branch_shift_rates, branch_table, coupling_growth, exchange_pressure,
    fixed_wavelength_crossover, identify_branches, monotonicity,
    raw_local_maxima, load_raw_spectra, zero_power_degeneracy,
)


@pytest.fixture(scope='module')
def branches():
    return identify_branches()


# --------------------------------------------------- both branches are real

def test_both_branches_exist_as_raw_maxima_at_every_pressure():
    """X1: neither branch is produced by interpolating between pressures."""
    spectra = load_raw_spectra()
    assert sorted(spectra) == [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0]
    for pressure, values in spectra.items():
        maxima = raw_local_maxima(values)
        low = [m for m in maxima if m[0] < 2.50]
        assert low, f'no ZPL maximum at {pressure} GPa'
        best_low = max(low, key=lambda m: m[1])
        high = [m for m in maxima if m[0] > best_low[0] + 0.20]
        assert high, f'no sideband maximum at {pressure} GPa'


def test_ambient_zpl_branch_lands_on_the_known_nv_zpl(branches):
    """X1: at 0 GPa the ZPL branch is at 637 nm -- the NV- ZPL.

    This is an independent check that the extraction is anchored correctly:
    the value is not fitted anywhere in the pipeline.
    """
    ambient = branches['zpl'][0]
    assert ambient.pressure == 0.0
    assert ambient.energy_ev == pytest.approx(1.945, abs=0.005)
    assert ambient.wavelength_nm == pytest.approx(637.1, abs=0.5)


# ------------------------------------------------------- Theorem X antecedent

def test_log_ratio_is_monotone_with_a_single_sign_change(branches):
    """X2: the antecedent of Theorem X holds, so the crossing is unique."""
    mono = monotonicity(branches)
    assert mono['monotone_increasing'] is True
    assert np.all(mono['increments'] > 0.0)
    assert mono['n_sign_changes'] == 1
    assert mono['total_change'] == pytest.approx(1.466, abs=0.01)


def test_exchange_pressure_and_jump(branches):
    """X3: P* and the size of the discontinuity, at fixed optical power."""
    star = exchange_pressure(branches)
    assert star['pressure'] == pytest.approx(87.9, abs=0.3)
    assert star['zpl_nm'] == pytest.approx(534.8, abs=1.0)
    assert star['sideband_nm'] == pytest.approx(463.5, abs=1.0)
    assert star['jump_nm'] == pytest.approx(71.3, abs=1.5)


def test_exchange_pressure_depends_on_the_power_convention(branches):
    """X3: fixed photon flux moves P* by ~12 GPa, so the convention matters."""
    weighted = exchange_pressure(branches, weighted=True)['pressure']
    flux = exchange_pressure(branches, weighted=False)['pressure']
    assert flux == pytest.approx(75.6, abs=0.5)
    assert weighted - flux == pytest.approx(12.3, abs=1.0)


# ------------------------------------------------------------ the driver

def test_electron_phonon_coupling_grows_monotonically(branches):
    """X4: S*hbar_omega rises 0.232 -> 0.404 eV.  This is what moves them."""
    growth = coupling_growth(branches)
    assert growth['monotone'] is True
    assert growth['gap_ev'][0] == pytest.approx(0.2319, abs=0.001)
    assert growth['gap_ev'][-1] == pytest.approx(0.4037, abs=0.001)
    assert growth['relative_growth'] == pytest.approx(0.74, abs=0.02)
    assert growth['growth_meV_per_GPa'] == pytest.approx(1.27, abs=0.05)


def test_the_two_branches_decay_at_very_different_rates(branches):
    """X5: the ZPL collapses (x6.0) while the sideband barely moves (x1.3).

    This asymmetry -- not the peak positions -- is what makes the ratio cross.
    """
    rates = branch_shift_rates(branches)
    assert rates['zpl']['sigma_falls_by'] == pytest.approx(6.04, rel=0.02)
    assert rates['sideband']['sigma_falls_by'] == pytest.approx(1.34, rel=0.02)
    assert rates['zpl']['sigma_falls_by'] > 4.0 * rates['sideband']['sigma_falls_by']


def test_branches_have_different_pressure_coefficients(branches):
    """X5: 3.8 vs 5.1 meV/GPa -- they are not one feature seen twice."""
    rates = branch_shift_rates(branches)
    assert rates['zpl']['mean_meV_per_GPa'] == pytest.approx(3.83, abs=0.1)
    assert rates['sideband']['mean_meV_per_GPa'] == pytest.approx(5.11, abs=0.1)


def test_zpl_shift_decelerates(branches):
    """X5: the ZPL pressure coefficient falls smoothly, as a real line should."""
    local = branch_shift_rates(branches)['zpl']['local_meV_per_GPa']
    assert np.all(np.diff(local) < 0.0)
    assert local[0] == pytest.approx(5.32, abs=0.1)
    assert local[-1] == pytest.approx(2.76, abs=0.1)


# --------------------------------------------------- distinct from other effects

def test_branch_exchange_is_not_the_fixed_wavelength_crossover(branches):
    """X6: within one kernel the two sit 36 GPa apart."""
    crossover = fixed_wavelength_crossover()
    star = exchange_pressure(branches)['pressure']
    assert crossover == pytest.approx(51.4, abs=0.5)
    assert star - crossover > 30.0


def test_a3_degeneracy_requires_no_power(branches):
    """X7: unlike A2's ladder, the branch pair is degenerate at zero power."""
    degenerate = zero_power_degeneracy(branches)
    assert degenerate['requires_power'] is False
    assert degenerate['separation_nm'] == pytest.approx(71.3, abs=1.5)
    assert len(degenerate['members_nm']) == 2


def test_branch_ordering_actually_inverts(branches):
    """X: the sideband starts below the ZPL and ends above it."""
    rows = branch_table(branches)
    assert rows[0]['A_ratio'] < 1.0
    assert rows[-1]['A_ratio'] > 1.0
    assert rows[0]['pressure'] == 0.0
    assert rows[-1]['pressure'] == 120.0
