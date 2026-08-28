"""Regression tests for the v1 failure diagnosis."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from repro_yield import HO_MAX_FRACTIONAL_RMS, HO_MAX_PEAK_ERROR_GPA  # noqa: E402
from v1_diagnosis import (  # noqa: E402
    V1_DE120_EV, V1_HW_EV, audit, phonon_energy, report, zpl_shift,
)


@pytest.fixture(scope='module')
def summary():
    return report()


def test_phonon_energy_is_derived_not_fitted():
    """hw comes from Fig. 1(b),(e) alone -- no Fig. 5(b) quantity enters."""
    pressure, curve, mean = phonon_energy()
    assert np.array_equal(pressure, np.arange(0.0, 121.0, 20.0))
    assert mean == pytest.approx(0.0879, abs=5e-4)


def test_phonon_energy_is_roughly_pressure_independent():
    """The correction is a wrong constant, not a missing pressure dependence."""
    _, curve, mean = phonon_energy()
    spread = curve[1:].max() - curve[1:].min()
    assert spread < 0.10 * mean


def test_v1_undersizes_the_phonon_energy_by_a_third():
    _, _, mean = phonon_energy()
    assert mean / V1_HW_EV == pytest.approx(1.35, abs=0.02)


def test_kernel_zpl_shift_exceeds_the_v1_anchor():
    """Ho's text quotes >400 meV as a bound; the kernel carries 464 meV."""
    shift = zpl_shift()
    assert shift == pytest.approx(0.464, abs=1e-3)
    assert shift > V1_DE120_EV


def test_frozen_v1_is_anticorrelated_at_532(summary):
    frozen = summary['results']['v1 as frozen']
    assert frozen['532']['correlation'] < 0.0
    assert frozen['pooled_fractional_rms'] > 0.35


def test_both_corrections_remove_the_anticorrelation(summary):
    corrected = summary['results']['both, from Ho panels']
    assert corrected['532']['correlation'] > 0.85
    assert corrected['457']['correlation'] > 0.85


def test_correction_is_a_factor_of_four_in_pooled_rms(summary):
    frozen = summary['results']['v1 as frozen']['pooled_fractional_rms']
    corrected = summary['results']['both, from Ho panels'][
        'pooled_fractional_rms']
    assert frozen / corrected > 3.5


def test_the_phonon_energy_carries_the_correction(summary):
    """dE120 alone does almost nothing; hw alone does almost all of it."""
    results = summary['results']
    frozen = results['v1 as frozen']['pooled_fractional_rms']
    hw_only = results['hw corrected only']['pooled_fractional_rms']
    dE_only = results['dE120 corrected only']['pooled_fractional_rms']
    assert dE_only > 0.9 * frozen
    assert hw_only < 0.4 * frozen


def test_corrected_v1_still_fails_both_gates(summary):
    """The trend is recovered; the peak position is not.  Do not oversell."""
    corrected = summary['results']['both, from Ho panels']
    assert corrected['pooled_fractional_rms'] > HO_MAX_FRACTIONAL_RMS
    assert any(corrected[line]['peak_error_GPa'] > HO_MAX_PEAK_ERROR_GPA
               for line in ('532', '457'))


def test_audit_is_pure_in_its_parameters():
    a = audit(V1_HW_EV, V1_DE120_EV)
    b = audit(V1_HW_EV, V1_DE120_EV)
    assert a['pooled_fractional_rms'] == b['pooled_fractional_rms']
