"""Regression tests for the v1 failure diagnosis."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from repro_yield import HO_MAX_FRACTIONAL_RMS, HO_MAX_PEAK_ERROR_GPA  # noqa: E402
from v1_diagnosis import (  # noqa: E402
    V1_DE120_EV, V1_HW_EV, audit, pekarian_peak, phonon_energy, report,
    wavelength_structure, zpl_shift,
)


@pytest.fixture(scope='module')
def summary():
    return report()


def test_phonon_energy_is_derived_not_fitted():
    """hw comes from Fig. 1(b),(e) alone -- no Fig. 5(b) quantity enters."""
    pressure, curve, mean = phonon_energy()
    assert np.array_equal(pressure, np.arange(0.0, 121.0, 20.0))
    assert mean == pytest.approx(0.1011, abs=5e-4)


def test_phonon_energy_is_roughly_pressure_independent():
    """The correction is a wrong constant, not a missing pressure dependence."""
    _, curve, mean = phonon_energy()
    spread = curve[1:].max() - curve[1:].min()
    assert spread < 0.15 * mean


def test_v1_undersizes_the_phonon_energy():
    _, _, mean = phonon_energy()
    assert mean / V1_HW_EV == pytest.approx(1.56, abs=0.02)


def test_pekarian_maximum_sits_half_a_phonon_below_S():
    """The continuum shortcut p* = S is what made the first attempt fail."""
    for s_abs in (3.0, 3.6, 4.1, 4.6):
        assert s_abs - pekarian_peak(s_abs) == pytest.approx(0.51, abs=0.01)


def test_the_continuum_shortcut_costs_the_gates(summary):
    shortcut = summary['results']['continuum shortcut p*=S']
    assert shortcut['pooled_fractional_rms'] > 0.10
    assert summary['hw_continuum'] == pytest.approx(0.0879, abs=5e-4)


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
    assert corrected['532']['correlation'] > 0.99
    assert corrected['457']['correlation'] > 0.99


def test_corrected_model_matches_the_reconstruction(summary):
    """1.0% pooled is what the reconstructed kernel scores on this test."""
    corrected = summary['results']['both, from Ho panels']
    assert corrected['pooled_fractional_rms'] == pytest.approx(0.010, abs=2e-3)


def test_the_phonon_energy_carries_the_correction(summary):
    """dE120 alone does almost nothing; hw alone does almost all of it."""
    results = summary['results']
    frozen = results['v1 as frozen']['pooled_fractional_rms']
    hw_only = results['hw corrected only']['pooled_fractional_rms']
    dE_only = results['dE120 corrected only']['pooled_fractional_rms']
    assert dE_only > 0.9 * frozen
    assert hw_only < 0.4 * frozen


def test_corrected_v1_passes_both_gates(summary):
    corrected = summary['results']['both, from Ho panels']
    assert summary['gates_pass']
    for line in ('532', '457'):
        assert corrected[line]['fractional_rms'] <= HO_MAX_FRACTIONAL_RMS
        assert corrected[line]['peak_error_GPa'] <= HO_MAX_PEAK_ERROR_GPA


def test_passing_fig5b_does_not_reinstate_v1_as_the_kernel(summary):
    """A two-wavelength slice cannot see the structure A2 is built on."""
    structure = summary['structure']
    assert structure['maxima_kernel'] == 4
    assert structure['maxima_model'] == 1
    assert structure['fractional_rms'] > 0.4
    assert structure['correlation'] > 0.99


def test_optical_limit_optimum_is_robust_to_the_envelope(summary):
    """440.65 nm survives the model change once the constants are right."""
    structure = summary['structure']
    assert abs(structure['lambda_opt_model']
               - structure['lambda_opt_kernel']) < 2.0


def test_audit_is_pure_in_its_parameters():
    a = audit(V1_HW_EV, V1_DE120_EV)
    b = audit(V1_HW_EV, V1_DE120_EV)
    assert a['pooled_fractional_rms'] == b['pooled_fractional_rms']
