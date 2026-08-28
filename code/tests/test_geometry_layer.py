"""Regression tests for the geometry layer: which geometry owns 440.65 nm."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geometry_layer import (  # noqa: E402
    FROZEN_OPTIMUM_NM, GeometryKernel, band_shift_eV, g_from_chapter7,
    g_from_fig63, g_from_hilberer, geometry_table, optimum,
    optimum_by_compression, report, structure_under_geometry,
)
from theory_a1_generalization import Kernel  # noqa: E402


@pytest.fixture(scope='module')
def summary():
    return report()


# --- the parameterisation reduces to the freeze ---------------------------

def test_g_of_one_is_the_frozen_kernel():
    """The whole layer must be a no-op at Ho's geometry."""
    assert band_shift_eV(120.0, 1.0) == 0.0
    assert optimum(1.0) == pytest.approx(FROZEN_OPTIMUM_NM, abs=0.02)


def test_g_of_one_reproduces_the_frozen_window_and_ladder():
    plain, shifted = Kernel(), GeometryKernel(1.0)
    assert shifted.lam_min == pytest.approx(plain.lam_min, abs=1e-6)
    assert shifted.lam_max == pytest.approx(plain.lam_max, abs=1e-6)
    assert shifted.a_max == pytest.approx(plain.a_max, rel=1e-12)
    structure = structure_under_geometry(1.0)
    assert structure['n_interior_maxima'] == 4
    # A2's frozen worked example: the lowest rung is the ZPL at I/Ic = 1.4414.
    assert structure['lowest_rung_I_over_Ic'] == pytest.approx(1.4414, abs=0.001)
    assert structure['lowest_rung_nm'] == pytest.approx(514.46, abs=0.1)


def test_the_shift_has_the_sign_its_name_says():
    assert band_shift_eV(120.0, 0.7) < 0.0        # red shift
    assert band_shift_eV(120.0, 1.3) > 0.0        # further blue shift


def test_the_optimum_moves_red_as_g_falls():
    values = [optimum(g) for g in (1.1, 1.0, 0.9, 0.8, 0.741)]
    assert values == sorted(values)


# --- estimate C: chapter 7 bounds the [111] group -------------------------

def test_chapter_seven_bounds_g_from_above(summary):
    ch7 = summary['chapter7']
    assert ch7['g_upper_bound'] == pytest.approx(0.741, abs=0.005)
    assert all(row['g_upper_bound'] < 1.0 for row in ch7['rows'])


def test_the_binding_bound_is_the_highest_pressure(summary):
    rows = summary['chapter7']['rows']
    assert min(row['g_upper_bound'] for row in rows) == rows[-1]['g_upper_bound']


def test_the_optimum_in_the_working_geometry_is_red_of_the_frozen_one(summary):
    ch7 = summary['chapter7']
    assert ch7['lambda_opt_lower_bound_nm'] == pytest.approx(460.3, abs=0.5)
    assert ch7['lambda_opt_lower_bound_nm'] > FROZEN_OPTIMUM_NM + 15.0


def test_457_nm_is_inside_the_band_the_bound_opens(summary):
    """Ho's practical line sits where the [111] optimum starts."""
    assert 455.0 < summary['chapter7']['lambda_opt_lower_bound_nm'] < 465.0


def test_the_rigid_shift_is_the_conservative_model(summary):
    """The compression model moves the optimum further red, not less."""
    ch7 = summary['chapter7']
    assert (ch7['lambda_opt_lower_bound_compression_nm']
            > ch7['lambda_opt_lower_bound_nm'])


def test_the_compression_model_maps_into_the_published_range():
    result = optimum_by_compression(0.741)
    assert 0.0 < result['effective_pressure_GPa'] < 120.0


# --- estimate B: Fig. 6.3(b) bounds the [100] culet -----------------------

def test_fig63_bounds_g_from_below(summary):
    fig63 = summary['fig63']
    assert fig63['g_lower_bound'] == pytest.approx(1.048, abs=0.01)
    assert fig63['g_lower_bound_1sigma'] > fig63['g_lower_bound']


def test_a_physical_anvil_excludes_a_red_shifted_100_culet(summary):
    """Every g below the boundary demands T(450)/T(532) > 1."""
    for row in summary['fig63']['rows']:
        if row['g'] < summary['fig63']['g_lower_bound']:
            assert row['transmission_ratio'] > 1.0


def test_the_two_anchors_straddle_unity(summary):
    """Which is why one scalar cannot describe both cuts."""
    assert summary['fig63']['g_lower_bound'] > 1.0
    assert summary['chapter7']['g_upper_bound'] < 1.0


# --- C-4 is not extended --------------------------------------------------

def test_c4s_factor_would_demand_an_unphysical_anvil(summary):
    hilberer = summary['hilberer']
    assert hilberer['g'] == pytest.approx(0.564, abs=0.002)
    assert hilberer['transmission_demanded'] > 2.5


def test_c4s_factor_lands_near_v1s_withdrawn_optimum(summary):
    """A coincidence, recorded so nobody reads it as a revival of v1's 475.5 nm."""
    assert summary['hilberer']['lambda_opt_nm'] == pytest.approx(474.9, abs=0.5)


# --- the table and the structure -----------------------------------------

def test_the_table_carries_all_three_geometries(summary):
    table = summary['table']
    assert table['quasi_hydrostatic']['lambda_opt_nm'] == pytest.approx(
        FROZEN_OPTIMUM_NM, abs=0.02)
    assert table['flat_culet_100']['lambda_opt_upper_bound_nm'] < FROZEN_OPTIMUM_NM
    assert table['flat_culet_111_axial']['lambda_opt_lower_bound_nm'] > 455.0


def test_the_structure_survives_the_geometry(summary):
    """Theorems M/G/X hold for any kernel; only the numbers on them move."""
    for row in summary['structure'].values():
        assert row['n_interior_maxima'] == 4
        assert row['n_critical_values'] == 8
        assert 1.4 < row['lowest_rung_I_over_Ic'] < 1.5


def test_but_the_rung_ratios_do_move(summary):
    """Not invariant: A = lambda sigma reweights when the band travels."""
    rungs = [row['lowest_rung_I_over_Ic'] for row in summary['structure'].values()]
    assert len(set(round(value, 4) for value in rungs)) == len(rungs)
