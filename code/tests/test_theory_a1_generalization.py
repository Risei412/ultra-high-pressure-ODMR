"""Regression tests for the numerical execution of Addendum A1.

These lock in both the propositions that survive the numerical check and the
three places where the analytic notes and the reconstructed kernel disagree.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from theory_a1_generalization import (  # noqa: E402
    DATA_WINDOW, MAIN_BAND, Kernel, MediatedResponse, eta_at,
    section1_geometry, section2_split_formula, section3_multiplet,
    section4_mechanisms, section4b_mechanism_degeneracy, section5_ordering,
    section6_t4_prediction, section7_hot_band, section8_kernel_provenance,
    sensitivity_optima,
)


@pytest.fixture(scope='module')
def kernel():
    return Kernel()


# ------------------------------------------------------------------ geometry

def test_band_maximum_matches_the_frozen_optical_optimum(kernel):
    """S1: the kernel maximum is the frozen 440.65 nm, within the grid step."""
    assert kernel.lam_abs == pytest.approx(440.65, abs=0.02)


def test_kernel_is_multimodal_and_the_zpl_is_a_serious_competitor(kernel):
    """S1: the note assumes one band; the reconstruction carries four maxima."""
    maxima = kernel.local_maxima()
    positions = [lam for lam, _ in maxima]
    assert len(maxima) == 4
    assert positions == pytest.approx([440.64, 475.55, 500.19, 514.46], abs=0.05)
    zpl_height = dict((round(lam, 2), value) for lam, value in maxima)[514.46]
    assert zpl_height == pytest.approx(0.694, abs=0.005)
    # The ZPL is as good an excitation line as 473 nm in the optical limit.
    assert 1.0 / np.sqrt(zpl_height) == pytest.approx(
        1.0 / np.sqrt(kernel.a(473.0)), rel=0.01)


def test_chord_curvature_reproduces_the_published_table(kernel):
    """S1: the note's kappa column is a chord curvature, and it reproduces."""
    expected = {405.0: 8.40e-4, 420.0: 8.92e-4, 430.0: 1.02e-3, 445.0: 8.00e-4,
                457.0: 6.58e-4, 473.0: 7.14e-4, 488.0: 7.68e-4}
    for lam, value in expected.items():
        assert float(kernel.kappa_chord(lam)) == pytest.approx(value, rel=0.01)


def test_pointwise_curvature_is_unusable_so_kappa_needs_a_stated_window(kernel):
    """S1: the local second derivative is interpolation noise, not curvature."""
    values = kernel.kappa_pointwise(np.arange(420.0, 500.0, 0.05))
    # Three orders of magnitude above the quoted 8e-4, and not even one-signed.
    assert np.max(np.abs(values)) > 1.0
    assert np.sum(np.diff(np.sign(values)) != 0) > 10


def test_kappa_depends_on_the_fit_window(kernel):
    """S1: 8e-4 corresponds to a +/-40 nm fit, not to the 5 % band."""
    assert kernel.kappa_fit(10.0) == pytest.approx(7.01e-4, rel=0.02)
    assert kernel.kappa_fit(40.0) == pytest.approx(8.07e-4, rel=0.02)
    assert kernel.kappa_fit(40.0) > kernel.kappa_fit(10.0)


def test_tolerance_band_is_asymmetric(kernel):
    """S1: the note quotes 440.65 +/- 15.6 nm; the real band is lopsided."""
    lo, hi = kernel.tolerance_band(1.05)
    assert lo == pytest.approx(426.43, abs=0.1)
    assert hi == pytest.approx(457.90, abs=0.1)
    assert (hi - kernel.lam_abs) - (kernel.lam_abs - lo) == pytest.approx(3.05, abs=0.2)


# ----------------------------------------------------------------- P2 (split)

def test_p2_is_accurate_in_the_mid_range():
    """S2: at l_G = 0.4 %/nm the linear formula is good to a few per cent."""
    rows = {row['l_G_per_nm']: row for row in section2_split_formula()['rows']}
    row = rows[0.004]
    assert row['local_shift_nm'] == pytest.approx(row['predicted_shift_nm'], rel=0.05)


def test_p2_fails_below_the_kernel_mesh():
    """S2: a 1.3 nm predicted shift cannot resolve on a 3-10 nm node mesh."""
    rows = {row['l_G_per_nm']: row for row in section2_split_formula()['rows']}
    row = rows[0.0005]
    assert row['predicted_shift_nm'] == pytest.approx(1.31, abs=0.1)
    # Pinned to the extraction node: the argmax does not leave the scan step.
    assert abs(row['local_shift_nm']) <= 0.01
    assert abs(row['local_shift_nm']) < row['predicted_shift_nm'] / 100.0


def test_global_optimum_jumps_to_the_zpl_well_before_the_quoted_threshold():
    """S2: the note's 0.62 %/nm is the band edge; the jump happens at 0.27."""
    result = section2_split_formula()
    assert result['zpl_jump_l_G'] == pytest.approx(0.00268, rel=0.05)
    assert result['zpl_jump_l_G'] < 0.00625
    jumped = [row for row in result['rows'] if row['jumped_to_zpl']]
    assert all(row['global_lambda_eta_nm'] == pytest.approx(514.46, abs=0.1)
               for row in jumped)


# ------------------------------------------------------------ P3 (degeneracy)

def test_doublet_degeneracy_is_exact(kernel):
    """S3: P3(c)'s central claim holds to machine precision."""
    response = MediatedResponse(gamma_contrast=1.0)
    gamma_max = response.gamma_star() / 0.9
    members = sensitivity_optima(kernel, response, gamma_max, MAIN_BAND)['optima']
    assert len(members) == 2
    etas = [float(eta_at(kernel, response, x, gamma_max)) for x in members]
    assert etas[0] == pytest.approx(etas[1], rel=1e-12)


def test_doublet_separation_matches_the_gaussian_formula_where_it_applies():
    """S3: the quoted 32 nm and 60 nm rows reproduce to ~1.5 %."""
    rows = {row['a_star_over_a_max']: row for row in section3_multiplet()['rows']}
    assert rows[0.9]['exact_separation_nm'] == pytest.approx(32.80, abs=0.2)
    assert rows[0.7]['exact_separation_nm'] == pytest.approx(60.55, abs=0.2)
    for ratio in (0.9, 0.7):
        assert rows[ratio]['exact_separation_nm'] == pytest.approx(
            rows[ratio]['gaussian_separation_nm'], rel=0.02)


def test_the_83_nm_row_is_not_verifiable_against_the_published_kernel():
    """S3: at a*/a_max = 0.5 the blue member lies below the data window."""
    rows = {row['a_star_over_a_max']: row for row in section3_multiplet()['rows']}
    assert rows[0.5]['truncated_blue'] is True
    assert np.isnan(rows[0.5]['exact_separation_nm'])


def test_doublet_is_degenerate_but_not_symmetric_in_wavelength():
    """S3: P6 calls the (M)-split symmetric; the red member sits further out."""
    rows = {row['a_star_over_a_max']: row for row in section3_multiplet()['rows']}
    for ratio in (0.9, 0.7):
        row = rows[ratio]
        assert row['red_offset_nm'] > row['blue_offset_nm'] + 2.0


# ---------------------------------------------------------- P4 (mechanisms)

def test_p4_mechanism_table_reproduces():
    """S4: only contrast collapse splits the optima on its own."""
    table = {row['mechanism']: row for row in section4_mechanisms()}
    assert table['rate saturation alone']['splits'] is False
    assert table['power broadening alone']['splits'] is False
    assert table['saturation + broadening']['splits'] is True
    assert table['contrast collapse alone']['splits'] is True
    # The note states Gamma_p* = Gamma_C for contrast collapse alone.
    assert table['contrast collapse alone']['gamma_star'] == pytest.approx(1.0, rel=1e-6)
    assert table['saturation + broadening']['gamma_star'] == pytest.approx(1.0, rel=1e-6)


def test_splitting_mechanisms_are_degenerate_in_eta():
    """S4b: eta cannot say which mechanism split the optima.

    With equal half-scales, saturation + power broadening and contrast collapse
    share Phi exactly, so the whole eta(lambda, I) surface is identical.  This
    is why the freeze's requirement to record contrast and linewidth separately
    is load-bearing rather than a convenience.
    """
    result = section4b_mechanism_degeneracy()
    for label in ('saturation+broadening vs contrast collapse',
                  'saturation alone vs broadening alone'):
        assert result[label]['max_slope_difference'] < 1e-12
        assert result[label]['max_shape_difference'] < 1e-12
    disc = result['discriminators_at_gamma_p_equals_1']
    # They are told apart only by C and dnu measured separately.
    assert disc['contrast: sat+broadening'] == pytest.approx(1.0)
    assert disc['contrast: contrast collapse'] == pytest.approx(0.5)
    assert disc['linewidth: sat+broadening'] == pytest.approx(np.sqrt(2.0))
    assert disc['linewidth: contrast collapse'] == pytest.approx(1.0)


def test_every_mechanism_starts_from_unit_logarithmic_slope():
    """S4: at vanishing power Phi grows like Gamma_p for all mechanisms."""
    for row in section4_mechanisms():
        assert row['slope_at_low_power'] == pytest.approx(1.0, abs=1e-4)


# ------------------------------------------------------------- P5 (ordering)

def test_intermediate_window_exists_but_full_containment_is_much_narrower():
    """S5: 1.43 decades of split, but only 0.27 decades with both members."""
    order = section5_ordering()
    assert order['window_decades'] == pytest.approx(1.43, abs=0.05)
    assert order['containment_decades'] == pytest.approx(0.27, abs=0.03)
    assert order['containment_decades'] < order['window_decades'] / 4.0


def test_lambda_pl_never_moves_under_mediation(kernel):
    """S5: under (M) with monotone R, lambda_PL is pinned to lambda_abs."""
    response = MediatedResponse(gamma_sat=1.0, gamma_contrast=1.0, gamma_width=1.0)
    grid = np.arange(*MAIN_BAND, 0.05)
    for gamma_max in (1e-2, 1.0, 1e2):
        rate = response.rate(gamma_max * kernel.a(grid))
        assert float(grid[int(np.argmax(rate))]) == pytest.approx(
            kernel.lam_abs, abs=0.05)


# -------------------------------------------------------------------- T4

def test_t4_recovers_gamma_star_from_a_single_wavelength_sweep():
    """S6: the inference half of T4 is robust to 2 % measurement noise."""
    result = section6_t4_prediction()
    assert result['gamma_star_fitted'] == pytest.approx(
        result['gamma_star_true'], rel=0.05)
    # The addendum's falsification threshold is a factor of two.
    assert result['ensemble']['within_factor_2'] == 1.0


def test_t4_closed_form_separation_breaks_once_the_doublet_is_truncated():
    """S6: the prediction survives; the Gaussian formula it feeds does not."""
    rows = {row['u_over_u_c']: row for row in section6_t4_prediction()['predictions']}
    contained = rows[1.43]
    assert contained['predicted_separation_nm'] == pytest.approx(
        contained['true_separation_nm'], abs=3.0)
    assert contained['gaussian_separation_nm'] == pytest.approx(
        contained['true_separation_nm'], rel=0.05)
    truncated = rows[2.0]
    assert truncated['truncated_blue'] is True
    assert truncated['gaussian_separation_nm'] > 2.0 * truncated['true_separation_nm']


# -------------------------------------------------------------------- E1, v1

def test_e1_hot_band_arithmetic():
    """S7: the erratum's Boltzmann numbers reproduce exactly."""
    hot = section7_hot_band()
    assert hot['gap_meV'] == pytest.approx(79.5, abs=0.5)
    assert hot['boltzmann_300K'] == pytest.approx(4.6e-2, rel=0.02)
    assert hot['boltzmann_90K'] == pytest.approx(3.7e-5, rel=0.05)
    # The interpolated value sits between the two limits, derived from neither.
    assert hot['boltzmann_90K'] < hot['interpolated_a_532'] < hot['boltzmann_300K']


def test_v1_and_ho_kernels_disagree_by_more_than_the_tolerance_band():
    """S8: 475.5 nm and 440.65 nm are not the same recommendation."""
    prov = section8_kernel_provenance()
    row = prov['per_pressure'][120.0]
    assert row['v1_franck_condon_nm'] == pytest.approx(475.51, abs=0.1)
    assert row['ho_kernel_nm'] == pytest.approx(440.64, abs=0.1)
    assert row['gap_nm'] == pytest.approx(34.87, abs=0.2)
    lo, hi = prov['band_5pct']
    assert not lo <= row['v1_franck_condon_nm'] <= hi


def test_v1_optimum_coincides_with_the_ho_kernel_second_local_maximum():
    """S8: 475.51 nm and 475.55 nm agree to 0.04 nm -- worth reconciling."""
    prov = section8_kernel_provenance()
    assert prov['second_local_max'][0] == pytest.approx(
        prov['per_pressure'][120.0]['v1_franck_condon_nm'], abs=0.1)


def test_planned_laser_lines_against_the_frozen_kernel():
    """S8: of 457/473/488 nm only 457 nm lies inside the frozen 5 % band."""
    prov = section8_kernel_provenance()
    lo, hi = prov['band_5pct']
    assert lo <= 457.0 <= hi
    assert not lo <= 473.0 <= hi
    assert not lo <= 488.0 <= hi
    assert prov['ho_penalty'][473.0] == pytest.approx(1.2054, rel=0.01)
