"""
test_freeze.py
--------------
Freeze tests for nv_model.NVModel.  Run:  python -m pytest tests/ -q
(or  python tests/test_freeze.py  for a plain report).

Three things are locked here:

  1. BACKWARD COMPATIBILITY.  With T = 0, the legacy ZPL parameters
     (Emax = 0.758 eV, P0 = 160 GPa) and the ISC prefactor off, the model must
     reproduce the pre-C-1/C-2/C-3 numbers to machine precision.  The revisions
     are extensions, not a rewrite.

  2. THE PHYSICS OF THE NEW ENVELOPE.  The finite-temperature Franck-Condon
     envelope must be normalised, have first moment S, agree with the discrete
     Struck-Fonger lineshape at integer phonon number, and reduce to the T = 0
     Pekarian.

  3. THE STRUCTURAL CLAIM THE PAPER RESTS ON.  lambda_opt must be invariant
     under every wavelength-independent factor -- the ISC contrast prefactor
     C0(P), the NV0 brightness w0, the linewidth -- and under the
     charge-transfer constants.  If a future edit makes lambda_opt depend on
     them, that claim is void and this test says so.

Plus the literature reproduction suite (repro_literature.py) must stay green.
"""

import os
import sys

import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.special import ive

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nv_model import NVModel, HW, KB, nm2eV, default_randomiser   # noqa: E402


# ---------------------------------------------------------------------------
# 1. backward compatibility with the pre-revision model
# ---------------------------------------------------------------------------

# sigma_abs and eta from the frozen T=0 model, keyed by (P/GPa, lambda/nm).
GOLDEN_T0 = {
    (0, 405): (3.358834733638e-07, 1.130622915648e+05),
    (0, 457): (1.577116702772e-03, 8.856832805116e+01),
    (0, 473): (1.026312625322e-02, 1.072543680499e+01),
    (0, 532): (1.000000000000e+00, 2.151657414560e+00),
    (50, 405): (4.402377472135e-04, 1.593847141112e+03),
    (50, 457): (2.115265262019e-01, 4.171866942156e+00),
    (50, 473): (6.770804811030e-01, 2.726553691247e+00),
    (50, 532): (3.566018325230e+00, 1.289286521577e+00),
    (100, 405): (2.811200244416e-02, 4.017594833794e+01),
    (100, 457): (1.788865839914e+00, 1.907529722228e+00),
    (100, 473): (2.934249237629e+00, 1.513553441330e+00),
    (100, 532): (5.433899882364e-01, 3.186211838357e+00),
    (120, 405): (8.519866283088e-02, 6.240703996900e+00),
    (120, 457): (2.575695198171e+00, 1.651733052410e+00),
    (120, 473): (3.261360955795e+00, 1.476901753598e+00),
    (120, 532): (3.954464125739e-02, 7.608462145806e+00),
    (140, 405): (1.513430263976e-01, 5.345939696024e+00),
    (140, 457): (3.002706009575e+00, 1.575184626197e+00),
    (140, 473): (3.179409084515e+00, 1.532926914209e+00),
    (140, 532): (5.221172928800e-05, 1.675349003914e+02),
}


def legacy():
    """The pre-revision model: T=0 envelope, legacy ZPL, no ISC prefactor,
    no collection efficiency.  Every extension must be off here."""
    return NVModel(Emax=0.758, P0=160.0, T=0.0, isc=False, collection=False)


@pytest.mark.parametrize('key', sorted(GOLDEN_T0))
def test_legacy_T0_is_bit_exact(key):
    P, L = key
    s_ref, e_ref = GOLDEN_T0[key]
    m = legacy()
    assert m.sigma_abs(nm2eV(float(L)), float(P)) == pytest.approx(s_ref, rel=1e-12)
    assert m.eta_lambda(float(L), float(P))[0] == pytest.approx(e_ref, rel=1e-12)


def test_low_temperature_falls_back_to_the_T0_branch():
    """Below ~50 K the thermal correction is negligible and the T=0 form is used."""
    lam = np.linspace(402, 700, 500)
    a = legacy().sigma_abs(nm2eV(lam), 120.)
    b = NVModel(Emax=0.758, P0=160.0, T=30.0, isc=False,
                collection=False).sigma_abs(nm2eV(lam), 120.)
    assert np.allclose(a, b, rtol=0, atol=0)


# ---------------------------------------------------------------------------
# 2. the finite-temperature Franck-Condon envelope
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('T', [0.0, 30.0, 90.0, 300.0])
@pytest.mark.parametrize('S', [3.08, 4.61])
def test_envelope_is_normalised_with_first_moment_S(T, S):
    m = NVModel(T=T)
    p = np.arange(-80.0, 120.0, 0.002)
    L = m._fc(p * HW, S)
    area = np.trapezoid(L, p)
    mean = np.trapezoid(L * p, p) / area
    # The continuum approximation of a discrete sum: exact to ~2% in area.
    assert area == pytest.approx(1.0, abs=0.02)
    assert mean == pytest.approx(S, rel=0.02)


@pytest.mark.parametrize('T', [90.0, 300.0, 600.0])
def test_envelope_matches_discrete_struck_fonger(T):
    """L_p = e^{-S(2n+1)} ((n+1)/n)^{p/2} I_|p|(2 S sqrt(n(n+1))) at integer p."""
    S = 4.61
    m = NVModel(T=T)
    nbar = 1.0 / np.expm1(HW / (KB * T))
    z = 2.0 * S * np.sqrt(nbar * (nbar + 1.0))
    for p in range(-8, 20):
        ref = (np.exp(-S * (2 * nbar + 1)) * ((nbar + 1) / nbar) ** (p / 2)
               * ive(abs(p), z) * np.exp(z))
        assert float(m._fc(np.array(p * HW), S)) == pytest.approx(ref, rel=1e-10)


def test_anti_stokes_absorption_exists_only_at_finite_T():
    """The T=0 hard cut below the ZPL is what issue C-1 removed."""
    x = -0.05                      # 50 meV BELOW the ZPL
    assert NVModel(T=0.0)._fc(np.array(x), 4.61) == 0.0
    assert NVModel(T=300.0)._fc(np.array(x), 4.61) > 0.0


def test_sigma_at_532_is_continuous_across_the_ZPL_crossing():
    """
    The frozen T=0 model dropped sigma_abs(532 nm) by 4.6x discontinuously where
    the ZPL crosses 532 nm (113.7 GPa with the legacy parameters).  At room
    temperature the crossing must leave no step.
    """
    m = NVModel(T=300.0)
    P = np.linspace(100., 130., 601)
    s = m.sigma_abs(nm2eV(532.), P)
    step = np.max(np.abs(np.diff(np.log(s))))
    assert step < 0.05, f'log-sigma jumps by {step:.3f} across the crossing'


# ---------------------------------------------------------------------------
# 3. lambda_opt must not depend on wavelength-independent factors
# ---------------------------------------------------------------------------

BASE_T = 300.0


def test_lambda_opt_is_invariant_under_the_ISC_prefactor():
    """
    C-3: C0(P) is the ISC contrast prefactor.  It is wavelength independent, so
    the paper's claim `lambda_opt = argmax sigma_abs` cannot depend on it.
    """
    ref = NVModel(T=BASE_T).lambda_opt(120.)
    for kw in ({'isc': False},
               {'C_amb': 0.05}, {'C_amb': 0.9},
               {'E_isc': 0.05}, {'E_isc': 1.0},
               {'C_floor': 0.05}):
        assert NVModel(T=BASE_T, **kw).lambda_opt(120.) == pytest.approx(ref, abs=1e-9)


def test_lambda_opt_is_invariant_under_the_charge_transfer_constants():
    ref = NVModel(T=BASE_T).lambda_opt(120.)
    for kw in ({'a_gs': 0.0}, {'a_gs': 20.0},
               {'a_es': 0.3}, {'a_es': 3.0},
               {'r0': 0.5}, {'r0': 4.0},
               {'rbg': 0.01}, {'rbg': 1.0},
               {'w0': 0.5}, {'w0': 2.0}):
        assert NVModel(T=BASE_T, **kw).lambda_opt(120.) == pytest.approx(ref, abs=0.05)


def test_lambda_opt_equals_the_absorption_maximum():
    m = NVModel(T=BASE_T)
    lam = np.arange(402., 640., 0.01)
    s = m.sigma_abs(nm2eV(lam), 120.)
    assert m.lambda_opt(120.) == pytest.approx(lam[s.argmax()], abs=0.05)


# ---------------------------------------------------------------------------
# 4. the C-2 ZPL parameterisation
# ---------------------------------------------------------------------------

def test_zpl_reproduces_its_two_measured_anchors():
    m = NVModel()
    slope0 = (m.ZPL(1e-5) - m.ZPL(0.0)) / 1e-5
    assert slope0 == pytest.approx(5.75e-3, rel=2e-3)        # Doherty et al. 2014
    assert m.dZPL(120.) == pytest.approx(0.400, rel=1e-6)    # Ho et al. 2026


@pytest.mark.parametrize('dE120,slope0', [(0.36, 5.0e-3), (0.40, 5.75e-3),
                                          (0.44, 6.5e-3), (0.40, 8.0e-3)])
def test_zpl_shape_solver_round_trips(dE120, slope0):
    m = NVModel(dE120=dE120, slope0=slope0)
    assert m.dZPL(120.) == pytest.approx(dE120, rel=1e-8)
    assert (m.ZPL(1e-6) - m.ZPL(0.)) / 1e-6 == pytest.approx(slope0, rel=1e-4)


def test_zpl_shift_saturates():
    """Doherty's linear extrapolation was wrong (Lyapin 2018, Dai 2022)."""
    m = NVModel()
    s0 = (m.ZPL(1e-5) - m.ZPL(0.)) / 1e-5
    s120 = (m.ZPL(120. + 1e-5) - m.ZPL(120.)) / 1e-5
    assert s120 < 0.6 * s0


def test_mc_band_is_tighter_than_the_legacy_randomisation():
    """C-2: randomising the measured quantities, not (Emax, P0), narrows the band."""
    rng = np.random.default_rng(0)
    new = np.array([default_randomiser(rng, T=BASE_T).lambda_opt(120.) for _ in range(120)])
    from nv_model import legacy_randomiser
    rng = np.random.default_rng(0)
    old = np.array([legacy_randomiser(rng).lambda_opt(120.) for _ in range(120)])
    assert new.std() < 0.8 * old.std()


# ---------------------------------------------------------------------------
# 5. headline numbers and the literature suite
# ---------------------------------------------------------------------------

def test_headline_numbers():
    m = NVModel(T=BASE_T)
    assert m.lambda_opt(120.) == pytest.approx(474.0, abs=1.0)
    # 473 nm DPSS must remain within 1% of the optimum
    pen = m.eta_lambda(473., 120.)[0] / m.eta_lambda(m.lambda_opt(120.), 120.)[0]
    assert pen < 1.01
    # green must survive to megabar (this is what C-1 fixed)
    for P, lim in [(120., 4.0), (140., 5.0), (150., 6.0)]:
        g = m.eta_lambda(532., P)[0] / m.eta_lambda(m.lambda_opt(P), P)[0]
        assert g < lim, f'eta(532)/eta_opt = {g:.1f} at {P} GPa'


def test_literature_reproduction_suite():
    import repro_literature
    repro_literature._rows.clear()
    npass, ntest = repro_literature.main()
    assert npass == ntest, f'{ntest - npass} literature target(s) not reproduced'


if __name__ == '__main__':
    raise SystemExit(pytest.main([os.path.abspath(__file__), '-q']))


# ---------------------------------------------------------------------------
# 6. C-4: anvil geometry (stress anisotropy)
# ---------------------------------------------------------------------------

def test_alpha_factor_is_calibrated_to_hilberer():
    from nv_model import _alpha_factor, ALPHA_REF
    assert _alpha_factor(ALPHA_REF) == pytest.approx(1.0, abs=1e-12)
    # -434 / -769 meV/(cm^3 mol^-1), standard culet vs micropillar
    assert _alpha_factor(0.56) == pytest.approx(434.0 / 769.0, abs=2e-3)


def test_deviatoric_stress_moves_lambda_opt_red():
    """
    Deviatoric stress red-shifts the 3E and so REDUCES the ZPL blue shift; the
    optimum must move back toward the green.  This is the C-4 scope condition:
    the 473 nm answer belongs to the quasi-hydrostatic geometry only.
    """
    lo_pillar = NVModel(T=BASE_T, alpha=0.95).lambda_opt(120.)
    lo_flat = NVModel(T=BASE_T, alpha=0.56).lambda_opt(120.)
    assert lo_flat > lo_pillar + 25.0
    # and in a flat culet the recommendation inverts: green beats 473 nm
    mflat = NVModel(T=BASE_T, alpha=0.56)
    assert (mflat.eta_lambda(532., 120.)[0]
            < mflat.eta_lambda(473., 120.)[0])


def test_alpha_is_monotone_in_lambda_opt():
    los = [NVModel(T=BASE_T, alpha=a).lambda_opt(120.)
           for a in (0.50, 0.60, 0.70, 0.80, 0.90, 1.00)]
    assert all(b < a for a, b in zip(los, los[1:]))


# ---------------------------------------------------------------------------
# 7. Sec. V: the intensity-explicit model
# ---------------------------------------------------------------------------

def _power_model(**kw):
    from nv_model_power import NVModelPower
    return NVModelPower(T=BASE_T, **kw)


def test_nv0_envelope_lies_BLUE_of_the_nv_minus_envelope():
    """
    NV0 ZPL 575 nm (2.156 eV) is 0.21 eV ABOVE NV- 637 nm (1.945 eV), so the NV0
    absorption envelope sits at SHORTER wavelength.  The paper's prose once said
    'displaced to the red', and that sign error is what produced the wrong
    physical expectation for the power dependence.  Guard it.
    """
    m = _power_model()
    lam = np.arange(402., 640., 0.5)
    s = m.sigma_abs(nm2eV(lam), 120.)
    s0 = m.sigma_NV0(nm2eV(lam), 120.)
    assert lam[s0.argmax()] < lam[s.argmax()] - 20.0


def test_power_model_returns_the_fixed_power_optimum_as_u_goes_to_zero():
    """The regression the paper demanded before any power claim could be made."""
    m = _power_model()
    ref = NVModel(T=BASE_T).lambda_opt(120.)
    lam = np.arange(402., 640., 0.1)
    for u in (1e-4, 1e-3, 1e-2):
        e = np.asarray(m.eta_u([(lam, 1.0)], 120., u)[0])
        assert lam[np.nanargmin(e)] == pytest.approx(ref, abs=0.5)


def test_saturation_reference_is_half_saturation_at_the_optimum():
    m = _power_model()
    nE = m.eta_u([(m.lam_sat_ref, 1.0)], m.sat_ref_P, 1.0)[4]
    assert float(nE) == pytest.approx(0.5, abs=1e-9)


def test_photon_rate_is_linear_in_power_at_low_u():
    """Dai et al. 2022 measured linear PL vs power at 100 GPa."""
    m = _power_model()
    u = np.logspace(-4, 2, 400)
    R = np.asarray(m.eta_u([(532., 1.0)], 100., u)[2])
    lin = (R / u) / (R[0] / u[0])
    u10 = u[np.argmax(lin < 0.9)]
    assert u10 > 0.05


def test_ridge_moves_blue_with_power():
    """
    With the NV0 envelope correctly placed blue of the NV- envelope, raising the
    power moves lambda_opt BLUE, pinned finally by the ground-state ionization
    edge -- not red as the paper's prose predicted.
    """
    m = _power_model()
    lam = np.arange(402., 640., 0.5)
    ridge = [lam[np.nanargmin(np.asarray(m.eta_u([(lam, 1.0)], 120., u)[0]))]
             for u in (1e-3, 0.1, 1.0, 10.0)]
    assert all(b <= a for a, b in zip(ridge, ridge[1:]))
    assert ridge[-1] < ridge[0] - 50.0
    # it cannot go blueward of the IP(3A2) edge
    edge = 1239.84 / NVModel(T=BASE_T).IP_A2(120.)
    assert ridge[-1] >= edge - 1.0


# ---------------------------------------------------------------------------
# 8. C-7: emission collection efficiency
# ---------------------------------------------------------------------------

def test_collection_efficiency_falls_monotonically_with_pressure():
    """
    The NV- emission band blue shifts with the ZPL (peak 719 nm at ambient,
    619 nm at 120 GPa), so a passband chosen at ambient pressure collects a
    shrinking fraction.  This is what a filtered count rate sees and a
    spectrally integrated PL yield does not.
    """
    m = NVModel(T=BASE_T)
    e = [m.eta_col(P) for P in (0., 25., 50., 75., 100., 120., 140.)]
    assert all(b < a for a, b in zip(e, e[1:]))
    assert e[0] > 0.75 and e[-1] < 0.30


def test_collection_efficiency_cannot_move_lambda_opt():
    """
    Emission follows Kasha's rule: it leaves the relaxed 3E whatever vibronic
    level was excited, so eta_col is wavelength independent.  Every frozen
    number must be identical with and without it.
    """
    on = NVModel(T=BASE_T, collection=True)
    off = NVModel(T=BASE_T, collection=False)
    assert on.lambda_opt(120.) == pytest.approx(off.lambda_opt(120.), abs=1e-9)
    for band in ((650., 800.), (600., 800.), (700., 900.)):
        m = NVModel(T=BASE_T, det_band=band)
        assert m.lambda_opt(120.) == pytest.approx(off.lambda_opt(120.), abs=1e-9)
        # and eta ratios at fixed pressure are untouched
        r_on = m.eta_lambda(532., 120.)[0] / m.eta_lambda(473., 120.)[0]
        r_off = off.eta_lambda(532., 120.)[0] / off.eta_lambda(473., 120.)[0]
        assert r_on == pytest.approx(r_off, rel=1e-12)


def test_emission_envelope_uses_its_own_huang_rhys_factor():
    m = NVModel()
    assert m.Sem(0.) == pytest.approx(3.39)       # Ho et al. 2026
    assert m.Sem(120.) == pytest.approx(5.25)
    # emission lies BELOW the ZPL, absorption above it
    E = np.linspace(1.2, 3.2, 2001)
    z = m.ZPL(120.)
    assert E[m.sigma_em(E, 120.).argmax()] < z
    assert E[m.sigma_abs(E, 120.).argmax()] > z


def test_widening_the_passband_recovers_sensitivity_at_high_pressure():
    """
    Moving the excitation to the blue is what frees the detection cut-on to move
    blue too; the two gains compound.
    """
    m = NVModel(T=BASE_T)
    gain = np.sqrt(m.eta_col(120., (600., 800.)) / m.eta_col(120., (650., 800.)))
    assert gain > 1.4


# ---------------------------------------------------------------------------
# 9. B-3: is the contrast wavelength independent?
# ---------------------------------------------------------------------------

def test_orbital_branch_analysis_leaves_lambda_opt_in_the_window():
    """
    The central claim lambda_opt = argmax sigma_abs requires C to be wavelength
    independent.  The threat is stress-split 3E orbital branches with different
    ISC rates.  analysis_C_lambda quantifies it: the splitting delta in the
    micropillar geometry must stay below the threshold at which a completely
    dark second branch would push lambda_opt out of the 5% window.
    """
    import analysis_C_lambda as a
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        thresh, d_lo, d_hi = a.main()
    assert d_hi < thresh, (f'orbital splitting {d_hi:.0f} meV exceeds the '
                           f'{thresh:.0f} meV threshold')


def test_worst_case_branch_weighting_stays_inside_the_tolerance_window():
    """
    Explicit worst case: no orbital relaxation at all, and the two branches'
    contrasts differing by any factor between 0 and 10.  At the micropillar
    splitting lambda_opt must not leave the 5% window.
    """
    import analysis_C_lambda as a
    m = NVModel(T=BASE_T)
    P, lam = 120.0, np.arange(402., 640., 0.05)
    E = nm2eV(lam)
    base = np.asarray(m.eta_lambda(lam, P)[0])
    w5 = lam[base <= 1.05 * base.min()]
    z = m.ZPL(P)
    for modulus in (1100.0, 578.0):
        delta = a.splitting_meV(P, 0.95, modulus)
        sA = m._sigma_raw(E, P, zpl=z - delta / 2e3) / m._norm
        sB = m._sigma_raw(E, P, zpl=z + delta / 2e3) / m._norm
        r = sA / (sA + sB)
        for rho in (0.0, 0.1, 0.5, 2.0, 10.0):
            lo = lam[(base / (r + (1.0 - r) * rho)).argmin()]
            assert w5.min() <= lo <= w5.max(), (
                f'delta={delta:.0f} meV, rho={rho}: lambda_opt={lo:.1f} nm '
                f'outside {w5.min():.1f}-{w5.max():.1f}')


def test_isc_is_far_slower_than_orbital_relaxation():
    """The timescale separation that makes C(lambda) flat in the first place."""
    import analysis_C_lambda as a
    assert 1.0 / (a.GAMMA_ISC * a.TAU_ORB) > 1e3


# ---------------------------------------------------------------------------
# 10. The threshold pressure of the blue advantage (Sec. IV G)
# ---------------------------------------------------------------------------

def _ratio(m, P):
    """eta(532)/eta(473) at pressure P; > 1 means blue wins."""
    return (float(np.asarray(m.eta_lambda(532., P)[0]))
            / float(np.asarray(m.eta_lambda(473., P)[0])))


def test_blue_advantage_reverses_below_the_crossover():
    """
    The central recommendation is bounded in pressure.  Below ~71 GPa the model
    must return GREEN as the better choice -- this is what makes the 50 GPa null
    result of the Bhattacharyya thesis a reproduction rather than a failure.
    """
    m = NVModel(T=BASE_T)
    assert _ratio(m, 20.) < 0.5           # green ahead by more than 2x
    assert _ratio(m, 50.) < 1.0           # green still ahead where it was tested
    assert _ratio(m, 120.) > 2.0          # blue clearly ahead at the target
    xo = brentq(lambda P: _ratio(m, P) - 1.0, 20., 145.)
    assert 65. < xo < 78.


def test_crossover_is_monotonic_in_the_blue_wavelength():
    """
    A line further from lambda_opt must wait for a larger shift of the envelope
    before it pays, so the crossover pressure increases as the blue line moves
    away from the optimum.
    """
    m = NVModel(T=BASE_T)

    def xo(lam_b):
        return brentq(lambda P: (float(np.asarray(m.eta_lambda(532., P)[0]))
                                 / float(np.asarray(m.eta_lambda(lam_b, P)[0]))) - 1.0,
                      20., 145.)
    assert xo(473.) < xo(457.) < xo(450.)


def test_crossover_survives_the_charge_transfer_constants():
    """
    The crossover is the prediction offered as falsifiable, so it must not
    depend on the phenomenological knobs any more than lambda_opt does.
    """
    ref = brentq(lambda P: _ratio(NVModel(T=BASE_T), P) - 1.0, 20., 145.)
    for kw in ({'a_gs': 0.0}, {'a_gs': 12.0}, {'r0': 0.5}, {'r0': 4.0},
               {'rbg': 1.0}, {'w0': 0.5}, {'w0': 2.0},
               {'C_amb': 0.12}, {'E_isc': 0.25}):
        m = NVModel(T=BASE_T, **kw)
        assert abs(brentq(lambda P: _ratio(m, P) - 1.0, 20., 145.) - ref) < 2.0


def test_collection_efficiency_is_vectorised_over_pressure():
    """
    eta_col used to accept only scalar pressures, which silently broke
    fig1_green_blue_mix.py.  Array and scalar calls must agree.
    """
    m = NVModel(T=BASE_T)
    P = np.array([0., 37.5, 75., 120.])
    arr = np.asarray(m.eta_col(P))
    assert arr.shape == P.shape
    for p, a in zip(P, arr):
        assert a == pytest.approx(m.eta_col(float(p)), rel=1e-12)


def test_power_penalty_of_the_headline_wavelength():
    """
    Sec. V: the 473 nm recommendation is a low-power statement.  The penalty for
    holding it at finite power is what the paper quotes, and it must remain
    negligible at u = 0.1 and substantial by u = 0.3.
    """
    mp = _power_model()
    lam = np.arange(402., 560., 0.1)

    def penalty(u):
        eta = np.asarray(mp.eta_lambda_u(lam, 120., u)[0])
        return float(np.asarray(mp.eta_lambda_u(473., 120., u)[0])) / eta.min()
    assert penalty(0.05) < 1.01
    assert penalty(0.10) < 1.05
    assert 1.2 < penalty(0.20) < 1.4
    assert 1.5 < penalty(0.30) < 1.8


def test_power_randomiser_keeps_the_saturation_reference_unit():
    """
    The MC randomiser must randomise the DIMENSIONLESS residual s_d, not the
    absolute Gamma_d.  Overriding Gamma_d puts every draw back on the arbitrary
    'green at ambient' cross-section unit -- the Sec. V bug -- which pushed the
    sampled eta(lambda,u) ridge blue of the central curve, so that the 16-84%
    band no longer bracketed the curve it was drawn around.
    """
    from nv_model_power import default_randomiser_power
    rng = np.random.default_rng(11)
    for _ in range(5):
        mm = default_randomiser_power(rng)
        assert mm.Gamma_d == pytest.approx(
            mm.s_d * mm.sigma_abs(nm2eV(mm.lam_sat_ref), mm.sat_ref_P), rel=1e-12)
        assert mm.Gamma_d0 == pytest.approx(
            mm.s_d0 * mm.sigma_NV0(nm2eV(mm.lam_sat_ref), mm.sat_ref_P), rel=1e-12)
        # u = 1 half-saturates (up to the residual s_d) at the reference point
        nE = mm.eta_u([(mm.lam_sat_ref, 1.0)], mm.sat_ref_P, 1.0)[4]
        assert float(nE) == pytest.approx(1.0 / (1.0 + mm.s_d), abs=1e-9)
