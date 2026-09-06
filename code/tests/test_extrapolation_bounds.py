"""
test_extrapolation_bounds.py
----------------------------
Pins the C-8 extrapolation policy and the claims of extrapolation_bounds.py.

Four things are locked here:

  1. THE DEFAULT IS UNCHANGED.  extrapolate='clip' must reproduce the inline
     np.clip(P, 0, 120) the model used to carry, so no frozen number moves.
  2. THE POLICY IS REAL.  'linear' must continue the anchors and 'error' must
     refuse, the way HoPublishedSpectrumModel does.
  3. THE SPREAD IS THE RESULT.  lambda_opt above the anchors is undetermined,
     and the width of that must grow with pressure.  If a future edit makes it
     collapse, the extrapolation has acquired information from somewhere and
     this test says so.
  4. THE BOUND SURVIVES, WITH ITS PRECONDITION.  lambda_opt >= hc/IP(3A2) in
     every variant that continues the ZPL and IP(3A2) under the same law, and
     once the optimum is pinned to that wall it stops depending on a_gs.  The
     mixed variants, which freeze one of that coupled pair, break it; they are
     flagged rather than quietly dropped.
"""

import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import extrapolation_bounds as eb                                # noqa: E402
from extrapolation_bounds import (DEFAULT_PRESSURES, HW_HARDENING,  # noqa: E402
                                  LinearZPLModel, build, lambda_ion,
                                  lambda_opt, spread, variant_table)
from nv_model import ANCHOR_MAX_P, HBARC, HW, NVModel, nm2eV      # noqa: E402


# ---------------------------------------------------------------------------
# 1. the default policy is byte-identical to the old inline clip
# ---------------------------------------------------------------------------

ANCHORED = [0.0, 37.5, 60.0, 119.9, 120.0]
EXTRAPOLATED = [120.1, 200.0, 300.0, 400.0]


@pytest.mark.parametrize('P', ANCHORED + EXTRAPOLATED)
def test_clip_policy_reproduces_the_inline_clip(P):
    m = NVModel()
    clipped = np.clip(P, 0, 120)
    assert m.Sabs(P) == pytest.approx(3.08 + m.S_slope * clipped / 120.0, rel=0, abs=0)
    assert m.Sem(P) == pytest.approx(3.39 + (5.25 - 3.39) * clipped / 120.0, rel=0, abs=0)
    assert m.IP_A2(P) == pytest.approx(2.68 + (3.06 - 2.68) * clipped / 120.0, rel=0, abs=0)
    assert m.IP_E(P) == pytest.approx(1.16 + (1.63 - 1.16) * clipped / 120.0, rel=0, abs=0)


@pytest.mark.parametrize('P', EXTRAPOLATED)
def test_clip_freezes_every_anchor_at_the_edge(P):
    """The complaint that motivated C-8: above 120 GPa nothing but dZPL moves."""
    m = NVModel()
    assert m.Sabs(P) == m.Sabs(ANCHOR_MAX_P)
    assert m.Sem(P) == m.Sem(ANCHOR_MAX_P)
    assert m.IP_A2(P) == m.IP_A2(ANCHOR_MAX_P)
    assert m.IP_E(P) == m.IP_E(ANCHOR_MAX_P)
    # ... and dZPL, the one that does move, is nearly saturated.
    assert m.dZPL(300.0) / m.Emax > 0.95
    assert m.dZPL(400.0) / m.Emax > 0.98


def test_the_two_dwf_factors_are_frozen_too():
    """Both Debye-Waller factors route through the policy (nv_model:_sigma_raw,
    :sigma_em), so the whole absorption and emission lineshape is frozen."""
    m = NVModel()
    for P in EXTRAPOLATED:
        frozen = m._anchor_P(P)
        assert frozen == pytest.approx(ANCHOR_MAX_P, rel=1e-12)
        abs_dwf = 0.022 * np.exp(-(frozen / 120.0) * np.log(0.022 / 0.0036))
        em_dwf = 0.049 * np.exp(-(frozen / 120.0) * np.log(0.049 / 0.008))
        assert abs_dwf == pytest.approx(0.0036, rel=1e-12)
        assert em_dwf == pytest.approx(0.008, rel=1e-12)


def test_golden_number_unmoved_by_the_refactor():
    """A frozen value from tests/test_freeze.py's GOLDEN_T0 grid."""
    m = NVModel(T=0.0, Emax=0.758, P0=160.0, isc=False, collection=False,
                intensity_basis='photon_flux')
    assert float(m.sigma_abs(nm2eV(532.0), 0.0)) == pytest.approx(1.0, rel=1e-12)


# ---------------------------------------------------------------------------
# 2. the other two policies do what they say
# ---------------------------------------------------------------------------

def test_linear_policy_continues_the_anchor_laws():
    m = NVModel(extrapolate='linear')
    assert m.Sabs(240.0) == pytest.approx(3.08 + 2 * (4.61 - 3.08), rel=1e-12)
    assert m.IP_A2(240.0) == pytest.approx(2.68 + 2 * (3.06 - 2.68), rel=1e-12)
    # still clamped at zero from below
    assert m.Sabs(-10.0) == pytest.approx(3.08, rel=1e-12)


@pytest.mark.parametrize('P', EXTRAPOLATED)
def test_error_policy_refuses_above_the_anchors(P):
    m = NVModel(extrapolate='error')
    with pytest.raises(ValueError, match='anchors stop at'):
        m.Sabs(P)


def test_error_policy_allows_the_anchored_range():
    m = NVModel(extrapolate='error')
    assert m.Sabs(120.0) == pytest.approx(4.61, rel=1e-12)


def test_unknown_policy_is_rejected():
    with pytest.raises(ValueError, match='extrapolate must be one of'):
        NVModel(extrapolate='extrapolate')


def test_boundary_pin_warns():
    """lambda_opt must not return a window edge silently."""
    m = build(zpl='linear', anchors='linear', hw_factor=HW_HARDENING[1])
    with pytest.warns(RuntimeWarning, match='pinned to the lower edge'):
        assert m.lambda_opt(400.0, lo=402.0, hi=640.0) == pytest.approx(402.0)


# ---------------------------------------------------------------------------
# 3. the spread is the result
# ---------------------------------------------------------------------------

def test_variants_agree_at_ambient():
    """Every variant is anchored at P = 0 by construction, so they must not
    differ there; a difference would mean a variant changed the anchor rather
    than the continuation."""
    values = [lambda_opt(build(**kw), 0.0) for _, kw, _ in eb.VARIANTS
              if 'hw_factor' not in kw]
    assert max(values) - min(values) == pytest.approx(0.0, abs=1e-9)


def test_linear_and_saturating_zpl_agree_at_the_anchor():
    a, b = NVModel(), LinearZPLModel()
    assert a.dZPL(ANCHOR_MAX_P) == pytest.approx(b.dZPL(ANCHOR_MAX_P), rel=1e-9)
    assert b.dZPL(240.0) == pytest.approx(2.0 * b.dZPL(120.0), rel=1e-12)


def test_spread_grows_with_pressure():
    widths = [w for _, _, _, w in spread()]
    assert widths == sorted(widths)
    assert all(a < b for a, b in zip(widths, widths[1:]))


def test_spread_magnitudes():
    """The numbers quoted in nv_model's C-8 note and in the strategy docs."""
    by_P = {P: w for P, _, _, w in spread()}
    assert by_P[120.0] == pytest.approx(21.4, abs=1.0)
    assert by_P[200.0] == pytest.approx(58.4, abs=2.0)
    assert by_P[400.0] == pytest.approx(133.8, abs=5.0)


def test_the_120_GPa_spread_is_not_zero():
    """hw has no pressure dependence in the model, so hardening it moves the
    answer INSIDE the anchored range.  This is a systematic on the published
    120 GPa number, not only on the extrapolated ones."""
    by_P = {P: w for P, _, _, w in spread()}
    assert by_P[120.0] > 15.0


def test_hw_hardening_moves_the_optimum_blue():
    """Larger hw pushes the sideband maximum further from the ZPL."""
    base = lambda_opt(build(), 120.0)
    for factor in HW_HARDENING:
        assert lambda_opt(build(hw_factor=factor), 120.0) < base


def test_baseline_flatlines_above_the_anchors():
    """The complaint, quantified: the shipped model's lambda_opt changes by
    only a few nm between 300 and 400 GPa because it is returning an
    asymptote."""
    a = lambda_opt(build(), 300.0)
    b = lambda_opt(build(), 400.0)
    assert abs(a - b) < 5.0


# ---------------------------------------------------------------------------
# 4. the one-sided bound
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('name,kw', list(eb.consistent_variants()))
@pytest.mark.parametrize('P', DEFAULT_PRESSURES)
def test_optimum_never_crosses_the_ionisation_wall(name, kw, P):
    """The bound, over the variants that continue ZPL and IP under one law."""
    m = build(**kw)
    assert lambda_opt(m, P) >= lambda_ion(m, P) - 1e-6


@pytest.mark.parametrize('name,kw', list(eb.mixed_variants()))
def test_mixed_continuation_is_flagged(name, kw):
    """Ho et al. Fig. 2(c): IP(3A2) closely follows the ZPL.  A variant that
    continues one and freezes the other is therefore not a physical
    alternative, and `mixed_variants` must keep flagging it -- silently
    treating these as candidates is what would let a bound be quoted that the
    model does not support."""
    assert (name, kw) not in eb.consistent_variants()


def test_freezing_IP_while_the_ZPL_moves_breaks_the_bound():
    """The concrete failure the flag exists for: the clip policy is safe
    applied to ALL anchors or to NONE, and unsafe applied selectively."""
    m = build(zpl='linear', anchors='clip')       # ZPL marches, IP frozen
    assert m.IP_A2(400.0) == m.IP_A2(ANCHOR_MAX_P)
    violation = lambda_ion(m, 400.0) - lambda_opt(m, 400.0)
    assert violation > 40.0                        # ~53 nm past the frozen wall

    consistent = build(zpl='linear', anchors='linear')
    assert lambda_opt(consistent, 400.0) >= lambda_ion(consistent, 400.0)


def test_wall_is_hc_over_IP():
    m = build(zpl='linear', anchors='linear')
    assert lambda_ion(m, 400.0) == pytest.approx(HBARC / m.IP_A2(400.0), rel=1e-12)


def test_margin_to_the_wall_closes_under_the_strongest_continuation():
    m = build(zpl='linear', anchors='linear', hw_factor=HW_HARDENING[1])
    margins = [abs(eb.wall_margin_eV(m, P)) for P in DEFAULT_PRESSURES]
    assert margins == sorted(margins, reverse=True)
    assert margins[-1] < 0.02          # pinned to the wall at 400 GPa


def test_pinned_optimum_is_insensitive_to_the_uncalibrated_a_gs():
    """a_gs has never been calibrated.  Three decades of it must not move the
    bound, which is why the bound is quotable and the value is not."""
    values = list(eb.a_gs_insensitivity().values())
    assert max(values) - min(values) < 1.0


def test_a_gs_still_matters_away_from_the_wall():
    """The insensitivity above is a property of being pinned, not of the model
    ignoring a_gs.  Deep in the blue, where the ReLU is active, it bites."""
    base = NVModel().a_gs
    lo = build(a_gs=0.01 * base, anchors='linear')
    hi = build(a_gs=100.0 * base, anchors='linear')
    E = nm2eV(360.0)
    f_lo, _ = lo.f_minus([(360.0, 1.0)], 200.0)
    f_hi, _ = hi.f_minus([(360.0, 1.0)], 200.0)
    assert float(f_lo) > float(f_hi)


# ---------------------------------------------------------------------------
# 5. the reconstruction still refuses, and must keep refusing
# ---------------------------------------------------------------------------

def test_published_kernel_still_refuses_above_the_anchors():
    from ho_spectrum_model import HoPublishedSpectrumModel
    k = HoPublishedSpectrumModel()
    with pytest.raises(ValueError, match='limited to 0--120 GPa'):
        k.sigma_abs(nm2eV(450.0), 300.0)


def test_report_helpers_run():
    assert len(eb.anchor_freeze_report()) == len(DEFAULT_PRESSURES)
    assert set(variant_table()) == {name for name, *_ in eb.VARIANTS}
    assert len(eb.consistent_variants()) + len(eb.mixed_variants()) == len(eb.VARIANTS)


# ---------------------------------------------------------------------------
# 6. the Lyapin (2018) curvature gate
# ---------------------------------------------------------------------------

def test_lyapin_coefficients_match_the_published_table():
    """Lyapin et al. (2018) Table 1, NV- row, verbatim."""
    assert eb.LYAPIN_2018[296.0]['alpha'] == pytest.approx(5.81e-3)
    assert eb.LYAPIN_2018[296.0]['beta'] == pytest.approx(-25e-6)
    assert eb.LYAPIN_2018[80.0]['alpha'] == pytest.approx(5.57e-3)
    assert eb.LYAPIN_2018[80.0]['beta'] == pytest.approx(-19e-6)
    assert eb.LYAPIN_P_MAX == pytest.approx(52.0)


def test_lyapin_alpha_reproduces_the_doherty_slope():
    """Their linear coefficient is an independent check on slope0 = 5.75
    meV/GPa, which nv_model takes from Doherty (2014).  Agreement to 3 %."""
    for T in (296.0, 80.0):
        assert eb.LYAPIN_2018[T]['alpha'] == pytest.approx(5.75e-3, rel=0.04)


def test_lyapin_quadratic_refuses_above_its_data():
    """It turns over at 116 GPa (296 K); evaluating it above 52 GPa would be
    reading a local expansion as an asymptotic law."""
    with pytest.raises(ValueError, match='turns over'):
        eb.lyapin_dZPL(120.0)
    assert eb.lyapin_dZPL(52.0) > eb.lyapin_dZPL(10.0)


def test_lyapin_quadratic_is_sublinear():
    """The whole point: beta < 0, so the shift falls below its own tangent."""
    for T in (296.0, 80.0):
        a = eb.LYAPIN_2018[T]['alpha']
        for P in (20.0, 40.0, 52.0):
            assert eb.lyapin_dZPL(P, T) < a * P


def test_gate_rejects_the_linear_zpl_and_keeps_the_saturating_one():
    gate = eb.lyapin_gate()
    for name, kw, _ in eb.VARIANTS:
        rms, passed = gate[name]
        assert passed == (kw.get('zpl') != 'linear'), name


def test_gate_separates_by_a_wide_margin():
    """6 meV against 40 meV: the rejection does not depend on where the
    threshold is put anywhere between them."""
    gate = eb.lyapin_gate()
    passing = [r for r, ok in gate.values() if ok]
    failing = [r for r, ok in gate.values() if not ok]
    assert max(passing) < 10.0
    assert min(failing) > 35.0
    assert min(failing) / max(passing) > 4.0


def test_gate_is_evaluated_at_lyapins_hydrostatic_condition():
    """Lyapin used a helium medium, alpha ~ 1, not the anchors' alpha = 0.95.
    Comparing at the wrong alpha flatters the model by ~3 meV, so the gate
    must fix alpha itself."""
    assert eb.LYAPIN_ALPHA == 1.0
    at_one = eb.zpl_rms_vs_lyapin(build(alpha=1.0))
    at_anchor = eb.zpl_rms_vs_lyapin(build(alpha=0.95))
    assert at_one > at_anchor
    assert eb.lyapin_gate()['baseline (code as shipped)'][0] == pytest.approx(
        at_one, rel=1e-9)


def test_gate_collapses_the_spread():
    before = {P: w for P, _, _, w in spread()}
    after = {P: w for P, _, _, w in spread(variants=eb.surviving_variants())}
    assert after[200.0] < 0.4 * before[200.0]
    assert after[400.0] < 0.3 * before[400.0]
    assert after[120.0] == pytest.approx(before[120.0], abs=1e-9)


def test_after_the_gate_hw_is_the_leading_uncertainty():
    """With the ZPL law settled, the residual spread must be reproducible from
    the hw variants alone -- which is the argument for measuring hw(P) next."""
    survivors = eb.surviving_variants()
    hw_only = tuple(v for v in survivors if 'hw_factor' in v[1]
                    or v[0].startswith('baseline'))
    full = {P: w for P, _, _, w in spread(variants=survivors)}
    hw = {P: w for P, _, _, w in spread(variants=hw_only)}
    assert hw[120.0] == pytest.approx(full[120.0], abs=1e-9)
    assert hw[200.0] > 0.5 * full[200.0]


# ---------------------------------------------------------------------------
# 7. the numbers quoted in docs/experiment/theory_limits_and_required_measurements.md
#    "120 GPa より上 -- 圧力外挿の限界".  Do not edit a number there without
#    editing the assertion here.
# ---------------------------------------------------------------------------

DOC_MINIMAX = {120.0: (466.5, 1.027), 200.0: (450.5, 1.027),
               300.0: (440.2, 1.047), 400.0: (432.2, 1.104)}
DOC_BAND5 = {120.0: (463.0, 488.0), 200.0: (448.0, 470.0),
             300.0: (440.0, 462.0), 400.0: (437.0, 459.0)}
DOC_WALL = {120.0: 405.2, 200.0: 374.2, 300.0: 341.6, 400.0: 314.1}
DOC_SPREAD_BEFORE = {120.0: 21.4, 200.0: 58.4, 300.0: 101.8, 400.0: 133.8}
DOC_SPREAD_AFTER = {120.0: 21.4, 200.0: 20.0, 300.0: 22.7, 400.0: 34.0}


@pytest.mark.parametrize('P', DEFAULT_PRESSURES)
def test_doc_minimax_wavelength(P):
    lam, pen = eb.minimax_wavelength(P)
    want_lam, want_pen = DOC_MINIMAX[P]
    assert lam == pytest.approx(want_lam, abs=0.3)
    assert pen == pytest.approx(want_pen, abs=0.002)


@pytest.mark.parametrize('P', DEFAULT_PRESSURES)
def test_doc_tolerance_band(P):
    lo, hi, width = eb.tolerance_band(P, 1.05)
    want_lo, want_hi = DOC_BAND5[P]
    assert lo == pytest.approx(want_lo, abs=1.0)
    assert hi == pytest.approx(want_hi, abs=1.0)
    assert 20.0 <= width <= 26.0        # "roughly constant in pressure"


def test_doc_tolerance_band_is_flat_in_pressure():
    """The argument that a 20-30 nm uncertainty costs only a few percent."""
    widths = [eb.tolerance_band(P, 1.05)[2] for P in DEFAULT_PRESSURES]
    assert max(widths) - min(widths) < 4.0


@pytest.mark.parametrize('P', DEFAULT_PRESSURES)
def test_doc_ionisation_wall(P):
    assert eb.lambda_ion(build(anchors='linear'), P) == pytest.approx(
        DOC_WALL[P], abs=0.2)


@pytest.mark.parametrize('P', DEFAULT_PRESSURES)
def test_doc_spreads(P):
    before = {p: w for p, _, _, w in spread()}
    after = {p: w for p, _, _, w in spread(variants=eb.surviving_variants())}
    assert before[P] == pytest.approx(DOC_SPREAD_BEFORE[P], abs=0.3)
    assert after[P] == pytest.approx(DOC_SPREAD_AFTER[P], abs=0.3)


def test_doc_model_discrepancy_at_the_anchor():
    """34.9 nm, and it exceeds the whole post-gate spread at 400 GPa -- the
    reason the document ranks resolving it first."""
    delta, single, kernel = eb.model_discrepancy_at_anchor()
    assert kernel == pytest.approx(440.65, abs=0.05)
    assert single == pytest.approx(475.5, abs=0.2)
    assert delta == pytest.approx(34.9, abs=0.2)
    after400 = {p: w for p, _, _, w in
                spread(variants=eb.surviving_variants())}[400.0]
    assert delta > after400


def test_doc_dzpl_saturation_fractions():
    m = NVModel()
    for P, frac in ((120.0, 0.702), (200.0, 0.867),
                    (300.0, 0.952), (400.0, 0.982)):
        assert float(m.dZPL(P) / m.Emax) == pytest.approx(frac, abs=0.002)


def test_doc_alpha_over_correction_flag():
    """Section 7, rank 4: the alpha = 0.95 model fits Lyapin's alpha ~ 1 data
    BETTER than the alpha = 1.0 model does, by ~3 meV."""
    assert eb.zpl_rms_vs_lyapin(build(alpha=0.95)) < eb.zpl_rms_vs_lyapin(
        build(alpha=1.0))
    from nv_model import _alpha_factor
    assert _alpha_factor(1.0) == pytest.approx(1.0559, abs=1e-3)
