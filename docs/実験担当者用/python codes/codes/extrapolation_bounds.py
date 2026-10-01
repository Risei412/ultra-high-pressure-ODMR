"""
extrapolation_bounds.py
-----------------------
What the published record does, and does not, say about the optimal excitation
wavelength above 120 GPa.

MOTIVATION.  Every pressure-dependent optical quantity in nv_model.NVModel is a
linear interpolation between the ambient and 120 GPa values of Ho et al.
(2026), and HoPublishedSpectrumModel -- the reconstruction that produces the
441 nm / 514.5 nm branch positions quoted in the manuscript -- refuses outright
above 120 GPa:

    ho_spectrum_model.py:  raise ValueError('published-curve interpolation is
                                             limited to 0--120 GPa')

NVModel does not refuse; under the default 'clip' policy it freezes each anchor
at its 120 GPa value and lets only dZPL keep moving.  dZPL is a saturating form
that is 95 % of the way to Emax by 300 GPa, so a call at 400 GPa returns an
asymptote.  That asymptote is not a prediction, and this module exists so that
it is never mistaken for one.

WHAT IS UNCONSTRAINED ABOVE 120 GPa.  Four independent choices, none of them
fixed by any published measurement:

  1. the ZPL law            saturating dE = Emax[1 - exp(-P/P0)]  vs  linear
  2. the anchor policy      S_abs, IP, DWF frozen  vs  continued linearly
  3. the phonon energy      hw is a CONSTANT 65 meV in the model, with no
                            pressure dependence at all, although compression
                            hardens the first-order diamond Raman line by
                            roughly 40-50 % over 0-400 GPa (the Akahama scale)
  4. the linewidth          nv_model.linewidth clips at 140 GPa; it is
                            wavelength independent, so it cannot move
                            lambda_opt, and it is not varied here

Choice 3 is not only an extrapolation problem: hardening hw moves lambda_opt at
120 GPa as well, by 21 nm, i.e. INSIDE the anchored range.  That is an
unevaluated systematic on the published 120 GPa numbers, not just on the
extrapolated ones.

WHAT SURVIVES.  One statement does, and it is one-sided.  The one-photon
ionisation threshold IP(3A2) acts as a hard blue wall: an excitation photon
above it ionises NV- out of the sensing charge state, and once the optimum is
pressed against the wall it stops depending on the ionisation rate constant.
`a_gs` is a phenomenological knob that has never been calibrated, and varying
it over three decades moves the pinned optimum by less than a nanometre.  So

    lambda_opt(P)  >=  lambda_ion(P) = hc / IP(3A2)(P)

is a bound that does not inherit the envelope's uncertainty.  It needs the
extrapolation of a SINGLE linear quantity rather than of the whole lineshape.

THE BOUND HAS A PRECONDITION, and finding it was worth the exercise.  It holds
only when the ZPL and IP(3A2) are continued under the SAME law.  Ho et al.
report that the ground-state ionisation energy closely follows the ZPL (their
Fig. 2(c); docs/references/notes/Ho2026.md), so a variant that marches the ZPL
blue while holding IP frozen at its 120 GPa value is not a physical
alternative -- it is an artifact of applying the clip policy to one of a pair
of coupled quantities.  That combination does break the inequality, by 53 nm at
400 GPa, and `ip_consistent` marks it.  The lesson generalises: the clip policy
is safe applied to ALL anchors or to NONE, and unsafe applied selectively.

WHAT THE LYAPIN GATE REMOVES.  Lyapin et al. (2018) measured the NV- ZPL in a
helium medium to 52 GPa and fitted E0 + alpha P + beta P^2 with beta < 0.  That
curvature rejects the linear-ZPL continuation -- 40 meV RMS against 6 meV for
the saturating one, a factor 6.5 -- WITHOUT any datum above 120 GPa, and their
alpha (5.81 meV/GPa at 296 K) independently reproduces the slope nv_model takes
from Doherty (2014).  The spread collapses accordingly: 58 -> 20 nm at
200 GPa, 134 -> 34 nm at 400 GPa.  What the gate does NOT do is select among
saturating laws, because their quadratic turns over at 116 GPa and is a local
expansion, not an asymptotic form.

After the gate the leading uncertainty is no longer the ZPL law.  It is hw,
which the model holds at its ambient 65 meV at every pressure, and which alone
accounts for the 21 nm that survives at 120 GPa.

WHAT THIS MODULE DOES NOT CLAIM.  It does not claim that any of the variants is
correct, that the linear continuation is more physical than the saturating one,
or that lambda_opt at 300-400 GPa is knowable from present data.  Its output is
a spread, and the spread is the result.

Provenance: nv_model.py (C-8), docs/references/notes/Ho2026.md,
docs/audits/positioning_in_high_pressure_sensing.md.
Tests: tests/test_extrapolation_bounds.py.

Run:  python extrapolation_bounds.py
"""

import warnings

import numpy as np

from nv_model import (ANCHOR_MAX_P, HBARC, HW, NVModel, nm2eV)

# The ZPL shift measured at the top anchor, Ho et al. (2026).  NVModel consumes
# this to solve for (Emax, P0) and does not retain it; the linear-ZPL variant
# needs it back.
DE120 = 0.400                       # eV

# Search window for the optimum.  Wider than the shipped [402, 640] nm, which
# was chosen for the anchored range and which the optimum leaves under the
# stronger continuations -- NVModel.lambda_opt warns when that happens.
LAM_LO, LAM_HI = 280.0, 700.0

# Pressures reported by the default table.  120 GPa is the anchor edge and is
# included so that the spread at the edge (which is not zero, because of hw)
# can be read next to the extrapolated ones.
DEFAULT_PRESSURES = (120.0, 200.0, 300.0, 400.0)

# hw multipliers bracketing the hardening of the first-order diamond Raman line
# over 0-400 GPa.  Unanchored: the model has no hw(P), so these are variants,
# not a correction.
HW_HARDENING = (1.30, 1.50)


class LinearZPLModel(NVModel):
    """NVModel with the saturating ZPL law replaced by a linear one.

    dE_ZPL(P) = (DE120 / 120 GPa) * P, carrying the same stress-anisotropy
    factor as the parent so that the two laws agree at the anchor by
    construction and differ only in how they continue.
    """

    def __init__(self, *args, dE120=DE120, **kwargs):
        self._dE120 = float(dE120)
        super().__init__(*args, **kwargs)

    def dZPL(self, P):
        P = np.clip(np.asarray(P, float), 0.0, None)
        return self._afac * self._dE120 * P / ANCHOR_MAX_P


def build(zpl='saturating', anchors='clip', hw_factor=1.0, **kwargs):
    """Construct one extrapolation variant.

    zpl        : 'saturating' (the fitted C-2 form) or 'linear'
    anchors    : NVModel `extrapolate` policy, 'clip' | 'linear' | 'error'
    hw_factor  : multiplier on the effective phonon energy (1.0 = the model's
                 ambient 65 meV, held at all pressures)
    """
    if zpl not in ('saturating', 'linear'):
        raise ValueError(f"zpl must be 'saturating' or 'linear', got {zpl!r}")
    cls = LinearZPLModel if zpl == 'linear' else NVModel
    return cls(hw=HW * float(hw_factor), extrapolate=anchors, **kwargs)


# The variants reported by default.  Each is one departure from the shipped
# model, then the three combined, so that a reader can attribute the spread.
#
# The third field is ip_consistent: whether the ZPL and IP(3A2) are continued
# under the same law.  The two mixed entries are kept in the spread -- they are
# legitimate statements about how much the POLICY CHOICE matters -- but they
# are not physical alternatives, and the ionisation bound is asserted only over
# the consistent ones.  See the module docstring.
VARIANTS = (
    ('baseline (code as shipped)', dict(), True),
    ('S_abs, IP, DWF linear',      dict(anchors='linear'), False),
    ('ZPL linear',                 dict(zpl='linear'), False),
    (f'hw x{HW_HARDENING[0]:.2f}', dict(hw_factor=HW_HARDENING[0]), True),
    (f'hw x{HW_HARDENING[1]:.2f}', dict(hw_factor=HW_HARDENING[1]), True),
    ('all linear',                 dict(zpl='linear', anchors='linear'), True),
    (f'all linear, hw x{HW_HARDENING[0]:.2f}',
     dict(zpl='linear', anchors='linear', hw_factor=HW_HARDENING[0]), True),
    (f'all linear, hw x{HW_HARDENING[1]:.2f}',
     dict(zpl='linear', anchors='linear', hw_factor=HW_HARDENING[1]), True),
)


# --- The Lyapin gate ------------------------------------------------------
# S. G. Lyapin, I. D. Ilichev, A. P. Novikov, V. A. Davydov, V. N. Agafonov,
# Nanosystems: Phys. Chem. Math. 9, 55 (2018), doi:10.17586/2220-8054-2018-9-1-55-57.
# PL of NV- in a DAC with a HELIUM pressure-transmitting medium, i.e. genuinely
# hydrostatic (alpha ~ 1), to 52 GPa, at 296 K and at 80 K.  Their Table 1 fits
#
#     E(P) = E0 + alpha P + beta P^2
#
# with a NEGATIVE quadratic coefficient for NV-.  That curvature is the point:
# it is an independent, room-temperature, hydrostatic measurement that the ZPL
# shift is measurably SUB-LINEAR well inside the anchored range, and it is what
# lets the linear-ZPL continuation be rejected without any data above 120 GPa.
#
# Their alpha, 5.81 meV/GPa at 296 K and 5.57 at 80 K, also independently
# reproduces the slope0 = 5.75 meV/GPa that nv_model takes from Doherty (2014).
#
# CAVEATS, both real.  (i) The quadratic turns over at -alpha/2beta = 116 GPa
# (296 K) and 147 GPa (80 K), so it is a local expansion and CANNOT be used as
# a law above ~50 GPa; it constrains curvature, not the asymptote.  (ii) It
# therefore rejects a linear continuation but does not select among saturating
# ones.
LYAPIN_2018 = {
    296.0: dict(E0=1.943, alpha=5.81e-3, beta=-25e-6),
    80.0:  dict(E0=1.946, alpha=5.57e-3, beta=-19e-6),
}
LYAPIN_P_MAX = 52.0      # GPa, the highest pressure they reached
LYAPIN_ALPHA = 1.0       # helium PTM: hydrostatic
LYAPIN_RMS_GATE = 20.0   # meV; see zpl_rms_vs_lyapin for how it separates


def lyapin_dZPL(P, temperature=296.0):
    """The measured NV- ZPL shift, Lyapin et al. (2018) Table 1.

    Raises above their maximum pressure rather than continuing a quadratic
    that turns over.
    """
    if temperature not in LYAPIN_2018:
        raise ValueError(f'Lyapin report 296 K and 80 K, not {temperature}')
    P = np.asarray(P, float)
    if np.any(P > LYAPIN_P_MAX):
        raise ValueError(
            f'Lyapin et al. reach {LYAPIN_P_MAX:g} GPa, and their quadratic '
            'turns over at 116 GPa (296 K) / 147 GPa (80 K); it constrains '
            'curvature, not the asymptote. Do not evaluate it above the data.')
    c = LYAPIN_2018[temperature]
    return c['alpha'] * P + c['beta'] * P ** 2


def zpl_rms_vs_lyapin(model, temperature=296.0, n=300):
    """RMS deviation of a model's dZPL from the Lyapin measurement, in meV.

    The comparison is made over 0-52 GPa.  A model built at the anvil geometry
    of the anchors (alpha = 0.95) is NOT at Lyapin's hydrostatic condition; to
    compare like with like, build it with alpha = LYAPIN_ALPHA.
    """
    P = np.linspace(0.0, LYAPIN_P_MAX, n)
    d = np.asarray(model.dZPL(P), float) - lyapin_dZPL(P, temperature)
    return 1e3 * float(np.sqrt(np.mean(d ** 2)))


def lyapin_gate(variants=None, temperature=296.0):
    """Which variants' ZPL law survives the Lyapin curvature measurement.

    Returns {name: (rms_meV, passed)}.  Evaluated at alpha = 1, Lyapin's
    condition, not at the anchors' alpha = 0.95.
    """
    variants = VARIANTS if variants is None else variants
    out = {}
    for name, kw, *_ in variants:
        kw = dict(kw, alpha=LYAPIN_ALPHA)
        rms = zpl_rms_vs_lyapin(build(**kw), temperature)
        out[name] = (rms, rms <= LYAPIN_RMS_GATE)
    return out


def surviving_variants(temperature=296.0):
    """The variants whose ZPL law is consistent with Lyapin et al. (2018)."""
    gate = lyapin_gate(temperature=temperature)
    return tuple((n, kw, ok) for n, kw, ok in VARIANTS if gate[n][1])


def consistent_variants():
    """The variants in which the ZPL and IP(3A2) follow the same law."""
    return tuple((name, kw) for name, kw, ok in VARIANTS if ok)


def mixed_variants():
    """The variants that continue one of a coupled pair and freeze the other."""
    return tuple((name, kw) for name, kw, ok in VARIANTS if not ok)


def lambda_opt(model, P, lo=LAM_LO, hi=LAM_HI):
    """Optimum over the widened window, without the boundary warning.

    The warning is suppressed because this function's callers report the
    boundary case through `pinned_at_wall` instead, which is the physically
    meaningful version of it.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return float(model.lambda_opt(P, lo=lo, hi=hi))


def lambda_ion(model, P):
    """The blue wall: the wavelength whose photon energy equals IP(3A2)(P)."""
    return float(HBARC / np.asarray(model.IP_A2(P), float))


def wall_margin_eV(model, P):
    """E(lambda_opt) - IP(3A2).  Negative means the optimum is below the wall."""
    return float(nm2eV(lambda_opt(model, P))) - float(np.asarray(model.IP_A2(P)))


def variant_table(pressures=DEFAULT_PRESSURES, variants=VARIANTS):
    """{variant name: [lambda_opt at each pressure]} for every variant."""
    return {name: [lambda_opt(build(**kw), P) for P in pressures]
            for name, kw, *_ in variants}


def spread(pressures=DEFAULT_PRESSURES, variants=VARIANTS):
    """Peak-to-peak spread of lambda_opt over the variants, per pressure.

    This is the headline: it is the width of what the published record leaves
    undetermined, in nanometres.
    """
    table = variant_table(pressures, variants)
    out = []
    for j, P in enumerate(pressures):
        col = [v[j] for v in table.values()]
        out.append((float(P), min(col), max(col), max(col) - min(col)))
    return out


def a_gs_insensitivity(P=400.0, factors=(0.1, 1.0, 10.0, 100.0), **kwargs):
    """lambda_opt against the uncalibrated ionisation rate constant a_gs.

    Evaluated in the strongest continuation, where the optimum is pressed
    against the wall.  A flat result is the point: once pinned, the bound no
    longer depends on the knob nobody has measured.
    """
    kwargs.setdefault('zpl', 'linear')
    kwargs.setdefault('anchors', 'linear')
    kwargs.setdefault('hw_factor', HW_HARDENING[1])
    base = NVModel().a_gs
    return {f * base: lambda_opt(build(a_gs=f * base, **kwargs), P)
            for f in factors}


def anchor_freeze_report(pressures=DEFAULT_PRESSURES):
    """How much of each anchored quantity is still moving above 120 GPa.

    Under the default 'clip' policy the answer is: only dZPL, and dZPL is
    nearly saturated.  Returned as a list of dicts, one per pressure.
    """
    m = NVModel()
    rows = []
    for P in pressures:
        rows.append(dict(
            P=float(P),
            dZPL=float(m.dZPL(P)),
            frac_of_Emax=float(m.dZPL(P) / m.Emax),
            Sabs=float(m.Sabs(P)),
            IP_A2=float(m.IP_A2(P)),
            linewidth=float(m.linewidth(P)),
        ))
    return rows


def minimax_wavelength(P, variants=None, lo=360.0, hi=620.0, step=0.25):
    """The excitation wavelength with the best worst-case penalty.

    For each surviving law, the penalty at lambda is eta(lambda) divided by
    that same law's own optimum, so a penalty of 1.00 means "optimal for this
    law".  The minimax choice is the wavelength whose LARGEST penalty across
    the laws is smallest, i.e. the wavelength that is safe whichever law turns
    out to hold.

    Returns (lambda_nm, worst_case_penalty).  This is a min-max over the
    variants considered, not a confidence interval: a law outside the family
    that still passes the Lyapin gate is not covered.
    """
    variants = surviving_variants() if variants is None else variants
    lam = np.arange(lo, hi, step)
    pen = []
    for _, kw, *_ in variants:
        e = np.asarray(build(**kw).eta_lambda(lam, P)[0], float)
        pen.append(e / e.min())
    worst = np.max(np.asarray(pen), axis=0)
    i = int(np.argmin(worst))
    return float(lam[i]), float(worst[i])


def tolerance_band(P, tol=1.05, model=None, lo=360.0, hi=620.0, step=0.25):
    """(lambda_lo, lambda_hi, width) of the band within `tol` of the best eta.

    The width is what makes the extrapolation spread tolerable: it is roughly
    constant in pressure, so a 20-30 nm uncertainty in the optimum costs only
    a few percent.
    """
    model = build() if model is None else model
    lam = np.arange(lo, hi, step)
    e = np.asarray(model.eta_lambda(lam, P)[0], float)
    inside = lam[e / e.min() <= tol]
    return float(inside.min()), float(inside.max()), float(np.ptp(inside))


def model_discrepancy_at_anchor():
    """Single-mode envelope minus reconstructed kernel, at 120 GPa, in nm.

    Both models are fully anchored here, so this is not an extrapolation
    error.  It is nonetheless larger than the whole post-gate extrapolation
    spread at 400 GPa, which is why resolving it comes first.
    """
    from ho_spectrum_model import HoPublishedSpectrumModel
    kernel = HoPublishedSpectrumModel().lambda_opt(ANCHOR_MAX_P)
    single = lambda_opt(build(), ANCHOR_MAX_P)
    return single - kernel, single, kernel


def main():
    pressures = DEFAULT_PRESSURES
    head = ' '.join(f'{P:>8.0f}' for P in pressures)

    print('Anchors stop at %.0f GPa. Under the default clip policy, what is '
          'still moving above it?' % ANCHOR_MAX_P)
    print(f"{'P (GPa)':>9s} {'dZPL':>8s} {'% Emax':>8s} {'S_abs':>8s} "
          f"{'IP(3A2)':>8s} {'linewidth':>10s}")
    for r in anchor_freeze_report(pressures):
        print(f"{r['P']:9.0f} {r['dZPL']:8.4f} {100 * r['frac_of_Emax']:7.1f}% "
              f"{r['Sabs']:8.3f} {r['IP_A2']:8.3f} {r['linewidth']:10.3f}")

    print(f'\nlambda_opt (nm) over the search window [{LAM_LO:.0f}, '
          f'{LAM_HI:.0f}] nm\n')
    print(f"{'variant':30s} {head}")
    print('-' * (30 + len(head) + 1))
    for name, values in variant_table(pressures).items():
        print(f'{name:30s} ' + ' '.join(f'{v:8.1f}' for v in values))

    print('-' * (30 + len(head) + 1))
    for P, lo, hi, width in spread(pressures):
        print(f'P = {P:3.0f} GPa   spread = {width:6.1f} nm   '
              f'[{lo:.1f} .. {hi:.1f}]')

    print('\nThe Lyapin (2018) curvature gate, at their hydrostatic condition')
    print(f'  RMS over 0-{LYAPIN_P_MAX:.0f} GPa, gate at {LYAPIN_RMS_GATE:.0f} meV')
    for name, (rms, ok) in lyapin_gate().items():
        print(f"  {'PASS' if ok else 'FAIL'}  {rms:6.1f} meV   {name}")

    survivors = surviving_variants()
    print('\nSpread after the gate (nm)')
    print(f"{'P (GPa)':>9s} {'before':>9s} {'after':>9s}   {'range after':>22s}")
    for (P, _, _, w0), (_, lo, hi, w) in zip(spread(pressures),
                                             spread(pressures, survivors)):
        print(f'{P:9.0f} {w0:9.1f} {w:9.1f}   [{lo:8.1f} .. {hi:6.1f}]')
    print('  The residual is now dominated by hw, which the model holds at its')
    print('  AMBIENT value at every pressure.  That is the next thing to measure.')

    print('\nThe one-sided bound: the ionisation wall lambda_ion = hc/IP(3A2)')
    strongest = build(zpl='linear', anchors='linear',
                      hw_factor=HW_HARDENING[1])
    print(f"{'P (GPa)':>9s} {'lambda_opt':>11s} {'lambda_ion':>11s} "
          f"{'E-IP (eV)':>10s}")
    for P in pressures:
        print(f'{P:9.0f} {lambda_opt(strongest, P):11.1f} '
              f'{lambda_ion(strongest, P):11.1f} '
              f'{wall_margin_eV(strongest, P):+10.3f}')

    print('\nlambda_opt(400 GPa) against the uncalibrated ionisation constant '
          'a_gs, at the wall')
    for value, lam in a_gs_insensitivity().items():
        print(f'  a_gs = {value:8.2f}   lambda_opt = {lam:6.1f} nm')


if __name__ == '__main__':
    main()

