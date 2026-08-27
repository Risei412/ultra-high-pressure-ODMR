"""
experiment_plan.py
------------------
Turning the frozen model into the two numbers a 120 GPa run actually needs:
which laser line, and how much power.

Both steps are written so that the experiment calibrates itself and does not
inherit assumptions it cannot check.

STEP 1 -- the excitation wavelength, from the measured ZPL
    lambda_opt is fixed by the absorption line shape, and the line shape is
    fixed by the ZPL position, the Huang-Rhys factor and the effective phonon
    energy.  Of those, the ZPL is the one that moves with anvil geometry, and
    it is also the one that is trivial to measure: it is a feature of the
    photoluminescence spectrum.

    This matters because the frozen anchors were measured on a micropillar
    (100) culet at alpha = 0.95, and a (111)-cut anvil selects the NV group
    aligned with the compression axis, whose 3E RED-shifts.  The anchors do not
    transfer.  Re-anchoring on the measured ZPL sidesteps the geometry problem
    entirely: whatever the stress state did, it did it to the ZPL, and the ZPL
    is what is read off.

    Residual assumption: S_abs(P) and hbar-omega are still taken from the
    near-hydrostatic reference.  A uniaxial stress that changed the
    electron-phonon coupling itself would not be caught this way.

STEP 2 -- the excitation power, from the shape of the saturation curve
    u = I/I_half is dimensionless and I_half is exactly what is not calibrated,
    so no absolute power can be quoted.  It does not have to be.  Both the
    onset of nonlinearity and the sensitivity optimum are at fixed values of u,
    so their RATIO is a model output that needs no calibration: measure the
    power at which the photoluminescence departs from linear, multiply by that
    ratio, and you are at the optimum without ever knowing I_half.

Run:  python experiment_plan.py            # the step 1 table, and step 2 on
                                           # synthetic data as a worked example
"""
import sys

import numpy as np
from scipy.optimize import brentq

from nv_model import NVModel, ALPHA_REF, eV2nm, nm2eV, _zpl_shape, _alpha_factor
from nv_model_power import NVModelPower

ZPL_AMBIENT = 1.945          # eV, the NV- zero-phonon line at ambient pressure
P_TARGET = 120.0             # GPa

# lines that exist as turn-key sources
COMMERCIAL = (405.0, 445.0, 457.0, 473.0, 488.0, 505.0, 532.0)


# ==========================================================================
# step 1 -- the excitation wavelength, measured directly at the working pressure
# ==========================================================================
SCAN_LINES = (445.0, 457.0, 473.0, 488.0, 505.0)


def lambda_opt_is_brightest(P=P_TARGET, u=None, lo=405.0, hi=560.0, n=621):
    """(argmax of the detected rate, argmin of eta) at this pressure.

    The two coincide, and that is the whole reason step 1 can be a measurement
    rather than a calculation: eta ~ dnu / (C sqrt(R)), and dnu and C carry no
    wavelength dependence, so minimising eta IS maximising the count rate.
    Whatever the absorption line shape turns out to be, THE BRIGHTEST LINE IS
    THE MOST SENSITIVE LINE -- provided the scan is taken below saturation,
    which is what `u` checks.
    """
    grid = np.linspace(lo, hi, n)
    if u is None:
        m = NVModel()
        R = np.asarray(m.eta_lambda(grid, P)[2])
        eta = np.asarray(m.eta_lambda(grid, P)[0])
    else:
        mp = NVModelPower()
        eta, _, R = (np.asarray(a) for a in mp.eta_lambda_u(grid, P, u)[:3])
    return float(grid[np.nanargmax(R)]), float(grid[np.nanargmin(eta)])


def excitation_scan(power_normalised_yield, lines=SCAN_LINES):
    """Reduce a measured excitation scan to the optimal wavelength.

    power_normalised_yield : detected PL, one value per line, each divided by
        the OPTICAL POWER actually delivered at that line (not the photon
        flux -- the lambda in R = f- (I/E_gamma) sigma_abs eta_col is already
        the thing being optimised, so dividing it out would remove the answer).

    Returns the peak of a log-quadratic through the points.  Three lines are
    the minimum; four or five bracketing the peak give a few nm.
    """
    lines = np.asarray(lines, float)
    y = np.asarray(power_normalised_yield, float)
    if len(y) != len(lines):
        raise ValueError('one yield per line')
    if len(y) < 3:
        raise ValueError('at least three lines are needed to locate a peak')
    c = np.polyfit(lines, np.log(np.clip(y, 1e-300, None)), 2)
    if c[0] >= 0:
        raise ValueError('the scan has no interior maximum -- the peak is '
                         'outside the lines used; add a bluer or redder line')
    return float(-c[1] / (2.0 * c[0]))


# ==========================================================================
# step 1b -- fallback: wavelength inferred from the measured ZPL
# ==========================================================================
def dE120_from_zpl(zpl_nm, P=P_TARGET, slope0=5.75e-3, alpha=ALPHA_REF):
    """Effective dE_ZPL(120 GPa) reproducing a ZPL measured at pressure P.

    'Effective' because it absorbs whatever the anvil geometry did: it is the
    number the frozen parameterisation needs in order to pass through the
    measured point, not a claim about the stress state.
    """
    shift = float(nm2eV(zpl_nm) - ZPL_AMBIENT)      # blue shift is positive
    if P == 120.0 and alpha == ALPHA_REF:
        return shift                                # dZPL(120) = dE120 exactly

    def residual(dE120):
        Emax, P0 = _zpl_shape(dE120, slope0)
        return _alpha_factor(alpha) * Emax * (1.0 - np.exp(-P / P0)) - shift

    # dE120 < 120*slope0 = 0.69 eV is the range in which the saturating form
    # has a solution at all; outside it the parameterisation degenerates
    lo, hi = 0.05, 120.0 * slope0 - 1e-4
    if residual(lo) > 0 or residual(hi) < 0:
        raise ValueError(
            f'a ZPL of {zpl_nm:.1f} nm at {P:.0f} GPa is outside the range this '
            f'parameterisation can represent (shift {shift*1e3:.0f} meV)')
    return brentq(residual, lo, hi)


def model_from_zpl(zpl_nm, P=P_TARGET, **kw):
    """The frozen model, re-anchored on a measured ZPL."""
    return NVModel(dE120=dE120_from_zpl(zpl_nm, P, **kw), **kw)


def lambda_opt_from_zpl(zpl_nm, P=P_TARGET, **kw):
    """Optimal excitation wavelength (nm) for a ZPL measured at pressure P."""
    return model_from_zpl(zpl_nm, P, **kw).lambda_opt(P)


def line_penalties(zpl_nm, P=P_TARGET, lines=COMMERCIAL, **kw):
    """[(wavelength, eta/eta_opt), ...] for each candidate line, best first."""
    m = model_from_zpl(zpl_nm, P, **kw)
    opt = m.lambda_opt(P)
    e_opt = float(np.asarray(m.eta_lambda(opt, P)[0]))
    out = [(lam, float(np.asarray(m.eta_lambda(lam, P)[0])) / e_opt)
           for lam in lines]
    return sorted(out, key=lambda r: r[1])


def zpl_range_for_line(lam, P=P_TARGET, tol=1.05, span=(0.15, 0.60), **kw):
    """ZPL window (nm) over which `lam` stays within `tol` of the optimum.

    This is the go/no-go for step 1: measure the ZPL, and if it falls inside
    the range returned here, the line is good enough and nothing needs redoing.
    """
    def penalty(dE):
        m = NVModel(dE120=float(dE), **kw)
        opt = m.lambda_opt(P)
        return (float(np.asarray(m.eta_lambda(lam, P)[0]))
                / float(np.asarray(m.eta_lambda(opt, P)[0])))

    grid = np.arange(span[0], span[1], 0.002)
    good = [d for d in grid if penalty(d) <= tol]
    if not good:
        return None
    lo, hi = min(good), max(good)
    # a larger ZPL shift is a shorter wavelength, so the ends swap over
    return (float(eV2nm(ZPL_AMBIENT + hi)), float(eV2nm(ZPL_AMBIENT + lo)))


# ==========================================================================
# step 2 -- power from the shape of the saturation curve
# ==========================================================================
def u_nonlinear(lam, P=P_TARGET, frac=0.10, model=None):
    """Normalised intensity at which R(u) falls `frac` below its linear slope.

    This is the feature the experiment can see without any calibration: the
    knee of the saturation curve.
    """
    mp = model or NVModelPower()
    R = lambda u: float(np.asarray(mp.eta_lambda_u(lam, P, u)[2]))
    u0 = 1e-4
    slope = R(u0) / u0
    return brentq(lambda u: R(u) / (slope * u) - (1.0 - frac), u0, 10.0)


def u_optimal(lam, P=P_TARGET, model=None):
    """Normalised intensity minimising eta at this wavelength and pressure."""
    mp = model or NVModelPower()
    u = np.logspace(-3, 0.7, 600)
    eta = np.array([float(np.asarray(mp.eta_lambda_u(lam, P, uu)[0]))
                    for uu in u])
    return float(u[int(eta.argmin())])


def power_working_point(lam=473.0, P=P_TARGET, frac=0.10):
    """The calibration-free recipe, as a dict.

    `ratio` is the whole point: the operating power is that many times the
    power at the knee, and neither number needs I_half to be known.
    """
    mp = NVModelPower()
    u_knee = u_nonlinear(lam, P, frac, mp)
    u_star = u_optimal(lam, P, mp)
    return dict(wavelength=lam, pressure=P, u_knee=u_knee, u_star=u_star,
                ratio=u_star / u_knee)


def analyse_power_sweep(power, pl, lam=473.0, P=P_TARGET, frac=0.10,
                        n_linear=3):
    """Reduce a measured saturation curve to a recommended operating power.

    power, pl : the measured excitation power and photoluminescence rate.
                Units are arbitrary and need not be absolute -- only the shape
                of the curve is used.
    n_linear  : how many of the lowest points define the linear slope.

    Returns the knee, the recommended operating power, and the model ratio
    that connects them.
    """
    power = np.asarray(power, float)
    pl = np.asarray(pl, float)
    order = np.argsort(power)
    power, pl = power[order], pl[order]
    if len(power) < n_linear + 2:
        raise ValueError('too few points to define a linear slope')

    slope = np.sum(pl[:n_linear] * power[:n_linear]) / np.sum(power[:n_linear] ** 2)
    ratio_curve = pl / (slope * power)

    below = np.flatnonzero(ratio_curve <= 1.0 - frac)
    if not len(below):
        knee = float('nan')                      # the sweep never reached it
    else:
        j = below[0]
        if j == 0:
            knee = float(power[0])
        else:                                    # linear interpolation in log P
            x0, x1 = np.log(power[j - 1]), np.log(power[j])
            y0, y1 = ratio_curve[j - 1], ratio_curve[j]
            t = (y0 - (1.0 - frac)) / (y0 - y1)
            knee = float(np.exp(x0 + t * (x1 - x0)))

    wp = power_working_point(lam, P, frac)
    return dict(knee_power=knee, operating_power=knee * wp['ratio'],
                linear_slope=float(slope), ratio_curve=ratio_curve, **wp)


# ==========================================================================
# report
# ==========================================================================
def _step1(P=P_TARGET):
    print(f'STEP 1  measured ZPL -> excitation wavelength   ({P:.0f} GPa)')
    print(f'  {"ZPL":>8}{"dE_ZPL":>10}{"lambda_opt":>12}   best commercial lines')
    for zpl in (565., 555., 545., 535., 529., 525., 518., 512.):
        dE = dE120_from_zpl(zpl, P)
        opt = lambda_opt_from_zpl(zpl, P)
        best = line_penalties(zpl, P)[:2]
        pretty = ', '.join(f'{l:.0f} nm x{p:.2f}' for l, p in best)
        print(f'  {zpl:6.0f} nm{dE*1e3:8.0f} meV{opt:10.1f} nm   {pretty}')

    print()
    for lam in (473.0, 488.0):
        rng = zpl_range_for_line(lam, P)
        print(f'  {lam:.0f} nm stays within 5% of the optimum for a ZPL of '
              f'{rng[0]:.0f}-{rng[1]:.0f} nm')
    print(f'  (the frozen micropillar anchor is a ZPL of '
          f'{eV2nm(ZPL_AMBIENT + 0.400):.0f} nm)')


def _step2(lam=473.0, P=P_TARGET):
    wp = power_working_point(lam, P)
    print(f'\nSTEP 2  saturation curve -> operating power   '
          f'({lam:.0f} nm, {P:.0f} GPa)')
    print(f'  knee of R(u), 10% below linear : u = {wp["u_knee"]:.3f}')
    print(f'  sensitivity optimum            : u = {wp["u_star"]:.3f}')
    print(f'  recipe: operate at {wp["ratio"]:.1f} x the knee power')

    mp = NVModelPower()
    e_star = float(np.asarray(mp.eta_lambda_u(lam, P, wp['u_star'])[0]))
    print(f'  {"u":>7}{"eta/eta*":>11}{"lambda_opt(u)":>16}')
    grid = np.linspace(402, 560, 700)
    for u in (0.03, 0.06, 0.10, 0.15, 0.20, 0.30, 0.50):
        e = float(np.asarray(mp.eta_lambda_u(lam, P, u)[0])) / e_star
        ridge = grid[int(np.nanargmin(
            np.asarray(mp.eta_lambda_u(grid, P, u)[0])))]
        print(f'  {u:7.2f}{e:11.3f}{ridge:14.0f} nm')

    # worked example on synthetic data, so the reduction can be exercised now
    truth = NVModelPower()
    I = np.logspace(-2, 0.3, 14)                 # arbitrary units
    R = np.array([float(np.asarray(truth.eta_lambda_u(lam, P, u)[2]))
                  for u in I]) * (1.0 + 0.0)     # noiseless
    got = analyse_power_sweep(I, R, lam, P)
    print(f'\n  worked example (synthetic sweep, I in units of I_half):')
    print(f'    knee found at I = {got["knee_power"]:.3f}  '
          f'(model says {wp["u_knee"]:.3f})')
    print(f'    operate at   I = {got["operating_power"]:.3f}  '
          f'(model says {wp["u_star"]:.3f})')


if __name__ == '__main__':
    P = float(sys.argv[1]) if len(sys.argv) > 1 else P_TARGET
    _step1(P)
    _step2(473.0, P)
