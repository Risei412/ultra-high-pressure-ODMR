"""
nv_bg.py
--------
Background-aware extension of nv_model.NVModel (Part E of PLAN.md).

WHY
    The baseline model scores an excitation choice with

        eta_0  ~  dnu / (C sqrt(R))

    where R counts NV photons only.  A DAC experiment also collects
    spin-independent light inside the same detection window (ruby R line,
    anvil luminescence, ...).  That light does TWO things, and getting only
    one of them right flips the sign of the effect:

      1. it adds shot noise            ->  noise ~ sqrt(R + B)
      2. it dilutes the contrast       ->  C_eff = C * R/(R + B)

    Putting both in:

        eta = dnu / (C_eff sqrt(R+B))
            = dnu sqrt(R+B) / (C R)
            = eta_0 * sqrt(1 + rho),        rho = B/R

    So background ALWAYS costs sensitivity (monotonic in rho), and the
    excitation-wavelength dependence enters through rho(lambda) = B/R, i.e.
    through the COMPETITION between the background shape g(lambda)
    (background.py) and the NV photon rate R(lambda) (nv_model.py).

    Note that a model which only wrote sqrt(R+B) in the denominator would
    predict eta_0/sqrt(1+rho) -- background IMPROVING the sensitivity, which
    is obviously wrong.  The contrast dilution is what fixes it.

CONVENTION
    rho0 is the background-to-signal ratio for 532 nm excitation at ambient
    pressure.  It is the single dial that carries the (experimentally
    unknown) absolute background level, and the analysis is reported as a
    function of it.
"""

import numpy as np

from nv_model import NVModel, default_randomiser
import background as bg


class NVModelBG(NVModel):
    """NVModel + spin-independent background inside the detection window."""

    def __init__(self, mix=None, c_P=0.0, **kw):
        super().__init__(**kw)
        self.mix = bg.default_mix() if mix is None else dict(mix)
        self.c_P = c_P                                   # pressure growth of `broad`
        self._R_ref = float(self.eta([(bg.LAM_REF, 1.0)], 0.0)[4])   # green @ ambient

    # ------------------------------------------------------------------
    def background(self, beams, P, rho0):
        """Background photon rate B in the same units as R.

        Beams are added incoherently, weighted by their relative intensity,
        exactly as `f_minus` adds their absorption.
        """
        total_I = sum(I for _, I in beams)
        if total_I == 0:
            return np.zeros_like(np.asarray(P, float))
        shape = sum(I * bg.g_total(lam, P, self.mix, self.c_P) for lam, I in beams)
        return rho0 * self._R_ref * shape

    def eta_bg(self, beams, P, rho0):
        """Background-aware sensitivity.

        Returns (eta, f_minus, R, C, B, rho).  rho0 = 0 reproduces
        NVModel.eta exactly.
        """
        eta0, f, s, C, R = self.eta(beams, P)
        B = self.background(beams, P, rho0)
        rho = B / R
        return eta0 * np.sqrt(1.0 + rho), f, R, C, B, rho

    def eta_lambda_bg(self, lam_nm, P, rho0):
        return self.eta_bg([(lam_nm, 1.0)], P, rho0)


# --------------------------------------------------------------- helpers ---
def optimum_wavelength(model, P, rho0, lam=None, refine=True):
    """argmin_lambda eta(lambda) with parabolic refinement around the grid min."""
    lam = np.linspace(402, 560, 1580) if lam is None else lam
    e = model.eta_lambda_bg(lam, P, rho0)[0]
    i = int(np.nanargmin(e))
    if not refine or i in (0, len(lam) - 1):
        return float(lam[i])
    y0, y1, y2 = e[i - 1], e[i], e[i + 1]
    denom = y0 - 2 * y1 + y2
    if denom <= 0:
        return float(lam[i])
    d = 0.5 * (y0 - y2) / denom
    return float(lam[i] + d * (lam[i + 1] - lam[i]))


def tolerance_band(model, P, rho0, lam=None, tol=1.05):
    """Wavelength interval over which eta stays within `tol` of the optimum."""
    lam = np.linspace(402, 560, 1580) if lam is None else lam
    e = model.eta_lambda_bg(lam, P, rho0)[0]
    ok = e <= tol * np.nanmin(e)
    idx = np.flatnonzero(ok)
    return float(lam[idx[0]]), float(lam[idx[-1]])


def crossover_pressure(model, rho0, blue=457.0, green=532.0, P=None):
    """Pressure at which blue overtakes green.  NaN if it never does."""
    P = np.linspace(0, 200, 2001) if P is None else P
    eB = model.eta_bg([(blue, 1.0)], P, rho0)[0]
    eG = model.eta_bg([(green, 1.0)], P, rho0)[0]
    better = eB < eG
    if not better.any():
        return float('nan')
    return float(P[int(np.argmax(better))])


def randomiser_bg(rng):
    """MC randomiser: baseline knobs + a random background composition.

    The three channel weights are drawn uniformly from the 2-simplex
    (Dirichlet(1,1,1)), so the Monte-Carlo band answers "does the conclusion
    survive ANY plausible background composition?" rather than assuming the
    equal-weight mix.
    """
    base = default_randomiser(rng)
    w = rng.dirichlet(np.ones(len(bg.CHANNELS)))
    return NVModelBG(
        mix=dict(zip(bg.CHANNELS, w)),
        c_P=rng.uniform(0.0, 1.5),
        Emax=base.Emax, P0=base.P0, S_slope=base.S_slope,
        a_gs=base.a_gs, r0=base.r0, rbg=base.rbg, w0=base.w0,
    )


if __name__ == '__main__':
    m = NVModelBG()
    print('rho0   lambda_opt@120GPa   crossover(457 vs 532)')
    for rho0 in (0.0, 0.1, 1.0, 10.0, 100.0):
        print(f'{rho0:6.2f} {optimum_wavelength(m, 120.0, rho0):12.1f} nm '
              f'{crossover_pressure(m, rho0):14.1f} GPa')
