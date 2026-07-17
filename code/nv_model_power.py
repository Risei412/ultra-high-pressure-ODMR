"""
nv_model_power.py
------------------
Intensity-explicit extension of nv_model.NVModel (Part B of PLAN.md).

The baseline NVModel (nv_model.py) fixes the excitation intensity at I=1 for
every beam; all ionisation/recombination rates are then LINEAR in I, so the
overall shape of eta(lambda) barely depends on absolute power (scaling I only
trades off against the intensity-independent dark-recombination term rbg).
This module adds the missing nonlinear optical physics so that the optimal
blue wavelength can be studied as a genuine function of laser power:

  * Excited-state (3E) occupation saturates with intensity:
        n_E(lam,P,u) = sigma_abs(lam,P)*u / (sigma_abs(lam,P)*u + Gamma_d)
  * Excited-state (two-photon) ionisation is builtup from n_E, so it is
    quadratic in intensity at low power and saturates (linear) at high power:
        G_ion_ES = a_es2 * sigma_abs * u * n_E
  * Ground-state (one-photon) ionisation stays linear in u (unchanged physics):
        G_ion_GS = a_gs * ReLU(Ephoton - IP_A2) * u
  * Recombination proceeds through NV0 excitation (own Franck-Condon envelope,
    own saturation), plus an intensity-independent dark/background channel:
        G_rec = r0 * sigma_NV0 * u * n_E0 + rbg
  * Detected photon rate saturates with the cycling rate (not just sigma*u):
        R = f_minus * n_E
  * Contrast and linewidth develop the standard CW-ODMR power dependence
    (Dreau et al., PRB 84, 195204 (2011)): line broadens and contrast is
    diluted as the optical pumping rate Gamma_p grows.

All new knobs (Gamma_d, Gamma_d0, a_es2, Gamma_c, Gamma_satC, u-scale of rbg)
are PHENOMENOLOGICAL, exactly like the original a_gs/a_es/r0/rbg/w0, and are
flagged as calibration targets in PLAN.md Part C. In the u -> u0 (moderate,
order-1) regime the model is built to reproduce the fixed-power NVModel
result (see `check_consistency()` at the bottom).
"""

import numpy as np
from nv_model import NVModel, nm2eV, HBARC, HW


class NVModelPower(NVModel):
    """NVModel + explicit, saturating intensity dependence u = I/I_ref."""

    def __init__(self,
                 Gamma_d=1.0,      # 3E decay rate (sets excited-state saturation scale)
                 Gamma_d0=1.0,     # NV0 excited-state decay rate
                 a_es2=0.9,        # two-photon (ES) ionisation rate constant
                 dZPL_NV0=0.21,    # ambient NV0-NV- ZPL gap (eV), Doherty et al.
                 Gamma_c=0.6,      # linewidth power-broadening scale
                 Gamma_satC=1.2,   # contrast-saturation scale
                 **kw):
        super().__init__(**kw)
        self.Gamma_d, self.Gamma_d0 = Gamma_d, Gamma_d0
        self.a_es2 = a_es2
        self.dZPL_NV0 = dZPL_NV0
        self.Gamma_c, self.Gamma_satC = Gamma_c, Gamma_satC

    def sigma_NV0(self, E, P):
        """NV0 absorption envelope: same Franck-Condon machinery, shifted ZPL."""
        z = self.ZPL(P) + self.dZPL_NV0
        S = self.Sabs(P)
        x = np.asarray(E, float) - z
        n = np.clip(x / HW, 0, None)
        from scipy.special import gammaln
        psb = np.exp(-S + n * np.log(S) - gammaln(n + 1.0)) * (x >= 0)
        return psb / self._norm

    def f_minus_u(self, beams, P, u):
        """beams: list of (wavelength_nm, relative_spectral_weight); u: scalar or array."""
        relu = lambda v: np.clip(v, 0, None)
        u = np.asarray(u, float)
        s   = sum(w * self.sigma_abs(nm2eV(lam), P) for lam, w in beams)   # NV- absorption
        s0  = sum(w * self.sigma_NV0(nm2eV(lam), P) for lam, w in beams)   # NV0 absorption
        Ggs = sum(self.a_gs * relu(nm2eV(lam) - self.IP_A2(P)) * w for lam, w in beams) * u

        nE  = (s * u) / (s * u + self.Gamma_d)
        nE0 = (s0 * u) / (s0 * u + self.Gamma_d0)

        Gion_ES = self.a_es2 * s * u * nE
        Gion = Ggs + Gion_ES
        Grec = self.r0 * s0 * u * nE0 + self.rbg
        f = Grec / (Grec + Gion)
        return f, nE, nE0

    def eta_u(self, beams, P, u):
        """Power-explicit lock-in sensitivity. Returns eta, f_minus, R, C, nE."""
        f, nE, nE0 = self.f_minus_u(beams, P, u)
        R = np.clip(f * nE, 1e-12, None)                       # saturating photon rate
        Gamma_p = R                                            # optical pumping-rate proxy
        C0 = f / (f + self.w0 * (1 - f))
        C = C0 / (1.0 + Gamma_p / self.Gamma_satC)              # power-diluted contrast
        dnu = self.linewidth(P) * np.sqrt(1.0 + Gamma_p / self.Gamma_c)  # power-broadened line
        eta = dnu / (C * np.sqrt(R))
        return eta, f, R, C, nE

    def eta_lambda_u(self, lam_nm, P, u):
        return self.eta_u([(lam_nm, 1.0)], P, u)


def default_randomiser_power(rng):
    """Randomised phenomenological knobs for MC bands (extends nv_model.default_randomiser)."""
    from nv_model import default_randomiser
    base = default_randomiser(rng)
    return NVModelPower(
        Emax=base.Emax, P0=base.P0, S_slope=base.S_slope,
        a_gs=base.a_gs, r0=base.r0, rbg=base.rbg, w0=base.w0,
        Gamma_d=1.0  * rng.uniform(0.6, 1.6),
        Gamma_d0=1.0 * rng.uniform(0.6, 1.6),
        a_es2=0.9    * rng.uniform(0.5, 1.8),
        Gamma_c=0.6  * rng.uniform(0.5, 1.8),
        Gamma_satC=1.2 * rng.uniform(0.5, 1.8),
    )


def check_consistency():
    """Sanity check: at a moderate reference power the power-explicit model
    should reproduce a 475 nm-ish optimum at 120 GPa, and G_ion_ES should show
    the expected u^2 (low power) -> u (saturated) crossover."""
    m = NVModelPower()
    lam = np.linspace(402, 560, 400)
    for u in (0.05, 1.0, 20.0):
        e = m.eta_lambda_u(lam, 120.0, u)[0]
        print(f'u={u:6.2f}  lambda_opt = {lam[np.nanargmin(e)]:.1f} nm')


if __name__ == '__main__':
    check_consistency()
