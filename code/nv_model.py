"""
nv_model.py
-----------
Phenomenological model of NV-center CW-ODMR lock-in sensitivity under pressure,
used to compare green / blue / mixed excitation and to optimise the blue
excitation wavelength for high-pressure ODMR.

Physical anchors (all taken from):
  K. O. Ho, C. Dailledouze, V. Zalandauskas, et al.,
  "Optical Stability and Photophysics of NV Centers in Diamond up to 120 GPa",
  arXiv:2606.02399 (2026).

Anchored quantities
  * NV- ZPL: 1.945 eV at ambient, blue-shifting by ~0.40 eV at 120 GPa.
  * Absorption Huang-Rhys factor (JT-inactive modes): S_abs = 3.08 -> 4.61 (0..120 GPa).
  * Ground-state ionisation threshold  IP(3A2): 2.68 -> 3.06 eV (0..120 GPa).
  * Excited-state ionisation threshold IP(3E):  1.16 -> 1.63 eV (0..120 GPa).
  * Effective phonon energy hw = 65 meV.

Everything else (rate constants a_gs, a_es, r0, rbg; NV0 background weight w0;
linewidth model) is PHENOMENOLOGICAL and flagged as such; these are the knobs
to be pinned down by the (I_405, I_457) intensity-sweep calibration.

Sensitivity model
  Lock-in CW-ODMR sensitivity:      eta  ~  d(nu) / ( C * sqrt(R) )   (lower = better)
    C = C0 * f- / (f- + w0 (1-f-))     contrast, diluted by NV0 background
    R = f- * sigma_abs                 detected photon rate (fixed excitation power)
    f- = G_rec / (G_rec + G_ion)       steady-state NV- fraction
      G_ion = a_gs * ReLU(E_photon - IP_A2)   one-photon ground-state ionisation
            + a_es * sigma_abs                two-photon (via 3E) ionisation
      G_rec = r0 * sigma_abs + rbg           recombination (NV0 excitation) + background
"""

import numpy as np
from scipy.special import gammaln

HBARC = 1239.84                      # eV * nm
def nm2eV(lam_nm): return HBARC / np.asarray(lam_nm, float)
def eV2nm(E):      return HBARC / np.asarray(E, float)

E532 = nm2eV(532.0)                  # reference green photon energy
HW   = 0.065                         # effective phonon energy (eV)


class NVModel:
    """Phenomenological NV optical / charge / sensitivity model vs pressure."""

    def __init__(self,
                 Emax=0.758, P0=160.0,          # ZPL shift:  dE = Emax (1 - e^{-P/P0})
                 S_slope=(4.61 - 3.08),         # S_abs(P) = 3.08 + S_slope * P/120
                 a_gs=6.0, a_es=0.9,            # ionisation rate constants (phenom.)
                 r0=1.2, rbg=0.15,              # recombination + background (phenom.)
                 w0=1.0):                       # NV0 background brightness weight (phenom.)
        self.Emax, self.P0 = Emax, P0
        self.S_slope = S_slope
        self.a_gs, self.a_es = a_gs, a_es
        self.r0, self.rbg = r0, rbg
        self.w0 = w0
        self._norm = self._sigma_raw(E532, 0.0)   # normalise sigma to green@ambient = 1

    # ---- pressure-dependent physical quantities (anchored to Ho et al.) ----
    def ZPL(self, P):    return 1.945 + self.Emax * (1.0 - np.exp(-np.asarray(P, float) / self.P0))
    def Sabs(self, P):   return 3.08 + self.S_slope * np.clip(P, 0, 120) / 120.0
    def IP_A2(self, P):  return 2.68 + (3.06 - 2.68) * np.clip(P, 0, 120) / 120.0
    def IP_E(self, P):   return 1.16 + (1.63 - 1.16) * np.clip(P, 0, 120) / 120.0

    # ---- absorption cross section: low-T Franck-Condon (Pekarian) envelope ----
    def _sigma_raw(self, E, P):
        z = self.ZPL(P); S = self.Sabs(P)
        x = np.asarray(E, float) - z                 # detuning above ZPL
        n = np.clip(x / HW, 0, None)                  # phonon quantum number
        psb = np.exp(-S + n * np.log(S) - gammaln(n + 1.0)) * (x >= 0)   # phonon sideband
        dwf = 0.022 * np.exp(-(np.clip(P, 0, 120) / 120.0) * np.log(0.022 / 0.0036))
        zpl = dwf * np.exp(-x ** 2 / (2 * 0.015 ** 2))                    # narrow ZPL line
        return psb + zpl

    def sigma_abs(self, E, P):
        """Normalised absorption cross section (green@ambient == 1)."""
        return self._sigma_raw(E, P) / self._norm

    # ---- steady-state NV- fraction ----
    def f_minus(self, beams, P):
        """beams: list of (wavelength_nm, relative_intensity)."""
        relu = lambda u: np.clip(u, 0, None)
        s   = sum(I * self.sigma_abs(nm2eV(lam), P) for lam, I in beams)
        Ggs = sum(self.a_gs * relu(nm2eV(lam) - self.IP_A2(P)) * I for lam, I in beams)
        Gion = Ggs + self.a_es * s
        Grec = self.r0 * s + self.rbg
        return Grec / (Grec + Gion), s

    # ---- lock-in ODMR sensitivity (arb. units; lower = better) ----
    def linewidth(self, P):
        return 1.0 * (1 + 0.15 * np.clip(P, 0, 140) / 100.0)   # mild strain broadening (phenom.)

    def eta(self, beams, P):
        f, s = self.f_minus(beams, P)
        C = f / (f + self.w0 * (1 - f))
        R = np.clip(f * s, 1e-12, None)
        return self.linewidth(P) / (C * np.sqrt(R)), f, s, C, R

    # convenience: single-wavelength sweep
    def eta_lambda(self, lam_nm, P):
        return self.eta([(lam_nm, 1.0)], P)


def mc_band(build_fn, evaluate_fn, n=250, seed=0, lo=16, hi=84):
    """
    Monte-Carlo 16-84% uncertainty band.
      build_fn(rng)   -> returns an NVModel with randomised phenomenological knobs
      evaluate_fn(m)  -> returns a 1-D array for model m
    """
    rng = np.random.default_rng(seed)
    stack = [evaluate_fn(build_fn(rng)) for _ in range(n)]
    stack = np.array(stack)
    return np.percentile(stack, lo, 0), np.percentile(stack, hi, 0)


def default_randomiser(rng):
    """Standard randomisation of the phenomenological knobs for the MC bands."""
    return NVModel(
        Emax=0.758 * rng.uniform(0.90, 1.10),
        P0=160.0  * rng.uniform(0.85, 1.15),
        S_slope=(4.61 - 3.08) * rng.uniform(0.85, 1.15),
        a_gs=6.0 * rng.uniform(0.50, 1.60),
        r0=1.2   * rng.uniform(0.70, 1.30),
        rbg=0.15 * rng.uniform(0.60, 1.50),
    )
