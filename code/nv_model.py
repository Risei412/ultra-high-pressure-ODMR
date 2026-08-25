"""
nv_model.py
-----------
Phenomenological model of NV-center CW-ODMR lock-in sensitivity under pressure,
used to compare green / blue / mixed excitation and to optimise the blue
excitation wavelength for high-pressure ODMR.

THEORY STATUS: FROZEN.  The three structural revisions C-1..C-3 below are in.
Any further change must keep `tests/test_freeze.py` green.

Physical anchors
  Optical / charge (Ho et al., arXiv:2606.02399, 2026 -- DFT + DAC to 120 GPa):
    * NV- ZPL shift            dE_ZPL(120 GPa) = 0.400 eV  (>400 meV reported)
    * Absorption Huang-Rhys    S_abs = 3.08 -> 4.61 (0..120 GPa, JT-inactive modes)
    * Ground-state ionisation  IP(3A2): 2.68 -> 3.06 eV
    * Excited-state ionisation IP(3E):  1.16 -> 1.63 eV
    * Absorption DWF           2.2% -> 0.36%
    * Effective phonon energy  hw = 65 meV  (Kehayias et al.)
  ZPL low-pressure slope (Doherty et al., PRL 112, 047601 (2014), 0..60 GPa):
    * dE_ZPL/dP|_{P->0} = 5.75 meV/GPa   (5.5 meV/GPa measured)
  ODMR contrast vs pressure (Dai et al., Chin. Phys. Lett. 39, 117601 (2022),
  532 nm, room temperature, microdiamonds to 146 GPa):
    * C = 14% at ambient -> ~1% plateau above ~50 GPa

Everything else (rate constants a_gs, a_es, r0, rbg; NV0 background weight w0;
linewidth model) is PHENOMENOLOGICAL and flagged as such; these are the knobs
to be pinned down by the (I_405, I_457) intensity-sweep calibration.

Sensitivity model
  Lock-in CW-ODMR sensitivity:      eta  ~  d(nu) / ( C * sqrt(R) )   (lower = better)
    C = C0(P) * f- / (f- + w0 (1-f-))  contrast: ISC prefactor x NV0 dilution
    R = f- * Phi * sigma_abs           detected rate at fixed optical power
      Phi = I / E_photon               photon flux (I denotes optical power density)
    f- = G_rec / (G_rec + G_ion)       steady-state NV- fraction
      G_ion = a_gs * ReLU(E_photon - IP_A2)   one-photon ground-state ionisation
            + a_es * sigma_abs                two-photon (via 3E) ionisation
      G_rec = r0 * sigma_abs + rbg           recombination (NV0 excitation) + background

---------------------------------------------------------------------------
C-1  FINITE-TEMPERATURE ABSORPTION ENVELOPE
---------------------------------------------------------------------------
The previous envelope was the T = 0 Pekarian, exp(-S) S^p / Gamma(p+1) with a
hard cut at p = (E - E_ZPL)/hw >= 0, i.e. *identically zero* absorption below
the ZPL.  At 120 GPa the 532 nm line sits 14 meV BELOW the ZPL (the crossing is
at 113.7 GPa), so that cut, not the photophysics, was setting sigma_abs(532).
It produced a 4.6x discontinuity in sigma_abs(532 nm) at 113.7 GPa and predicted
eta(532)/eta_opt = 110 at 140 GPa -- i.e. "green is dead" -- which is refuted by
Dai et al. (532 nm cwODMR to 140 GPa at ROOM temperature).

The fix is the standard single-mode multiphonon lineshape at temperature T
(Huang & Rhys 1950; Struck & Fonger 1975), with nbar = 1/(exp(hw/kT) - 1):

    L(p) = exp[-S(2 nbar + 1)] * ((nbar+1)/nbar)^(p/2) * I_{|p|}( 2 S sqrt(nbar(nbar+1)) )

continued to real p.  It is nonzero for p < 0 (anti-Stokes / hot-band
absorption) and reduces EXACTLY to exp(-S) S^p / Gamma(p+1) for p >= 0 and to 0
for p < 0 as T -> 0.  T is a measured experimental condition, not a fitted knob.

---------------------------------------------------------------------------
C-2  ZPL PARAMETERISED BY THE MEASURED QUANTITIES
---------------------------------------------------------------------------
The saturating form dE(P) = Emax [1 - exp(-P/P0)] is retained, but (Emax, P0)
are no longer independent inputs.  They are DERIVED from the two quantities that
were actually measured:

    dE120  = dE_ZPL(120 GPa)        = 0.400 eV   (Ho et al.)
    slope0 = dE_ZPL/dP |_{P->0}     = 5.75 meV/GPa  (Doherty et al., 0..60 GPa)

The old defaults (Emax = 0.758 eV, P0 = 160 GPa) reproduce dE120 = 0.400 eV but
give slope0 = 4.74 meV/GPa, ~20% below the measured low-pressure slope, and --
more importantly -- randomising Emax and P0 independently let dE_ZPL(120 GPa)
wander far outside its measured uncertainty, which is what inflated the Monte
Carlo band on lambda_opt.  Randomising (dE120, slope0) instead propagates
measurement uncertainty rather than parameterisation slack.

---------------------------------------------------------------------------
C-3  CONTRAST: EXPLICIT ISC PREFACTOR
---------------------------------------------------------------------------
Contrast at high pressure is set by the intersystem crossing, not by NV0
dilution: the 3E -> 1A1 rate follows the phonon density of states resonant with
the 3E-1A1 gap, which opens under hydrostatic compression (Goldman et al., PRB
91, 165201 (2015); Bhattacharyya, PhD thesis, UC Berkeley 2022).  Contrast is
therefore lost under hydrostatic / [100] stress and RETAINED for the [111] NV
group of a [111]-cut culet, whose 3E red-shifts (Davies & Hamer 1976).

C0(P) is added as an explicit, WAVELENGTH-INDEPENDENT prefactor calibrated to
the Dai et al. contrast data.  Being wavelength independent it cannot move
lambda_opt -- test_freeze.py asserts this -- but it makes the assumption visible
and lets absolute contrast be quoted.  The condition under which it would fail
(orbital-branch-selective excitation making the ISC rate wavelength dependent)
is a stated limitation, not something this model covers.
"""

import numpy as np
from scipy.special import gammaln, ive
from scipy.optimize import brentq

HBARC = 1239.84                      # eV * nm
KB    = 8.617333262e-5               # eV / K
def nm2eV(lam_nm): return HBARC / np.asarray(lam_nm, float)
def eV2nm(E):      return HBARC / np.asarray(E, float)

E532 = nm2eV(532.0)                  # reference green photon energy
HW   = 0.065                         # effective phonon energy (eV)

# Temperature below which the T=0 Pekarian is used verbatim (nbar < 1e-8).
_T0_TOL = 1e-8
_COL_GRID = 1201      # energy-grid points for the collection-efficiency integral

# --- C-4: stress anisotropy of the ZPL shift -------------------------------
# At an anvil tip the normal stress sigma_perp equals the chamber pressure P,
# while the tangential component is reduced, sigma_par = alpha * P.  Hilberer
# et al., arXiv:2301.05094 / PRB 107, L220102 (2023) MEASURED alpha from the
# ODMR field response:
#     alpha = 0.56  standard flat culet   (theory estimate 0.6)
#     alpha = 0.95  FIB-machined micropillar (quasi-hydrostatic)
# and measured the ZPL shift per unit compressed volume in the two geometries:
#     -769 +/- 4 meV/(cm^3 mol^-1)  micropillar
#     -434 +/- 2 meV/(cm^3 mol^-1)  standard anvil
# i.e. at equal compression the deviatorically stressed culet shifts the ZPL by
# only 434/769 = 0.564 of the quasi-hydrostatic value: deviatoric stress RED
# shifts the 3E, opposing the hydrostatic blue shift (Davies & Hamer 1976).
#
# Writing dE_ZPL ~ c1 * sigma_mean + c2 * sigma_dev with
#     sigma_mean(alpha) = (1 + 2 alpha)/3,   sigma_dev(alpha) = 1 - alpha
# the two measurements fix c2/c1; ALPHA_REF is the geometry the anchors were
# taken in (Ho et al. also quote alpha ~ 0.95), so g(ALPHA_REF) = 1 by
# construction and the frozen 120 GPa answer is unchanged.
ALPHA_REF = 0.95
ZPL_DEV_K = -0.39211        # c2/c1, from the 434/769 ratio at alpha = 0.56


def _alpha_factor(alpha):
    """Multiplier on the hydrostatic ZPL shift for stress anisotropy alpha."""
    f = lambda a: (1.0 + 2.0 * a) / 3.0 + ZPL_DEV_K * (1.0 - a)
    return f(alpha) / f(ALPHA_REF)


def _zpl_shape(dE120, slope0):
    """
    Solve  dE(P) = Emax [1 - exp(-P/P0)]  for (Emax, P0) from the two measured
    quantities dE120 = dE(120 GPa) and slope0 = dE'(0).

    With x = 120/P0 the constraint is  (1 - e^-x)/x = dE120 / (120 * slope0).
    The RHS is in (0,1) whenever the curve saturates; RHS -> 1 is the linear
    limit, where P0 -> inf and the exponential degenerates.
    """
    rhs = dE120 / (120.0 * slope0)
    if not (0.0 < rhs < 1.0):
        # No saturating solution (slope0 too small for dE120): fall back to the
        # linear law dE = slope0 * P, represented by a very large P0.
        P0 = 1e6
        return slope0 * P0, P0
    g = lambda x: (1.0 - np.exp(-x)) / x - rhs
    x = brentq(g, 1e-9, 700.0)
    P0 = 120.0 / x
    return slope0 * P0, P0


class NVModel:
    """Phenomenological NV optical / charge / sensitivity model vs pressure."""

    def __init__(self,
                 dE120=0.400, slope0=5.75e-3,   # ZPL, from measurement (C-2)
                 Emax=None, P0=None,            # legacy override of the above
                 T=300.0,                       # measurement temperature, K (C-1)
                 hw=HW,                         # effective phonon energy (eV)
                 intensity_basis='optical_power', # beam weights are equal optical powers
                 alpha=ALPHA_REF,               # stress anisotropy sigma_par/sigma_perp (C-4)
                 S_slope=(4.61 - 3.08),         # S_abs(P) = 3.08 + S_slope * P/120
                 a_gs=6.0, a_es=0.9,            # ionisation rate constants (phenom.)
                 r0=1.2, rbg=0.15,              # recombination + background (phenom.)
                 w0=1.0,                        # NV0 background brightness weight (phenom.)
                 zpl_width=0.015,               # ZPL gaussian width (eV), phenom.
                 C_amb=0.2469, C_floor=0.0,     # ISC contrast prefactor (C-3),
                 E_isc=0.1807, isc=True,        #   calibrated to Dai et al. 2022
                 det_band=(650.0, 800.0),       # detection passband, nm (C-7)
                 collection=True):              # apply the collection efficiency
        if Emax is not None or P0 is not None:
            # Legacy path: (Emax, P0) given directly.
            self.Emax = 0.758 if Emax is None else Emax
            self.P0   = 160.0 if P0   is None else P0
        else:
            self.Emax, self.P0 = _zpl_shape(dE120, slope0)
        self.T = float(T)
        self.hw = float(hw)
        if intensity_basis not in ('optical_power', 'photon_flux'):
            raise ValueError("intensity_basis must be 'optical_power' or 'photon_flux'")
        self.intensity_basis = intensity_basis
        self.alpha = float(alpha)
        self._afac = _alpha_factor(self.alpha)
        self.S_slope = S_slope
        self.a_gs, self.a_es = a_gs, a_es
        self.r0, self.rbg = r0, rbg
        self.w0 = w0
        self.zpl_width = zpl_width
        self.C_amb, self.C_floor, self.E_isc, self.isc = C_amb, C_floor, E_isc, isc
        self.det_band, self.collection = det_band, collection
        self._col_cache = {}
        self._norm = self._sigma_raw(E532, 0.0)   # normalise sigma to green@ambient = 1

    # ---- pressure-dependent physical quantities (anchored to Ho et al.) ----
    def dZPL(self, P):
        """ZPL blue shift at pressure P for this anvil geometry (C-2, C-4)."""
        return self._afac * self.Emax * (1.0 - np.exp(-np.asarray(P, float) / self.P0))
    def ZPL(self, P):    return 1.945 + self.dZPL(P)
    def Sabs(self, P):   return 3.08 + self.S_slope * np.clip(P, 0, 120) / 120.0
    def Sem(self, P):    return 3.39 + (5.25 - 3.39) * np.clip(P, 0, 120) / 120.0
    def IP_A2(self, P):  return 2.68 + (3.06 - 2.68) * np.clip(P, 0, 120) / 120.0
    def IP_E(self, P):   return 1.16 + (1.63 - 1.16) * np.clip(P, 0, 120) / 120.0

    # ---- finite-temperature Franck-Condon envelope (C-1) ----
    def _fc(self, x, S):
        """
        Single-mode multiphonon absorption envelope at temperature self.T.

        x : E - E_ZPL (eV), may be negative (anti-Stokes side)
        S : Huang-Rhys factor

        T -> 0 limit is exactly exp(-S) S^p / Gamma(p+1) for p >= 0, else 0.
        """
        x = np.asarray(x, float)
        p = x / self.hw
        S = np.asarray(S, float)

        if self.T <= 0.0:
            nbar = 0.0
        else:
            r = self.hw / (KB * self.T)
            nbar = 1.0 / np.expm1(r) if r < 700.0 else 0.0

        if nbar < _T0_TOL:                                  # frozen T=0 branch
            n = np.clip(p, 0, None)
            return np.exp(-S + n * np.log(S) - gammaln(n + 1.0)) * (x >= 0)

        z = 2.0 * S * np.sqrt(nbar * (nbar + 1.0))
        v = np.abs(p)
        # Work with the EXPONENTIALLY SCALED Bessel ive(v,z) = I_v(z) exp(-z),
        # whose exp(-z) is exactly what the prefactor below supplies.  Where ive
        # underflows, fall back to the small-argument series
        # log I_v(z) ~ v log(z/2) - lgamma(v+1), converted to the scaled form.
        iv = ive(v, z)
        with np.errstate(divide='ignore'):
            log_ive = np.where(iv > 0.0, np.log(np.where(iv > 0.0, iv, 1.0)),
                               v * np.log(z / 2.0) - gammaln(v + 1.0) - z)
        # exp(-S(2nbar+1)) * I_v(z) = exp(-S (sqrt(nbar+1) - sqrt(nbar))^2) * ive(v,z)
        log_pref = -S * (np.sqrt(nbar + 1.0) - np.sqrt(nbar)) ** 2
        log_therm = x / (2.0 * KB * self.T)                 # = p * log((nbar+1)/nbar)/2
        return np.exp(log_pref + log_therm + log_ive)

    # ---- absorption cross section ----
    def _sigma_raw(self, E, P, zpl=None):
        z = self.ZPL(P) if zpl is None else zpl
        S = self.Sabs(P)
        x = np.asarray(E, float) - z                 # detuning above ZPL
        psb = self._fc(x, S)                          # phonon sideband
        dwf = 0.022 * np.exp(-(np.clip(P, 0, 120) / 120.0) * np.log(0.022 / 0.0036))
        zpl_line = dwf * np.exp(-x ** 2 / (2 * self.zpl_width ** 2))
        return psb + zpl_line

    def sigma_abs(self, E, P):
        """Normalised absorption cross section (green@ambient == 1)."""
        return self._sigma_raw(E, P) / self._norm

    # ---- emission envelope and collection efficiency (C-7) ----
    def sigma_em(self, E, P):
        """
        NV- emission lineshape, the mirror image of the absorption envelope:
        emission is displaced BELOW the ZPL, and carries its own Huang-Rhys
        factor S_em = 3.39 -> 5.25 and DWF_em = 4.9% -> <1% (Ho et al. 2026).
        Not normalised; only ratios are used.
        """
        z = self.ZPL(P)
        x = z - np.asarray(E, float)                 # detuning BELOW the ZPL
        psb = self._fc(x, self.Sem(P))
        dwf = 0.049 * np.exp(-(np.clip(P, 0, 120) / 120.0) * np.log(0.049 / 0.008))
        return psb + dwf * np.exp(-x ** 2 / (2 * self.zpl_width ** 2))

    def eta_col(self, P, band=None):
        """
        Fraction of NV- emission falling inside the detection passband.

        R = f- * photon_flux * sigma_abs counts ABSORBED photons.  What a confocal setup
        records is the emission that survives a long-pass filter chosen at
        ambient pressure to reject the laser and NV0 -- typically 650-800 nm.
        The emission band blue shifts with the ZPL (peak 719 nm at ambient,
        619 nm at 120 GPa), so a fixed passband collects a rapidly shrinking
        fraction of it.  This factor is monotonic in pressure and is
        WAVELENGTH INDEPENDENT (emission follows Kasha's rule: it leaves the
        relaxed 3E whatever vibronic level was excited), so it cannot move
        lambda_opt -- but it is the missing piece in the pressure dependence of
        the detected count rate.
        """
        lo, hi = self.det_band if band is None else band
        Parr = np.asarray(P, float)
        scalar = Parr.ndim == 0
        if scalar:
            key = (float(Parr), lo, hi)
            if key in self._col_cache:
                return self._col_cache[key]
        # broadcast pressure against the energy grid in one vectorised pass:
        # a Python loop over pressure made the Monte Carlo figures unusably slow.
        E = np.linspace(0.6, 3.2, _COL_GRID)         # eV, spans the whole band
        sem = self.sigma_em(E, Parr[..., None])
        inside = (E >= HBARC / hi) & (E <= HBARC / lo)
        tot = np.trapezoid(sem, E, axis=-1)
        num = np.trapezoid(np.where(inside, sem, 0.0), E, axis=-1)
        out = np.where(tot > 0.0, num / np.where(tot > 0.0, tot, 1.0), 0.0)
        if scalar:
            out = float(out.reshape(()) if out.ndim == 0 else out.ravel()[0])
            self._col_cache[key] = out
        return out

    # ---- steady-state NV- fraction ----
    def photon_flux_factor(self, lam_nm):
        """Photon flux per beam-weight, normalised to one at 532 nm.

        Experimental comparisons in this project hold optical power fixed.  A
        beam of power density I contains I/E_photon photons, hence its optical
        rates carry a factor E_532/E_photon = lambda/532.  The legacy
        ``photon_flux`` basis is retained only to reproduce historical outputs.
        """
        lam_nm = np.asarray(lam_nm, float)
        if self.intensity_basis == 'photon_flux':
            return np.ones_like(lam_nm)
        return lam_nm / 532.0

    def f_minus(self, beams, P):
        """beams: ``(wavelength_nm, relative optical power)`` by default."""
        relu = lambda u: np.clip(u, 0, None)
        s = sum(I * self.photon_flux_factor(lam) * self.sigma_abs(nm2eV(lam), P)
                for lam, I in beams)
        Ggs = sum(self.a_gs * relu(nm2eV(lam) - self.IP_A2(P)) * I
                  * self.photon_flux_factor(lam) for lam, I in beams)
        Gion = Ggs + self.a_es * s
        Grec = self.r0 * s + self.rbg
        return Grec / (Grec + Gion), s

    # ---- contrast: ISC prefactor (C-3) x NV0 dilution ----
    def C0(self, P):
        """
        Wavelength-independent ISC contrast prefactor.  The 3E->1A1 rate falls as
        the 3E-1A1 gap opens with the ZPL blue shift, so C0 = C_amb exp(-dE_ZPL/E_isc);
        the plateau observed above ~50 GPa is produced by the saturation of dE_ZPL
        itself, so no separate floor is needed (C_floor = 0).

        Calibrated to Dai et al. 2022, 532 nm, room temperature: 14% at ambient,
        ~3% at 102 GPa, ~1.2% at 138 GPa -> C_amb = 0.247, E_isc = 0.181 eV.
        The calibration is order-of-magnitude: the two microdiamonds in that work
        differ by a factor ~2 in contrast at the same pressure.

        This prefactor is WAVELENGTH INDEPENDENT and therefore cannot move
        lambda_opt; test_freeze.py asserts that.  It applies to the hydrostatic /
        randomly-oriented case.  For the [111] NV group of a [111]-cut culet the
        3E RED-shifts under the deviatoric stress and contrast is retained to
        megabar pressures (Bhattacharyya thesis Fig. 6.10); no quantitative
        contrast-vs-pressure data exists for that geometry, so it is not modelled.
        """
        if not self.isc:
            return 1.0
        return self.C_floor + (self.C_amb - self.C_floor) * np.exp(-self.dZPL(P) / self.E_isc)

    # ---- lock-in ODMR sensitivity (arb. units; lower = better) ----
    def linewidth(self, P):
        return 1.0 * (1 + 0.15 * np.clip(P, 0, 140) / 100.0)   # mild strain broadening (phenom.)

    def eta(self, beams, P):
        f, s = self.f_minus(beams, P)
        C = self.C0(P) * f / (f + self.w0 * (1 - f))
        R = f * s
        if self.collection:
            R = R * self.eta_col(P)
        R = np.clip(R, 1e-12, None)
        return self.linewidth(P) / (C * np.sqrt(R)), f, s, C, R

    # convenience: single-wavelength sweep
    def eta_lambda(self, lam_nm, P):
        return self.eta([(lam_nm, 1.0)], P)

    # convenience: optimum of the sweep, parabolically refined
    def lambda_opt(self, P, lo=402.0, hi=640.0, step=0.05):
        lam = np.arange(lo, hi + step, step)
        e = np.asarray(self.eta_lambda(lam, P)[0])     # vectorised over lam
        i = int(e.argmin())
        if 0 < i < len(lam) - 1:                    # parabolic refinement
            y0, y1, y2 = e[i - 1], e[i], e[i + 1]
            d = y0 - 2 * y1 + y2
            if d != 0:
                return lam[i] + 0.5 * step * (y0 - y2) / d
        return lam[i]


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


def default_randomiser(rng, T=300.0):
    """
    Standard randomisation for the MC bands.

    C-2: the ZPL is randomised through the two MEASURED quantities, while the
    single-effective-mode approximation is exposed through a +/-15% range on hw.
      dE120  = 0.400 +/- 0.020 eV   (Ho et al., ">400 meV at 120 GPa")
      slope0 = 5.75 +/- 0.25 meV/GPa (Doherty et al.; 5.5 measured / 5.75 theory)
    rather than through (Emax, P0) independently.  dE_ZPL(120 GPa) is the only
    ZPL quantity lambda_opt(120 GPa) responds to, so this is the randomisation
    that actually propagates measurement uncertainty.
    """
    return NVModel(
        dE120=0.400 + 0.020 * rng.normal(),
        slope0=5.75e-3 + 0.25e-3 * rng.normal(),
        T=T,
        hw=0.065 * rng.uniform(0.85, 1.15),
        S_slope=(4.61 - 3.08) * rng.uniform(0.85, 1.15),
        zpl_width=0.015 * rng.uniform(0.7, 1.4),
        a_gs=6.0 * rng.uniform(0.50, 1.60),
        r0=1.2   * rng.uniform(0.70, 1.30),
        rbg=0.15 * rng.uniform(0.60, 1.50),
    )


def legacy_randomiser(rng):
    """The pre-C-2 randomisation, kept so the old band can be reproduced."""
    return NVModel(
        Emax=0.758 * rng.uniform(0.90, 1.10),
        P0=160.0  * rng.uniform(0.85, 1.15),
        T=0.0,
        intensity_basis='photon_flux',
        S_slope=(4.61 - 3.08) * rng.uniform(0.85, 1.15),
        a_gs=6.0 * rng.uniform(0.50, 1.60),
        r0=1.2   * rng.uniform(0.70, 1.30),
        rbg=0.15 * rng.uniform(0.60, 1.50),
        isc=False, collection=False,
    )
