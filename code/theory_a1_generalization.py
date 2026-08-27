"""Numerical execution of Addendum A1 (coincidence/divergence theory).

`theory_freeze_v3_ho_integrated.md` Addendum A1 and `docs/theory_optima_coincidence.md`
state propositions P1-P7 analytically and quote numbers (kappa_A, the 2.5e3 nm^2
conversion factor, the 32/60/83 nm doublet separations, the mechanism table)
without an executable derivation.  This module executes them against the frozen
Ho kernel and reports where the analytic shortcuts hold and where they do not.

Nothing here changes a v3 frozen value.  The optical-limit table is reproduced
by `ho_odmr_sensitivity.py` and is used here only as an input.

Two structural properties of the reconstructed kernel drive most of the results
below, and neither is acknowledged in the analytic notes:

* A(lambda) at 120 GPa is **not unimodal**.  Besides the band maximum at
  440.64 nm it carries local maxima at 475.55, 500.19 and -- decisively -- the
  zero-phonon line at 514.46 nm, which reaches 69 % of the band maximum.
* In the blue window the extracted samples are 3-10 nm apart, so A is a
  piecewise-linear interpolant on a mesh coarser than the displacements P2
  predicts.

Structure
---------
S1  kernel geometry: A(lambda), lambda_abs, the multimodal structure, and the
    scale dependence of the curvature kappa.
S2  P2: the split formula lambda_eta - lambda_PL ~= 2 l_G / kappa_R, tested
    against exact argmax both locally and globally.
S3  P3: the (M)-mediation theorem.  The degenerate doublet generalises to a
    degenerate multiplet on a multimodal kernel.
S4  P4: which mechanism can split the optima, by direct evaluation.
S5  P5: the ordering theorem and the operational power window.
S6  T4: predicting the wavelength scan from a single-wavelength power sweep.
S7  E1: the sub-ZPL hot-band arithmetic.
S8  kernel provenance: the v1 Franck-Condon envelope against the Ho kernel.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from ho_spectrum_model import HBARC, HoPublishedSpectrumModel

PRESSURE = 120.0
KB_EV = 8.617333262e-5  # eV/K

# The main absorption band, excluding the ZPL resonance and the structure that
# section 7 of the analytic note places outside the theory's domain.
MAIN_BAND = (402.0, 500.0)
# The full window over which extracted data exist (E >= 2.396 eV).
DATA_WINDOW = (402.0, 517.0)


# ---------------------------------------------------------------- S1 geometry

class Kernel:
    """A(lambda) = lambda * sigma_abs(hc/lambda, P), normalised to its maximum."""

    def __init__(self, pressure=PRESSURE, window=DATA_WINDOW):
        self.pressure = float(pressure)
        self.optical = HoPublishedSpectrumModel()
        self.lam_min, self.lam_max = (float(x) for x in window)
        grid = np.arange(self.lam_min, self.lam_max + 1e-9, 0.01)
        values = self._raw(grid)
        self.lam_abs = float(grid[int(np.argmax(values))])
        self.a_max = float(np.max(values))

    def _raw(self, lam):
        lam = np.asarray(lam, float)
        return lam * self.optical.sigma_abs(HBARC / lam, self.pressure)

    def a(self, lam):
        """Normalised absorbed-photon proxy, a(lambda_abs) = 1."""
        return self._raw(lam) / self.a_max

    def local_maxima(self, threshold=0.01, step=0.01):
        """All interior local maxima of A above `threshold`, refined."""
        grid = np.arange(self.lam_min, self.lam_max + 1e-9, step)
        values = self.a(grid)
        rising = np.diff(values)
        idx = np.where((rising[:-1] > 0.0) & (rising[1:] <= 0.0))[0] + 1
        return [(float(grid[i]), float(values[i]))
                for i in idx if values[i] > threshold]

    def level_set(self, target, window=None, step=0.01):
        """Every wavelength in `window` at which a(lambda) equals `target`.

        A multimodal kernel can cross a level more than twice, so this returns
        the full crossing set rather than a bracketed pair.
        """
        lo, hi = window or (self.lam_min, self.lam_max)
        grid = np.arange(lo, hi + 1e-9, step)
        values = self.a(grid) - target
        crossings = np.where(np.sign(values[:-1]) * np.sign(values[1:]) < 0.0)[0]
        roots = []
        for i in crossings:
            roots.append(float(brentq(lambda x: self.a(x) - target,
                                      grid[i], grid[i + 1], xtol=1e-10)))
        return roots

    def kappa_chord(self, lam):
        """Chord curvature: ln(a) = -kappa (lam - lam_abs)^2 / 2 through one point.

        This is the definition that reproduces the table in
        `docs/theory_optima_coincidence.md` section 2.  It is a band-averaged
        quantity, not a local second derivative.
        """
        lam = np.asarray(lam, float)
        delta = lam - self.lam_abs
        with np.errstate(divide='ignore', invalid='ignore'):
            return -2.0 * np.log(self.a(lam)) / delta ** 2

    def kappa_fit(self, half_width):
        """Local curvature from a parabolic least-squares fit of ln a.

        Reported because the pointwise second derivative of the reconstructed
        kernel is dominated by interpolation artefacts.
        """
        lam = np.arange(self.lam_abs - half_width,
                        self.lam_abs + half_width + 1e-9, 0.05)
        lam = lam[(lam >= self.lam_min) & (lam <= self.lam_max)]
        coeffs = np.polyfit(lam - self.lam_abs, np.log(self.a(lam)), 2)
        return -2.0 * coeffs[0]

    def kappa_pointwise(self, lam, step=0.05):
        """Second derivative of ln a by central differences (diagnostic only)."""
        lam = np.asarray(lam, float)
        f = lambda x: np.log(self.a(x))
        return -(f(lam + step) - 2.0 * f(lam) + f(lam - step)) / step ** 2

    def tolerance_band(self, penalty=1.05):
        """Contiguous interval around lambda_abs whose penalty is within limit."""
        target = penalty ** -2.0
        roots = self.level_set(target, MAIN_BAND)
        below = [r for r in roots if r < self.lam_abs]
        above = [r for r in roots if r > self.lam_abs]
        return (max(below) if below else self.lam_min,
                min(above) if above else self.lam_max)


# ------------------------------------------------------- S3/S4 response model

@dataclass(frozen=True)
class MediatedResponse:
    """(M)-compliant response: every factor depends on lambda only via Gamma_p.

    `gamma_sat`, `gamma_contrast`, `gamma_width` are the half-scales of rate
    saturation, contrast collapse and linewidth power broadening, in units of
    Gamma_p.  `np.inf` switches a mechanism off.
    """

    gamma_sat: float = np.inf
    gamma_contrast: float = np.inf
    gamma_width: float = np.inf

    def rate(self, gamma_p):
        return gamma_p / (1.0 + gamma_p / self.gamma_sat)

    def contrast(self, gamma_p):
        return 1.0 / (1.0 + gamma_p / self.gamma_contrast)

    def linewidth(self, gamma_p):
        return np.sqrt(1.0 + gamma_p / self.gamma_width)

    def phi(self, gamma_p):
        """Phi = (C/dnu)^2 R.  eta = 1/sqrt(Phi), so argmax Phi == argmin eta."""
        gamma_p = np.asarray(gamma_p, float)
        return ((self.contrast(gamma_p) / self.linewidth(gamma_p)) ** 2
                * self.rate(gamma_p))

    def dlnphi_dlngamma(self, gamma_p):
        """Closed form of d ln Phi / d ln Gamma_p (P4's left-hand side)."""
        gamma_p = np.asarray(gamma_p, float)
        rho_s = gamma_p / self.gamma_sat
        rho_c = gamma_p / self.gamma_contrast
        rho_w = gamma_p / self.gamma_width
        return (1.0
                - rho_s / (1.0 + rho_s)
                - 2.0 * rho_c / (1.0 + rho_c)
                - rho_w / (1.0 + rho_w))

    def dlnr_dlngamma(self, gamma_p):
        gamma_p = np.asarray(gamma_p, float)
        return 1.0 / (1.0 + gamma_p / self.gamma_sat)

    def gamma_star(self, lo=1e-8, hi=1e8):
        """Interior maximiser of Phi, or nan when Phi is monotone increasing."""
        f = self.dlnphi_dlngamma
        if f(hi) > 0.0 or f(lo) < 0.0:
            return float('nan')
        root = brentq(lambda x: f(np.exp(x)), np.log(lo), np.log(hi), xtol=1e-12)
        return float(np.exp(root))


# ------------------------------------------------------- optima on the kernel

def eta_at(kernel, response, lam, gamma_max):
    """eta(lambda) at a power set by Gamma_max = gamma I A_max, up to a constant."""
    return 1.0 / np.sqrt(response.phi(gamma_max * kernel.a(lam)))


def sensitivity_optima(kernel, response, gamma_max, window=None):
    """Minimisers of eta: one point, or an exactly degenerate multiplet.

    Under (M) with an interior maximum of Phi at Gamma_p*, every wavelength
    solving A(lambda) = Gamma_p*/(gamma I) is an exactly equivalent optimum.
    On a multimodal kernel that level set can hold more than two members, and
    on a kernel with a finite window a member can fall off the edge entirely.

    Returns a dict so that window truncation is never silently reported as
    coincidence: `split` is decided by Gamma_max vs Gamma_p*, not by how many
    roots survive inside the window.
    """
    lo, hi = window or (kernel.lam_min, kernel.lam_max)
    star = response.gamma_star()
    if not np.isfinite(star) or gamma_max <= star:
        return {'optima': (kernel.lam_abs,), 'split': False, 'a_star': float('nan'),
                'truncated_blue': False, 'truncated_red': False}
    a_star = star / gamma_max
    return {
        'optima': tuple(kernel.level_set(a_star, (lo, hi))),
        'split': True,
        'a_star': float(a_star),
        # A member is lost off an edge when the kernel never falls to a_star
        # before the window ends.
        'truncated_blue': bool(kernel.a(lo) > a_star),
        'truncated_red': bool(kernel.a(hi) > a_star),
    }


# ------------------------------------------------------------------ sections

def section1_geometry():
    k = Kernel()
    lines = (405.0, 420.0, 430.0, 445.0, 457.0, 473.0, 488.0)
    lam = np.arange(420.0, 500.0, 0.05)
    pointwise = k.kappa_pointwise(lam)
    return {
        'lambda_abs_nm': k.lam_abs,
        'local_maxima': k.local_maxima(),
        'chord_kappa': {lam_i: float(k.kappa_chord(lam_i)) for lam_i in lines},
        'fit_kappa': {hw: float(k.kappa_fit(hw)) for hw in (10.0, 15.0, 20.0, 30.0, 40.0)},
        'pointwise': {
            'min': float(np.min(pointwise)), 'max': float(np.max(pointwise)),
            'median_abs': float(np.median(np.abs(pointwise))),
            'sign_changes': int(np.sum(np.diff(np.sign(pointwise)) != 0)),
        },
        'band_5pct': k.tolerance_band(1.05),
        'band_10pct': k.tolerance_band(1.10),
    }


def section2_split_formula():
    """Test lambda_eta - lambda_PL ~= 2 l_G / kappa_R against exact argmax.

    A constant logarithmic slope l_G is imposed on G = C/dnu and the exact
    minimiser of eta = 1/(G sqrt(R)) is located on the real kernel, once inside
    the main band and once over the whole data window.  In the low-power limit
    R is proportional to A, so ln Phi = 2 l_G (lam - lam_abs) + ln a.
    """
    k = Kernel()
    kappa = k.kappa_fit(15.0)
    fine = np.arange(k.lam_min, k.lam_max + 1e-9, 0.005)
    main = fine[(fine >= MAIN_BAND[0]) & (fine <= MAIN_BAND[1])]
    log_a_fine, log_a_main = np.log(k.a(fine)), np.log(k.a(main))
    rows = []
    for l_g in (0.0005, 0.001, 0.002, 0.004, 0.00625, 0.01, 0.02,
                -0.001, -0.004, -0.00625, -0.02):
        local = float(main[int(np.argmax(2.0 * l_g * (main - k.lam_abs) + log_a_main))])
        glob = float(fine[int(np.argmax(2.0 * l_g * (fine - k.lam_abs) + log_a_fine))])
        rows.append({
            'l_G_per_nm': l_g,
            'predicted_shift_nm': 2.0 * l_g / kappa,
            'local_shift_nm': local - k.lam_abs,
            'global_lambda_eta_nm': glob,
            'jumped_to_zpl': glob > 505.0,
        })
    return {'kappa_used': kappa, 'conversion_factor_nm2': 2.0 / kappa,
            'rows': rows, 'zpl_jump_l_G': _zpl_jump_threshold(k, fine, log_a_fine)}


def _zpl_jump_threshold(kernel, fine, log_a_fine):
    """Smallest positive l_G at which the global optimum jumps to the ZPL."""
    def outside(l_g):
        best = fine[int(np.argmax(2.0 * l_g * (fine - kernel.lam_abs) + log_a_fine))]
        return 1.0 if best > 505.0 else -1.0

    if outside(0.02) < 0.0:
        return float('nan')
    return float(brentq(outside, 1e-5, 0.02, xtol=1e-7))


def section3_multiplet():
    """Solve P3(c) on the real kernel: degeneracy is exact, geometry is not."""
    k = Kernel()
    response = MediatedResponse(gamma_contrast=1.0)  # contrast collapse alone
    star = response.gamma_star()
    kappa = k.kappa_fit(15.0)

    # The blue flank rises monotonically to lambda_abs, so a doublet keeps both
    # members inside the reconstructed window only while a* exceeds the value
    # the kernel takes at the blue edge.
    a_blue_edge = float(k.a(k.lam_min))
    widest = k.level_set(a_blue_edge * 1.0001, DATA_WINDOW)
    max_verifiable = (widest[-1] - widest[0]) if len(widest) > 1 else 0.0

    rows = []
    for ratio in (0.9, 0.7, 0.5):
        gamma_max = star / ratio
        main = sensitivity_optima(k, response, gamma_max, MAIN_BAND)
        full = sensitivity_optima(k, response, gamma_max, DATA_WINDOW)
        members = main['optima']
        etas = [float(eta_at(k, response, x, gamma_max)) for x in full['optima']]
        centre = float(eta_at(k, response, k.lam_abs, gamma_max))
        rows.append({
            'a_star_over_a_max': ratio,
            'main_band_members': tuple(round(x, 2) for x in members),
            'full_window_multiplet': tuple(round(x, 2) for x in full['optima']),
            'multiplicity': len(full['optima']),
            'truncated_blue': main['truncated_blue'],
            'exact_separation_nm': (members[-1] - members[0]) if len(members) > 1 else float('nan'),
            'gaussian_separation_nm': 2.0 * np.sqrt(2.0 * np.log(1.0 / ratio) / kappa),
            'blue_offset_nm': (k.lam_abs - members[0]) if len(members) > 1 else float('nan'),
            'red_offset_nm': (members[-1] - k.lam_abs) if len(members) > 1 else float('nan'),
            'degeneracy_spread': float(np.ptp(etas) / np.mean(etas)),
            'centre_penalty': centre / etas[0],
        })
    return {'gamma_star': star, 'kappa': kappa, 'rows': rows,
            'a_at_blue_edge': a_blue_edge, 'blue_edge_nm': k.lam_min,
            'max_verifiable_separation_nm': max_verifiable}


def section4_mechanisms():
    """P4: evaluate which single mechanism produces an interior maximum of Phi."""
    cases = {
        'rate saturation alone': MediatedResponse(gamma_sat=1.0),
        'power broadening alone': MediatedResponse(gamma_width=1.0),
        'saturation + broadening': MediatedResponse(gamma_sat=1.0, gamma_width=1.0),
        'contrast collapse alone': MediatedResponse(gamma_contrast=1.0),
        'all three': MediatedResponse(gamma_sat=1.0, gamma_contrast=1.0,
                                      gamma_width=1.0),
    }
    k = Kernel()
    out = []
    for name, response in cases.items():
        star = response.gamma_star()
        detail = {'mechanism': name, 'gamma_star': star,
                  'splits': bool(np.isfinite(star)),
                  'slope_at_low_power': float(response.dlnphi_dlngamma(1e-6)),
                  'slope_at_high_power': float(response.dlnphi_dlngamma(1e6))}
        if np.isfinite(star):
            detail['doublet'] = tuple(
                round(x, 2) for x in
                sensitivity_optima(k, response, star / 0.9, MAIN_BAND)['optima'])
        out.append(detail)
    return out


def section4b_mechanism_degeneracy():
    """Can eta alone say *which* mechanism split the optima?  No.

    With equal half-scales, saturation + power broadening and contrast collapse
    give d ln Phi/d ln Gamma_p = 1 - 2 rho/(1 + rho) identically, so they share
    the same Phi, the same Gamma_p*, and the same doublet at every power.  The
    pair (saturation alone, broadening alone) is degenerate the same way.
    """
    gamma = np.logspace(-3, 3, 400)
    pairs = {
        'saturation+broadening vs contrast collapse': (
            MediatedResponse(gamma_sat=1.0, gamma_width=1.0),
            MediatedResponse(gamma_contrast=1.0)),
        'saturation alone vs broadening alone': (
            MediatedResponse(gamma_sat=1.0),
            MediatedResponse(gamma_width=1.0)),
    }
    out = {}
    for label, (first, second) in pairs.items():
        phi_a, phi_b = first.phi(gamma), second.phi(gamma)
        out[label] = {
            'max_slope_difference': float(np.max(np.abs(
                first.dlnphi_dlngamma(gamma) - second.dlnphi_dlngamma(gamma)))),
            'max_shape_difference': float(np.max(np.abs(
                phi_a / phi_a.max() - phi_b / phi_b.max()))),
        }
    # What *does* separate them: the contrast and linewidth spectra themselves.
    probe = 1.0
    out['discriminators_at_gamma_p_equals_1'] = {
        'contrast: sat+broadening': float(
            MediatedResponse(gamma_sat=1.0, gamma_width=1.0).contrast(probe)),
        'contrast: contrast collapse': float(
            MediatedResponse(gamma_contrast=1.0).contrast(probe)),
        'linewidth: sat+broadening': float(
            MediatedResponse(gamma_sat=1.0, gamma_width=1.0).linewidth(probe)),
        'linewidth: contrast collapse': float(
            MediatedResponse(gamma_contrast=1.0).linewidth(probe)),
    }
    return out


def section5_ordering():
    """P5: Phi turns over below the power at which the PL spectrum flattens."""
    response = MediatedResponse(gamma_sat=1.0, gamma_contrast=1.0, gamma_width=1.0)
    star = response.gamma_star()
    k = Kernel()
    # The PL spectrum visibly flattens once R is within 10 % of its asymptote.
    gamma_r_flat = 9.0 * response.gamma_sat
    rows = []
    for gamma_max in np.logspace(-2, 3, 11):
        result = sensitivity_optima(k, response, gamma_max, DATA_WINDOW)
        if not result['split']:
            regime = 'coincident with lambda_PL'
        elif result['truncated_blue']:
            regime = 'blue member off the data window'
        elif gamma_max < gamma_r_flat:
            regime = 'split, PL still single-peaked'
        else:
            regime = 'split, PL flat'
        rows.append({
            'gamma_max': float(gamma_max),
            'n_optima': len(result['optima']),
            'a_star': result['a_star'],
            'lambda_eta_nm': tuple(round(x, 2) for x in result['optima']),
            'R_over_R_sat': float(response.rate(gamma_max) / response.gamma_sat),
            'regime': regime,
        })
    # The doublet is fully contained only while a* exceeds the blue-edge value.
    a_blue_edge = float(k.a(k.lam_min))
    return {'gamma_star_phi': star, 'gamma_r_flat': gamma_r_flat,
            'window_decades': float(np.log10(gamma_r_flat / star)),
            'lambda_PL_nm': k.lam_abs, 'rows': rows,
            'gamma_max_contained_max': float(star / a_blue_edge),
            'containment_decades': float(np.log10(1.0 / a_blue_edge))}


def section6_t4_prediction(seed=0, noise=0.02, n_powers=14, probe_nm=457.0,
                           n_trials=200):
    """T4: predict the wavelength scan from one single-wavelength power sweep.

    Ground truth is a mediated response with all three mechanisms.  A synthetic
    power sweep at `probe_nm` is generated with multiplicative noise, the three
    half-scales are recovered by least squares, Gamma_p* follows from P4, and
    the critical power and multiplet geometry are predicted without ever
    evaluating the kernel away from `probe_nm`.
    """
    k = Kernel()
    truth = MediatedResponse(gamma_sat=1.4, gamma_contrast=0.8, gamma_width=2.2)
    a_probe = float(k.a(probe_nm))
    star_true = truth.gamma_star()
    u = np.logspace(-1.5, 1.5, n_powers)
    gamma_p = u * a_probe

    def one_trial(rng):
        obs_r = truth.rate(gamma_p) * (1.0 + noise * rng.standard_normal(n_powers))
        obs_c = truth.contrast(gamma_p) * (1.0 + noise * rng.standard_normal(n_powers))
        obs_w = truth.linewidth(gamma_p) * (1.0 + noise * rng.standard_normal(n_powers))
        fitted = MediatedResponse(
            gamma_sat=_fit_scale(u, obs_r / u, lambda x: 1.0 / (1.0 + x)) * a_probe,
            gamma_contrast=_fit_scale(u, obs_c, lambda x: 1.0 / (1.0 + x)) * a_probe,
            gamma_width=_fit_scale(u, obs_w, lambda x: np.sqrt(1.0 + x)) * a_probe)
        return fitted

    rng = np.random.default_rng(seed)
    reference = one_trial(rng)
    kappa = k.kappa_fit(15.0)
    predictions = []
    for factor in (1.11, 1.43, 2.0, 4.0):
        u_probe = factor * star_true
        pred = sensitivity_optima(k, reference, u_probe, DATA_WINDOW)
        true = sensitivity_optima(k, truth, u_probe, DATA_WINDOW)
        p, t = pred['optima'], true['optima']
        predictions.append({
            'u_over_u_c': factor,
            'predicted_nm': tuple(round(x, 2) for x in p),
            'true_nm': tuple(round(x, 2) for x in t),
            'predicted_separation_nm': float(p[-1] - p[0]) if len(p) > 1 else float('nan'),
            'true_separation_nm': float(t[-1] - t[0]) if len(t) > 1 else float('nan'),
            # The formula T4 actually quotes, evaluated at the same power.
            'gaussian_separation_nm': float(
                2.0 * np.sqrt(2.0 * np.log(factor) / kappa)),
            'truncated_blue': true['truncated_blue'],
        })

    rng = np.random.default_rng(seed + 1)
    stars = np.array([one_trial(rng).gamma_star() for _ in range(n_trials)])
    stars = stars[np.isfinite(stars)]
    return {
        'probe_nm': probe_nm, 'a_probe': a_probe, 'noise': noise,
        'n_powers': n_powers,
        'recovered': {'gamma_sat': reference.gamma_sat,
                      'gamma_contrast': reference.gamma_contrast,
                      'gamma_width': reference.gamma_width},
        'truth': {'gamma_sat': truth.gamma_sat,
                  'gamma_contrast': truth.gamma_contrast,
                  'gamma_width': truth.gamma_width},
        'gamma_star_true': star_true,
        'gamma_star_fitted': reference.gamma_star(),
        'ensemble': {
            'n_valid': int(stars.size), 'n_trials': n_trials,
            'median': float(np.median(stars)),
            'p16': float(np.percentile(stars, 16)),
            'p84': float(np.percentile(stars, 84)),
            'within_factor_2': float(np.mean((stars > star_true / 2.0)
                                             & (stars < star_true * 2.0))),
        },
        'predictions': predictions,
    }


def _fit_scale(u, observed, shape):
    """Least-squares fit of a single half-scale, amplitude profiled out."""
    def cost(log_scale):
        model = shape(u / np.exp(log_scale))
        amplitude = np.sum(model * observed) / np.sum(model * model)
        return float(np.sum((observed - amplitude * model) ** 2))

    result = minimize_scalar(cost, bounds=(np.log(1e-3), np.log(1e3)),
                             method='bounded')
    return float(np.exp(result.x))


def section7_hot_band():
    zpl_ev = 2.410
    e532 = HBARC / 532.0
    gap = zpl_ev - e532
    wide = Kernel(window=(402.0, 560.0))
    return {
        'zpl_eV': zpl_ev, 'zpl_nm': HBARC / zpl_ev, 'E_532_eV': e532,
        'gap_meV': gap * 1e3,
        'boltzmann_300K': float(np.exp(-gap / (KB_EV * 300.0))),
        'boltzmann_90K': float(np.exp(-gap / (KB_EV * 90.0))),
        'interpolated_a_532': float(wide.a(532.0)),
        'penalty_300K': float(np.sqrt(1.0 / np.exp(-gap / (KB_EV * 300.0)))),
        'penalty_90K': float(np.sqrt(1.0 / np.exp(-gap / (KB_EV * 90.0)))),
    }


def section8_kernel_provenance():
    from nv_model import NVModel
    v1 = NVModel()
    k = Kernel()
    out = {'per_pressure': {}}
    for pressure in (100.0, 120.0):
        kernel = Kernel(pressure=pressure)
        v1_opt = float(v1.lambda_opt(pressure))
        out['per_pressure'][pressure] = {
            'v1_franck_condon_nm': v1_opt,
            'ho_kernel_nm': kernel.lam_abs,
            'gap_nm': v1_opt - kernel.lam_abs,
        }
    out['ho_penalty'] = {lam: float(1.0 / np.sqrt(k.a(lam)))
                         for lam in (440.65, 445.0, 457.0, 473.0, 475.55,
                                     488.0, 514.46)}
    out['band_5pct'] = k.tolerance_band(1.05)
    out['second_local_max'] = k.local_maxima()[1] if len(k.local_maxima()) > 1 else None
    return out


# ----------------------------------------------------------------------- main

def main():
    bar = '=' * 78
    print(bar)
    print('Addendum A1 - numerical execution against the frozen Ho kernel (120 GPa)')
    print(bar)

    geo = section1_geometry()
    print('\n[S1] Kernel geometry')
    print(f"  lambda_abs = {geo['lambda_abs_nm']:.2f} nm")
    print('  local maxima of A (the note assumes a single band):')
    for lam, value in geo['local_maxima']:
        print(f'    {lam:7.2f} nm   a = {value:.4f}   penalty = '
              f'x{1.0 / np.sqrt(value):.4f}')
    print('  chord curvature (reproduces the note\'s table):')
    for lam, value in geo['chord_kappa'].items():
        print(f'    {lam:7.2f} nm : {value:.3e} nm^-2')
    print('  parabolic-fit curvature vs fit half-width:')
    for hw, value in geo['fit_kappa'].items():
        print(f'    +/- {hw:4.0f} nm : {value:.3e} nm^-2   '
              f'(2/kappa = {2.0 / value:.0f} nm^2)')
    p = geo['pointwise']
    print('  pointwise d2 ln a/dlam2 (420-500 nm, diagnostic):')
    print(f"    range [{p['min']:.2e}, {p['max']:.2e}] nm^-2, "
          f"median |.| = {p['median_abs']:.2e}, {p['sign_changes']} sign changes")
    lo, hi = geo['band_5pct']
    print(f'  5 % band  = [{lo:.2f}, {hi:.2f}] nm  '
          f'(blue {geo["lambda_abs_nm"] - lo:.2f}, red {hi - geo["lambda_abs_nm"]:.2f})')
    lo10, hi10 = geo['band_10pct']
    print(f'  10 % band = [{lo10:.2f}, {hi10:.2f}] nm')

    split = section2_split_formula()
    print('\n[S2] P2 split formula vs exact argmax')
    print(f"  kappa (+/-15 nm fit) = {split['kappa_used']:.3e} nm^-2, "
          f"2/kappa = {split['conversion_factor_nm2']:.0f} nm^2")
    print(f"  global optimum jumps to the ZPL once l_G > "
          f"{split['zpl_jump_l_G'] * 100:.3f} %/nm")
    print('    l_G [%/nm]  predicted [nm]  exact local [nm]  global lam_eta [nm]')
    for row in split['rows']:
        flag = '  <-- ZPL' if row['jumped_to_zpl'] else ''
        print(f"    {row['l_G_per_nm'] * 100:9.3f}  {row['predicted_shift_nm']:13.2f}  "
              f"{row['local_shift_nm']:15.2f}  {row['global_lambda_eta_nm']:18.2f}{flag}")

    mult = section3_multiplet()
    print('\n[S3] P3(c) degenerate optima on the real kernel (contrast collapse alone)')
    print(f"  Gamma_p* = {mult['gamma_star']:.4g} Gamma_C")
    print(f"  a at the blue window edge ({mult['blue_edge_nm']:.0f} nm) = "
          f"{mult['a_at_blue_edge']:.4f}; a doublet keeps both members inside the "
          f"reconstructed kernel only while a* exceeds it")
    print(f"  widest verifiable separation = "
          f"{mult['max_verifiable_separation_nm']:.1f} nm")
    for row in mult['rows']:
        if row['truncated_blue']:
            print(f"    a*/a_max = {row['a_star_over_a_max']:.1f}: blue member lies "
                  f"BELOW the data window -- gaussian formula would give "
                  f"{row['gaussian_separation_nm']:.2f} nm, unverifiable")
        else:
            print(f"    a*/a_max = {row['a_star_over_a_max']:.1f}: "
                  f"{row['main_band_members']}, "
                  f"sep {row['exact_separation_nm']:.2f} nm "
                  f"(gaussian {row['gaussian_separation_nm']:.2f} nm), "
                  f"blue/red {row['blue_offset_nm']:.1f}/{row['red_offset_nm']:.1f}")
        print(f"      full window multiplicity {row['multiplicity']}: "
              f"{row['full_window_multiplet']}, "
              f"degeneracy spread {row['degeneracy_spread']:.1e}, "
              f"centre penalty x{row['centre_penalty']:.4f}")

    print('\n[S4] P4 mechanism classification')
    print('    mechanism                   splits   Gamma_p*      slope(low)  slope(high)')
    for row in section4_mechanisms():
        star = row['gamma_star']
        text = '     -   ' if not np.isfinite(star) else f'{star:9.4g}'
        print(f"    {row['mechanism']:26s} {str(row['splits']):6s} {text}   "
              f"{row['slope_at_low_power']:+9.4f}  {row['slope_at_high_power']:+9.4f}")

    degeneracy = section4b_mechanism_degeneracy()
    print('\n[S4b] Can eta identify which mechanism split the optima?')
    for label, row in degeneracy.items():
        if label.startswith('discriminators'):
            continue
        print(f"    {label}:")
        print(f"      max |slope difference| = {row['max_slope_difference']:.2e}, "
              f"max |shape difference| = {row['max_shape_difference']:.2e}")
    disc = degeneracy['discriminators_at_gamma_p_equals_1']
    print('    -> degenerate in eta; separated only by C and dnu measured apart:')
    print(f"       C  : {disc['contrast: sat+broadening']:.4f} (sat+broad) vs "
          f"{disc['contrast: contrast collapse']:.4f} (contrast collapse)")
    print(f"       dnu: {disc['linewidth: sat+broadening']:.4f} (sat+broad) vs "
          f"{disc['linewidth: contrast collapse']:.4f} (contrast collapse)")

    order = section5_ordering()
    print('\n[S5] P5 ordering and the operational window')
    print(f"  Phi peaks at Gamma_p*       = {order['gamma_star_phi']:.4g}")
    print(f"  PL within 10 % of saturation = {order['gamma_r_flat']:.4g}")
    print(f"  intermediate window          = {order['window_decades']:.2f} decades")
    print(f"  lambda_PL stays at           = {order['lambda_PL_nm']:.2f} nm throughout")
    print(f"  doublet fully inside the data window only up to Gamma_max = "
          f"{order['gamma_max_contained_max']:.3g} "
          f"({order['containment_decades']:.2f} decades above Gamma_p*)")
    print('    Gamma_max    n   lambda_eta                R/R_sat  regime')
    for row in order['rows']:
        print(f"    {row['gamma_max']:9.3g}    {row['n_optima']}   "
              f"{str(row['lambda_eta_nm']):25s} {row['R_over_R_sat']:6.3f}  "
              f"{row['regime']}")

    t4 = section6_t4_prediction()
    print(f"\n[S6] T4 prediction from a {t4['probe_nm']:.0f} nm power sweep "
          f"({t4['n_powers']} powers, {t4['noise'] * 100:.0f} % noise)")
    print(f"  recovered Gamma scales: sat={t4['recovered']['gamma_sat']:.4g}, "
          f"contrast={t4['recovered']['gamma_contrast']:.4g}, "
          f"width={t4['recovered']['gamma_width']:.4g}")
    print(f"  true      Gamma scales: sat={t4['truth']['gamma_sat']:.4g}, "
          f"contrast={t4['truth']['gamma_contrast']:.4g}, "
          f"width={t4['truth']['gamma_width']:.4g}")
    print(f"  Gamma_p*: fitted {t4['gamma_star_fitted']:.4g} vs true "
          f"{t4['gamma_star_true']:.4g}")
    ens = t4['ensemble']
    print(f"  ensemble of {ens['n_trials']} noisy sweeps: Gamma_p* median "
          f"{ens['median']:.4g} [16-84 %: {ens['p16']:.4g}, {ens['p84']:.4g}]")
    print(f"  fraction within a factor of two of the truth: "
          f"{ens['within_factor_2'] * 100:.1f} %")
    print('    I/I_c   predicted [nm]                true [nm]                '
          'sep pred/true/gaussian')
    for row in t4['predictions']:
        flag = '  (blue member off window)' if row['truncated_blue'] else ''
        print(f"    {row['u_over_u_c']:5.2f}   {str(row['predicted_nm']):28s}  "
              f"{str(row['true_nm']):22s}  "
              f"{row['predicted_separation_nm']:.1f}/"
              f"{row['true_separation_nm']:.1f}/"
              f"{row['gaussian_separation_nm']:.1f}{flag}")

    hot = section7_hot_band()
    print('\n[S7] E1 sub-ZPL hot band at 532 nm')
    print(f"  ZPL {hot['zpl_eV']:.3f} eV = {hot['zpl_nm']:.1f} nm; "
          f"gap to 532 nm = {hot['gap_meV']:.1f} meV")
    print(f"  Boltzmann: 300 K = {hot['boltzmann_300K']:.3e} "
          f"(penalty x{hot['penalty_300K']:.1f}), "
          f"90 K = {hot['boltzmann_90K']:.3e} (penalty x{hot['penalty_90K']:.0f})")
    print(f"  interpolated a(532) = {hot['interpolated_a_532']:.3e} "
          '-- between the two, derived from neither')

    prov = section8_kernel_provenance()
    print('\n[S8] Kernel provenance: v1 Franck-Condon envelope vs Ho kernel')
    for pressure, row in prov['per_pressure'].items():
        print(f"  {pressure:5.0f} GPa: v1 = {row['v1_franck_condon_nm']:.2f} nm, "
              f"Ho = {row['ho_kernel_nm']:.2f} nm, gap = {row['gap_nm']:+.2f} nm")
    print(f"  Ho kernel second local maximum: {prov['second_local_max'][0]:.2f} nm "
          f"(a = {prov['second_local_max'][1]:.4f})")
    print('  Ho-kernel optical-limit penalty of the candidate lines:')
    for lam, penalty in prov['ho_penalty'].items():
        print(f'    {lam:7.2f} nm : x{penalty:.4f}')
    lo, hi = prov['band_5pct']
    print(f'  Ho-kernel 5 % band = [{lo:.2f}, {hi:.2f}] nm')


if __name__ == '__main__':
    main()
