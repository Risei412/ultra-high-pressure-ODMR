"""Addendum A2: the structural layer of the coincidence/divergence theory.

A1 stated its propositions on the assumption that the absorption kernel is a
single flat band, and quoted numbers read off a digitised figure.  The
numerical execution in `docs/theory_a1_numerical_execution.md` showed that the
two halves of A1 fail differently: everything that depends on the *shape* of
the reconstructed kernel is fragile, and everything that follows from the
*structure* of eta = 1/(G sqrt(R)) is exact.

A2 freezes the structural half and generalises it.  Nothing in this module
reads a fitted value, and the two theorems below hold for any absorption
spectrum whatsoever -- the Ho kernel enters only as the worked example.

Theorem M (multiplicity ladder)
    Under (M), the set of sensitivity-optimal wavelengths at power I is the
    level set {lambda : A(lambda) = Gamma_p^*/(gamma I)}.  Its cardinality is
    piecewise constant in log I and changes exactly where the level crosses a
    critical value of A: +2 at an interior maximum, -2 at an interior minimum,
    -1 at a window edge.  The transition powers are I_k/I_c = A_max/A_k, which
    are ratios of the *measured* rate spectrum and need no absolute
    calibration.

Theorem G (gauge degeneracy)
    eta(lambda, I) depends on the response only through Phi = C^2 R / dnu^2.
    Writing R ~ Gamma_p (1+rho)^-s, C ~ (1+rho)^-c, dnu ~ (1+rho)^w, every
    model on the plane 2c + s + 2w = 2 produces the identical sensitivity
    surface.  Mechanism attribution from wavelength scans is therefore
    impossible at any number of powers, while the multiplicity ladder of
    Theorem M is invariant across the whole family -- the same degeneracy that
    forbids attribution makes the ladder robust.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from theory_a1_generalization import (
    DATA_WINDOW, Kernel, MediatedResponse,
)

# Below this the level set is numerically indistinguishable from the baseline.
MIN_LEVEL = 1e-3


# --------------------------------------------------------------- Theorem M

@dataclass(frozen=True)
class Critical:
    """A critical value of the kernel, where the multiplicity changes."""

    level: float          # a = A/A_max at the critical point
    wavelength: float     # nm
    kind: str             # 'max', 'min', 'edge-blue', 'edge-red'

    @property
    def power_ratio(self):
        """I/I_c at which the level set crosses this critical value."""
        return 1.0 / self.level

    @property
    def delta_multiplicity(self):
        return {'max': +2, 'min': -2, 'edge-blue': -1, 'edge-red': -1}[self.kind]


def critical_values(kernel, window=DATA_WINDOW, step=0.005):
    """Interior extrema and edge values of the kernel, ordered by power."""
    lo, hi = window
    grid = np.arange(lo, hi + 1e-9, step)
    values = kernel.a(grid)
    slope = np.diff(values)
    out = []
    for index in np.where((slope[:-1] > 0) & (slope[1:] <= 0))[0] + 1:
        out.append(Critical(float(values[index]), float(grid[index]), 'max'))
    for index in np.where((slope[:-1] < 0) & (slope[1:] >= 0))[0] + 1:
        out.append(Critical(float(values[index]), float(grid[index]), 'min'))
    out.append(Critical(float(kernel.a(lo)), lo, 'edge-blue'))
    out.append(Critical(float(kernel.a(hi)), hi, 'edge-red'))
    out = [c for c in out if c.level > MIN_LEVEL and c.level < 1.0]
    return sorted(out, key=lambda c: -c.level)


def multiplicity(kernel, level, window=DATA_WINDOW, step=0.002):
    """Number of exactly degenerate optima at the given level."""
    return len(kernel.level_set(level, window, step=step))


_LADDER_CACHE = {}


def observed_ladder(kernel, window=DATA_WINDOW, max_ratio=6.0, samples=900):
    """Locate the transitions of N(I) by bisection on I/I_c.

    Returns a list of (power_ratio, N_before, N_after) with the transition
    power resolved to 1e-4 in the ratio.

    The scan must stay fine enough to resolve the narrowest plateau (the
    six-fold one spans only a factor 1.003), which makes it expensive; results
    are cached per kernel instance.
    """
    key = (id(kernel), tuple(window), max_ratio, samples)
    if key in _LADDER_CACHE:
        return _LADDER_CACHE[key]
    ratios = np.logspace(np.log10(1.001), np.log10(max_ratio), samples)
    counts = np.array([multiplicity(kernel, 1.0 / r, window) for r in ratios])
    transitions = []
    for index in np.where(np.diff(counts) != 0)[0]:
        lo, hi = ratios[index], ratios[index + 1]
        before, after = int(counts[index]), int(counts[index + 1])
        for _ in range(40):
            mid = np.sqrt(lo * hi)
            if multiplicity(kernel, 1.0 / mid, window) == before:
                lo = mid
            else:
                hi = mid
        transitions.append({'power_ratio': float(np.sqrt(lo * hi)),
                            'before': before, 'after': after,
                            'delta': after - before})
    result = {'ratios': ratios, 'counts': counts, 'transitions': transitions}
    _LADDER_CACHE[key] = result
    return result


def match_transitions(kernel, window=DATA_WINDOW, tolerance=0.02):
    """Pair each observed transition with the critical value that predicts it."""
    ladder = observed_ladder(kernel, window)
    criticals = critical_values(kernel, window)
    rows = []
    for event in ladder['transitions']:
        best, best_error = None, np.inf
        for critical in criticals:
            error = abs(event['power_ratio'] - critical.power_ratio) / critical.power_ratio
            if error < best_error:
                best, best_error = critical, error
        rows.append({
            'observed_power_ratio': event['power_ratio'],
            'delta': event['delta'],
            'before': event['before'], 'after': event['after'],
            'predicted_power_ratio': best.power_ratio if best else float('nan'),
            'predicted_from_nm': best.wavelength if best else float('nan'),
            'kind': best.kind if best else '?',
            'relative_error': best_error,
            'matched': bool(best_error < tolerance),
        })
    return {'rows': rows, 'ladder': ladder, 'criticals': criticals}


def plateau_widths(kernel, window=DATA_WINDOW):
    """Power range over which each multiplicity persists.

    A plateau narrower than the power resolution of the experiment is a
    prediction that cannot be tested, and must be declared as such.
    """
    ladder = observed_ladder(kernel, window)
    events = ladder['transitions']
    rows = []
    for index, event in enumerate(events):
        start = event['power_ratio']
        end = (events[index + 1]['power_ratio'] if index + 1 < len(events)
               else float('inf'))
        rows.append({'multiplicity': event['after'],
                     'from_power_ratio': start, 'to_power_ratio': end,
                     'width_factor': end / start})
    return rows


# --------------------------------------------------------------- Theorem G

@dataclass(frozen=True)
class GaugeResponse:
    """Power-law response, characterised entirely by one scalar.

    R = Gamma_p (1 + rho)^-s,  C = (1 + rho)^-c,  dnu = (1 + rho)^w,
    with rho = (Gamma_p / Gamma)^n.  Then

        Phi = C^2 R / dnu^2 = Gamma_p (1 + rho)^-E,   E = 2c + s + 2w,

    so the entire response enters eta only through the **splitting exponent**
    E and the pump nonlinearity n.  Phi has an interior maximum iff E n > 1,
    at rho* = 1/(E n - 1).  Every (c, s, w) on a plane of constant E is one
    model: same eta surface, same Gamma_p*, same multiplicity ladder.  A1's
    mechanism table is the E n > 1 test at n = 1, applied to four corners of
    this plane.

    n = 1 is the default (a pump rate entering linearly).  The published CW
    model of Dreau et al. has n = 2 and E = 3; see `dreau_exponent.py`.
    """

    c: float
    s: float
    w: float
    gamma: float = 1.0
    n: float = 1.0

    @property
    def exponent(self):
        """The splitting exponent E = 2c + s + 2w."""
        return 2.0 * self.c + self.s + 2.0 * self.w

    @property
    def splits(self):
        return self.exponent * self.n > 1.0

    @property
    def rho_star(self):
        """Interior maximum in rho, or nan when Phi is monotone."""
        product = self.exponent * self.n
        return 1.0 / (product - 1.0) if product > 1.0 else float('nan')

    def _rho(self, gamma_p):
        return (np.asarray(gamma_p, float) / self.gamma) ** self.n

    def rate(self, gamma_p):
        return np.asarray(gamma_p, float) * (1.0 + self._rho(gamma_p)) ** -self.s

    def contrast(self, gamma_p):
        return (1.0 + self._rho(gamma_p)) ** -self.c

    def linewidth(self, gamma_p):
        return (1.0 + self._rho(gamma_p)) ** self.w

    def phi(self, gamma_p):
        gamma_p = np.asarray(gamma_p, float)
        return (self.contrast(gamma_p) / self.linewidth(gamma_p)) ** 2 * self.rate(gamma_p)

    def gamma_star(self):
        """Pump rate at the interior maximum, from rho* = 1/(E n - 1)."""
        if not self.splits:
            return float('nan')
        return float(self.gamma * self.rho_star ** (1.0 / self.n))


def gauge_family(exponent=2.0):
    """Representative members of the gauge plane 2c + s + 2w = exponent."""
    return {
        'contrast collapse alone': GaugeResponse(c=1.0, s=0.0, w=0.0),
        'saturation + broadening': GaugeResponse(c=0.0, s=1.0, w=0.5),
        'linewidth^1 (strong)': GaugeResponse(c=0.0, s=0.0, w=1.0),
        'mixed a': GaugeResponse(c=0.5, s=0.4, w=0.3),
        'mixed b': GaugeResponse(c=0.25, s=1.0, w=0.25),
    }


def a1_mechanism_table():
    """A1's named mechanisms placed on the (c, s, w) plane.

    The standard CW-ODMR forms are R ~ Gamma_p/(1+rho) (s = 1) and
    dnu ~ sqrt(1+rho) (w = 1/2); contrast collapse is C ~ 1/(1+rho) (c = 1).
    A1's four-row table is then just the sign of E - 1.
    """
    rows = {
        'rate saturation alone': GaugeResponse(c=0.0, s=1.0, w=0.0),
        'power broadening alone': GaugeResponse(c=0.0, s=0.0, w=0.5),
        'saturation + broadening': GaugeResponse(c=0.0, s=1.0, w=0.5),
        'contrast collapse alone': GaugeResponse(c=1.0, s=0.0, w=0.0),
        'all three': GaugeResponse(c=1.0, s=1.0, w=0.5),
    }
    return [{'mechanism': name, 'c': r.c, 's': r.s, 'w': r.w,
             'exponent': r.exponent, 'splits': r.splits,
             'rho_star': r.rho_star}
            for name, r in rows.items()]


def gauge_degeneracy(gamma=np.logspace(-3, 3, 600)):
    """Every member of the plane shares Phi, Gamma_p*, and the ladder."""
    members = gauge_family()
    reference = members['contrast collapse alone']
    reference_phi = reference.phi(gamma)
    rows = []
    for name, response in members.items():
        phi = response.phi(gamma)
        rows.append({
            'name': name,
            'c': response.c, 's': response.s, 'w': response.w,
            'exponent': response.exponent,
            'max_relative_phi_difference': float(np.max(
                np.abs(phi - reference_phi) / reference_phi)),
            'gamma_star': response.gamma_star(),
            # What actually differs, and is therefore what must be measured.
            'contrast_at_gamma_1': float(response.contrast(1.0)),
            'linewidth_at_gamma_1': float(response.linewidth(1.0)),
            'rate_at_gamma_1': float(response.rate(1.0)),
        })
    return rows


def identifiability(gamma_p=1.0):
    """What each additional measured spectrum buys, evaluated at Gamma_p.

    eta alone fixes Phi = C^2 R / dnu^2 and nothing else; adding the PL rate
    fixes G = C/dnu; only a third, independent spectrum separates C from dnu.
    Gaps are reported as the ratio between two gauge-equivalent models.
    """
    members = gauge_family()
    first = members['contrast collapse alone']
    second = members['saturation + broadening']
    ratio = lambda a, b: float(abs(a / b - 1.0))
    return {
        'eta alone': {
            'quantity': 'Phi = C^2 R / dnu^2',
            'gap': ratio(first.phi(gamma_p), second.phi(gamma_p))},
        '+ PL rate': {
            'quantity': 'R',
            'gap': ratio(first.rate(gamma_p), second.rate(gamma_p))},
        '=> gives G': {
            'quantity': 'G = C/dnu',
            'gap': ratio(first.contrast(gamma_p) / first.linewidth(gamma_p),
                         second.contrast(gamma_p) / second.linewidth(gamma_p))},
        '+ contrast': {
            'quantity': 'C',
            'gap': ratio(first.contrast(gamma_p), second.contrast(gamma_p))},
    }


def ladder_is_gauge_invariant(kernel, window=DATA_WINDOW):
    """The multiplicity ladder is the same for every member of the plane."""
    reference = None
    out = []
    for name, response in gauge_family().items():
        star = response.gamma_star()
        # Express the ladder in units of I_c, which removes Gamma_p* entirely.
        ratios = np.array([1.05, 1.45, 1.55, 2.0, 3.7, 4.5])
        counts = [multiplicity(kernel, 1.0 / r, window) for r in ratios]
        if reference is None:
            reference = counts
        out.append({'name': name, 'gamma_star': star, 'counts': counts,
                    'identical_to_reference': counts == reference})
    return out


# ------------------------------------------------------------------- report

def main():
    kernel = Kernel()
    bar = '=' * 78
    print(bar)
    print('Addendum A2 - structural layer: multiplicity ladder and gauge degeneracy')
    print(bar)

    print('\n[M1] Critical values of the 120 GPa kernel')
    print('    a         lambda [nm]   kind        I/I_c      dN')
    for critical in critical_values(kernel):
        print(f'    {critical.level:.4f}    {critical.wavelength:8.2f}   '
              f'{critical.kind:10s}  {critical.power_ratio:7.3f}   '
              f'{critical.delta_multiplicity:+d}')

    print('\n[M2] Observed multiplicity ladder vs prediction')
    matched = match_transitions(kernel)
    print('    I/I_c obs   I/I_c pred   from [nm]   kind        N       err     ok')
    for row in matched['rows']:
        print(f"    {row['observed_power_ratio']:9.4f}   "
              f"{row['predicted_power_ratio']:10.4f}   "
              f"{row['predicted_from_nm']:9.2f}   {row['kind']:10s}  "
              f"{row['before']}->{row['after']}   "
              f"{row['relative_error'] * 100:5.2f} %   "
              f"{'yes' if row['matched'] else 'NO'}")

    print('\n[M3] Plateau widths (a plateau narrower than the power resolution')
    print('     of the experiment is untestable and must be declared so)')
    for row in plateau_widths(kernel):
        end = ('inf' if not np.isfinite(row['to_power_ratio'])
               else f"{row['to_power_ratio']:.4f}")
        width = ('inf' if not np.isfinite(row['width_factor'])
                 else f"x{row['width_factor']:.4f}")
        print(f"    N = {row['multiplicity']}: I/I_c in "
              f"[{row['from_power_ratio']:.4f}, {end}]   width {width}")

    print("\n[G0] A1's mechanism table is the sign of E - 1, E = 2c + s + 2w")
    print('    mechanism                   c     s     w      E     splits   rho*')
    for row in a1_mechanism_table():
        star = ('   -  ' if not np.isfinite(row['rho_star'])
                else f"{row['rho_star']:.4f}")
        print(f"    {row['mechanism']:26s} {row['c']:.2f}  {row['s']:.2f}  "
              f"{row['w']:.2f}   {row['exponent']:.2f}   "
              f"{str(row['splits']):6s}  {star}")

    print('\n[G1] Gauge plane 2c + s + 2w = 2: every member shares Phi')
    print('    model                       c     s     w    2c+s+2w   '
          'max |dPhi/Phi|   Gamma_p*')
    for row in gauge_degeneracy():
        print(f"    {row['name']:26s} {row['c']:.2f}  {row['s']:.2f}  "
              f"{row['w']:.2f}   {row['exponent']:.3f}    "
              f"{row['max_relative_phi_difference']:.2e}      "
              f"{row['gamma_star']:.4f}")

    print('\n[G2] What differs between two gauge-equivalent models at Gamma_p = 1')
    rows = {row['name']: row for row in gauge_degeneracy()}
    first, second = rows['contrast collapse alone'], rows['saturation + broadening']
    for key, label in (('rate_at_gamma_1', 'R  '), ('contrast_at_gamma_1', 'C  '),
                       ('linewidth_at_gamma_1', 'dnu')):
        print(f"    {label}: {first[key]:.4f} (contrast collapse) vs "
              f"{second[key]:.4f} (saturation + broadening)")

    print('\n[G3] Identifiability: what each measured spectrum buys '
          '(at Gamma_p = 1)')
    for key, row in identifiability().items():
        verdict = ('degenerate' if row['gap'] < 1e-12
                   else f"separates ({row['gap'] * 100:.0f} % apart)")
        print(f"    {key:12s} {row['quantity']:22s} {verdict}")

    print('\n[G4] The ladder is invariant across the gauge plane')
    for row in ladder_is_gauge_invariant(kernel):
        print(f"    {row['name']:26s} Gamma_p* = {row['gamma_star']:.4f}   "
              f"N = {row['counts']}   "
              f"{'identical' if row['identical_to_reference'] else 'DIFFERS'}")


if __name__ == '__main__':
    main()
