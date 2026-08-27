"""Anchor A2's splitting exponent in published CW-ODMR data.

A2's Theorem M needs two antecedents: assumption (M), and a unique interior
maximum of Phi = C^2 R / dnu^2 in the pump rate.  The second one is not a
conjecture -- it is what Dreau et al. measured and fitted in 2011.

    A. Dreau, M. Lesik, L. Rondin, P. Spinicelli, O. Arcizet, J.-F. Roch and
    V. Jacques, "Avoiding power broadening in optically detected magnetic
    resonance of single NV defects for enhanced dc magnetic field
    sensitivity", Phys. Rev. B 84, 195204 (2011); arXiv:1108.0178.

Their fitted CW model, valid for s > 1e-2 where intrinsic relaxation is
negligible (their Eqs. 10-12), is

    x = s / (1 + s),                     s = P_opt / P_sat
    R    = R_inf * x                                        (Eq. 12)
    C    = Theta * Omega^2 / (Omega^2 + Gp_inf * Gc_inf * x^2)   (Eq. 10)
    dnu  = (Gc_inf / 2pi) * sqrt(x^2 + Omega^2/(Gp_inf*Gc_inf))  (Eq. 11)

Writing kappa = Omega^2 / (Gp_inf * Gc_inf) and rho = x^2 / kappa,

    R ∝ x,        C ∝ (1+rho)^-1,        dnu ∝ (1+rho)^(1/2),

so in A2's language c = 1, w = 1/2, and the rate carries no extra
(1+rho) factor, s_exp = 0:

    E = 2c + s_exp + 2w = 3.

The one wrinkle is that Dreau's rho is quadratic in the pump (rho ∝ x^2),
not linear.  A2's criterion generalises to E*n > 1 with rho ∝ Gamma_p^n and

    rho* = 1 / (E n - 1),

giving E n = 6 and rho* = 1/5 for the published model.  Splitting is not
marginal there: it is comfortably satisfied by the fitted parameters.

Theorem M survives the saturating pump map.  Under (M) the response depends on
lambda through the excitation parameter xi = I A(lambda) only, and Gamma_p is a
monotone function of xi, so the optimal set is still {lambda : A(lambda) =
xi*/I} and the rung ratios I_k/I_j = a_j/a_k are unchanged.  Only the
conversion from xi to laser power is nonlinear.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

# Dreau et al. fitted values (their Fig. 3(g) caption and Sec. II C).
GP_INF = 5.0e6    # s^-1, polarisation rate at saturation
GC_INF = 8.0e7    # s^-1, optical cycling rate at saturation
THETA = 0.2       # contrast normalisation
R_INF = 250.0e3   # counts/s at saturation
OMEGA_R = 3.0e6   # rad/s, representative Rabi frequency from their Fig. 3


def kappa(omega_r=OMEGA_R, gp_inf=GP_INF, gc_inf=GC_INF):
    """The dimensionless microwave parameter Omega^2/(Gp_inf Gc_inf)."""
    return omega_r ** 2 / (gp_inf * gc_inf)


def x_of_s(s):
    """Dreau's saturating pump variable x = s/(1+s)."""
    s = np.asarray(s, float)
    return s / (1.0 + s)


def rate(s, r_inf=R_INF):
    return r_inf * x_of_s(s)


def contrast(s, omega_r=OMEGA_R, theta=THETA):
    x = x_of_s(s)
    return theta * omega_r ** 2 / (omega_r ** 2 + GP_INF * GC_INF * x ** 2)


def linewidth(s, omega_r=OMEGA_R):
    x = x_of_s(s)
    return (GC_INF / (2.0 * np.pi)) * np.sqrt(x ** 2 + kappa(omega_r))


def phi(s, omega_r=OMEGA_R):
    """Phi = C^2 R / dnu^2; eta = 1/sqrt(Phi), so argmax Phi = argmin eta."""
    return (contrast(s, omega_r) ** 2 * rate(s)
            / linewidth(s, omega_r) ** 2)


def sensitivity(s, omega_r=OMEGA_R, profile='lorentzian'):
    """Dreau's Eq. (3) up to the constant h/(g mu_B), in arbitrary units."""
    p_f = 4.0 / (3.0 * np.sqrt(3.0)) if profile == 'lorentzian' else np.sqrt(
        np.e / (8.0 * np.log(2.0)))
    return p_f * linewidth(s, omega_r) / (contrast(s, omega_r)
                                          * np.sqrt(rate(s)))


def exponents():
    """Read (c, s_exp, w, n) straight off the published functional forms."""
    return {'c': 1.0, 's_exp': 0.0, 'w': 0.5, 'n': 2.0}


def splitting_exponent():
    e = exponents()
    exponent = 2.0 * e['c'] + e['s_exp'] + 2.0 * e['w']
    return {'E': exponent, 'n': e['n'], 'E_times_n': exponent * e['n'],
            'splits': exponent * e['n'] > 1.0,
            'rho_star': 1.0 / (exponent * e['n'] - 1.0)}


def optimum_pump(omega_r=OMEGA_R):
    """Locate the interior maximum of Phi in Dreau's own variables."""
    result = minimize_scalar(lambda ls: -np.log(phi(np.exp(ls), omega_r)),
                             bounds=(np.log(1e-4), np.log(1e4)),
                             method='bounded')
    s_star = float(np.exp(result.x))
    x_star = float(x_of_s(s_star))
    return {
        's_star': s_star,
        'x_star': x_star,
        'rho_star_measured': float(x_star ** 2 / kappa(omega_r)),
        'x_star_closed_form': float(np.sqrt(kappa(omega_r) / 5.0)),
        'kappa': float(kappa(omega_r)),
    }


def first_rung_power(a_first_rung, omega_r=OMEGA_R):
    """Optical power at the ladder's first rung, in units of P_sat.

    Under (M) the excitation parameter xi ∝ I A(lambda) sets everything, and
    Dreau's s is proportional to xi, so I/I_c = s/s*.  The first rung sits at
    I/I_c = 1/a_first_rung.
    """
    star = optimum_pump(omega_r)
    ratio = 1.0 / a_first_rung
    return {'I_over_Ic': ratio,
            's_at_first_rung': star['s_star'] * ratio,
            's_star': star['s_star'],
            'fraction_of_saturation_power': star['s_star'] * ratio}


def headroom(omega_r=OMEGA_R, s_max=10.0):
    """How many rungs are reachable before the practical power ceiling."""
    star = optimum_pump(omega_r)
    return {'max_I_over_Ic': s_max / star['s_star'], 's_max': s_max}


def robustness_over_rabi(omega_values=None):
    """E n > 1 does not depend on the microwave setting."""
    omega_values = omega_values or (1.1e6, 2.5e6, 3.0e6, 3.6e6, 5.5e6)
    rows = []
    for omega in omega_values:
        star = optimum_pump(omega)
        rows.append({'omega_R': omega, 's_star': star['s_star'],
                     'rho_star': star['rho_star_measured'],
                     'x_star': star['x_star']})
    return rows


def main():
    print('=' * 74)
    print('Splitting exponent from Dreau et al., PRB 84, 195204 (2011)')
    print('=' * 74)

    e = exponents()
    split = splitting_exponent()
    print('\n[D1] Exponents read off the published closed forms')
    print(f"  R    ∝ x                      -> s_exp = {e['s_exp']:.1f}")
    print(f"  C    ∝ (1+rho)^-1             -> c     = {e['c']:.1f}")
    print(f"  dnu  ∝ (1+rho)^(1/2)          -> w     = {e['w']:.1f}")
    print(f"  rho  ∝ x^2                    -> n     = {e['n']:.1f}")
    print(f"  E = 2c + s_exp + 2w = {split['E']:.1f}, "
          f"E*n = {split['E_times_n']:.1f} > 1  -> splits: {split['splits']}")
    print(f"  predicted rho* = 1/(E n - 1) = {split['rho_star']:.4f}")

    star = optimum_pump()
    print('\n[D2] Interior maximum of Phi, located numerically')
    print(f"  kappa = Omega^2/(Gp_inf Gc_inf) = {star['kappa']:.4e}")
    print(f"  s*    = {star['s_star']:.5f}  (fraction of P_sat)")
    print(f"  x*    = {star['x_star']:.5f}  "
          f"(closed form sqrt(kappa/5) = {star['x_star_closed_form']:.5f})")
    print(f"  rho*  = {star['rho_star_measured']:.5f}  "
          f"(predicted {split['rho_star']:.5f})")

    print('\n[D3] Is the ladder reachable?  First rung of the 120 GPa example')
    rung = first_rung_power(0.6938)  # the ZPL level from A2's worked example
    print(f"  first rung at I/I_c = {rung['I_over_Ic']:.3f}")
    print(f"  optimal power  s* = {rung['s_star']:.4f} P_sat")
    print(f"  first rung at  s  = {rung['s_at_first_rung']:.4f} P_sat "
          f"({rung['fraction_of_saturation_power'] * 100:.1f} % of saturation)")
    head = headroom()
    print(f"  with s up to {head['s_max']:.0f} P_sat, reachable "
          f"I/I_c <= {head['max_I_over_Ic']:.0f}")

    print('\n[D4] Robustness across the microwave settings Dreau reports')
    print('    Omega_R [rad/s]     s*        rho*      x*')
    for row in robustness_over_rabi():
        print(f"    {row['omega_R']:.2e}    {row['s_star']:.5f}   "
              f"{row['rho_star']:.5f}   {row['x_star']:.5f}")


if __name__ == '__main__':
    main()
