"""
anvil_orientation.py
--------------------
Which culet orientation, [111] or [100], gives the better CW-ODMR sensitivity
in a DAC at 120 GPa?

The usual argument for a [111] culet is that one of the four NV families then
lies exactly along the load axis, so that family keeps its C3v symmetry and
its transverse strain splitting E vanishes.  This module checks that argument
quantitatively instead of assuming it.

SPIN-STRESS HAMILTONIAN
    For an NV oriented along [111], with the stress tensor written in the
    CUBIC CRYSTAL frame (X,Y,Z = [100],[010],[001]):

        M_z = a1 (s_XX + s_YY + s_ZZ) + 2 a2 (s_YZ + s_ZX + s_XY)
        M_x = b (2 s_ZZ - s_XX - s_YY) + c (2 s_XY - s_YZ - s_ZX)
        M_y = sqrt(3) [ b (s_XX - s_YY) + c (s_YZ - s_ZX) ]

    D_eff = D0 + M_z,   E_eff = sqrt(M_x^2 + M_y^2),
    zero-field ODMR lines at  D_eff +/- E_eff.

    Coupling constants measured by Barson et al., Nano Lett. 17, 1496 (2017):
        a1 = +4.86, a2 = -3.70, b = -2.30, c = +3.50   MHz/GPa

    The other three NV families are obtained by the cubic symmetry operations
    S = diag(s1,s2,s3) with s1*s2*s3 = +1 that map [111] onto each <111>
    direction; under those, s_ij -> s_i s_j s_ij.

DAC STRESS STATE
    Load along the culet normal n.  Mean (hydrostatic) pressure P plus a
    uniaxial deviatoric stress t = s_radial - s_axial > 0 (compression is
    negative):

        s_axial  = -(P + 2t/3),   s_radial = -(P - t/3)
        s = s_radial * I + (s_axial - s_radial) * n (x) n

    t is what a non-hydrostatic medium adds at 120 GPa (a few GPa), and it is
    not uniform across the culet: its spread dt is what broadens the lines.
"""

import numpy as np

# Barson et al. (2017) spin-stress coupling constants, MHz/GPa
A1, A2, B_, C_ = 4.86, -3.70, -2.30, 3.50
D0 = 2870.0                                  # MHz, ambient zero-field splitting

# the four NV axes, as sign patterns with s1*s2*s3 = +1
NV_SIGNS = ((1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1))

CULETS = {'[111]': np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
          '[100]': np.array([1.0, 0.0, 0.0]),
          '[110]': np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)}


# ------------------------------------------------------------------ stress --
def stress_tensor(P, t, n):
    """DAC stress tensor in the cubic crystal frame (GPa, compression < 0)."""
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    s_ax, s_rad = -(P + 2.0 * t / 3.0), -(P - t / 3.0)
    return s_rad * np.eye(3) + (s_ax - s_rad) * np.outer(n, n)


# ----------------------------------------------------------- spin response --
def _MzMxMy(s, a1=A1, a2=A2, b=B_, c=C_):
    """Barson form, valid for an NV along [111] with `s` in the crystal frame."""
    sXX, sYY, sZZ = s[0, 0], s[1, 1], s[2, 2]
    sYZ, sZX, sXY = s[1, 2], s[2, 0], s[0, 1]
    Mz = a1 * (sXX + sYY + sZZ) + 2.0 * a2 * (sYZ + sZX + sXY)
    Mx = b * (2.0 * sZZ - sXX - sYY) + c * (2.0 * sXY - sYZ - sZX)
    My = np.sqrt(3.0) * (b * (sXX - sYY) + c * (sYZ - sZX))
    return Mz, Mx, My


def family_response(s, signs, **kw):
    """(D, E) in MHz for one NV family, given the crystal-frame stress tensor."""
    S = np.diag(np.asarray(signs, float))
    Mz, Mx, My = _MzMxMy(S @ s @ S, **kw)
    return D0 + Mz, np.hypot(Mx, My)


def families(P, t, culet, **kw):
    """List of dicts (one per NV family) for a given culet orientation.

    Each entry carries the family axis, its angle to the load/optical axis,
    D and E, the optical weight, and dD/dt, dE/dt (numerical derivatives used
    for the deviatoric-stress broadening).
    """
    n = CULETS[culet] if isinstance(culet, str) else np.asarray(culet, float)
    n = n / np.linalg.norm(n)
    s = stress_tensor(P, t, n)
    h = 1e-4
    s_p, s_m = stress_tensor(P, t + h, n), stress_tensor(P, t - h, n)

    out = []
    for sg in NV_SIGNS:
        axis = np.array(sg, float) / np.sqrt(3.0)
        D, E = family_response(s, sg, **kw)
        Dp, Ep = family_response(s_p, sg, **kw)
        Dm, Em = family_response(s_m, sg, **kw)
        cos2 = float(np.dot(axis, n)) ** 2
        out.append(dict(axis=axis, signs=sg,
                        theta_deg=np.degrees(np.arccos(min(1.0, abs(np.dot(axis, n))))),
                        D=D, E=E,
                        dDdt=(Dp - Dm) / (2 * h), dEdt=(Ep - Em) / (2 * h),
                        weight=(1.0 + cos2) / 2.0))     # dipole factor for k || n
    return out


# ------------------------------------------------------------------- lines --
def odmr_lines(P, t, culet, dt=0.0, dnu0=5.0, **kw):
    """Group the four families into distinguishable ODMR lines.

    Families that share (D, E) to within a tolerance are merged, since their
    transitions coincide.  Each line gets:
      nu       : transition frequency (MHz)
      weight   : summed optical weight of the contributing families
      contrast : weight / (total PL weight), i.e. contrast dilution by the
                 NV population that does NOT dip at this frequency
      width    : sqrt(dnu0^2 + (|dnu/dt| dt)^2), inhomogeneous deviatoric
                 broadening added in quadrature to the intrinsic width
    """
    fam = families(P, t, culet, **kw)
    total_w = sum(f['weight'] for f in fam)

    raw = []
    for f in fam:
        for sign in (-1.0, +1.0):
            if f['E'] < 1e-9 and sign < 0:
                continue                       # E = 0: the two lines coincide
            raw.append(dict(nu=f['D'] + sign * f['E'],
                            dnudt=f['dDdt'] + sign * f['dEdt'],
                            w=f['weight'] * (1.0 if f['E'] < 1e-9 else 0.5)))

    lines = []
    for r in sorted(raw, key=lambda r: r['nu']):
        if lines and abs(r['nu'] - lines[-1]['nu']) < 1e-6:
            lines[-1]['weight'] += r['w']
        else:
            lines.append(dict(nu=r['nu'], weight=r['w'], dnudt=r['dnudt']))

    for L in lines:
        L['contrast'] = L['weight'] / total_w
        L['width'] = np.hypot(dnu0, abs(L['dnudt']) * dt)
        L['total_pl'] = total_w
    return lines


def best_line(P, t, culet, dt=0.0, dnu0=5.0, **kw):
    """The line with the best sensitivity eta ~ width / (contrast sqrt(R)).

    R is the total detected photon rate, which is the SAME for [111] and [100]
    (sum of dipole weights = 8/3 in both cases), so it does not drive the
    comparison -- contrast dilution and linewidth do.
    """
    lines = odmr_lines(P, t, culet, dt, dnu0, **kw)
    scored = [(L['width'] / (L['contrast'] * np.sqrt(L['total_pl'])), L) for L in lines]
    eta, L = min(scored, key=lambda x: x[0])
    return eta, L


def eta_ratio(P, t, dt, dnu0=5.0, culet_a='[111]', culet_b='[100]', **kw):
    """eta(culet_a)/eta(culet_b); > 1 means culet_b is the better choice."""
    return best_line(P, t, culet_a, dt, dnu0, **kw)[0] / \
           best_line(P, t, culet_b, dt, dnu0, **kw)[0]


def randomiser(rng):
    """Randomised coupling constants for Monte-Carlo bands (Barson uncertainties)."""
    return dict(a1=A1 * rng.uniform(0.92, 1.08),
                a2=A2 * rng.uniform(0.85, 1.15),
                b=B_ * rng.uniform(0.85, 1.15),
                c=C_ * rng.uniform(0.85, 1.15))


if __name__ == '__main__':
    P, t, dt = 120.0, 5.0, 1.0
    for culet in ('[111]', '[100]'):
        print(f'--- culet {culet},  P={P} GPa, t={t} GPa ---')
        for f in families(P, t, culet):
            print(f"  NV{f['signs']}  theta={f['theta_deg']:5.1f} deg  "
                  f"D={f['D']:9.1f}  E={f['E']:7.2f} MHz  "
                  f"dD/dt={f['dDdt']:+6.2f}  dE/dt={f['dEdt']:+6.2f} MHz/GPa")
        for L in odmr_lines(P, t, culet, dt):
            print(f"   line {L['nu']:9.1f} MHz  contrast={L['contrast']:.3f}  "
                  f"dnu/dt={L['dnudt']:+6.2f}  width={L['width']:6.2f} MHz")
        e, L = best_line(P, t, culet, dt)
        print(f'   best eta = {e:.4f}  at {L["nu"]:.1f} MHz')
    print(f'\neta([111])/eta([100]) = {eta_ratio(P, t, dt):.3f}')
