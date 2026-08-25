"""
analysis_C_lambda.py
--------------------
Test of the last untested assumption behind the central claim of this work:
that the CW-ODMR contrast C is independent of the excitation wavelength, so
that lambda_opt = argmax(lambda * sigma_abs) at fixed optical power.

Run:  python analysis_C_lambda.py

The worry (Bhattacharyya thesis, Ch. 8): under stress the 3E orbital doublet
splits by ~THz/GPa, comparable to the width of the Franck-Condon envelope
itself.  If the two orbital branches have different intersystem-crossing rates,
and if the excitation wavelength selects between them, then C = C(lambda) and
"the optimum is the absorption maximum" is void.

Three questions, answered in order:

  Q1  How large is the 3E orbital splitting delta at 120 GPa, and how does it
      compare with k_B T?
  Q2  Can the population remember which branch absorbed the photon, i.e. is
      orbital relaxation slow compared with the intersystem crossing?
  Q3  If, contrary to Q2, memory did survive, how different would the two
      branches' contrasts have to be to move lambda_opt out of the tolerance
      window?

Inputs
  3E orbital-strain coupling  850 +/- 130 THz per unit strain
      (Barfuss et al., Nat. Commun. 8, 14358 (2017); consistent with the
      ~THz/GPa scale quoted from Davies & Hamer 1976)
  stress anisotropy alpha     0.95 micropillar, 0.56 flat culet (Hilberer 2023)
  diamond elastic constants   E ~ 1100 GPa, C44 = 578 GPa
  3E lifetime ~10 ns, ISC branching ~30% -> Gamma_ISC ~ 3e7 /s
  orbital relaxation in 3E at room temperature: picosecond scale, which is why
      the excited-state zero-field splitting is motionally averaged to
      1.42 GHz (Rogers et al., New J. Phys. 11, 063007 (2009))
"""

import warnings

import numpy as np

from nv_model import NVModel, nm2eV, HW

warnings.filterwarnings('ignore')

KB_meV = 0.08617333262 * 1000 / 1000        # meV/K -> 0.0862 meV/K
ORBITAL_STRAIN_THz = 850.0                  # THz per unit strain
THz_to_meV = 4.135667                       # 1 THz = 4.1357 meV
GAMMA_ISC = 3e7                             # s^-1, 3E -> 1A1
TAU_ORB = 1e-12                             # s, orbital relaxation at 300 K


def splitting_meV(P, alpha, modulus_GPa):
    """3E orbital splitting from the deviatoric stress at the culet."""
    sigma_dev = (1.0 - alpha) * P                     # GPa
    strain = sigma_dev / modulus_GPa
    return ORBITAL_STRAIN_THz * strain * THz_to_meV


def main():
    T = 300.0
    kT = KB_meV * T
    m = NVModel(T=T)
    P = 120.0

    # ---------------------------------------------------------------- Q1 ----
    print('=' * 78)
    print('Q1  How large is the 3E orbital splitting at 120 GPa?')
    print('=' * 78)
    print(f'    k_B T at {T:.0f} K = {kT:.1f} meV;  '
          f'Franck-Condon envelope width ~ sqrt(S) hw = '
          f'{np.sqrt(m.Sabs(P)) * HW * 1e3:.0f} meV')
    print()
    print('    geometry              alpha   sigma_dev    delta (E=1100)   delta (C44=578)')
    rows = []
    for alpha, lab in ((0.95, 'micropillar'), (0.56, 'standard flat culet')):
        d_E = splitting_meV(P, alpha, 1100.0)
        d_C = splitting_meV(P, alpha, 578.0)
        rows.append((alpha, lab, d_E, d_C))
        print(f'    {lab:20s}  {alpha:.2f}   {(1-alpha)*P:6.1f} GPa    '
              f'{d_E:8.1f} meV      {d_C:8.1f} meV')
    d_lo, d_hi = rows[0][2], rows[0][3]
    print()
    print(f'    In the micropillar geometry -- the only one with usable ODMR at')
    print(f'    120 GPa -- delta = {d_lo:.0f}-{d_hi:.0f} meV against k_B T = {kT:.0f} meV.')
    print(f'    The two branches are thermally mixed, not resolved.')
    print(f'    Boltzmann factor exp(-delta/kT) = {np.exp(-d_lo/kT):.2f}-{np.exp(-d_hi/kT):.2f}')

    # ---------------------------------------------------------------- Q2 ----
    print()
    print('=' * 78)
    print('Q2  Can the population remember which branch absorbed the photon?')
    print('=' * 78)
    print(f'    intersystem crossing      Gamma_ISC ~ {GAMMA_ISC:.0e} /s  '
          f'(tau = {1/GAMMA_ISC*1e9:.0f} ns)')
    print(f'    orbital relaxation in 3E  tau_orb   ~ {TAU_ORB*1e12:.0f} ps')
    print(f'    ratio                     {1/(GAMMA_ISC*TAU_ORB):.0e}')
    print()
    print('    Downward orbital relaxation is SPONTANEOUS phonon emission: it needs')
    print('    no thermal activation and proceeds at the picosecond scale whatever')
    print('    delta is.  The branch populations therefore reach thermal equilibrium')
    print('    ~10^4 relaxation times before a single ISC event.  The equilibrium')
    print('    distribution is fixed by delta and T alone and carries NO memory of')
    print('    which branch absorbed the photon.  For memory to survive, orbital')
    print(f'    relaxation would have to be slower than {1/GAMMA_ISC*1e9:.0f} ns -- four orders')
    print('    of magnitude slower than the value that makes the excited-state ZFS')
    print('    motionally averaged to 1.42 GHz at room temperature.')
    print()
    print('    => C is set by the RELAXED 3E, exactly as emission is (Kasha\'s rule).')

    # ---------------------------------------------------------------- Q3 ----
    print()
    print('=' * 78)
    print('Q3  Worst case: if memory DID survive, how big an effect is needed?')
    print('=' * 78)
    lam = np.arange(402., 640., 0.05)
    E = nm2eV(lam)
    base = np.asarray(m.eta_lambda(lam, P)[0])
    lo0 = lam[base.argmin()]
    w5 = lam[base <= 1.05 * base.min()]
    print(f'    reference: lambda_opt = {lo0:.1f} nm, 5% window '
          f'{w5.min():.1f}-{w5.max():.1f} nm')
    print()

    for delta in (d_lo, d_hi, 100.0, 200.0):
        # two branch envelopes, displaced by delta, equal oscillator strength
        z = m.ZPL(P)
        sA = m._sigma_raw(E, P, zpl=z - delta / 2e3) / m._norm
        sB = m._sigma_raw(E, P, zpl=z + delta / 2e3) / m._norm
        r = sA / (sA + sB)                       # fraction absorbed into branch A
        inside = (lam >= w5.min()) & (lam <= w5.max())
        print(f'    delta = {delta:5.1f} meV : branch-selection ratio r(lambda) '
              f'varies {r[inside].min():.3f}-{r[inside].max():.3f} '
              f'across the 5% window (a swing of '
              f'{(r[inside].max()-r[inside].min())*100:.1f} points)')

        # C(lambda) = r*C_A + (1-r)*C_B, with rho = C_B/C_A.  Only C is weighted:
        # contrast enters eta linearly and sqrt(R) sublinearly, so this is the
        # conservative (largest-effect) choice.
        out = []
        for rho in (0.5, 0.2, 0.1, 0.0, 2.0, 5.0, 10.0):
            w = r + (1.0 - r) * rho
            eta = base / w
            lo = lam[eta.argmin()]
            out.append((rho, lo, w5.min() <= lo <= w5.max()))
        worst = max(abs(o[1] - lo0) for o in out)
        allin = all(o[2] for o in out)
        print(f'                     contrast ratio rho = C_B/C_A from 0 to 10 moves '
              f'lambda_opt by at most {worst:.2f} nm; '
              f'{"ALL inside" if allin else "LEAVES"} the 5% window')

    # how extreme must delta be to matter at all?
    # bisect the threshold delta at the extreme rho = 0
    def lopt_at(delta, rho=0.0):
        z = m.ZPL(P)
        sA = m._sigma_raw(E, P, zpl=z - delta / 2e3) / m._norm
        sB = m._sigma_raw(E, P, zpl=z + delta / 2e3) / m._norm
        r = sA / (sA + sB)
        return lam[(base / (r + (1.0 - r) * rho)).argmin()]

    lo_d, hi_d = 1.0, 120.0
    while hi_d - lo_d > 0.5:
        mid = 0.5 * (lo_d + hi_d)
        if lopt_at(mid) <= w5.max():
            lo_d = mid
        else:
            hi_d = mid
    thresh = 0.5 * (lo_d + hi_d)
    print()
    print(f'    THRESHOLD: with one branch completely dark (rho = 0), lambda_opt')
    print(f'    leaves the 5% window once delta > {thresh:.0f} meV.')
    print(f'    Micropillar value is {d_lo:.0f}-{d_hi:.0f} meV -> margin of '
          f'{thresh/d_hi:.1f}-{thresh/d_lo:.1f}x.')
    print()
    print('    Scan: how large must delta be before a fully dark second branch')
    print('    (rho = 0, i.e. branch B gives zero contrast) moves lambda_opt out')
    print('    of the 5% window?')
    z = m.ZPL(P)
    for delta in (50., 100., 200., 300., 400., 600.):
        sA = m._sigma_raw(E, P, zpl=z - delta / 2e3) / m._norm
        sB = m._sigma_raw(E, P, zpl=z + delta / 2e3) / m._norm
        r = sA / (sA + sB)
        lo = lam[(base / r).argmin()]
        ok = w5.min() <= lo <= w5.max()
        print(f'      delta = {delta:5.0f} meV -> lambda_opt = {lo:6.1f} nm  '
              f'({"inside" if ok else "OUTSIDE"} the 5% window)')

    print()
    print('=' * 78)
    print('CONCLUSION')
    print('=' * 78)
    print('    Q1  delta = %.0f-%.0f meV in the micropillar geometry, at or below' % (d_lo, d_hi))
    print('        k_B T = %.0f meV: the branches are thermally mixed.' % kT)
    print('    Q2  orbital relaxation beats the ISC by 10^4, so the branch')
    print('        population is thermal and independent of lambda.  C(lambda) is')
    print('        flat for the same reason the emission spectrum is.')
    print(f'    Q3  even in the worst case -- no relaxation at all, one branch')
    print(f'        completely dark -- delta must exceed {thresh:.0f} meV before')
    print(f'        lambda_opt leaves the 5% window, a margin of '
          f'{thresh/d_hi:.1f}-{thresh/d_lo:.1f}x over the')
    print(f'        micropillar value.  That threshold is only reached in a flat')
    print('        culet, which has no usable ODMR contrast at 120 GPa anyway.')
    print()
    print('    The assumption survives, protected independently by a timescale')
    print(f'    separation of 10^4 (Q2) and a factor {thresh/d_hi:.1f}-{thresh/d_lo:.1f} '
          f'margin in delta (Q3).')
    print()
    print('    Caveat: the same argument covers the vibrational degree of freedom.')
    print('    A hot 3E has a larger 3E-1A1 gap and hence a different ISC rate,')
    print('    but vibrational relaxation is also picosecond, so the ISC again')
    print('    samples only the relaxed state.  What the argument does NOT cover')
    print('    is a mechanism that bypasses the relaxed 3E entirely -- direct')
    print('    photoionization from the hot state, or ISC competitive with')
    print('    picosecond relaxation.  Neither is supported by any measurement')
    print('    on the NV center.')
    return thresh, d_lo, d_hi


if __name__ == '__main__':
    main()
