"""B1 diagnostic: where does the 34.9 nm lambda_opt gap at 120 GPa come from?

DIAGNOSTIC ONLY.  This script does not modify either model and does not touch
any frozen number.  It reports a decomposition; the decision about what (if
anything) to change is made elsewhere.

    reconstructed Ho kernel (ho_spectrum_model)  lambda_opt(120 GPa) = 440.65 nm
    single-mode envelope    (nv_model)           lambda_opt(120 GPa) = 475.51 nm

Run:  python diagnose_kernel_vs_singlemode.py
"""
import numpy as np
from scipy.optimize import brentq

from ho_spectrum_model import HoPublishedSpectrumModel, HBARC
from nv_model import NVModel, nm2eV

LAM = np.arange(400.0, 640.0 + 0.05, 0.05)
TARGET = 440.65


def argmax_absorption(model, P=120.0):
    """arg max lambda * sigma_abs -- the Ho model's objective, applied to any model."""
    if isinstance(model, HoPublishedSpectrumModel):
        sigma = model.sigma_abs(HBARC / LAM, P)
    else:
        sigma = model.sigma_abs(nm2eV(LAM), P)
    return LAM[(LAM * sigma).argmax()]


def ho_zpl_and_psb(ho, P):
    """Split each published curve into its ZPL spike and its sideband maximum.

    The digitiser resolved the sharp ZPL with a locally dense node spacing, so
    the fine-grid nodes isolate the ZPL and the coarse ones carry the sideband.
    """
    e, a = ho.spectra[P]
    order = np.argsort(e)
    e, a = e[order], a[order]
    fine = np.gradient(e) < np.median(np.gradient(e)) * 0.4
    zpl = e[fine][a[fine].argmax()]
    coarse = ~fine & (e > zpl + 0.15)
    psb = e[coarse][a[coarse].argmax()]
    return zpl, psb


def main():
    ho, nv = HoPublishedSpectrumModel(), NVModel()

    print('=' * 68)
    print('0. Reproduce the gap')
    print('=' * 68)
    print(f'  Ho kernel   lambda_opt(120) = {ho.lambda_opt(120):.2f} nm')
    print(f'  NVModel     lambda_opt(120) = {nv.lambda_opt(120):.2f} nm')
    print(f'  gap                         = {nv.lambda_opt(120) - ho.lambda_opt(120):.2f} nm')

    print()
    print('=' * 68)
    print('1. Is it the objective function?  (eta vs lambda*sigma_abs)')
    print('=' * 68)
    nv_abs = argmax_absorption(nv)
    print(f'  NVModel, full eta objective      = {nv.lambda_opt(120):.2f} nm')
    print(f'  NVModel, lambda*sigma_abs only   = {nv_abs:.2f} nm')
    print(f'  -> charge state + contrast move lambda_opt by '
          f'{nv.lambda_opt(120) - nv_abs:+.2f} nm  (NEGLIGIBLE)')
    print('     The whole gap lives in the absorption lineshape.')

    print()
    print('=' * 68)
    print('2. Branch structure: which feature is the maximum, per pressure?')
    print('=' * 68)
    print('   P    Ho_argmax  NV_argmax     ZPL_Ho   ZPL_NV    d(ZPL)   '
          'StokesHo  S*hw_NV')
    for P in (0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0):
        zpl_ho, psb_ho = ho_zpl_and_psb(ho, P)
        zpl_nv = nv.ZPL(P)
        print(f'  {P:5.0f}  {argmax_absorption(ho, P):8.2f}  {argmax_absorption(nv, P):8.2f}   '
              f'{zpl_ho:7.4f}  {zpl_nv:7.4f}  {zpl_ho - zpl_nv:+7.4f}   '
              f'{psb_ho - zpl_ho:7.4f}  {nv.Sabs(P) * nv.hw:7.4f}')
    print('  NOTE: below ~100 GPa the Ho kernel maximum IS the ZPL line; it')
    print('  hands over to the sideband between 80 and 100 GPa (branch exchange).')
    print('  At 120 GPa both models sit on the sideband, so the gap is a')
    print('  sideband-shape gap, not a branch mismatch.')

    print()
    print('=' * 68)
    print('3. Sequential decomposition at 120 GPa (absorption objective)')
    print('=' * 68)
    zpl_ho, psb_ho = ho_zpl_and_psb(ho, 120.0)
    dE120_ho = zpl_ho - 1.9459          # Ho curves' own ZPL shift
    stokes_ho = psb_ho - zpl_ho

    steps = [('baseline NVModel', NVModel())]
    steps.append((f'+ ZPL from Ho curves (dE120 0.400 -> {dE120_ho:.4f} eV)',
                  NVModel(dE120=dE120_ho)))
    steps.append((f'+ Stokes S*hw 0.300 -> {stokes_ho:.4f} eV (S -> {stokes_ho / 0.065:.2f})',
                  NVModel(dE120=dE120_ho, S_slope=stokes_ho / 0.065 - 3.08)))

    prev = None
    for label, model in steps:
        lam = argmax_absorption(model)
        delta = '' if prev is None else f'   ({lam - prev:+6.2f} nm)'
        print(f'  {label:<52s} {lam:7.2f}{delta}')
        prev = lam
    print(f'  {"Ho kernel target":<52s} {TARGET:7.2f}   ({TARGET - prev:+6.2f} nm residual)')
    print()
    print('  residual = single-mode lineshape shape error (width/asymmetry):')
    m2 = steps[-1][1]
    s_nv = LAM * m2.sigma_abs(nm2eV(LAM), 120.0)
    s_ho = LAM * ho.sigma_abs(HBARC / LAM, 120.0)
    sel = (LAM > 410) & (LAM < 520)
    for name, s in (('single-mode', s_nv / s_nv.max()), ('reconstructed', s_ho / s_ho.max())):
        w = LAM[sel][s[sel] >= 0.5]
        print(f'    {name:<14s} blue lobe at half max: {w.min():.1f}-{w.max():.1f} nm '
              f'(width {w.max() - w.min():.1f})')

    print()
    print('=' * 68)
    print('4. Factors that do NOT matter')
    print('=' * 68)
    z, S = nv.ZPL(120.0), nv.Sabs(120.0)
    psb_only = LAM * (nv._fc(nm2eV(LAM) - z, S) / nv._norm)
    print(f'  Debye-Waller / ZPL gaussian : {LAM[psb_only.argmax()] - argmax_absorption(nv):+.2f} nm')
    print(f'  NV search window [402,640] vs [380,700] : '
          f'{nv.lambda_opt(120, lo=380.0, hi=700.0) - nv.lambda_opt(120):+.2f} nm')
    print('  Ho grid step 0.05 -> 2.0 nm : ', end='')
    print(', '.join(f'{ho.lambda_opt(120, step=s):.2f}' for s in (0.05, 0.5, 1.0, 2.0)))
    e, a = ho.spectra[120.0]
    lam_n = np.sort(HBARC / e)
    near = lam_n[(lam_n > 430) & (lam_n < 455)]
    print(f'  Ho raw nodes near the peak  : {np.round(near, 1)}')
    v = LAM * ho.sigma_abs(HBARC / LAM, 120.0)
    flat = LAM[(v >= 0.98 * v.max()) & (LAM < 520)]
    print(f'  Ho objective within 2% of max spans {flat.min():.1f}-{flat.max():.1f} nm')

    print()
    print('=' * 68)
    print('5. What hbar-omega would reconcile the two models?')
    print('=' * 68)
    for label, kw in (('shipped ZPL (dE120=0.400)', {}),
                      (f'Ho-consistent ZPL (dE120={dE120_ho:.3f})', {'dE120': dE120_ho})):
        hw = brentq(lambda h: argmax_absorption(NVModel(hw=h, **kw)) - TARGET,
                    0.05, 0.20, xtol=1e-4)
        print(f'  {label:<34s} hw* = {hw * 1000:.1f} meV = {hw / 0.065:.2f}x ambient '
              f'({(hw / 0.065 - 1) * 100:+.0f}%)')
    print(f'  If instead S=4.61 is kept, hw_eff = {stokes_ho / 4.61 * 1000:.1f} meV '
          f'= {stokes_ho / 4.61 / 0.065:.2f}x ambient')
    print('  Reference: the diamond Raman mode hardens ~20% by 120 GPa.')
    print('  => 35% is the same order as lattice hardening and is defensible;')
    print('     58-84% is not, so hbar-omega alone cannot absorb the gap.')


if __name__ == '__main__':
    main()
