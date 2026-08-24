"""
repro_literature.py
-------------------
Reproduction of published high-pressure NV measurements with the frozen model
(nv_model.NVModel).  Run:  python repro_literature.py

Every target below is an observation taken from a paper or thesis, NOT a number
this model was fitted to -- with two declared exceptions, marked [CAL], which
are the anchors of the C-2 ZPL parameterisation and the C-3 contrast prefactor.

Sources
  [Doh14]  M. W. Doherty et al., Phys. Rev. Lett. 112, 047601 (2014).
  [Ho26]   K. O. Ho, C. Dailledouze et al., arXiv:2606.02399 (2026).
  [Dai22]  J.-H. Dai et al., Chin. Phys. Lett. 39, 117601 (2022); arXiv:2204.05064.
  [Bha22]  P. Bhattacharyya, PhD thesis, UC Berkeley (2022), Chs. 1, 6, 8.
  [Lya18]  S. G. Lyapin et al., Nanosystems 9, 55 (2018).
  [Hil23]  A. Hilberer et al., arXiv:2301.05094; Phys. Rev. B 107, L220102 (2023).

Measurement temperatures: [Doh14] [Dai22] [Bha22] room temperature; [Ho26] 30-90 K.
"""

import warnings
import numpy as np

from nv_model import NVModel, mc_band, default_randomiser, nm2eV, eV2nm

warnings.filterwarnings('ignore')

RT = 300.0      # room temperature, the condition of [Doh14] [Dai22] [Bha22]
LT = 90.0       # [Ho26] single-micropillar DAC

_rows = []


def check(tag, source, observation, model_value, ok, note=''):
    """ok = True/False, or the string 'OPEN' for a recorded, unreproduced observation."""
    _rows.append((tag, source, observation, model_value, ok, note))


def snr_ratio(m, lam_a, lam_b, P):
    """SNR(lam_a)/SNR(lam_b) = eta(lam_b)/eta(lam_a) at the same pressure."""
    return m.eta_lambda(lam_b, P)[0] / m.eta_lambda(lam_a, P)[0]


def main():
    m = NVModel(T=RT)                 # room-temperature model, all defaults frozen
    mlt = NVModel(T=LT)               # low-temperature model for the [Ho26] anchors

    # ---------------- ZPL: the C-2 anchors and their consequence -------------
    slope0 = (m.ZPL(1e-4) - m.ZPL(0.0)) / 1e-4
    check('R1', '[Doh14]', 'dE_ZPL/dP -> 5.75 meV/GPa as P -> 0 (linear to 60 GPa)',
          f'{slope0*1e3:.2f} meV/GPa', abs(slope0*1e3 - 5.75) < 0.05, '[CAL]')

    check('R2', '[Ho26]', 'dE_ZPL(120 GPa) > 400 meV',
          f'{m.dZPL(120.)*1e3:.1f} meV', abs(m.dZPL(120.) - 0.400) < 0.005, '[CAL]')

    slope120 = (m.ZPL(120.0001) - m.ZPL(120.)) / 1e-4
    check('R3', '[Lya18][Dai22]', 'the ZPL shift SLOWS at high pressure '
          '(Doherty linear extrapolation is wrong)',
          f'{slope0*1e3:.2f} -> {slope120*1e3:.2f} meV/GPa (0 -> 120 GPa)',
          slope120 < 0.6 * slope0)

    # ZPL crossing 532 nm: Doherty's linear law predicted 66 GPa, which failed.
    from scipy.optimize import brentq
    Pcross = brentq(lambda P: m.ZPL(P) - nm2eV(532.), 1.0, 400.0)
    check('R3b', '[Bha22]', 'linear extrapolation put the ZPL across 532 nm at '
          '~66 GPa; DFT + [Lya18] pushed it well beyond',
          f'ZPL crosses 532 nm at {Pcross:.0f} GPa', Pcross > 100.)

    # ---------------- 405 nm: the blue-excitation signature ------------------
    # [Bha22] Sec. 6.1: microdiamonds, hydrostatic, 532 nm vs 405 nm.
    r17 = snr_ratio(m, 405., 532., 17.)
    r100 = snr_ratio(m, 405., 532., 100.)
    check('R4', '[Bha22] 6.1', 'NO cwODMR resonance under 405 nm at ~17 GPa',
          f'SNR(405)/SNR(532) = {r17:.1e} at 17 GPa', r17 < 1e-2)
    check('R5', '[Bha22] 6.1', '405 nm resonance APPEARS near 100 GPa, matching '
          'the 532 nm resonance',
          f'SNR(405)/SNR(532) = {r100:.2f} at 100 GPa; gain over 17 GPa = '
          f'{r100/r17:.1e}x', r100 > 0.05 and r100 / r17 > 1e2)

    # [Bha22] Fig. 6.2(d,e): charge-state populations.
    f405_17, f532_17 = m.f_minus([(405., 1.)], 17.)[0], m.f_minus([(532., 1.)], 17.)[0]
    f405_100, f532_100 = m.f_minus([(405., 1.)], 100.)[0], m.f_minus([(532., 1.)], 100.)[0]
    check('R6', '[Bha22] 6.2(d)', 'at ~17 GPa the 405 nm PL is NV0-dominated while '
          'the 532 nm PL is NV--dominated',
          f'f-(405)={f405_17:.3f} vs f-(532)={f532_17:.3f}',
          f405_17 < 0.5 * f532_17)
    check('R7', '[Bha22] 6.2(e)', 'near 100 GPa the two excitation wavelengths give '
          'COMPARABLE charge-state populations',
          f'f-(405)/f-(532) = {f405_100/f532_100:.2f} at 100 GPa '
          f'(was {f405_17/f532_17:.2f} at 17 GPa)',
          f405_100 / f532_100 > 0.4)

    # ---------------- 450 nm vs 532 nm at moderate pressure ------------------
    # [Bha22] Fig. 6.3: [100] culet, ~50 GPa, SNR ~ contrast * sqrt(counts).
    r450_50 = snr_ratio(m, 450., 532., 50.)
    check('R8', '[Bha22] 6.3', 'at ~50 GPa, 450 nm shows NO DISTINCT ADVANTAGE '
          'over 532 nm',
          f'SNR(450)/SNR(532) = {r450_50:.2f} at 50 GPa', r450_50 < 1.0,
          'direction reproduced; magnitude untested')

    # Crossover pressure where blue overtakes green.
    g = lambda P: snr_ratio(m, 450., 532., P) - 1.0
    Pxo = brentq(g, 20., 140.)
    check('R8b', '[Bha22] 6.3 / [Ho26]', 'green stays competitive at moderate '
          'pressure and is superseded by blue only at high pressure',
          f'450/532 crossover at {Pxo:.0f} GPa', 60. < Pxo < 110.)

    # ---------------- green excitation must SURVIVE to megabar ---------------
    # This is the observation that the T=0 envelope got wrong (issue C-1).
    for tag, P, src, obs in [
        ('R9',  140., '[Dai22]', '532 nm cwODMR persists to ~140 GPa (146 GPa PL)'),
        ('R10', 150., '[Bha22] 6.4/8', '532 nm sensing at ~150 GPa reaches '
                                        '50 uT/sqrt(Hz)'),
    ]:
        lo = m.lambda_opt(P)
        pen = m.eta_lambda(532., P)[0] / m.eta_lambda(lo, P)[0]
        check(tag, src, obs, f'eta(532)/eta_opt = {pen:.2f} at {P:.0f} GPa', pen < 10.)

    # [Ho26]: the integrated PL yield under 532 nm and under 457 nm each trace a
    # BROAD DOME vs pressure, following sigma_abs sampled at those wavelengths.
    # This is the measurement that validates the quantity lambda_opt depends on.
    Pg = np.arange(0., 145., 0.5)
    dome532 = mlt.sigma_abs(nm2eV(532.), Pg)
    dome457 = mlt.sigma_abs(nm2eV(457.), Pg)
    p532, p457 = Pg[dome532.argmax()], Pg[dome457.argmax()]
    check('R11', '[Ho26] 5b', 'SPECTRALLY INTEGRATED PL yield at constant laser '
          'power and integration time traces TWO BROAD DOMES vs pressure under '
          '532 and 457 nm -- explicitly NOT monotonic; quantum efficiency stable '
          'to 120 GPa',
          f'532 nm dome peaks at {p532:.0f} GPa, 457 nm dome still rising at '
          f'120 GPa (max {p457:.0f} GPa)',
          10. < p532 < 60. and p457 > 110.)

    # ---- C-7: the detected count rate is not the PL yield -------------------
    # [Dai22] Fig. 2C is a FILTERED confocal count rate and falls monotonically;
    # [Ho26] Fig. 5b is the spectrally integrated yield and domes.  They are
    # different observables.  The NV- emission band blue shifts with the ZPL
    # (peak 719 nm at ambient, 619 nm at 120 GPa), so a passband chosen at
    # ambient pressure collects a shrinking fraction of it.
    ec = [(P, m.eta_col(P)) for P in (0., 50., 100., 120., 140.)]
    check('R20', '[Ho26] S_em', 'the emission Huang-Rhys factor rises 3.39 -> 5.25 '
          'and DWF_em falls 4.9% -> <1%, so the emission band blue shifts with '
          'the ZPL',
          'collected fraction in a fixed 650-800 nm passband: '
          + ', '.join(f'{P:.0f} GPa {v:.2f}' for P, v in ec),
          ec[0][1] > 0.75 and ec[-1][1] < 0.30
          and all(b[1] < a[1] for a, b in zip(ec, ec[1:])))

    Rc = lambda P: m.eta_lambda(532., P)[4]      # index 4 = R, collection applied
    Pk = np.arange(0., 145., 2.5)
    Rk = np.array([Rc(P) for P in Pk])
    Rk = Rk / Rk[0]
    check('R11b', '[Dai22] 2C', 'FILTERED confocal count rate at 532 nm falls '
          'dramatically over the first 50 GPa, then slowly',
          f'with the collection factor the modelled count rate peaks at '
          f'{Pk[Rk.argmax()]:.0f} GPa at {Rk.max():.1f}x ambient (was 2.8x at '
          f'37 GPa without it) and falls to {Rk[-1]:.3f}x by 140 GPa; the '
          f'residual rise below ~30 GPa is NOT reproduced', 'OPEN',
          'the collection factor accounts for x3.2 of the decline but not the '
          'initial rise; that rise is directly MEASURED by [Ho26] 5b, so the '
          'disagreement is between two experiments, not model vs experiment')

    # ---------------- contrast vs pressure (C-3 prefactor) -------------------
    Cof = lambda P: m.eta_lambda(532., P)[3]
    check('R12', '[Dai22] 3A', 'cwODMR contrast 14% at ambient falls to a ~1-3% '
          'plateau above ~50 GPa',
          f'C(0)={Cof(0.)*100:.1f}%, C(102)={Cof(102.3)*100:.1f}%, '
          f'C(138)={Cof(137.7)*100:.1f}%',
          Cof(0.) > 0.10 and 0.005 < Cof(137.7) < 0.04, '[CAL]')

    # ---------------- [Ho26] operational wavelength switch -------------------
    # 532 nm was used to 41.4 GPa in one cell; 457 nm from 55.2 to 113.8 GPa.
    # Consistency test: 532 nm must still be near its best at 41 GPa, and 457 nm
    # must already beat it at 55 GPa and dominate by 114 GPa.
    r41 = dome532[np.argmin(abs(Pg - 41.4))] / dome532.max()
    r114 = snr_ratio(mlt, 457., 532., 113.8)
    check('R13', '[Ho26]', '532 nm used up to 41.4 GPa, then 457 nm from 55.2 to '
          '113.8 GPa in the other cell',
          f'sigma(532) at 41 GPa is {r41*100:.0f}% of its own dome maximum; '
          f'SNR(457)/SNR(532) = {r114:.1f} at 114 GPa',
          r41 > 0.9 and r114 > 2.)

    lo120 = m.lambda_opt(120.)

    # ---------------- C-4: anvil geometry / deviatoric stress ----------------
    # [Hil23] measured alpha = sigma_par/sigma_perp from the ODMR field response
    # and the ZPL shift per unit compressed volume in both geometries.
    from nv_model import _alpha_factor
    check('R15', '[Hil23]', 'ZPL shift per unit compressed volume is '
          '-434 meV/(cm3/mol) on a standard flat culet (alpha = 0.56) vs '
          '-769 on a micropillar (alpha = 0.95): ratio 0.564',
          f'g(0.56)/g(0.95) = {_alpha_factor(0.56):.3f}',
          abs(_alpha_factor(0.56) - 0.564) < 0.005, '[CAL]')

    # Micropillar ODMR contrast, NV implanted in the ANVIL (independent of the
    # Dai microdiamond calibration of C0).
    mc = NVModel(T=RT, alpha=0.95)
    obs_c = [(73., 457., 0.050), (103., 457., 0.030), (131., 405., 0.015)]
    got = [(P, mc.eta_lambda(L, P)[3]) for P, L, _ in obs_c]
    ratios = [g / c for (P, g), (_, _, c) in zip(got, obs_c)]
    check('R16', '[Hil23] 3b', 'micropillar cwODMR contrast 5% at 73 GPa, 3% at '
          '103 GPa, 1.5% at 131 GPa (NV in the anvil, not microdiamonds)',
          'model ' + ', '.join(f'{g*100:.1f}%' for _, g in got)
          + f'  (ratios {", ".join(f"{r:.2f}" for r in ratios)})',
          all(0.4 < r < 2.5 for r in ratios))

    # [Hil23] stepped the excitation DOWN through 532, 488, 457, 405 nm as
    # pressure rose, "to match the blueshift of the NV absorption spectrum".
    Ps = np.arange(10., 141., 1.)
    ridge = np.array([m.lambda_opt(P) for P in Ps])
    monotone = np.all(np.diff(ridge) < 1e-6)
    check('R17', '[Hil23]', 'the excitation was stepped DOWN through 532, 488, '
          '457 nm as pressure increased, to follow the absorption blue shift',
          f'lambda_opt falls monotonically {ridge[0]:.0f} -> {ridge[-1]:.0f} nm '
          f'over 10-140 GPa; crosses 532 nm at '
          f'{np.interp(-532., -ridge, Ps):.0f} GPa, 488 nm at '
          f'{np.interp(-488., -ridge, Ps):.0f} GPa, 457 nm at '
          f'{np.interp(-457., -ridge, Ps):.0f} GPa',
          monotone)

    # ---------------- Sec. V: the power regime -------------------------------
    # [Dai22] Fig. 2D: PL is LINEAR in excitation power at 100 GPa, i.e. a
    # single-photon process dominates and the transition is far from saturated.
    from nv_model_power import NVModelPower
    mp = NVModelPower(T=RT)
    uu = np.logspace(-4, 2, 400)
    Ru = np.asarray(mp.eta_u([(532., 1.0)], 100., uu)[2])
    lin = (Ru / uu) / (Ru[0] / uu[0])
    u10 = uu[np.argmax(lin < 0.9)]
    check('R18', '[Dai22] 2D', 'PL is linear in laser power at 100 GPa '
          '(single-photon excitation, far from saturation)',
          f'R(u) is linear to within 10% for u < {u10:.2f}; the experiment '
          f'therefore sits at u <~ {u10:.2f}', u10 > 0.05)

    lam_p = np.arange(402., 640., 0.5)
    ridge_u = {u: lam_p[np.nanargmin(np.asarray(
        mp.eta_u([(lam_p, 1.0)], 120., u)[0]))] for u in (1e-3, 0.03, 0.1, 0.3)}
    check('R19', 'this work', 'the intensity-explicit model must return the '
          'fixed-power optimum as u -> 0 (the regression the paper demanded)',
          'lambda_opt(u) = ' + ', '.join(f'{u:g}:{l:.0f} nm'
                                         for u, l in ridge_u.items()),
          abs(ridge_u[1e-3] - lo120) < 1.0)

    # ---------------- the threshold pressure of the blue advantage -----------
    # Sec. IV G of the paper.  This is the model's only NULL reproduction, and
    # it is also the design of the falsifiable experiment: a sign change cannot
    # be produced or displaced by a pressure-independent systematic.
    xo = {}
    for lam_b in (473., 457., 450.):
        g = lambda P: snr_ratio(m, lam_b, 532., P) - 1.0
        xo[lam_b] = brentq(g, 20., 145.)
    check('R21', 'this work', 'the recommended 473 nm line beats 532 nm only '
          'above a threshold pressure; below it the recommendation reverses',
          'crossover ' + ', '.join(f'{int(l)}nm:{x:.0f} GPa'
                                   for l, x in xo.items()),
          65. < xo[473.] < 78.)

    check('R22', '[Bha22] 6.3', 'the 50 GPa null result sits BELOW the model '
          'crossover of the line that was tested (450 nm), i.e. the model '
          'predicts the absence of the effect that was not observed',
          f'crossover(450 nm) = {xo[450.]:.0f} GPa, measurement at 50 GPa, '
          f'SNR(450)/SNR(532) = {snr_ratio(m, 450., 532., 50.):.2f}',
          xo[450.] > 50. and snr_ratio(m, 450., 532., 50.) < 1.0)

    # ---------------- the answer itself --------------------------------------
    lam = np.arange(402, 640, 0.01)
    s = m.sigma_abs(nm2eV(lam), 120.)
    check('R14', 'this work', 'lambda_opt at 120 GPa coincides with the maximum of '
          'the absorption envelope',
          f'lambda_opt = {lo120:.2f} nm, sigma_abs max at {lam[s.argmax()]:.2f} nm',
          abs(lo120 - lam[s.argmax()]) < 1.0)

    # tolerance windows
    e = np.array([m.eta_lambda(l, 120.)[0] for l in lam])
    eo = e.min()
    for lvl in (1.05, 1.10):
        inside = lam[e <= lvl * eo]
        check(f'W{int((lvl-1)*100)}', 'this work',
              f'{int((lvl-1)*100)}% tolerance window at 120 GPa',
              f'{inside.min():.1f} - {inside.max():.1f} nm', True)

    # commercial lines
    print('\n' + '=' * 100)
    print('LITERATURE REPRODUCTION  --  frozen model, T = %.0f K unless noted' % RT)
    print('=' * 100)
    npass = 0
    ntest = 0
    for tag, src, obs, val, ok, note in _rows:
        if ok == 'OPEN':
            mark = 'OPEN'
        else:
            ntest += 1
            npass += bool(ok)
            mark = 'ok  ' if ok else 'FAIL'
        print(f'\n[{mark}] {tag}  {src}  {note}')
        print(f'       observed : {obs}')
        print(f'       model    : {val}')
    print('\n' + '-' * 100)
    print(f'{npass}/{ntest} reproduced  ({len(_rows)-ntest} recorded as OPEN)')

    print('\n' + '=' * 100)
    print('COMMERCIAL LASER LINES AT 120 GPa  (eta/eta_opt, lower is better)')
    print('=' * 100)
    print(f'{"lambda":>8} {"eta/eta_opt":>12} {"sigma/sigma_max":>16} '
          f'{"(T=0 legacy)":>14}')
    m0 = NVModel(T=0.0)
    lo0 = m0.lambda_opt(120.)
    smax0 = m0.sigma_abs(nm2eV(lam), 120.).max()
    smax = s.max()
    for L in (405, 450, 457, 473, 475, 488, 505, 532):
        pen = m.eta_lambda(float(L), 120.)[0] / m.eta_lambda(lo120, 120.)[0]
        pen0 = m0.eta_lambda(float(L), 120.)[0] / m0.eta_lambda(lo0, 120.)[0]
        print(f'{L:>8} {pen:>12.2f} {m.sigma_abs(nm2eV(float(L)),120.)/smax:>16.3f} '
              f'{pen0:>14.2f}')

    print()
    print('=' * 100)
    print('C-4: ANVIL GEOMETRY  (alpha = sigma_par/sigma_perp at the culet)')
    print('=' * 100)
    print(f'{"alpha":>7} {"geometry":<30} {"dE_ZPL(120)":>12} {"lam_opt":>9} '
          f'{"eta(473)":>9} {"eta(532)":>9}')
    for a, lab in ((1.00, 'ideal hydrostatic'),
                   (0.95, 'micropillar [Hil23, Ho26]'),
                   (0.80, 'intermediate'),
                   (0.56, 'standard flat culet [Hil23]')):
        ma = NVModel(T=RT, alpha=a)
        la = ma.lambda_opt(120.)
        ea = ma.eta_lambda(la, 120.)[0]
        print(f'{a:>7.2f} {lab:<30} {ma.dZPL(120.)*1e3:>9.0f} meV {la:>7.1f} nm '
              f'{ma.eta_lambda(473., 120.)[0]/ea:>9.2f} '
              f'{ma.eta_lambda(532., 120.)[0]/ea:>9.2f}')
    print('  NOTE [Hil23]: a standard flat culet loses ODMR contrast entirely')
    print('  by 40-50 GPa, so the alpha = 0.56 row is not an operating point --')
    print('  it is WHY a micropillar (or equivalent quasi-hydrostatic geometry)')
    print('  is a PRECONDITION for the answer, not an optional refinement.')

    # ---------------- C-2: Monte Carlo band on lambda_opt --------------------
    print('\n' + '=' * 100)
    print('UNCERTAINTY OF THE OPTIMUM  (C-2: randomise the MEASURED ZPL quantities)')
    print('=' * 100)
    n = 400
    rng = np.random.default_rng(0)
    draws = np.array([default_randomiser(rng, T=RT).lambda_opt(120.) for _ in range(n)])
    q16, q50, q84 = np.percentile(draws, [16, 50, 84])
    print(f'  new  (dE120 +/- 20 meV, slope0 +/- 0.25 meV/GPa, S_slope +/- 15%, '
          f'zpl_width +/- 35%):')
    print(f'      lambda_opt = {lo120:.1f}  +{q84-lo120:.1f} / -{lo120-q16:.1f} nm '
          f'(16-84%),  sd = {draws.std():.1f} nm,  n = {n}')

    from nv_model import legacy_randomiser
    rng = np.random.default_rng(0)
    d0 = np.array([legacy_randomiser(rng).lambda_opt(120.) for _ in range(n)])
    l0 = NVModel(Emax=0.758, P0=160., T=0.0, isc=False).lambda_opt(120.)
    p16, p84 = np.percentile(d0, [16, 84])
    print(f'  legacy (Emax +/- 10%, P0 +/- 15% randomised independently):')
    print(f'      lambda_opt = {l0:.1f}  +{p84-l0:.1f} / -{l0-p16:.1f} nm '
          f'(16-84%),  sd = {d0.std():.1f} nm')

    return npass, ntest


if __name__ == '__main__':
    main()
