"""
talk_figs.py
============
Figures for the 8-minute internal talk, generated from the frozen model in
``code/nv_model.py`` / ``code/nv_model_power.py`` so that every number on a
slide comes from the same source as the manuscript.

These are deliberately NOT the paper figures.  A paper figure is read for a
minute; a talk figure is read in ten seconds from the back of a room.  So each
one carries a single message, has Japanese annotation baked in, and uses large
type and few elements.

Run:
    python talk_figs.py            # all figures, 120 GPa
    python talk_figs.py 4 7        # only the figures for slides 4 and 7
Out:
    slides/figs/fig4_absorption.png   two edges close in on the blue window
    slides/figs/fig5_tornado.png      only the optical inputs move the optimum
    slides/figs/fig6_window.png       the answer is a window, not a point
    slides/figs/fig7_crossover.png    the sign change at 73 GPa
    slides/figs/fig8_power.png        the optimum walks blue with power
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'code'))
from nv_model import NVModel, nm2eV, mc_band, default_randomiser   # noqa: E402
from nv_model_power import (NVModelPower, default_randomiser_power)  # noqa: E402

# --------------------------------------------------------------------------
# house style
# --------------------------------------------------------------------------
NAVY   = '#1e3a6e'
PURPLE = '#6b2d91'
GREY   = '#8a94a6'
DARK   = '#2b3038'
BAND   = '#c9d2e6'

JP = ['IPAGothic', 'IPAPGothic', 'DejaVu Sans']
rcParams.update({
    'font.family': JP,
    'font.size': 19,
    'axes.linewidth': 1.3,
    'axes.edgecolor': DARK,
    'axes.labelcolor': DARK,
    'xtick.color': DARK, 'ytick.color': DARK,
    'xtick.major.width': 1.3, 'ytick.major.width': 1.3,
    'xtick.major.size': 6, 'ytick.major.size': 6,
    'mathtext.fontset': 'dejavusans',
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.06,
})

P0 = 120.0
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs')
DPI = 200


def _frame(ax):
    """Two-sided frame: talk figures read better without a box."""
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.grid(True, color='#e8ebf2', lw=1.0)
    ax.set_axisbelow(True)


def _save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print('wrote', path)
    return path


# --------------------------------------------------------------------------
# slide 4 -- the two edges close in on the blue window
# --------------------------------------------------------------------------
def fig4_absorption():
    m = NVModel()
    lam = np.linspace(380, 700, 1200)
    E = nm2eV(lam)
    s0 = m.sigma_abs(E, 0.0)
    s120 = m.sigma_abs(E, P0)

    ip0 = 1239.841984 / m.IP_A2(0.0)      # ionisation edge, ambient (nm)
    ip120 = 1239.841984 / m.IP_A2(P0)     # ionisation edge, 120 GPa (nm)
    pk0 = lam[s0.argmax()]
    pk120 = lam[s120.argmax()]

    fig, ax = plt.subplots(figsize=(11.2, 6.1))
    _frame(ax)

    ax.fill_between(lam, 0, s0, color=GREY, alpha=0.20, lw=0)
    ax.plot(lam, s0, color=GREY, lw=3.0, label='0 GPa')
    ax.fill_between(lam, 0, s120, color=NAVY, alpha=0.20, lw=0)
    ax.plot(lam, s120, color=NAVY, lw=4.0, label='120 GPa')

    top = float(max(s0.max(), s120.max())) * 1.16
    ax.set_ylim(0, top)
    ax.set_xlim(380, 700)

    # the hard blue wall
    ax.axvspan(380, ip120, color=PURPLE, alpha=0.10, lw=0)
    ax.axvline(ip120, color=PURPLE, lw=2.6)
    ax.text(ip120 + 7, top * 0.55, 'これより青は\nNV$^-$ を直接壊す', color=PURPLE,
            ha='left', va='center', fontsize=17, linespacing=1.3)

    # the two commercial lines
    for x, c, lab in ((473.0, NAVY, '473 nm'), (532.0, GREY, '532 nm')):
        ax.axvline(x, color=c, lw=2.0, ls=(0, (5, 4)))
        ax.text(x + 3, top * 0.995, lab, color=c, ha='left', va='top', fontsize=18)

    # arrow 1: the absorption maximum walks blue
    y1 = top * 0.80
    ax.annotate('', xy=(pk120, y1), xytext=(pk0, y1),
                arrowprops=dict(arrowstyle='-|>,head_width=0.32,head_length=0.7',
                                color=NAVY, lw=2.6, shrinkA=0, shrinkB=0))
    ax.text((pk0 + pk120) / 2, y1 + top * 0.022,
            f'吸収極大  {pk0:.0f} → {pk120:.0f} nm',
            color=NAVY, ha='center', va='bottom', fontsize=19)

    # arrow 2: the ionisation edge walks blue too
    y2 = top * 0.135
    ax.annotate('', xy=(ip120, y2), xytext=(ip0, y2),
                arrowprops=dict(arrowstyle='-|>,head_width=0.32,head_length=0.7',
                                color=PURPLE, lw=2.6, shrinkA=0, shrinkB=0))
    ax.text(ip0 + 6, y2, f'イオン化端  {ip0:.0f} → {ip120:.0f} nm',
            color=PURPLE, ha='left', va='center', fontsize=19)

    ax.set_xlabel('励起波長  $\\lambda$  (nm)')
    ax.set_ylabel('吸収断面積  $\\sigma_{\\mathrm{abs}}$\n(常圧 532 nm = 1)')
    ax.legend(loc='upper right', frameon=False, fontsize=19,
              bbox_to_anchor=(1.0, 0.94))
    return _save(fig, 'fig4_absorption.png')


# --------------------------------------------------------------------------
# slide 5 -- only the optical inputs move the optimum
# --------------------------------------------------------------------------
# (label, kwargs at the low end, kwargs at the high end, "is this an optical
#  input rather than a phenomenological charge-transfer constant?")
TORNADO = [
    ('$\\hbar\\omega$  有効フォノン  ±15%',
     dict(hw=0.065 * 0.85), dict(hw=0.065 * 1.15), True),
    ('$\\Delta E_{\\mathrm{ZPL}}$(120 GPa)  ±20 meV',
     dict(dE120=0.380), dict(dE120=0.420), True),
    ('$S_{\\mathrm{abs}}$ slope  ±15%',
     dict(S_slope=(4.61 - 3.08) * 0.85), dict(S_slope=(4.61 - 3.08) * 1.15), True),
    ('$\\sigma_{\\mathrm{ZPL}}$  0.7–1.4×',
     dict(zpl_width=0.015 * 0.7), dict(zpl_width=0.015 * 1.4), False),
    ('d$E_{\\mathrm{ZPL}}$/d$P|_0$  ±0.25 meV/GPa',
     dict(slope0=5.50e-3), dict(slope0=6.00e-3), False),
    ('$a_{\\mathrm{gs}}$  ±50%',
     dict(a_gs=3.0), dict(a_gs=9.0), False),
    ('$a_{\\mathrm{es}}$  ±50%',
     dict(a_es=0.45), dict(a_es=1.35), False),
    ('$r_0$  ±30%',
     dict(r0=0.84), dict(r0=1.56), False),
    ('$r_{\\mathrm{bg}}$  ±40%',
     dict(rbg=0.09), dict(rbg=0.21), False),
    ('$w_0$  ±30%',
     dict(w0=0.7), dict(w0=1.3), False),
    ('$C_{\\mathrm{amb}}$  ±50%',
     dict(C_amb=0.1235), dict(C_amb=0.3704), False),
    ('$E_{\\mathrm{ISC}}$  ±30%',
     dict(E_isc=0.1265), dict(E_isc=0.2349), False),
]


def _tornado_rows():
    ref = NVModel().lambda_opt(P0)
    rows = []
    for label, lo_kw, hi_kw, optical in TORNADO:
        lo = NVModel(**lo_kw).lambda_opt(P0) - ref
        hi = NVModel(**hi_kw).lambda_opt(P0) - ref
        rows.append((label, lo, hi, optical))
    rows.sort(key=lambda r: max(abs(r[1]), abs(r[2])), reverse=True)
    return ref, rows


def fig5_tornado():
    ref, rows = _tornado_rows()
    n = len(rows)
    fig, ax = plt.subplots(figsize=(12.6, 6.6))

    ys = np.arange(n)[::-1]
    span = max(max(abs(r[1]), abs(r[2])) for r in rows)
    for y, (label, lo, hi, optical) in zip(ys, rows):
        col = NAVY if optical else GREY
        neg, pos = min(lo, hi), max(lo, hi)
        if max(abs(neg), abs(pos)) < 0.02:
            ax.plot([0], [y], marker='o', ms=9, color=GREY, zorder=3)
            continue
        ax.barh(y, neg, color=col, height=0.62, zorder=2)
        ax.barh(y, pos, color=col, height=0.62, zorder=2)
        ax.text(pos + 0.30, y, f'±{max(abs(neg), abs(pos)):.1f} nm',
                va='center', ha='left', color=col, fontsize=19)

    n_optical = sum(1 for r in rows if r[3])
    ax.axhline(ys[n_optical - 1] - 0.5, color='#c8cede', lw=1.4)
    ax.axvline(0, color=DARK, lw=2.0, zorder=4)

    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=18)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xlabel(f'120 GPa における $\\lambda_{{\\mathrm{{opt}}}}$ の変化  (nm)'
                  f'    [基準 {ref:.1f} nm]')
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.grid(True, axis='x', color='#e8ebf2', lw=1.0)
    ax.set_axisbelow(True)
    ax.set_xlim(-span - 1.0, span + 3.6)

    ax.text(span * 0.30, ys[n_optical] - 0.7, 'すべて 0.00 nm',
            color=GREY, fontsize=20, va='center')
    handles = [plt.Rectangle((0, 0), 1, 1, color=NAVY),
               plt.Rectangle((0, 0), 1, 1, color=GREY)]
    ax.legend(handles, ['光学入力 — 答えを動かす', '現象論定数 — 動かさない'],
              loc='lower right', frameon=True, framealpha=1.0,
              edgecolor='#dde2ec', fontsize=18)
    return _save(fig, 'fig5_tornado.png')


# --------------------------------------------------------------------------
# slide 6 -- the answer is a window, not a point
# --------------------------------------------------------------------------
def fig6_window(n_mc=200):
    m = NVModel()
    lam = np.arange(400.0, 561.0, 0.25)
    eta = np.asarray(m.eta_lambda(lam, P0)[0])
    lopt = m.lambda_opt(P0)
    eopt = float(np.asarray(m.eta_lambda(lopt, P0)[0]))
    rel = eta / eopt

    # 5% tolerance window
    inside = np.flatnonzero(rel <= 1.05)
    lo_w, hi_w = lam[inside[0]], lam[inside[-1]]

    # MC band on the curve, and on the optimum itself
    lo_b, hi_b = mc_band(default_randomiser,
                         lambda mm: (np.asarray(mm.eta_lambda(lam, P0)[0])
                                     / float(np.asarray(
                                         mm.eta_lambda(mm.lambda_opt(P0), P0)[0]))),
                         n=n_mc, seed=3)
    rng = np.random.default_rng(11)
    opts = np.array([default_randomiser(rng).lambda_opt(P0) for _ in range(n_mc)])
    q16, q84 = np.percentile(opts, [16, 84])

    fig, ax = plt.subplots(figsize=(11.6, 6.2))
    _frame(ax)
    ax.set_ylim(0.9, 4.2)
    ax.set_xlim(400, 560)

    ax.axvspan(lo_w, hi_w, color=BAND, alpha=0.55, lw=0)
    ax.fill_between(lam, lo_b, hi_b, color=NAVY, alpha=0.13, lw=0)
    ax.plot(lam, rel, color=NAVY, lw=4.2, zorder=3)
    ax.axhline(1.05, color=NAVY, lw=1.4, ls=':')

    # the MC interval of the optimum itself: stated, not drawn on the crowded floor
    ax.text(494, 4.05,
            f'$\\lambda_{{\\mathrm{{opt}}}}$ = {lopt:.1f} nm\n'
            f'MC 16–84%: $^{{+{q84 - lopt:.1f}}}_{{-{lopt - q16:.1f}}}$ nm',
            color=PURPLE, ha='left', va='top', fontsize=18, linespacing=1.5)

    ax.text((lo_w + hi_w) / 2, 4.12, f'5% 許容窓\n{lo_w:.0f}–{hi_w:.0f} nm',
            color=NAVY, ha='center', va='top', fontsize=19, linespacing=1.2)

    for x, lab, xt, yt, ha in ((405.0, '405 nm', 410, 3.55, 'left'),
                               (532.0, '532 nm', 528, 2.95, 'right')):
        y = float(np.asarray(m.eta_lambda(x, P0)[0])) / eopt
        ax.plot([x], [y], marker='o', ms=11, color=GREY, zorder=5)
        ax.axvline(x, color=GREY, lw=1.6, ls=(0, (5, 4)))
        ax.text(xt, yt, f'{lab}\n×{y:.1f}', color=GREY,
                ha=ha, va='center', fontsize=19, linespacing=1.2)

    y473 = float(np.asarray(m.eta_lambda(473.0, P0)[0])) / eopt
    ax.plot([473.0], [y473], marker='o', ms=13, color=PURPLE, zorder=6)
    ax.annotate(f'473 nm DPSS\n最適の {100 * (y473 - 1):.1f}% 落ち',
                xy=(473.0, y473), xytext=(414, 1.62),
                color=PURPLE, fontsize=19, ha='left', linespacing=1.2,
                arrowprops=dict(arrowstyle='-', color=PURPLE, lw=1.8,
                                shrinkA=2, shrinkB=6))

    ax.set_xlabel('励起波長  $\\lambda$  (nm)')
    ax.set_ylabel('ロックイン感度  $\\eta/\\eta_{\\mathrm{opt}}$\n(小さいほど良い)')
    return _save(fig, 'fig6_window.png')


# --------------------------------------------------------------------------
# slide 7 -- the sign change near 73 GPa
# --------------------------------------------------------------------------
def _ratio(m, P, blue=473.0):
    return (float(np.asarray(m.eta_lambda(532.0, P)[0]))
            / float(np.asarray(m.eta_lambda(blue, P)[0])))


def fig7_crossover(n_mc=120):
    m = NVModel()
    P = np.linspace(15, 150, 400)
    r = np.array([_ratio(m, p) for p in P])
    xover = brentq(lambda p: _ratio(m, p) - 1.0, 20., 145.)

    lo_b, hi_b = mc_band(default_randomiser,
                         lambda mm: np.array([_ratio(mm, p) for p in P]),
                         n=n_mc, seed=5)

    fig, ax = plt.subplots(figsize=(11.6, 6.2))
    _frame(ax)
    ax.set_yscale('log')
    ax.set_xlim(15, 150)
    ax.set_ylim(0.22, 5.0)

    ax.axhspan(0.22, 1.0, color=BAND, alpha=0.45, lw=0)
    ax.fill_between(P, lo_b, hi_b, color=NAVY, alpha=0.13, lw=0)
    ax.plot(P, r, color=NAVY, lw=4.2, zorder=3)
    ax.axhline(1.0, color=DARK, lw=1.8)

    ax.axvline(xover, color=PURPLE, lw=2.0, ls=':')
    ax.plot([xover], [1.0], marker='o', ms=13, color=PURPLE, zorder=5)
    ax.annotate(f'符号が変わる  {xover:.0f} GPa',
                xy=(xover, 1.0), xytext=(xover + 9, 0.42),
                color=PURPLE, fontsize=20, ha='left',
                arrowprops=dict(arrowstyle='-|>,head_width=0.28,head_length=0.6',
                                color=PURPLE, lw=1.8))

    r50 = _ratio(m, 50.0)
    ax.plot([50.0], [r50], marker='s', ms=13, color=GREY, zorder=5)
    ax.annotate(f'50 GPa 実測: 優位なし [Bha22]\nモデルは {r50:.2f} — 緑が有利',
                xy=(50.0, r50), xytext=(20, 2.3),
                color=GREY, fontsize=19, ha='left', linespacing=1.25,
                bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'),
                arrowprops=dict(arrowstyle='-|>,head_width=0.28,head_length=0.6',
                                color=GREY, lw=1.8))

    r120 = _ratio(m, 120.0)
    ax.plot([120.0], [r120], marker='o', ms=13, color=NAVY, zorder=5)
    ax.text(117.0, r120 * 1.16, f'120 GPa  ×{r120:.1f}', color=NAVY,
            ha='right', va='bottom', fontsize=20)

    ax.set_yticks([0.25, 0.5, 1, 2, 4])
    ax.set_yticklabels(['0.25', '0.5', '1', '2', '4'])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel('圧力  $P$  (GPa)')
    ax.set_ylabel('$\\eta$(532 nm) / $\\eta$(473 nm)\n(>1 なら青が有利)')
    return _save(fig, 'fig7_crossover.png')


# --------------------------------------------------------------------------
# slide 8 -- the optimum walks blue with power (uncalibrated)
# --------------------------------------------------------------------------
def fig8_power(n_mc=60):
    m = NVModelPower()
    lam = np.linspace(402, 560, 700)
    u = np.logspace(-2.3, 1.0, 90)

    def ridge(mm):
        e = np.array([np.asarray(mm.eta_lambda_u(lam, P0, uu)[0]) for uu in u])
        return lam[np.nanargmin(e, axis=1)]

    lo_r = ridge(m)
    lo_b, hi_b = mc_band(default_randomiser_power, ridge, n=n_mc, seed=11)

    fig, ax = plt.subplots(figsize=(11.6, 6.2))
    _frame(ax)
    ax.set_xscale('log')
    ax.set_xlim(u[0], u[-1])
    ax.set_ylim(398, 500)

    # the regime in which the fixed-power recommendation is the answer
    ax.axvspan(u[0], 0.1, color=BAND, alpha=0.5, lw=0)
    ax.text(np.sqrt(u[0] * 0.1), 494, '低励起  $u \\lesssim 0.1$\n473 nm が最適',
            color=NAVY, ha='center', va='top', fontsize=19, linespacing=1.25)

    ax.fill_between(u, lo_b, hi_b, color=PURPLE, alpha=0.16, lw=0)
    ax.plot(u, lo_r, color=PURPLE, lw=4.2, zorder=3)

    ax.axhline(473.0, color=NAVY, lw=2.0, ls=(0, (5, 4)))
    ax.text(u[-1], 475, '473 nm', color=NAVY, ha='right', va='bottom', fontsize=19)
    ax.axhline(405.2, color=DARK, lw=2.0, ls=':')
    ax.text(u[0] * 1.15, 406.5, 'イオン化端 405 nm — 稜線はここで止まる',
            color=DARK, ha='left', va='bottom', fontsize=18)

    # where the existing high-pressure literature sits
    ax.axvspan(0.1, 0.3, color='#f2dede', alpha=0.6, lw=0)
    ax.text(np.sqrt(0.1 * 0.3), 494, '既存実験は\nこの範囲 [Dai22]',
            color='#a03b3b', ha='center', va='top', fontsize=18, linespacing=1.25)

    mp = NVModelPower()
    anno = {0.1: (0.092, 461, 'right'), 0.3: (0.27, 436, 'right')}
    for uu in (0.1, 0.3):
        e = np.asarray(mp.eta_lambda_u(lam, P0, uu)[0])
        best = lam[int(np.nanargmin(e))]
        pen = (float(np.asarray(mp.eta_lambda_u(473.0, P0, uu)[0]))
               / float(np.nanmin(e)))
        xt, yt, ha = anno[uu]
        ax.plot([uu], [best], marker='o', ms=11, color=PURPLE, zorder=5)
        ax.text(xt, yt, f'$u$={uu}: {best:.0f} nm\n473 固定は ×{pen:.2f}',
                color=PURPLE, ha=ha, va='top', fontsize=18, linespacing=1.25)

    ax.set_xlabel('規格化強度  $u = I/I_{1/2}$   ($u$=1 で NV$^-$ 遷移が半飽和)')
    ax.set_ylabel('最適励起波長  (nm)')
    ax.set_title('未校正 — 向きは全ドローで一致、量は校正前のシナリオ',
                 fontsize=19, color='#a03b3b', pad=12)
    return _save(fig, 'fig8_power.png')


FIGS = {4: fig4_absorption, 5: fig5_tornado, 6: fig6_window,
        7: fig7_crossover, 8: fig8_power}

if __name__ == '__main__':
    want = [int(a) for a in sys.argv[1:]] or sorted(FIGS)
    for k in want:
        FIGS[k]()
