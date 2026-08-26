"""
talk_figs.py
============
Figures for the 8-minute internal talk, generated from the frozen model in
``code/nv_model.py`` / ``code/nv_model_power.py`` so that every number on a
slide comes from the same source as the manuscript.

These are deliberately BARE: axes, ticks, axis labels and the data, and nothing
else.  No call-outs, no reference lines, no shaded regions, no legends, no
titles.  Everything that explains a plot is added in PowerPoint afterwards,
where it can be moved and reworded without regenerating anything.  (The earlier
annotated versions are in git history, at commit dc17a85.)

Each figure is drawn at exactly the size it occupies on the slide -- FIGW
inches wide, with no tight-bbox trim -- so an 18 pt tick label here is 18 pt
when projected, and any text added in PowerPoint at 18 pt will match it.

Run:
    python talk_figs.py            # all figures, 120 GPa
    python talk_figs.py 4 7        # only the figures for slides 4 and 7
Out:
    slides/figs/fig4_absorption.png   sigma_abs at 0 and 120 GPa
    slides/figs/fig5_tornado.png      d(lambda_opt) for each input
    slides/figs/fig6_window.png       eta / eta_opt vs wavelength
    slides/figs/fig7_crossover.png    eta(532) / eta(473) vs pressure
    slides/figs/fig8_power.png        optimal wavelength vs intensity
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'code'))
from nv_model import NVModel, nm2eV, mc_band, default_randomiser   # noqa: E402
from nv_model_power import (NVModelPower, default_randomiser_power)  # noqa: E402

# --------------------------------------------------------------------------
# house style -- Institute of Science Tokyo navy, one neutral grey, nothing else
# --------------------------------------------------------------------------
NAVY = '#1c3177'
GREY = '#8a94a6'
DARK = '#2b3038'

JP = ['IPAGothic', 'IPAPGothic', 'DejaVu Sans']
BASE = 18.0                 # the floor: nothing on a slide may be smaller
FIGW = 11.60                # inches; the figure spans the slide
FIGH = 4.55
rcParams.update({
    'font.family': JP,
    'font.size': BASE,
    'axes.linewidth': 1.3,
    'axes.edgecolor': DARK,
    'axes.labelcolor': DARK,
    'xtick.color': DARK, 'ytick.color': DARK,
    'xtick.major.width': 1.3, 'ytick.major.width': 1.3,
    'xtick.major.size': 6, 'ytick.major.size': 6,
    'mathtext.fontset': 'dejavusans',
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
})

P0 = 120.0
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs')
DPI = 200


def _axes():
    ax = plt.subplots(figsize=(FIGW, FIGH), layout='constrained')[1]
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.grid(True, color='#e8ebf2', lw=1.0)
    ax.set_axisbelow(True)
    return ax


def _save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    # no bbox='tight': the figure must keep the size it was created at, or the
    # deck rescales it and the point sizes stop matching the slide
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print('wrote', path)
    return path


# --------------------------------------------------------------------------
# slide 4 -- the absorption envelope at ambient pressure and at 120 GPa
# --------------------------------------------------------------------------
def fig4_absorption():
    m = NVModel()
    lam = np.linspace(380, 700, 1200)
    E = nm2eV(lam)

    ax = _axes()
    for P, colour, lw in ((0.0, GREY, 3.0), (P0, NAVY, 4.0)):
        s = m.sigma_abs(E, P)
        ax.fill_between(lam, 0, s, color=colour, alpha=0.20, lw=0)
        ax.plot(lam, s, color=colour, lw=lw)

    ax.set_xlim(380, 700)
    ax.set_ylim(0, None)
    ax.set_xlabel('励起波長  $\\lambda$  (nm)')
    ax.set_ylabel('$\\sigma_{\\mathrm{abs}}$  (常圧 532 nm = 1)')
    return _save(ax.figure, 'fig4_absorption.png')


# --------------------------------------------------------------------------
# slide 5 -- how far each input moves the optimum
# --------------------------------------------------------------------------
# (tick label, kwargs at the low end, kwargs at the high end, optical input?)
TORNADO = [
    ('$\\hbar\\omega$  有効フォノン  ±15%',
     dict(hw=0.065 * 0.85), dict(hw=0.065 * 1.15), True),
    ('$\\Delta E_{\\mathrm{ZPL}}$(120)  ±20 meV',
     dict(dE120=0.380), dict(dE120=0.420), True),
    ('$S_{\\mathrm{abs}}$ 勾配  ±15%',
     dict(S_slope=(4.61 - 3.08) * 0.85), dict(S_slope=(4.61 - 3.08) * 1.15), True),
    ('$\\sigma_{\\mathrm{ZPL}}$  0.7–1.4×',
     dict(zpl_width=0.015 * 0.7), dict(zpl_width=0.015 * 1.4), False),
    ('d$E_{\\mathrm{ZPL}}$/d$P|_0$  ±0.25',
     dict(slope0=5.50e-3), dict(slope0=6.00e-3), False),
    ('$a_{\\mathrm{gs}}$  ±50%', dict(a_gs=3.0), dict(a_gs=9.0), False),
    ('$a_{\\mathrm{es}}$  ±50%', dict(a_es=0.45), dict(a_es=1.35), False),
    ('$r_0$  ±30%', dict(r0=0.84), dict(r0=1.56), False),
    ('$r_{\\mathrm{bg}}$  ±40%', dict(rbg=0.09), dict(rbg=0.21), False),
    ('$w_0$  ±30%', dict(w0=0.7), dict(w0=1.3), False),
    ('$C_{\\mathrm{amb}}$  ±50%', dict(C_amb=0.1235), dict(C_amb=0.3704), False),
    ('$E_{\\mathrm{ISC}}$  ±30%', dict(E_isc=0.1265), dict(E_isc=0.2349), False),
]


def fig5_tornado():
    ref = NVModel().lambda_opt(P0)
    rows = []
    for label, lo_kw, hi_kw, optical in TORNADO:
        lo = NVModel(**lo_kw).lambda_opt(P0) - ref
        hi = NVModel(**hi_kw).lambda_opt(P0) - ref
        rows.append((label, min(lo, hi), max(lo, hi), optical))
    rows.sort(key=lambda r: max(abs(r[1]), abs(r[2])), reverse=True)

    ax = _axes()
    ys = np.arange(len(rows))[::-1]
    for y, (_, neg, pos, optical) in zip(ys, rows):
        colour = NAVY if optical else GREY
        if max(abs(neg), abs(pos)) < 0.02:
            ax.plot([0], [y], marker='o', ms=9, color=GREY, zorder=3)
            continue
        ax.barh(y, neg, color=colour, height=0.62, zorder=2)
        ax.barh(y, pos, color=colour, height=0.62, zorder=2)

    ax.axvline(0, color=DARK, lw=2.0, zorder=4)
    span = max(max(abs(r[1]), abs(r[2])) for r in rows)
    ax.set_xlim(-span - 1.0, span + 1.0)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=BASE)
    ax.spines['left'].set_visible(False)
    ax.grid(False, axis='y')
    ax.set_xlabel(f'$\\lambda_{{\\mathrm{{opt}}}}$ の変化  (nm)   '
                  f'[基準 {ref:.1f} nm]')
    return _save(ax.figure, 'fig5_tornado.png')


# --------------------------------------------------------------------------
# slide 6 -- sensitivity against excitation wavelength
# --------------------------------------------------------------------------
def fig6_window(n_mc=200):
    m = NVModel()
    lam = np.arange(400.0, 561.0, 0.25)
    eopt = float(np.asarray(m.eta_lambda(m.lambda_opt(P0), P0)[0]))
    rel = np.asarray(m.eta_lambda(lam, P0)[0]) / eopt

    def norm(mm):
        e0 = float(np.asarray(mm.eta_lambda(mm.lambda_opt(P0), P0)[0]))
        return np.asarray(mm.eta_lambda(lam, P0)[0]) / e0
    lo, hi = mc_band(default_randomiser, norm, n=n_mc, seed=3)

    ax = _axes()
    ax.fill_between(lam, lo, hi, color=NAVY, alpha=0.13, lw=0)
    ax.plot(lam, rel, color=NAVY, lw=4.2)
    ax.set_xlim(400, 560)
    ax.set_ylim(0.95, 4.2)
    ax.set_xlabel('励起波長  $\\lambda$  (nm)')
    ax.set_ylabel('$\\eta/\\eta_{\\mathrm{opt}}$')
    return _save(ax.figure, 'fig6_window.png')


# --------------------------------------------------------------------------
# slide 7 -- green against blue, as a function of pressure
# --------------------------------------------------------------------------
def _ratio(m, P, blue=473.0):
    return (float(np.asarray(m.eta_lambda(532.0, P)[0]))
            / float(np.asarray(m.eta_lambda(blue, P)[0])))


def fig7_crossover(n_mc=120):
    m = NVModel()
    P = np.linspace(15, 150, 400)
    r = np.array([_ratio(m, p) for p in P])
    lo, hi = mc_band(default_randomiser,
                     lambda mm: np.array([_ratio(mm, p) for p in P]),
                     n=n_mc, seed=5)

    ax = _axes()
    ax.set_yscale('log')
    ax.fill_between(P, lo, hi, color=NAVY, alpha=0.13, lw=0)
    ax.plot(P, r, color=NAVY, lw=4.2)
    ax.set_xlim(15, 150)
    ax.set_ylim(0.22, 5.0)
    ax.set_yticks([0.25, 0.5, 1, 2, 4])
    ax.set_yticklabels(['0.25', '0.5', '1', '2', '4'])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel('圧力  $P$  (GPa)')
    ax.set_ylabel('$\\eta$(532) / $\\eta$(473)')
    return _save(ax.figure, 'fig7_crossover.png')


# --------------------------------------------------------------------------
# slide 8 -- the optimal wavelength against excitation intensity
# --------------------------------------------------------------------------
def fig8_power(n_mc=60):
    m = NVModelPower()
    lam = np.linspace(402, 560, 700)
    u = np.logspace(-2.3, 1.0, 90)

    def ridge(mm):
        e = np.array([np.asarray(mm.eta_lambda_u(lam, P0, uu)[0]) for uu in u])
        return lam[np.nanargmin(e, axis=1)]
    lo, hi = mc_band(default_randomiser_power, ridge, n=n_mc, seed=11)

    ax = _axes()
    ax.set_xscale('log')
    ax.fill_between(u, lo, hi, color=NAVY, alpha=0.16, lw=0)
    ax.plot(u, ridge(m), color=NAVY, lw=4.2)
    ax.set_xlim(u[0], u[-1])
    ax.set_ylim(398, 500)
    ax.set_xlabel('規格化強度  $u = I/I_{1/2}$')
    ax.set_ylabel('最適励起波長  (nm)')
    return _save(ax.figure, 'fig8_power.png')


FIGS = {4: fig4_absorption, 5: fig5_tornado, 6: fig6_window,
        7: fig7_crossover, 8: fig8_power}

if __name__ == '__main__':
    want = [int(a) for a in sys.argv[1:]] or sorted(FIGS)
    for k in want:
        FIGS[k]()
