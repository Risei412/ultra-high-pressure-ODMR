"""
report_figures.py
Figures used by the weekly reports, drawn in the deck palette
(see ``presentaion/slide_style_v3.md``).

These are deliberately narrower than the four-panel exploration figures:
a report figure carries the evidence the report argues from, and nothing else.

  report/image/anvil_culet_decision.png   sensitivity ratio + stress
                                          susceptibility   (report 3)
  report/image/lockin_decision_curve.png  sigma_D vs technical-noise knee
                                          (report 4)

Colour grammar: colour encodes argumentative role, not category.
  Science Blue  = the option under test ([111] culet / modulated detection)
  Deep Indigo   = exceptional structure (the E=0 sector, the break-even knee)
  muted ink     = the reference it is measured against ([100] / DC sweep)
  Mist Lavender = uncertainty band

Run:  python report_figures.py
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

import anvil_orientation as ao
import lockin_sim as ls
import style_v3 as sty

sty.use_style()

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'report', 'image')

BLUE, INDIGO, INK, LAV = (sty.SCIENCE_BLUE, sty.DEEP_INDIGO,
                          sty.NEUTRAL, sty.MIST_LAVENDER)


# ============================================================ report 3 ======
def anvil_culet_decision():
    P, T_DEV, DNU0 = 120.0, 5.0, 5.0
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.3))

    def ratio_best(dt, **kw):
        return ao.eta_ratio(P, T_DEV, dt, dnu0=DNU0, **kw)

    def ratio_axial(dt, **kw):
        lines = ao.odmr_lines(P, T_DEV, '[111]', dt, DNU0, **kw)
        axl = max(lines, key=lambda L: L['dnudt'])       # the E=0 axial line
        e_ax = axl['width'] / (axl['contrast'] * np.sqrt(axl['total_pl']))
        return e_ax / ao.best_line(P, T_DEV, '[100]', dt, DNU0, **kw)[0]

    # ----------------------------------------- left: eta ratio vs dt --------
    ax = axs[0]
    dts = np.linspace(0.0, 8.0, 220)
    rb = np.array([ratio_best(d) for d in dts])
    ra = np.array([ratio_axial(d) for d in dts])

    rng = np.random.default_rng(3)
    mc = np.array([[ratio_best(d, **kw) for d in dts]
                   for kw in (ao.randomiser(rng) for _ in range(160))])
    ax.fill_between(dts, np.percentile(mc, 16, 0), np.percentile(mc, 84, 0),
                    color=LAV, lw=0, label='MC 16-84% (coupling consts.)')
    ax.plot(dts, rb, color=BLUE, lw=2.6, label='[111], best line (off-axis)')
    ax.plot(dts, ra, color=INDIGO, lw=1.8, ls='--',
            label='[111], axial $E{=}0$ line only')
    ax.axhline(1.0, color=INK, lw=1.6)

    dtc = brentq(lambda d: ratio_best(d) - 1.0, 0.1, 20.0)
    ax.axvline(dtc, color=INDIGO, lw=1.0, ls=':')
    ax.plot(dtc, 1.0, 'o', color=INDIGO, ms=6)
    ax.annotate(f'crossover\n$\\delta t$ = {dtc:.2f} GPa',
                (dtc, 1.0), (dtc + 0.7, 1.45), fontsize=9, color=INDIGO,
                arrowprops=dict(arrowstyle='->', color=INDIGO, lw=1))
    ax.set_xlim(0, 8); ax.set_ylim(0.6, 2.3)
    ax.set_xlabel('Deviatoric-stress spread across the culet  $\\delta t$  (GPa)')
    ax.set_ylabel('$\\eta_{[111]}\\,/\\,\\eta_{[100]}$   ($>1$: [100] wins)')
    sty.panel_title(ax, f'Sensitivity ratio vs. stress inhomogeneity '
                        f'($\\Delta\\nu_0$ = {DNU0:.0f} MHz)')
    ax.legend(fontsize=8.5, loc='upper left'); ax.grid(True, alpha=0.35)

    # ----------------------------------------- right: |dnu/dt| per line -----
    ax = axs[1]
    lines111 = ao.odmr_lines(P, T_DEV, '[111]', dt=0.0, dnu0=DNU0)
    axial_nu = max(lines111, key=lambda L: abs(L['dnudt']))['nu']
    labels, vals, cols = [], [], []
    for culet, col in (('[111]', BLUE), ('[100]', INK)):
        for L in ao.odmr_lines(P, T_DEV, culet, dt=0.0, dnu0=DNU0):
            labels.append(f"{culet}\n{L['nu']:.0f} MHz")
            vals.append(abs(L['dnudt']))
            cols.append(INDIGO if (culet == '[111]' and L['nu'] == axial_nu) else col)
    xs = np.arange(len(vals))
    ax.bar(xs, vals, color=cols, width=0.62)
    for x, v in zip(xs, vals):
        ax.text(x, v + 0.14, f'{v:.2f}', ha='center', fontsize=9)
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylim(0, 8.6)
    ax.set_ylabel('$|d\\nu/dt|$  (MHz/GPa)   lower = less stress-broadened')
    sty.panel_title(ax, 'Stress susceptibility of the ODMR transitions')
    ax.text(0.98, 0.96,
            '[100]: $dD/dt = 0$ exactly $\\Rightarrow$\n'
            '$D$ reads mean stress, $E$ reads $t$\n(clean $P$/$t$ decoupling)',
            transform=ax.transAxes, ha='right', va='top', fontsize=8.5, color=INK)
    ax.grid(True, axis='y', alpha=0.35)

    fig.tight_layout()
    out = os.path.join(OUT, 'anvil_culet_decision.png')
    fig.savefig(out, dpi=185, bbox_inches='tight')
    plt.close(fig)
    print(f'saved {out}   (crossover dt = {dtc:.2f} GPa)')


# ============================================================ report 4 ======
def lockin_decision_curve(n_mc=200, seed=1):
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    scheme = {'dc': ('DC sweep', INK, '-', 'o'),
              'onoff': ('MW on-off', BLUE, '--', 's'),
              'fm': ('FM lock-in', BLUE, '-', 'o')}

    knees = np.concatenate([[0.0], np.logspace(0, 3.7, 16)])
    res = {s: [] for s in scheme}
    for fk in knees:
        sd = ls.sigma_D(fk, n_mc=n_mc, seed=seed)
        for s in scheme:
            res[s].append(sd[s])
    kx = np.where(knees == 0, 0.3, knees)      # plot f_knee = 0 at the left edge
    for s, (lbl, c, lsty, mk) in scheme.items():
        ax.loglog(kx, res[s], color=c, lw=2.4, ls=lsty, marker=mk, ms=3.5, label=lbl)

    dc = np.array(res['dc'])
    be = {}
    for s in ('onoff', 'fm'):
        flat = float(np.median(res[s][:6]))
        i = int(np.argmax(dc > flat))
        be[s] = float(np.interp(flat, dc[i - 1:i + 1], kx[i - 1:i + 1]))
        ax.plot(be[s], flat, 'v', color=INDIGO, ms=9, zorder=5)
    ax.annotate(f'break-even\n$f_{{knee}}$ $\\approx$ {be["fm"]:.0f} Hz',
                (be['fm'], np.median(res['fm'][:6])), (0.45, 0.55), fontsize=9,
                color=INDIGO, arrowprops=dict(arrowstyle='->', color=INDIGO, lw=1))
    ax.set_xlim(0.25, 6e3)
    ax.set_xlabel('Technical-noise knee  $f_{knee}$  (Hz)     '
                  '[leftmost point: $f_{knee}=0$]')
    ax.set_ylabel('$\\sigma_D$  (MHz)   lower = better')
    sty.panel_title(ax, 'Line-centre precision vs. technical-noise knee')
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, which='both', alpha=0.35)

    fig.tight_layout()
    out = os.path.join(OUT, 'lockin_decision_curve.png')
    fig.savefig(out, dpi=185, bbox_inches='tight')
    plt.close(fig)
    print(f'saved {out}   (break-even: on-off {be["onoff"]:.0f} Hz, '
          f'FM {be["fm"]:.0f} Hz)')


if __name__ == '__main__':
    anvil_culet_decision()
    lockin_decision_curve()
