"""
fig6_anvil_orientation.py
[111] vs [100] culet orientation for CW-ODMR at 120 GPa.

  (a) zero-field ODMR line structure (stick spectrum, height = contrast)
  (b) simulated spectra for two deviatoric-stress spreads dt
  (c) eta([111])/eta([100]) vs dt: best line, axial line only, MC band
  (d) sensitivity of each line to deviatoric stress, |dnu/dt|, and the
      P/t decoupling that [100] provides

Run:  python fig6_anvil_orientation.py
Out:  anvil_orientation_120GPa.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.optimize import brentq

import anvil_orientation as ao

rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                 'axes.linewidth': 0.9, 'mathtext.fontset': 'dejavusans'})

P, T_DEV, DNU0 = 120.0, 5.0, 5.0
c111, c100, cMC = '#c1272d', '#2660c4', '#c1272d'

fig, axs = plt.subplots(2, 2, figsize=(12.6, 8.6))


def lorentz(nu, nu0, w):
    return 1.0 / (1.0 + ((nu - nu0) / (w / 2.0)) ** 2)


# --------------------------------------------------- (a) line structure ----
ax = axs[0, 0]
for culet, col, off in (('[111]', c111, +0.0), ('[100]', c100, -0.0)):
    lines = ao.odmr_lines(P, T_DEV, culet, dt=0.0, dnu0=DNU0)
    for L in lines:
        ax.vlines(L['nu'], 0, L['contrast'], color=col, lw=3.0,
                  alpha=0.9 if culet == '[111]' else 0.7)
fam111 = ao.families(P, T_DEV, '[111]')
ax.annotate('axial NV\n$E=0$ exactly\n(contrast 0.375)',
            (fam111[0]['D'], 0.375), (fam111[0]['D'] - 12, 0.47),
            color=c111, fontsize=8.5, ha='center',
            arrowprops=dict(arrowstyle='->', color=c111, lw=1))
ax.annotate('3 off-axis NV,\nsplit by $E$', (1131.4, 0.312), (1131.4 - 34, 0.40),
            color=c111, fontsize=8.5, ha='center',
            arrowprops=dict(arrowstyle='->', color=c111, lw=1))
ax.annotate('all 4 NV coincide,\nsplit by $E=2|b|t$', (1143.4, 0.50), (1163, 0.44),
            color=c100, fontsize=8.5, ha='left',
            arrowprops=dict(arrowstyle='->', color=c100, lw=1))
ax.plot([], [], color=c111, lw=3, label='[111] culet')
ax.plot([], [], color=c100, lw=3, label='[100] culet')
ax.set_xlim(1070, 1190); ax.set_ylim(0, 0.58)
ax.set_xlabel('ODMR frequency  (MHz)')
ax.set_ylabel('Contrast weight (fraction of total PL)')
ax.set_title(f'(a)  Line structure at {P:.0f} GPa, $t$ = {T_DEV:.0f} GPa',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=9, loc='upper left'); ax.grid(True, alpha=0.18)

# --------------------------------------------------- (b) simulated spectra --
ax = axs[0, 1]
nu = np.linspace(1060, 1200, 2400)
for dt, ls in ((1.0, '-'), (5.0, '--')):
    for culet, col in (('[111]', c111), ('[100]', c100)):
        y = np.ones_like(nu)
        for L in ao.odmr_lines(P, T_DEV, culet, dt=dt, dnu0=DNU0):
            y -= L['contrast'] * lorentz(nu, L['nu'], L['width'])
        ax.plot(nu, y, color=col, lw=2.0, ls=ls, alpha=0.95 if dt == 1 else 0.55,
                label=f'{culet}, $\\delta t$={dt:.0f} GPa')
ax.set_xlim(1060, 1200); ax.set_ylim(0.35, 1.03)
ax.set_xlabel('ODMR frequency  (MHz)'); ax.set_ylabel('Normalised PL')
ax.set_title('(b)  Deviatoric-stress inhomogeneity washes the lines out',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='lower left'); ax.grid(True, alpha=0.18)

# --------------------------------------------------- (c) eta ratio vs dt ---
ax = axs[1, 0]
dts = np.linspace(0.0, 8.0, 220)


def ratio_best(dt, **kw):
    return ao.eta_ratio(P, T_DEV, dt, dnu0=DNU0, **kw)


def ratio_axial(dt, **kw):
    lines = ao.odmr_lines(P, T_DEV, '[111]', dt, DNU0, **kw)
    axl = max(lines, key=lambda L: L['dnudt'])          # the E=0 axial line
    e_ax = axl['width'] / (axl['contrast'] * np.sqrt(axl['total_pl']))
    return e_ax / ao.best_line(P, T_DEV, '[100]', dt, DNU0, **kw)[0]


rb = np.array([ratio_best(d) for d in dts])
ra = np.array([ratio_axial(d) for d in dts])

rng = np.random.default_rng(3)
mc = np.array([[ratio_best(d, **kw) for d in dts]
               for kw in (ao.randomiser(rng) for _ in range(160))])
ax.fill_between(dts, np.percentile(mc, 16, 0), np.percentile(mc, 84, 0),
                color=cMC, alpha=0.20, lw=0, label='MC 16-84% (coupling consts.)')
ax.plot(dts, rb, color=c111, lw=2.6, label='[111], best line (off-axis)')
ax.plot(dts, ra, color=c111, lw=1.8, ls='--', label='[111], axial $E{=}0$ line only')
ax.axhline(1.0, color=c100, lw=1.6)
dtc = brentq(lambda d: ratio_best(d) - 1.0, 0.1, 20.0)
ax.axvline(dtc, color='0.5', lw=1.0, ls=':')
ax.plot(dtc, 1.0, 'o', color='k', ms=6)
ax.annotate(f'crossover\n$\\delta t$ = {dtc:.2f} GPa\n($=0.42\\,\\Delta\\nu_0$ per MHz)',
            (dtc, 1.0), (dtc + 0.7, 1.42), fontsize=8.5,
            arrowprops=dict(arrowstyle='->', color='k', lw=1))
ax.text(7.8, 0.90, '[111] better', color=c111, fontsize=9, ha='right', va='top')
ax.text(7.8, 1.06, '[100] better', color=c100, fontsize=9, ha='right', va='bottom')
ax.set_xlim(0, 8); ax.set_ylim(0.6, 2.3)
ax.set_xlabel('Deviatoric-stress spread across the culet  $\\delta t$  (GPa)')
ax.set_ylabel('$\\eta_{[111]}\\,/\\,\\eta_{[100]}$   ($>1$: [100] wins)')
ax.set_title(f'(c)  Which culet wins  ($\\Delta\\nu_0$ = {DNU0:.0f} MHz)',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='upper left'); ax.grid(True, alpha=0.18)

# --------------------------------------------------- (d) dnu/dt per line ---
ax = axs[1, 1]
labels, vals, cols = [], [], []
for culet, col in (('[111]', c111), ('[100]', c100)):
    for L in ao.odmr_lines(P, T_DEV, culet, dt=0.0, dnu0=DNU0):
        labels.append(f"{culet}\n{L['nu']:.0f} MHz")
        vals.append(abs(L['dnudt'])); cols.append(col)
xs = np.arange(len(vals))
bars = ax.bar(xs, vals, color=cols, alpha=0.85, width=0.62)
for x, v in zip(xs, vals):
    ax.text(x, v + 0.12, f'{v:.2f}', ha='center', fontsize=9)
ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylim(0, 8.6)
ax.set_ylabel('$|d\\nu/dt|$  (MHz/GPa)   lower = less stress-broadened')
ax.set_title('(d)  The $E{=}0$ line is the MOST stress-sensitive one',
             loc='left', weight='bold', fontsize=11)
ax.text(0.98, 0.96,
        '[100]: $dD/dt = 0$ exactly $\\Rightarrow$\n'
        '$D$ reads mean stress, $E$ reads $t$\n(clean $P$/$t$ decoupling)',
        transform=ax.transAxes, ha='right', va='top', fontsize=8.5, color=c100)
ax.grid(True, axis='y', alpha=0.18)

plt.tight_layout()
out = f'anvil_orientation_{int(P)}GPa.png'
plt.savefig(out, dpi=185, bbox_inches='tight')

# ------------------------------------------------------------- summary -----
print(f'saved {out}')
for culet in ('[111]', '[100]'):
    print(f'--- {culet} ---')
    for L in ao.odmr_lines(P, T_DEV, culet, dt=0.0, dnu0=DNU0):
        print(f"   nu={L['nu']:8.1f} MHz  contrast={L['contrast']:.3f}  "
              f"|dnu/dt|={abs(L['dnudt']):5.2f} MHz/GPa")
print(f'\neta ratio at dt=0      : {ratio_best(0.0):.3f}   (axial line: {ratio_axial(0.0):.3f})')
print(f'eta ratio at dt=1 GPa  : {ratio_best(1.0):.3f}   (axial line: {ratio_axial(1.0):.3f})')
print(f'eta ratio at dt=5 GPa  : {ratio_best(5.0):.3f}   (axial line: {ratio_axial(5.0):.3f})')
print(f'crossover dt           : {dtc:.2f} GPa  = 0.42 * dnu0[MHz]')
print(f'axial line never wins  : min over dt = {ra.min():.3f}')
