"""
fig5_threshold.py
-----------------
The blue advantage has a threshold pressure (Sec. IV G of the manuscript, where
it appears as Table tab:threshold).  Also TALK FIGURE 4 -- the falsifiable
prediction.

eta(532 nm)/eta(473 nm) at the same pressure and the same excitation: above 1
blue wins, below 1 green wins.  The model crosses at ~71 GPa, so it predicts
that a comparison made at 50 GPa finds NO advantage for blue -- which is what
Ref. [Bha22] Sec. 6.3 reports.  A sign change is the one signature a
wavelength-dependent systematic of a single-pressure comparison cannot produce
or displace, which is why this curve, and not the 474 nm number, is the
measurement that would falsify this work.

Run:  python fig5_threshold.py
Out:  threshold_green_blue.png
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
from scipy.optimize import brentq

import talk_style as ts
from nv_model import NVModel

ts.use()

m = NVModel(T=300.0)


def ratio(P):
    return float(np.asarray(m.eta_lambda(532., P)[0])
                 / np.asarray(m.eta_lambda(473., P)[0]))


P = np.linspace(15., 150., 271)
r = np.array([ratio(p) for p in P])
xo = brentq(lambda p: ratio(p) - 1.0, 20., 145.)

fig, ax = plt.subplots(figsize=(10.8, 5.9))

ax.axhspan(1.0, 6.0, color=ts.W473, alpha=0.07, lw=0)
ax.axhspan(0.2, 1.0, color=ts.W532, alpha=0.07, lw=0)
ax.axhline(1.0, color=ts.INK2, lw=1.6)
ax.plot(P, r, color=ts.ACCENT, lw=3.0, zorder=4)

ax.axvline(xo, color=ts.INK2, lw=1.2, ls=':')
ax.plot(xo, 1.0, 'o', color=ts.INK, ms=11, mec='white', mew=2, zorder=6)
ax.annotate(ts.t(f'符号が変わる  {xo:.0f} GPa', f'sign change  {xo:.0f} GPa'),
            xy=(xo, 1.0), xytext=(xo + 8, 0.47), color=ts.INK, fontsize=14,
            weight='bold',
            arrowprops=dict(arrowstyle='->', color=ts.INK, lw=1.6))

# the published null result the model has to reproduce
ax.plot(50., ratio(50.), 's', color=ts.W532, ms=12, mec='white', mew=1.8,
        zorder=6)
ax.annotate(ts.t('50 GPa の実測:\n「青の明確な優位なし」\n[Bhattacharyya 2022]',
                 'measured at 50 GPa:\n"no distinct advantage"\n[Bhattacharyya 2022]'),
            xy=(50., ratio(50.)), xytext=(20., 1.35), color=ts.W532,
            fontsize=12.5,
            arrowprops=dict(arrowstyle='->', color=ts.W532, lw=1.5))

# the target
ax.plot(120., ratio(120.), 'o', color=ts.W473, ms=12, mec='white', mew=1.8,
        zorder=6)
ax.annotate(ts.t(f'120 GPa: 青が ×{ratio(120.):.1f}',
                 f'120 GPa: blue by x{ratio(120.):.1f}'),
            xy=(120., ratio(120.)), xytext=(88., 3.9), color=ts.W473,
            fontsize=13.5, weight='bold',
            arrowprops=dict(arrowstyle='->', color=ts.W473, lw=1.6))

ax.text(147, 1.15, ts.t('青(473 nm)が有利', 'blue (473 nm) wins'),
        color=ts.W473, fontsize=13, ha='right', va='bottom', weight='bold')
ax.text(147, 0.87, ts.t('緑(532 nm)が有利', 'green (532 nm) wins'),
        color=ts.W532, fontsize=13, ha='right', va='top', weight='bold')

ax.set_yscale('log')
ax.set_xlim(15, 150)
ax.set_ylim(0.22, 5.2)
ax.yaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4]))
ax.yaxis.set_minor_locator(NullLocator())
ax.yaxis.set_major_formatter(FixedFormatter(['0.25', '0.5', '1', '2', '4']))
ax.set_xlabel('Pressure  $P$  (GPa)')
ax.set_ylabel(ts.t('$\\eta$(532 nm) / $\\eta$(473 nm)\n(>1 なら青が有利)',
                   '$\\eta$(532 nm) / $\\eta$(473 nm)\n(>1: blue is better)'))
ax.grid(alpha=0.22)

plt.savefig('threshold_green_blue.png')
print(f'saved threshold_green_blue.png ; crossover {xo:.1f} GPa ; '
      + ', '.join(f'{p:.0f} GPa: {ratio(float(p)):.2f}'
                  for p in (20, 30, 50, 100, 120, 150)))
