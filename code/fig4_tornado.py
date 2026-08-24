"""
fig4_tornado.py
---------------
Anchor sensitivity of lambda_opt at 120 GPa (PLAN.md Part A3 / Part D).
Also TALK FIGURE 3 -- why the answer is not a fit.

Each input is displaced to the edges of its Monte-Carlo range with all others
held nominal, and the resulting shift of lambda_opt is recorded.  Two measured
optical quantities move it.  Every phenomenological charge-transfer constant --
the ones a calibration experiment would fit -- moves it by exactly nothing,
because excited-state ionisation and recombination both scale with sigma_abs and
cancel from f_minus.  That cancellation is what makes lambda_opt = argmax
sigma_abs a structural result rather than a fitted number; tests/test_freeze.py
locks it.

Run:  python fig4_tornado.py
Out:  tornado_lambda_opt_120GPa.png  (+ the table, on stdout)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import talk_style as ts
from nv_model import NVModel

ts.use()

P = 120.0
BASE = dict(T=300.0)
base = NVModel(**BASE).lambda_opt(P)

# (label, kwarg, low, high, range caption, is_measured)
ROWS = [
    ('$\\Delta E_{\\mathrm{ZPL}}$(120 GPa)', 'dE120', 0.380, 0.420,
     '±20 meV', True),
    ('$S_{\\mathrm{abs}}$ slope', 'S_slope', 0.85 * (4.61 - 3.08),
     1.15 * (4.61 - 3.08), '±15%', True),
    ('d$E_{\\mathrm{ZPL}}$/d$P|_0$', 'slope0', 5.50e-3, 6.00e-3,
     '±0.25 meV/GPa', False),
    ('$\\sigma_{\\mathrm{ZPL}}$', 'zpl_width', 0.7 * 0.015, 1.4 * 0.015,
     '0.7–1.4×', False),
    ('$a_{\\mathrm{gs}}$', 'a_gs', 3.0, 9.0, '±50%', False),
    ('$a_{\\mathrm{es}}$', 'a_es', 0.45, 1.35, '±50%', False),
    ('$r_0$', 'r0', 0.84, 1.56, '±30%', False),
    ('$r_{\\mathrm{bg}}$', 'rbg', 0.09, 0.21, '±40%', False),
    ('$w_0$', 'w0', 0.7, 1.3, '±30%', False),
    ('$C_{\\mathrm{amb}}$', 'C_amb', 0.5 * 0.2469, 1.5 * 0.2469, '±50%', False),
    ('$E_{\\mathrm{ISC}}$', 'E_isc', 0.7 * 0.1807, 1.3 * 0.1807, '±30%', False),
]

rows = []
for label, key, lo, hi, rng, measured in ROWS:
    d_lo = NVModel(**BASE, **{key: lo}).lambda_opt(P) - base
    d_hi = NVModel(**BASE, **{key: hi}).lambda_opt(P) - base
    rows.append((label, rng, d_lo, d_hi, measured))
    print(f'{label:34s} {rng:16s} {d_lo:+7.2f}  {d_hi:+7.2f} nm')

# measured movers first, then the flat ones in the order given
rows.sort(key=lambda r: (-max(abs(r[2]), abs(r[3]))))

fig, ax = plt.subplots(figsize=(8.4, 5.2) if ts.LEAN else (10.6, 6.2))
y = np.arange(len(rows))[::-1]

for yi, (label, rng, d_lo, d_hi, measured) in zip(y, rows):
    col = ts.ACCENT if measured else ts.MUTED
    span = max(abs(d_lo), abs(d_hi))
    if span < 0.01:
        ax.plot(0, yi, 'o', color=col, ms=7, zorder=3)
        if not ts.LEAN:      # one repeated label per row is redundant on a slide
            ax.text(0.35, yi, ts.t('0.00 nm — 動かない', '0.00 nm — no effect'),
                    va='center', ha='left', fontsize=11.5, color=ts.MUTED)
    else:
        ax.barh(yi, d_lo, height=0.55, color=col, alpha=0.95, zorder=3)
        ax.barh(yi, d_hi, height=0.55, color=col, alpha=0.55, zorder=3)
        ax.text(max(d_lo, d_hi) + 0.25, yi, f'±{span:.1f} nm',
                va='center', ha='left', fontsize=14 if ts.LEAN else 12, color=ts.INK,
                weight='bold')

ax.axvline(0, color=ts.INK2, lw=1.4, zorder=4)
ax.set_yticks(y)
ax.set_yticklabels([f'{r[0]}   {r[1]}' for r in rows],
                   fontsize=14 if ts.LEAN else 12.5)
ax.set_xlim(-5.2, 6.8)
ax.set_xlabel(ts.t('120 GPa における $\\lambda_{\\mathrm{opt}}$ の変化  (nm)',
                   'shift of $\\lambda_{opt}$ at 120 GPa  (nm)'))
ax.grid(axis='x', alpha=0.25)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)

n_meas = sum(1 for r in rows if r[4])
ax.axhline(y[n_meas] + 0.5, color=ts.RULE, lw=1.2)
ax.set_ylim(y.min() - 0.7, y.max() + 0.7)

handles = [Patch(facecolor=ts.ACCENT, label=ts.t('測定量 — 答えを動かす',
                                                 'measured — moves it')),
           Patch(facecolor=ts.MUTED, label=ts.t('現象論定数 — 動かさない',
                                                'fitted — moves nothing'))]
if ts.LEAN:
    zero_rows = [yi for yi, r in zip(y, rows) if max(abs(r[2]), abs(r[3])) < 0.01]
    ax.text(0.5, float(np.mean(zero_rows)), ts.t('すべて 0.00 nm', 'all 0.00 nm'),
            va='center', ha='left', fontsize=13, color=ts.MUTED)

leg = ax.legend(handles=handles, loc='center left', bbox_to_anchor=(0.0, 0.40),
                fontsize=12.5, handlelength=1.1, handleheight=1.1,
                borderpad=0.9, labelspacing=0.8, frameon=True)
leg.get_frame().set_facecolor('white')
leg.get_frame().set_edgecolor(ts.RULE)
leg.get_frame().set_alpha(1.0)
leg.set_zorder(10)

plt.savefig(f'tornado_lambda_opt_120GPa{ts.SUF}.png')
print(f'\nsaved tornado_lambda_opt_120GPa{ts.SUF}.png ; baseline lambda_opt = {base:.2f} nm')
