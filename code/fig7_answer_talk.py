"""
fig7_answer_talk.py
-------------------
TALK FIGURE 2 -- the answer.

Single-panel version of Fig. fig:sweep(a) for projection: one pressure, a linear
ordinate, the 5% tolerance window drawn as a band, and the three laser lines
labelled with what they cost.  The pressure tracking of panel (c) survives as an
inset, because "choose it for your target pressure" is one sentence of the talk,
not a panel of its own.

Run:  python fig7_answer_talk.py [P]      (default 120 GPa)
Out:  talk_answer_<P>GPa.png
"""

import sys

import numpy as np
import matplotlib.pyplot as plt

import talk_style as ts
from nv_model import NVModel, mc_band, default_randomiser

ts.use()

P = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
Pi = int(round(P))

m = NVModel(T=300.0)
lam = np.linspace(400., 560., 1200)
eta = np.asarray(m.eta_lambda(lam, P)[0])
eta = eta / np.nanmin(eta)
lam_opt = m.lambda_opt(P)

# 5% tolerance window
inside = lam[eta <= 1.05]
w_lo, w_hi = inside.min(), inside.max()

bL, bH = mc_band(lambda rng: default_randomiser(rng, T=300.0),
                 lambda mm: (lambda e: e / np.nanmin(e))(
                     np.asarray(mm.eta_lambda(lam, P)[0])),
                 n=250, seed=5)

fig, ax = plt.subplots(figsize=(11.0, 5.8))

ax.axvspan(w_lo, w_hi, color=ts.BAND, alpha=0.85, lw=0, zorder=0)
ax.axhline(1.05, color=ts.ACCENT, lw=1.0, ls=':', alpha=0.8, zorder=1)

ax.fill_between(lam, bL, bH, color=ts.ACCENT, alpha=0.16, lw=0, zorder=2)
ax.plot(lam, eta, color=ts.ACCENT, lw=3.0, zorder=3)

# --- the three lines, placed by hand so nothing collides -------------------
PLACE = {405.: (ts.W405, ts.t('イオン化端', 'ionisation edge'), 407., 3.95, 'left'),
         473.: (ts.W473, ts.t('市販 DPSS', 'commercial DPSS'), 492., 2.62, 'left'),
         532.: (ts.W532, ts.t('従来の既定', 'inherited default'), 527., 3.95, 'right')}
for lam0, (col, note, tx, ty, ha) in PLACE.items():
    pen = float(np.asarray(m.eta_lambda(lam0, P)[0])
                / np.asarray(m.eta_lambda(lam_opt, P)[0]))
    ax.axvline(lam0, color=col, lw=1.8, ls='--', alpha=0.9, zorder=4)
    if pen < 4.0:
        ax.plot(lam0, pen, 'o', color=col, ms=10, mec='white', mew=1.6, zorder=6)
    ax.text(tx, ty, f'{lam0:.0f} nm', color=col, fontsize=13.5, weight='bold',
            ha=ha, va='top')
    ax.text(tx, ty - 0.22, ts.t(f'×{pen:.1f}  {note}', f'x{pen:.1f}  {note}'),
            color=col, fontsize=12, ha=ha, va='top')

# window label: inside the band, clear of everything else
ax.text(0.5 * (w_lo + w_hi), 3.95,
        ts.t(f'5% 許容窓\n{w_lo:.0f}–{w_hi:.0f} nm',
             f'5% window\n{w_lo:.0f}-{w_hi:.0f} nm'),
        color=ts.ACCENT, ha='center', va='top', fontsize=13, weight='bold')

# the optimum
ax.plot(lam_opt, 1.0, 'o', color=ts.ACCENT, ms=13, mec='white', mew=2, zorder=6)
ax.annotate(ts.t(f'最適  {lam_opt:.0f} nm', f'optimum  {lam_opt:.0f} nm'),
            xy=(lam_opt, 1.0), xytext=(492., 1.45),
            color=ts.ACCENT, fontsize=14, weight='bold',
            arrowprops=dict(arrowstyle='->', color=ts.ACCENT, lw=1.6))

# pressure tracking, as one line of text rather than an inset panel
ax.text(403., 1.28,
        ts.t(f'$\\lambda_{{\\mathrm{{opt}}}}$ は圧力で動く: '
             f'{m.lambda_opt(100.):.0f} nm @100 GPa → {lam_opt:.0f} nm @120 GPa',
             f'$\\lambda_{{opt}}$ moves with pressure: '
             f'{m.lambda_opt(100.):.0f} nm @100 GPa -> {lam_opt:.0f} nm @120 GPa'),
        color=ts.INK2, fontsize=12.5, ha='left', va='center')

ax.set_xlim(400, 560)
ax.set_ylim(0.85, 4.05)
ax.set_xlabel('Excitation wavelength  $\\lambda$  (nm)')
ax.set_ylabel(ts.t('ロックイン感度  $\\eta/\\eta_{\\mathrm{opt}}$\n(小さいほど良い)',
                   'Lock-in sensitivity  $\\eta/\\eta_{\\mathrm{opt}}$\n(lower is better)'))
ax.grid(axis='y', alpha=0.25)
ax.set_title(ts.t(f'{Pi} GPa の答え:474 nm、市販 473 nm でよい',
                  f'The answer at {Pi} GPa: 474 nm, and 473 nm buys it'),
             loc='left', weight='bold', pad=14)

out = f'talk_answer_{Pi}GPa.png'
plt.savefig(out)
print(f'saved {out} ; lambda_opt = {lam_opt:.1f} nm ; '
      f'5% window {w_lo:.1f}-{w_hi:.1f} nm')
