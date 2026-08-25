"""
fig8_repro_contrast.py
----------------------
Reproduction of published cwODMR contrast against pressure — the quantitative
half of `repro_literature.py`, drawn.

Two independent experiments, one model curve:

  [Dai22]  J.-H. Dai et al., Chin. Phys. Lett. 39, 117601 (2022).  Microdiamonds,
           532 nm, room temperature: 14% at ambient falling to a ~1-3% plateau.
           These three points CALIBRATE the ISC prefactor C0(P) (two constants,
           C_amb and E_isc) and are marked as such — they are not a test.
  [Hil23]  A. Hilberer et al., PRB 107, L220102 (2023).  NV implanted in the
           ANVIL of a FIB-machined micropillar, excited at 457 nm and (at the
           highest pressure) 405 nm: 5% at 73 GPa, 3% at 103, 1.5% at 131.
           A different sample geometry, a different excitation, no free
           parameters left — this is the test, and the model lands within a
           factor 0.6-1.2 of it.

Contrast is diluted by the NV0 background, which is wavelength dependent, so the
model is drawn at both excitation wavelengths rather than at one.

Run:  python fig8_repro_contrast.py [--style st]
Out:  repro_contrast[_st].png
"""

import numpy as np
import matplotlib.pyplot as plt

import talk_style as ts
from nv_model import NVModel

ts.use()

m = NVModel(T=300.0)
P = np.linspace(0., 150., 301)
C532 = np.array([float(m.eta_lambda(532., p)[3]) for p in P]) * 100.
P457 = np.linspace(40., 150., 221)      # 457 nm is the high-pressure line;
C457 = np.array([float(m.eta_lambda(457., p)[3])                # below ~12 GPa
                 for p in P457]) * 100.  # its one-photon ionisation edge kinks

# observations, as published
DAI = [(0.0, 14.0), (102.3, 3.0), (137.7, 1.2)]                 # 532 nm  [CAL]
HIL = [(73.0, 5.0, 457.), (103.0, 3.0, 457.), (131.0, 1.5, 405.)]

ratios = [float(m.eta_lambda(lam, p)[3]) * 100. / obs for p, obs, lam in HIL]

fig, ax = plt.subplots(figsize=(8.4, 4.6) if ts.LEAN else (10.8, 5.9))

ax.plot(P, C532, color=ts.ACCENT, lw=2.8, label='model, 532 nm')
ax.plot(P457, C457, color=ts.ACCENT, lw=1.8, ls='--', label='model, 457 nm')

ax.plot([p for p, _ in DAI], [c for _, c in DAI], 's', ms=15, mfc='none',
        mec=ts.MUTED, mew=2.2, ls='none', zorder=7,
        label=ts.t('Dai 2022 — 較正に使用', 'Dai 2022 — used for calibration'))
ax.plot([p for p, _, _ in HIL], [c for _, c, _ in HIL], 'o', ms=10,
        color=ts.CRIT, mec='white', mew=1.4, ls='none', zorder=6,
        label=ts.t('Hilberer 2023 micropillar — 独立',
                   'Hilberer 2023 micropillar — independent'))

ax.set_yscale('log')
ax.set_xlim(-3, 152)
ax.set_ylim(0.8, 25)
ax.set_yticks([1, 2, 5, 10, 20])
ax.set_yticklabels(['1', '2', '5', '10', '20'])
ax.set_xlabel('Pressure  $P$  (GPa)')
ax.set_ylabel(ts.t('cwODMR コントラスト  $C$  (%)',
                   'cwODMR contrast  $C$  (%)'))
ax.grid(alpha=0.22)
ax.legend(fontsize=12 if ts.LEAN else 11.5, loc='lower left')

# The interpretation — agreement within x0.6-1.2 with the data the model was
# NOT calibrated on — belongs on the slide, not inside the frame.

plt.savefig(f'repro_contrast{ts.SUF}.png')
print(f'saved repro_contrast{ts.SUF}.png ; '
      + ', '.join(f'{p:.0f} GPa obs {obs:.1f}% model '
                  f'{float(m.eta_lambda(lam, p)[3]) * 100:.1f}% '
                  f'(x{r:.2f})' for (p, obs, lam), r in zip(HIL, ratios)))
