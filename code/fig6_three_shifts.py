"""
fig6_three_shifts.py
--------------------
TALK FIGURE 1 -- why the question is not trivial.

At 120 GPa two edges close in on the excitation window from opposite sides:

  * the Franck-Condon ABSORPTION envelope blue shifts and sharpens, its maximum
    moving 586 -> 474 nm (ZPL +0.40 eV, S_abs 3.08 -> 4.61);
  * the ground-state PHOTOIONISATION edge blue shifts too, 463 -> 405 nm
    (IP(3A2) 2.68 -> 3.06 eV), and everything bluer than it converts NV- to NV0.

The gap they leave is the answer.  The two conventional choices sit outside it
for opposite reasons: 532 nm is left behind by the envelope, 405 nm lands on the
ionisation edge.

The ambient curve is drawn for contrast; note that the single effective phonon
mode of the model is a poor approximation below ~50 GPa (its ambient maximum
falls at 586 nm, redder than the ~560 nm of the real NV- absorption), so the
0 GPa trace is an illustration of the SHIFT, not a prediction at ambient.

Run:  python fig6_three_shifts.py
Out:  talk_three_shifts.png
"""

import numpy as np
import matplotlib.pyplot as plt

import talk_style as ts
from nv_model import NVModel, nm2eV, HBARC

ts.use()

m = NVModel(T=300.0)
lam = np.linspace(380., 700., 1400)
s0 = m.sigma_abs(nm2eV(lam), 0.)
s120 = m.sigma_abs(nm2eV(lam), 120.)
pk0, pk120 = lam[s0.argmax()], lam[s120.argmax()]
ed0, ed120 = HBARC / m.IP_A2(0.), HBARC / m.IP_A2(120.)

fig, ax = plt.subplots(figsize=(11.0, 5.6))

# --- forbidden region: one-photon ionisation out of 3A2 at 120 GPa ---------
ax.axvspan(lam[0], ed120, color=ts.W405, alpha=0.10, lw=0)
ax.axvline(ed120, color=ts.W405, lw=2.0)
ax.axvline(ed0, color=ts.W405, lw=1.4, ls=':')

# --- the two envelopes -----------------------------------------------------
ax.fill_between(lam, s0, color=ts.PAST, alpha=0.16, lw=0)
ax.plot(lam, s0, color=ts.PAST, lw=2.2, label='0 GPa')
ax.fill_between(lam, s120, color=ts.ACCENT, alpha=0.16, lw=0)
ax.plot(lam, s120, color=ts.ACCENT, lw=2.8, label='120 GPa')

ax.set_xlim(380, 700)
ax.set_ylim(0, 3.5)
ax.set_xlabel('Excitation wavelength  $\\lambda$  (nm)')
ax.set_ylabel('Absorption cross section  $\\sigma_{\\mathrm{abs}}$\n'
              '(normalised to 532 nm at ambient)')
ax.grid(axis='y', alpha=0.25)

# --- the two shifts, as arrows --------------------------------------------
ax.annotate('', xy=(pk120, 2.92), xytext=(pk0, 2.92),
            arrowprops=dict(arrowstyle='-|>', color=ts.ACCENT, lw=2.2,
                            shrinkA=0, shrinkB=0))
ax.text(0.5 * (pk0 + pk120), 2.99,
        ts.t(f'吸収極大  {pk0:.0f} → {pk120:.0f} nm',
           f'absorption maximum  {pk0:.0f} -> {pk120:.0f} nm'),
        color=ts.ACCENT, ha='center', va='bottom', fontsize=13.5, weight='bold')

ax.annotate('', xy=(ed120, 0.42), xytext=(ed0, 0.42),
            arrowprops=dict(arrowstyle='-|>', color=ts.W405, lw=2.2,
                            shrinkA=0, shrinkB=0))
ax.text(0.5 * (ed0 + ed120), 0.50,
        ts.t(f'イオン化端  {ed0:.0f} → {ed120:.0f} nm',
           f'ionisation edge  {ed0:.0f} -> {ed120:.0f} nm'),
        color=ts.W405, ha='center', va='bottom', fontsize=13.5, weight='bold')
ax.text(0.5 * (lam[0] + ed120) - 1.5, 1.9,
        ts.t('NV$^-$ → NV$^0$\n1 光子イオン化', 'NV$^-$ -> NV$^0$\n1-photon'),
        color=ts.W405, fontsize=11, va='center', ha='center', rotation=90)
ax.text(ed0 + 3, 0.07, ts.t('463 nm (0 GPa)', '463 nm (0 GPa)'),
        color=ts.W405, fontsize=10.5, ha='left')

# --- the three laser lines -------------------------------------------------
for lam0, col, lab, ha in ((532., ts.W532, '532 nm', 'left'),
                           (473., ts.W473, '473 nm', 'right')):
    ax.axvline(lam0, color=col, lw=1.8, ls='--', alpha=0.9)
    ax.text(lam0 + (4 if ha == 'left' else -4), 3.44, lab, color=col,
            fontsize=13.5, weight='bold', ha=ha, va='top')

# the ratio that carries the message, taken AT ONE PRESSURE so that it does not
# inherit the ambient-envelope approximation
b473 = float(m.sigma_abs(nm2eV(473.), 120.))
b532 = float(m.sigma_abs(nm2eV(532.), 120.))
for lam0, col, val in ((473., ts.W473, b473), (532., ts.W532, b532)):
    ax.plot(lam0, val, 'o', color=col, ms=11, mec='white', mew=1.6, zorder=5)
ax.annotate('', xy=(532, b532), xytext=(532, b473),
            arrowprops=dict(arrowstyle='<->', color=ts.INK2, lw=1.6))
ax.plot([473, 532], [b473, b473], color=ts.INK2, lw=1.0, ls=':')
ax.text(537, 0.5 * (b473 + b532),
        ts.t(f'120 GPa で\n$\\sigma$(473) は\n$\\sigma$(532) の {b473 / b532:.0f} 倍',
             f'at 120 GPa\n$\\sigma$(473) is {b473 / b532:.0f}x\n$\\sigma$(532)'),
        color=ts.INK2, fontsize=12.5, va='center', ha='left')

ax.legend(loc='upper right', fontsize=13.5)

plt.savefig('talk_three_shifts.png')
print(f'saved talk_three_shifts.png ; peak {pk0:.0f} -> {pk120:.0f} nm ; '
      f'edge {ed0:.0f} -> {ed120:.0f} nm ; '
      f'sigma(473)/sigma(532) at 120 GPa = {b473 / b532:.1f}')
