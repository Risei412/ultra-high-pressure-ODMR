"""
fig1_green_blue_mix.py
Sensitivity of green (532 nm) / blue (457 nm) / mixture vs pressure (0-140 GPa),
plus the underlying absorption cross section and NV- fraction.

Run:  python fig1_green_blue_mix.py
Out:  sensitivity_green_blue_mix.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from nv_model import NVModel, mc_band, default_randomiser

rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                 'axes.linewidth': 0.9, 'mathtext.fontset': 'dejavusans'})

P = np.linspace(0, 140, 561)

# collection=False on purpose.  eta_col (C-7) is WAVELENGTH INDEPENDENT, so it
# is a factor common to all three excitation schemes and cancels exactly from
# every ratio this figure shows -- including the green/blue crossover pressure.
# Carrying it here would only fold the pressure dependence of the DETECTION
# band into a figure about the EXCITATION choice; it is treated on its own in
# Sec. IV E.  (It is also ~10^3 times more expensive than the rest of the
# model, which made the Monte Carlo bands unusable.)
m = NVModel(collection=False)

green = [(532, 1.0)]
blue  = [(457, 1.0)]
mix   = [(532, 0.5), (457, 0.5)]          # equal TOTAL optical power

etaG, fG, sG, CG, RG = m.eta(green, P)
etaB, fB, sB, CB, RB = m.eta(blue,  P)
etaM, fM, sM, CM, RM = m.eta(mix,   P)
ref = np.min(etaG)                        # normalise to best green value
etaG, etaB, etaM = etaG / ref, etaB / ref, etaM / ref

# Monte-Carlo 16-84% bands
def band(beams):
    def rnd(rng):
        mm = default_randomiser(rng)
        mm.collection = False
        return mm
    return mc_band(rnd,
                   lambda mm: mm.eta(beams, P)[0] / np.min(mm.eta(green, P)[0]),
                   n=300, seed=1)
gL, gH = band(green); bL, bH = band(blue); mL, mH = band(mix)

# crossover pressure (blue overtakes green)
xover = P[np.argmax(etaB < etaG)]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.5))
cg, cb, cm = '#2e8b3f', '#2660c4', '#c1272d'

ax.axvspan(0, xover, color=cg, alpha=0.05); ax.axvspan(xover, 140, color=cb, alpha=0.05)
ax.axvline(70, color='0.45', ls=':', lw=1.2)
for lo, hi, c in [(gL, gH, cg), (bL, bH, cb), (mL, mH, cm)]:
    ax.fill_between(P, lo, hi, color=c, alpha=0.18, lw=0)
ax.plot(P, etaG, color=cg, lw=2.4, label='Green 532 nm')
ax.plot(P, etaB, color=cb, lw=2.4, label='Blue 457 nm')
ax.plot(P, etaM, color=cm, lw=2.4, ls='--', label='Mix 532+457 nm')
ax.set_yscale('log'); ax.set_xlim(0, 140); ax.set_ylim(0.8, 30)
ax.set_xlabel('Pressure  $P$  (GPa)')
ax.set_ylabel(r'Lock-in ODMR sensitivity  $\eta \propto \Delta\nu/(C\sqrt{R})$'
              + '\n(normalised, lower = better)')
ax.text(43, 0.93, 'green\noptimal', ha='center', va='bottom', color=cg, fontsize=9)
ax.text(112, 0.93, 'blue\noptimal', ha='center', va='bottom', color=cb, fontsize=9)
ax.text(71, 18, 'ref. Ho et al.\n70 GPa', fontsize=8, color='0.4', ha='left')
ax.axhline(1, color='0.7', lw=0.7)
ax.legend(frameon=False, loc='upper center', fontsize=9.5)
ax.set_title('(a)  Sensitivity vs pressure', loc='left', fontsize=11, weight='bold')
ax.grid(True, which='both', alpha=0.18)

axr = ax2.twinx()
ax2.plot(P, sG, color=cg, lw=2.2, label=r'$\sigma_{abs}$ 532')
ax2.plot(P, sB, color=cb, lw=2.2, label=r'$\sigma_{abs}$ 457')
axr.plot(P, fG, color=cg, lw=1.6, ls=':')
axr.plot(P, fB, color=cb, lw=1.6, ls=':')
ax2.set_xlim(0, 140); ax2.set_ylim(0, 4.0); axr.set_ylim(0, 1.02)
ax2.set_xlabel('Pressure  $P$  (GPa)')
ax2.set_ylabel(r'Absorption cross section  $\sigma_{abs}$  (norm., solid)')
axr.set_ylabel(r'NV$^-$ fraction  $f_-$  (dotted)')
ax2.axvline(10.4, color=cb, ls='--', lw=1, alpha=0.6)
ax2.text(11, 3.6, '457 nm 1-photon\nionisation off', color=cb, fontsize=8)
ax2.set_title('(b)  Underlying mechanism', loc='left', fontsize=11, weight='bold')
ax2.legend(frameon=False, loc='upper right', fontsize=9)
ax2.grid(True, alpha=0.18)

plt.tight_layout()
plt.savefig('sensitivity_green_blue_mix.png', dpi=190, bbox_inches='tight')
print('saved sensitivity_green_blue_mix.png ; crossover ~%.0f GPa' % xover)
