"""
fig2_blue_wavelength_sweep.py
Blue-excitation-wavelength sweep of the lock-in sensitivity at a fixed pressure:
  (a) eta vs wavelength (main pressure + two comparison pressures)
  (b) mechanism (sigma_abs and f-) at the main pressure
  (c) optimal blue wavelength vs pressure (tracks ZPL / sideband edge)

Run:
  python fig2_blue_wavelength_sweep.py            # 120 GPa (compare 100, 140)
  python fig2_blue_wavelength_sweep.py 100 75 125 # 100 GPa (compare 75, 125)
Out:
  blue_wavelength_sensitivity_<P>GPa.png
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from nv_model import NVModel, HBARC, HW, mc_band, default_randomiser

rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                 'axes.linewidth': 0.9, 'mathtext.fontset': 'dejavusans'})

# ---- pressures (main + two comparison) from CLI, defaults = 120,100,140 ----
args = [float(a) for a in sys.argv[1:]]
PMAIN = args[0] if len(args) >= 1 else 120.0
PC1   = args[1] if len(args) >= 2 else 100.0
PC2   = args[2] if len(args) >= 3 else 140.0

m = NVModel()
lam = np.linspace(402, 560, 900)

def curve(P):                                   # eta(lambda) normalised to its own optimum
    e = m.eta_lambda(lam, P)[0]
    return e / np.nanmin(e)

cM, cC1, cC2 = curve(PMAIN), curve(PC1), curve(PC2)
iM = np.nanargmin(cM); lam_opt = lam[iM]

# MC band for the main-pressure curve
bL, bH = mc_band(default_randomiser,
                 lambda mm: (lambda e: e / np.nanmin(e))(mm.eta_lambda(lam, PMAIN)[0]),
                 n=250, seed=5)

# mechanism at main pressure
_, fM, sM, _, _ = m.eta_lambda(lam, PMAIN)

# optimal lambda vs pressure (+ band)
Pg = np.linspace(40, 140, 120)
def optlam(mm):
    return np.array([lam[np.nanargmin(mm.eta_lambda(lam, p)[0])] for p in Pg])
lopt = optlam(m)
loL, loH = mc_band(default_randomiser, optlam, n=150, seed=7)
zpl_nm  = HBARC / m.ZPL(Pg)
edge_nm = HBARC / (m.ZPL(Pg) + 0.5 * m.Sabs(Pg) * HW)

Pi = int(round(PMAIN))
cO, c457, c487, cZ, cIP = '#c1272d', '#2660c4', '#7a3fb0', '0.45', '#d98a00'
fig, axs = plt.subplots(1, 3, figsize=(15.2, 4.4))

# (a) eta vs wavelength
ax = axs[0]
ax.fill_between(lam, bL, bH, color=cO, alpha=0.16, lw=0)
ax.plot(lam, cM,  color=cO,    lw=2.6, label=f'{Pi} GPa')
ax.plot(lam, cC1, color='0.55', lw=1.4, ls='--', label=f'{int(PC1)} GPa')
ax.plot(lam, cC2, color='0.2',  lw=1.4, ls=':',  label=f'{int(PC2)} GPa')
ax.axvline(457, color=c457, lw=1.4, ls='--')
ax.axvline(487, color=c487, lw=1.2, ls='-.')
ax.plot(lam_opt, cM[iM], 'o', color=cO, ms=7)
ax.annotate(f'optimum\n{lam_opt:.0f} nm', (lam_opt, 1.0), (lam_opt + 4, 1.7),
            color=cO, fontsize=9, ha='left')
ax.text(455, 4.3, '457 nm', color=c457, fontsize=8.5, ha='right')
ax.text(489, 5.0, '487 nm', color=c487, fontsize=8, ha='left')
ax.set_yscale('log'); ax.set_ylim(0.9, 7); ax.set_xlim(400, 560)
ax.set_xlabel('Blue excitation wavelength  $\\lambda$  (nm)')
ax.set_ylabel(r'Lock-in sensitivity  $\eta$  (norm. to each optimum, lower=better)')
ax.set_title('(a)  $\\eta$ vs blue wavelength', loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=9, loc='upper right'); ax.grid(True, which='both', alpha=0.18)

# (b) mechanism
ax = axs[1]; axr = ax.twinx()
ax.plot(lam, sM, color=cO, lw=2.4, label=r'$\sigma_{abs}$')
axr.plot(lam, fM, color=cIP, lw=2.0, ls=':', label=r'$f_-$')
ax.axvline(487, color=c487, lw=1.0, ls='-.')
ax.axvline(HBARC / m.IP_A2(PMAIN), color=cIP, lw=1.2, ls='-.')
ax.axvline(HBARC / m.ZPL(PMAIN),  color=cZ,  lw=1.2, ls='-')
ax.plot(lam_opt, np.interp(lam_opt, lam, sM), 'o', color=cO, ms=6)
ax.text(HBARC / m.IP_A2(PMAIN) + 2, 3.1,
        f'IP($^3A_2$) {HBARC/m.IP_A2(PMAIN):.0f} nm\n(ionisation edge)', color=cIP, fontsize=8, ha='left')
ax.text(HBARC / m.ZPL(PMAIN) + 2, 2.0, f'ZPL\n{HBARC/m.ZPL(PMAIN):.0f} nm', color=cZ, fontsize=8, ha='left')
ax.text(lam_opt - 3, np.interp(lam_opt, lam, sM) + 0.1, f'opt {lam_opt:.0f} nm',
        color=cO, fontsize=8, ha='right')
ax.set_xlim(400, 560); ax.set_ylim(0, 3.9); axr.set_ylim(0, 1.0)
ax.set_xlabel('Blue excitation wavelength  $\\lambda$  (nm)')
ax.set_ylabel(r'Absorption cross section  $\sigma_{abs}$ (norm.)')
axr.set_ylabel(r'NV$^-$ fraction  $f_-$')
ax.set_title(f'(b)  Mechanism at {Pi} GPa', loc='left', weight='bold', fontsize=11)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axr.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=9, loc='upper left'); ax.grid(True, alpha=0.18)

# (c) optimal lambda vs pressure
ax = axs[2]
ax.fill_between(Pg, loL, loH, color=cO, alpha=0.18, lw=0)
ax.plot(Pg, lopt,    color=cO, lw=2.6, label='optimal blue $\\lambda$')
ax.plot(Pg, zpl_nm,  color=cZ, lw=1.6, ls='-',  label='ZPL')
ax.plot(Pg, edge_nm, color='#1a9850', lw=1.6, ls='--', label='sideband edge')
ax.axhline(457, color=c457, lw=1.4, ls='--'); ax.text(42, 459, '457 nm fixed', color=c457, fontsize=8.5)
ax.axhline(487, color=c487, lw=1.2, ls='-.'); ax.text(42, 489, '487 nm', color=c487, fontsize=8.5)
lP = np.interp(PMAIN, Pg, lopt)
ax.axvline(PMAIN, color='0.6', lw=1.0, ls=':')
ax.plot(PMAIN, lP, 'o', color=cO, ms=8, zorder=5)
ax.annotate(f'{lP:.0f} nm @{Pi} GPa', (PMAIN, lP), (PMAIN - 24, lP - 20),
            color=cO, fontsize=9, arrowprops=dict(arrowstyle='->', color=cO, lw=1))
ax.set_xlim(40, 140); ax.set_ylim(450, 590)
ax.set_xlabel('Pressure  $P$  (GPa)'); ax.set_ylabel('Wavelength  (nm)')
ax.set_title('(c)  Optimal blue $\\lambda$ shifts with $P$', loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=9, loc='upper right'); ax.grid(True, alpha=0.18)

plt.tight_layout()
out = f'blue_wavelength_sensitivity_{Pi}GPa.png'
plt.savefig(out, dpi=185, bbox_inches='tight')
print(f'saved {out} ; optimum@{Pi} = {lam_opt:.0f} nm ; '
      f'eta(487)/opt={np.interp(487,lam,cM):.2f} ; eta(457)/opt={np.interp(457,lam,cM):.2f}')
