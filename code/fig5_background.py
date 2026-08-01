"""
fig5_background.py
Green-vs-blue re-assessment with spin-independent background light (Part E).

  (a) background shapes g_k(lambda_exc) and their sum, with the ruby
      transmission gap and the N3 penalty on 405 nm marked;
      diamond-Raman check printed/annotated (never enters the >650 nm window)
  (b) eta(lambda) at the main pressure for several rho0
  (c) optimal blue lambda vs rho0, with MC band and the 5% tolerance band
  (d) green/blue crossover pressure vs rho0, and eta(blue)/eta(green) at the
      main pressure -- i.e. "how bad may the background be before blue loses?"

Run:  python fig5_background.py            # 120 GPa
      python fig5_background.py 100        # other pressure
Out:  background_green_blue_<P>GPa.png
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

import background as bg
from nv_model import mc_band
from nv_bg import (NVModelBG, optimum_wavelength, tolerance_band,
                   crossover_pressure, randomiser_bg)

rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                 'axes.linewidth': 0.9, 'mathtext.fontset': 'dejavusans'})

PMAIN = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
Pi = int(round(PMAIN))

m = NVModelBG()
lam = np.linspace(402, 560, 1580)
rho_grid = np.logspace(-2, 3, 60)
RHOS_SHOWN = (0.0, 0.1, 1.0, 10.0, 100.0)
LINES = (405, 457, 473, 488, 532)

cO, cG, cB, cZ = '#c1272d', '#2e8b3f', '#2660c4', '0.45'
cRUBY, cN3, cBROAD = '#d1495b', '#8a6fbf', '#5b8c5a'

fig, axs = plt.subplots(2, 2, figsize=(12.6, 9.0))

# ------------------------------------------------- (a) background shapes ---
ax = axs[0, 0]
ax.plot(lam, bg.g_channel('ruby', lam),  color=cRUBY,  lw=2.0, label='ruby $R$ line (Cr$^{3+}$ U/Y bands)')
ax.plot(lam, bg.g_channel('n3', lam),    color=cN3,    lw=2.0, ls='--', label='N3 / A-band anvil')
ax.plot(lam, bg.g_channel('broad', lam), color=cBROAD, lw=2.0, ls=':',  label='broad deformation lum.')
ax.plot(lam, bg.g_total(lam),            color='k',    lw=2.6, label='total (equal mix)')
for x in LINES:
    ax.axvline(x, color='0.75', lw=0.8, ls='-')
    ax.text(x, 22, f'{x}', rotation=90, fontsize=7.5, color='0.4', ha='right', va='top')
i_gap = int(np.argmin(bg.g_channel('ruby', lam)))
ax.plot(lam[i_gap], bg.g_channel('ruby', lam)[i_gap], 'o', color=cRUBY, ms=6)
ax.annotate(f'Cr$^{{3+}}$ transmission gap\n{lam[i_gap]:.0f} nm',
            (lam[i_gap], bg.g_channel('ruby', lam)[i_gap]), (lam[i_gap] + 8, 0.09),
            color=cRUBY, fontsize=8.5,
            arrowprops=dict(arrowstyle='->', color=cRUBY, lw=1))
ax.annotate('N3 punishes 405 nm', (405, bg.g_channel('n3', 405)), (424, 3.4),
            color=cN3, fontsize=8.5, arrowprops=dict(arrowstyle='->', color=cN3, lw=1))
r1 = bg.raman_nm(np.array(LINES, float), 1)
ax.text(0.02, 0.03,
        f'diamond Raman (1332 cm$^{{-1}}$): {r1.min():.0f}-{r1.max():.0f} nm\n'
        f'all below the {bg.WINDOW_NM:.0f} nm window edge $\\Rightarrow$ excluded',
        transform=ax.transAxes, fontsize=8, color='0.35', va='bottom')
ax.set_yscale('log'); ax.set_xlim(402, 560); ax.set_ylim(0.03, 40)
ax.set_xlabel('Excitation wavelength  $\\lambda$  (nm)')
ax.set_ylabel('Background shape  $g_k(\\lambda)$   (norm. to 532 nm)')
ax.set_title('(a)  What the detection window collects besides NV',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='upper right')
ax.grid(True, which='both', alpha=0.18)

# ------------------------------------------------------ (b) eta(lambda) ----
ax = axs[0, 1]
shades = plt.cm.viridis(np.linspace(0.05, 0.85, len(RHOS_SHOWN)))
for rho0, c in zip(RHOS_SHOWN, shades):
    e = m.eta_lambda_bg(lam, PMAIN, rho0)[0]
    e = e / np.nanmin(e)
    lo = optimum_wavelength(m, PMAIN, rho0, lam)
    ax.plot(lam, e, color=c, lw=2.2,
            label=f'$\\rho_0={rho0:g}$  ($\\lambda_{{opt}}={lo:.0f}$ nm)')
    ax.plot(lo, 1.0, 'o', color=c, ms=6)
lo0 = optimum_wavelength(m, PMAIN, 0.0, lam)
ax.axvline(lo0, color='0.7', lw=0.9, ls=':')
ax.set_yscale('log'); ax.set_xlim(430, 545); ax.set_ylim(0.95, 4)
ax.set_xlabel('Blue excitation wavelength  $\\lambda$  (nm)')
ax.set_ylabel(r'$\eta$  (norm. to each curve''\'''s optimum, lower = better)')
ax.set_title(f'(b)  Optimum barely moves at {Pi} GPa',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='upper right')
ax.grid(True, which='both', alpha=0.18)

# ------------------------------------------- (c) lambda_opt vs rho0 + MC ---
ax = axs[1, 0]
lopt = np.array([optimum_wavelength(m, PMAIN, r, lam) for r in rho_grid])
tb = np.array([tolerance_band(m, PMAIN, r, lam) for r in rho_grid])
mcL, mcH = mc_band(randomiser_bg,
                   lambda mm: np.array([optimum_wavelength(mm, PMAIN, r, lam)
                                        for r in rho_grid]),
                   n=120, seed=11)
ax.fill_between(rho_grid, tb[:, 0], tb[:, 1], color='0.6', alpha=0.20, lw=0,
                label='$\\eta \\leq 1.05\\,\\eta_{opt}$ tolerance band')
ax.fill_between(rho_grid, mcL, mcH, color=cO, alpha=0.22, lw=0,
                label='MC 16-84% (random background mix)')
ax.plot(rho_grid, lopt, color=cO, lw=2.6, label='$\\lambda_{opt}$')
for x, c, ls in ((473, cB, '--'), (488, '#7a3fb0', '-.')):
    ax.axhline(x, color=c, lw=1.2, ls=ls)
    ax.text(rho_grid[1], x + 1, f'{x} nm (commercial)', color=c, fontsize=8)
ax.set_xscale('log'); ax.set_xlim(rho_grid[0], rho_grid[-1]); ax.set_ylim(455, 500)
ax.set_xlabel('Background-to-signal ratio at 532 nm, ambient   $\\rho_0$')
ax.set_ylabel('Optimal blue wavelength  (nm)')
ax.set_title(f'(c)  $\\lambda_{{opt}}$ vs background level at {Pi} GPa',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='lower right')
ax.grid(True, which='both', alpha=0.18)

# ------------------------------------------------ (d) crossover vs rho0 ----
ax = axs[1, 1]
xo457 = np.array([crossover_pressure(m, r, blue=457.0) for r in rho_grid])
xoOpt = np.array([crossover_pressure(m, r, blue=optimum_wavelength(m, PMAIN, r, lam))
                  for r in rho_grid])
ax.plot(rho_grid, xo457, color=cB, lw=2.6, label='blue = 457 nm (fixed line)')
ax.plot(rho_grid, xoOpt, color=cO, lw=2.6, ls='--', label='blue = $\\lambda_{opt}$')
ax.axhline(PMAIN, color='0.5', lw=1.0, ls=':')
ax.text(rho_grid[1], PMAIN + 1.5, f'target {Pi} GPa', color='0.4', fontsize=8.5)
ax.set_xscale('log'); ax.set_xlim(rho_grid[0], rho_grid[-1])
ax.set_ylim(60, max(PMAIN + 15, 110))
ax.set_xlabel('Background-to-signal ratio at 532 nm, ambient   $\\rho_0$')
ax.set_ylabel('Green $\\rightarrow$ blue crossover pressure  (GPa)')
ax.set_title('(d)  Crossover, and how much background blue can absorb',
             loc='left', weight='bold', fontsize=11)

axr = ax.twinx()
ratio = np.array([float(m.eta_bg([(457.0, 1.0)], PMAIN, r)[0]
                        / m.eta_bg([(532.0, 1.0)], PMAIN, r)[0]) for r in rho_grid])
ratioOpt = np.array([float(m.eta_bg([(optimum_wavelength(m, PMAIN, r, lam), 1.0)], PMAIN, r)[0]
                           / m.eta_bg([(532.0, 1.0)], PMAIN, r)[0]) for r in rho_grid])
axr.plot(rho_grid, ratio, color=cB, lw=1.4, ls=':')
axr.plot(rho_grid, ratioOpt, color=cO, lw=1.4, ls=':')
axr.axhline(1.0, color=cG, lw=1.2)
axr.text(rho_grid[0] * 1.4, 0.72, 'green wins above this line', color=cG,
         fontsize=8.5, ha='left', va='top')
axr.set_yscale('log'); axr.set_ylim(0.02, 2.0)
axr.set_ylabel(f'$\\eta_{{blue}}/\\eta_{{green}}$ at {Pi} GPa  (dotted)')
ax.legend(frameon=False, fontsize=8.5, loc='lower right')
ax.grid(True, which='both', alpha=0.18)

plt.tight_layout()
out = f'background_green_blue_{Pi}GPa.png'
plt.savefig(out, dpi=185, bbox_inches='tight')

# --------------------------------------------------------------- summary ---
print(f'saved {out}')
print(f'{"rho0":>8} {"lam_opt":>9} {"tol band":>18} {"xover(457)":>11} '
      f'{"xover(opt)":>11} {"eta_b/eta_g":>12}')
for r in RHOS_SHOWN:
    lo = optimum_wavelength(m, PMAIN, r, lam)
    t = tolerance_band(m, PMAIN, r, lam)
    rb = float(m.eta_bg([(lo, 1.0)], PMAIN, r)[0] / m.eta_bg([(532.0, 1.0)], PMAIN, r)[0])
    print(f'{r:8.2f} {lo:8.1f}nm  [{t[0]:6.1f},{t[1]:6.1f}] nm '
          f'{crossover_pressure(m, r, blue=457.0):10.1f} '
          f'{crossover_pressure(m, r, blue=lo):11.1f} {rb:12.3f}')
print(f'saturated (rho0 -> inf): lam_opt = {optimum_wavelength(m, PMAIN, 1e6, lam):.1f} nm, '
      f'eta_blue/eta_green = '
      f'{float(m.eta_bg([(optimum_wavelength(m, PMAIN, 1e6, lam), 1.0)], PMAIN, 1e6)[0] / m.eta_bg([(532.0, 1.0)], PMAIN, 1e6)[0]):.3f} '
      f'-> blue never loses at {Pi} GPa')
