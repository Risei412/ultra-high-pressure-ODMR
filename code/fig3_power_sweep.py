"""
fig3_power_sweep.py
Power-explicit blue-wavelength optimisation at a fixed pressure (Part B of PLAN.md):
  (a) eta(lambda, u) heatmap with the optimal-wavelength ridge overlaid
  (b) optimal blue wavelength vs power u (+ MC band), showing whether/where it
      jumps between the two competing local optima
  (c) eta(u) for real laser lines (405/457/473/475/488/532 nm) so the crossover
      power (which line is best) is visible directly

Run:
  python fig3_power_sweep.py            # 120 GPa
  python fig3_power_sweep.py 100        # any other pressure
Out:
  power_wavelength_sensitivity_<P>GPa.png
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from nv_model import mc_band
from nv_model_power import NVModelPower, default_randomiser_power

rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                 'axes.linewidth': 0.9, 'mathtext.fontset': 'dejavusans'})

PMAIN = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
Pi = int(round(PMAIN))

m = NVModelPower()
lam = np.linspace(402, 560, 900)
u_grid = np.logspace(-2, 2, 81)

# ---- (a) full eta(lambda,u) map ----
ETA = np.array([m.eta_lambda_u(lam, PMAIN, u)[0] for u in u_grid])   # shape (u, lam)
lam_opt = lam[np.nanargmin(ETA, axis=1)]
eta_opt = np.nanmin(ETA, axis=1)

# ---- (b) MC band on the optimal-wavelength ridge ----
def ridge(mm):
    e = np.array([mm.eta_lambda_u(lam, PMAIN, u)[0] for u in u_grid])
    return lam[np.nanargmin(e, axis=1)]
loL, loH = mc_band(default_randomiser_power, ridge, n=80, seed=11)

# ---- (c) candidate laser lines ----
lines = [(405, '#7a3fb0'), (457, '#2660c4'), (473, '#1a9850'),
         (475, '#c1272d'), (488, '#d98a00'), (532, '#555555')]
eta_lines = {lam0: np.array([m.eta_lambda_u(lam0, PMAIN, u)[0] for u in u_grid])
             for lam0, _ in lines}

fig, axs = plt.subplots(1, 3, figsize=(15.6, 4.6))

# (a) heatmap
ax = axs[0]
LOG_ETA = np.log10(np.clip(ETA, 1e-3, 1e4))
pc = ax.pcolormesh(lam, u_grid, LOG_ETA, shading='auto', cmap='viridis_r',
                    vmin=np.log10(1), vmax=np.log10(50))
ax.plot(lam_opt, u_grid, color='w', lw=2.2)
ax.plot(lam_opt, u_grid, color='#c1272d', lw=1.2, label='optimal $\\lambda$')
ax.set_yscale('log')
ax.set_xlabel('Blue excitation wavelength  $\\lambda$  (nm)')
ax.set_ylabel('Normalised intensity  $u = I/I_{ref}$')
ax.set_title(f'(a)  $\\log_{{10}}\\eta(\\lambda,u)$ @ {Pi} GPa', loc='left', weight='bold', fontsize=11)
cb = fig.colorbar(pc, ax=ax, pad=0.02); cb.set_label(r'$\log_{10}\eta$ (lower=better)')
ax.legend(frameon=False, fontsize=9, loc='upper right')

# (b) ridge vs power
ax = axs[1]
ax.fill_between(u_grid, loL, loH, color='#c1272d', alpha=0.18, lw=0)
ax.plot(u_grid, lam_opt, color='#c1272d', lw=2.4, label='optimal blue $\\lambda$')
for lam0, c in lines:
    ax.axhline(lam0, color=c, lw=1.0, ls='--', alpha=0.7)
    ax.text(u_grid[-1] * 1.05, lam0, f'{lam0}', color=c, fontsize=7.5, va='center')
ax.set_xscale('log')
ax.set_xlabel('Normalised intensity  $u = I/I_{ref}$')
ax.set_ylabel('Optimal blue wavelength  (nm)')
ax.set_title(f'(b)  Optimal $\\lambda$ vs power @ {Pi} GPa', loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=9, loc='upper left'); ax.grid(True, which='both', alpha=0.18)

# (c) eta(u) for candidate lines
ax = axs[2]
ax.plot(u_grid, eta_opt, color='k', lw=2.0, ls=':', label='global optimum')
for lam0, c in lines:
    ax.plot(u_grid, eta_lines[lam0], color=c, lw=1.8, label=f'{lam0} nm')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Normalised intensity  $u = I/I_{ref}$')
ax.set_ylabel(r'Lock-in sensitivity  $\eta$  (lower=better)')
ax.set_title(f'(c)  $\\eta(u)$ for candidate lines @ {Pi} GPa', loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, ncol=2); ax.grid(True, which='both', alpha=0.18)

plt.tight_layout()
out = f'power_wavelength_sensitivity_{Pi}GPa.png'
plt.savefig(out, dpi=185, bbox_inches='tight')

# ---- text summary ----
print(f'saved {out}')
print(f'{"u":>10} {"lambda_opt(nm)":>16} {"eta_opt":>10}')
for u, l, e in zip(u_grid[::10], lam_opt[::10], eta_opt[::10]):
    print(f'{u:10.3f} {l:16.1f} {e:10.3f}')
best_line_idx = np.argmin([eta_lines[l] for l, _ in lines], axis=0)
print('\nbest real laser line vs power (index into', [l for l, _ in lines], '):')
for u, idx in zip(u_grid[::10], best_line_idx[::10]):
    print(f'  u={u:8.3f}  best line = {lines[idx][0]} nm')
