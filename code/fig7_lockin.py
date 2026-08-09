"""
fig7_lockin.py
Is modulated detection (a lock-in) worth it for our CW-ODMR setup?

  (a) assumed noise PSD S(f) = S0 (1 + f_knee/f) and where each scheme measures
  (b) example simulated sweeps at a large f_knee -- DC baseline wander
  (c) sigma_D vs f_knee for DC / MW on-off / FM, with the break-even knees
  (d) sigma_D vs modulation frequency, and the FM depth optimum

Run:  python fig7_lockin.py
Out:  lockin_decision.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

import lockin_sim as ls

rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                 'axes.linewidth': 0.9, 'mathtext.fontset': 'dejavusans'})

N_MC = 200
SEED = 1
cDC, cOO, cFM = '#c1272d', '#2660c4', '#1a9850'
SCHEME = {'dc': ('DC sweep', cDC), 'onoff': ('MW on-off', cOO), 'fm': ('FM lock-in', cFM)}

fig, axs = plt.subplots(2, 2, figsize=(12.6, 8.8))

# ------------------------------------------------------------- (a) PSD -----
ax = axs[0, 0]
f = np.logspace(-1, 4, 500)
S0 = ls.shot_floor_psd(ls.R_PHOT)
for fk, c, lw in ((10.0, '0.65', 1.6), (100.0, '0.4', 1.8), (1000.0, 'k', 2.2)):
    ax.loglog(f, S0 * (1 + fk / f), color=c, lw=lw, label=f'$f_{{knee}}$ = {fk:.0f} Hz')
ax.axhline(S0, color=cFM, lw=1.4, ls='--')
ax.text(9e3, S0 * 1.3, 'shot-noise floor  $2/R$', color=cFM, fontsize=8.5, ha='right')
T_sweep = ls.T_TOTAL
ax.axvspan(1.0 / T_sweep / 3, 3.0 / T_sweep, color=cDC, alpha=0.13, lw=0)
ax.text(1.0 / T_sweep, S0 * 4e3, 'DC sweep\nmeasures here\n($\\sim 1/T_{sweep}$)',
        color=cDC, fontsize=8.5, ha='center')
ax.axvline(ls.F_MOD, color=cOO, lw=1.6, ls=':')
ax.text(ls.F_MOD * 0.85, S0 * 4e3, f'modulated\nschemes\n$f_{{mod}}$={ls.F_MOD/1e3:.0f} kHz',
        color=cOO, fontsize=8.5, ha='right')
ax.set_xlim(0.1, 1e4); ax.set_ylim(S0 * 0.4, S0 * 2e4)
ax.set_xlabel('Frequency  $f$  (Hz)')
ax.set_ylabel('Relative intensity noise PSD  (1/Hz)')
ax.set_title('(a)  The whole decision is set by $f_{knee}$',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='lower left'); ax.grid(True, which='both', alpha=0.18)

# --------------------------------------------------- (b) example sweeps ----
ax = axs[0, 1]
rng = np.random.default_rng(7)
fk_demo = 300.0
m = ls.dwell_samples(ls.FS, ls.T_TOTAL, ls.N_PTS, ls.F_MOD)[0]
noise = ls.technical_noise(m * ls.N_PTS, ls.FS, S0, fk_demo, rng)
clean = {}
for scheme in ('dc', 'onoff', 'fm'):
    nu, sig = ls.simulate_sweep(fk_demo, rng, scheme, noise=noise)
    _, ideal = ls.simulate_sweep(0.0, rng, scheme, tech_gain=0.0, shot=False)
    scale = np.max(np.abs(ideal - np.median(ideal)))
    lbl, c = SCHEME[scheme]
    ax.plot(nu - ls.D_TRUE, (sig - np.median(sig)) / scale,
            color=c, lw=1.7, marker='o', ms=3, label=lbl)
    clean[scheme] = (ideal - np.median(ideal)) / scale
for scheme, c in (('dc', cDC), ('onoff', cOO), ('fm', cFM)):
    ax.plot(nu - ls.D_TRUE, clean[scheme], color=c, lw=1.2, ls='--', alpha=0.45)
ax.axvline(0.0, color='0.6', lw=1.0, ls=':')
ax.text(0.6, -2.6, 'true line centre', color='0.4', fontsize=8.5)
ax.plot([], [], color='0.5', lw=1.2, ls='--', label='noise-free shape')
ax.set_ylim(-3.0, 3.0)
ax.set_xlabel('MW detuning from $D$  (MHz)')
ax.set_ylabel('Signal / noise-free peak amplitude')
ax.set_title(f'(b)  One sweep at $f_{{knee}}$ = {fk_demo:.0f} Hz',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='upper left'); ax.grid(True, alpha=0.18)

# --------------------------------------------------- (c) sigma_D vs knee ---
ax = axs[1, 0]
knees = np.concatenate([[0.0], np.logspace(0, 3.7, 16)])
res = {s: [] for s in SCHEME}
for fk in knees:
    sd = ls.sigma_D(fk, n_mc=N_MC, seed=SEED)
    for s in SCHEME:
        res[s].append(sd[s])
kx = np.where(knees == 0, 0.3, knees)          # plot f_knee = 0 at the left edge
for s, (lbl, c) in SCHEME.items():
    ax.loglog(kx, res[s], color=c, lw=2.4, marker='o', ms=3.5, label=lbl)

dc = np.array(res['dc'])
be = {}
for s in ('onoff', 'fm'):
    flat = float(np.median(res[s][:6]))
    i = int(np.argmax(dc > flat))
    be[s] = float(np.interp(flat, dc[i - 1:i + 1], kx[i - 1:i + 1]))
    ax.plot(be[s], flat, 'v', color=SCHEME[s][1], ms=9, zorder=5)
ax.annotate(f'break-even\n$f_{{knee}}$ $\\approx$ {be["fm"]:.0f} Hz',
            (be['fm'], np.median(res['fm'][:6])), (0.45, 0.55),
            fontsize=8.5, arrowprops=dict(arrowstyle='->', color='k', lw=1))
ax.set_xlim(0.25, 6e3)
ax.set_xlabel('Technical-noise knee  $f_{knee}$  (Hz)     [leftmost point: $f_{knee}=0$]')
ax.set_ylabel('$\\sigma_D$  (MHz)   lower = better')
ax.set_title('(c)  Modulation pays off only above the break-even knee',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='upper left')
ax.grid(True, which='both', alpha=0.18)

# ------------------------------------------------- (d) f_mod and FM depth --
ax = axs[1, 1]
ps = np.array([196, 98, 64, 48, 32, 24, 16, 12, 8, 6, 4])
fmods = ls.FS / ps
fk_d = 500.0
for s in ('onoff', 'fm'):
    y = [ls.sigma_D(fk_d, n_mc=250, seed=4, schemes=(s,), f_mod=fm)[s] for fm in fmods]
    ax.loglog(fmods, y, color=SCHEME[s][1], lw=2.4, marker='o', ms=4, label=SCHEME[s][0])
dcref = ls.sigma_D(fk_d, n_mc=250, seed=4, schemes=('dc',))['dc']
ax.axhline(dcref, color=cDC, lw=1.8, ls='--', label='DC sweep (no modulation)')
ax.axvline(fk_d, color='0.5', lw=1.2, ls=':')
ax.text(fk_d * 0.9, dcref * 0.93, f'$f_{{knee}}$ = {fk_d:.0f} Hz', color='0.4',
        fontsize=8.5, va='top', ha='right')
ax.set_xlabel('Modulation frequency  $f_{mod}$  (Hz)')
ax.set_ylabel('$\\sigma_D$  (MHz)')
ax.set_title('(d)  Push $f_{mod}$ above $f_{knee}$; gain then saturates',
             loc='left', weight='bold', fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc='lower left')
ax.grid(True, which='both', alpha=0.18)

axi = ax.inset_axes([0.60, 0.58, 0.37, 0.38])
depths = np.array([0.25, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0])
yd = [ls.sigma_D(fk_d, n_mc=250, seed=5, schemes=('fm',), fm_depth=d)['fm'] for d in depths]
axi.plot(depths, yd, color=cFM, lw=1.8, marker='o', ms=3.5)
axi.plot(depths[int(np.argmin(yd))], min(yd), 'o', color='k', ms=6)
axi.set_xlabel('FM depth / $\\Delta\\nu$', fontsize=8)
axi.set_ylabel('$\\sigma_D$ (MHz)', fontsize=8)
axi.tick_params(labelsize=7); axi.grid(True, alpha=0.2)

plt.tight_layout()
out = 'lockin_decision.png'
plt.savefig(out, dpi=185, bbox_inches='tight')

# ------------------------------------------------------------- summary -----
print(f'saved {out}')
print(f'{"f_knee":>9} {"DC":>9} {"ON-OFF":>9} {"FM":>9}   sigma_D (MHz)')
for i, fk in enumerate(knees):
    print(f'{fk:9.1f} {res["dc"][i]:9.4f} {res["onoff"][i]:9.4f} {res["fm"][i]:9.4f}')
print(f'\nbreak-even f_knee : on-off {be["onoff"]:.0f} Hz,  FM {be["fm"]:.0f} Hz')
print(f'gain at f_knee=1 kHz : DC/FM = '
      f'{np.interp(1000., kx, dc) / np.interp(1000., kx, res["fm"]):.1f}x')
print(f'best FM depth        : {depths[int(np.argmin(yd))]:.2f} x linewidth')
