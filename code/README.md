# NV high-pressure ODMR — lock-in sensitivity model

Numerical code behind the excitation-wavelength / sensitivity analyses
(green vs blue vs mix; blue-wavelength optimisation at 100 and 120 GPa).

Physics anchored to: K. O. Ho, C. Dailledouze, V. Žalandauskas, *et al.*,
"Optical Stability and Photophysics of NV Centers in Diamond up to 120 GPa",
arXiv:2606.02399 (2026). Rate constants (`a_gs, a_es, r0, rbg, w0`) are
**phenomenological** and are the quantities to be calibrated from the
`(I_405, I_457)` intensity sweep.

## Files
- `nv_model.py` — shared model: absorption cross section (low-T Franck–Condon /
  Pekarian envelope), steady-state NV⁻ fraction `f_minus`, and lock-in
  sensitivity `eta ∝ Δν/(C√R)`; also Monte-Carlo band helper.
- `fig1_green_blue_mix.py` — sensitivity of green(532)/blue(457)/mix vs pressure
  → `sensitivity_green_blue_mix.png`.
- `experiment_plan.py` — the two laser parameters for a run, from measurements
  the cell can supply: the excitation line from a measured ZPL, and the
  operating power from the knee of the saturation curve. Neither needs the
  uncalibrated intensity scale. See `docs/experiment_120GPa_laser_parameters.md`.
- `fig2_blue_wavelength_sweep.py` — blue-wavelength sweep at a chosen pressure
  (optimum, mechanism, optimum-vs-pressure)
  → `blue_wavelength_sensitivity_<P>GPa.png`.

## Requirements
Python ≥ 3.9 with `numpy`, `scipy`, `matplotlib` (see `requirements.txt`).

```bash
pip install -r requirements.txt
```

## Reproduce the figures
```bash
python fig1_green_blue_mix.py                     # green/blue/mix vs pressure
python fig2_blue_wavelength_sweep.py              # blue sweep @120 GPa (compare 100,140)
python fig2_blue_wavelength_sweep.py 100 75 125   # blue sweep @100 GPa (compare 75,125)
```

## Key numbers reproduced
| Quantity | Value |
|---|---|
| green/blue crossover (457 nm fixed) | ~86 GPa |
| optimal blue λ @100 GPa | 487 nm (η(457)/opt = 1.35) |
| optimal blue λ @120 GPa | 475.5 nm at equal optical power (473 nm penalty = 0.2%) |
| optimal blue λ tracking | ~0.6–0.7 nm/GPa (follows ZPL / sideband edge) |

## Model summary
```
eta ∝ Δν / (C √R)                      # lower = better
  C = C0 · f₋/(f₋ + w0(1−f₋))          # contrast, diluted by NV0 background
  R = f₋ · (I/Eγ) · σ_abs               # detected rate at fixed optical power
  f₋ = G_rec/(G_rec+G_ion)             # steady-state NV⁻ fraction
    G_ion = a_gs·ReLU(Eγ − IP(³A₂)) + a_es·σ_abs
    G_rec = r0·σ_abs + rbg
σ_abs : low-T Franck–Condon envelope with ZPL(P), S_abs(P) from the reference
```
