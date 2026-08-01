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
- `fig2_blue_wavelength_sweep.py` — blue-wavelength sweep at a chosen pressure
  (optimum, mechanism, optimum-vs-pressure)
  → `blue_wavelength_sensitivity_<P>GPa.png`.
- `background.py` — excitation-wavelength shapes of the spin-independent
  background collected inside the NV detection window (ruby R line, N3 /
  A-band, broad deformation luminescence; diamond-Raman exclusion check).
- `nv_bg.py` — `NVModel` + background: `eta = eta_0 · √(1+ρ)`, ρ = B/R
  (shot noise **and** contrast dilution), plus helpers for λ_opt, tolerance
  band and green/blue crossover.
- `fig5_background.py` — green-vs-blue re-assessment as a function of the
  background level ρ₀ → `background_green_blue_<P>GPa.png`.
- `tests/test_bg.py` — ρ₀=0 ↔ baseline equivalence, monotonicity, regressions.

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
python fig5_background.py                         # background re-assessment @120 GPa
python -m pytest tests/ -q                        # regressions
```

## Key numbers reproduced
| Quantity | Value |
|---|---|
| green/blue crossover (457 nm fixed) | ~86 GPa |
| optimal blue λ @100 GPa | 487 nm (η(457)/opt = 1.35) |
| optimal blue λ @120 GPa | 475 nm (η(487)/opt = 1.06, η(457)/opt = 1.12) |
| optimal blue λ tracking | ~0.6–0.7 nm/GPa (follows ZPL / sideband edge) |

### With spin-independent background (`fig5_background.py`, ρ₀ = B/R at 532 nm, ambient)
| Quantity | ρ₀ = 0 | ρ₀ = 1 | ρ₀ → ∞ |
|---|---|---|---|
| optimal blue λ @120 GPa | 474.8 nm | 478.3 nm | 480.9 nm |
| 5% tolerance band | 463–486 nm | 469–488 nm | 473–488 nm |
| crossover, blue = 457 nm | 86.2 GPa | 91.4 GPa | 93.0 GPa |
| crossover, blue = λ_opt | 75.2 GPa | 76.1 GPa | 76.5 GPa |
| η(blue λ_opt)/η(532) @120 GPa | 0.194 | 0.056 | 0.035 |

λ_opt moves ≲6 nm to the red and stays inside the zero-background 5% band;
the green/blue decision **never flips at 120 GPa** for any background level.

## Model summary
```
eta ∝ Δν / (C √R)                      # lower = better
  C = C0 · f₋/(f₋ + w0(1−f₋))          # contrast, diluted by NV0 background
  R = f₋ · σ_abs                        # detected photon rate (fixed power)
  f₋ = G_rec/(G_rec+G_ion)             # steady-state NV⁻ fraction
    G_ion = a_gs·ReLU(Eγ − IP(³A₂)) + a_es·σ_abs
    G_rec = r0·σ_abs + rbg
σ_abs : low-T Franck–Condon envelope with ZPL(P), S_abs(P) from the reference
```
