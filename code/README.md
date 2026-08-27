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
- `ho_spectrum_model.py`, `ho_odmr_sensitivity.py`,
  `report_120gpa_sensitivity.py` — the v3 external-kernel chain
  (`σ_abs^Ho → R_det → η`) behind `theory_freeze_v3_ho_integrated.md`.
- `theory_a1_generalization.py` — numerical execution of Addendum A1
  (the coincidence/divergence propositions P1–P7 and tests T1–T4) against the
  frozen Ho kernel.  Findings written up in
  `docs/theory_a1_numerical_execution.md`.
- `fig5_a1_generalization.py` — figures for the above
  → `a1_generalization_120GPa.png`.
- `theory_a2_multiplicity.py` — Addendum A2, the kernel-independent structural
  layer: the multiplicity-ladder theorem (the optimal set is a level set of `A`,
  whose size steps at the critical values of `A`, at power ratios needing no
  absolute calibration) and the gauge-degeneracy theorem (η sees the response
  only through `E = 2c + s + 2w`).  Write-up in `docs/theory_a2_multiplicity.md`.
- `fig6_a2_multiplicity.py` — figures for A2 → `a2_multiplicity_120GPa.png`.
- `theory_a3_branch_exchange.py` — Addendum A3, pressure-driven branch exchange:
  the ZPL and phonon-sideband branches scale differently with the Huang-Rhys
  factor, so pressure exchanges which one carries the global optimum.  Branches
  are identified in the **raw extracted samples**, not in the pressure
  interpolation.  Write-up in `docs/theory_a3_branch_exchange.md`.
- `fig7_a3_branch_exchange.py` — figures for A3 → `a3_branch_exchange.png`.
- `figure_validation.py` — pixel-level check of the extracted kernel against the
  source figure.  The sideband branch reproduces it to better than 1 %; the
  zero-phonon-line peaks turn out to be **clipped by the axis**, so A3's
  ×6.04 collapse and P* = 87.9 GPa are withdrawn (erratum E3).  Supplies the
  corrected, bandwidth-dependent treatment using the published Debye–Waller
  factor from `data/ho_fig1_panels_bc.csv`.
- `dreau_exponent.py` — reads A2's splitting exponent off the published CW-ODMR
  model of Dréau et al., PRB 84, 195204 (2011): `E = 3`, pump nonlinearity
  `n = 2`, so `E n = 6 > 1` and `rho* = 1/5`, independent of the microwave
  setting.  A2's splitting antecedent is therefore measured, not assumed.
  Audit in `docs/novelty_and_exponent_audit.md`.

Note that the v1 Franck–Condon envelope (`nv_model.py`) and the v3 Ho kernel
(`ho_spectrum_model.py`) do **not** agree on the optimum: 475.51 nm against
440.64 nm at 120 GPa.  `theory_a1_generalization.py` section S8 quantifies the
gap; do not mix numbers from the two chains without saying which one produced
them.

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
