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
- `thesis_crosscheck.py` — T1-T3, the retrodictions against chapter 6 of
  Bhattacharyya's thesis (`docs/ref/Principle_and_Applications_of_.pdf`).
- `thesis_from_the_core.py` — N1-N4, the same chapter read out of the core
  (eta, assumption (M), Theorems M/G/X) rather than sorted against it.
- `anvil_transmission.py` — bounds the anvil transmission from Fig. 6.3(b).
- `external_audit.py` — X1-X5, what the theory explains in the two external
  sources and what it does not: the power-convention erratum E5, how far
  erratum E1's sub-ZPL gap actually reaches, chapter 7's cryogenic 532 nm
  megabar data against the hydrostatic kernel, and the band-position/anvil
  degeneracy.  Write-ups in `docs/ho_results_audit.md` and
  `docs/bhattacharyya_thesis_scope.md`.
- `geometry_layer.py` — the optimum as a function of culet geometry.  A single
  scalar `g` scales the band shift relative to Ho's quasi-hydrostatic
  micropillar; `g = 1` reduces the layer to the frozen kernel exactly (window,
  a_max, the four maxima and the 1.4414 ZPL rung).  Chapter 7 of the thesis
  bounds `g <= 0.741` for the [111] group in a [111] culet, so the optimum
  there is `>= 460.33 nm`, not 440.65 nm; Fig. 6.3(b) with `T <= 1` bounds
  `g >= 1.048` for a [100] culet.  The two signs are why `nv_model`'s C-4
  factor is NOT extended per NV group.  Write-up in
  `docs/geometry_of_the_optimum.md`.
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

## Key numbers

**Two chains, two answers.**  Every number below is labelled with the chain that
produced it.  Do not quote a v1 number as a v3 result -- they differ by 34.87 nm
at 120 GPa, more than twice the half-width of the 5 % tolerance band.

### v1 chain — phenomenological Franck–Condon envelope (`nv_model.py`)

Internally consistent, and pinned by `tests/test_freeze.py`, but **superseded as
a statement about the physics**: this envelope is unimodal at every pressure, so
it cannot represent the zero-phonon-line / sideband branch structure.

| Quantity | Value |
|---|---|
| green/blue crossover (457 nm fixed) | ~86 GPa |
| optimal blue λ @100 GPa | 487 nm (η(457)/opt = 1.35) |
| optimal blue λ @120 GPa | 475.5 nm (473 nm penalty = 0.2%) |
| optimal blue λ tracking | ~0.6–0.7 nm/GPa |

The ~86 GPa crossover is **not** the branch exchange of Addendum A3, despite
the numerical proximity to 87.9 GPa.  It is the geometry of one unimodal peak
sweeping past the midpoint of 457 and 532 nm (494.5 nm); the same test on the
Ho kernel lands at 51.4 GPa.  See `docs/novelty_and_exponent_audit.md` §4.

### v3 chain — Ho published kernel (`ho_spectrum_model.py`), frozen

| Quantity | Value |
|---|---|
| optimal λ @120 GPa | **440.65 nm** |
| 5 % tolerance band @120 GPa | **[426.43, 457.90] nm** (asymmetric) |
| penalty at 457 nm | ×1.04494 |
| penalty at 473 nm | ×1.2054 (**not** 0.2 %) |
| penalty at 514.5 nm (ZPL) | ×1.2006 |
| penalty at 532 nm | ×12.53 — but see erratum E1, it is an interpolation artefact |

### Anchors taken from published measurements

| Quantity | Value | Source |
|---|---|---|
| splitting exponent | E = 3, n = 2, ρ\* = 1/5 | Dréau et al., PRB 84, 195204 |
| S_abs vs pressure | 3.023 → 4.554 (+51 %, monotone) | source Fig. 1(b) |
| DWF_abs vs pressure | 0.0205 → 0.00226 (×9.07, monotone) | source Fig. 1(c) |

## Model summary (v1 chain)
```
eta ∝ Δν / (C √R)                      # lower = better
  C = C0 · f₋/(f₋ + w0(1−f₋))          # contrast, diluted by NV0 background
  R = f₋ · (I/Eγ) · σ_abs               # detected rate at fixed optical power
  f₋ = G_rec/(G_rec+G_ion)             # steady-state NV⁻ fraction
    G_ion = a_gs·ReLU(Eγ − IP(³A₂)) + a_es·σ_abs
    G_rec = r0·σ_abs + rbg
σ_abs : low-T Franck–Condon envelope with ZPL(P), S_abs(P) from the reference
```
