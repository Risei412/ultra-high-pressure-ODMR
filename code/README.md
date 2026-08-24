# NV high-pressure ODMR — lock-in sensitivity model

Numerical code behind the excitation-wavelength / sensitivity analyses
(green vs blue vs mix; blue-wavelength optimisation at 100 and 120 GPa;
laser-power dependence; literature reproduction).

Physics anchored to: K. O. Ho, C. Dailledouze, V. Žalandauskas, *et al.*,
"Optical Stability and Photophysics of NV Centers in Diamond up to 120 GPa",
arXiv:2606.02399 (2026). Rate constants (`a_gs, a_es, r0, rbg, w0`) are
**phenomenological** and are the quantities to be calibrated from the
`(I_405, I_457)` intensity sweep.

## Files
- `nv_model.py` — shared model: absorption cross section (finite-temperature
  Franck–Condon / Struck–Fonger envelope, C-1), ZPL from its two measured
  anchors (C-2), ISC contrast prefactor (C-3), stress anisotropy of the anvil
  geometry (C-4), emission collection efficiency (C-7), steady-state NV⁻
  fraction `f_minus`, lock-in sensitivity `eta ∝ Δν/(C√R)`, Monte-Carlo helper.
- `nv_model_power.py` — intensity-explicit extension `NVModelPower`
  (3E saturation, two-photon ionisation, power broadening); `u = 1` is defined
  as the intensity that half-saturates NV⁻ at `lambda_opt(120 GPa)`.
- `fig1_green_blue_mix.py` — sensitivity of green(532)/blue(457)/mix vs pressure
  → `sensitivity_green_blue_mix.png`.
- `fig2_blue_wavelength_sweep.py` — blue-wavelength sweep at a chosen pressure
  (optimum, mechanism, optimum-vs-pressure)
  → `blue_wavelength_sensitivity_<P>GPa.png`.
- `fig3_power_sweep.py` — power-explicit sweep: `eta(lambda,u)` map, ridge
  `lambda_opt(u)` with MC band, `eta(u)` for real laser lines
  → `power_wavelength_sensitivity_<P>GPa.png`.
- `analysis_C_lambda.py` — test of the assumption `C` is wavelength independent
  (stress-split 3E orbital branches); prints Q1–Q3 and the threshold splitting.
- `repro_literature.py` — reproduction of 26 published observations with the
  frozen model (plus 1 recorded as OPEN); prints the commercial-line table,
  the anvil-geometry table and the MC uncertainty of `lambda_opt`.
- `fig4_tornado.py` — anchor sensitivity of `lambda_opt` at 120 GPa (PLAN A3):
  two measured optical quantities move it, the nine phenomenological constants
  move it by exactly zero → `tornado_lambda_opt_120GPa.png`.
- `fig5_threshold.py` — `eta(532)/eta(473)` vs pressure: the recommendation
  reverses at ~71 GPa, which is the falsifiable prediction and the reason the
  50 GPa null result of [Bha22] is a reproduction → `threshold_green_blue.png`.
- `fig6_three_shifts.py` — talk figure: the absorption envelope (586 → 474 nm)
  and the ionisation edge (463 → 405 nm) closing on the window
  → `talk_three_shifts.png`.
- `fig7_answer_talk.py` — talk figure: single-panel `eta(lambda)` with the 5%
  tolerance window and the cost of each commercial line
  → `talk_answer_<P>GPa.png`.
- `talk_style.py` — shared style for the four presentation figures (one message
  per figure; only 405/473/532 nm are ever drawn together, because 457 and
  473 nm are perceptually the same colour).  Two styles: `--style talk`
  (default, wavelength-as-colour, for `slides/talk_deck.html`) and
  `--style st` (Institute of Science Tokyo deck style — colour encodes
  argumentative role, annotation budget 0–3, larger type; writes `*_st.png`
  for `slides/build_pptx.py`).
- `tests/test_freeze.py` — freeze tests (bit-exact legacy T=0 golden values,
  envelope physics, invariance of `lambda_opt`, C-4/C-7 and power-model
  regressions) plus the literature suite.

## Requirements
Python ≥ 3.9 with `numpy` (≥ 2.0 — `np.trapezoid` is used), `scipy`,
`matplotlib`, and `pytest` to run the tests (see `requirements.txt`).

```bash
pip install -r requirements.txt
```

## Reproduce the figures and checks
```bash
python fig1_green_blue_mix.py                     # green/blue/mix vs pressure
python fig2_blue_wavelength_sweep.py              # blue sweep @120 GPa (compare 100,140)
python fig2_blue_wavelength_sweep.py 100 75 125   # blue sweep @100 GPa (compare 75,125)
python fig3_power_sweep.py                        # power dependence @120 GPa
python analysis_C_lambda.py                       # C(lambda) assumption
python fig4_tornado.py                            # anchor sensitivity of lambda_opt
python fig5_threshold.py                          # green/blue crossover pressure
python fig6_three_shifts.py                       # talk: the two closing edges
python fig7_answer_talk.py                        # talk: the answer, single panel
python fig6_three_shifts.py --style st            # Science Tokyo variants (*_st.png)
python repro_literature.py                        # literature reproduction table
python -m pytest tests/ -q                        # freeze tests
```

## Key numbers reproduced
(room temperature, micropillar geometry `alpha = 0.95`, fixed excitation)

| Quantity | Value |
|---|---|
| optimal blue λ @120 GPa | 474.0 nm (MC 16–84%: +3.7 / −4.3 nm) |
| optimal blue λ @100 GPa | 484.1 nm |
| 5% / 10% tolerance window @120 GPa | 461.6–486.2 nm / 456.6–491.1 nm |
| η/η_opt @120 GPa | 473 nm: 1.00, 475 nm: 1.00, 457 nm: 1.10, 487 nm: 1.06, 532 nm: 2.66 |
| η/η_opt @100 GPa | 487 nm: 1.00, 473 nm: 1.04, 457 nm: 1.25 |
| green→blue crossover (473 nm) | 71 GPa (457 nm: 82 GPa, 450 nm: 87 GPa) |
| optimal blue λ tracking | ~0.5 nm/GPa near 100–120 GPa (0.62 nm/GPa averaged 40–140 GPa) |
| literature suite | 26/26 reproduced, 1 OPEN |

## Model summary
```
eta ∝ Δν / (C √R)                      # lower = better
  C = C0(P) · f₋/(f₋ + w0(1−f₋))       # ISC prefactor x NV0 dilution
  R = f₋ · σ_abs · η_col(P)             # detected photon rate (fixed excitation)
  f₋ = G_rec/(G_rec+G_ion)             # steady-state NV⁻ fraction
    G_ion = a_gs·ReLU(Eγ − IP(³A₂)) + a_es·σ_abs
    G_rec = r0·σ_abs + rbg
σ_abs : finite-T Franck–Condon envelope with ZPL(P), S_abs(P) from the reference
η_col : fraction of NV⁻ emission inside the fixed detection passband
```

## Known limitations
The claim these figures support is frozen in `docs/claim_freeze.md`.
See `docs/code_audit_2026-08.md` for the full audit. The open items that touch
numbers quoted in the paper are: the excitation weight `I` is a photon flux, not
an optical power (a fixed-power convention moves `lambda_opt` +1.6 nm and the
crossover +2.3 GPa); the ZPL Debye–Waller line is added with the tabulated DWF
as its PEAK, not its AREA, so the effective DWF is 0.58x the nominal value; and
all pressure-dependent anchors except the ZPL are clipped at 120 GPa, so results
above 120 GPa are extrapolations.
