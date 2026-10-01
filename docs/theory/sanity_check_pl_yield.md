# Ho et al. Fig. 5(b) reproduction audit

## Scope

This audit keeps three statements separate:

1. reproduction of Ho et al.'s calculated absorption curves;
2. conditional calibration to the measured PL-yield curves;
3. the engineering choice of a laser line for a future 120 GPa experiment.

Passing item 2 does not imply item 1. Reconstructing published curves also does
not constitute an independent first-principles calculation.

## Observables

The CSV contains experimental PL series (`expt532`, `expt457`) and
digitised Ho-theory series (`theory532_ho`, `theory457_ho`).

- Ho theory is compared with `sigma_abs` only.
- Experimental PL is reported against both `sigma_abs` and
  `sigma_abs * eta_col`.
- The experimental detection proxy may be selected only from the documented
  optical path or supplement, not from which choice gives a lower residual.
- Each wavelength receives one arbitrary multiplicative scale because the
  source series are separately normalised.

The Ho-theory score is evaluated only over the pressure interval sampled by the
corresponding experiment. This prevents an unmeasured extrapolation from
dominating the score.

## Reproduction gates

The executable audit uses the following deliberately generous limits:

| gate | shape requirement | peak requirement |
|---|---:|---:|
| Ho calculated absorption | fractional RMS ≤ 10% per branch | error ≤ 5 GPa |
| experimental PL shape | fractional RMS ≤ 16% per branch | error ≤ 10 GPa |

The experimental peak tolerance reflects the 16 GPa sampling gap around the
blue dome. The smooth theoretical curve retains the tighter 5 GPa tolerance.

Overall reproduction requires both gates. Run:

```bash
cd code
python repro_yield.py
python -m pytest tests/test_repro_yield.py tests/test_ho_spectrum_model.py -q
```

## Current result

At 300 K, the frozen model fails the 457 nm experimental branch. A conditional
fit using `sigma_abs * eta_col` gives approximately:

- `dE120 = 0.550 eV`;
- experimental fractional RMS ≈ 13% on both branches;
- experimental peak pressures ≈ 26 and 83 GPa;
- `lambda_opt(120 GPa) ≈ 450 nm`.

The same calibrated model does not reproduce Ho's calculated absorption:
the branch-wise discrepancies and peak offsets exceed the independent gate.
A one-parameter fit of `dE120` directly to the Ho curves also remains above
20% pooled fractional RMS. Therefore 550 meV is a conditional fit parameter,
not a measured ZPL shift, and the model-form question remains open.

Removing `eta_col` changes the best experimental fit materially (the fit runs
to the edge of the original parameter grid and the inferred optimum moves by
more than 15 nm). The detection-band assumption is therefore load-bearing and
must be verified before quoting 443–455 nm as a validated prediction.

### Published-curve reconstruction

`extract_ho_fig1e.py` extracts the absorption half of all seven vector curves
in the arXiv source figure `theory1ev_V4.pdf` at 0, 20, ..., 120 GPa. The
extracted spectra are stored in `data/ho_fig1e_absorption.csv`.
`HoPublishedSpectrumModel` interpolates energy on each published curve and
uses a cubic spline in pressure, matching the interpolation described for
Fig. 5(b).

Fig. 5(b) is not used to construct this model and therefore supplies a
cross-figure validation:

| comparison | 532 nm | 457 nm |
|---|---:|---:|
| fractional RMS against Ho Fig. 5(b) theory | 1.5% | 0.4% |
| reconstructed / reference dome peak | 17 / 17 GPa | 88 / 88 GPa |
| fractional RMS against experimental PL | 10.0% | 10.2% |

The published-curve model passes both gates and gives
`lambda_opt(120 GPa) = 440.7 nm`. At fixed optical power, 457 nm retains
91.6% of the maximum absorbed-photon rate, corresponding to only ×1.045 in
shot-noise-limited sensitivity.

## Interpretation

- **Cleared:** v1 is falsified on the measured 457 nm pressure dependence.
- **Cleared:** Ho's published Fig. 1(e) spectra reproduce the independently
  digitised Fig. 5(b) theory curves and experimental dome shapes.
- **Supported:** 457 nm is a practical line at 120 GPa; the reconstructed
  published spectrum puts the optimum at 440.7 nm with a ×1.045 penalty.
- **Not cleared:** independent reproduction from mode-resolved
  electron--phonon and dynamical-Jahn--Teller inputs.
- **Not identified:** which of `dE120`, effective phonon energy, Huang–Rhys
  evolution, or stress anisotropy accounts for the experimental shape.

The next microscopic extension requires Ho's supplemental mode-resolved
`S_k` and `K_k^2` data (or equivalent raw DFT outputs), which are absent from
the arXiv source archive. It should remain a separate branch. The frozen v1
output remains the historical pre-comparison prediction.

