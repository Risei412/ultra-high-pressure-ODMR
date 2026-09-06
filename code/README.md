# Numerical code

Numerical models, reconstructed spectral data, tests, and figure generators
for the high-pressure NV ODMR manuscript.

The optical inputs are anchored to K. O. Ho et al., "Optical Stability and
Photophysics of NV Centers in Diamond up to 120 GPa," arXiv:2606.02399 (2026).
Phenomenological rate constants remain calibration targets rather than measured
high-pressure inputs.

## Contents

- `ho_spectrum_model.py` — reconstructed published absorption kernel.
- `ho_odmr_sensitivity.py` — optical-limit ODMR sensitivity model.
- `nv_model.py`, `nv_model_power.py` — legacy phenomenological and
  intensity-explicit models.
- `extrapolation_bounds.py` — what the published record leaves undetermined
  above 120 GPa.
- `theory_a1_generalization.py` — coincidence and divergence analysis.
- `theory_a2_multiplicity.py` — level-set multiplicity and gauge degeneracy.
- `theory_a3_branch_exchange.py` — pressure-driven branch exchange.
- `data/` — digitized or reconstructed numerical inputs.
- `tests/` — regression and figure tests.
- `outputs/` — regenerable previews and reports; ignored by Git.

## [111]-aligned anisotropic-stress calculation

The literature-anchored extension is intentionally split into four parts:

- `nv_model_111.py`: stress tensor, low-power optical model, and normalized
  finite-power optimum sets.
- `report_111_anisotropic.py`: command-line report for `sigma_zz > sigma_xx =
  sigma_yy` with the DAC axis parallel to [111].
- `fig_alpha_invariance.py`: experimental-handoff plot showing that the stress
  tensor changes with `alpha` while the current hydrostatic optical baseline
  remains approximately 441 nm.
- `tests/test_nv_model_111.py`: tensor, approximately 441-nm hydrostatic
  baseline, and power-split
  regression tests.
- `../docs/theory/nv111_literature_reconstruction.md`: literature provenance
  and the boundary between reproduced and conditional claims.

Run the standard 120-GPa report from `code/`:

```powershell
python report_111_anisotropic.py --pressure 120 --alpha 0.56
```

`alpha = sigma_xx/sigma_zz = sigma_yy/sigma_zz`; the finite-power axis is
reported as `I/Ic` because the absolute `Ic` requires calibration in the actual
[111] DAC.

The source modules remain in one directory because the analyses import each
other as flat modules. Do not split them into subpackages without updating the
imports and tests together.

## Above 120 GPa

Every optical anchor stops at 120 GPa. `ho_spectrum_model.py` refuses outright
past that edge; `nv_model.py` does not, and used to continue by freezing each
anchor inline, so a call at 400 GPa returned an asymptote that looked like a
prediction. That clip is now the default value of an explicit policy (C-8):

```python
NVModel(extrapolate='clip')    # default; freezes the anchors, as before
NVModel(extrapolate='linear')  # continues the fitted laws -- an extrapolation
NVModel(extrapolate='error')   # refuses, as ho_spectrum_model.py does
```

Every frozen number is unchanged, because `'clip'` is what the model already
did. `NVModel.lambda_opt` now also warns when the optimum is pinned to an edge
of its search window instead of returning the edge silently.

`extrapolation_bounds.py` quantifies what the choice costs:

```powershell
python extrapolation_bounds.py
```

The optimum is undetermined by 21 nm at 120 GPa, 58 nm at 200 GPa and 134 nm at
400 GPa, and the spread at 120 GPa is not zero because the effective phonon
energy carries no pressure dependence at all. One statement survives: the
one-photon ionisation threshold is a one-sided bound,
`lambda_opt >= hc/IP(3A2)`, insensitive to the uncalibrated `a_gs` over three
decades — but only where the ZPL and `IP(3A2)` are continued under the same
law. Do not quote a single optimal wavelength above 120 GPa.

## Environment and tests

Run from this directory:

```powershell
python -m pip install -r requirements.txt
python -m pytest . -q
```

## Current manuscript figures

```powershell
python fig_m_level_set_ladder.py
python fig_g_gauge_degeneracy.py
python fig_x_branch_exchange.py
python esa_figv_kernel_sanity.py
```

The four scripts produce five figures: `esa_figv_kernel_sanity.py` draws the
cross-figure check and the internal checks K1–K4 as two separate files.

Vector PDFs are written to `../paper/figures/`; PNG previews are written to
`outputs/figures/`.

All five are full-width `figure*` floats, included at `\textwidth` ≈ 7.06 in,
so each is drawn at that width and nothing is rescaled on the page. The shared
style lives in `fig_style.py` — 8.5 pt axis labels, 8 pt ticks, 7 pt legends
and in-axes text. **Set font sizes there, not in a figure script**: a figure
drawn at some other canvas width is rescaled by the ratio, lettering included,
which is how the kernel figure once reached the page at 4.6 pt.

Numbers belong in the captions. Each script prints the numbers its caption
quotes; none is rendered into the raster.

## Legacy exploratory figures

```powershell
python fig1_green_blue_mix.py
python fig2_blue_wavelength_sweep.py
python fig3_power_sweep.py
```

These scripts write PNG files to `outputs/figures/`. Their v1 numerical claims
are historical and should not replace the frozen Ho-integrated results.
