# Theory freeze v3 — Ho-integrated 120 GPa ODMR sensitivity

Date: 2026-08-27

## Frozen scope

This freeze adopts the pressure-dependent absorption spectra published by Ho
et al. as an **external optical kernel** and extends them into an ODMR
sensitivity framework:

\[
\sigma_{\rm abs}^{\rm Ho}(\lambda,P)
\rightarrow R_{\rm det}(\lambda,P,I)
\rightarrow \eta(\lambda,P,I)
=\frac{\Delta\nu}{C\sqrt{R_{\rm det}}}.
\]

It does not claim an independent DFT or dynamical-Jahn--Teller calculation.
The Ho kernel is frozen because its Fig. 1(e) reconstruction independently
reproduces both branches of Fig. 5(b). The ODMR layer is frozen only at the
level of equations, factor definitions, limiting case, and falsification
criteria. Unmeasured wavelength-dependent charge, contrast, linewidth, and
saturation factors are not assigned fitted values.

## Optical-limit result at 120 GPa

Assumptions: fixed incident optical power, low saturation, and wavelength-
independent charge yield, ODMR contrast, linewidth, and collection within the
scan. Under exactly these assumptions:

| quantity | frozen result |
|---|---:|
| optimal wavelength | 440.65 nm |
| sensitivity penalty at 457 nm | ×1.04494 |
| sensitivity penalty at 532 nm | ×12.53 |
| sensitivity advantage, 457 over 532 | ×11.99 |
| minimum C(457)/C(532) for 457 to remain better, equal linewidth/yield | 0.0834 |

Thus 457 nm is a robust practical line in the optical limit: it retains 91.6%
of the maximum absorbed-photon rate, and its contrast could be about twelve
times lower than the 532 nm contrast before the optical advantage is erased
(if linewidth and charge yield are equal).

The ×11.99 blue/green sensitivity advantage is **not yet an experimental
prediction**. It is the optical-kernel limit. The 532 nm point lies in the very
weak tail of the reconstructed 120 GPa absorption curve, making that ratio
sensitive to unmeasured hot-band absorption and non-optical response.

## Conditional ODMR layer

The executable model defines

\[
A(\lambda,P)=\lambda\sigma_{\rm abs}^{\rm Ho}(hc/\lambda,P),
\]

\[
R_{\rm det}=q_-(\lambda,P)
\frac{x}{1+x/s_{\rm sat}(\lambda,P)},\qquad
x=I_{\rm rel}\frac{A(\lambda,P)}{A_{\max}(P)},
\]

\[
\eta=\frac{\Delta\nu(\lambda,P)}
{C(\lambda,P)\sqrt{R_{\rm det}}}.
\]

Here `q_-` collects NV- population, radiative yield, and wavelength-dependent
detection effects. `C`, `Delta nu`, and `s_sat` are replaceable response
functions. They are deliberately not inferred from Ho's PL-yield curves,
because those curves do not identify ODMR contrast or linewidth.

## What is now cleared

- Ho Fig. 1(e) → Fig. 5(b) cross-figure absorption reproduction passes.
- The optical-limit optimum at 120 GPa is 440.65 nm.
- 457 nm is within 4.5% of the optical-limit best sensitivity.
- The theory provides a measurable robustness condition, rather than assuming
  wavelength-independent ODMR contrast silently.

## What is not cleared

- Absolute sensitivity in T/sqrt(Hz).
- The actual 457/532 sensitivity gain at 120 GPa.
- A power-dependent shift of the optimum. With flat non-optical factors,
  saturation broadens/flattens the optimum but does not move it; a shift needs
  measured wavelength-dependent charge, contrast, linewidth, or saturation.
- Independent microscopic reproduction of Ho's mode-resolved calculation.

## Experimental freeze: minimum decisive measurements

1. At 120 GPa, measure 457 and 532 nm at matched incident power and record
   detected rate, ODMR contrast, and linewidth separately.
2. At 120 GPa, scan at least 405/445/457/473/488/532 nm below saturation to
   locate the sensitivity optimum, not merely the PL maximum.
3. Repeat the wavelength scan at two or more powers to test whether the
   optimum moves. A shift is evidence for wavelength-dependent non-optical
   response; no shift supports the optical-kernel limit.
4. Record the emission spectrum and NV charge-state proxy concurrently, so a
   change in detected counts is not misassigned to absorption.

Primary pre-registered outcomes are: `lambda_opt`, `eta(457)/eta_opt`,
`eta(532)/eta(457)`, and the separate wavelength dependences of rate, contrast,
and linewidth.

## Novelty and publication gate

| result | novelty assessment | publication role |
|---|---|---|
| Reproducing Ho's curves | validation, not novelty | methods/audit |
| 440.65 nm optical optimum from Ho's curve | derivative result | design prior |
| Absorption-to-ODMR sensitivity framework with explicit nuisance factors | moderate | theory core |
| Measured separation of PL optimum and sensitivity optimum | strong if observed | main result |
| Measured 457/532 sensitivity gain at 120 GPa | strong | main result |
| Power-dependent movement of sensitivity optimum | strong if reproduced and explained | main result |

For PRA or a specialist journal, the minimum credible package is the frozen
framework plus a 120 GPa 457/532 ODMR comparison reporting rate, contrast, and
linewidth. A multi-wavelength, multi-power scan makes the novelty materially
stronger. The most defensible central proposition is:

> Under megabar pressure, the excitation wavelength that maximizes absorbed
> or detected photons need not minimize ODMR sensitivity once charge
> conversion, contrast, linewidth, and saturation are included.

## Falsification and update rule

The optical-limit reduction is falsified if, below saturation, the measured
rate spectrum cannot be reconciled with the Ho absorption kernel after
independently measured charge/detection corrections. The full sensitivity
ranking is falsified if the measured `C sqrt(R) / Delta nu` ordering disagrees
with the model after those measured factors are inserted.

This v3 file is immutable as the pre-experiment record. Bug fixes receive an
erratum note. New measurements generate `v4`; they do not overwrite v3.

## Reproduction

```bash
cd code
python ho_odmr_sensitivity.py
python report_120gpa_sensitivity.py
python -m pytest tests/test_repro_yield.py \
  tests/test_ho_spectrum_model.py tests/test_ho_odmr_sensitivity.py -q
```
