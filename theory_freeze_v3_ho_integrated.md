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

---

# Addendum A1 — general coincidence/divergence theory

Date: 2026-08-27. Appended under the v3 immutability rule as a **theory-only**
extension. It introduces no new measurement, no fitted value, and no change to
the frozen optical-limit table above. Every proposition below is derived from
the equations already frozen in *Conditional ODMR layer*; full derivations are
in `docs/theory_optima_coincidence.md`.

The addendum answers one question the frozen body left open: **under what
condition does the wavelength that maximises excitation coincide with the
wavelength that optimises ODMR sensitivity, and how far apart are they when it
does not?**

## A1.0 Definitions

At fixed pressure and fixed incident optical power,

\[
A(\lambda)=\lambda\,\sigma_{\rm abs}^{\rm Ho}(hc/\lambda,P),\qquad
\Gamma_p=\gamma\,I\,A(\lambda),\qquad
G\equiv\frac{C}{\Delta\nu},\qquad
\Phi\equiv G^2R,
\]

\[
\lambda_{\rm abs}=\arg\max A,\qquad
\lambda_{\rm PL}=\arg\max R,\qquad
\lambda_\eta=\arg\min \eta=\arg\max \Phi .
\]

Write \(\ell_f=\mathrm{d}\ln f/\mathrm{d}\lambda\) and
\(\kappa_R=-\mathrm{d}^2\ln R/\mathrm{d}\lambda^2\big|_{\lambda_{\rm PL}}>0\).
\(\lambda_{\rm abs}\) is a theoretical quantity; only \(\lambda_{\rm PL}\) and
\(\lambda_\eta\) are measurable, and they differ in what must be recorded:
\(\lambda_{\rm PL}\) needs the count rate alone, \(\lambda_\eta\) needs rate,
contrast and linewidth separately.

## A1.1 Frozen propositions

**P1 (coincidence condition).**
\(\lambda_\eta=\lambda_{\rm PL}\iff\ell_G(\lambda_{\rm PL})=0\).
Since \(\eta=1/(G\sqrt R)\), \(\lambda_\eta\) maximises \(G^2R\) while
\(\lambda_{\rm PL}\) maximises \(R\); at \(\lambda_{\rm PL}\) the stationarity of
\(\Phi\) reduces to that of \(G\).
*Corollary.* If \(C/\Delta\nu\) is wavelength independent, the two optima
coincide **exactly, at every pressure and every power**, no matter how strongly
\(q_-\) or \(s_{\rm sat}\) depend on wavelength. Conversely, an observed
separation is direct evidence of a wavelength-dependent \(C/\Delta\nu\) and
cannot be produced by the detected rate alone.

**P2 (split magnitude).**
\(\lambda_\eta-\lambda_{\rm PL}\simeq 2\,\ell_G(\lambda_{\rm PL})/\kappa_R\).
The 120 GPa published curve gives, in the low-power limit,
\(\kappa_R\to\kappa_A\simeq 8\times10^{-4}\,\mathrm{nm^{-2}}\)
(\(9\times10^{-4}\) on the blue flank, \(7\times10^{-4}\) on the red flank),
hence a frozen conversion factor

\[
\lambda_\eta-\lambda_{\rm PL}\simeq 2.5\times10^{3}\,[\mathrm{nm^2}]\times\ell_G .
\]

The 5% penalty band is \(440.65\pm15.6\) nm, so \(\lambda_\eta\) leaves that band
once \(|\ell_G|>0.62\,\%/\mathrm{nm}\) — equivalently, once \(C/\Delta\nu\) varies
by more than about \(\pm28\%\) across a 430–470 nm scan. **The flatness of the
absorption band that makes the optical optimum soft is the same flatness that
makes the sensitivity optimum easy to displace.**

**P3 (mediation theorem).** Assume (M): \(q_-,C,\Delta\nu,s_{\rm sat}\) depend on
wavelength **only** through \(\Gamma_p\propto I\,A(\lambda)\), i.e. only on how
many photons are absorbed and not on which photon. Then \(\eta=h(\Gamma_p)\) and

- (a) \(R\) increasing in \(\Gamma_p\) \(\Rightarrow\lambda_{\rm PL}=\lambda_{\rm abs}\) at all powers;
- (b) \(\Phi\) increasing in \(\Gamma_p\) \(\Rightarrow\lambda_\eta=\lambda_{\rm PL}=\lambda_{\rm abs}\);
- (c) \(\Phi\) with an interior maximum at \(\Gamma_p^\*<\gamma I A_{\max}\)
  \(\Rightarrow\lambda_\eta\) is **not a point but an exactly degenerate doublet**
  \(\lambda_\pm\) solving \(A(\lambda_\pm)=\Gamma_p^\*/(\gamma I)\), straddling
  \(\lambda_{\rm abs}\) with \(\eta(\lambda_-)=\eta(\lambda_+)\); \(\lambda_{\rm abs}\)
  is then the locally *worst* choice between them;
- (d) the doublet opens at \(I_c=\Gamma_p^\*/(\gamma A_{\max})\) and closes as \(I\to0\),
  with separation \(2\sqrt{2\ln(A_{\max}/A^\*)/\kappa_A}\) — 32 nm already at
  \(A^\*/A_{\max}=0.9\), hence resolvable with the 445/457/473 nm lines.

**P4 (what can and cannot split the optima).** An interior maximum of \(\Phi\)
requires

\[
\frac{\mathrm{d}\ln(C/\Delta\nu)}{\mathrm{d}\ln\Gamma_p}
<-\frac12\frac{\mathrm{d}\ln R}{\mathrm{d}\ln\Gamma_p}.
\]

With \(\rho_i=\Gamma_p/\Gamma_i\): rate saturation alone gives
\(\mathrm{d}\ln\Phi/\mathrm{d}\ln\Gamma_p=1/(1+\rho_s)>0\) and **cannot** split;
linewidth power broadening alone gives \(1/(1+\rho_c)>0\) and **cannot** split;
saturation together with power broadening **can**; contrast collapse
\(C=C_0/(1+\rho_C)\) alone **can**, at \(\Gamma_p^\*=\Gamma_C\), because \(\eta\)
is first order in \(C\) and only half order in \(R\).

> This upgrades the frozen statement *"with flat non-optical factors, saturation
> broadens/flattens the optimum but does not move it"* from a numerical
> observation to a theorem, and identifies its cause: **saturation does not move
> the optimum; contrast does.**

**P5 (ordering).** While \(G\) decreases with \(\Gamma_p\),
\(\mathrm{d}\ln\Phi/\mathrm{d}\ln\Gamma_p<\mathrm{d}\ln R/\mathrm{d}\ln\Gamma_p\),
so \(\Phi\) turns over at lower power than \(R\). Raising the power therefore
produces, in order: (i) coincidence at low power; (ii) an **intermediate window
in which the PL spectrum is still single-peaked at \(\lambda_{\rm abs}\) while the
sensitivity optimum has already split**; (iii) turnover of the PL itself. Window
(ii) is where the central proposition of this work is observable.

**P6 (mechanism discrimination).** A split produced under (M) is symmetric,
degenerate, and closes as \(I\to0\). A split produced by a violation of (M) — a
channel that depends on photon energy rather than on absorbed photon number — is
one-sided and survives \(I\to0\). The (M)-violating channels identifiable at
120 GPa are: the ground-state ionisation edge IP(\(^3A_2\)) = 3.06 eV =
405.2 nm (active only below 405 nm, so it sets the blue wall but does not move a
440 nm optimum); the **NV\(^0\) absorption band** entering recombination, whose
envelope differs from that of NV\(^-\) (the only candidate active inside the blue
window); the detection passband; and sub-ZPL hot-band absorption above 514.5 nm
(see erratum E1).

**P7 (the charge state alone cannot split the optima).** In the frozen charge
model both ionisation and recombination scale with the same
\(\sigma_{\rm abs}\), so \(f_-\) is (M)-mediated and, by P3(b), cannot separate
\(\lambda_\eta\) from \(\lambda_{\rm PL}\). This is the structural reason for the
previously reported flatness of \(f_-\) to 0.2% across the blue window. Charge
conversion can move the sensitivity optimum **only** if
\(\sigma_{\rm abs}^{\rm NV^0}\) has a spectral shape different from
\(\sigma_{\rm abs}^{\rm NV^-}\).

## A1.2 Pre-registered tests

None requires absolute calibration.

| test | procedure | discriminates |
|---|---|---|
| T1 collapse | plot \(\eta\) against measured \(R\) for all \((\lambda,I)\), discarding \(\lambda\) | collapse onto one curve ⟹ (M) holds; scatter measures \(\ell_G\) |
| T2 power closure | measure \(\lambda_\eta(I)-\lambda_{\rm PL}(I)\) | closes as \(I\to0\) ⟹ intensity-mediated; persistent offset ⟹ (M) violated |
| T3 degeneracy | compare \(\eta(\lambda_-)\) and \(\eta(\lambda_+)\) if a doublet appears | equal ⟹ pure (M); imbalance gives \(\ell_G\) |
| T4 single-line prediction | from a 457 nm power sweep alone extract \(\Gamma_C,\Gamma_{\rm sat},\Gamma_c\Rightarrow\Gamma_p^\*\) | predicts \(I_c\) and the doublet separation **without any wavelength scan** |

T4 is the strongest pre-registered item: it converts a single-wavelength
measurement into a quantitative prediction of the multi-wavelength sensitivity
scan.

## A1.3 Falsification of this addendum

- If T1 collapses onto a single curve **and** \(\lambda_\eta=\lambda_{\rm PL}\)
  holds within \(\pm3\) nm at every power, then no divergence exists at 120 GPa,
  the conventional practice of maximising PL is correct, and the central
  proposition of this work fails.
- If the T4 prediction misses the measured \(\lambda_\eta(I)\) by more than a
  factor of two, the functional form \(\Phi=G^2R\) itself is in question, not
  merely its parameters.
- P2's conversion factor is falsified if the measured \(\kappa_R\) at low power
  differs from \(8\times10^{-4}\,\mathrm{nm^{-2}}\) by more than a factor of two
  after detection corrections.

---

# Erratum E1 — sub-ZPL interpolation artefact at 532 nm

Date: 2026-08-27. Filed under the v3 bug-fix/erratum rule. It does not change
any equation or any value at \(\lambda\le500\) nm.

The 120 GPa published curve carries 50 extracted points. Its lowest-energy real
sample is at 2.396 eV (517.4 nm), just below the reconstructed ZPL at
2.410 eV (514.5 nm); the only other point below that is a single baseline sample
at 1.400 eV (885.6 nm, the figure's left axis edge). The interval between them
is empty.

**532 nm (2.3305 eV) falls inside that empty interval.** The frozen values
`relative rate 0.00637`, `penalty 12.53`, and `advantage 11.99` are therefore
produced by linear interpolation between the axis-edge baseline point and the
ZPL onset, not by extracted data. Physically, absorption below the ZPL vanishes
at 0 K and is thermally activated at finite temperature,

\[
\sigma_{\rm abs}(E<E_{\rm ZPL})\propto
\exp\!\left(-\frac{E_{\rm ZPL}-E}{k_BT}\right),\qquad
E_{\rm ZPL}-E_{532}=79.5\ \mathrm{meV},
\]

giving \(4.6\times10^{-2}\) at 300 K and \(3.7\times10^{-5}\) at 90 K. The
interpolated 6.4e-3 lies between these two without deriving from either.

**Consequence.** The 120 GPa 457/532 comparison is governed by hot-band
absorption and temperature, not by the optical kernel. The ratio ×11.99 must not
be quoted as a single number; it is to be treated as a temperature-dependent
quantity, or the discussion restricted to \(\lambda\le517\) nm where extracted
data exist. This strengthens, and supersedes, the qualitative caveat already
recorded in *Optical-limit result at 120 GPa*. The optimum 440.65 nm, the
×1.04494 penalty at 457 nm, and all of Addendum A1 are unaffected.

