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

---

# Erratum E2 — corrections to Addendum A1 from its numerical execution

Date: 2026-08-27. Filed under the v3 bug-fix/erratum rule. It corrects
Addendum A1 only. **No value in the frozen body changes**: the optical-limit
table was re-run and reproduces exactly (440.65 nm, ×1.04494, ×12.53, ×11.99,
0.0834), and the existing test suite passes.

Full working: `docs/theory_a1_numerical_execution.md`
(`code/theory_a1_generalization.py`, 23 regression tests).

**Provenance note.** `code/ho_odmr_sensitivity.py`,
`code/report_120gpa_sensitivity.py`, `code/tests/test_ho_odmr_sensitivity.py`
and `code/data/ho_120gpa_wavelength_scan.csv` — all required by the
*Reproduction* section above — were absent from the repository and present only
inside `ho_integrated_odmr_v3_20260827.zip`. They have been restored.

## E2.1 The reconstructed kernel is not a single band

A1 treats \(A(\lambda)\) as one flat band peaked at 440.65 nm. The
reconstruction carries **four** local maxima at 120 GPa:

| λ [nm] | \(a/a_{\max}\) | optical-limit penalty |
|---:|---:|---:|
| 440.64 | 1.000 | ×1.0000 |
| 475.55 | 0.660 | ×1.2306 |
| 500.19 | 0.277 | ×1.8994 |
| **514.46 (ZPL)** | **0.694** | **×1.2006** |

The ZPL penalty ×1.2006 is indistinguishable from that of 473 nm (×1.2054), and
514.46 nm lies inside the region §7 declares to hold real data (λ ≤ 517.4 nm),
sampled there at 0.1–0.5 nm spacing. It is absent from every candidate-line
list in the project.

## E2.2 P2's threshold and destination

A1 states that \(\lambda_\eta\) leaves the 5 % band once
\(|\ell_G|>0.62\,\%/\mathrm{nm}\). Exact argmax over the reconstruction gives a
**jump to the ZPL at \(\ell_G>0.268\,\%/\mathrm{nm}\)**, a displacement of
73.8 nm rather than 15.6 nm. P2 is a local expansion and cannot represent this.

P2's usable range is bounded on both sides: below
\(\ell_G\approx0.1\,\%/\mathrm{nm}\) the predicted shift is smaller than the
3–10 nm node spacing of the extraction and the argmax does not move at all;
above \(0.268\,\%/\mathrm{nm}\) the global optimum has already left the band.
In between the formula is good to a few per cent (at 0.4 %/nm, predicted
10.51 nm against exact 10.21 nm).

## E2.3 κ has two incompatible definitions

P2 defines \(\kappa_R=-\mathrm{d}^2\ln R/\mathrm{d}\lambda^2\) as a **local**
second derivative, but the quoted \(8\times10^{-4}\,\mathrm{nm^{-2}}\) is
reproduced only by the **chord** curvature
\(-2\ln a/(\lambda-\lambda_{\rm abs})^2\) (all seven rows of
`theory_optima_coincidence.md` §2 reproduce to 1 %). The true pointwise second
derivative spans \([-1.12,+1.14]\,\mathrm{nm^{-2}}\) with 20 sign changes and is
unusable.

A parabolic fit is window dependent: 7.01e-4 (±10 nm), 7.61e-4 (±15 nm),
8.07e-4 (±40 nm). **The quoted 8e-4, hence the 2.5e3 nm² conversion factor,
corresponds to a ±40 nm fit**; over the 5 % band it is 2628 nm². A1.3's
falsification clause is not well posed until the fit window is stated; read it
as "±40 nm parabolic fit".

## E2.4 Three geometric corrections

- **The 83 nm doublet row is unverifiable.** At \(a^*/a_{\max}=0.5\) the blue
  member lies below 402 nm, outside the reconstruction and below the
  ground-state ionisation edge where (M) fails by A1's own P6. The widest
  verifiable separation is 112.6 nm. The 0.9 and 0.7 rows reproduce well
  (32.80 and 60.55 nm against the Gaussian 33.28 and 61.24 nm).
- **The doublet is degenerate but not symmetric.** η-degeneracy is exact
  (2.2e-16), yet the red member sits ~3 nm further out at both levels. The 5 %
  band is likewise asymmetric: **[426.43, 457.90] nm**, not 440.65 ± 15.6 nm.
  T3 must be phrased as a search for equal η, not for symmetric placement.
- **P5's window is narrower than advertised.** The split spans 1.43 decades, but
  both members stay inside the data window over only **0.27 decades** (a factor
  1.9 in power).

## E2.5 P4's conclusion is too strong

A1 concludes "saturation does not move the optimum; contrast does." The
single-mechanism table is correct, but **saturation together with power
broadening, and contrast collapse alone, give identical \(\Phi(\Gamma_p)\)**
(difference 1.1e-16) and hence identical η surfaces. An observed split therefore
does *not* identify contrast as its cause. Read P4 as: splitting requires
\(\mathrm{d}\ln G/\mathrm{d}\ln\Gamma_p<-\tfrac12\,\mathrm{d}\ln R/\mathrm{d}\ln\Gamma_p\),
which either mechanism achieves and which η alone cannot tell apart. This is why
the experimental freeze's requirement to record rate, contrast and linewidth
*separately* is load-bearing rather than a convenience. Addendum A2 turns this
into a theorem.

## E2.6 What is unaffected

P1, P3's exact degeneracy, P4's algebra, P5's ordering, T4's inference procedure
(\(\Gamma_p^*\) recovered to 2.7 % from a 457 nm sweep at 2 % noise; 200/200
trials inside A1.3's factor-of-two bound), and all of E1's arithmetic reproduce
exactly. The frozen body is untouched.

---

# Addendum A2 — multiplicity ladder and gauge degeneracy

Date: 2026-08-27. Appended under the v3 immutability rule as a **theory-only**
extension, like A1. It introduces no measurement and no fitted value, and
changes nothing in the frozen body. Full derivations, numerical verification and
pre-registered tests: `docs/theory_a2_multiplicity.md`
(`code/theory_a2_multiplicity.py`, 15 regression tests).

A2 exists because E2 showed that A1 fails along a clean seam: **every broken
claim depends on the shape of the reconstructed kernel, and every surviving
claim follows from the structure of \(\eta=1/(G\sqrt R)\).** A2 freezes the
structural half and generalises it. No proposition below assumes anything about
the absorption spectrum; the Ho kernel appears only as a worked example.

**Theorem M (multiplicity ladder).** Under (M) with a unique interior maximum of
\(\Phi\) at \(\Gamma_p^*\), the optimal set at power \(I\) is the level set

\[
\Lambda_\eta(I)=\{\lambda:\ a(\lambda)=I_c/I\},
\qquad I_c=\Gamma_p^*/(\gamma A_{\max}).
\]

Its cardinality is piecewise constant in \(\ln I\) and changes only where the
level crosses a critical value of \(A\): **+2** at an interior maximum, **−2** at
a minimum, **−1** at a window edge. The transition powers are
\(I_k/I_c=A_{\max}/A_k\).

*Corollary (calibration-free).* The ratios \(I_k/I_j=a_j/a_k\) contain no
response parameter, and under (M) the \(a_k\) are read directly off the low-power
PL scan. The whole ladder follows from a relative spectrum plus the single scale
\(I_c\), which T4 already supplies. A1's doublet is the bottom rung: a unimodal
kernel has no interior critical values and admits only \(2\to1\).

**Theorem G (gauge degeneracy).** η depends on \((C,R,\Delta\nu)\) only through
\(\Phi=C^2R/\Delta\nu^2\). For power-law responses \(R=\Gamma_p(1+\rho)^{-s}\),
\(C=(1+\rho)^{-c}\), \(\Delta\nu=(1+\rho)^{w}\),

\[
\Phi=\Gamma_p(1+\rho)^{-E},\qquad E\equiv 2c+s+2w,
\qquad \rho^*=\frac{1}{E-1}\ \ (E>1).
\]

The response enters η through the single scalar \(E\), and splitting is exactly
\(E>1\). A1's four-row mechanism table is the sign of \(E-1\); its "all three"
case \(E=4\) gives \(\rho^*=1/3\), matching A1's numerical 0.3333.

*Corollary.* Responses sharing \(E\) share the entire η surface, \(\Gamma_p^*\)
and the ladder, so mechanism attribution from wavelength scans is impossible at
any number of powers — only \(C\) and \(\Delta\nu\) measured apart separate them.
*Conversely*, the ladder is invariant across the whole gauge plane, so it can be
predicted and tested **without knowing which mechanism operates**.

**Pre-registered tests.** T5 (ladder: predict \(I_k/I_c=1/a_k\) from the
low-power scan, verify the steps); T6 (universal degeneracy: all \(N\) members
share η exactly — an arbitrary-precision null test for (M), strengthening T3);
T7 (exponent closure: \(E\) from separately measured \(C,R,\Delta\nu\) must
predict the observed splitting power).

**Worked example, conditional on the reconstruction.** At 120 GPa the predicted
rungs fall at \(I/I_c=\) 1.4414 (ZPL, 514.46 nm), 1.5145, 1.5191, 1.8646, 3.6078
and 4.2989, giving \(N=2\to4\to6\to4\to3\to5\to3\); all six match their predicted
critical values to 0.00 %. **The lowest rung is the ZPL**, so a candidate-line
list without it cannot see the effect at all. The \(N=6\) plateau spans only
×1.003 in power and must not be claimed; the resolvable plateaux are ×1.05,
×1.23, ×1.93 and ×1.19. These *values* are conditional on the Fig. 1(e)
reconstruction — the *structure* is what A2 freezes.

---

# Addendum A3 — pressure-driven branch exchange of the optimum

Date: 2026-08-27. Appended under the v3 immutability rule as a **theory-only**
extension, like A1 and A2. It introduces no measurement and no fitted value.
Full derivations, numerical verification and pre-registered tests:
`docs/theory_a3_branch_exchange.md` (`code/theory_a3_branch_exchange.py`,
12 regression tests).

A2 asked what *power* does to the optimum at one pressure. A3 asks the
orthogonal question, and finds something A2's ladder cannot produce. An
absorption spectrum built on a Franck–Condon progression carries two
structurally distinct branches — the zero-phonon line and the phonon sideband —
whose peak heights scale differently with the Huang–Rhys factor \(S\). Pressure
raises \(S\), the branches move relative to one another, and where their
weighted peaks cross the global optimum **switches branch discontinuously**.

**Theorem X (branch exchange).** With \(A=\lambda\sigma_{\rm abs}\) the
fixed-power figure of merit and \(D(P)=\ln(A_{\rm SB}/A_{\rm ZPL})\), the global
optimum exchanges branch wherever \(D\) changes sign, and the optimal wavelength
jumps from \(\lambda_{\rm ZPL}(P^*)\) to \(\lambda_{\rm SB}(P^*)\) taking no
intermediate value. For a Franck–Condon pair
\(A_{\rm ZPL}\sim e^{-S}/\Gamma_{\rm ZPL}\) and
\(A_{\rm SB}\sim(1-e^{-S})/\Gamma_{\rm SB}\), so

\[
\frac{\mathrm{d}D}{\mathrm{d}P}
=\frac{\mathrm{d}S}{\mathrm{d}P}
+\frac{\mathrm{d}}{\mathrm{d}P}\ln\frac{\Gamma_{\rm ZPL}}{\Gamma_{\rm SB}}+\cdots>0
\]

whenever pressure strengthens the electron–phonon coupling and broadens the ZPL
at least as fast as the sideband. \(D\) is then monotone and **the crossing is
unique**: a generic, once-only exchange rather than an accident of one material.
Both antecedents are generic for a colour centre under compression, so what
varies between materials is only where \(P^*\) falls inside the accessible
range — not whether an exchange occurs. That is the general condition for
high-pressure optical sensing.

**Two degeneracies, not one.** A2's ladder is driven by power and needs
\(I>I_c\); A3's pair is driven by pressure and exists at **zero power**. At
\(P^*\) the low-power optimum is already twofold, before any power-induced
structure appears — which is what makes the two separable in practice. Together
they give the full \((P,I)\) phase diagram: A3 selects the branch, A2 splits it.

**Verification, in the raw samples rather than the interpolation.** Because E2
established that this kernel's pressure interpolation is unreliable, A3
identifies both branches as local maxima of the *extracted samples themselves*.
Both are present at all seven published pressures. The ZPL branch at 0 GPa lands
at **637.1 nm (1.945 eV)** — the NV⁻ ZPL, a value fitted nowhere in the
pipeline and therefore an independent check that the extraction is anchored
correctly. \(\ln(A_{\rm SB}/A_{\rm ZPL})\) rises monotonically at all seven
points, by ×4.33 in total, with exactly one sign change.

**Worked example, conditional on the reconstruction.** At fixed incident optical
power \(P^*=87.9\) GPa, where the ZPL branch sits at 534.8 nm and the sideband at
463.5 nm: the optimum jumps **71.3 nm**. At fixed incident photon flux
\(P^*=75.6\) GPa instead, so the power convention must be stated whenever
\(P^*\) is quoted. The driver is visible directly: the Franck–Condon
displacement \(S\hbar\omega\) grows from 0.232 to 0.404 eV (+74 %, monotone,
1.27 meV/GPa), the ZPL cross-section falls by ×6.04 over 0–120 GPa while the
sideband falls by only ×1.34, and the two branches carry different pressure
coefficients (3.83 against 5.11 meV/GPa).

**What this is not.** The 457/532 fixed-wavelength ranking crossover sits at
51.4 GPa in the same kernel, 36 GPa away: that crossover is the geometry of one
unimodal peak sweeping past the midpoint of two probes and needs no branch
structure at all. The v1 Franck–Condon envelope has exactly one local maximum at
every pressure and so cannot show this effect; its "green/blue crossover
≈ 86 GPa" is unrelated, and its numerical closeness to 87.9 GPa is a coincidence
(`docs/novelty_and_exponent_audit.md` §4).

**Pre-registered tests.** T8 (track \(\lambda_{\rm PL}(P)\) at low power through
\(P^*\); a jump confirms Theorem X), T9 (check both branches coexist either side
of \(P^*\), distinguishing exchange from branch disappearance), T10 (compare the
two wavelengths at low power, separating A3's zero-power degeneracy from A2's
ladder), T11 (measure the branch separation \(S\hbar\omega\) against pressure).
T11 is the cheapest and most decisive: it needs neither absolute nor intensity
calibration, only two peak positions per pressure, and if the separation does
not grow the antecedent fails and no exchange is predicted.

---

# Erratum E3 — the panel (e) zero-phonon-line peaks are clipped

Date: 2026-08-27. Filed under the v3 bug-fix/erratum rule. It corrects
Addendum A3 only. **The frozen body and Addendum A2's structural propositions
are untouched.**

The extraction in `code/data/ho_fig1e_absorption.csv` was checked pixel by pixel
against the source Fig. 1. Working: `code/figure_validation.py`, 12 regression
tests.

## E3.1 The sideband branch is exact

Tracing each coloured curve in panel (e) and taking its maximum on the
absorption side reproduces the extracted CSV to **better than 1 % in height at
all seven pressures**, and to 0.015 eV in position at six of them (20 GPa picks
a shoulder 0.064 eV away). Everything A3 says about the sideband stands, as does
its 0 GPa anchor at 637.1 nm.

## E3.2 The ZPL peak heights are an artefact

Every ZPL spike in panel (e) rises from its own baseline to within one unit of
the top of the axis. Measured tops, on an axis ending near 15:

| P [GPa] | 0 | 20 | 40 | 60 | 80 | 100 | 120 |
|---|---:|---:|---:|---:|---:|---:|---:|
| spike top | 14.93 | 14.95 | 14.67 | 14.03 | 14.09 | 14.62 | 14.94 |

The spread across all seven is 0.92 units. They are drawn as near-delta lines
and **the plot simply cuts them off**. The CSV heights (6.364 falling to 1.054)
are a digitisation artefact of a clipped feature.

**Withdrawn:** the "ZPL cross-section falls by ×6.04 over 0–120 GPa" claim, and
with it \(P^*=87.9\) GPa (fixed optical power), 75.6 GPa (fixed photon flux),
and the 71.3 nm jump. `theory_a3_branch_exchange.exchange_pressure` is retained
as a faithful computation of what the CSV contains, and marked superseded.

## E3.3 Theorem X survives, on better evidence

Panels (b) and (c) publish precisely what is needed, as theory curves with
markers, calibrated from their own axis ticks and not clipped
(`code/data/ho_fig1_panels_bc.csv`):

- **\(S_{\rm abs}\) rises 3.023 → 4.554** over 0–120 GPa, +51 %, monotone,
  \(\mathrm{d}S/\mathrm{d}P = 12.7\) milli/GPa. Theorem X's **first driver is
  confirmed directly from the published theory**, not inferred.
- **DWF\(_{\rm abs}\) falls 0.0205 → 0.00226**, a factor **9.07**, monotone.
  The ZPL weight collapse is real — steeper than the clipped spikes suggested.
  (A single effective mode, \(e^{-S}\), would give only ×4.62, so the published
  DWF is genuinely multi-mode.)

## E3.4 What changes: P\(^*\) is bandwidth dependent

The clipped-peak version silently compared two peak heights. It should not have:
the sideband enters through a **lineshape density** (per eV) while the ZPL
enters through a **dimensionless weight**, which must be divided by whatever
bandwidth samples it — the laser linewidth, or the ZPL's own width, whichever is
larger. Writing

\[
r(P)=\frac{\lambda_{\rm SB}\,\sigma_{\rm SB}}{\lambda_{\rm ZPL}\,{\rm DWF}_{\rm abs}}
\quad[\mathrm{eV^{-1}}],
\qquad
\frac{A_{\rm SB}}{A_{\rm ZPL}} = r(P)\,W ,
\]

the exchange occurs where \(r(P)=1/W\). **\(r\) is monotone increasing, ×6.51
over 0–120 GPa, so the crossing is still unique** and Theorem X is intact. But
\(P^*\) now moves with \(W\):

| \(W\) [meV] | 1.5 | 2 | 3 | 5 | 7 | 9 |
|---|---:|---:|---:|---:|---:|---:|
| \(P^*\) [GPa] | 119.6 | 104.0 | 77.7 | 45.2 | 24.1 | 6.7 |

Outside roughly 1.5–9.7 meV there is no crossing inside 0–120 GPa. Equivalently,
the ZPL only competes at all for an excitation bandwidth below
\(1/r(P)\) — 9.7 meV at ambient falling to 1.5 meV at 120 GPa.

**This is a strengthening, not merely a retraction.** The bandwidth dependence
is a new, sharply testable prediction: sweeping the laser linewidth moves the
exchange pressure monotonically, and no other mechanism in the framework does
that. It also makes T11 (measure \(S\hbar\omega\) against pressure) more
important, since it tests the antecedent without touching the bandwidth
question at all.

**Reporting rule.** Never quote \(P^*\) without stating the excitation
bandwidth, alongside the existing rule about the fixed-power versus fixed-flux
convention.

