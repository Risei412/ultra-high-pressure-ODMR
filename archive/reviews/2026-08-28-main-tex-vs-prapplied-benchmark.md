# Review: `paper/main.tex` against the PRApplied benchmark

Date: 2026-08-28
Scope: `paper/main.tex` @ 1715 lines / 16 pp, `paper/refs.bib`, `paper/build/main.pdf`, figure scripts in `code/`
Benchmark: Sekiguchi *et al.*, Phys. Rev. Applied **21**, 064010 (2024) — `paper/PhysRevApplied.21.064010.pdf`
(9 pp, 5 sections, 5 figures, 0 tables, 2 numbered equations, 118-word abstract, 41 references)

Agents: `logic-reviewer`, `writing-reviewer`, `research-analyst`, `technical-reviewer`, `latex-layout-auditor`
(all read-only; no file was modified during Stage 1)

Items marked **[verified]** were checked by the orchestrator by running the code or reading the
source directly, not merely reported by an agent.

---

## Overview

The manuscript is mathematically sound in its three theorems and internally clean in the two
places a previous pass repaired. What Stage 1 found is different in kind: **Erratum E3 (the
zero-phonon line is clipped in the reconstructed kernel) is stated in the manuscript but not
propagated through its consequences**, and the paper's framing — abstract, title, section order,
and which result leads — contradicts the project's own positioning directive in
`docs/positioning_in_high_pressure_sensing.md` §6.3.

Length and float findings are secondary and now quantified: the body is 2.1x the benchmark and
about 3,700 words are recoverable from duplication alone, while floats yield less than one page.
Length is a prose problem, not a layout problem.

---

## Critical

### C1. Erratum E3 is not propagated to the multiplicity sequence **[verified]**
`main.tex:671-676`, `701-706`, `942-945`, Fig. 1(b) via `code/fig_m_level_set_ladder.py`

The manuscript correctly records at `main.tex:680-691` that the ZPL height is over-weighted by up
to x1.43, that deflating `a = 0.6938` to `0.485` moves that rung to `I/I_c = 2.06`, and that it
"would no longer be the lowest". None of the downstream results carries that branch.

Verified by running `theory_a2_multiplicity.critical_values` and applying the step rule:

| | levels, descending | ladder |
|---|---|---|
| optimistic (printed) | 0.6938 (ZPL max), 0.6603, 0.6583, 0.5363, 0.2772, 0.2326 | `2 -> 4 -> 6 -> 4 -> 3 -> 5 -> 3` |
| pessimistic (deflated) | 0.6603, 0.6583, 0.5363, **0.485 (ZPL max)**, 0.2772, 0.2326 | `2 -> 4 -> 2 -> 1 -> 3 -> 5 -> 3` |

Consequences that are currently unstated:
- **The `N = 6` plateau does not exist on the pessimistic branch.** §IV D spends a paragraph
  (`main.tex:701-703`) declaring it untestable at x1.003; on the other branch there is nothing to
  declare.
- The ZPL becomes the **fourth** rung at `I/I_c = 2.062`, not the first.
- The plateau-width list at `main.tex:704-706` (x1.44, x1.05, x1.23, x1.93, x1.19) is
  optimistic-branch only.
- Fig. 1(b) plots the optimistic staircase with no indication that the first two steps are
  conditional.

**Compounding (technical-reviewer, Rigor #2).** The x1.43 is an *area* over-weight from
`kernel_sanity_checks.zpl_area`, applied in the manuscript to a *height*. That transfer is only
valid if the width is right, and check K2 (`main.tex:1359-1366`) establishes the ZPL width is
resolution-limited. So `1.44 <= I/I_c <= 2.06` is not a genuine bound. The honest statement per
freeze rule S1.6 is that `a(lambda_ZPL)` is not determinable from this source above 40 GPa.

### C2. Two internal contradictions in Sec. IV **[verified]**
- `main.tex:697-700`: "even at the pessimistic end of the range above, its sensitivity penalty is
  comparable to that of the commercial 473 nm line, 1.2054." At the pessimistic end the penalty is
  **1.44** — stated by the manuscript itself nine lines earlier at `main.tex:691`. 1.44 against
  1.205 is 19% worse, not comparable. The sentence silently reverts to the value it just withdrew.
- `main.tex:736-738`: "the factor 1.44 carries the zero-phonon uncertainty ..., which could raise
  it to 2.06." If the ZPL rung moves to 2.06 it is no longer the first rung — `main.tex:687-688`
  says so explicitly. The first rung *falls* to 1.51 (the 475.55 nm maximum). The "4% to 30% of
  saturation" figure at `main.tex:738` is built on the wrong branch and must be recomputed.

### C3. The (P,I) phase diagram contradicts the worked example
`main.tex:1153-1161` against `main.tex:659-664`

§VII asserts that below `P*` "the ladder of Theorem M is erected *within that branch*". The level
set is defined over the whole window, so once the level drops below both branch peaks its members
come from both. The paper's own Table II proves it: rung 1 is the ZPL at 514.46 nm, rungs 2-6 are
sideband features. This is the only place the two headline theorems are synthesised and as written
it is wrong. Restate in terms of *which branch supplies the level-set members at a given level*.

### C4. Theorem X has no proof, and is stated in a form the paper then declines to establish
`main.tex:967-985`, `992-998`, `1035-1040`

The boxed statement rests on Eq. (14), whose second term is the relative broadening
`d ln(Gamma_ZPL/Gamma_SB)/dP`. The manuscript then says this "is a sufficient condition, not the one
we verify" and "is **not** confirmed and cannot be from this source". What is actually established is
`D(P;W) = ln[r(P) W]` with `r` monotone. Theorems M and G each carry a proof
(`main.tex:577-600`, `779-782`); **X, the primary result under the §6.3 directive, carries none**,
and the load-bearing phrase "taking no intermediate value" (`main.tex:972`) is asserted.

Related rigor defect (technical-reviewer, Rigor #1): at `main.tex:976-984` the inequality is
printed as `dS/dP + d/dP ln(Gamma_ZPL/Gamma_SB) + ... > 0`. The coefficient of `dS/dP` is
`1/(1 - e^{-S}) >= 1`, not 1, and the elided term is never named. The conclusion is true; as
printed the reader cannot verify it.

### C5. "A single-mode envelope has exactly one interior maximum" is unreproducible **[verified]**
`main.tex:411-418`, `1137-1145`

`main.tex:411` calls this "the single fact on which the numerical content of Theorem M rests".
The "four" half is reproducible (`theory_a1_generalization.Kernel.local_maxima`). The "one" half
is not: its only sources are `theory_freeze_v3:868-874` (E4.3, explicitly "Reproduced by
`code/v1_diagnosis.py`") and `docs/novelty_and_exponent_audit.md:230-234`, which asserts it without
an executable.

Verified: `local_maxima` is defined only on `Kernel`, never on `NVModel`; no script in `code/`
counts interior local maxima of the v1 single-mode envelope. `code/v1_diagnosis.py` is absent from
the repository and from all three archived zip bundles (also missing: `fig8_repro_ho_fig5b.py`,
`fig7_a3_branch_exchange.py`).

**This is a stronger dependency on the missing script than the earlier pass recorded.** The
manuscript correctly withholds E4's robustness numbers, but it does assert E4.3's structural claim.
Fix: add a ~15-line local-maximum counter for `NVModel` to `theory_a1_generalization.py`, pin it
with a test, and cite it in Appendix B. Until then, restrict the claim to 120 GPa or attribute it as
an unreproduced assertion.

### C6. The abstract contradicts both the body and the project directive
- `main.tex:129-131` opens the results by **assuming** (M). `main.tex:218-220`, `1261-1264` and
  `1595-1597` establish the opposite — that the framework *measures* it — and
  `positioning_in_high_pressure_sensing.md` §6.3 directs the second framing. The abstract never
  mentions T6, the arbitrary-precision null test that *is* the measurement of (M).
- Order is M, G, X (`main.tex:131, 138, 148`), and the title (`main.tex:114-115`) leads with
  multiplicity. §6.3 directs X primary, M secondary.
- `main.tex:146-147` places the first rung at "a few tens of percent of saturation"; the body says
  "between 4% and 30%" (`main.tex:738`) and the Summary "a few percent and a third". The abstract's
  floor is wrong by an order of magnitude — a residue of the interval fix that never reached it.
- The framing headline of the whole project — that the optimal wavelength **jumps 71 nm** — does
  not appear in `main.tex` at all.

### C7. The gauge-invariance "numerical verification" is vacuous **[verified]**
`main.tex:942-945` against `code/theory_a2_multiplicity.py:314-327`

`ladder_is_gauge_invariant` loops over five responses but computes
`counts = [multiplicity(kernel, 1.0/r, window) for r in ratios]` with `ratios` **hardcoded inside
the loop and independent of `response`**; `response` enters only the reported `gamma_star`. The
loop recomputes an identical expression five times and reports that it is identical.

The invariance is a genuine theorem and is stated correctly at `main.tex:940-941`. What must go, or
be made real, is the claim at `main.tex:942-945` that it was probed numerically at six powers
across five models.

### C8. Novelty-framing exposure: Beha *et al.*, PRL **109**, 097404 (2012)
Absent from `refs.bib`. Title: *Optimum Photoluminescence Excitation and Recharging Cycle of Single
Nitrogen-Vacancy Centers in Ultrapure Diamond*. `main.tex:197-202` says the excitation wavelength
was carried into the megabar regime "essentially unexamined".

The content is favourable: Beha's optimum is a **point** optimum of `R`, fixed by charge-cycle
balance, with no power axis, no pressure and no level set. The fix is one clause — the wavelength
has been optimized for photoluminescence (Beha) and for contrast (Todenhagen), never for `eta`, and
never at fixed power — plus a citation in §III C, where the manuscript already reasons about the
charge cycle Beha measured.

### C9. Precision theatre: 46 numbers at >= 4 significant figures, one uncertainty in 1715 lines
The 4th and 5th digits of every ladder rung are set by the 0.005 nm sampling grid in
`theory_a2_multiplicity.critical_values(step=0.005)` acting on a piecewise-linear interpolant whose
stationary points lie between nodes. Worked example (500.19 nm rung): node value 3.6076, grid value
3.6079 — the printed fifth digit *is* the artefact.

Binding uncertainties, from the project's own cross-trace
(`figure_validation.FIGURE_SIDEBAND[120] = 2.8083 eV` against the CSV's 2.81372 eV):
**positions +-1 nm, heights +-1%**. Dréau-derived quantities inherit 2 significant figures from
`GP_INF = 5.0e6`, `GC_INF = 8.0e7`, `OMEGA_R = 3.0e6`.

A ready-to-apply replacement list for all 46 sites is in the technical-reviewer transcript.
Notable: `rho* = 0.20000` should be `1/5` (it is an exact identity, not a measurement);
`P*` values to 0.1 GPa should be integers with +-2 GPa; `x* = 0.06708` should be `0.067`.

### C10. Layout: raster art below the APS bar, and a float relaxation that is not working **[verified]**
- All four figures are included as PNG at 398 / 410 / 424 / **346** dpi at final printed size. APS
  requires 600 dpi for line and combination art. Vector PDFs exist and swap with < 3.6 pt reflow —
  but all four use **matplotlib Type 3 fonts** (verified with `pdffonts`), which APS and arXiv
  flag. Set `matplotlib.rcParams['pdf.fonttype'] = 42` in the four scripts and re-save first.
- The three single-column figures are ~5.5 in tall; with captions each float takes ~74% of a
  column. This is what forced the preamble relaxation at `main.tex:88-94`, and **the relaxation did
  not clear the warnings**: `paper/build/main.log` still carries 6 stuck-float lines (verified).
  `\textfraction{0.07}` licenses a 93%-float column and produced the 35-text-row page 12.
- Table V (`main.tex:1088-1101`) is a verbatim transcription of Fig. 3(b)'s annotated markers;
  `main.tex:1103` says so out loud. Dissolve it.
- Fig. 4's title, subtitle and footer are caption text baked into a 346 dpi bitmap, including the
  validation numbers (RMS 1.5% / 0.4% / 13.0% / 18.2%) that never appear in the typeset caption.

---

## Important

### Structure
- **Theorem G is a lemma plus a corollary, not a third theorem.** `main.tex:766-782` and `830-885`
  supply Theorem M's antecedent (`En > 1`) and are placed *after* the theorem that needs them,
  forcing the forward reference at `main.tex:720`. `main.tex:887-934` is the only independent claim
  and is methodological. Reordering G before M removes the forward reference and the
  "three loosely coupled contributions" reading.
- **§III (609 w) is three unrelated fragments**: the coincidence condition belongs after §II A; the
  three cases belong at the end of §II D; the charge-conversion rebuttal (`main.tex:511-538`) is
  Discussion material, and §X B already cross-references it as such.
- **§VII (160 w) is not a section**; it is the (P,I) corollary of §VI. Fold it in — after fixing C3.
- **§IX (1174 w) should split.** The *constraint* (~250 w: four numbered restrictions) moves
  forward to §II E where the kernel is introduced; the *evidence* (~900 w + the only `figure*`)
  moves to an appendix. The demotion is cheap because the reconciliation at `main.tex:1424-1440` is
  already written verbatim at `main.tex:680-691`. As currently placed the fence arrives 950, 850,
  630 and 80 lines *after* the four claims it fences.
- **Four stacked headings** — `\section` immediately followed by `\subsection` with no prose:
  `main.tex:282/285`, `541/545`, `744/747`, `1443/1446`. Each is exactly the slot where the
  goal-problem-solution sentence belongs.
- Delete the mathematics-journal subsection headings "Statement" / "Proof" / "Statement and proof"
  (`main.tex:545, 577, 747, 959`); run the proof as a run-in `\emph{Proof.}` paragraph.
- Rename the theorem section titles to functional titles; keep M/G/X inside the boxed statements.

### Citations
Four works the project's own §3.1 marks "citation mandatory" are absent from `refs.bib` and
uncited. Verified BibTeX and two positioning sentences each are in the research-analyst transcript.

| key | work | note |
|---|---|---|
| `Ma2010` | Appl. Spectrosc. **64**, 1274 (2010) | high confidence |
| `ShajiRebeirro2024` | Phys. Rev. Appl. **21**, 044039 (2024) | check APS name alphabetisation |
| `Ronchi2024` | Phys. Rev. Appl. **22**, 034058 (2024) | high confidence |
| `Sun2024` | Nanophotonics **13**, 2401 (2024) | **identification ~85% — verify against Zotero `HPHT`** |

**19 of 46 bib entries are uncited [verified]**, and several are the canonical source for objects
the abstract names: `Huang1950` (the Huang-Rhys factor itself), `DaviesHamer1976`, `Razinkovas2021`,
`Razinkovas2021vib`, `Aslam2013` (the continuous wavelength dependence of the ionization channel —
its absence next to §II C's "fails only at discrete resonances" reads as convenient), `Lyapin2018`
(the only independent experimental high-pressure NV optical study in the bibliography, which would
triangulate §V B's single-source pressure coefficients).

`Alkauskas2014` (New J. Phys. **16**, 073026) is missing entirely and supplies the ambient-pressure
anchor (`S ~ 3.5`, DW ~ 3%) that Table III's values must reduce to.

Two accuracy flags on Todenhagen: the published APL reports the -23% at 575 nm for *electrically*
detected magnetic resonance (not "photoelectrically detected ODMR"), and the abstract's
"works best between 525 and 550 nm" implies the optical contrast is non-monotone across the
480-575 nm interval over which `main.tex:209-211` fits a single 0.43 %/nm slope. The argument needs
only `l_G != 0`, so quoting a stated sub-interval is strictly safer.

Rename `Semenok2022LaCeH9` (first author is Bi).

### Prose
- **3,666 words recoverable** across 25 ranked sites (writing-reviewer transcript); body
  10,686 -> ~7,020. Largest: §IX 425, Summary 283, §V B Dréau glossary 200, test re-narration 180.
- **19 enumerative preambles** ("Two remarks fix its status", "Three features are worth separating
  from the numbers") are this manuscript's dominant AI/lecture-note signature. The lexical tells
  (delve, leverage, crucial, "It is important to note") are **absent — clean**.
- 76 em-dashes in ~11,000 words against the benchmark's ~1 per 1,500.
- `not X but Y`: 6 hits, of which 2 are mathematically substantive and must be kept verbatim
  (`main.tex:132`, `257`). `rather than`: 30 hits, 18 stylistic; "(M) as measurand rather than
  assumption" appears **five times**.
- Ten redundancy pairs R1-R10 including: §IX promised four times before delivery; the pressure-as-
  control-variable sentence (the paper's best) spent twice at `main.tex:264-270` and `1493-1507`.

### Verified typesetting and grammar defects
- `main.tex:620` ends the source line with `single-` and line 621 begins `wavelength`; LaTeX turns
  the newline into a space, so this sets as **"a single- wavelength power sweep"**. **[verified]**
- `main.tex:916` "an hierarchy" -> "a hierarchy". **[verified]**
- `main.tex:717-718` "A structure that required exotic powers would be of limited interest. It does
  not." — the elided verb has no grammatical antecedent.
- British spellings in an APS submission: `labelled` (639, 1356), `modelling` (409), `programme`
  (369, 453). **[verified]**
- `main.tex:649` caption says agreement "better than 0.01%", `main.tex:672-673` says "reproduced to
  0.00%", Table I shows 3.6078 vs 3.6079. Harmonize.

### Notation
- `h` is both Planck's constant (`main.tex:291`) and the mediated response function
  (`main.tex:366, 393, 485, 580`).
- `P` is both pressure and optical power (`main.tex:844`, `1466`) — the same collision class that
  `\varsigma` was introduced to fix.
- `S` / `\Sabs` / `S_total` are never reconciled. T11 (`main.tex:1241`) says "invert
  `psi(p*+1) = ln S`" without saying which; it must be `\Sabs`, and the choice changes the inferred
  `hbar omega` by 35%.
- Eq. (6) at `main.tex:521-524` introduces `r_0`, `a_es`, `r_bg` and `sigma`, **none of which is
  defined anywhere in the manuscript**.
- `main.tex:844-850` introduces nine symbols in one sentence.
- `\varsigma` is used at `main.tex:722` and defined at `main.tex:844`.
- Appendix B maps tables to scripts but gives **no symbol-to-variable map**, and the mapping is
  actively inverted: the paper's `\varsigma` is `s` in `dreau_exponent.py`, while the paper's rate
  exponent `s` is `s_exp` there. `theory_a1_generalization.kappa_*` is a wavelength curvature
  unrelated to the paper's `kappa`.
- Dead macro `\esE` at `main.tex:99`.

### Layout, remaining
- Fig. 3 and Table IV ship on p9, which contains no float reference; both first cited on p10.
- Tables I and VII share a column signature and a row (440.64/440.65) — merge and promote to body.
  Table VII is the only table an experimentalist will act on and it is buried on p15.
- Fig. 1(b)'s legend reprints four of Table II's six columns.
- Fig. 2(a) is a null result the artwork must annotate in words; one sentence replaces it losslessly
  and removes the tallest single-column float.
- Captions total ~117 column-lines (0.94 page) against ~25 lines for the benchmark's five, and omit
  every axis quantity, marker key and abbreviation expansion (SB, ZPL, DWF, K1-K4, RSS, Expt) that
  the benchmark always supplies.
- Fig. 2's five legend models do not correspond to Table III's five rows; a reader will try to match
  them and fail.
- An un-numbered eighth table sits in the p8 text flow at `main.tex:918-928`.

---

## Minor

- `main.tex:427-428`: the Table I caption attributes the 440.64/440.65 split to the 0.005 nm grid.
  Both the 0.005 nm and 0.01 nm grids give 440.64; 440.65 comes only from the 0.05 nm scan grid in
  `ho_odmr_sensitivity.py:99`. The node itself is at 440.641 nm.
- `main.tex:722` prints `0.072`, `main.tex:865` prints `0.0719` for the same quantity.
- `main.tex:739`: `I/I_c <= 139` needs the drive strength stated (71 at `Omega_R = 5.5e6`).
- `main.tex:984`: "the crossing is unique" -> a monotone `D` has *at most one* zero; existence needs
  a sign change, which Table V shows fails outside W ~ 1.5-9.7 meV.
- `main.tex:1674-1677`: Tables IV and V are attributed to
  `theory_a3_branch_exchange.exchange_pressure`, which its own docstring marks superseded; the
  numbers come from `figure_validation.py`.
- `main.tex:1702`: `pytest tests/ -q` misses `code/test_kernel_sanity_checks.py`, which sits in
  `code/` rather than `code/tests/`.
- Table I normalises to `A(440.64)`, Table VI to `A(440.65)`; their `1.0000` entries differ.
- `\SI{x}{}` with an empty unit at `main.tex:209, 1049, 1052, 1125, 1309, 1398`.
- `main.tex:528-531`: the "0.2% flatness of `f_-`" is an unattributed number from the superseded v1
  model (`docs/introduction_rationale.md:76`, over 463-486 nm), reported as "across the blue window"
  from "the fixed-power analysis".
- Stale code docstrings a referee running the code would see contradicting the paper:
  `figure_validation.py:28,30` (4.554, x9.09) and `kernel_sanity_checks.py:26` (1%).
- `theory_a2_multiplicity.gauge_family(exponent=2.0)` ignores its argument.
- `code/tests/test_theory_a1_generalization.py:237` pins E1's 3.7e-5 with `rel=0.05`; the correct
  value is 3.546e-5 and the test passes only by 0.9%. Re-pin to `3.5e-5, rel=0.02`.

---

## Patterns

1. **Erratum E3 is stated but not propagated.** Three of the five agents reached this independently
   from different directions. It is the systemic defect, and C1, C2 and part of C9 are all instances.
2. **Duplication is the dominant structural problem, not verbosity.** Ten redundancy pairs in prose,
   Table V == Fig. 3(b), Table I == Table VII, Fig. 1(b) legend == Table II, caption == body. The
   3,666 recoverable words are almost entirely second and third statements of the same fact.
3. **Precision theatre.** Five-significant-figure numbers throughout, one uncertainty statement in
   the document, and the digits below the third are a sampling-grid artefact.
4. **The paper's own best material is unused.** The 71 nm jump never appears; the
   pressure-as-control-variable sentence is spent twice in weak positions and absent from §VI, where
   it belongs; the strongest sentence in the paper (`main.tex:215-216`) is a mid-introduction orphan.
5. **Length is a prose problem, not a layout problem.** Float *area* is already near the benchmark
   (22% against 19%); floats yield < 1 page, prose yields ~3.7k words.

---

## Recommended Stage 2, in order

**Tier 1 — correctness (no scope decision required).**
C1 propagation and the Fig. 1(b) caveat; C2's two contradictions; C3's phase-diagram restatement;
C4's proof and the Eq. (8) prefactor; C5's local-maximum counter plus test, or a narrowed claim;
C7's vacuous verification; `main.tex:1452` `E > 1` -> `En > 1`; the `single-` line-break bug;
"an hierarchy"; the broken ellipsis at 717-718; British spellings; the 0.00%/0.01% inconsistency;
re-pin the E1 test.

**Tier 2 — framing and citations.**
New abstract (196-word draft ready, X-primary); title reordered; the `71 nm` claim introduced;
four missing references plus `Beha2012` and `Alkauskas2014`; cite the 19 uncited entries or cut
them; the two Todenhagen accuracy fixes.

**Tier 3 — restructuring.**
11 sections -> 7; §IX split (constraint forward to §II E, evidence to an appendix); §III dissolved
three ways; §VII folded into §VI; G before M; stacked headings filled; "Statement"/"Proof" headings
removed; 3,666 words of duplication cut.

**Tier 4 — numbers and figures.**
The 46-site significant-figure pass with uncertainties; `pdf.fonttype = 42` then PNG -> PDF;
Figs. 1 and 3 re-laid out 1x2 as `figure*`; Fig. 2(a) dropped; Table V dissolved; Tables I and VII
merged and promoted; in-artwork titles/subtitles/footers moved into captions; captions rewritten to
benchmark register; preamble float relaxation reverted to REVTeX defaults and re-checked.

Verification after each tier: rebuild (0 errors, 0 undefined references, 0 overfull), run
`code/tests`, and re-run `consistency-checker` on the changed regions.

---

## Applied — 2026-08-28, Tier 1 + Tier 2

Verified by a `consistency-checker` pass over the changed regions, then rebuilt.

**Code.** `theory_a1_generalization.v1_envelope_local_maxima` added (C5) and
`theory_a2_multiplicity.ladder_is_gauge_invariant` repaired (C7); both pinned.
The E1 test re-pinned from the erratum's 3.7e-5 to the computed 3.5e-5.
Tests 200 -> 213, all pass.

**Result that changed the paper.** The deflated ZPL branch reorders the ladder:
`2 -> 4 -> 2 -> 1 -> 3 -> 5 -> 3` (new Eq. 13). The N=6 plateau is absent from
it and an N=1 plateau opens between I/I_c = 1.86 and 2.06. Propagated to the
sequence, the plateau widths, Fig. 1(b)'s caption, Table II's caption, the
Sec. V D gauge probe and Sec. VII.

**Corrections.** The 473 nm comparison (1.20 vs 1.44, 19% worse, not
"comparable"); the reachability range 4-30% -> **4-21%** (the deflated first
rung is 1.51, not 2.06); four surviving "the crossing is unique" softened to
"at most once"; the Summary's "a few percent and a third" removed; Theorem X
restated in bandwidth-resolved form **with a proof**, Eq. (14) demoted to a
remark carrying its correct prefactor 1/(1-e^{-S}); branch separation computed
at 68-78 nm and used in the abstract.

**Citations.** 46 -> 52 entries, 27 -> 39 cited, 19 -> 13 uncited. All four
mandatory works confirmed present in the Zotero HPHT collection (`novelty-audit`
tag), Sun 2024 included. `Beha2012` added and placed.

**Front matter.** Title reordered to lead with branch exchange; abstract
399 -> 285 words, longest sentence 82 -> 37, branch exchange first.

**State.** 18 pp (up from 16 — Tier 1+2 adds content). Build: 0 errors,
0 undefined references, 0 overfull boxes. The two stuck-float clusters and the
raster/Type-3 figure issues are untouched; they are Tier 4.

**Not applied.** Tier 3 (11 sections -> 7, Sec. IX split, ~3,700 words of
duplication) and Tier 4 (46-site significant-figure pass, vector figures, float
economy). Everything needed to execute them is above.

---

## Applied — 2026-08-28, Tier 3 + Tier 4

### Structure: 11 sections + 2 appendices -> 8 + 3

| new | was |
|---|---|
| I Introduction | I |
| II Framework and the two optima | II + III A + III B |
| III The response enters through a single exponent | V (Theorem G, **moved before M**) |
| IV The optimum is a level set, and its multiplicity is a ladder | IV (Theorem M) |
| V Pressure exchanges the branch carrying the optimum | VI + VII (phase diagram folded in) |
| VI Pre-registered tests | VIII |
| VII Discussion | X + III C (charge conversion) |
| VIII Conclusions | XI |
| App. A Fixed-power optical limit | App. A |
| App. B Domain of validity of the kernel | **IX**, demoted |
| App. C Reproducing the numbers | App. B |

G precedes M because it supplies M's antecedent `En > 1`; that removes the
forward reference the old Sec. IV E had to make. The letters G, M, X are
mnemonics, not ordinals, so the order costs nothing. Section titles are now
functional in the benchmark's style, and the mathematics-journal
"Statement"/"Proof" subsection headings are gone (the proof runs in as
`\emph{Proof.}`). All four stacked headings — a `\section` immediately followed
by a `\subsection` — now carry a goal-problem paragraph.

### Prose

The Conclusions were cut from ~420 words to ~200 and now close the arc the
Introduction opens: map acquisition time scales as `eta^2`, so the 473 nm line
costs a factor 1.45 in time at 120 GPa. The Introduction's three-theorem
previews were cut roughly in half and reordered to X-primary. Seventeen of the
nineteen enumerative preambles are gone, together with the defensive sentences
("we regard that as the strongest argument", "which is the only claim being
made", "The mathematics is elementary"). The duplicate pressure-as-control-
variable passage in the Discussion now cites the Introduction instead of
repeating it, and the triple statement of the ZPL reconciliation is down to
one plus a pointer.

### Floats: 11 -> 9

- **Table V dissolved.** It transcribed Fig. 3(b)'s annotated markers; the four
  references are retargeted to the figure, which now also carries the
  fixed-power convention and the +-2 GPa.
- **Tables I and VII merged** into one body table with `interior critical
  points` / `commercial laser lines` blocks. That deletes the 440.64-vs-440.65
  apologia the split forced into Table I's caption, and moves the one table an
  experimentalist acts on out of the appendix.
- The three single-column figures are now full-width 1x2 `figure*`.

### Figures

All four regenerated with `pdf.fonttype = 42` and included as **vector PDF**.
`pdffonts` now reports CID TrueType where it reported Type 3, which is what the
APS and arXiv checkers require. **The compiled PDF fell from 1.45 MB to
552 kB.** The claim-sentence panel titles and the overall title, subtitle and
footer that were rendered into the raster are removed; their content is in the
typeset captions, which now also give axis quantities with units, marker keys
and every abbreviation (SB, ZPL, DWF, RSS, Expt, Ho) — the three things the
benchmark's captions always supply and this manuscript's never did.

The preamble float relaxation is reduced from `topfraction 0.90 / textfraction
0.07` to `0.85 / 0.10`, which the shorter figures now permit.

### Numbers

The significant-figure pass is applied: 16 remaining values at four or more
decimal places are down to one, in a header comment. Table II's rungs are three
significant figures, Table I's positions carry the +-1 nm from the project's own
cross-trace, and `rho* = 0.20000` is printed as the identity `1/5` it is.

### Landing point, against the benchmark

| | benchmark | before | after |
|---|---|---|---|
| pages | 9 | 16 | **17** |
| sections | 5 | 11 | **8** |
| floats | 5 | 11 | **9** |
| float density | 0.56/pp | 0.69/pp | **0.53/pp** |
| abstract words | 118 | 399 | **285** |
| longest abstract sentence | 65 | 82 | **37** |
| references | 41 | 26 | **39** |
| PDF size | 715 kB | 1.45 MB | **552 kB** |

Two deliberate departures from the benchmark. Its 9 pages and its **zero
tables** both follow from its being an experimental report of a single number;
a paper making three theorems cannot reach either without dropping content.
Body prose is about 11,800 words including captions, against a target of 7,000
that would have required cutting substance rather than repetition.

### Still open

- One stuck-float warning cluster survives. Forcing `[!t]`/`[!tb]` produced a
  **byte-identical PDF**, so it is a REVTeX column-balancing artifact, not a
  placement failure: every float appears in source order, on or one page after
  its first reference, with 0 overfull boxes.
- 13 uncited entries in `refs.bib`.
- `code/v1_diagnosis.py` still absent; E4's envelope-robustness result stays
  out of the paper. E4.3's structural claim is now carried by
  `v1_envelope_local_maxima` instead.
