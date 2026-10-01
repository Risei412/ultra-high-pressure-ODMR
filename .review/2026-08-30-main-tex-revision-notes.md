# Aggregated review — `paper/main.tex` against `PRA_high_pressure_NV_theory_revision_notes.md`

Date: 2026-08-30
Scope: full manuscript (1858 lines), `refs.bib` (51 entries), `build/main.pdf` (18 pp), the five figure
scripts and `code/tests/`.
Agents: technical-reviewer, logic-reviewer, consistency-checker, writing-reviewer,
latex-layout-auditor, bibliography-auditor (6, parallel).
Prior reviews treated as settled: `archive/reviews/2026-08-28-…`, `archive/reviews/2026-08-29-…`.

Orchestrator note: the memo's items were assigned across agents with a mandatory
CONFIRMED / ALREADY ADDRESSED / DISPUTE verdict so the memo itself was audited rather
than executed. Three memo items came back partly disputed, and the disputes are sound.

---

## 0. Verdicts on the revision memo

| Memo item | Verdict | Note |
|---|---|---|
| **A1** Eq. (6) `⇔` too strong | **CONFIRMED — understated** | `⇐` is not "missing a global condition"; it is **false**, and the paper's own case (c) is the counterexample. See §2.1. |
| **A2a** Morse regularity mid-proof | **CONFIRMED — understated** | Also: nondegenerate critical *values* never stated (and nearly violated); `A` is a piecewise-linear interpolant, so it is not Morse and the Milnor citation does not apply. |
| **A2b** `I<I_c` monotonicity | **CONFIRMED in substance, DISPUTED in form** | No new axiom needed: monotonicity below `Γ_p*` follows from `Φ(0⁺)=0` + *unique interior local maximum*. The real defect is the ambiguity of "taken here to be unique". |
| **A2c** `R ∝ A` needs a low-power condition | **CONFIRMED** | Plus a sharper gap the memo missed: T5 states no fixed-power convention; a fixed-photon-flux scan biases rung ratios by up to 29 %. |
| **A3** `W → W_eff` | **CONFIRMED** | The fix is **number-neutral**. Upgrade available: the knee in `P*(W_laser)` *locates* `Γ_ZPL(P)`. |
| **B1** window edge vs physical rung | **CONFIRMED — worsened today** | But excluding the edge would contradict Eq. (14), which makes the edge case theorem content. The illegitimate part is only that 402 nm is set by where the source figure was cropped. |
| **B2** ZPL height out of main claim | **ALREADY ADDRESSED for the sequence; CONFIRMED for the height** | Memo aimed at the wrong object: the residual dependence is on the height being usable at all, and §IV D contradicts App. B about it. |
| **B3** falsification-oriented tests | **CONFIRMED** (T11 DISPUTED) | 9 cells + 3 body leaks. T11 measures the antecedent directly, so "consistent with" is wrong there. |
| **B4** Theorem G domain | **PARTIALLY ALREADY ADDRESSED** | Conditions *are* stated at the theorem and in the Table I caption. Every restatement (Abstract, Intro, §III B italic, §IV C, Conclusions) drops them. Two further hidden conditions: same `n`, same functional class. |
| **B5** `E=3, n=2` scope | **CONFIRMED** | Table I caption localizes correctly; Abstract, §III A and Conclusions do not. |
| **C** keep G → M → X | **ENDORSE, premise DISPUTED** | G ⇒ M is a real dependency and is cashed. M ∥ X is juxtaposition, not chaining. §II C already delivers Theorem M's conclusion 300 lines early, undercutting §IV's "has never been asked". |
| **D** `pre-registered` | **CONFIRMED** | 5 occurrences. Recommended: `pre-specified` (keeps the temporal claim, which §VI's opening depends on). |
| **E** acquisition time | **CONFIRMED — Introduction is wrong** | Confirmed by 3 agents independently. |

---

## 1. Cross-validated findings (independent agreement)

| Finding | Agents |
|---|---|
| `\ref` destroyed into body text, `main.tex:1030` | 5 of 6 + orchestrator (bytes + PDF text) |
| Acquisition-time scaling inverted in Introduction | 3 |
| `ℓ_C` reported as `ℓ_G`; Abstract/Conclusions overclaim | 3 |
| Live `TODO(author)` on a load-bearing uncited claim | 3 |
| Fig. 2 caption's "every multiplicity below `I/I_c=3.6` changes" is false | 2 |
| Fig. 5 silently rescaled on the page | 2 (build log; CTM matrix) |
| `E=3, n=2` unlocalized outside Table I | 2 |

---

## 2. Critical

### 2.1 Eq. (6)'s `⇐` is false, and the stronger true result is already computed and pinned
`main.tex:395–406`, and App. A `:1631` leans on the false direction.

`ℓ_G(λ_PL)=0` gives only that `λ_PL` is a critical point of `Φ`. What the paper needs is
`ℓ_G ≡ 0` on the window, which is the *next* sentence (`:401–403`) — a different, sufficient
condition. Case (c) at `:477–481` is an explicit counterexample: under (M), `λ_PL = λ_abs` is an
interior maximum of `A`, so `ℓ_G(λ_PL)=0` automatically, yet the paper says `λ_abs` is then the
locally *worst* member of the optimal set.

**The repository already contains the quantitative result the paper does not report.**
`theory_a1_generalization._zpl_jump_threshold` → `zpl_jump_l_G = 0.268 %/nm`, pinned by
`test_theory_a1_generalization.py:101`. Above it the global optimum moves 440.64 → 514.46 nm,
a **73.8 nm** jump. Verified by the orchestrator by direct execution.

*Orchestrator caveat:* the computation is parameterised by `ℓ_G`, but the measurement is of `ℓ_C`.
"The cited 0.43 %/nm already exceeds the threshold" holds only if the linewidth does not cancel
the contrast slope — the very assumption A1 is about. State it conditionally:
*if the linewidth does not carry the same slope, the separation is not a perturbation but a
branch change of 74 nm.*

Fix: split into (6a) stationarity (`⇒` only, with interiority and differentiability stated) and
(6b) sufficiency (`ℓ_G ≡ 0` on the window). Use `6a`/`6b` subnumbering so Eqs. (7)–(23) do not shift.

### 2.2 Table IV's absorption Debye-Waller factor is 37 % off the source's own stated value
`main.tex:1119–1149`; data `code/data/ho_fig1_panels_bc.csv`.

`Ho2026` states in running text: DWF_abs **2.2 %** at ambient and **0.36 %** at 120 GPa, and
describes the collapse as "more than a fivefold reduction" (= 6.11×).

| | ours | source text | Δ |
|---|---|---|---|
| `S_abs` 0 → 120 GPa | 3.023 → 4.553 | 3.08 → 4.61 | −1.9 %, −1.2 % (harmless) |
| `DWF_abs` 0 GPa | 0.0205 | 0.022 | −6.8 % |
| **`DWF_abs` 120 GPa** | **0.00226** | **0.0036** | **−37 %** |

Consequence: the quoted "factor 9.07" should be **6.11**, and the multi-mode margin against
`e^{ΔS}=4.62` falls from 2.0× to **1.3×**. The conclusion survives; the number does not.

Two diagnoses, both pointing at panel (c) only (`S_abs` from panel (b) is clean):
- *Agent's hypothesis:* ~20 GPa pressure-axis offset — our 100 GPa value (0.336 %) is within 7 %
  of their stated 120 GPa value (0.36 %).
- *Orchestrator's alternative, which fits both endpoints better:* a log-axis calibration stretch.
  A pure 20 GPa shift does not explain the ambient point (2.05 vs 2.2 is only −7 %), whereas
  `ln 9.07 / ln 6.11 = 1.22` — a ~22 % stretch of the decade calibration anchored near ambient —
  reproduces both ends.

Either way panel (c) must be re-extracted and anchored to the two values the source states in text.
Then re-run App. B's K3 check, which draws on the low-pressure end of the same row.

Referee exposure: PRApplied will plausibly send this to Roch, Razinkovas or Loubeyre — the authors
of `Ho2026`. Appendix B exists to survive exactly this audit.

### 2.3 `\ref` destroyed into body text — visible in the compiled PDF
`main.tex:1030–1031`. Bytes are `S e c . ~ \n e f {` — the backslash and `r` are gone.
`build/main.pdf` p. 9 prints **"the gauge degeneracy of Sec. efsec:G:identifiability,"**.
No LaTeX error is raised because `\ref` is never called, which is why the 08-29 "0 undefined
references" check passed. Not in `HEAD`; introduced by the 08-29 restructuring pass, i.e. a
regression on the sentence that review deliberately added.
**Grep for other `\r`-initial macros (`\rho`, `\right`, `\raggedright`) before the next build.**

### 2.4 Eq. (14)'s window-edge rule is not generally true
`main.tex:791–799, 825–828`; `code/theory_a2_multiplicity.py:62` hard-codes `-1`.
The sign depends on the endpoint slope: `ΔN = −1` where `A` increases into the interior,
`+1` where it decreases. **No value changes for this kernel** (both edges increase inward), but
the rule is stated as a general result in a boxed theorem.

### 2.5 Table II's fractions are cut by the `\colrule`
`main.tex:642–652`, p. 5, confirmed at 400 dpi: the `\tfrac12` denominator of the *all three* row is
struck through by the rule, and the `\tfrac13` in the `ρ*` column is cut in half. `ruledtabular`
rows are too tight for two-level fractions. Fix: slashed form `$1/2$`, `$1/3$`, `$1/5$`.

### 2.6 Two placeholders in the typeset output
- `ACKNOWLEDGMENTS — [To be completed.]`, p. 14.
- `% TODO(author)` at `main.tex:962–967`, on the claim "the ZPL is absent from the excitation
  wavelengths conventionally considered for high-pressure work" — which supports the paper's most
  concrete experimental recommendation and carries no source, as the comment itself says.

---

## 3. Important

### Theory / logic
1. **A2 restructure** → Theorem M1 (level set) + Corollary M2 (staircase) + Corollary M3
   (calibration-free rungs). Keep the existing `\label`s so no `\eqref` breaks. Replace "Morse
   function" + `\cite[Thm. 3.1]{Milnor1963}` with the piecewise-`C¹` argument the code actually
   implements. Note the naming clash: `code/tests/` already uses M1/M2/M3 for other things.
2. **New hypothesis (L)**: `R(Γ_p) = R'(0)Γ_p + o(Γ_p)`, `R'(0) > 0`, background subtracted.
   Independent of (M) and load-bearing for "calibration free".
3. **T5 must state fixed incident optical power.** A fixed-photon-flux scan returns `σ_abs`, not
   `λσ_abs`; rung ratios are then wrong by `λ_k/λ_j`, up to **29 %** over 402–517 nm, against the
   ±2 % the manuscript assigns to `a_k`.
4. **Theorem M has no high-power boundary case.** Once `I_c/I < min(a(λ_b), a(λ_r))` the level set
   is empty. `theory_a1_generalization.py:222–228` tracks this (`truncated_blue/red`); the theorem
   does not mention it.
5. **`S` in Eq. (20) is `S_total`, not `S_abs`** — K4 establishes `DWF = e^{−S_total}` with
   `S_total = 3.887 → 6.092`, while §V A "confirms the driver" with `dS_abs/dP`. Three objects
   wear the letter `S`. The conclusion survives (both monotone).
6. **`I/I_c` maps to `ς`, not `x`.** §III A sets `Γ_p ≡ x` (so `n=2`), but the ladder's power ratios
   are ratios of `ς`, with `x = ς/(1+ς)`. Never stated; a reader applying Eq. (15) with `x` gets the
   wrong power.
7. **Hypothesis (M) is not currently falsifiable as stated** — no domain, no tolerance, and the
   manuscript itself calls it "an approximation everywhere". State it as (M|Ω,ε) and let T6 report
   the smallest ε. This is the same move B3 asks for in Table V.
8. **§II C delivers Theorem M's conclusion 300 lines early**, using `Γ_p*` before it is defined,
   which undercuts §IV's "What has never been asked is what the solution set looks like".
9. **§IV D vs App. B contradiction on the ZPL height**: `:937–941` says the range "indicates a scale,
   not a limit"; `:968–971` then says "on either reading it belongs on the list", treating it as
   exhaustive, and the interval propagates unconditioned into the Conclusions' 4–21 % claim.
10. **`0.01 %` in the text is pinned at `0.1 %` by the test** (`test_theory_a2_multiplicity.py:42`) —
    a violation of the manuscript's own header rule.
11. **"about 10 GPa at fixed photon flux"** (`:1228–1230`) is supported by no live code. The only
    fixed-flux number, 12.3 GPa, belongs to the withdrawn clipped-ZPL path.

### Framing
12. Introduction leans on the weaker evidence (Todenhagen contrast slope, ambient, 480–575 nm) and
    withholds the stronger (Dréau `En = 6 > 1`, independent of `Δν`). The Todenhagen range is also
    **disjoint** from the 426.4–457.9 nm blue window, and contains both known (M)-violating
    resonances (521, 575 nm). Not acknowledged.
13. Abstract states a conditional number unconditionally ("jumping roughly 70 nm").
14. "Both antecedents are generic for a color center under compression" (`:1101–1105`) generalizes
    from n = 1 host. The sentence a referee will use to argue PRA-2 rather than PRA-3.
15. §VII E does not deliver the "role of pressure is inverted" claim the Introduction promises.
16. Confidence gradient is **inverted**: own-reconstruction claims meticulously hedged;
    published-data claims flatly asserted ("confirmed", "real", "This is not an assumption").

### Figures / floats (most of these are from today's rework)
17. **Fig. 5 rescaled 0.96762** — `bbox_inches='tight'` expands the canvas to hold the protruding
    twin label, then the page shrinks it. Legend/annotation lands at **6.77 pt**, below the 7 pt
    floor `fig_style.py` declares, and Fig. 4/Fig. 5 differ by 4 % on facing pages.
    Fix `fig_style.save` to `bbox_inches=Bbox([[0,0], fig.get_size_inches()])`.
    **The guard test is also wrong** — it measures the canvas before tight-bbox, not the saved width.
18. **Fig. 2(b) rung numeral "1" is overprinted by the rung-2/3 dashed rule** (the +2.5 pt offset
    pushes it onto the neighbouring line).
19. **Fig. 3(b) legend is struck by the `W = 9.7 meV` window rule** and sits unframed on the grey band.
20. **Fig. 5 caption mislabels `S_total`** as the Jahn-Teller-active remainder; that is `S_JT`.
21. **Fig. 4's in-figure legend "Ho Fig. 5(b)" now collides with this paper's own Fig. 5**, created
    by today's split. Neither caption nor body introduces Ho's Fig. 5(b).
22. **Fig. 3(b) caption** says the square markers are "the six bandwidths annotated in the panel";
    the code annotates `P*` (120, 104, 78, 45, 24, 7 GPa). The caption also omits the two green
    window verticals and says "the shaded region" where two `axhspan`s are drawn.
23. **Fig. 2 caption** implies six numerals; `rung_labels()` emits five (`2,3` merged). Explanation
    lives only in the Table III caption.
24. **Table IV prints a page before its first citation**, and first-citation order is I → III → II → IV → V.
25. **Unnumbered sixth table** at `main.tex:749–759` on the same spread as Table II.
26. Table IV discussed in §V A without a `\ref`; Fig. 5 referenced once, bare; **Fig. 2(a) never
    discussed in the body**.

### Bibliography
27. `Hao2025` is missing an author (**Xiaobing Liu**, between Gang-Qin Liu and Yanming Ma).
28. `Goldman2015` carries PRB 91's title/DOI with **PRL 114's author list**. (Uncited, but wrong in the file.)
29. `Todenhagen2023`: the "525–550 nm" claim at `:253` **is not in the preprint**; and "over its
    steepest stretch" misdescribes 0.43 %/nm, which is the whole-span mean. The entry cites the
    **APL 2025** version under a changed title — every number was verified against the 2023 preprint.
    Re-verify against APL 126, 194003 before submission.
30. **Taylor et al., Nat. Phys. 4, 810 (2008) is absent from the file entirely** — the paper that
    established the sensitivity scaling this manuscript optimises.
31. `Doherty2013` (Phys. Rep. NV review) and `Huang2025ISC` (now PRL 137, 093801) are in the file and
    uncited; the latter closes a limitation the author's own notes flag (`main.tex` has zero
    occurrences of "intersystem").

**Bibliography positives:** 45/45 DOIs resolve with matching metadata; 44/45 author lists match
Crossref exactly; abbreviations uniform; brace protection complete and `longbibliography`-safe;
no duplicates, no undefined citations; the 3 BibTeX warnings are cosmetic artefacts of correct
preprint entries. All three arXiv-only entries are genuinely still preprints.

---

## 4. Minor (selected)

Prose: `pre-registered` → `pre-specified` (5×); `gauge` never defined and the metaphor is inexact
(motion in the plane *is* observable — Fig. 1(b) shows a factor 2 — only `η` is invariant); `rung`
never defined though it is formal terminology and appears in the Abstract; `knob` → `control
variable`; dangling "Conditional on the reconstruction" openers in the Table II and Table III
captions; ODMR never expanded, ZPL and DWF used before expansion; Oxford comma ~50/50;
`\subsection{What this is not}`; B7 offenders at `:1420–1426`, `:1431`, `:1591–1599`, `:979–989`,
`:239–244`.

Figures: `linewidth^1` prints a literal caret; `(normalised)` vs `normalized`; Fig. 5(b) `ylim(0,4)`
on data at 2.53–2.59 meV wastes 60 % of the panel and hides the 2.4 % the caption claims;
Appendix C's "Figs. 2, 1, 3, 4 and 5" ordering; `colorlinks=true,allcolors=blue`; missing space
"unremarkable.Its" on p. 7; four badness-10000 underfull lines in the Appendix C bullet list.

Notation collisions: `Γ` (half-scale / pump rate / linewidth / rates), `W` vs `w`, `E` vs `E_γ`,
`r` (branch ratio / correlation / rate coefficients), `p*` vs `P*`, `a` vs `a_es`, `s` vs `S`.

Naming: "window edge" / "blue edge" / "edge of the analysis window" / `'edge-blue'`;
"analysis interval" / "analysis window" / "plotted window" / `DATA_WINDOW`;
"rung" / "step" / "transition power"; "level set" / "optimal set" / "solution set".

**PASS, cleanly:** negation-contrast audit (zero `not X but Y`), AI-writing tells (one "Moreover" in
1858 lines), float definition order, all floats referenced, section structure vs the header comment,
`kernel`/`envelope` term separation, Type 3 fonts (none), overfull boxes (none), missing characters
(none), and — apart from §2.2 — the full numeric cross-check against `code/tests/`.

---

## 5. The nugget

*Once you optimize for sensitivity instead of brightness, the optimal excitation wavelength stops
being a wavelength: under rate mediation it is a level set whose size steps with power and whose
carrying branch swaps with pressure — so what is to be measured is the optimum's structure, not
its location.*

The Conclusions state this in nine words (`:1588–1589`, "Relaxing it does not displace the optimum;
it changes what kind of object the optimum is"). The Introduction has it at `:260–261`. **The
Abstract does not** — a reader who stops there comes away with "blue is better at high pressure and
there is a jump", the PRA-1 reading of a PRA-3 paper. Highest-leverage single revision available:
lift that sentence into the Abstract.

Second nugget, unresolved: Theorem G's non-identifiability is the most reusable result here and is
currently framed as scaffolding. Title leads with X, Abstract leads with X, body leads with G,
Conclusions lead with M. Choose one.

---

## 6. Applied on 2026-08-30 (mechanical pass)

Scope chosen by the author: mechanical corrections only, plus the `Ho2026`
panel-(c) re-extraction. Theorem statements (A1, A2, A3, Eq. (14)) were **not**
touched.

### Critical #1 — the DWF discrepancy: diagnosed and measured, not yet applied

Re-extracted from the source's own vector paths (`pdftocairo -svg`, page 2),
not by pixel tracing. Panel (c) is a **linear** axis (0.00 / 0.02 / 0.04), so
neither reviewer hypothesis holds: there is no ~20 GPa pressure-axis offset and
no log-calibration stretch. Both rows carry a **constant baseline offset** —
`S_abs` low by 0.060, `DWF_abs` low by 0.00138 — invisible at ambient (−6 %)
and severe at 120 GPa (−60 %), because 0.0036 is 7 % of full scale on a linear
axis.

| | re-extracted | committed CSV | source running text |
|---|---|---|---|
| `S_abs` 0 → 120 GPa | 3.081 → 4.613 | 3.023 → 4.553 | 3.08 → 4.61 |
| `DWF_abs` 0 → 120 GPa | 0.0218 → 0.00363 | 0.0205 → 0.00226 | 0.022 → 0.0036 |
| collapse | ×6.01 | ×9.07 | "more than a fivefold" |

Calibration validated independently on the emission curve: 0.0489 against the
stated 0.049, and 0.0091 against "below 1 %".

Unchanged by the correction: `dS_abs/dP = 1.27×10⁻²`, `e^{ΔS} = 4.62`.
Changed: the collapse factor, hence the multi-mode margin 2.0× → **1.30×**; and,
because `r ∝ 1/DWF_abs`, the whole Theorem X worked example — `r` 103–672 →
**97–419**, growth ×6.51 → **×4.31**, critical bandwidths 9.68/1.49 →
**10.29/2.39 meV**, the window 1.5–9.7 → **≈2.4–10.3 meV**, `P*(W)`, Fig. 3,
Table IV, and the pinned values in `test_figure_validation.py` and
`test_theory_figures.py`. Held pending an explicit decision, since this exceeds
a mechanical fix.

### Applied

`paper/main.tex`
- restored the destroyed `\ref` (the file held a bare CR; `Sec. efsec:G:…`
  was printing on p. 9) and verified no stray CR remains anywhere
- Introduction acquisition-time scaling: `divided by` → `times the square of`
- Fig. 2 caption: "every multiplicity below `I/I_c = 3.6`" → "between 1.44 and
  2.06"; noted the merged `2,3` numeral; noted that rung 4 is drawn in gray
- Table III: rung 4 relabelled `edge` with a footnote saying it is not a
  critical point of `A`; the sequence at first quotation now flags the
  censoring step in place rather than 77 lines later
- Fig. 5 caption: `S_total` is the sum, `S_JT` the Jahn-Teller remainder
- Fig. 1 caption: the bars are raw values, not normalized to the first model
- Table I and Table III captions: dangling "Conditional on the reconstruction"
- Appendix C: singular `\ref`, correct target, figures in numerical order
- Table II: `\tfrac` → slashed fractions (the `\colrule` was cutting them)
- the uncited `TODO(author)` claim restated as "we are not aware of"
- ODMR, ZPL and DWF expanded at first use

`code/`
- `fig_style.FULL_WIDTH_IN` = 510 pt exactly; `save()` writes at `figsize` and
  `verify_width()` re-reads the PDF and refuses any other width. All five
  figures are now 510.0 pt, scale 1.000
- the guard test measured `get_size_inches()` — the canvas *before* saving —
  so it could not have caught the 7.29 in file. Replaced, in both test files,
  by an assertion on the written PDF
- Fig. 2(b): rung 4 drawn gray and dotted with its own legend entry (B1); the
  numeral `1` moved left of its rule instead of onto rung 2/3's
- Fig. 3(b): legend framed and opaque, clear of the window rule
- Fig. 1: `linewidth^1` no longer prints a literal caret; `normalised` →
  `normalized`
- Fig. 4: legend reads "Ho *et al.* Fig. 5(b)", no longer colliding with this
  paper's own Fig. 5
- margins reset for the fixed bounding box; nothing clipped

`paper/refs.bib`
- `Hao2025`: restored the dropped author Xiaobing Liu
- `Goldman2015`: PRL 114's author list replaced by PRB 91's

Build: 18 pages, **0 overfull boxes**, 0 undefined references.

### Not applied — needs the author

- `ACKNOWLEDGMENTS [To be completed.]`
- the unnumbered sixth table at `main.tex:749` (promote or dissolve)
- `Todenhagen2023`: the 525–550 nm claim is absent from the preprint, and
  "over its steepest stretch" misdescribes a whole-span mean. Needs APL 126,
  194003 (2025) to settle
- Taylor *et al.*, Nat. Phys. **4**, 810 (2008) is absent from the file
- Table IV printed a page before its first citation; Table III cited before
  Table II
- everything in §2.1, §2.4 and §3 above (the theorem-level items)

---

## 7. The DWF propagation (option 1), and what it uncovered

The correction was applied. Its consequences are larger than the estimate in §6,
and one of them is a finding in its own right.

### What the data now says

`code/extract_ho_fig1_panels_bc.py` (new) reads panels (b) and (c) from the
source's vector paths and checks itself against the four values Ho *et al.*
state in their running text. All four agree to ≤0.9 %. The committed CSV is
regenerated from it.

### The blast radius, measured

| | before | after |
|---|---|---|
| `S_abs` 0 → 120 GPa | 3.023 → 4.553 | **3.081 → 4.613** |
| `DWF_abs` | 0.0205 → 0.00226 | **0.0218 → 0.00363** |
| DWF collapse | ×9.07 | **×6.01** |
| multi-mode margin vs `e^{ΔS}` | 2.0× | **1.30×** |
| `r(P)` | 103.3 → 672.0 | **97.1 → 418.4** |
| growth of `r` | ×6.51 | **×4.31** |
| critical bandwidths | 9.68 / 1.49 meV | **10.30 / 2.39 meV** |
| accessible window | 1.5–9.7 meV | **2.4–10.3 meV** |
| `S_total` | 3.887 → 6.092 | **3.826 → 5.619** |
| `dS_abs/dP` | 1.27×10⁻² | **1.27×10⁻² (unchanged)** |
| **K3 ratio, kernel ZPL / published DWF** | **0.96 → 1.43, monotone** | **0.90 → 0.89, flat** |

Unaffected: everything fed by `ho_fig1e_absorption.csv` — the kernel `A(λ)`,
Table I, the 5 % band, the 440.6 nm optimum, and Theorem M's face-value ladder.
The two data files are separate, and the paper's headline number never touched
panels (b) or (c).

### The finding: the ×1.43 over-weight was an artefact

K3 compares the kernel's own zero-phonon area with the published Debye-Waller
factor. On the old extraction that ratio climbed 0.96 → 1.43 and was read as
*the reconstruction over-weights the ZPL, monotonically, and may not be used
above 40 GPa*. On the corrected extraction it is **0.88–0.91 at every pressure
with no trend**, and the two collapses agree to 1.3 % (kernel ×6.08, published
×6.01). What was a pressure-dependent defect is a constant 10 % offset.

The direction inverts: the kernel's ZPL is not 43 % too *high* at 120 GPa, it is
about 11 % too *low* everywhere. So the ZPL height correction is an inflation by
1.12, not a deflation by 1.43; `a(λ_ZPL)` goes 0.6938 → 0.7796, still above the
475.55 nm maximum, so **the rung ordering does not change and the deflated
branch does not arise**. The uncertainty interval collapses from
`1.44 ≤ I/I_c ≤ 2.06` to roughly `1.28 ≤ I/I_c ≤ 1.44`.

This is good news for the manuscript — the largest caveat in Sec. IV D mostly
evaporates — but it is a re-derivation, not an edit, so it was not applied.

### Applied

- `code/extract_ho_fig1_panels_bc.py` (new, self-auditing against the four
  published anchors) and the regenerated `code/data/ho_fig1_panels_bc.csv`
- `fig_x_branch_exchange.py`: marked bandwidths moved into the new window
  (2.5, 3, 4, 5, 7, 9 meV → `P*` = 117, 101, 74, 57, 31, 13 GPa)
- `main.tex`: Table IV all four rows; Sec. V A drivers paragraph; Fig. 3
  caption (growth, endpoints, markers, window); Sec. V B growth, `P*` values
  and the narrowing bandwidth; App. B K4; Fig. 5 caption
- `tests/`: `test_figure_validation.py`, `test_kernel_sanity_checks.py`,
  `test_theory_figures.py` re-pinned, with the K3 test renamed to state what
  it now measures

### Deliberately not applied — these are arguments, not numbers

Each of these currently rests on the withdrawn ×1.43:

| location | what it asserts |
|---|---|
| Table I caption | "the overweight of up to ×1.43" |
| Sec. IV D | the deflation of the ZPL rung; Eq. (18); the `N = 1` plateau; "one plateau must not be claimed" |
| Fig. 2 caption, Table III caption | the `1.44 ≤ I/I_c ≤ 2.06` interval |
| Sec. IV C | the deflated invariance probe sequence `2,2,2,1,5,3` |
| Conclusions | the 4 %–21 % first-rung reachability range |
| Limitations | "the ×1.43 overweight" |
| App. B K3 | "over-weights the ZPL monotonically, reaching 1.43 at 120 GPa"; "usable to 40 GPa" |

Until these are re-derived the manuscript is internally inconsistent: Table IV
and App. B K4 now carry the corrected data while the passages above still
carry the artefact.

---

## 8. Items 1–5 of the work order

### Item 1 — the excitation-power convention (T5): **the convention is right, the target is stale**

The convention is specified correctly in the experiment plan.
`docs/experiment/next_step_power_dependence_experiment.md:8` premises the frozen
claim on "同一入射光学パワー" — **fixed incident optical power** — which is the
convention `A = λσ_abs` requires. T5 itself stated no convention; it now does, and
so does the calibration-free subsection.

**But that document targets the superseded v1 prediction.** It states the optimum
as 475.5 nm with 473 nm at "感度損失約 0.2%", and its Stage 2 laser lines are
457 / 473 / 488 nm. Against the frozen v3 kernel the optimum is 440.6 nm, and
Table I gives:

| line | penalty | status |
|---|---|---|
| 445 nm | 1.004 | inside the 5 % band, essentially optimal |
| 457 nm | 1.045 | inside |
| 473 nm | 1.205 | **outside — 21 %, not 0.2 %** |
| 488 nm | 1.538 | far outside |

A Stage-1 power sweep at 473 nm and a Stage-2 scan at 457/473/488 would place one
point inside the band and would never visit 440–445 nm. **Before any of this data
is taken or used, that plan has to be re-cut against the v3 freeze.** This is a
larger risk to data usability than the power convention was.

### Item 2 — Sec. IV D re-derived under the corrected K3

`a(ZPL)` 0.694 → 0.780 (an inflation by 1.12, not a deflation by 1.43); the rung
moves 1.44 → 1.28 and **stays the lowest**, so the ordering and the sequence
`2→4→6→4→3→5→3` hold at both ends. Removed: Eq. (18), the reordered ladder, the
`N = 1` plateau, and the `[1.44, 2.06]` interval wherever it appeared (Fig. 2
caption, Table III caption, Sec. IV C's probe, Table I's ZPL row and caption,
Limitations). Reachability 4–21 % → **3–20 %**. The invariance probe now returns
the same six multiplicities at both ends, which it did not before. App. B's K3
paragraph rewritten, with a note recording that the ×1.43 trend was an extraction
artifact.

### Item 3 — Eq. (6), Theorem M, `W_eff`

- **Eq. (6)** split into (6a) stationarity (`⇔`, with interiority and
  differentiability stated) and (6b) sufficiency (`ℓ_G ≡ 0` on the window ⇒
  coincidence). Subequations, so (7)–(23) keep their numbers. Case (c) is named
  as the counterexample to the converse of (6a), and the 74 nm jump above
  `ℓ_G = 0.27 %/nm` is now reported — conditionally on the reconstruction, and
  with the `ℓ_C` vs `ℓ_G` proviso stated. App. A retargeted to (6b).
- **Theorem M** now assumes: `Φ` continuous, `Φ(0⁺) = 0`, a **unique interior
  local maximum**; `A` continuous and piecewise `C¹` on the closed window with
  finitely many strict local extrema, **pairwise distinct critical values** also
  distinct from the endpoint values, and monotone near each endpoint. The proof's
  Morse framing is replaced by the elementary piecewise-`C¹` count the code
  implements — `A` is a piecewise-linear interpolant and is not Morse — with
  Milnor kept as the smooth-case reference. **Eq. (14)'s edge rule now carries
  its sign condition** (`−1` inward-increasing, `+1` inward-decreasing); both
  edges here are inward-increasing, so no number changes.
- **Hypothesis (L)** stated: `R(Γ_p) = R'(0)Γ_p + o(Γ_p)`, background subtracted,
  independent of (M) and load-bearing for the calibration-free claim.
- **`W → W_eff`** throughout, with `W_eff ≃ max(W_laser, Γ_ZPL(P))` as a new
  equation. Turned into a prediction rather than a caveat: `P*(W_laser)` falls
  monotonically while `W_laser > Γ_ZPL` and flattens below it, so **the knee
  locates `Γ_ZPL(P)`** — which K2 says cannot be read from the reconstruction.
  Number-neutral, as promised.

### Item 4 — Results in G→M→X order, kernel demoted: **blocked, correctly**

Not attempted. Writing a Results section without results would be fabrication.
The scaffolding is what the framework already specifies (T1–T3/T6 → Result 1,
T4/T5 → Result 2, T8–T11 → Result 3), and the demotion is the same edit in
reverse: once `a_k` are read from a measured spectrum, App. B stops defending the
reconstruction and becomes a note on where the pre-registration came from. Both
wait on data — and, per item 1, on a re-cut experiment plan.

### Item 5 — Table V, scoping, housekeeping

- Table V: all ten decision cells rewritten from proof language to falsification
  language; the epistemic rule stated once in the caption at no cost in row
  height. T11's first clause left as an assertion, because it measures Theorem X's
  antecedent directly rather than inferring a hypothesis from agreement.
- Three body leaks fixed ("(M) holds broadly" → "consistent with (M)"; T7
  "verifies" → "constrains"; Theorem X "intact" → its antecedent holds).
- `pre-registered` → `pre-specified`, 5 occurrences.
- **B4**: "within the common-scale power-law response class" added to the Abstract,
  the Introduction preview, the Sec. III B italic, and the Conclusions; the two
  further hidden conditions (fixed `n`, one functional class) stated once.
- **B5**: `E = 3, n = 2` localized to the ambient-pressure single-NV fit in the
  Abstract, Sec. III A and the Conclusions. "This is not an assumption: it has
  been measured" replaced by what was actually done.
- Appendix A's "only 457 nm falls inside the 5 % band" corrected — 445 nm is
  inside at a penalty of 1.004.
- **Taylor et al., Nat. Phys. 4, 810 (2008)** added and cited at Eq. (1).
- Acknowledgments: the `[To be completed.]` placeholder removed from the output,
  with a commented skeleton left in the source. Funding is a disclosure and
  cannot be invented.

Build: **19 pages** (up from 18), 0 overfull boxes, 0 undefined references,
3 cosmetic BibTeX preprint warnings.
