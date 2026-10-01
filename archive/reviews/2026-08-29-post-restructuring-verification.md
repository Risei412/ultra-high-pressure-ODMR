# Post-restructuring verification — `paper/main.tex`

Date: 2026-08-29
Scope: `paper/main.tex` (1825 → 1843 lines), `paper/refs.bib`, `paper/build/main.pdf`, `code/`
Prior review: `archive/reviews/2026-08-28-main-tex-vs-prapplied-benchmark.md` (all four tiers applied 2026-08-28)

Agents: `consistency-checker`, `logic-reviewer`, `technical-reviewer`, `bibliography-auditor`
(all read-only; the orchestrator applied the fixes recorded below)

Items marked **[verified]** were checked by the orchestrator against the source or by running
code, not merely reported by an agent.

---

## Overview

The 2026-08-28 restructuring (11 → 8 sections, G moved before M, ~3,700 words cut) landed its
numerical content correctly: **C1, C2 and C3 are fully propagated [verified]**. What it left behind
is splice damage — contradictions between merged lead paragraphs and the subsections they
introduce, one subsection stranded 182 lines ahead of its own dependencies, and captions describing
figures that no longer match.

Two findings are new and more serious than anything in the prior review, and **both are in `code/`,
not the prose**: the C7 gauge-invariance repair is cosmetic, and Eq. (18) is not pinned by any test.
Together they falsify the manuscript's own claim at `main.tex:1758` that every number is
code-produced and pinned.

---

## Critical

### V1. The gauge-invariance check is still an identity **[verified]**
`code/theory_a2_multiplicity.py:314-347`, asserted at `main.tex:750` and `1776-1780`

The C7 repair gave each model its own `gamma`, so `i_c` differs. But the level passed to
`multiplicity` is `i_c / intensity` where `intensity = probe * i_c` — algebraically `1/probe`,
with `gamma` cancelling identically. The five models cannot disagree. `i_c = star / gamma_a_max`
never uses the model's `gamma` either (`gamma_a_max` defaults to 1.0).

Appendix C's claim that "the agreement of their ladders is a result rather than an identity" is
false, and `main.tex:750` ("We check the collapse numerically rather than assume it") is
unsupported.

**Not fixed** — deciding what the honest numerical content of this test is, is a scientific call.
A real check probes at a *shared absolute intensity grid* and normalises by each model's own
`I_c` afterwards, so the levels genuinely differ across models.

### V2. Eq. (18) and the 4–21% range are not produced by code and not pinned **[verified]**
`main.tex:1758` against `code/`

Grep of `code/theory_a2_multiplicity.py` for `0.485` / `deflat` returns zero hits. The deflated
ladder `2→4→2→1→3→5→3`, the five deflated plateau widths, and the 4–21% reachability figure exist
only in the manuscript. Worse, `code/tests/test_theory_a2_multiplicity.py:82` is still named
`test_zpl_is_the_first_rung_so_omitting_it_hides_the_effect` and pins `power_ratio == 1.441` as
*the first rung* — asserting the claim the erratum withdrew.

A referee who runs the suite finds the repository contradicting Sec. IV B.

**Not fixed** — needs a `deflated_ladder()` / `critical_values(zpl_scale=...)` path, tests pinning
Eq. (18) and both plateau-width lists, and that test renamed.

### V3. §III C was parasitic on §IV — **fixed**
It cited `eq:ladder-deflated` (defined 182 lines later) and quoted multiplicity sequences the
reader had not met. The G-before-M reorder was performed to remove a forward reference and created
a worse one. Moved to the end of §IV, after "The ladder is reachable", where all three objects
exist. Label `sec:G:invariance` was unreferenced, so no cross-reference broke.

### V4. Antecedent count contradicted itself three ways — **fixed**
`main.tex:1378` said four antecedents / three measured; §VII A (`1385`) and the Conclusions
(`1561`) both said three / two. The §VII lead was the odd one out; corrected, and extended to
announce all six subsections rather than two.

---

## Important — fixed

| Finding | Sites | Fix |
|---|---|---|
| Table I captioned "interior critical points … on which every rung sits" but lists only the four maxima; rungs also sit on two interior minima and the window edge | caption, block header, `§II F` | Relabelled "interior maxima"; claim corrected and pointed at `tab:ladder` |
| Fig. 1 caption: "the five response models of Table II" — Table II has only two `E=2` rows | Fig. 1 caption | "five representatives of the plane `E=2`, two of which appear in Table II" |
| Fig. 1(a) caption says "log axes"; the ordinate is linear (`fig_g_gauge_degeneracy.py` sets no `set_yscale`) | Fig. 1 caption | "on a logarithmic abscissa" |
| Fig. 3(b) caption: "the six bandwidths tabulated in the text" — text gives three; dissolved Table V's content never reached the caption | Fig. 3 caption | Six values listed **[verified against `fig_x_branch_exchange.py:68`]**, plus the fixed-power convention and ±2 GPa; "accessible" → "published" to match the in-figure annotation **[verified]** |
| Fig. 2 caption self-references its own section and says "only the zero-phonon step moves" — four of seven multiplicities change | Fig. 2 caption | Points to Eq. (18); "that step moves past the blue edge and every multiplicity below `I/I_c = 3.6` changes" |
| Appendix B: "Every number in Secs. IV and V descends from one reconstruction" — contradicted by Table IV's caption and by `1696-1697` | App. B lead | Exceptions named (reachability estimate, Table IV drivers) |
| Duplicated sentence inside one paragraph ("… is gauge fixing" twice) | §III | Second deleted |
| "Three checks support the two-branch identification. **Both** appear as …" | §V A | "First, the two branches appear as …" |
| "by case (b) above" — pointer ~1000 lines back, and the logic was wrong (case (c) is precisely when the optima separate) | §VII C | Restated as "enters the three cases of Sec. II E through `A` alone and cannot itself be the source of a separation", with a real `\ref` |
| Abstract's "First / Second / Third" signalled an order the body reverses | abstract | Ordinals removed; the deliberate X-primary framing preserved |
| Negation-contrast (B2/F2): "the optimum is **not** a point **but** the level set" | abstract | Rephrased positively |
| Four `\SI{…}{}` calls with an empty unit argument | 4 sites | `\SIrange` / `\SIlist` / explicit units |
| Precision contradicted the stated ±1 nm: `440.65^{+17.25}_{-14.22}`, band to 0.01 nm | 3 sites | One decimal throughout; `440.6` everywhere |
| Table I's ZPL row: `a = 0.49–0.69` with penalty `1.20–1.44`, but `1/√0.49 = 1.43` | Table I | `0.485–0.694`, which reproduces 1.20–1.44 |
| Orphan labels `app:optical` and `sec:phase`, never referenced | 3 sites | Referenced from §II E, the Conclusions and the Introduction |
| Conclusions quoted "a factor 1.45" with no conditionality flag, the manuscript's own rule | Conclusions | Flagged and pointed at Appendix A |
| "Weaker, because the whole response enters through one number" — a simplification, not a weakening | §III lead | "Easier to meet … so a single measured exponent decides it" |
| `I_c` used in §III, defined in §IV | §III | Glossed inline with a forward `\eqref` |

### Bibliography — fixed
- **`Barfuss2017` removed.** Its author list (Barfuss, Kasperczyk, Kölbl, Maletinsky) did not match
  its DOI/title/volume, which Crossref returns as MacQuarrie *et al.* — two papers conflated into
  one entry. It was uncited, and `main.tex` contains no occurrence of "strain" or "stress", so
  nothing could have cited it.
- **Four preprints rendered their identifier twice** in the PDF
  (`arXiv:2510.26605 10.48550/arXiv.2510.26605`). Converted to `eprint`/`archivePrefix`.
  apsrev4-2 requires an `@article` to *have* a `journal` field, so an empty one is retained —
  this yields three benign "empty journal" BibTeX warnings and correct output **[verified in the
  rendered PDF]**.
- **`Huang2025ISC` updated** to its published version: Phys. Rev. Lett. **137**, 093801 (2026),
  DOI `10.1103/hxtk-vmjq`, published 26 Aug 2026 — three days before this review.
- `J. Phys. D` → `J. Phys. D: Appl. Phys.`, the only venue-abbreviation outlier.

---

## Still open

**Scientific calls, deliberately not made here:**

1. **V1** — repair the gauge check, or withdraw the claim that it is a result.
2. **V2** — pin Eq. (18) in `code/`, or soften `main.tex:1758`.
3. **Theorem X's proof premise is false on the deflated branch.** `main.tex:1036-1038` asserts that
   between the two branch maxima `A` is strictly smaller than at either; on the deflated branch
   `a(475.6) = 0.660 > a_ZPL = 0.485`. The conclusion survives (`0.660 < 1`), the justification
   does not. Replace with: let `λ_SB` maximise `A` over the sideband branch and `λ_ZPL` over the
   zero-phonon branch; the global maximiser is then in `{λ_ZPL, λ_SB}` by definition.
4. **Theorem X lacks a hypothesis on `W`.** "At most one zero" needs `W` fixed in `P`, but `W` is
   defined as `max(laser linewidth, Γ_ZPL)` and check K2 declares `Γ_ZPL` unreadable from this
   source. State it as an explicit assumption.
5. **`D` names two functions** — `D(P;W) = ln[r(P)W]` in the theorem, and
   `ln[(1-e^{-S})Γ_ZPL/(e^{-S}Γ_SB)]` in Eq. (20). Subscript them.
6. **`h` names two things** — Planck's constant in Eq. (1), the mediated response `η = h(Γ_p)` in
   §II D and §IV.
7. **The ZPL-omission claim carries no source.** Marked with a `TODO(author)` in the file; it is
   load-bearing for "it belongs on the list". An earlier draft of this pass inserted three
   plausible citations and they were withdrawn as unverified.
8. **12 uncited entries** remain (13 minus `Barfuss2017`). The auditor's priorities: cite
   `Doherty2013` (canonical NV review, currently absent entirely) and `Huang2025ISC`; cite
   `StruckFonger1975` at the hot-band equation, which is an uncited source for a numerical claim.
   `Occelli2003`, `Akahama2004`, `Bhattacharyya2022` lean toward removal.
9. **`\begin{acknowledgments}\emph{[To be completed.]}\end{acknowledgments}`** still renders in the
   PDF. Submission blocker; only the authors can write it.
10. **Duplication the merge exposed**: §VII D repeats §VII C point for point; §VII E restates §V
    and the Introduction; the Todenhagen measurement appears with the same three numbers in §I and
    §II C.

**Header comment navigation is stale** — it cites "Sec. IV D", "Sec. VI B", "Sec. II D" and
"Eq. (13)" for material that moved. Documentation only, but it is the first thing a co-author reads.

---

## Build state

| | before | after |
|---|---|---|
| pages | 17 | **18** |
| errors | 0 | **0** |
| undefined references | 0 | **0** |
| overfull boxes | 0 | **0** |
| stuck floats | 1 cluster | **0** |
| BibTeX warnings | 0 | 3 (benign empty `journal` on preprints) |
| tests | 213 pass | 213 pass |

The stuck-float cluster carried over from the prior review resolved on its own once the captions
were re-cut. The page count rose by one: the corrections add content (the §VII roadmap, the six
bandwidths, the Table I and Appendix B qualifiers).

**Note on running the tests:** the suite must be invoked from `code/`. From the repository root,
12 tests error and 1 fails on missing data paths — a working-directory artifact, not a regression.

---

## Second pass — breaks created by the §III C relocation (V3), then repaired

The logic reviewer re-ran against the edited file and confirmed C1–C4, I9–I12 and the "Weaker,
because" non sequitur as resolved. It also caught three breaks that **the V3 move itself created**:

- **The anaphor was stranded.** §IV D opened "The same degeneracy that obstructs mechanism
  attribution…", whose referent (§III B) was now ~260 lines back with two subsections in between.
  Named explicitly: "The gauge degeneracy of Sec.~\ref{sec:G:identifiability}…".
- **§III lost its positive closer**, which had travelled with the moved subsection, leaving the
  section ending on pure negation with no hand-off. Added one sentence: the ladder's transition
  powers depend on the response only through `Γ_p*`, so the degeneracy that defeats mechanism
  attribution leaves the next section's prediction untouched.
- **`\label{sec:G:invariance}` sat inside `sec:M`.** Renamed `sec:M:invariance` (unreferenced, so
  nothing broke).

Also taken while there: **Theorem M now invokes Theorem G** instead of silently re-assuming its
conclusion — the entire point of the G-before-M reorder, which the theorem statement had not
cashed.

Rebuild after: 18 pp, 0 errors, 0 undefined, 0 overfull, 0 stuck floats.

**Lesson for the next structural move:** relocating a subsection strands (a) anaphora pointing at
its old neighbour, (b) whatever closing work it was doing for the section it left, and (c) its
label's namespace. Check all three at both seams.
