# Theory freeze v2 — 2026-08-27

Computed from the tree at commit `4e60c83` on branch
`claude/repository-theory-presentation-rm3txv`.

This is the second freeze. The first (`docs/theory_freeze.md`, 2026-08-26) was
made *before* any comparison with measurement, tested against Ho *et al.* 2026
Fig. 5(b), and **failed on one of its two branches**. This one is made *after*
that test, which makes it weaker as a pre-registration, and the document is
written so that the weakness is visible rather than hidden:

- every claim is stated as a **range over all four model forms that fit the
  yield data**, not as a single number from the one that was fitted;
- the claims are sorted into what survives that spread, what does not, and what
  is not claimed at all;
- the headline claim is a **theorem about the measurement**, not a predicted
  number, and therefore does not depend on the part of the model the comparison
  called into question.

Reproduce with `python repro_yield.py`, `python experiment_plan.py`,
`python -m pytest tests/ -q` (107 tests).

---

## What v1 got wrong, and how it was caught

v1 predicted λ_opt(120 GPa) = 475.5 nm and recommended 473 nm. Compared with 61
measured points of integrated PL yield versus pressure at 532 and 457 nm, the
532 nm branch was reproduced to 13 % and **the 457 nm branch was not**: v1 has
no yield maximum below 120 GPa, the data peak at 82 GPa.

The defect was the anchor, not the structure. Ho *et al.* write "a pronounced
shift **exceeding** 400 meV at 120 GPa"; v1 took the lower edge of a one-sided
statement as a point estimate and carried it with a ±20 meV bar. **The
mis-specified error bar was the real error**, and ħω's ±15 % has the same
problem: the value that fits is 103 meV against an assumed 65 ± 10.

Four different single inputs each reconcile both branches on their own, and all
four put λ_opt(120 GPa) between 443 and 455 nm. The displacement of the optimum
is robust; the attribution is not. Full account, with caveats:
`docs/sanity_check_pl_yield.md`.

---

## Tier 1 — what holds in every model form that fits

These are the claims the talk and the paper should rest on. Each is asserted by
`code/tests/`.

### 1. The brightest line is the most sensitive line

η ∝ Δν / (C √R), and neither Δν nor C carries wavelength dependence, so
minimising η **is** maximising the detected count rate. At equal optical power
and below saturation:

> **argmax_λ R(λ) = λ_opt, to under 1 nm in every model form.**

This is why λ_opt is measurable rather than predictable, and it is the one
claim that is untouched by the failure above. It **fails above the saturation
knee**: at u = 0.3 the brightest line sits 13 nm red of the best one, costing
×1.11. That failure is itself a prediction.

### 2. 457 nm is the line to load

| line | η(λ)/η_opt at 120 GPa, 300 K | at 77 K |
|---|---|---|
| 405 nm | ×1.40–2.37 | |
| **445 nm** | **×1.00–1.04** | |
| **457 nm** | **×1.00–1.05** | **×1.00–1.06** |
| 473 nm | ×1.08–1.25 | ×1.08–1.30 |
| 488 nm | ×1.24–1.81 | |
| 505 nm | ×1.56–3.11 | |
| 532 nm | ×2.44–11.71 | ×3.5 to unbounded |

457 nm is within ×1.05 of optimal under every surviving form, and within ×1.11
even under the rejected v1. **Load 457 nm. 473 nm is usable, not preferred.**

### 3. Blue beats green at 120 GPa, by a lot, and crosses over well below it

The green/blue crossing is at **54–67 GPa** across the surviving forms (v1 said
84 GPa). At 120 GPa green costs at least ×2.4 in η — a factor of 5.9 in
measurement time — and possibly far more.

### 4. The optimum moves with pressure, and nobody has followed it

| P | λ_opt range |
|---|---|
| 40 GPa | 509–527 nm |
| 60 | 490–503 |
| 100 | 457–466 |
| 120 | **443–455** |

Ho *et al.* used 532 nm to 51 GPa and 457 nm from 55 GPa — a two-step choice.
457 nm costs ×1.21–2.75 at 55 GPa and only becomes the right line above about
100 GPa. **Exactly one pressure in the published data (51.0 GPa) carries both
lines, and the two series are separately normalised there**, so no pressure
anywhere in the literature compares two excitation wavelengths on a common
scale. σ_abs(λ) under pressure has been calculated, never measured.

### 5. The optimum does not move with temperature

λ_opt(4 K) − λ_opt(300 K) = **+0.4 to +0.6 nm** in every form, against a
12 nm form-spread and a ~20 nm window. Temperature enters only through the
Franck–Condon envelope, which broadens symmetrically. Cooling *strengthens* the
blue case: 457 nm stays within ×1.06 while green gets worse.

### 6. The optimal wavelength moves blue with excitation power

The sensitivity ridge runs **25–30 nm blue between u → 0 and u = 0.3** in every
form. Sign and magnitude survive; the calibration-free recipe ratio (u*/u_knee)
does **not** — it is 1.8–2.2 across the surviving forms and 5.1 under v1.

---

## Tier 2 — central model, quoted with its spread

Central model: `NVModel(dE120=0.550)`. It is singled out only because ΔE_ZPL is
a *measured input* and re-anchoring it is the sanctioned exception to a freeze;
the other three forms change model-form parameters and are equally good fits.

| quantity | central | surviving spread |
|---|---|---|
| λ_opt(120 GPa) | 449.6 nm | 443–455 nm |
| 5 % window | 438–460 nm | 429–468 nm |
| ZPL at 120 GPa | 497 nm | 497–529 nm |
| η(532)/η_opt, 300 K | ×11.7 | ×2.4–11.7 |
| green/blue crossing | 67 GPa | 54–67 GPa |
| power ridge, u → 0 / 0.1 / 0.3 | 450 / 441 / 424 nm | — |
| u_knee / u* at 457 nm | 0.063 / 0.116 | ratio 1.8–2.2 |

**A prediction worth singling out.** In the forms with the larger ZPL shift, the
ZPL at 120 GPa lands at 497–504 nm — *red* of 532 nm. If that is right, green
excitation at 120 GPa is a purely thermally-activated anti-Stokes process and
**collapses on cooling** rather than degrading by ×2.5. One measurement decides
which: the ZPL position at 120 GPa against 532 nm.

---

## Tier 3 — explicitly not claimed

- **The size of the win over green.** ×2.4 to ×11.7 is not a result. Only the
  A/B measurement settles it.
- **Which input was wrong.** ΔE_ZPL, ħω, the S_abs slope and α are not separable
  by the available data.
- **That η behaves as modelled at all.** Nothing in the literature measures ODMR
  sensitivity versus excitation wavelength under pressure. σ_abs is one factor
  of it; f₋ flatness across the window is untested and untestable by the yield
  data.
- **λ_opt at 120 GPa as a computed number.** The data stop at 113.8 GPa and
  sample σ_abs at two photon energies; the forms spread ~30 nm at 120 GPa.
  Measuring the ZPL does not close this — the ZPL fixes one end of the
  absorption band and S·ħω the other, and S·ħω is what the forms disagree about.
- **Absolute power in mW**, MW power, modulation depth. By construction.
- **The (111) geometry.** All anchors are near-hydrostatic micropillar (100).

---

## What tests what

| claim | tested by |
|---|---|
| argmax R = λ_opt below saturation | step ①, and step ② by violating it above the knee |
| λ_opt(120 GPa) itself | **step ① — the five-line excitation scan.** ±2–5 nm at 5 % photometry, against ~30 nm from theory |
| the 5 % window width | step ①, from the curvature of the same scan |
| the power ridge and the recipe ratio | step ② — the 457 nm power sweep |
| the size of the win over green | step ③ — A/B at the operating point |
| the temperature invariance | step ③ repeated warm and cold |
| ZPL red of 532 nm ⇒ green collapses cold | the PL spectrum taken during step ① |
| the 54–67 GPa crossing, the ×5.9+ swing | step ④ — a 40 → 120 GPa sweep with both lines |
| f₋ flatness across the window | **nothing.** No measurement in this protocol touches it |

Steps ① – ④ are `docs/experiment_120GPa_laser_parameters.md`.

---

## Falsification

What would show the framework wrong, rather than a constant mis-set:

- **The excitation scan has no interior maximum** in 445–505 nm at 120 GPa —
  the band is not where any surviving form puts it, by more than the whole
  spread.
- **The brightest line is not the most sensitive line** below the saturation
  knee. This would break the η ∝ Δν/(C√R) reduction, i.e. show Δν or C carrying
  wavelength dependence, which is the load-bearing assumption of the whole
  argument.
- **The A/B at 120 GPa shows green ≥ blue.** Every form says blue wins by ×2.4
  or more; no admissible parameter produces a tie.
- **The power ridge does not move blue** as intensity rises.

What would *not* falsify it, and must not be presented as doing so: the scan
peak landing anywhere in 429–468 nm, the green advantage coming out anywhere in
×2.4–11.7, or the crossing pressure landing anywhere in 54–67 GPa. Those are the
freeze's own stated ignorance.

---

## What "frozen" means, and what would unfreeze it

Unchanged from v1: bug fixes only; extensions on a branch; re-anchoring on a
**measured** input is the sanctioned exception, and `experiment_plan.py` does
that without touching `nv_model.py`.

What would legitimately unfreeze it:

- **The excitation scan lands outside 429–468 nm.** No single-input re-anchoring
  covers that; the line shape itself would be wrong. The single-effective-phonon
  Pekarian with no Jahn–Teller treatment is the known caricature — Ho *et al.*
  compute mode-resolved lineshapes with explicit dynamical JT, and their curve
  fits the yield data where ours needed refitting.
- **The measured line shape disagrees with the fitted one** even when the peak
  agrees. Width and asymmetry are what the JT treatment changes.
- **Above ~100 GPa the measured ZPL is split** by residual non-hydrostatic
  stress; `nv_model.py` has one ZPL, and even Ho's theory models only the lower
  branch. On a (111) anvil this is the expected regime, not an edge case.

`code/tests/test_freeze.py` still enforces the v1 numbers as the record of what
was predicted before the comparison. It is not a statement that those numbers
are right — `docs/sanity_check_pl_yield.md` shows they are not — and it must not
be edited to agree with v2.
