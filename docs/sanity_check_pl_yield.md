# The sanity check, and what it did to the freeze

## Why this comparison and not the other one

`fig7b_contrast.png` puts measured ODMR contrast against modelled contrast.
It is a real check, but it is a check of the **absolute sensitivity scale**,
not of the wavelength recommendation: C is wavelength independent in this model
by construction, and `test_freeze.py` asserts that λ_opt does not move when
`C_amb` or `E_isc` are changed. A figure that cannot move the answer cannot
test it.

λ_opt is made of σ_abs(λ, P) and almost nothing else — the charge-state
fraction is flat to 0.2 % across the blue window. So the comparison that tests
it has to be a measurement of σ_abs(λ, P). That measurement exists: **Ho *et
al.* 2026, Fig. 5(b)** — integrated photoluminescence yield versus pressure at
a fixed excitation line, at constant laser power and camera integration time.
532 nm covers 4.7–51 GPa and 457 nm covers 51–114 GPa; 61 points in total.

As the ZPL blue shifts, the absorption band sweeps **past** the fixed laser
line, so the yield rises, peaks and falls. The pressure at which it peaks is a
direct read-out of where the absorption maximum is, and — this is the point —
it does not depend on the arbitrary vertical scale of the measurement. Each
series is allowed its own free scale factor and the test still has teeth.

Data: `code/data/pl_yield_vs_pressure.csv`, recovered from the vector content
of the figure rather than eyeballed, so the coordinates are exact to the
plotting precision. Analysis: `code/repro_yield.py`. Figure:
`slides/figs/fig9_yield.png`.

## The result

```
                            532 nm    457 nm    pooled
  frozen  dE120=0.400        13.5%     35.6%     26.4%
  refit   dE120=0.550        13.0%     13.5%     13.2%

  pressure of maximum yield (GPa)     532 nm    457 nm
  frozen                                  28       120
  refit                                   26        83
  measured                                26        82
```

**The green branch is reproduced. The blue branch is not.** The frozen model's
yield at 457 nm rises monotonically to 120 GPa; the measurement peaks at
82 GPa and falls. That is not a scale error or a noise question — it is a
qualitative feature the frozen model does not have.

Re-fitting the one anchor the comparison constrains, the ZPL shift at 120 GPa,
fixes both branches at once: 400 → **550 meV**, pooled residual 26 % → 13 %,
and the 457 nm peak lands at 83 GPa against a measured 82. The fit has a single
clean minimum and gives the same answer whether the model is evaluated at the
measurement temperature (90 K) or at 300 K.

Dropping `eta_col` — the fraction of emission still inside the fixed 650–800 nm
passband, which falls steeply as the emission band blue shifts — makes both
branches worse, so it belongs in the comparison. It is wavelength independent
and cannot move λ_opt.

## What it costs

| at 120 GPa | frozen | re-anchored |
|---|---|---|
| λ_opt | 475.5 nm | **449.6 nm** |
| 5 % window | 463–488 nm | **438–460 nm** |
| best commercial line | 473 nm ×1.00 | **445 nm ×1.01, 457 nm ×1.02** |
| 473 nm | ×1.00 | ×1.25 |
| η(532)/η_opt | ×2.5 | **×11.7** |
| 532/473 crossover | 73 GPa | 60 GPa |
| ZPL at 120 GPa | 529 nm | 497 nm |

The recommendation moves from 473 nm to **457 nm**, and the advantage over
green gets much larger, not smaller. The internal consistency is worth noting:
the re-anchored ZPL of 497 nm is outside the 512–541 nm range over which
`experiment_plan.py` says 473 nm holds, which is exactly the go/no-go that
step ① was built to apply.

## How much of this to believe

### The conclusion is robust; the attribution is not

The obvious worry is that 550 meV is just a fitted number and the whole thing is
one knob turned until it agreed. It is worth checking which knobs *can* do it,
because that is the difference between a wrong constant and a wrong theory.
Fitting each input on its own, against both branches at once:

| refit | pooled residual | λ_opt(120 GPa) |
|---|---|---|
| nothing (frozen) | 26.4 % | 475.5 nm |
| ΔE_ZPL(120) → 550 meV | 13.2 % | 449.6 nm |
| ħω → 103 meV | 9.8 % | 450.7 nm |
| S_abs slope → 4.55 | 12.3 % | 442.9 nm |
| α → 1.21 | 14.4 % | 455.2 nm |
| dE_ZPL/dP\|₀ → 8.3 meV/GPa | 23.0 % | 475.5 nm (fails) |

**Four different inputs each reconcile both branches, and every one of them
puts λ_opt between 443 and 455 nm.** So "the optimum is near 450 nm, not
475 nm" does not depend on which knob is blamed; "the ZPL shift is 550 meV"
very much does. The best commercial line is 445 or 457 nm in all four cases,
and the 473 nm penalty is ×1.08–×1.25 — the recommendation degrades, it does
not collapse.

What is **not** pinned by these data is the size of the advantage over green:
η(532)/η_opt comes out ×2.4 to ×11.7 depending on the knob. The ×11.7 quoted
below is the dE120 case alone and should not be used as a result.

### Why this is a wrong constant and not a wrong structure

Three things say so.

- **No trade-off.** Re-anchoring to fix the blue branch made the green branch
  *better* (13.5 % → 13.0 %), not worse. A model with the wrong structure buys
  agreement on one branch by giving it up on the other.
- **The mechanism is what was confirmed.** The test only exists because the
  absorption band sweeps *past* a fixed line as the ZPL shifts, which is the
  same statement as "λ_opt is set by the absorption line shape alone". The data
  show both domes, in the right order, with the green one at the predicted
  pressure. The layer of the argument that the talk is about passed.
- **400 meV was a bound, not a value.** Ho *et al.* write "a pronounced shift
  **exceeding** 400 meV at 120 GPa". The freeze took the lower edge of a
  one-sided statement and used it as a point estimate with a ±20 meV bar. That
  is a defect in how the anchor was read, not in the physics.

### Where the model form really is a caricature

The paper's own Fig. 5(b) theory curve reproduces both domes, and its σ_abs is
a mode-resolved calculation with explicit dynamical Jahn–Teller treatment of
the JT-active modes (the K² constants). `nv_model.py` replaces all of that with
a single 65 meV effective phonon and no JT. The paper notes that absorption is
*not* the mirror of emission precisely because of the JT effect, and that the
sideband broadening under pressure "is most pronounced in absorption" — which
is exactly the regime this comparison probes. The ħω row above, which gives the
best fit of the four, is that story: the effective phonon needs to be 103 meV,
far outside the ±15 % that was carried as the model-form uncertainty.

A second structural simplification: above ~100 GPa the measured ZPL is
**split** by residual non-hydrostatic stress, and even the paper's theory
models only the lower branch. `nv_model.py` has one ZPL.

### What genuinely broke: the uncertainty budget

The tornado figure carries ΔE_ZPL(120) at ±20 meV, a measurement-precision
number for a quantity that was quoted as a bound. The data displace it by
150 meV. The central value being wrong is recoverable; claiming ±20 meV on it
was the actual error, and the same question should be asked of ħω, whose ±15 %
also fails to cover the value that fits.

Other caveats:

- The measurements are on a near-hydrostatic micropillar sample at ≈90 K. The
  planned run is a (111)-cut anvil. Step ① re-anchors on the measured ZPL
  either way, which is what makes the protocol robust to exactly this.
- The two series are separately normalised in the source figure, so nothing
  here tests the *ratio* between green and blue at a given pressure. Only the
  shape of each dome is tested. The ×11.7 above is a model output, not a
  measurement.
- Sample degradation over a long pressure run would also depress the high-
  pressure end of the 457 nm series. It would not move the peak by 40 GPa.

## What is decided, and what the theory is still needed for

Run the same five models — the frozen one and the four single-knob refits, each
of which fits the data — and ask what they agree on.

**Penalty η(λ)/η_opt at 120 GPa, across all five:**

| line | 405 | 445 | 457 | 473 | 488 | 505 | 532 |
|---|---|---|---|---|---|---|---|
| range | ×1.4–3.8 | **×1.00–1.34** | **×1.00–1.11** | ×1.00–1.25 | ×1.05–1.81 | ×1.33–3.11 | ×2.4–11.7 |

**457 nm is within ×1.11 of optimal under every model considered, including the
one the data rejected.** So as an engineering decision the wavelength question
*is* closed: use 457 nm (or 445), and no plausible revision of the model moves
that. The pressure at which blue overtakes green is 54–84 GPa in all five, so
"blue wins at 120 GPa" is robust too; "it crosses at 73 GPa" is not.

What is *not* decided, and what no amount of further modelling will decide:

- **The size of the win.** η(532)/η_opt spans ×2.4 to ×11.7. Only an A/B
  measurement in the actual apparatus settles it — step ③.
- **Whether η behaves as modelled at all.** Nothing here measures ODMR
  sensitivity versus excitation wavelength under pressure. That measurement
  does not exist in the literature. This comparison tests σ_abs, which is one
  factor of it.
- **f₋ flatness across the window**, which is the reason the answer reduces to
  the absorption shape. Untested, and untestable by this figure.
- **The power dependence.** The ridge still runs ~25–30 nm blue between u → 0
  and u = 0.3 in all five models, so the effect is robust in sign and rough
  size; the calibration-free recipe ratio is not (1.8–5.1 depending on knob and
  line). The procedure survives, the number does not — step ② measures it.

The honest summary is that this comparison closed the question the model was
*least* needed for and left open the ones it was built for.

## What this means for the freeze

The freeze survives in form and fails in content. Re-anchoring on a measured
input is the sanctioned exception (`docs/theory_freeze.md`), and this is that:
one anchor, fitted to data, no change to `nv_model.py`. But the headline number
the deck is built on moves 26 nm, and the recommended line changes.

The honest reading is that the pre-registration worked. A prediction was frozen,
a comparison was made that could kill it, and it was killed on one of its two
branches. That is a better talk than one where everything agreed.

**Nothing in this changes the structure of the argument** — λ_opt is still set
by the absorption line shape alone, the optimum is still a window rather than a
point, and the protocol still measures the ZPL first and reads the line off it.
What changes is which line comes out of it.

---

## The gap this leaves: nobody has swept the wavelength

It is worth being precise about what the state of the art did, because it is
less than it looks.

Ho *et al.* used **532 nm from 4.7 to 51 GPa and 457 nm from 55 to 114 GPa**.
That is a two-step function, chosen once, not a scan. Consequences:

- **Exactly one pressure in the whole dataset carries both lines** — 51.0 GPa —
  and the two series are separately normalised in the figure, so even that one
  is unusable as a ratio. **The literature contains no pressure at which two
  excitation wavelengths were compared on a common scale.**
- Each dome therefore samples σ_abs at a *single* photon energy per pressure.
  Two scale-free numbers come out of 61 points — the two peak pressures — and
  that is why the model forms that all fit still spread ~30 nm at 120 GPa.
- **σ_abs(λ) under pressure has never been measured.** It has been calculated,
  and the calculation has been checked at two energies.

And the fixed choice was not the optimum. Across every model form that fits,
evaluated at the measurement temperature:

| P (GPa) | λ_opt range | η(457)/η_opt |
|---|---|---|
| 55.2 | 494–518 nm | ×1.21–2.75 |
| 74.5 | 478–502 | ×1.07–1.89 |
| 85.8 | 469–495 | ×1.03–1.59 |
| 99.1 | 458–487 | ×1.00–1.35 |
| 113.8 | 447–479 | ×1.00–1.19 |

457 nm only becomes the right line above about 100 GPa. At the bottom of its
own blue run the state of the art was leaving up to ×2.75 in η — a factor of
**7.6 in measurement time** — on the table. Symmetrically, 532 nm was fine to
about 41 GPa and cost ×1.04–1.27 by 51 GPa.

`slides/figs/fig10_tracking.png` draws this: λ_opt(P) as a band over the
surviving model forms, with the two fixed lines and the pressures they were
used at.

So the earlier statement that "the answer is already out" was too pessimistic
and should be narrowed. What is out is a *calculated* absorption lineshape from
which a number could be read. What is not out, anywhere:

1. any **measurement** of the excitation-wavelength dependence at pressure;
2. the statement in terms of **sensitivity** rather than brightness — which is
   what makes the measurement worth making, and which is a theorem here, not an
   assumption (below saturation, argmax R = λ_opt, to under 1 nm in every model
   form);
3. **λ_opt(P) as a law** rather than a choice made twice;
4. the **power dependence**, which moves the ridge 25–30 nm by u = 0.3 in every
   surviving form.

A single excitation scan at one high pressure supplies (1) and tests the
calculated lineshape peak — the quantity this whole re-analysis showed to be
the fragile one — for the first time.
