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

The fitted 550 meV is an **effective** ZPL shift inside this model form, not
independently a measurement of the ZPL. What the data actually say is that the
absorption maximum at 120 GPa is about 20 nm bluer than the frozen model puts
it. Within this parameterisation the only knob that does that is dE120; the
same displacement could come from the single-effective-phonon envelope being
the wrong shape, which is the ±15 % ħω model-form error already named as the
dominant uncertainty. The two are not separable by these data.

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
