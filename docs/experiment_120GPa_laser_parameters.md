# Laser parameters for the 120 GPa run

> **2026-08-27 — step ① has been rewritten; read `docs/sanity_check_pl_yield.md`
> first.** Measured PL yield versus pressure puts the absorption maximum about
> 20 nm bluer at 120 GPa than the anchor these tables were built on, and the
> model forms that fit those data still disagree by ~30 nm about λ_opt at
> 120 GPa. Reading the line off a measured ZPL does not resolve that. Step ① is
> therefore now a **direct measurement of λ_opt in the cell** — a five-line
> excitation scan — which is possible because the brightest line is the most
> sensitive line below saturation. Load **457 nm** as the working line: it is
> within ×1.11 of optimal under every model form that fits, including the
> rejected one. The ZPL route is kept below as a fallback.

Two numbers have to be fixed before the cell is loaded: **which line** and
**how much power**. Both are decided by measurements made in the cell itself,
not by the frozen model alone — the model supplies the rule, the cell supplies
the number.

Everything below is produced by `code/experiment_plan.py`:

```bash
cd code
python experiment_plan.py          # both tables, at 120 GPa
python experiment_plan.py 100      # at another pressure
python -m pytest tests/test_experiment_plan.py -q
```

---

## Why step 1 exists: the (111) anvil is outside the frozen anchors

The anchors — ΔE_ZPL(120 GPa) = 400 meV, α = 0.95 — were measured on a
micropillar (100) culet under near-hydrostatic stress. A **(111)-cut anvil
selects the NV group aligned with the compression axis, and that group's ³E
red-shifts** rather than blue-shifts (Davies & Hamer 1976; Bhattacharyya thesis
2022). `paper/main.tex` states plainly that the anchors would have to be
re-derived for that geometry, and `nv_model.py` does not model it.

So **473 nm does not transfer automatically**. It has to be re-derived — but
only from one measurement, because λ_opt is set by the absorption line shape
and the only part of that shape the geometry moves is the ZPL position.

There is a bonus: Huang *et al.* (arXiv:2511.20750) find *ab initio* that
near-uniaxial stress along [111] is not merely tolerable but **genuinely
optimal for the contrast** of that NV subensemble. The cost is that only one of
the four orientations is usable, which rules out vector magnetometry but is
harmless — arguably favourable — for flux exclusion along the compression axis.

---

## Step 1 — the excitation wavelength, measured at 120 GPa

**The theory does not have to supply this number, and should not be asked to.**

Ho *et al.*'s yield data stop at 113.8 GPa and probe σ_abs at only two photon
energies. Getting λ_opt(120 GPa) out of them means extrapolating the absorption
peak past the last measurement, and `docs/sanity_check_pl_yield.md` shows what
that costs: the model forms that all fit those data disagree by ~30 nm about
λ_opt at 120 GPa. Measuring the ZPL does **not** rescue it either — the ZPL
fixes one end of the band and the Stokes shift S·ħω fixes the other, and it is
S·ħω that the forms disagree about.

The way out is that λ_opt is directly measurable in the cell, in an afternoon.

### Why it is measurable

η ∝ Δν / (C √R), and neither Δν nor C carries any wavelength dependence. So
minimising η **is** maximising the detected count rate:

> **At equal optical power and below saturation, the brightest line is the most
> sensitive line.**

This holds identically in all five model forms — `test_experiment_plan.py`
asserts it to under 1 nm — so it does not depend on the part of the theory the
yield comparison called into question. The theory's job here is to have proved
the equivalence and to say how to use it, not to predict the number.

### The measurement

1. At the working pressure, illuminate at **445, 457, 473, 488 and 505 nm**
   (three is the minimum; the peak must be bracketed).
2. Set each line to the **same optical power**, well below the saturation knee.
3. Record the detected PL through the usual filter, same integration time.
4. Divide by the delivered power and fit the peak.

```python
from experiment_plan import excitation_scan, SCAN_LINES
excitation_scan(power_normalised_counts, SCAN_LINES)   # -> lambda_opt, nm
```

At 5 % photometry this returns the optimum to **±2–5 nm** whichever model form
is true — against the ~30 nm spread the theory alone leaves. Nothing about the
absorption line shape has to be assumed; only that it has a maximum.

### The one thing the theory must be believed on

**Stay below the knee.** The equivalence is exact at u ≤ 0.06 and fails above
it: at u = 0.3 the brightest line sits 13–24 nm red of the most sensitive one,
because saturation flattens the bright end of the scan. Run the scan at low
power, then set the power by step 2. If the scan has to be run hot, the peak it
returns is not λ_opt.

### Go/no-go

The lines to bring are settled without any of this: **457 nm is within ×1.11 of
optimal under every model form that fits the data, including the rejected one**,
and 445 nm is within ×1.34. So load 457 nm as the working line and let the scan
confirm or refine it. 473 nm is ×1.00–1.25 — usable, not preferred.

---

## Step 2 — the excitation power, from the saturation curve

The model's intensity is `u = I/I_half`, and `I_half` is exactly what is not
calibrated, so **no absolute power can be quoted**. It does not have to be:
the knee of the saturation curve and the sensitivity optimum both sit at fixed
values of `u`, so their **ratio** needs no calibration.

At 473 nm and 120 GPa:

| | u |
|---|---|
| knee of R(u), 10% below the linear slope | 0.060 |
| sensitivity optimum | 0.114 |
| **ratio** | **1.9** |

**Measure:** sweep the 473 nm power, record photoluminescence, find the power
at which it falls 10% below the low-power linear slope.
**Then operate at 1.9 × that power.** Units cancel; absolute intensity is never
needed.

```python
from experiment_plan import analyse_power_sweep
r = analyse_power_sweep(power_mW, pl_counts, lam=473.0, P=120.0)
r['knee_power'], r['operating_power']
```

η is flat to 2% over u = 0.10–0.15, so the setting does not need to be precise.

**Guard rails** — past the optimum, two losses stack:

| u | η/η* | λ_opt(u) |
|---|---|---|
| 0.03 | 1.40 | 475 nm |
| 0.06 | 1.10 | 475 |
| **0.114** | **1.00** | **~465** |
| 0.20 | 1.09 | 451 |
| 0.30 | 1.25 | 447 |
| 0.50 | 1.53 | 442 |

Above u ≈ 0.3 the sensitivity is 25% worse **and** the optimal wavelength has
moved 28 nm blue of the line you are using. Note that even at the optimum the
ridge has already left 475 nm; the penalty for staying at 473 nm there is only
×1.006, which is why the recommendation survives.

**Take these at the same time**, because they cost nothing extra and they fix
the two remaining uncalibrated constants:

- ODMR contrast **C(I)** → fixes Γ_sat
- ODMR linewidth **Δν(I)** → fixes Γ_c

With those, the power-dependence result stops being a pre-calibration scenario.

---

## Cooling does not change the answer

The run will be cold, not at room temperature. λ_opt moves **0.6 nm between
4 K and 300 K**, so nothing in steps 1 and 2 has to be redone on cooling: pick
the line warm, and it is still the right line cold.

The margin gets better, not worse. The sideband narrows as the anti-Stokes wing
freezes out, so 532 nm on the far red flank falls further while the peak stays:

| T | η(532)/η_opt |
|---|---|
| 300 K | 2.53 |
| 150 K | 3.39 |
| 77 K | 4.10 |
| 4 K | 4.89 |

**The ×2.5 quoted at room temperature is conservative for the measurement that
will actually be made.** If step ③ is run both warm and cold, the growth of
that ratio is itself a prediction being tested.

## What this does not give

- **Absolute mW.** By construction. Step 2 returns a multiple of a measured
  power, never an intensity.
- **MW power and modulation depth.** Γ_c and Γ_sat are uncalibrated, so the
  model has no opinion; the framework is Dréau *et al.*, PRB **84**, 195204
  (2011). Step 2's sweep is what supplies them.
- **A check on ħω.** The ±15% model-form range on the single effective phonon
  energy is the largest single contributor to the uncertainty in λ_opt
  (±6.9 nm) and no measurement in this protocol touches it. Only a direct
  measurement of σ_abs(λ) at pressure would.

---

## Order of operations

1. **PL spectrum at 120 GPa** → ZPL → pick the line (this document, step 1).
2. **473 nm power sweep** → knee → operating power, plus C(I) and Δν(I)
   (step 2).
3. **A/B against 532 nm at that operating point** → the achieved gain in *this*
   apparatus, systematics included. This is the result the paper is about.
4. *Optional:* a 40 → 120 GPa sweep with both lines, which is the falsification
   test — the ratio must swing ×5.8 and cross unity. Piggyback it on a run that
   sweeps pressure anyway; it does not deserve a dedicated load.

Steps 1 and 2 replace every predicted condition with a measured one. Step 3 is
the measurement the theory was written to enable.
