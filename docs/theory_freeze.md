# Theory freeze — 2026-08-26

The model is frozen here. Everything below is a **prediction made before the
measurement**, which is the only thing that makes the comparison worth
anything: a model improved after seeing the data can always be made to agree.

Frozen at commit `ab5b4f3` on branch `claude/repository-theory-presentation-rm3txv`.
(A `theory-freeze-2026-08-26` tag was cut locally but the remote refused the tag
ref, so the commit hash is the marker. Re-cut the tag from a machine with direct
push access if a tag is wanted.)

## Why freeze rather than develop

The largest single contributor to the uncertainty in λ_opt is ħω at ±6.9 nm,
and that is **model-form error** — the single-effective-phonon approximation —
not measurement error. More modelling cannot reduce it; only a direct
measurement of σ_abs(λ) at pressure can.

The one genuinely novel result, the power dependence, is explicitly
uncalibrated: `main.tex` states that its magnitude depends on a′_es, s_d, s_d0
and on using the NV⁻ Huang–Rhys factor for NV⁰. Extending the theory without
data means adding unconstrained parameters to a model whose existing
unconstrained parameters already dominate — the model gets more elaborate and
less falsifiable.

Meanwhile a single cheap measurement (a 473 nm power sweep at the working
pressure) collapses a whole class of unknowns and converts the strongest claim
from a scenario into a result.

## What "frozen" means

- **Bug fixes only.** A defect that makes the code disagree with its own
  documented model is a fix. Anything that changes the model is not.
- **Model extensions go on a branch** and are not merged before the experiment.
- **The one sanctioned exception** is re-anchoring on measured inputs, which is
  not a change to the model but the calibration the model asked for:
  `experiment_plan.py` re-anchors the ZPL on a measurement without touching
  `nv_model.py`.
- `code/tests/test_freeze.py` enforces the numbers. If it fails, either a bug
  was fixed and the golden values need updating with a stated reason, or the
  freeze was broken.

## The frozen predictions

At 120 GPa, room temperature, near-hydrostatic (α = 0.95) unless stated.

| quantity | prediction |
|---|---|
| λ_opt | **475.5 nm** |
| 5% tolerance window | **463–488 nm** |
| η(405)/η_opt | 3.76 |
| η(457)/η_opt | 1.11 |
| η(473)/η_opt | **1.002** |
| η(488)/η_opt | 1.05 |
| η(532)/η_opt | **2.53** → map time ×6.4 |
| f₋ across the window | 0.5835–0.5848 (0.2%) |
| ZPL at 120 GPa | 529 nm |
| IP(³A₂) edge | 405.2 nm |
| λ_opt at 100 GPa | 486 nm |
| green/blue crossover | **73 GPa** |
| η(532)/η(473) swing, 40 → 120 GPa | **×5.8**, crossing unity |
| power ridge, u → 0 / 0.1 / 0.3 | 475.5 / 465 / 447 nm |
| 473 nm penalty at u = 0.2 / 0.3 | ×1.22 / ×1.51 |
| ZPL range for which 473 nm holds | **512–541 nm** |
| power recipe | operate at **1.9 ×** the 10% knee (u: 0.060 → 0.114) |
| λ_opt(4 K) − λ_opt(300 K) | **+0.6 nm** |
| η(532)/η_opt at 4 / 77 / 150 / 300 K | **4.89 / 4.10 / 3.39 / 2.53** |
| 5% window at 4 K | 465–487 nm (300 K: 463–488) |

Reproduce with `python repro_literature.py` and `python experiment_plan.py`.

## Temperature

The Meissner measurements happen below T_c, not at room temperature, so the
recommendation has to survive being cooled. It does, and not by assumption:
temperature enters the model only through the Franck–Condon envelope, which
broadens symmetrically and leaves its peak where it is. **λ_opt moves 0.6 nm
between 4 K and 300 K**, against a ±5.5 nm uncertainty and a 25 nm window.

Cooling in fact *strengthens* the case. The anti-Stokes wing disappears, the
sideband narrows, and 532 nm — which sits on the far red flank — falls further
while the peak does not move: the green penalty rises from ×2.53 at 300 K to
×4.89 at 4 K. **The headline ×2.5 is the room-temperature value and is
conservative for the experiment that will actually be run.**

Two things the model does not carry here. The real NV⁻ ZPL shifts by roughly
15 meV over 4–300 K through thermal expansion and electron–phonon coupling,
which the model's ZPL does not; at −0.18 nm of λ_opt per meV that is about
3 nm, larger than the 0.6 nm above but still well inside the window. And C₀ is
temperature dependent in reality, which moves the absolute sensitivity — but
C₀ is wavelength independent and therefore cannot move λ_opt at all, which
`test_freeze.py` asserts.

## What each prediction is exposed to

| prediction | tested by |
|---|---|
| λ_opt, the window, the line penalties | step ③ — A/B against 532 nm at the operating point |
| the ZPL range for 473 nm | step ① — one PL spectrum |
| the power ridge and its penalties | step ② — the 473 nm power sweep |
| the 73 GPa crossover and the ×5.8 swing | step ④ — a 40 → 120 GPa sweep with both lines |
| the temperature invariance | step ③, if the A/B is repeated warm and cold |
| ħω, the dominant uncertainty | **nothing in this protocol.** Only a direct σ_abs(λ) measurement at pressure |

Steps ① – ④ are `docs/experiment_120GPa_laser_parameters.md`.

## Caveat carried into the experiment

The run is planned on a **(111)-cut anvil**, and the frozen anchors were
measured on a micropillar (100) culet. A (111) culet selects the NV group whose
³E red-shifts; the anchors do not transfer. This is not a reason to unfreeze —
step ① re-anchors on the measured ZPL, which absorbs whatever the stress state
did, and 473 nm survives for any ZPL between 512 and 541 nm against an anchor
of 529 nm.

What step ① cannot check is whether the uniaxial stress also changed S_abs or
ħω. If the measured ZPL lands in range but the observed line shape does not
match, that is the signal to unfreeze.
