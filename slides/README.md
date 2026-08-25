# 8-minute internal talk

`odmr_473nm_talk.pptx` is the talk deck. It is generated, not hand-edited:

```bash
python talk_figs.py          # figures -> figs/
python build_talk_deck.py    # odmr_474nm_talk.pptx -> odmr_473nm_talk.pptx
```

`odmr_474nm_talk.pptx` is the earlier draft and is kept only as the layout
source that `build_talk_deck.py` transforms; edit the build script rather than
the output deck, or the next run will overwrite the change.

## What changed from the draft

The draft was written against the pre-revision numbers. Three things moved:

| | draft | now |
|---|---|---|
| mechanism | optimum = absorption maximum | optimum = maximum of `lambda * sigma_abs` (equal optical power) |
| optimum at 120 GPa | 474.0 nm | 475.5 nm, MC 16-84% `+5.6 / -5.3` nm |
| 5% window | 462-486 nm | 463-488 nm |
| 473 nm penalty | 0.03% | 0.2% |
| gain over 532 nm | x2.7 (map x7) | x2.5 (map x6) |
| green/blue crossover | 71 GPa | 73 GPa |
| largest uncertainty | ZPL shift | the single-effective-phonon approximation |

Because the Monte Carlo band now spans about half the tolerance window, the
claim is stated as a window with the commercial 473 nm line inside it, which is
what `paper/main.tex` says as well.

## Structure

| # | slide | figure |
|---|---|---|
| 1 | headline: measure at 473 nm | - |
| 2 | the sample sits between two anvils | literature schematic |
| 3 | NV inside the anvil images flux exclusion | literature data |
| 4 | both edges move blue at 120 GPa | `fig4_absorption.png` |
| 5 | only the optical inputs move the optimum | `fig5_tornado.png` |
| 6 | the answer is a window, 463-488 nm | `fig6_window.png` |
| 7 | reproduction, and the sign change at 73 GPa | `fig7_crossover.png` |
| 8 | the optimum walks blue with power - uncalibrated | `fig8_power.png` |
| 9 | where the recommendation holds and where it fails | - |

Slides 4-8 are regenerated from `code/nv_model.py` and
`code/nv_model_power.py`, so the numbers on them cannot drift from the
manuscript. Slides 2 and 3 carry scanned literature figures and are untouched.
