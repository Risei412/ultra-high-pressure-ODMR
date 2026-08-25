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

## Type size

Nothing on a slide is below 18 pt. `build_talk_deck.py` enforces the floor on
the deck (`raise_font_sizes`, which also covers table cells and end-paragraph
properties), and `talk_figs.py` enforces it on the figures by drawing each one
at exactly the size it occupies on the slide -- 7.70 in wide -- so an 18 pt
label in matplotlib is 18 pt when projected. `replace_picture` places the PNGs
1:1 for the same reason; rescaling them would quietly break the floor.

Citations are the one exemption, at 12 pt: nobody reads a reference list from
the back of a room, and dropping them frees the foot of the slide.

## Colour

Three values, nothing else: Institute of Science Tokyo navy `1C3177` for all
text, white `FFFFFF` where the table header reverses it, and
grey `6E7488` on the citation strip, which is not meant to compete. The figures
follow the same rule -- navy for what the slide is about, one neutral grey for
the option being argued against, and no third colour.

## Where text goes

On a figure slide the right-hand column had been narrating the plot. It now
carries only what the figure cannot show, in 50-90 Japanese characters; the
fact the figure establishes is
**one sentence** across the foot (`BOTTOM`), in navy, at the largest of
28/26/24/22/20 pt that still keeps it on one line -- which is what makes the
sentences short. The title already carries the headline, so the foot carries
the number. The citation sits at 12 pt beneath it. One sentence is the rule, not a target: anything the sentence
does not need in order to carry the claim -- a legend gloss, a value already
labelled on the plot, a definition the axis already gives -- comes out of the
foot. The grey strip is references and nothing else -- annotations that used
to hide there are said out loud instead. Raising the body to 18 pt costs room, so phrasing was
tightened throughout -- but no claim was dropped. The qualifiers that carry one
(`micropillar は前提条件`, `1 件は OPEN と申告`, the 50 GPa caveat on the
ambient absorption curve, the power penalty in the conditions table) are all
still on the slides.

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
