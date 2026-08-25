# Figure panels reproduced from the literature

These are crops of published figures, used in `slides/odmr_474nm_talk.pptx`
with the citation printed on the slide.  They are NOT this work's results.

| file | source |
|---|---|
| `lit_dac_schematic.png` | P. Bhattacharyya, *Principle and applications of quantum metrology in diamond anvil cells using the Nitrogen Vacancy color center*, PhD thesis, UC Berkeley (2022), Fig. 2.2(a,b) — the anvil and a DAC loading. |
| `lit_nv_in_anvil.png` | ibid., Fig. 7.2(a) — the NV layer implanted 50 nm below the culet of the upper anvil, with the sample and transport leads below. |
| `lit_meissner_map.png` | ibid., Fig. 7.8(b,c) — confocal image of the CeH₉ sample and the field ratio *B/H* along the marked line cut; the same measurement campaign as P. Bhattacharyya *et al.*, Nature **627**, 73 (2024). |

The thesis PDF itself is in `docs/ref/`.  Regenerate the crops with:

```bash
pdftoppm -png -r 300 -f 52  -l 52  docs/ref/Principle_and_Applications_of_.pdf p52
# then crop (left, top, right, bottom) as fractions of the page:
#   Fig 2.2(a,b)  0.170 0.135 0.415 0.485
#   Fig 7.2(a)    0.222 0.283 0.383 0.424   (page 115)
#   Fig 7.8(b,c)  0.375 0.325 0.845 0.462   (page 123)
```
