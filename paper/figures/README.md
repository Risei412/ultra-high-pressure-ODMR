# Manuscript figures

The PDF files in this directory are the vector figures included by
`../main.tex`. Regenerate them by running the matching scripts in `../../code/`.
PNG previews are written to `../../code/outputs/figures/` and are not versioned.

| file | script | manuscript |
| --- | --- | --- |
| `fig_g_gauge_degeneracy.pdf` | `fig_g_gauge_degeneracy.py` | Theorem G |
| `fig_m_level_set_ladder.pdf` | `fig_m_level_set_ladder.py` | Theorem M |
| `fig_x_branch_exchange.pdf` | `fig_x_branch_exchange.py` | Theorem X |
| `esa_figv_kernel_cross.pdf` | `esa_figv_kernel_sanity.py` | App. B, cross-figure check |
| `esa_figv_kernel_internal.pdf` | `esa_figv_kernel_sanity.py` | App. B, checks K1–K4 |

All five are full-width `figure*` floats included at `\textwidth` ≈ 7.06 in,
and each is drawn at that width so the page rescales nothing. Font sizes come
from `../../code/fig_style.py`; do not override them per script.
