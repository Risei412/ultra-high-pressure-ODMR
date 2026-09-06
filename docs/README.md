# Documentation

- `theory/` — derivations, numerical execution records, validation notes, and
  the frozen theory record.
- `audits/` — novelty and journal-positioning audits.
- `experiment/` — model-development and future experimental plans.
- `manuscript/` — notes used to develop the manuscript narrative.
- `references/notes/` — literature notes.
- `references/papers/` — locally retained reference PDFs.

Historical review reports are under `../archive/reviews/`.

## Where the limits are recorded

- `experiment/theory_limits_and_required_measurements.md` is the single
  statement of what the theory does and does not determine. It covers two
  axes: anisotropic stress at 120 GPa (the α axis) and extrapolation above
  120 GPa (the pressure axis, added 2026-09-02).
- `references/notes/audit_2026-09-02_six_additions.md` records which published
  papers were checked for pressure-resolved optical spectroscopy and what they
  contain. Six of seven contain none.

Two rules follow from those documents and are repeated here because they are
easy to violate by accident:

1. **Do not quote a single optimal excitation wavelength above 120 GPa.** Every
   optical anchor stops there; `ho_spectrum_model.py` refuses past it and
   `nv_model.py` returns an asymptote under its default clip policy. Quote a
   band and a worst-case penalty instead — `code/extrapolation_bounds.py`
   produces both.
2. **The largest single uncertainty is not the extrapolation.** The
   reconstructed kernel and the single-mode envelope disagree by 34.9 nm at
   120 GPa, inside the fully anchored range. Resolving that needs no new
   experiment and no new literature, and comes before any extrapolation work.
