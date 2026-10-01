# Merge resolution review (2026-10-01)

Scope: merge main (4e75f60) into integrate-local-work-20260906 (333ef89).

## User-selected resolution

The user selected the main excitation-design manuscript, with explicit dependence
of the optimal blue laser wavelength on alpha. paper/main.tex retains that structure.
The branch theory manuscript is preserved verbatim in
archive/manuscripts/main_theory_GMX_20260906.tex.

## Scientific integration

- Preserve the branch corrected vector extraction of Ho Fig. 1(b,c), matching
  scientific code and tests, stronger gauge checks, and executable peak counting.
- Preserve main figure-script layout and working-directory-independent data paths.
- Clarify alpha = sigma_xx/sigma_zz = sigma_yy/sigma_zz, axial stress 120 GPa,
  and temperature 300 K. Show the differential shift anchored at the published
  hydrostatic Ho optimum 440.65 nm, with both Hilberer normalization assumptions.
- Distinguish legacy single-mode 475.5 nm from the conditional hydrostatically
  anchored estimates. The interval between normalizations is not a confidence
  interval and anisotropic absorption has not been experimentally validated.
- Preserve frozen theory records, with integration notes about corrected vector
  inputs (DWF decline 6.01; K3 ratio 0.89; density growth 4.31).
- Update v1 diagnostic snapshots independently from corrected vector data,
  without changing algorithms, tolerances, or test gates.
- Correct the inherited conflated Barfuss2017 citation metadata to MacQuarrie,
  Otten, Gray, and Fuchs, Nat. Commun. 8, 14358 (2017), DOI
  https://doi.org/10.1038/ncomms14358. Publisher text explicitly reports the
  same 850 +/- 130 THz/strain coupling; the scientific value is unchanged.
- Keep canonical docs/experiment, docs/manuscript, and docs/theory locations.
  Obsolete root copies were backed up to temporary storage before removal.

## Verification

- Full integrated suite: 409 passed (188.12 s).
- Alpha optimum/tolerance figure generation succeeded; PNG visually checked.
- Standalone introductory slides generation succeeded.
- Final manuscript build and Git conflict-index checks are recorded below after
  completion. Draft PR #5 remains the review destination; main is not merged.

- Final paper/build.cmd: exit 0, 18-page PDF; no undefined citations or references.
- Three inherited table-width warnings remain in unrelated baseline tables.
- Working tree scan: no conflict markers; all 30 original conflicts resolved.
