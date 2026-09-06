# Ultra-high-pressure ODMR theory

[![tests](https://github.com/Risei412/ultra-high-pressure-ODMR/actions/workflows/ci.yml/badge.svg)](https://github.com/Risei412/ultra-high-pressure-ODMR/actions/workflows/ci.yml)

Theory, numerical analysis, and manuscript sources for excitation-wavelength
optimization in nitrogen-vacancy ODMR under pressure.

## Canonical files

- Manuscript: `paper/main.tex`
- Bibliography: `paper/refs.bib`
- Manuscript build: run `build.cmd` from `paper/`
- Numerical code and tests: `code/`
- Frozen theory record: `docs/theory/theory_freeze_v3_ho_integrated.md`

Files under `archive/` are historical snapshots. They are retained for
provenance and should not be used as the current implementation.

## Reproduce

The numerical record is frozen. `code/tests/test_freeze.py` and
`code/tests/test_freeze_a4_handoff.py` assert published and derived values to a
stated tolerance, and `code/repro_literature.py` checks the frozen model against
six published sources, with the two fitted calibration anchors marked `[CAL]` so
they are never mistaken for independent predictions.

Install the exact environment those values were verified in, and run the tests
from `code/`:

```powershell
python -m pip install -r requirements.lock
python -m pytest . -q
```

`requirements.lock` pins every package, including the test runner, to the
versions the frozen values were checked against; `requirements.txt` holds the
looser development constraints. Python version: see `.python-version`.

The full suite is 357 tests and takes about four minutes, roughly half of it the
one level-set scan described in `code/tests/conftest.py`. For a fast check of the
frozen numerical record alone (90 tests, about 10 seconds):

```powershell
python -m pytest tests/test_freeze.py tests/test_freeze_a4_handoff.py -q
```

Both are run in a clean pinned environment by CI on every change; the badge above
reports the result.

Generate the four manuscript figures from `code/`:

```powershell
python fig_m_level_set_ladder.py
python fig_g_gauge_degeneracy.py
python fig_x_branch_exchange.py
python esa_figv_kernel_sanity.py
```

The scripts write submission PDFs to `paper/figures/` and regenerable PNG
previews to `code/outputs/figures/`. Then build the manuscript:

```powershell
cd ..\paper
.\build.cmd
```

## Layout

- `paper/` — current manuscript, bibliography, submitted figure PDFs, build
  script, and ignored LaTeX output.
- `code/` — models, analyses, figure generators, digitized data, and tests.
- `docs/` — theory derivations, audits, experimental plans, manuscript notes,
  and reference material.
- `archive/` — immutable historical bundles, superseded manuscripts, and past
  review records.
- `slides/` — presentation material.
