"""
conftest.py
-----------
One shared kernel for the whole test session.

WHY THIS EXISTS.  theory_a2_multiplicity.observed_ladder scans 900 power
ratios and then bisects 40 times per transition, each step evaluating a level
set on a 0.002 nm wavelength grid.  That is about 1100 level-set evaluations
and it costs roughly 95 s.  The module caches the result -- but keyed on
``id(kernel)``:

    key = (id(kernel), tuple(window), max_ratio, samples)

tests/test_theory_a2_multiplicity.py and tests/test_theory_figures.py each
used to define their own ``@pytest.fixture(scope='module') def kernel()``.
Module scope means two distinct Kernel objects, two distinct ids, and a cache
miss, so the 95 s scan ran TWICE.  Those two tests were 192 s of a 299 s
suite.

Hoisting the fixture here at session scope makes both files share one Kernel,
so the scan is paid once.  Nothing in the suite mutates the kernel -- it wraps
the digitised Ho spectra and exposes only readers -- so sharing is safe.

Holding the instance for the whole session also removes a latent hazard in the
``id()`` key: CPython reuses ids after garbage collection, so a Kernel created
after another was freed could collide with the freed one's cache entry and get
its results.  A session-scoped fixture never dies, so the collision cannot
arise for the suite.  The hazard is still there for library callers who create
and drop kernels in a loop; that is noted in theory_a2_multiplicity, not fixed
here, because fixing it means changing a frozen module's cache semantics.

DO NOT re-add a module-scoped ``kernel`` fixture to an individual test file.
pytest resolves the nearest definition, so a local one silently shadows this
and the 95 s comes straight back.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from theory_a1_generalization import Kernel  # noqa: E402


@pytest.fixture(scope='session')
def kernel():
    """The reconstructed Ho kernel, built once per test session."""
    return Kernel()
