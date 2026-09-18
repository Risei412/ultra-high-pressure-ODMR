"""Print the numerical environment this run is using.

Cross-environment agreement is reported at a stated tolerance rather than
bit-for-bit, because different CPUs and BLAS kernels legitimately differ in
the last bits. That claim is only checkable if each run says which CPU and
which BLAS it used, so CI records it here.
"""

import platform
import sys

import numpy as np


def main() -> None:
    print("python      :", sys.version.replace("\n", " "))
    print("platform    :", platform.platform())
    print("machine     :", platform.machine())
    print("processor   :", platform.processor() or "(not reported)")
    print("numpy       :", np.__version__)

    try:
        import scipy

        print("scipy       :", scipy.__version__)
    except ImportError:  # pragma: no cover - scipy is pinned in the lock file
        print("scipy       : (not installed)")

    print("--- numpy build configuration ---")
    try:
        np.show_config()
    except Exception as exc:  # pragma: no cover - show_config is best-effort
        print("(numpy.show_config() unavailable:", exc, ")")


if __name__ == "__main__":
    main()
