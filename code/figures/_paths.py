"""Shared locations for figure generation, independent of the working directory."""
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parents[1]
FIGURE_DIR = CODE_DIR / 'outputs' / 'figures'


def figure_output(filename):
    """Create the output directory and return the path for a generated figure."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    return FIGURE_DIR / filename
