"""One place for the publication figure style.

Every figure in `paper/main.tex` is a full-width `figure*` included at
`width=\textwidth`, and in REVTeX 4.2 `aps,prapplied,reprint` that width is
about 510 pt = 7.06 in.  A figure drawn at some other width is *rescaled* by
the ratio, and its lettering is rescaled with it: the six-panel kernel figure
used to be drawn 15.6 in wide, so its 10 pt labels reached the page at 4.6 pt.

So the rule this module enforces is: draw at the width the page will use, and
set the font sizes to what they must be on the page.  Nothing downstream
rescales them.

    FULL_WIDTH_IN   draw a `figure*` this wide
    COLUMN_WIDTH_IN draw a single-column `figure` this wide
    ANNOT           font size for text drawn inside the axes

`save()` writes at exactly `figsize`, and `verify_width()` re-reads the PDF and
refuses a file that is not the width it was drawn at.  Nothing downstream is
allowed to rescale a figure, because rescaling a figure rescales its lettering.

APS asks that lettering be legible at final size.  8 pt for axis labels and
tick labels, 7 pt for legends and in-axes annotation, is the floor used here.
"""

import re

import matplotlib.pyplot as plt

# REVTeX 4.2, aps/prapplied/reprint, letterpaper.  Measured from the class:
# \textwidth = 510.0 pt, \columnwidth = 246.0 pt.  A figure drawn at exactly
# these widths is placed on the page at a scale of 1.
FULL_WIDTH_IN = 510.0 / 72.0
COLUMN_WIDTH_IN = 246.0 / 72.0

# In-axes annotation.  Never smaller than this.
ANNOT = 7.0

RCPARAMS = {
    # TrueType, not Type 3: the APS and arXiv submission checkers require it.
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.size': 8.5,
    'axes.titlesize': 9.0,
    'axes.labelsize': 8.5,
    'xtick.labelsize': 8.0,
    'ytick.labelsize': 8.0,
    'legend.fontsize': 7.0,
    'axes.linewidth': 0.7,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
    'figure.dpi': 200,
    'savefig.dpi': 400,
    'hatch.linewidth': 0.4,
}


def use(**overrides):
    """Apply the shared style; `overrides` is for the rare per-figure tweak."""
    plt.rcParams.update(RCPARAMS)
    if overrides:
        plt.rcParams.update(overrides)


def save(fig, png_path, pdf_path):
    """Write both formats at exactly `figsize`, so the page rescales nothing.

    `bbox_inches='tight'`, which this used to pass, is not a trim: it is a
    recomputed bounding box that *grows* when an artist sticks out past the
    canvas.  The six-panel kernel figure did exactly that -- one twin-axis
    ordinate label protruded, tight bbox widened the canvas to 7.29 in, and
    `\\includegraphics[width=\\textwidth]` then shrank the whole figure by
    3.2 %, putting its 7 pt annotation on the page at 6.8 pt.  Saving at the
    figure's own bounds removes the freedom: what is drawn 7.083 in wide is
    printed 7.083 in wide.

    The cost is that anything outside the canvas is clipped rather than
    accommodated, so margins have to be set by `subplots_adjust` in each
    script.  `verify_width` below is the guard.
    """
    for path in (png_path, pdf_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    verify_width(pdf_path, fig.get_size_inches()[0])
    return '%s\n%s' % (png_path, pdf_path)


def verify_width(pdf_path, expected_in, tolerance_pt=0.75):
    """Fail loudly if the written PDF is not the width it was drawn at."""
    data = pdf_path.read_bytes()
    match = re.search(rb'MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+'
                      rb'([\d.]+)\s+([\d.]+)', data)
    if match is None:                       # pragma: no cover - defensive
        raise AssertionError('no MediaBox in %s' % pdf_path)
    width_pt = float(match.group(3)) - float(match.group(1))
    expected_pt = expected_in * 72.0
    if abs(width_pt - expected_pt) > tolerance_pt:
        raise AssertionError(
            '%s is %.2f pt wide, drawn at %.2f pt: the page would rescale it '
            'by %.3f and shrink its lettering with it'
            % (pdf_path.name, width_pt, expected_pt, expected_pt / width_pt))
    return width_pt
