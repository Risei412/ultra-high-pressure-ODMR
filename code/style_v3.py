"""
style_v3.py
Figure style for the ultra-high-pressure ODMR deck / reports.

Implements the fixed palette and the colour grammar of
``presentaion/slide_style_v3.md``: colour encodes the *argumentative role*
of an object, not its category.

    SCIENCE_BLUE   main scientific subject being followed through the argument
    DEEP_INDIGO    exceptional / critical / symmetry-protected structure
    MIST_LAVENDER  support field: uncertainty band, reference regime
    MIDNIGHT_INK   text, axes, neutral geometry (replaces pure black)
    NEUTRAL        muted derivative of Midnight Ink: comparison / baseline
    OPTIC_WHITE    background, active whitespace

Two objects that share an argumentative role are separated by line style or
marker shape, never by adding a new hue.
"""
from matplotlib import rcParams

# --------------------------------------------------------- fixed palette ----
SCIENCE_BLUE = '#1C3177'
OPTIC_WHITE = '#FFFFFF'
DEEP_INDIGO = '#4B0082'
MIST_LAVENDER = '#E8EAF1'
MIDNIGHT_INK = '#101426'

# neutral comparison series: a muted derivative of Midnight Ink, deliberately
# subordinate to Science Blue and *not* a sixth brand colour.
NEUTRAL = '#6B7080'
NEUTRAL_FAINT = '#9AA0AC'


def use_style():
    """Apply the deck typography and neutral chrome to Matplotlib."""
    rcParams.update({
        'font.size': 11,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'figure.facecolor': OPTIC_WHITE,
        'axes.facecolor': OPTIC_WHITE,
        'savefig.facecolor': OPTIC_WHITE,
        'axes.linewidth': 0.9,
        'axes.edgecolor': MIDNIGHT_INK,
        'axes.labelcolor': MIDNIGHT_INK,
        'axes.titlecolor': MIDNIGHT_INK,
        'text.color': MIDNIGHT_INK,
        'xtick.color': MIDNIGHT_INK,
        'ytick.color': MIDNIGHT_INK,
        'grid.color': NEUTRAL_FAINT,
        'grid.linewidth': 0.6,
        'legend.frameon': False,
    })


def panel_title(ax, text):
    """Descriptive, non-argumentative panel title (see style guide, sec. 9)."""
    ax.set_title(text, loc='left', weight='bold', fontsize=11,
                 color=MIDNIGHT_INK)
