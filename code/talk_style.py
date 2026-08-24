"""
talk_style.py
-------------
Shared style for the presentation figures (fig4_tornado, fig5_threshold,
fig6_three_shifts, fig7_answer_talk).

Two things differ from the manuscript figures:

  * ONE MESSAGE PER FIGURE.  No twin axes, no three-panel composites, no
    comparison curves that are not the point.  A figure projected for 60 s in
    an eight-minute talk can carry one idea.
  * THE PALETTE IS THE SUBJECT.  Each categorical colour is the appearance of
    the laser line it labels, as in slides/intro_deck.html.  Only THREE lines
    are ever drawn together -- 405 (the boundary), 473 (the answer), 532 (the
    inherited default) -- because 457 and 473 nm are perceptually the same
    colour: as a fourth simultaneous series they fall to dE 12.2 against 473 nm
    for normal vision, below the legibility floor.  457 nm belongs in the table
    and in the speech, not as a fourth coloured curve.
"""

import matplotlib
from matplotlib import font_manager, rcParams

# --- palette (light surface; validated as a 3-slot categorical set) ---------
W405   = '#6D3AA8'      # violet   -- the ionisation boundary
W473   = '#0E86C4'      # cyan     -- the answer
W532   = '#3E9B45'      # green    -- the inherited default
ACCENT = W473
PAST   = '#8A9BA6'      # ambient-pressure reference curves
INK    = '#111A21'
INK2   = '#3D505C'
MUTED  = '#6D8391'
RULE   = '#CBD7DE'
BAD    = '#B4432E'
BAND   = '#D6EDF9'      # tolerance / accepted region

_JP_FONTS = ('IPAPGothic', 'IPAGothic', 'Noto Sans CJK JP', 'Noto Sans JP',
             'Hiragino Kaku Gothic ProN', 'Yu Gothic', 'TakaoPGothic')


def _jp_font():
    have = {f.name for f in font_manager.fontManager.ttflist}
    for name in _JP_FONTS:
        if name in have:
            return name
    return None


JP = _jp_font()


def t(ja, en):
    """Japanese annotation when a CJK font is installed, English otherwise."""
    return ja if JP else en


def use(base=15.0):
    """Apply the talk rcParams.  Axis labels stay in English (physics
    convention and font-safe); only the annotations switch language."""
    family = [JP, 'DejaVu Sans'] if JP else ['DejaVu Sans']
    rcParams.update({
        'font.family': family,
        'font.size': base,
        'axes.titlesize': base + 1,
        'axes.labelsize': base,
        'axes.labelcolor': INK,
        'axes.edgecolor': RULE,
        'axes.linewidth': 1.1,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.unicode_minus': False,
        'xtick.color': INK2, 'ytick.color': INK2,
        'xtick.labelsize': base - 1, 'ytick.labelsize': base - 1,
        'text.color': INK,
        'grid.color': RULE, 'grid.alpha': 0.45, 'grid.linewidth': 0.8,
        'legend.frameon': False,
        'legend.fontsize': base - 1.5,
        'mathtext.fontset': 'dejavusans',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'savefig.dpi': 200,
        'savefig.bbox': 'tight',
    })
    matplotlib.use('Agg')
