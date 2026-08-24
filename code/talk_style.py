"""
talk_style.py
-------------
Shared style for the presentation figures (fig4_tornado, fig5_threshold,
fig6_three_shifts, fig7_answer_talk).

TWO STYLES.  `--style talk` (default) is the wavelength-coloured house style of
slides/intro_deck.html.  `--style st` is the Science Tokyo deck style of
presentaion/slide_style_v3.md: colour encodes ARGUMENTATIVE ROLE, not category
(main subject = Science Blue, critical structure = Deep Indigo, support field =
Mist Lavender, text and axes = Midnight Ink), and the in-plot annotation budget
drops to 0-3 because the slide, not the figure, carries the interpretation.

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

import os
import sys

import matplotlib
from matplotlib import font_manager, rcParams

# --- palettes ---------------------------------------------------------------
# 'talk': the wavelength-as-colour house style (validated as a 3-slot
#         categorical set: 405/473/532 only -- 457 nm is perceptually the same
#         colour as 473 nm and cannot join them).
# 'st'  : Science Tokyo, colour by argumentative role (slide_style_v3.md S3).
PALETTES = {
    'talk': dict(
        W405='#6D3AA8', W473='#0E86C4', W532='#3E9B45',
        CRIT='#6D3AA8', OP='#0E86C4', PAST='#8A9BA6', INK='#111A21', INK2='#3D505C', MUTED='#6D8391',
        RULE='#CBD7DE', BAD='#B4432E', BAND='#D6EDF9', LEAN=False,
    ),
    'st': dict(
        # Science Blue = the subject; Deep Indigo = the critical structure
        # (the ionisation edge, the operating line, the crossover); the two
        # inherited laser lines are neutral references, not team colours.
        W405='#6E7488', W473='#1C3177', W532='#6E7488',
        CRIT='#4B0082', OP='#4B0082', PAST='#6E7488', INK='#101426', INK2='#101426', MUTED='#6E7488',
        RULE='#C7CBD8', BAD='#4B0082', BAND='#E8EAF1', LEAN=True,
    ),
}
SCIENCE_BLUE = '#1C3177'

W405, W473, W532 = (PALETTES['talk'][k] for k in ('W405', 'W473', 'W532'))
PAST, INK, INK2 = (PALETTES['talk'][k] for k in ('PAST', 'INK', 'INK2'))
MUTED, RULE, BAD = (PALETTES['talk'][k] for k in ('MUTED', 'RULE', 'BAD'))
BAND = PALETTES['talk']['BAND']
CRIT, OP = PALETTES['talk']['CRIT'], PALETTES['talk']['OP']
LEAN = False
STYLE = 'talk'
SUF = ''
ACCENT = W473

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


def _requested_style(explicit=None):
    if explicit:
        return explicit
    if '--style' in sys.argv:
        return sys.argv[sys.argv.index('--style') + 1]
    return os.environ.get('TALK_STYLE', 'talk')


def positional():
    """sys.argv[1:] with the '--style <name>' pair removed."""
    a = sys.argv[1:]
    if '--style' in a:
        i = a.index('--style')
        a = a[:i] + a[i + 2:]
    return a


def use(base=None, style=None):
    """Apply the rcParams for the requested style.  Axis labels stay in English
    (physics convention and font-safe); only the annotations switch language."""
    g = globals()
    name = _requested_style(style)
    if name not in PALETTES:
        raise SystemExit(f'unknown style {name!r}; choose from {sorted(PALETTES)}')
    g.update(PALETTES[name])
    # a slide-native figure is scaled to ~0.7 on the page, so the Science Tokyo
    # variant is typeset larger to stay legible from presentation distance
    base = (18.0 if name == 'st' else 15.0) if base is None else base
    g['STYLE'] = name
    g['SUF'] = '_st' if name == 'st' else ''
    g['ACCENT'] = SCIENCE_BLUE if name == 'st' else g['W473']
    # Science Tokyo asks for Arial; Liberation Sans is its metric-compatible
    # stand-in where Arial is not installed.
    sans = (['Arial', 'Liberation Sans', 'DejaVu Sans'] if name == 'st'
            else ['DejaVu Sans'])
    family = ([JP] + sans) if JP else sans
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
