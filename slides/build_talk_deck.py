"""
build_talk_deck.py
==================
Rebuild the 8-minute internal talk deck from the previous draft.

The draft (``odmr_474nm_talk.pptx``) was written against the pre-revision
numbers, when the fixed-power optimum was believed to be the absorption maximum
at 474.0 nm.  Two things changed in the model and one thing changed in the
argument:

  * at equal optical power the photon flux is I/E_gamma, so the optimum is the
    maximum of lambda*sigma_abs (475.5 nm), not of sigma_abs (474.0 nm);
  * the +/-15% model-form range on the single effective phonon energy is now
    propagated, which widens the Monte Carlo band to about +/-5.5 nm and makes
    hbar-omega the largest single contributor;
  * consequently the claim is stated as a window (463-488 nm, with the
    commercial 473 nm line inside it) rather than as a point.

This script applies those changes to the deck: it reorders the slides so the
mechanism precedes the answer, merges the two validation slides into one, adds
a slide for the uncalibrated power dependence, replaces every generated figure
with the output of ``talk_figs.py``, and rewrites the text.

Run:
    python talk_figs.py                       # figures first
    python build_talk_deck.py SRC.pptx        # -> odmr_473nm_talk.pptx

Slide map (draft -> talk):
    1 title            -> 1   headline becomes 473 nm, x2.5, x6
    2 DAC motivation   -> 2   unchanged
    3 NV imaging       -> 3   x2.7/x7 -> x2.5/x6
    4 two edges        -> 4   fig4
    6 cancellation     -> 5   fig5, moved ahead of the answer
    5 the optimum      -> 6   fig6, restated as a window
    8 crossover        -> 7   fig7, absorbs the 26/26 line from draft slide 7
    -                  -> 8   NEW: power dependence, fig8
    9 conditions       -> 9   numbers updated
    7 reproduction     -> (dropped; folded into 7)
"""
import copy
import os
import re
import shutil
import sys
import zipfile

from lxml import etree

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'figs')
DEFAULT_SRC = os.path.join(HERE, 'odmr_474nm_talk.pptx')
OUT = os.path.join(HERE, 'odmr_473nm_talk.pptx')

A = 'http://schemas.openxmlformats.org/drawingml/2006/main'
P = 'http://schemas.openxmlformats.org/presentationml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
PKG = 'http://schemas.openxmlformats.org/package/2006/relationships'
CT = 'http://schemas.openxmlformats.org/package/2006/content-types'
NS = {'a': A, 'p': P, 'r': R}

ACCENT = '1C3177'          # the deck's navy, used for the one line that matters


# ==========================================================================
# content
# ==========================================================================
# Markup inside a paragraph:
#   **bold**      bold, body colour
#   ##accent##    bold, deck navy -- reserved for the sentence of the slide
#   _{x}          subscript
# A shape's runs inherit their size/typeface from the run already in the deck,
# so only weight and colour are specified here.

TITLES = {
    1: '120 GPa の ODMR は 473 nm で測るべきである',
    4: '120 GPa では吸収帯とイオン化端が同時に青へ動く',
    5: '最適波長を動かすのは光学入力だけである',
    6: '答えは点ではなく窓である — 463–488 nm',
    7: '先行実験を再現し、73 GPa で符号が変わる',
    8: '最適波長はパワーで青へ動く — ここは未校正である',
    9: '473 nm が成り立つ条件、成り立たない条件',
}

BODY = {
    (1, 'Content Placeholder 2'): [
        '緑 532 nm 比で感度 ×2.5、磁気マップ取得時間 ×6',
    ],
    (4, 'TextBox 3'): [
        'ZPL が +0.40 eV 青方偏移し、Huang–Rhys 因子が 3.08 → 4.61 に増える。'
        '吸収極大は 586 → 474 nm。',
        '同時に IP(³A₂) が 2.68 → 3.06 eV に上がる。3.06 eV = 405 nm。'
        'これより青は NV⁻ を直接 NV⁰ に変える。',
        '120 GPa では σ(473) は σ(532) の **10 倍**。',
        'η ∝ Δν/(C√R) はコントラストに線形、光子レートに √ で効く。'
        '釣り合いの位置は量的問題である。',
    ],
    (5, 'TextBox 3'): [
        '励起状態イオン化と再結合はどちらも σ_{abs} に比例するので、'
        '定常 NV⁻ 分率 f₋ から相殺する。f₋ は 463–488 nm で 0.582–0.583、'
        '**0.2% しか動かない**。',
        '校正でフィットする現象論定数は **1 つも λ_{opt} を動かさない**。',
        '効くのは 3 つだけ — ZPL 位置、電子格子結合、そして'
        '**単一有効フォノン近似そのもの**。最大の ±6.9 nm は他人の測定値ではなく'
        '我々のモデル形である。',
        '##だから答えは一点ではなく窓で述べる。##',
    ],
    (5, 'TextBox 4'): [
        '各入力をモンテカルロ範囲の端に振り、他を公称値に固定したときの λopt の変化。'
        '等パワー比較では光子束が λ に比例するので、最適は σabs 極大 474.0 nm ではなく '
        'λσabs 極大 475.5 nm。',
    ],
    (6, 'TextBox 4'): [
        '**473 nm DPSS は最適の 0.2% 落ち** — 窓の中に市販線がある。'
        '457 nm で ×1.11、488 nm で ×1.05。',
        'λ_{opt} は ZPL に追随して 0.4–0.6 nm/GPa で動く。100 GPa では 486 nm。',
        '97–172 GPa を通す実験なら 488 nm と 473 nm の 2 本立てが有利。',
    ],
    (6, 'TextBox 5'): [
        '帯は入力のモンテカルロ 16–84%。λopt = 475.5 +5.6 / −5.3 nm。'
        '単一有効フォノン近似の ±15% を含む。',
    ],
    (7, 'TextBox 3'): [
        '較正に使ったのは Dai 2022 の 3 点だけ（定数 2 個）。'
        'その上で 0–150 GPa・4 波長の文献値を **26/26 再現**する（1 件は OPEN と申告）。',
        '50 GPa で青と緑を比較した既報は、**青の明確な優位を見つけなかった**。'
        'モデルはそこで 0.56、つまり緑が有利と答える。'
        '見つからなかったことが再現になっている。',
        '単一圧力の比較に乗る波長依存の系統誤差は圧力に依らないので、'
        '比が 1 を横切る位置を作ることも動かすこともできない。',
        '##これが本研究を反証する測定である。##',
    ],
    (7, 'TextBox 4'): [
        'P. Bhattacharyya, PhD thesis, UC Berkeley (2022), Sec. 6.3;  '
        'J.-H. Dai et al., Chin. Phys. Lett. 39, 117601 (2022);  '
        'A. Hilberer et al., PRB 107, L220102 (2023).',
    ],
    (8, 'TextBox 3'): [
        '稜線は 475.5 nm (u→0) → 465 nm (u = 0.1) → 447 nm (u = 0.3)。'
        '405 nm のイオン化端で止まる。',
        '473 nm 固定のコストは u = 0.2 で ×1.22、u = 0.3 で ×1.51。'
        '**緑に対する ×2.5 の大半を食う。**',
        'Dai 2022 の線形 PL は既存の高圧実験を u ≲ 0.3 に置く。'
        '移動の**向き**は全ドローで一致するが、**量は校正前のシナリオである**。',
        '##次の一手：473 nm 単色でパワー掃引 1 本。##'
        'R(I) の線形離脱点が「低励起」を実パワーに換算する。波長比較より先にやる。',
    ],
    (8, 'TextBox 4'): [
        'u = 1 は λopt における NV⁻ 遷移の半飽和強度。'
        '未校正の量は docs/next_step_power_dependence_experiment.md の Stage 1–3 で決める。',
    ],
    (9, 'TextBox 3'): [
        '検出帯も ZPL とともに青へ動く。常圧で選んだパスバンドは 120 GPa で'
        '発光の 26% しか拾わない。青励起にして初めて検出カットオンも動かせ、'
        '×1.6 が上乗せされる — 合計 ×4。',
        '##120 GPa では 532 nm を 473 nm DPSS に替える。'
        '同じ磁気マップが 1/6 の時間で撮れる。##',
    ],
}

# Whole-string substitutions for shapes that only need a number changed.
PATCHES = {
    3: [('×2.7', '×2.5'), ('×7', '×6')],
}

# table cell -> new text  {slide: {(row, col): text}}
TABLES = {
    6: {
        (1, 1): '3.76',
        (2, 1): '1.002', (2, 2): '市販 DPSS — 窓の中',
        (3, 1): '2.53',
    },
    9: {
        (1, 2): '平坦キュレット (α = 0.56) なら 510 nm。その形状は 40–50 GPa で '
                'ODMR コントラストを失う — micropillar は前提条件',
        (2, 1): '約 73 GPa 以上',
        (3, 2): '未校正（前ページ）。u = 0.3 では最適が 447 nm へ動き、'
                '473 nm 固定は ×1.51',
    },
}

NOTES = {
    1: '0:00–0:40  結論を先に置く。答えは窓で、その中に市販の 473 nm DPSS がある。',
    2: '0:40–1:35  DAC とは何かを図で見せる。狙う圧力帯 → 試料が小さすぎる → '
       '抵抗では決着しない、の 3 段。表は H₃S と (La,Ce)H₉ だけ読む。',
    3: '1:35–2:20  探針はアンビル内。既に撮れている実測を見せ、'
       'マップ時間 ∝ 画素数/η² で感度が通貨だと言う。',
    4: '2:20–3:20  2 つの端が窓を挟み込む。常圧曲線は単一有効フォノン近似が粗い旨を'
       '先に断る（次ページの伏線）。',
    5: '3:20–3:55  相殺 → 現象論非依存。ゼロが並ぶ絵が主張。'
       'ただし上の 3 本のうち最大は自分のモデル形であることを言い、'
       '「点ではなく窓」へ渡す。λσ の 1.5 nm は聞かれたら答える程度に留める。',
    6: '3:55–4:45  答え。窓の広さと 473 nm の位置。'
       'MC 帯が窓の半分を占めることを隠さない。',
    7: '4:45–5:55  再現と反証。較正は 2 定数だけ、26/26、'
       '50 GPa の null も再現。符号変化は単一圧力の系統誤差では作れない。',
    8: '5:55–6:50  未解決。向きは頑健、量は未校正。'
       '既存実験がいる u ≲ 0.3 でコストが ×1.5 に達することを言い、'
       '473 nm 単色パワー掃引 1 本という次の一手に落とす。',
    9: '6:50–7:35  適用条件とまとめ。校正は 473 nm の検証ではなく、'
       '絶対感度・パワー依存・f₋ の直接測定のために行う。',
}

# figure for each slide position
PICTURES = {4: 'fig4_absorption.png', 5: 'fig5_tornado.png',
            6: 'fig6_window.png', 7: 'fig7_crossover.png',
            8: 'fig8_power.png'}


# ==========================================================================
# package helpers
# ==========================================================================
def unpack(src, work):
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    with zipfile.ZipFile(src) as z:
        z.extractall(work)


def repack(work, out):
    if os.path.exists(out):
        os.remove(out)
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(work):
            for f in files:
                full = os.path.join(root, f)
                z.write(full, os.path.relpath(full, work))


def load(path):
    return etree.parse(path)


def save(tree, path):
    tree.write(path, xml_declaration=True, encoding='UTF-8', standalone=True)


def rels_path(part):
    d, f = os.path.split(part)
    return os.path.join(d, '_rels', f + '.rels')


def add_rel(rels_file, rid, type_, target):
    tree = load(rels_file)
    root = tree.getroot()
    el = etree.SubElement(root, f'{{{PKG}}}Relationship')
    el.set('Id', rid)
    el.set('Type', type_)
    el.set('Target', target)
    save(tree, rels_file)


def next_rid(rels_file):
    tree = load(rels_file)
    used = {int(r.get('Id')[3:]) for r in tree.getroot()
            if r.get('Id', '').startswith('rId')}
    return f'rId{max(used) + 1}'


def add_override(work, partname, content_type):
    p = os.path.join(work, '[Content_Types].xml')
    tree = load(p)
    root = tree.getroot()
    el = etree.SubElement(root, f'{{{CT}}}Override')
    el.set('PartName', partname)
    el.set('ContentType', content_type)
    save(tree, p)


# ==========================================================================
# structural edits
# ==========================================================================
def duplicate_slide(work, src_name, new_name):
    """Copy ppt/slides/<src_name> to <new_name> with all package bookkeeping.

    The copy deliberately does NOT inherit the source's notes slide: a notes
    slide points back at exactly one slide, so the new slide gets its own.
    """
    sd = os.path.join(work, 'ppt', 'slides')
    shutil.copy(os.path.join(sd, src_name), os.path.join(sd, new_name))

    src_rels = load(os.path.join(sd, '_rels', src_name + '.rels'))
    for rel in list(src_rels.getroot()):
        if rel.get('Type').endswith('/notesSlide'):
            src_rels.getroot().remove(rel)
    save(src_rels, os.path.join(sd, '_rels', new_name + '.rels'))

    add_override(work, f'/ppt/slides/{new_name}',
                 'application/vnd.openxmlformats-officedocument.'
                 'presentationml.slide+xml')

    pres_rels = os.path.join(work, 'ppt', '_rels', 'presentation.xml.rels')
    rid = next_rid(pres_rels)
    add_rel(pres_rels, rid,
            'http://schemas.openxmlformats.org/officeDocument/2006/'
            'relationships/slide', f'slides/{new_name}')
    return rid


def duplicate_notes(work, src_notes, new_notes, slide_name):
    nd = os.path.join(work, 'ppt', 'notesSlides')
    shutil.copy(os.path.join(nd, src_notes), os.path.join(nd, new_notes))

    tree = load(os.path.join(nd, '_rels', src_notes + '.rels'))
    for rel in tree.getroot():
        if rel.get('Type').endswith('/slide'):
            rel.set('Target', f'../slides/{slide_name}')
    save(tree, os.path.join(nd, '_rels', new_notes + '.rels'))

    add_override(work, f'/ppt/notesSlides/{new_notes}',
                 'application/vnd.openxmlformats-officedocument.'
                 'presentationml.notesSlide+xml')

    slide_rels = os.path.join(work, 'ppt', 'slides', '_rels',
                              slide_name + '.rels')
    add_rel(slide_rels, next_rid(slide_rels),
            'http://schemas.openxmlformats.org/officeDocument/2006/'
            'relationships/notesSlide', f'../notesSlides/{new_notes}')


def slide_order(work):
    """Map presentation order -> (slide part name, notes part name or None)."""
    pres = load(os.path.join(work, 'ppt', 'presentation.xml'))
    rels = load(os.path.join(work, 'ppt', '_rels', 'presentation.xml.rels'))
    target = {r.get('Id'): r.get('Target') for r in rels.getroot()}
    out = []
    for sld in pres.getroot().find(f'{{{P}}}sldIdLst'):
        name = os.path.basename(target[sld.get(f'{{{R}}}id')])
        sr = os.path.join(work, 'ppt', 'slides', '_rels', name + '.rels')
        notes = None
        for rel in load(sr).getroot():
            if rel.get('Type').endswith('/notesSlide'):
                notes = os.path.basename(rel.get('Target'))
        out.append((name, notes))
    return out


def reorder(work, rids):
    p = os.path.join(work, 'ppt', 'presentation.xml')
    tree = load(p)
    lst = tree.getroot().find(f'{{{P}}}sldIdLst')
    keep = {}
    for sld in list(lst):
        keep[sld.get(f'{{{R}}}id')] = sld
        lst.remove(sld)
    for i, rid in enumerate(rids):
        el = keep.get(rid)
        if el is None:                       # a slide added after the draft
            el = etree.Element(f'{{{P}}}sldId')
            el.set(f'{{{R}}}id', rid)
        el.set('id', str(256 + i))
        lst.append(el)
    save(tree, p)


def drop_unreferenced(work):
    """Delete slides (and their notes/rels) no longer in <p:sldIdLst>."""
    live = {name for name, _ in slide_order(work)}
    sd = os.path.join(work, 'ppt', 'slides')
    removed = []
    for name in sorted(os.listdir(sd)):
        if not name.endswith('.xml') or name in live:
            continue
        notes = None
        rp = os.path.join(sd, '_rels', name + '.rels')
        for rel in load(rp).getroot():
            if rel.get('Type').endswith('/notesSlide'):
                notes = os.path.basename(rel.get('Target'))
        os.remove(os.path.join(sd, name))
        os.remove(rp)
        removed.append(f'ppt/slides/{name}')
        if notes:
            nd = os.path.join(work, 'ppt', 'notesSlides')
            for f in (os.path.join(nd, notes),
                      os.path.join(nd, '_rels', notes + '.rels')):
                if os.path.exists(f):
                    os.remove(f)
            removed.append(f'ppt/notesSlides/{notes}')

    # drop the matching content-type overrides and presentation relationships
    ctp = os.path.join(work, '[Content_Types].xml')
    tree = load(ctp)
    for ov in list(tree.getroot()):
        pn = ov.get('PartName', '').lstrip('/')
        if pn in removed:
            tree.getroot().remove(ov)
    save(tree, ctp)

    prp = os.path.join(work, 'ppt', '_rels', 'presentation.xml.rels')
    tree = load(prp)
    used = {s.get(f'{{{R}}}id') for s in
            load(os.path.join(work, 'ppt', 'presentation.xml'))
            .getroot().find(f'{{{P}}}sldIdLst')}
    for rel in list(tree.getroot()):
        if rel.get('Type').endswith('/slide') and rel.get('Id') not in used:
            tree.getroot().remove(rel)
    save(tree, prp)

    # finally any media nobody points at
    referenced = set()
    for root, _, files in os.walk(os.path.join(work, 'ppt')):
        for f in files:
            if f.endswith('.rels'):
                for rel in load(os.path.join(root, f)).getroot():
                    referenced.add(os.path.basename(rel.get('Target')))
    md = os.path.join(work, 'ppt', 'media')
    for f in sorted(os.listdir(md)):
        if f not in referenced:
            os.remove(os.path.join(md, f))
            removed.append(f'ppt/media/{f}')
    return removed


# ==========================================================================
# text edits
# ==========================================================================
TOKEN = re.compile(r'(\*\*.+?\*\*|##.+?##|_\{.+?\})')


def _runs(text, bold=False, accent=False, sub=False):
    """Flatten the markup into (text, bold, accent, subscript) runs.

    Spans nest: ``**1 つも λ_{opt} を動かさない**`` has to come out as a bold
    run, a bold subscript run and a bold run, not as one literal string.
    """
    out = []
    for piece in TOKEN.split(text):
        if not piece:
            continue
        if piece.startswith('**'):
            out += _runs(piece[2:-2], True, accent, sub)
        elif piece.startswith('##'):
            out += _runs(piece[2:-2], True, True, sub)
        elif piece.startswith('_{'):
            out += _runs(piece[2:-1], bold, accent, True)
        else:
            out.append((piece, bold, accent, sub))
    return out


def _mk_rPr(base, bold=False, accent=False, sub=False):
    rPr = copy.deepcopy(base)
    rPr.set('b', '1' if bold else '0')
    if sub:
        rPr.set('baseline', '-25000')
    else:
        rPr.attrib.pop('baseline', None)
    if accent:
        clr = rPr.find(f'{{{A}}}solidFill/{{{A}}}srgbClr')
        if clr is not None:
            clr.set('val', ACCENT)
    return rPr


def set_text(sp, paragraphs):
    """Replace a shape's paragraphs, inheriting run/paragraph formatting."""
    tx = sp.find(f'{{{P}}}txBody')
    if tx is None:
        tx = sp.find(f'{{{A}}}txBody')
    old = tx.findall(f'{{{A}}}p')
    base_p = old[0]
    base_pPr = base_p.find(f'{{{A}}}pPr')
    base_rPr = base_p.find(f'{{{A}}}r/{{{A}}}rPr')   # None on notes slides,
    # whose runs carry no explicit formatting and inherit from the notes master
    for p in old:
        tx.remove(p)

    for text in paragraphs:
        para = etree.SubElement(tx, f'{{{A}}}p')
        if base_pPr is not None:
            para.append(copy.deepcopy(base_pPr))
        for piece, bold, accent, sub in _runs(text):
            run = etree.SubElement(para, f'{{{A}}}r')
            if base_rPr is not None:
                run.append(_mk_rPr(base_rPr, bold, accent, sub))
            t = etree.SubElement(run, f'{{{A}}}t')
            t.text = piece
            if piece != piece.strip():
                t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')


def shapes_by_name(tree):
    out = {}
    for sp in tree.getroot().iter():
        if sp.tag in (f'{{{P}}}sp', f'{{{P}}}pic', f'{{{P}}}graphicFrame'):
            cNvPr = sp.find(f'.//{{{P}}}cNvPr')
            if cNvPr is not None:
                out[cNvPr.get('name')] = sp
    return out


def patch_runs(tree, pairs):
    for t in tree.getroot().iter(f'{{{A}}}t'):
        if t.text:
            for old, new in pairs:
                if t.text == old:
                    t.text = new


def set_table(tree, cells):
    tbl = tree.getroot().find(f'.//{{{A}}}tbl')
    rows = tbl.findall(f'{{{A}}}tr')
    for (ri, ci), text in cells.items():
        tc = rows[ri].findall(f'{{{A}}}tc')[ci]
        set_text(tc, [text])


# ==========================================================================
# picture edits
# ==========================================================================
def replace_picture(work, slide_name, tree, png, new_media=None):
    """Swap the picture's image and refit it inside its original box."""
    from PIL import Image

    pic = tree.getroot().find(f'.//{{{P}}}pic')
    blip = pic.find(f'.//{{{A}}}blip')
    rid = blip.get(f'{{{R}}}embed')

    rp = os.path.join(work, 'ppt', 'slides', '_rels', slide_name + '.rels')
    rels = load(rp)
    rel = [r for r in rels.getroot() if r.get('Id') == rid][0]
    if new_media:
        rel.set('Target', f'../media/{new_media}')
        save(rels, rp)
        target = new_media
    else:
        target = os.path.basename(rel.get('Target'))
    shutil.copy(png, os.path.join(work, 'ppt', 'media', target))

    # refit: keep the original box, preserve the new aspect ratio, centre it
    xfrm = pic.find(f'{{{P}}}spPr/{{{A}}}xfrm')
    off, ext = xfrm.find(f'{{{A}}}off'), xfrm.find(f'{{{A}}}ext')
    x, y = int(off.get('x')), int(off.get('y'))
    bw, bh = int(ext.get('cx')), int(ext.get('cy'))
    iw, ih = Image.open(png).size
    scale = min(bw / iw, bh / ih)
    nw, nh = int(iw * scale), int(ih * scale)
    off.set('x', str(x + (bw - nw) // 2))
    off.set('y', str(y + (bh - nh) // 2))
    ext.set('cx', str(nw))
    ext.set('cy', str(nh))

    descr = pic.find(f'.//{{{P}}}cNvPr')
    descr.set('descr', os.path.basename(png))


# ==========================================================================
# main
# ==========================================================================
def main(src=DEFAULT_SRC, out=OUT):
    work = os.path.join(os.path.dirname(out), '_deck_build')
    unpack(src, work)

    order = slide_order(work)
    assert len(order) == 9, f'expected the 9-slide draft, found {len(order)}'

    # --- structure: duplicate slide 8 for the new power slide -------------
    new_rid = duplicate_slide(work, 'slide8.xml', 'slide10.xml')
    duplicate_notes(work, 'notesSlide8.xml', 'notesSlide10.xml', 'slide10.xml')

    rels = load(os.path.join(work, 'ppt', '_rels', 'presentation.xml.rels'))
    rid_of = {os.path.basename(r.get('Target')): r.get('Id')
              for r in rels.getroot() if r.get('Type').endswith('/slide')}
    reorder(work, [rid_of['slide1.xml'], rid_of['slide2.xml'],
                   rid_of['slide3.xml'], rid_of['slide4.xml'],
                   rid_of['slide6.xml'],          # cancellation, moved up
                   rid_of['slide5.xml'],          # the answer, moved down
                   rid_of['slide8.xml'],          # crossover
                   new_rid,                       # NEW power slide
                   rid_of['slide9.xml']])
    for gone in drop_unreferenced(work):
        print('dropped', gone)

    # --- content ----------------------------------------------------------
    for pos, (name, notes) in enumerate(slide_order(work), 1):
        path = os.path.join(work, 'ppt', 'slides', name)
        tree = load(path)
        shapes = shapes_by_name(tree)

        if pos in TITLES:
            set_text(shapes['Title 1'], [TITLES[pos]])
        for (sp_pos, sp_name), paras in BODY.items():
            if sp_pos == pos:
                set_text(shapes[sp_name], paras)
        if pos in PATCHES:
            patch_runs(tree, PATCHES[pos])
        if pos in TABLES:
            set_table(tree, TABLES[pos])
        if pos in PICTURES:
            fig = os.path.join(FIGS, PICTURES[pos])
            replace_picture(work, name, tree, fig,
                            new_media='image23.png' if pos == 8 else None)
        save(tree, path)

        if notes and pos in NOTES:
            np_ = os.path.join(work, 'ppt', 'notesSlides', notes)
            ntree = load(np_)
            body = None
            for sp in ntree.getroot().iter(f'{{{P}}}sp'):
                if sp.find(f'.//{{{A}}}t') is not None and \
                        sp.find(f'.//{{{P}}}ph[@type="body"]') is not None:
                    body = sp
            if body is not None:
                set_text(body, [NOTES[pos]])
                save(ntree, np_)
        print(f'slide {pos}: {name}')

    repack(work, out)
    shutil.rmtree(work)
    print('wrote', out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SRC,
         sys.argv[2] if len(sys.argv) > 2 else OUT)
