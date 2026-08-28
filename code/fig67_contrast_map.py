"""Test N1 against Bhattacharyya's Fig. 6.7(b).  It does not survive.

N1 (`thesis_from_the_core.py`) observes that Phi = (C/dnu)^2 R is EVEN in the
ODMR contrast C, so a resonance with C = -5% is exactly as good a magnetometer
as one with C = +5%, and the thesis's "positive contrast" is not a loss of
sensitivity.  What the core does force is that C, being continuous in stress,
must pass through ZERO between the negative and positive regimes -- and there
eta diverges.  In a widefield map the boundary between positive- and
negative-contrast regions is therefore a CONTOUR OF ZERO CONTRAST: a dark ODMR
line across an otherwise bright field.

STATUS listed this as possibly confirmable in existing data, and Fig. 6.7(b) --
a line cut across a [110] cut culet at ~30 GPa carrying both signs at once --
is the obvious place to look.  This module looks.  **The answer is no.**

## What the figure turns out to contain

Panel (b) is a raster, 3875 x 1491, in an INDEXED colour space whose 254-entry
palette is embedded in the PDF.  Every pixel therefore inverts to its palette
index exactly -- 100% of them, with no nearest-colour guessing -- and the
palette is reversed viridis, monotone in the plotted normalised fluorescence.
So the published figure yields a quantitative, signed contrast map:

    v     = (253 - index) / 253          colourbar coordinate, high = bright
    D     = v - baseline(v)              proportional to true contrast C

D is proportional to C and not equal to it: the colourbar maps fluorescence
affinely onto [0, 1] with limits the figure does not state, so only the SIGN
and the RATIOS of D are meaningful.  The sign is what N1 needs.

## What it shows, and why that is not N1

Thirteen resonance tracks survive ridge-linking over >= 40 rows, nine negative
and four positive.  Both signs are present, as the caption says.  But:

  * No track changes sign along the cut.  Every candidate turned out to be one
    line crossing another IN FREQUENCY, not one line's contrast passing
    through zero in POSITION.  The clearest candidate -- a positive track
    ending at y = 735 and another beginning at y = 1035 at nearly the same
    frequency -- is the deep negative line sweeping through that frequency in
    between, which the bridge test shows directly.

  * The negative lines do fade at BOTH ends of the cut -- summed negative
    weight falls from 285 in the middle to 102 at the ends -- which looks
    like N1 until you ask what is there instead.  Summed positive weight at
    the ends is 52 against a no-resonance floor of 24, i.e. baseline residual
    rather than lines, and NO positive resonance is identified beyond
    y = 1410 or below y = 280.  The one positive track that reaches into the
    collapsed zone is already present where the negatives are at full depth,
    so it is coexistence, not emergence.  Both signs die together: that is
    the edge of the useful field, not a zero-contrast contour.

    This is also why the summed positive weight is NOT the statistic to test
    on -- it counts noise and baseline residual, so it rises wherever the map
    empties out.  Scored on it, the figure would falsely "confirm" N1.

So the line cut of Fig. 6.7(b) does not cross a zero-contrast boundary within
its useful range, or the boundary is finer than this figure resolves.  N1
remains untested, and STATUS's "possibly confirmable in existing data" should
be struck.

## One methodological result worth keeping

The baseline matters more than it looks.  With a 401-px running median the
brightest positive feature in the figure reads D = +0.494; widening the window
to 801 px drops it to +0.090 and to 1201 px, +0.021.  It is mostly the median
being dragged down between two deep adjacent dips, not a positive resonance.
Two other positive features are stable across all three windows (+0.17 and
+0.13, moving by under 10%) and are real.  **Any contrast read off a published
waterfall has to carry this check**, and the default here is the 801-px window
with the 401-px value reported alongside so the difference stays visible.

Run for the tables.  Needs the thesis PDF at ../docs/ref/.
"""
import os

import numpy as np
from scipy.ndimage import median_filter, uniform_filter, uniform_filter1d

THESIS_PDF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          'docs', 'ref', 'Principle_and_Applications_of_.pdf')
# Fig. 6.7 sits on PDF page 105 (0-indexed 104); panel (b) is the larger raster.
FIG67_PAGE_INDEX = 104
PANEL_B_SHAPE = (1491, 3875)
PALETTE_ENTRIES = 254
# Baseline windows in pixels along the frequency axis.  801 is the default;
# 401 is kept because it is what a careless reading would use.
BASELINE_WINDOWS = (401, 801, 1201)
DEFAULT_WINDOW = 801
# Columns dropped at each edge, where a running baseline has no room.
EDGE_MARGIN = 150


def load_panel_b(path=THESIS_PDF):
    """Return (indices, palette) for Fig. 6.7(b), inverted exactly.

    The image is an /Indexed colour space, so the palette in the PDF is the
    colour map itself and the inversion is a lookup, not a fit.
    """
    from pypdf import PdfReader
    from PIL import Image
    import io

    reader = PdfReader(path)
    page = reader.pages[FIG67_PAGE_INDEX]

    def images(resources, depth=0, found=None):
        found = [] if found is None else found
        xobjects = resources.get('/XObject')
        if xobjects is None:
            return found
        for _, ref in xobjects.get_object().items():
            obj = ref.get_object()
            if obj.get('/Subtype') == '/Image':
                found.append(obj)
            elif depth < 4 and obj.get('/Resources') is not None:
                images(obj['/Resources'].get_object(), depth + 1, found)
        return found

    target = None
    for obj in images(page['/Resources'].get_object()):
        if (int(obj['/Height']), int(obj['/Width'])) == PANEL_B_SHAPE:
            target = obj
            break
    if target is None:
        raise ValueError('panel (b) not found; is this the right PDF?')

    space = target['/ColorSpace']
    palette = np.frombuffer(bytes(space[3]), dtype=np.uint8).reshape(-1, 3)

    data = None
    for image in page.images:
        candidate = np.asarray(Image.open(io.BytesIO(image.data)).convert('RGB'))
        if candidate.shape[:2] == PANEL_B_SHAPE:
            data = candidate.astype(np.int32)
            break
    if data is None:
        raise ValueError('panel (b) pixels not recoverable')

    keys = (palette[:, 0].astype(np.int32) << 16 | palette[:, 1].astype(np.int32) << 8
            | palette[:, 2].astype(np.int32))
    order = np.argsort(keys)
    pixel_keys = (data[..., 0] << 16 | data[..., 1] << 8 | data[..., 2]).ravel()
    position = np.clip(np.searchsorted(keys[order], pixel_keys), 0, len(keys) - 1)
    exact = float((keys[order][position] == pixel_keys).mean())
    return order[position].reshape(PANEL_B_SHAPE), palette, exact


def contrast_map(indices, window=DEFAULT_WINDOW, decimate=4):
    """Signed contrast proxy D = v - baseline(v), proportional to true C."""
    v = (float(PALETTE_ENTRIES - 1) - indices.astype(float)) / (PALETTE_ENTRIES - 1)
    # Baseline on a column-decimated copy, then interpolated back: a 801-px
    # median over 3875 columns is the slow step and decimation is exact enough
    # for a background that varies on the scale of the whole axis.
    small = v[:, ::decimate]
    base_small = median_filter(small, size=(1, max(window // decimate, 3)),
                               mode='nearest')
    columns = np.arange(v.shape[1])
    base = np.empty_like(v)
    source = np.arange(small.shape[1]) * decimate
    for row in range(v.shape[0]):
        base[row] = np.interp(columns, source, base_small[row])
    base = uniform_filter1d(base, 81, axis=1, mode='nearest')
    return uniform_filter(v - base, size=(9, 9))


def baseline_sensitivity(indices, probes=None):
    """How much each named feature depends on the baseline window."""
    probes = probes or {
        'bright top feature': (130, 2145),
        'positive line, mid': (600, 2110),
        'positive line, right': (700, 3239),
        'deep negative line': (600, 2159),
        'background': (600, 1700),
    }
    out = {name: {} for name in probes}
    for window in BASELINE_WINDOWS:
        D = contrast_map(indices, window)
        for name, (y, x) in probes.items():
            out[name][window] = float(D[y, x])
    for name, values in out.items():
        span = max(values.values()) - min(values.values())
        reference = abs(values[DEFAULT_WINDOW])
        out[name]['stable'] = bool(reference > 0.05 and span < 0.5 * reference)
    return out


def tracks(D, step=5, prominence=0.05, distance=40, tolerance=25, min_rows=40):
    """Ridge-link the resonances, keeping the sign of each."""
    from scipy.signal import find_peaks

    rows = range(0, D.shape[0], step)
    found = []
    for y in rows:
        row = D[y]
        detections = []
        for sign in (+1, -1):
            peaks, _ = find_peaks(sign * row, prominence=prominence,
                                  distance=distance)
            detections += [(int(x), float(row[x])) for x in peaks]
        detections.sort()
        used = set()
        for track in found:
            if track['last'] != y - step:
                continue
            best = None
            for i, (x, value) in enumerate(detections):
                if i in used:
                    continue
                gap = abs(x - track['x'][-1])
                if gap <= tolerance and (best is None or gap < best[0]):
                    best = (gap, i, x, value)
            if best:
                _, i, x, value = best
                used.add(i)
                track['y'].append(y)
                track['x'].append(x)
                track['v'].append(value)
                track['last'] = y
        for i, (x, value) in enumerate(detections):
            if i not in used:
                found.append({'y': [y], 'x': [x], 'v': [value], 'last': y})
    found = [t for t in found if len(t['y']) >= min_rows]
    found.sort(key=lambda t: -len(t['y']))
    for track in found:
        values = np.array(track['v'])
        track['sign'] = ('mixed' if values.min() < 0 < values.max()
                         else ('positive' if values.min() > 0 else 'negative'))
    return found


def bridge(D, lower, upper, half=45, step=20):
    """Read D between the ends of two tracks along the interpolated frequency.

    The test that kills the best N1 candidate: if the gap between a positive
    track ending and another beginning is a zero crossing, D stays small
    through it.  If it is one line crossing another, D dives.
    """
    y0, x0 = lower['y'][-1], lower['x'][-1]
    y1, x1 = upper['y'][0], upper['x'][0]
    values = []
    for y in range(y0, y1 + 1, step):
        fraction = (y - y0) / max(y1 - y0, 1)
        x = int(round(x0 + fraction * (x1 - x0)))
        window = D[y, max(x - half, 0):x + half + 1]
        values.append(float(window.max() if abs(window.max()) > abs(window.min())
                            else window.min()))
    values = np.array(values)
    return {'y_from': int(y0), 'y_to': int(y1),
            'D_from': float(lower['v'][-1]), 'D_to': float(upper['v'][0]),
            'min': float(values.min()), 'max': float(values.max()),
            'is_zero_crossing': bool(abs(values).max() < 3.0 * max(
                abs(lower['v'][-1]), abs(upper['v'][0])))}


def sign_weight_profile(D, bins=25, margin=EDGE_MARGIN):
    """Negative and positive contrast weight against position along the cut.

    N1 predicts that where the negative lines vanish the positive ones appear.
    If both fall together, the cut has simply run off the useful field.
    """
    core = D[:, margin:D.shape[1] - margin]
    negative = np.where(core < 0.0, -core, 0.0).sum(1)
    positive = np.where(core > 0.0, core, 0.0).sum(1)
    edges = np.linspace(0, D.shape[0], bins + 1).astype(int)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        rows.append({'y_from': int(lo), 'y_to': int(hi),
                     'negative': float(negative[lo:hi].mean()),
                     'positive': float(positive[lo:hi].mean())})
    middle = [r for r in rows if 0.25 * D.shape[0] < r['y_from'] < 0.75 * D.shape[0]]
    ends = [rows[0], rows[-1]]
    return {
        'rows': rows,
        'negative_middle': float(np.mean([r['negative'] for r in middle])),
        'negative_ends': float(np.mean([r['negative'] for r in ends])),
        'positive_middle': float(np.mean([r['positive'] for r in middle])),
        'positive_ends': float(np.mean([r['positive'] for r in ends])),
    }


def verdict(profile, found, D=None):
    """Does the figure show a zero-contrast contour?

    The summed positive weight is NOT the statistic to test on: it counts the
    noise floor and the baseline residual in the wings of the deep negative
    lines, so it rises wherever the map is empty.  The statistic that means
    something is whether an identified positive RESONANCE lives in the zone
    where the negatives have collapsed.  It does not.
    """
    rows = profile['rows']
    peak = max(r['negative'] for r in rows)
    collapsed = [r for r in rows if r['negative'] < 0.5 * peak]
    zone = [(r['y_from'], r['y_to']) for r in collapsed]
    positives_in_zone = 0
    for track in found:
        if track['sign'] != 'positive':
            continue
        if any(lo <= y < hi for y in track['y'] for lo, hi in zone):
            positives_in_zone += 1
    noise = None
    if D is not None:
        # positive weight the map would carry with no resonances at all
        core = D[:, EDGE_MARGIN:D.shape[1] - EDGE_MARGIN]
        scale = float(np.median(np.abs(core)))
        noise = float(np.where(core > 0.0, core, 0.0).sum(1).min())
    return {
        'tracks': len(found),
        'signs': {s: sum(1 for t in found if t['sign'] == s)
                  for s in ('negative', 'positive', 'mixed')},
        'any_track_changes_sign': any(t['sign'] == 'mixed' for t in found),
        'collapsed_zone': zone,
        'negative_collapses_at_the_ends': bool(collapsed),
        'positive_tracks_in_the_collapsed_zone': positives_in_zone,
        'positive_weight_noise_floor': noise,
        # N1 names ONE observable: a resonance whose contrast changes sign
        # along the cut.  The collapsed-zone counts above are context, not the
        # test -- a positive line that merely persists into the zone while
        # fading is coexistence, which the figure shows everywhere.
        'n1_confirmed': bool(any(t['sign'] == 'mixed' for t in found)),
    }


def report():
    indices, palette, exact = load_panel_b()
    out = {'palette_entries': len(palette), 'exact_inversion': exact}
    print('Fig. 6.7(b) as data')
    print(f'      panel {indices.shape[1]} x {indices.shape[0]}, '
          f'{len(palette)}-entry indexed palette embedded in the PDF')
    print(f'      exact palette inversion: {exact * 100:.2f}% of pixels')

    print('\nbaseline sensitivity (D at named features, by window)')
    out['baseline'] = sensitivity = baseline_sensitivity(indices)
    print(f'      {"feature":24}' + ''.join(f'{w:>9}' for w in BASELINE_WINDOWS)
          + '   verdict')
    for name, values in sensitivity.items():
        print(f'      {name:24}'
              + ''.join(f'{values[w]:+9.3f}' for w in BASELINE_WINDOWS)
              + ('   stable' if values['stable'] else '   BASELINE ARTEFACT'))

    D = contrast_map(indices, DEFAULT_WINDOW)
    out['found'] = found = tracks(D)
    print(f'\nresonance tracks (>= 40 rows, {DEFAULT_WINDOW}-px baseline)')
    print(f'      {"#":>3}{"rows":>6}{"y range":>14}{"x range":>14}'
          f'{"D min":>9}{"D max":>9}  sign')
    for i, track in enumerate(found[:12]):
        values = np.array(track['v'])
        span_y = '{}-{}'.format(track['y'][0], track['y'][-1])
        span_x = '{}-{}'.format(min(track['x']), max(track['x']))
        print(f'      {i:3d}{len(track["y"]):6d}{span_y:>14}{span_x:>14}'
              f'{values.min():9.3f}{values.max():9.3f}  {track["sign"]}')

    positives = [t for t in found if t['sign'] == 'positive']
    if len(positives) >= 2:
        pair = sorted(positives, key=lambda t: t['y'][0])
        lower = max((t for t in pair if t['y'][-1] < pair[-1]['y'][0]),
                    key=lambda t: t['y'][-1], default=None)
        if lower is not None:
            out['bridge'] = result = bridge(D, lower, pair[-1])
            print('\nbridge test on the best N1 candidate')
            print(f'      positive track ends y={result["y_from"]} '
                  f'(D={result["D_from"]:+.3f}), next begins y={result["y_to"]} '
                  f'(D={result["D_to"]:+.3f})')
            print(f'      D between them spans {result["min"]:+.3f} .. '
                  f'{result["max"]:+.3f}')
            print('      -> ' + ('a zero crossing' if result['is_zero_crossing']
                                 else 'the deep negative line sweeping through: '
                                      'a frequency crossing, not a sign change'))

    out['profile'] = profile = sign_weight_profile(D)
    print('\ncontrast weight against position along the cut')
    print(f'      {"y":>12}{"sum |neg|":>12}{"sum pos":>12}')
    for row in profile['rows'][::4]:
        span = '{}-{}'.format(row['y_from'], row['y_to'])
        print(f'      {span:>12}{row["negative"]:12.1f}'
              f'{row["positive"]:12.1f}')
    print(f'      middle: neg {profile["negative_middle"]:.0f}, '
          f'pos {profile["positive_middle"]:.0f}')
    print(f'      ends  : neg {profile["negative_ends"]:.0f}, '
          f'pos {profile["positive_ends"]:.0f}')

    out['verdict'] = check = verdict(profile, found, D)
    print('\nverdict')
    print(f'      tracks {check["tracks"]}, signs {check["signs"]}')
    print(f'      any track changes sign along the cut: '
          f'{check["any_track_changes_sign"]}')
    print(f'      zone where the negatives have collapsed: '
          f'{check["collapsed_zone"]}')
    print(f'      identified POSITIVE resonances in that zone: '
          f'{check["positive_tracks_in_the_collapsed_zone"]}')
    print(f'      (summed positive weight there is {profile["positive_ends"]:.0f} '
          f'against a no-resonance floor of '
          f'{check["positive_weight_noise_floor"]:.0f}, so the sum is baseline '
          f'residual,\n       not lines -- which is why the sum is not the '
          f'statistic to test on)')
    print(f'      N1 confirmed: {check["n1_confirmed"]}  '
          '(criterion: a resonance changes sign along the cut)')
    print('      Both signs die together at the ends, so that is the edge of '
          'the useful\n      field, not a zero-contrast contour.  The one '
          'positive line reaching into\n      the collapsed zone is already '
          'there while the negatives are strong, so it\n      is coexistence, '
          'not an emergence.  N1 remains untested.')
    return out


if __name__ == '__main__':
    report()
