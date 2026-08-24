"""
build_standalone.py
Inline every <img src="..."> of a deck as a data: URI so it is a single
self-contained file (needed for hosting, e-mail, or offline viewing).

Run from the repository root or from slides/:
    python slides/build_standalone.py                  # intro_deck.html
    python slides/build_standalone.py talk_deck.html   # any other deck
Out:
    slides/<name>.standalone.html
"""
import base64
import mimetypes
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
NAME = sys.argv[1] if len(sys.argv) > 1 else 'intro_deck.html'
SRC = (HERE / NAME) if not pathlib.Path(NAME).is_absolute() else pathlib.Path(NAME)
OUT = SRC.with_suffix('').with_suffix('.standalone.html')


def inline(match):
    src = match.group(1)
    if src.startswith(('data:', 'http:', 'https:')):
        return match.group(0)
    path = (SRC.parent / src).resolve()
    if not path.is_file():
        raise FileNotFoundError(f'referenced asset not found: {path}')
    mime = mimetypes.guess_type(path.name)[0] or 'application/octet-stream'
    b64 = base64.b64encode(path.read_bytes()).decode('ascii')
    return f'src="data:{mime};base64,{b64}"'


html = SRC.read_text(encoding='utf-8')
html = re.sub(r'src="([^"]+)"', inline, html)
OUT.write_text(html, encoding='utf-8')
print(f'wrote {OUT}  ({OUT.stat().st_size / 1024:.0f} kB)')
