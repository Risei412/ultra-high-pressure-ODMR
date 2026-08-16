"""
build_standalone.py
Inline every <img src="..."> of intro_deck.html as a data: URI so the deck is a
single self-contained file (needed for hosting, e-mail, or offline viewing).

Run from the repository root or from slides/:
    python slides/build_standalone.py
Out:
    slides/intro_deck.standalone.html
"""
import base64
import mimetypes
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / 'intro_deck.html'
OUT = HERE / 'intro_deck.standalone.html'


def inline(match):
    src = match.group(1)
    if src.startswith(('data:', 'http:', 'https:')):
        return match.group(0)
    path = (HERE / src).resolve()
    if not path.is_file():
        raise FileNotFoundError(f'referenced asset not found: {path}')
    mime = mimetypes.guess_type(path.name)[0] or 'application/octet-stream'
    b64 = base64.b64encode(path.read_bytes()).decode('ascii')
    return f'src="data:{mime};base64,{b64}"'


html = SRC.read_text(encoding='utf-8')
html = re.sub(r'src="([^"]+)"', inline, html)
OUT.write_text(html, encoding='utf-8')
print(f'wrote {OUT}  ({OUT.stat().st_size / 1024:.0f} kB)')
