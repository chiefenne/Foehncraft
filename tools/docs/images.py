"""Gives every image of the gallery pages and the README its display size.

    python -m tools.docs.images           rewrite the pages
    python -m tools.docs.images --check   exit 1 when an image has no size

Why: a link to an entry far down a gallery page (gallery.md#mrf opens
gallery/flow.md#mrf) makes the browser jump to the entry while the images
above it are still loading. An image without a size takes no room until
it has arrived, so every image that arrives after the jump pushes the
entry further down, and the reader lands thousands of pixels above it.
With width and height in the page, the browser reserves the room from the
start and the jump lands on the entry.

The script turns ``![alt](src)`` into
``<img src="src" alt="alt" width="W" height="H">`` and adds the size to an
``<img>`` that has none. The size is the size of the file (PNG or SVG),
scaled down to at most MAX_WIDTH pixels (MAX_WIDTH_TABLE inside a table
row), with the aspect ratio kept. Images inside code blocks and remote
images are left alone.
"""

from __future__ import annotations

import re
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAGES = (sorted((ROOT / "docs" / "gallery").glob("*.md"))
         + sorted((ROOT / "docs" / "applications").glob("*.md")) + [ROOT / "README.md"])

MAX_WIDTH = 800
MAX_WIDTH_TABLE = 200

MD_IMAGE = re.compile(r"!\[([^\]]*)\]\(([^)\s]+)\)")
HTML_IMAGE = re.compile(r"<img\s[^>]*>")
ATTR = re.compile(r'(\w+)="([^"]*)"')
FENCE = re.compile(r"^\s*(```|~~~)")


def svg_length(value: str) -> float:
    m = re.match(r"\s*([0-9.]+)\s*(px|pt)?\s*$", value)
    if not m:
        raise ValueError(value)
    return float(m.group(1)) * (4.0 / 3.0 if m.group(2) == "pt" else 1.0)


def file_size(path: Path) -> tuple[float, float]:
    """(width, height) in CSS pixels of a PNG or SVG file."""
    data = path.read_bytes()
    if data[:8] == b"\x89PNG\r\n\x1a\n":
        w, h = struct.unpack(">II", data[16:24])
        return float(w), float(h)
    head = data[:4000].decode("utf-8", "replace")
    tag = re.search(r"<svg\b[^>]*>", head, re.S)
    if tag:
        attrs = dict(ATTR.findall(tag.group(0)))
        if "width" in attrs and "height" in attrs and "%" not in attrs["width"]:
            return svg_length(attrs["width"]), svg_length(attrs["height"])
        if "viewBox" in attrs:
            _, _, w, h = (float(v) for v in attrs["viewBox"].replace(",", " ").split())
            return w, h
    raise SystemExit(f"{path}: cannot read the image size (PNG or SVG with width/height or viewBox)")


def display(path: Path, in_table: bool) -> tuple[int, int]:
    w, h = file_size(path)
    cap = MAX_WIDTH_TABLE if in_table else MAX_WIDTH
    s = min(1.0, cap / w)
    return round(w * s), round(h * s)


def esc(text: str) -> str:
    return (
        text.replace("&", "&amp;").replace('"', "&quot;").replace("<", "&lt;").replace(">", "&gt;")
    )


def fix_line(line: str, page: Path) -> tuple[str, list[str]]:
    in_table = line.lstrip().startswith("|")
    problems: list[str] = []

    def local(src: str) -> Path | None:
        if re.match(r"^[a-z]+:", src) or src.startswith("/"):
            return None
        return (page.parent / src.split("#")[0]).resolve()

    def md(m: re.Match) -> str:
        alt, src = m.group(1), m.group(2)
        p = local(src)
        if p is None:
            return m.group(0)
        if not p.exists():
            problems.append(f"missing file {src}")
            return m.group(0)
        w, h = display(p, in_table)
        problems.append(f"no size: {src}")
        return f'<img src="{src}" alt="{esc(alt)}" width="{w}" height="{h}">'

    def html(m: re.Match) -> str:
        tag = m.group(0)
        attrs = dict(ATTR.findall(tag))
        if "width" in attrs and "height" in attrs:
            return tag
        src = attrs.get("src", "")
        p = local(src)
        if p is None or not p.exists():
            return tag
        problems.append(f"no size: {src}")
        if "width" in attrs:
            # keep the chosen width, add the matching height
            fw, fh = file_size(p)
            w = int(attrs["width"])
            h = round(fh * w / fw)
            return tag[:-1].rstrip(" /") + f' height="{h}">'
        w, h = display(p, in_table)
        return tag[:-1].rstrip(" /") + f' width="{w}" height="{h}">'

    line = HTML_IMAGE.sub(html, line)
    line = MD_IMAGE.sub(md, line)
    return line, problems


def fix_page(page: Path) -> tuple[str, list[str]]:
    out, problems, fenced = [], [], False
    for n, line in enumerate(page.read_text(encoding="utf-8").split("\n"), 1):
        if FENCE.match(line):
            fenced = not fenced
        if not fenced:
            line, p = fix_line(line, page)
            problems += [f"{page.relative_to(ROOT)}:{n}: {x}" for x in p]
        out.append(line)
    return "\n".join(out), problems


def main(argv: list[str]) -> int:
    check = "--check" in argv
    all_problems = []
    for page in PAGES:
        new, problems = fix_page(page)
        all_problems += problems
        if not check and new != page.read_text(encoding="utf-8"):
            page.write_text(new, encoding="utf-8")
    if check and all_problems:
        print("\n".join(all_problems))
        print("images without a size (run python -m tools.docs.images)")
        return 1
    print(f"sized {len(all_problems)} images" if all_problems else "every image has a size")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
