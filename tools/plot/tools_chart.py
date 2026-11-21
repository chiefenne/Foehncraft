"""The tool catalogue of `foehn commands` (D196) as a chart: the tools of
each category, coloured by action kind (read, write files, start a run),
each bar labelled with its tool ids. The gallery image of the registry.

    python3 -m tools.plot.tools_chart docs/figures/gallery/tool_catalogue.svg [FOEHN]

FOEHN defaults to target/release/foehn.
"""
import json
import subprocess
import sys
from xml.sax.saxutils import escape

COLOUR = {"read": "#1f6fb4", "write-files": "#d9822b", "start-run": "#c23b3b"}
LABEL = {"read": "read", "write-files": "writes files", "start-run": "starts a run"}


def catalogue(foehn):
    out = subprocess.run([foehn, "commands"], capture_output=True, check=True, text=True).stdout
    return json.loads(out)


def chart(cat, dst):
    cats = []
    for t in cat["tools"]:
        if t["category"] not in cats:
            cats.append(t["category"])
    rows = {c: [t for t in cat["tools"] if t["category"] == c] for c in cats}
    most = max(len(v) for v in rows.values())
    w, left, right, top, bar, gap = 760, 110, 20, 64, 30, 14
    cell = (w - left - right) / most
    h = top + len(cats) * (bar + gap) + 56
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        'font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
        f'<text x="{w / 2}" y="26" text-anchor="middle" font-size="15">'
        f'The tool registry: {len(cat["tools"])} tools in {len(cats)} categories '
        f'(foehn commands, Foehncraft {escape(cat["foehncraft_version"])})</text>',
    ]
    for i, c in enumerate(cats):
        y = top + i * (bar + gap)
        out.append(f'<text x="{left - 10}" y="{y + bar / 2 + 4}" text-anchor="end" '
                   f'font-size="13">{escape(c)}</text>')
        for j, t in enumerate(rows[c]):
            x = left + j * cell
            verb = t["id"].split(".", 1)[1]
            out.append(f'<rect x="{x + 1:.1f}" y="{y}" width="{cell - 2:.1f}" height="{bar}" rx="3" '
                       f'fill="{COLOUR[t["action"]]}"><title>{escape(t["id"])}: '
                       f'{escape(t["summary"])}</title></rect>')
            out.append(f'<text x="{x + cell / 2:.1f}" y="{y + bar / 2 + 4}" text-anchor="middle" '
                       f'fill="#ffffff" font-size="12">{escape(verb)}</text>')
    ly = h - 26
    x = left
    for k in ("read", "write-files", "start-run"):
        n = sum(1 for t in cat["tools"] if t["action"] == k)
        out.append(f'<rect x="{x}" y="{ly - 11}" width="14" height="14" rx="2" fill="{COLOUR[k]}"/>')
        out.append(f'<text x="{x + 20}" y="{ly}">{LABEL[k]} ({n})</text>')
        x += 170
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    return w, h


def main():
    dst = sys.argv[1]
    foehn = sys.argv[2] if len(sys.argv) > 2 else "target/release/foehn"
    cat = catalogue(foehn)
    w, h = chart(cat, dst)
    print(f"wrote {dst} ({w} x {h}): {len(cat['tools'])} tools")


if __name__ == "__main__":
    main()
