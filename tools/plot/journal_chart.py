"""A journaled session and its replay (D197) as a chart: the calls of the
session in order, each with the time it took when it was recorded and when
`foehn run --script` ran it again, and whether the replay ended as recorded.
The gallery image of the journal.

    python3 -m tools.plot.journal_chart docs/figures/gallery/journal_replay.svg [FOEHN]

Runs the script examples/cases/tools/cavity_session.jsonl with a journal
into output/tools/journal/, then replays the journal. FOEHN defaults to
target/release/foehn.
"""
import json
import math
import os
import shutil
import subprocess
import sys
from xml.sax.saxutils import escape

SCRIPT = "examples/cases/tools/cavity_session.jsonl"
OUT = "output/tools/journal"
RECORDED, REPLAYED = "#1f6fb4", "#2e9e5b"


def session(foehn):
    shutil.rmtree(OUT, ignore_errors=True)
    shutil.rmtree("output/tools/cavity", ignore_errors=True)
    os.makedirs(OUT)
    journal = os.path.join(OUT, "session.jsonl")
    # record: each call of the script as its own command, with a journal
    for c in (json.loads(line) for line in open(SCRIPT)):
        subprocess.run([foehn, "--quiet", "--journal", journal, "call", c["tool"], "--input", "-"],
                       input=json.dumps(c["input"]), text=True, check=True,
                       stdout=subprocess.DEVNULL)
    calls = [r for r in (json.loads(line) for line in open(journal)) if r["type"] == "call"]
    shutil.rmtree("output/tools/cavity", ignore_errors=True)
    replay = subprocess.run([foehn, "--json", "--quiet", "run", "--script", journal],
                            capture_output=True, text=True, check=True).stdout
    return calls, json.loads(replay)


def chart(calls, replay, dst):
    w, left, right, top, row = 760, 150, 130, 70, 64
    h = top + len(calls) * row + 60
    times = [c["seconds"] for c in calls] + [c.get("seconds", 0) for c in replay["calls"]]
    lo = 1e-4
    hi = max(max(times), 1.0) * 1.5
    span = w - left - right
    sx = lambda v: left + (math.log10(max(v, lo)) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * span  # noqa: E731
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        'font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
        f'<text x="{w / 2}" y="26" text-anchor="middle" font-size="15">'
        f'A journaled session and its replay: {replay["replayed"]} calls, '
        f'{replay["different"]} different (foehn run --script)</text>',
    ]
    # time axis
    ay = top - 14
    for k in range(int(math.log10(lo)), int(math.ceil(math.log10(hi))) + 1):
        x = sx(10.0**k)
        if x > w - right + 1:
            break
        label = f"{10.0**k:g} s"
        out.append(f'<line x1="{x:.1f}" y1="{ay}" x2="{x:.1f}" y2="{h - 50}" stroke="#e3e3e3"/>')
        out.append(f'<text x="{x:.1f}" y="{ay - 4}" text-anchor="middle" fill="#666">{label}</text>')
    for i, (c, r) in enumerate(zip(calls, replay["calls"])):
        y = top + i * row
        out.append(f'<text x="{left - 10}" y="{y + 20}" text-anchor="end" font-size="13">'
                   f'{c["seq"]}. {escape(c["tool"])}</text>')
        for j, (secs, colour) in enumerate(((c["seconds"], RECORDED), (r.get("seconds", 0), REPLAYED))):
            by = y + 6 + j * 22
            out.append(f'<rect x="{left}" y="{by}" width="{max(sx(secs) - left, 3):.1f}" height="18" rx="2" '
                       f'fill="{colour}"><title>{secs:.3f} s</title></rect>')
        mark = "as recorded" if r["same"] else "different"
        colour = "#2e7d32" if r["same"] else "#c23b3b"
        out.append(f'<text x="{w - right + 14}" y="{y + 20}" fill="{colour}">'
                   f'{escape(r["status"])}</text>')
        out.append(f'<text x="{w - right + 14}" y="{y + 38}" fill="#666">{mark}</text>')
    ly = h - 22
    for k, (name, colour) in enumerate((("recorded (journal)", RECORDED),
                                        ("replayed (run --script)", REPLAYED))):
        x = left + k * 220
        out.append(f'<rect x="{x}" y="{ly - 11}" width="14" height="14" rx="2" fill="{colour}"/>')
        out.append(f'<text x="{x + 20}" y="{ly}">{name}</text>')
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {dst} ({w} x {h}): {len(calls)} calls")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    foehn = sys.argv[2] if len(sys.argv) > 2 else "target/release/foehn"
    calls, replay = session(foehn)
    chart(calls, replay, sys.argv[1])


if __name__ == "__main__":
    main()
