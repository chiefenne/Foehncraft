"""The diagnosis of a steady run (D198 (d), gallery 166): the residual
histories against the outer iteration up to the iteration budget, the
window of the diagnosis, the fitted trends and the verdict.

    python3 -m tools.plot.diagnosis_chart CHART OUTDIR

OUTDIR is the output directory of a steady flow run with a diagnosis
(its residuals.csv and run.jsonl). No dependencies beyond the standard
library.
"""
import csv
import json
import math
import os
import sys

COLOURS = {"momentum_x": "#d0582f", "momentum_y": "#2f6fd0", "continuity": "#2e9a5a"}
EXTRA = ["#8a4fc8", "#b08a1e", "#1e9fb0"]


def read(d):
    with open(os.path.join(d, "residuals.csv")) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    diag = budget = None
    with open(os.path.join(d, "run.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            if r.get("type") == "diagnosis":
                diag = r
            elif r.get("type") == "budget":
                budget = r
    if diag is None:
        sys.exit(f"{d}: no diagnosis row in run.jsonl")
    return rows, diag, budget


def main(argv):
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    dst, d = argv
    rows, diag, budget = read(d)
    names = [e["equation"] for e in diag["equations"]]
    it = [int(r["iteration"]) for r in rows]
    series = {n: [float(r[n]) for r in rows] for n in names if n in rows[0]}
    vals = [v for s in series.values() for v in s if v > 0]
    lo = math.floor(math.log10(min(vals)))
    hi = max(math.ceil(math.log10(max(vals))), lo + 1)
    n = diag["budget"]
    x0, y0, w, h = 80, 46, 560, 320
    fx = lambda i: x0 + (i - 1) / max(n - 1, 1) * w
    fy = lambda v: y0 + (hi - math.log10(v)) / (hi - lo) * h
    at, win = diag["iteration"], diag["window"]
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{x0 + w + 200}" height="{y0 + h + 60}" '
         f'viewBox="0 0 {x0 + w + 200} {y0 + h + 60}" font-family="sans-serif" font-size="12">',
         f'<rect width="100%" height="100%" fill="white"/>']
    # window of the diagnosis and the remaining budget
    o.append(f'<rect x="{fx(at - win + 1):.1f}" y="{y0}" width="{fx(at) - fx(at - win + 1):.1f}" height="{h}" fill="#f2e6c9"/>')
    o.append(f'<rect x="{fx(at):.1f}" y="{y0}" width="{fx(n) - fx(at):.1f}" height="{h}" fill="#eeeeee"/>')
    for k in range(lo, hi + 1):
        y = fy(10.0 ** k)
        o.append(f'<line x1="{x0}" x2="{x0 + w}" y1="{y:.1f}" y2="{y:.1f}" stroke="#ddd"/>')
        o.append(f'<text x="{x0 - 8}" y="{y + 4:.1f}" text-anchor="end">1e{k}</text>')
    step = 500
    for i in range(0, n + 1, step):
        x = fx(max(i, 1))
        o.append(f'<text x="{x:.1f}" y="{y0 + h + 18}" text-anchor="middle">{i}</text>')
    o.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="none" stroke="#444"/>')
    extra = iter(EXTRA)
    legend_y = y0 + 10
    for name, s in series.items():
        c = COLOURS.get(name) or next(extra, "#555")
        pts = " ".join(f"{fx(i):.1f},{fy(v):.1f}" for i, v in zip(it, s) if v > 0)
        o.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="1.4"/>')
        o.append(f'<line x1="{x0 + w + 16}" x2="{x0 + w + 36}" y1="{legend_y}" y2="{legend_y}" stroke="{c}" stroke-width="2"/>')
        o.append(f'<text x="{x0 + w + 42}" y="{legend_y + 4}">{name}</text>')
        legend_y += 20
    # fitted trends over the window, extended to the budget
    for e in diag["equations"]:
        s = e["decades_per_1000"] / 1000.0
        r = e["residual"]
        a, b = at - win + 1, n
        ya, yb = r * 10 ** (s * (a - at)), r * 10 ** (s * (b - at))
        if all(10.0 ** lo <= v <= 10.0 ** hi for v in (ya, yb)):
            o.append(f'<line x1="{fx(a):.1f}" x2="{fx(b):.1f}" y1="{fy(ya):.1f}" y2="{fy(yb):.1f}" '
                     f'stroke="#222" stroke-width="1" stroke-dasharray="5,4"/>')
    o.append(f'<line x1="{fx(at):.1f}" x2="{fx(at):.1f}" y1="{y0}" y2="{y0 + h}" stroke="#a33" stroke-width="2"/>')
    verdict = diag["verdict"].upper()
    o.append(f'<text x="{fx(at) - 6:.1f}" y="{y0 + h - 10}" text-anchor="end" fill="#a33">diagnosis at {at}: {verdict}</text>')
    stopped = budget and budget.get("iterations", {}).get("diagnosis") == "stop"
    o.append(f'<text x="{(fx(at) + fx(n)) / 2:.1f}" y="{y0 + 18}" text-anchor="middle" fill="#666">'
             f'{"not run (stopped)" if stopped else "rest of the budget"}</text>')
    o.append(f'<text x="{(fx(at - win + 1) + fx(at)) / 2:.1f}" y="{y0 + 18}" text-anchor="middle" fill="#8a6a1e">window {win}</text>')
    o.append(f'<text x="{x0 + w / 2}" y="{y0 + h + 40}" text-anchor="middle">outer iteration (budget {n})</text>')
    o.append(f'<text x="{x0}" y="24" font-size="14">Scaled residuals; dashed: trend over the window, extended to the budget</text>')
    o.append(f'<line x1="{x0 + w + 16}" x2="{x0 + w + 36}" y1="{legend_y}" y2="{legend_y}" stroke="#222" stroke-dasharray="5,4"/>')
    o.append(f'<text x="{x0 + w + 42}" y="{legend_y + 4}">trend</text>')
    o.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(o) + "\n")
    print(f"{dst}: {verdict} at {at} of {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
