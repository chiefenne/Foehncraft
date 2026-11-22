"""An assistant session through the MCP server (L18) as a figure: the
calls of a scripted session on the finned heat sink, the answers of
`foehn mcp` (a refusal with its fix, the person's approval, the result
summary with its limitations, the residual history paged by cursor), and
the residual history it read. The gallery image of the MCP server.

    python3 -m tools.plot.mcp_session docs/figures/gallery/mcp_session.svg [FOEHN]

A scripted client stands in for the assistant: it sends the calls of
examples/cases/tools/mcp/heat_sink_session.jsonl to `foehn mcp`, answers
the approval question as the line says, and keeps every message. The run
writes output/tools/mcp/; the journal goes to
output/tools/mcp-journal.jsonl, which `foehn run --script` replays. FOEHN
defaults to target/release/foehn.
"""
import json
import math
import os
import shutil
import subprocess
import sys
from xml.sax.saxutils import escape

SCRIPT = "examples/cases/tools/mcp/heat_sink_session.jsonl"
POLICY = "examples/cases/tools/mcp/foehn-mcp.toml"
OUT = "output/tools/mcp"
JOURNAL = "output/tools/mcp-journal.jsonl"


def session(foehn):
    shutil.rmtree(OUT, ignore_errors=True)
    if os.path.exists(JOURNAL):
        os.remove(JOURNAL)
    os.makedirs(OUT)
    p = subprocess.Popen(
        [foehn, "--quiet", "--journal", JOURNAL, "mcp", "--read-root", "examples/cases",
         "--read-root", OUT, "--write-root", OUT, "--policy", POLICY],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
    n = [0]

    def send(m):
        p.stdin.write(json.dumps(m) + "\n")
        p.stdin.flush()

    def request(method, params, answer=None, asked=None):
        n[0] += 1
        send({"jsonrpc": "2.0", "id": n[0], "method": method, "params": params})
        while True:
            m = json.loads(p.stdout.readline())
            if m.get("method") == "elicitation/create":
                asked.append(m["params"]["message"])
                send({"jsonrpc": "2.0", "id": m["id"],
                      "result": {"action": "accept", "content": {"answer": answer or "refuse-once"}}})
            elif m.get("id") == n[0]:
                return m

    request("initialize", {"protocolVersion": "2025-06-18", "capabilities": {"elicitation": {}},
                           "clientInfo": {"name": "scripted-assistant", "version": "1"}})
    send({"jsonrpc": "2.0", "method": "notifications/initialized"})
    steps = []
    for line in open(SCRIPT):
        c = json.loads(line)
        asked = []
        args = dict(c["arguments"])
        while True:
            r = request("tools/call", {"name": c["tool"], "arguments": args}, c.get("answer"), asked)
            env = r["result"]["structuredContent"]
            steps.append({"tool": c["tool"], "args": dict(args), "asked": asked,
                          "answer": c.get("answer"), "env": env,
                          "headline": r["result"]["content"][0]["text"].split("\n")[0]})
            print(steps[-1]["headline"][:160])
            # the assistant follows the cursor to the last page
            nxt = (env["result"] or {}).get("next_cursor")
            if not (c.get("follow_cursor") and nxt):
                break
            args["cursor"] = nxt
            asked = []
    p.stdin.close()
    p.wait()
    return steps


def wrap(text, width):
    out, line = [], ""
    for w in text.split():
        if len(line) + len(w) + 1 > width and line:
            out.append(line)
            line = w
        else:
            line = f"{line} {w}".strip()
    if line:
        out.append(line)
    return out


def bubbles(steps):
    """The transcript: (who, text) in order."""
    b = []
    for s in steps:
        a = s["args"]
        extra = ""
        if a.get("set"):
            extra = " with set " + ", ".join(a["set"])
        if a.get("series"):
            extra = f" series {a['series']}, limit {a.get('limit', 200)}"
            if a.get("cursor"):
                extra += f", cursor {a['cursor']}"
        b.append(("assistant", f"{s['tool']}{extra}"))
        for q in s["asked"]:
            b.append(("server asks", q.replace("The assistant wants to call", "Call")))
            b.append(("person", s["answer"] or "refuse-once"))
        env = s["env"]
        res = env["result"] or {}
        if env["diagnostics"]:
            d = env["diagnostics"][0]
            fix = d.get("fix") or {}
            t = f"{d['code']}: {d['message']}"
            if fix.get("set"):
                t += " Fix: set " + ", ".join(fix["set"])
            b.append(("server", t))
        elif s["tool"] == "case.check":
            e = res["estimate"]
            b.append(("server", f"case can be run: {e['cells']} cells, up to {e['steps'][0]} iterations, "
                                f"host memory {e['host_memory'] / 1e6:.0f} MB"))
        elif s["tool"] == "flow.run":
            sv = res["solves"][0]
            lims = res.get("limitations", [])
            t = (f"{env['status']} (exit {env['exit_code']}): {sv['iterations']} iterations, "
                 f"converged {str(sv['converged']).lower()}; {len(lims)} limitation(s)")
            if lims:
                t += ": " + "; ".join(f"{x['code']} {x['message']}" for x in lims)
            b.append(("server", t))
        elif s["tool"] == "result.page" and "available" in res:
            names = [f"{x['name']} ({x['rows']})" for x in res["available"]]
            b.append(("server", "series: " + ", ".join(names)))
        elif s["tool"] == "result.page":
            b.append(("server", f"{res['series']}: rows {res['first']} to "
                                f"{res['first'] + len(res['rows'])} of {res['total']}, "
                                f"next cursor {res['next_cursor'] or 'none (last page)'}"))
    return b


COLOURS = {"assistant": "#1f6fb4", "server": "#4a4a4a", "server asks": "#b86e00", "person": "#2e7d32"}


def chart(steps, dst):
    w, tw, gap = 1000, 560, 30
    out = []
    y = 56
    for who, text in bubbles(steps):
        lines = wrap(text, 76)[:5]
        hgt = 22 + 15 * len(lines)
        x = 20 if who == "assistant" or who == "person" else 60
        colour = COLOURS[who]
        out.append(f'<rect x="{x}" y="{y}" width="{tw - 40}" height="{hgt}" rx="6" '
                   f'fill="{colour}" fill-opacity="0.09" stroke="{colour}"/>')
        out.append(f'<text x="{x + 10}" y="{y + 15}" font-weight="bold" fill="{colour}">{escape(who)}</text>')
        for i, ln in enumerate(lines):
            out.append(f'<text x="{x + 10}" y="{y + 31 + 15 * i}" font-family="Menlo, Consolas, monospace" '
                       f'font-size="10.5">{escape(ln)}</text>')
        y += hgt + 8
    h = max(y + 20, 470)
    # the residual history the assistant paged
    pages = [s["env"]["result"] for s in steps if s["tool"] == "result.page" and s["args"].get("series")]
    page = pages[0]
    cols = page["columns"]
    rows = [r for p in pages for r in p["rows"]]
    px, py, pw, ph = tw + gap + 50, 70, w - tw - gap - 80, 300
    it = cols.index("iteration")
    names = [c for c in ("continuity", "momentum_x", "momentum_y", "energy", "k", "omega") if c in cols]
    palette = ["#1f6fb4", "#d1495b", "#edae49", "#2e9e5b", "#7a5195", "#00798c"]
    vals = [r[cols.index(c)] for r in rows for c in names if isinstance(r[cols.index(c)], float) and r[cols.index(c)] > 0]
    vals.sort()
    lo = 10 ** math.floor(math.log10(vals[len(vals) // 100]))
    hi = 10 ** math.ceil(math.log10(max(vals)))
    n = rows[-1][it]
    sx = lambda v: px + v / n * pw  # noqa: E731
    sy = lambda v: py + (math.log10(hi) - math.log10(max(v, lo))) / (math.log10(hi) - math.log10(lo)) * ph  # noqa: E731
    out.append(f'<text x="{px + pw / 2}" y="{py - 26}" text-anchor="middle" font-size="13">'
               f'Residuals read with result.page: {page["total"]} rows in {len(pages)} pages</text>')
    k = math.log10(lo)
    while k <= math.log10(hi) + 1e-9:
        yy = sy(10 ** k)
        out.append(f'<line x1="{px}" y1="{yy:.1f}" x2="{px + pw}" y2="{yy:.1f}" stroke="#e6e6e6"/>')
        out.append(f'<text x="{px - 6}" y="{yy + 4:.1f}" text-anchor="end" fill="#666">1e{int(k)}</text>')
        k += 1
    for t in range(0, int(n) + 1, max(int(n) // 4 // 100 * 100, 100)):
        out.append(f'<text x="{sx(t):.1f}" y="{py + ph + 16}" text-anchor="middle" fill="#666">{t}</text>')
    out.append(f'<text x="{px + pw / 2}" y="{py + ph + 34}" text-anchor="middle" fill="#444">outer iteration</text>')
    out.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" fill="none" stroke="#999"/>')
    for j, c in enumerate(names):
        ci = cols.index(c)
        pts = " ".join(f"{sx(r[it]):.1f},{sy(r[ci]):.1f}" for r in rows if isinstance(r[ci], float) and r[ci] > 0)
        out.append(f'<polyline points="{pts}" fill="none" stroke="{palette[j]}" stroke-width="1.4"/>')
        out.append(f'<rect x="{px + 10 + (j % 3) * 120}" y="{py + ph + 46 + (j // 3) * 18}" width="12" height="12" fill="{palette[j]}"/>')
        out.append(f'<text x="{px + 28 + (j % 3) * 120}" y="{py + ph + 56 + (j // 3) * 18}">{escape(c)}</text>')
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            'font-family="Helvetica, Arial, sans-serif" font-size="12">',
            f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
            f'<text x="{tw / 2}" y="30" text-anchor="middle" font-size="15">'
            'An assistant session through foehn mcp: the finned heat sink</text>']
    with open(dst, "w") as f:
        f.write("\n".join(head + out + ["</svg>"]) + "\n")
    print(f"wrote {dst} ({w} x {h})")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    foehn = sys.argv[2] if len(sys.argv) > 2 else "target/release/foehn"
    chart(session(foehn), sys.argv[1])


if __name__ == "__main__":
    main()
