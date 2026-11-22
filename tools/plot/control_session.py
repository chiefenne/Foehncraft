"""An assistant steers a running case through the live control port (L19) as
a figure: the calls of a scripted session on the transient lid-driven cavity
(`run.status` while `flow.run` is going, `run.set` that is applied, one that
the case refuses, `run.stop`), and the continuity residual of the run with
the three events marked. The gallery image of the control port.

    python3 -m tools.plot.control_session docs/figures/gallery/control_port.svg [FOEHN]

A scripted client stands in for the assistant: it starts `flow.run`, does
not wait for it, and sends the control calls on the same connection while
the run is going. The run writes output/tools/control/live; the journal goes
to output/tools/control-journal.jsonl. FOEHN defaults to target/release/foehn.
"""
import json
import math
import os
import shutil
import subprocess
import sys
import time
from xml.sax.saxutils import escape

CASE = "examples/cases/flow/cavity_re100_transient.toml"
OUT = "output/tools/control"
RUN = OUT + "/live"
JOURNAL = "output/tools/control-journal.jsonl"
POLICY = OUT + "-policy/foehn-mcp.toml"


class Client:
    def __init__(self, foehn):
        shutil.rmtree(OUT, ignore_errors=True)
        shutil.rmtree(os.path.dirname(POLICY), ignore_errors=True)
        if os.path.exists(JOURNAL):
            os.remove(JOURNAL)
        os.makedirs(OUT)
        os.makedirs(os.path.dirname(POLICY))
        with open(POLICY, "w") as f:
            f.write('[actions]\nstart-run = "allow"\n')
        self.p = subprocess.Popen(
            [foehn, "--quiet", "--journal", JOURNAL, "mcp", "--read-root", "examples/cases",
             "--read-root", OUT, "--write-root", OUT, "--policy", POLICY],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
        self.n = 0
        self.inbox = {}
        self.request("initialize", {"protocolVersion": "2025-06-18", "capabilities": {},
                                    "clientInfo": {"name": "scripted-assistant", "version": "1"}})
        self.send({"jsonrpc": "2.0", "method": "notifications/initialized"})

    def send(self, m):
        self.p.stdin.write(json.dumps(m) + "\n")
        self.p.stdin.flush()

    def post(self, method, params):
        self.n += 1
        self.send({"jsonrpc": "2.0", "id": self.n, "method": method, "params": params})
        return self.n

    def wait(self, i):
        while i not in self.inbox:
            m = json.loads(self.p.stdout.readline())
            if "id" in m and "method" not in m:
                self.inbox[m["id"]] = m
        return self.inbox.pop(i)

    def request(self, method, params):
        return self.wait(self.post(method, params))

    def call(self, tool, args):
        r = self.wait(self.post("tools/call", {"name": tool, "arguments": args}))
        return r["result"]["structuredContent"], r["result"]["content"][0]["text"].split("\n")[0]


def session(foehn):
    c = Client(foehn)
    events, talk = [], []
    run = c.post("tools/call", {"name": "flow.run", "arguments": {
        "case": CASE, "out": RUN,
        "set": ["time.time_step=0.005", "time.steps=1000000",
                "run.limits.max_wall_time=120", "run.limits.max_output=2e9"]}})
    talk.append(("assistant", "flow.run cavity_re100_transient.toml, time step 0.005 s, "
                              "run limits 120 s and 2 GB (does not wait for the result)"))

    def status(wait_for=None):
        for _ in range(400):
            if os.path.exists(RUN + "/control.json"):
                env, head = c.call("run.status", {"dir": RUN})
                r = env["result"]
                if env["status"] == "finished" and r["last"] and (wait_for is None or wait_for(r)):
                    return r
            time.sleep(0.05)
        sys.exit("the run did not come up")

    r = status(lambda r: r["last"]["step"] > 400)
    talk.append(("assistant", "run.status"))
    talk.append(("server", f"running for {r['wall_time']:.1f} s, step {r['last']['step']}, "
                           f"t = {r['last']['time']:.3f} s, continuity {r['last']['continuity']:.2e}"))
    # the time step: applied at the end of the current step
    env, _ = c.call("run.set", {"dir": RUN, "set": {"time.time_step": 0.01}})
    talk.append(("assistant", "run.set time.time_step = 0.01"))
    talk.append(("server", f"queued {', '.join(env['result']['queued'])}; applied at the end of the current step"))
    # a setting the case refuses: relaxation below 1 with one outer iteration
    c.call("run.set", {"dir": RUN, "set": {"solver.momentum.relaxation": 0.7}})
    talk.append(("assistant", "run.set solver.momentum.relaxation = 0.7"))

    def settled(r):
        done = {e["key"]: e for e in r["control"]}
        return len(done) == 2 and r["last"]["step"] > 0
    r = status(settled)
    for e in r["control"]:
        if e["status"] == "applied":
            talk.append(("server", f"step {e['at']}: {e['key']} {e['old']} -> {e['new']} applied "
                                   f"(equivalent: {e['schedule']})"))
        else:
            talk.append(("server", f"step {e['at']}: {e['key']} = {e['value']} not applied: "
                                   + e["error"].split(" (")[0]))
        events.append(e)
    r = status(lambda r: r["last"]["step"] > events[-1]["at"] + 600)
    talk.append(("assistant", "run.status"))
    talk.append(("server", f"step {r['last']['step']}, continuity {r['last']['continuity']:.2e}; "
                           "the new time step is in force"))
    c.call("run.stop", {"dir": RUN})
    talk.append(("assistant", "run.stop"))
    res = c.wait(run)["result"]["structuredContent"]
    s = res["result"]
    talk.append(("server", f"flow.run ends: {res['status']} (exit {res['exit_code']}), "
                           f"{s['stop']['code']} {s['stop']['message']}; restart file and results written; "
                           f"{len(s['limitations'])} limitation(s)"))
    events.append({"status": "stop", "at": int(s["solves"][0]["steps"] if "steps" in s["solves"][0] else s["solves"][0]["iterations"])})
    c.p.stdin.close()
    c.p.wait()
    return talk, events


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


COLOURS = {"assistant": "#1f6fb4", "server": "#4a4a4a"}


def chart(talk, events, dst):
    w, tw, gap = 1000, 520, 30
    out, y = [], 56
    for who, text in talk:
        lines = wrap(text, 70)[:5]
        hgt = 22 + 15 * len(lines)
        x = 20 if who == "assistant" else 50
        colour = COLOURS[who]
        out.append(f'<rect x="{x}" y="{y}" width="{tw - 70}" height="{hgt}" rx="6" '
                   f'fill="{colour}" fill-opacity="0.09" stroke="{colour}"/>')
        out.append(f'<text x="{x + 10}" y="{y + 15}" font-weight="bold" fill="{colour}">{who}</text>')
        for i, ln in enumerate(lines):
            out.append(f'<text x="{x + 10}" y="{y + 31 + 15 * i}" font-family="Menlo, Consolas, monospace" '
                       f'font-size="10.5">{escape(ln)}</text>')
        y += hgt + 8
    h = max(y + 20, 420)
    rows = []
    for line in open(RUN + "/history.csv"):
        if line.startswith("#") or line.startswith("step"):
            continue
        rows.append(line.strip().split(","))
    head = [l for l in open(RUN + "/history.csv") if not l.startswith("#")][0].strip().split(",")
    si, ci = head.index("step"), head.index("continuity")
    data = [(int(r[si]), float(r[ci])) for r in rows if float(r[ci]) > 0]
    data = data[::max(len(data) // 800, 1)]
    px, py, pw, ph = tw + gap + 50, 70, w - tw - gap - 80, 260
    n = data[-1][0]
    lo = 10 ** math.floor(math.log10(min(v for _, v in data)))
    hi = 10 ** math.ceil(math.log10(max(v for _, v in data)))
    sx = lambda v: px + v / n * pw  # noqa: E731
    sy = lambda v: py + (math.log10(hi) - math.log10(max(v, lo))) / (math.log10(hi) - math.log10(lo)) * ph  # noqa: E731
    out.append(f'<text x="{px + pw / 2}" y="{py - 26}" text-anchor="middle" font-size="13">'
               f'Continuity residual of the steered run</text>')
    k = math.log10(lo)
    while k <= math.log10(hi) + 1e-9:
        yy = sy(10 ** k)
        out.append(f'<line x1="{px}" y1="{yy:.1f}" x2="{px + pw}" y2="{yy:.1f}" stroke="#e6e6e6"/>')
        out.append(f'<text x="{px - 6}" y="{yy + 4:.1f}" text-anchor="end" fill="#666">1e{int(k)}</text>')
        k += 1
    step = max(int(n) // 4 // 100 * 100, 100)
    for t in range(0, int(n) + 1, step):
        out.append(f'<text x="{sx(t):.1f}" y="{py + ph + 16}" text-anchor="middle" fill="#666">{t}</text>')
    out.append(f'<text x="{px + pw / 2}" y="{py + ph + 34}" text-anchor="middle" fill="#444">time step</text>')
    out.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" fill="none" stroke="#999"/>')
    pts = " ".join(f"{sx(s):.1f},{sy(v):.1f}" for s, v in data)
    out.append(f'<polyline points="{pts}" fill="none" stroke="#1f6fb4" stroke-width="1.5"/>')
    marks = {"applied": ("#2e7d32", "time step 0.005 to 0.01"),
             "rejected": ("#b86e00", "relaxation 0.7 refused"),
             "stop": ("#b3261e", "run.stop")}
    for j, e in enumerate(events):
        colour, label = marks[e["status"]]
        xx = sx(e["at"])
        out.append(f'<line x1="{xx:.1f}" y1="{py}" x2="{xx:.1f}" y2="{py + ph}" stroke="{colour}" '
                   f'stroke-dasharray="4 3" stroke-width="1.4"/>')
        out.append(f'<rect x="{px + 10}" y="{py + ph + 46 + 18 * j}" width="12" height="12" fill="{colour}"/>')
        out.append(f'<text x="{px + 28}" y="{py + ph + 56 + 18 * j}">step {e["at"]}: {escape(label)}</text>')
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            'font-family="Helvetica, Arial, sans-serif" font-size="12">',
            f'<rect width="{w}" height="{h}" fill="#ffffff"/>',
            f'<text x="{tw / 2}" y="30" text-anchor="middle" font-size="15">'
            'An assistant steers a running case: the cavity at Re = 100</text>']
    with open(dst, "w") as f:
        f.write("\n".join(head + out + ["</svg>"]) + "\n")
    print(f"wrote {dst} ({w} x {h})")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    foehn = sys.argv[2] if len(sys.argv) > 2 else "target/release/foehn"
    talk, events = session(foehn)
    chart(talk, events, sys.argv[1])


if __name__ == "__main__":
    main()
