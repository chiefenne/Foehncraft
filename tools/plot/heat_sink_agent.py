"""What the assistant did in the agent-driven heat sink study (gallery 181),
as two figures.

heat_sink_agent.png: the session, call by call. Every call the assistant
made through `foehn mcp` is a mark on its own time line (case checks and
estimates, the runs with their fin pitch, the status polls while a run goes,
the pages and result reads), the three refused calls with their codes
(`MCP002`, `MCP003` twice), and the cell-hours of every run against the
study budget the assistant set itself (`study_budget.jsonl`).

heat_sink_flow.png: the speed of the four designs on one colour scale and
the temperature along the top of the base the assistant read by page
(`result.page`, series `line_base_top`).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.heat_sink_agent summary TRANSCRIPT.jsonl SESSION.csv
    python3 -m tools.plot.heat_sink_agent figures SESSION.csv STUDY_DIR OUT_DIR [LEDGER_DIR]

`summary` condenses the transcript of `tools/agent/stdio_host.py` (every
call, its arguments and its answer) to one row per call; the recorded
summary is examples/cases/tools/mcp/heat_sink_study/session.csv. STUDY_DIR
holds the runs `pitch_<mm>` (output/agent/study); LEDGER_DIR the
`study_budget.jsonl` of the recorded session (default STUDY_DIR; the replay
writes its own ledger with its own cell-hours).
"""
import csv
import glob
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.cht_turbulent import outline, panel
from tools.plot.energy_fields import grid
from tools.plot.turbulence_fields import read

FIELDS = ["ticket", "tool", "start_s", "seconds", "outcome", "code", "pitch_mm", "note"]


def summary(transcript, out):
    rows = [json.loads(ln) for ln in open(transcript) if ln.strip()]
    t0 = min(r["started"] for r in rows)
    result = []
    for r in sorted(rows, key=lambda r: r["started"]):
        env = r["envelope"]
        diags = env.get("diagnostics") or []
        code = diags[0].get("code", "") if diags else ""
        m = re.search(r"pitch_(\d+)", json.dumps(r["arguments"]))
        outcome = "refused" if env.get("status") == "error" else "ok"
        note = r["arguments"].get("series", "")
        if r["tool"] == "run.status":
            note = "polled while the run goes"
        if code == "MCP002":
            note = "no run limits in the case"
        elif code == "MCP003":
            note = "outside the write root" if "/tmp/elsewhere" in json.dumps(r["arguments"]) else (
                "60 000 cells, cap 40 000")
        result.append({
            "ticket": r["ticket"], "tool": r["tool"], "start_s": f"{r['started'] - t0:.1f}",
            "seconds": f"{r['seconds']:.2f}", "outcome": outcome, "code": code,
            "pitch_mm": m.group(1) if m else "", "note": note})
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        w.writerows(result)
    print(f"wrote {out}: {len(result)} calls")


def session_rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def budget(study):
    spent = []
    with open(os.path.join(study, "study_budget.jsonl")) as f:
        for ln in f:
            d = json.loads(ln)
            spent.append((d["directory"], d["spend"]["cell_hours"]))
    return spent


def timeline(rows, study, out):
    lanes = ["case.check", "case.estimate", "flow.run", "run.status", "result.page", "result.read"]
    colour = {"case.check": "#6b7a8f", "case.estimate": "#8c6bb1", "flow.run": "#1f6fb4",
              "run.status": "#2a9d8f", "result.page": "#e9a23b", "result.read": "#d1495b"}
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(12.0, 7.4), constrained_layout=True,
                                 gridspec_kw={"height_ratios": [3.2, 1.35]})
    ax.set_yticks(range(len(lanes)), lanes)
    ax.invert_yaxis()
    for r in rows:
        y = lanes.index(r["tool"])
        t, d = float(r["start_s"]), max(float(r["seconds"]), 1.8)
        if r["outcome"] == "refused":
            ax.plot(t + 0.9, y, "X", color="#c1121f", ms=11, zorder=5)
            ypos = {"MCP002": 0.85, "outside the write root": 1.1}.get(r["code"] if r["code"] == "MCP002" else r["note"], 1.35)
            ax.annotate(f"{r['code']}: {r['note']}", (t + 0.9, y), xytext=(30, ypos), textcoords="data",
                        fontsize=8.5, color="#c1121f", va="center",
                        bbox={"fc": "white", "ec": "none", "pad": 1.5},
                        arrowprops={"arrowstyle": "-", "color": "#c1121f", "lw": 0.6})
        else:
            ax.barh(y, d, left=t, height=0.55, color=colour[r["tool"]], zorder=3,
                    edgecolor="white", linewidth=0.4)
            if r["tool"] == "flow.run":
                ax.text(t + d / 2, y, f"pitch {r['pitch_mm']} mm", ha="center", va="center",
                        fontsize=8.5, color="white", zorder=4)
    ax.set_xlabel("seconds since the first call")
    ax.grid(axis="x", alpha=0.3)
    ax.set_title("What the assistant did through foehn mcp: 27 calls, three refused (red crosses), "
                 "four runs of two minutes each, the histories read after the last run", fontsize=10)
    ax.set_xlim(-5, max(float(r["start_s"]) + float(r["seconds"]) for r in rows) + 8)
    spent = budget(study)
    cum = np.cumsum([c for _, c in spent])
    names = [re.sub("pitch_", "", n) + " mm" for n, _ in spent]
    bx.barh(range(len(spent)), [c for _, c in spent], left=np.concatenate([[0], cum[:-1]]),
            color="#1f6fb4", height=0.6, edgecolor="white")
    for k, (n, c) in enumerate(spent):
        bx.text(cum[k] - c / 2, k, f"{c:.0f}", ha="center", va="center", color="white", fontsize=8.5)
    bx.axvline(3500, color="#c1121f", ls="--", lw=1.2)
    bx.text(3490, 0, "study budget 3500 cell-hours (set by the assistant)", ha="right",
            va="center", fontsize=8.5, color="#c1121f")
    bx.set_yticks(range(len(spent)), names)
    bx.invert_yaxis()
    bx.set_xlim(0, 3700)
    bx.set_xlabel(f"cell-hours spent (ledger study_budget.jsonl); {cum[-1]:.0f} of 3500 used")
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


def flow(study, out):
    runs = sorted((r for r in glob.glob(os.path.join(study, "pitch_*"))
                   if re.fullmatch(r"pitch_\d+", os.path.basename(r))
                   and os.path.exists(os.path.join(r, "fields.vtkhdf"))),
                  key=lambda r: int(os.path.basename(r)[6:]))
    pitch = [int(os.path.basename(r)[6:]) for r in runs]
    data = [read(os.path.join(r, "fields.vtkhdf")) for r in runs]
    hi = max(float(np.max(d[2]["speed"])) for d in data)
    fig = plt.figure(figsize=(12.0, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.0])
    xlim, ylim = (0.07, 0.19), (0.0, 0.04)
    for k, (p, (quads, _, d)) in enumerate(zip(pitch, data)):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        xs, ys, solid = grid(quads, (d["zone"] > 0).astype(float))
        panel(ax, quads, d["speed"], "viridis", "speed (m/s)", fig, xlim, ylim, 0.0, hi)
        outline(ax, xs, ys, solid, "w")
        ax.set_aspect("equal")
        ax.set_title(f"pitch {p} mm", fontsize=10)
        ax.set_ylabel("y (m)")
        if k // 2 == 1:
            ax.set_xlabel("x (m)")
    ax = fig.add_subplot(gs[2, :])
    cmap = plt.get_cmap("plasma")
    for k, r in enumerate(runs):
        with open(os.path.join(r, "line_base_top.csv")) as f:
            lines = [ln for ln in f if not ln.startswith("#")]
        t = list(csv.DictReader(lines))
        ax.plot([float(a["x"]) * 1e3 for a in t], [float(a["temperature"]) for a in t], lw=1.8,
                color=cmap(0.15 + 0.7 * k / max(len(runs) - 1, 1)), label=f"pitch {pitch[k]} mm")
    ax.set_xlabel("x along the top of the base (mm)")
    ax.set_ylabel("temperature (K)")
    ax.grid(alpha=0.3)
    ax.set_ylim(top=ax.get_ylim()[1] + 1.4)
    ax.legend(fontsize=9, ncol=4, loc="upper center")
    ax.set_title("Read by the assistant by page (result.page, series line_base_top, 100 cells per run): "
                 "the temperature along the base top", fontsize=10)
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


def main():
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == "summary":
        summary(a[1], a[2])
    elif len(a) in (4, 5) and a[0] == "figures":
        rows = session_rows(a[1])
        timeline(rows, a[4] if len(a) == 5 else a[2], os.path.join(a[3], "heat_sink_agent.png"))
        flow(a[2], os.path.join(a[3], "heat_sink_flow.png"))
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
