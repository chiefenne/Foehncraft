"""Grid convergence of the U-duct of Monson et al. (1990; D45, stage 8):
the separated region on the inner wall behind the turn and the wall
pressure, k-omega SST low-Re at Re_b = 1e5 on three nested structured
meshes (examples/meshes/ubend.py, levels 0 to 2), observed order,
extrapolated value and GCI_fine after Celik et al. (2008), ratio 2. No
reference value is used (D201); the comparison with the published
measurements is pending.

    python3 examples/studies/ubend_gci.py [OUT] [--levels 0,1,2] [--jobs 3] [--report-only] [--sharp]

--sharp runs ubend_sst_sharp.toml on the meshes with sharp edges
(ubend.py --sharp: a square divider end and a square end wall).

Positions along the inner wall: in the bend the angle theta (0 deg at the
bend entry, 180 deg at its exit), behind it the distance s = -x from the
bend exit (sharp: from the divider end) in H. The separated region is the longest stretch of the inner
wall where the shear stress along the flow direction is negative.
C_p = (p - p_ref) / (1/2 rho U_b^2) with p_ref the outer-wall pressure at
x = -2 H in the inlet leg; taken on both walls in the middle of the turn
(theta 90 deg; sharp: y = 0 on the divider end and the end wall) and on the
outer wall at s = 10 H (the pressure loss of the turn).
The levels run at the same time (--jobs). Needs the `foehn` executable
(FOEHN=target/release/foehn by default).
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plate_channel_gci import celik, gci_txt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SHARP = "--sharp" in sys.argv
CASE = ROOT / ("examples/cases/flow/ubend_sst_sharp.toml" if SHARP else "examples/cases/flow/ubend_sst.toml")
FOEHN = os.environ.get("FOEHN", str(ROOT / "target/release/foehn"))
MESHER = ROOT / "examples/meshes/ubend.py"


def run(base: Path, level: int) -> None:
    d = base / f"run_{level}"
    if (d / "flow.json").exists():
        return
    d.mkdir(parents=True, exist_ok=True)
    su2 = base / f"ubend_{level}.su2"
    subprocess.run([sys.executable, str(MESHER), str(su2), "--level", str(level)] + (["--sharp"] if SHARP else []),
                   check=True)
    cmd = [FOEHN, "flow", str(CASE), "--out", str(d), "--threads", "1",
           "--set", f'mesh.file="{su2}"']
    with open(base / f"run_{level}.log", "w") as log:
        r = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    print(f"level {level}: exit {r.returncode}", flush=True)


def wall(d: Path, name: str):
    rows = [x for x in (d / f"wall_{name}.csv").read_text().splitlines() if x and not x.startswith("#")]
    head = rows[0].split(",")
    return [dict(zip(head, map(float, r.split(",")))) for r in rows[1:]]


def station(x: float, y: float) -> float:
    """Monotone position along a wall in flow direction (H, radians in the
    rounded bend). Rounded: the inlet leg x (< 0), the bend theta (0 to pi),
    the outlet leg pi - x. Sharp: the inlet leg x, then the divider end or
    the end wall by y (0 to 1 or 3), then the outlet leg."""
    if not SHARP:
        if x > 0.0:
            return math.atan2(y, x) + 0.5 * math.pi
        return x if y < 0.0 else math.pi - x
    if abs(y) < 0.5 - 1e-12 and abs(x) < 1e-12:      # divider end
        return y + 0.5
    if x > 1.0 - 1e-12:                              # end wall
        return 1.0 + (y + 1.5)
    if y < 0.0:                                      # inlet leg and gap, lower side
        return x
    return (1.0 if y < 1.0 else 5.0) - x             # outlet leg (and gap, upper side)


def tangent(x: float, y: float):
    if not SHARP and x > 0.0:
        phi = math.atan2(y, x)
        return (-math.sin(phi), math.cos(phi))
    if SHARP and ((abs(y) < 0.5 - 1e-12 and abs(x) < 1e-12) or x > 1.0 - 1e-12):
        return (0.0, 1.0)
    return (1.0, 0.0) if y < 0.0 else (-1.0, 0.0)


def along(d: Path, name: str):
    out = []
    for r in wall(d, name):
        t = tangent(r["x"], r["y"])
        out.append((station(r["x"], r["y"]), r["tau_x"] * t[0] + r["tau_y"] * t[1], r["pressure"], r["x"], r["y"]))
    return sorted(out)


def bubble(w) -> tuple[float, float]:
    """Start and end (stations) of the longest stretch with negative shear."""
    cuts, sign = [(w[0][0], w[0][1] < 0.0)], w[0][1] < 0.0
    for a, b in zip(w, w[1:]):
        if (b[1] < 0.0) != sign:
            t = a[1] / (a[1] - b[1])
            cuts.append((a[0] + t * (b[0] - a[0]), not sign))
            sign = not sign
    cuts.append((w[-1][0], not sign))
    runs = [(x0, x1) for (x0, neg), (x1, _) in zip(cuts, cuts[1:]) if neg]
    return max(runs, key=lambda r: r[1] - r[0])


def p_at(w, s: float) -> float:
    for a, b in zip(w, w[1:]):
        if a[0] <= s <= b[0]:
            t = (s - a[0]) / (b[0] - a[0]) if b[0] > a[0] else 0.0
            return a[2] + t * (b[2] - a[2])
    raise ValueError(s)


def values(d: Path) -> dict | None:
    if not (d / "flow.json").exists():
        return None
    s = json.loads((d / "flow.json").read_text())["solves"][-1]
    inner, outer = along(d, "inner"), along(d, "outer")
    s0, s1 = bubble(inner)
    exit_ = station(0.0, 0.5)                        # the inner wall leaves the turn
    p0 = p_at(outer, station(-2.0, -1.5))
    mid_in, mid_out = ((0.0, 0.0), (1.0, 0.0)) if SHARP else ((0.5, 0.0), (1.5, 0.0))
    v = {"cells": s["cells"], "iterations": s["iterations"], "converged": s.get("converged"),
         "s_reat": s1 - exit_,
         "cp_in_90": 2 * (p_at(inner, station(*mid_in)) - p0),
         "cp_out_90": 2 * (p_at(outer, station(*mid_out)) - p0),
         "cp_out_10": 2 * (p_at(outer, station(-10.0, 1.5)) - p0)}
    # rounded: the separation angle in the bend; sharp: the start behind the divider end (s/H)
    v["sep"] = s0 - exit_ if SHARP or s0 >= exit_ else math.degrees(s0)
    return v


KEYS = [("sep", "separation (sharp: s/H behind the divider end; rounded: theta in deg)"), ("s_reat", "reattachment s/H behind the bend"),
        ("cp_in_90", "C_p inner wall, 90 deg"), ("cp_out_90", "C_p outer wall, 90 deg"),
        ("cp_out_10", "C_p outer wall, 10 H behind the bend")]


def main(argv: list[str]) -> int:
    levels, jobs, skip = [0, 1, 2], 3, set()
    for flag in ("--levels", "--jobs"):
        if flag in argv:
            i = argv.index(flag)
            skip |= {i, i + 1}
            if flag == "--levels":
                levels = [int(x) for x in argv[i + 1].split(",")]
            else:
                jobs = int(argv[i + 1])
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and i not in skip]
    base = Path(args[0] if args else ("output/ubend_gci_sharp" if SHARP else "output/ubend_gci"))
    base.mkdir(parents=True, exist_ok=True)
    if "--report-only" not in argv:
        with ThreadPoolExecutor(jobs) as ex:
            list(ex.map(lambda lv: run(base, lv), levels))
    vs = [values(base / f"run_{lv}") for lv in levels]
    if any(v is None for v in vs):
        print("runs missing")
        return 1
    lines = [f"# Grid convergence of the U-duct{' with sharp edges' if SHARP else ''}, Rc/H = 1, Re_b 1e5", "",
             "k-omega SST low-Re; GCI_fine after Celik et al. (2008), ratio 2. Pending: comparison with published measurements.", "",
             "| cells | separation | s_reat/H | C_p in 90 | C_p out 90 | C_p out 10 H | iterations |",
             "|---|---|---|---|---|---|---|"]
    for v in vs:
        lines.append(f"| {v['cells']} | {v['sep']:.3f} | {v['s_reat']:.3f} | {v['cp_in_90']:.4f} | "
                     f"{v['cp_out_90']:.4f} | {v['cp_out_10']:.4f} | {v['iterations']}{'' if v['converged'] else '*'} |")
    lines.append("")
    res = {"values": vs}
    if len(vs) >= 3:
        for key, label in KEYS:
            g = celik(vs[-1][key], vs[-2][key], vs[-3][key])
            res[key] = g
            p = "–" if g["p"] is None else f"{g['p']:.2f}"
            ext = "–" if g["ext"] is None else f"{g['ext']:.4f}"
            lines.append(f"- {label}: p = {p} ({g['verdict']}), GCI_fine {gci_txt(g.get('gci'))}, extrapolated {ext}")
    lines += ["", "\\* the run stopped at its iteration limit before the tolerance."]
    text = "\n".join(lines) + "\n"
    (base / "report.md").write_text(text)
    (base / "results.json").write_text(json.dumps(res, indent=1))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
