"""Grid convergence of the laminar T-junction mixer (D45, stage 8): three
tetrahedral Gmsh meshes converted into polyhedra, the mean pressures of the
two inlets and the scalar difference across the pipe 0.6 m behind the
junction, observed order, extrapolated value and GCI_fine (Celik et al.
2008) for non-constant refinement ratios r = (N_coarse / N_fine)^(1/3).
No reference value of another solver is used (D201); a comparison with
published measurements is pending.

    python3 examples/studies/tjunction_gci.py [OUT] [--report-only]

Base case: examples/cases/flow/tjunction_mixing.toml; the meshes come from
examples/meshes/tjunction_gci.geo with h = 0.0126, 0.0079 and 0.0055 m
(needs `gmsh`, `pip install gmsh`). Writes OUT/report.md and results.json.
Needs the `foehn` executable (FOEHN=target/release/foehn by default).
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "examples/cases/flow/tjunction_mixing.toml"
GEO = ROOT / "examples/meshes/tjunction_gci.geo"
FOEHN = os.environ.get("FOEHN", str(ROOT / "target/release/foehn"))
SIZES = [0.0126, 0.0079, 0.0055]      # coarse to fine


def make_mesh(base: Path, h: float) -> Path:
    msh = base / f"tj_{h}.msh"
    if not msh.exists():
        subprocess.run(["gmsh", "-3", str(GEO), "-setnumber", "h", str(h), "-format", "msh41",
                        "-o", str(msh), "-v", "1"], check=True, stdout=subprocess.DEVNULL)
    return msh


def run(base: Path, h: float) -> None:
    d = base / f"run_{h}"
    if (d / "flow.json").exists():
        return
    d.mkdir(parents=True, exist_ok=True)
    msh = make_mesh(base, h)
    cmd = [FOEHN, "flow", str(CASE), "--out", str(d), "--set", f'mesh.file="{msh}"']
    with open(d.parent / (d.name + ".log"), "w") as log:
        r = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    print(f"h {h}: exit {r.returncode}", flush=True)


def line(d: Path):
    rows = [l for l in (d / "line_cross.csv").read_text().splitlines() if l and not l.startswith("#")]
    head = rows[0].split(",")
    return [dict(zip(head, map(float, r.split(",")))) for r in rows[1:]]


def at(pts, y: float) -> float:
    for a, b in zip(pts, pts[1:]):
        if a["y"] <= y <= b["y"]:
            t = (y - a["y"]) / (b["y"] - a["y"])
            return a["theta"] + t * (b["theta"] - a["theta"])
    raise ValueError(y)


def values(d: Path) -> dict | None:
    if not (d / "flow.json").exists():
        return None
    s = json.loads((d / "flow.json").read_text())["solves"][-1]
    b = {x["name"]: x for x in s["boundaries"]}
    pts = line(d)
    return {"cells": s["cells"], "iterations": s["iterations"], "converged": s.get("converged"),
            "p_main": b["inlet_main"]["mean_pressure"], "p_branch": b["inlet_branch"]["mean_pressure"],
            "dtheta": at(pts, 0.035) - at(pts, -0.035)}


def celik3(f: list[float], n: list[int]) -> dict:
    """Celik et al. (2008), coarse to fine f1<-fine: f = [fine, medium, coarse]."""
    f1, f2, f3 = f
    r21 = (n[0] / n[1]) ** (1 / 3)
    r32 = (n[1] / n[2]) ** (1 / 3)
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0:
        return {"verdict": "unchanged", "p": None, "ext": f1, "gci": 0.0}
    s = math.copysign(1.0, e32 / e21)
    if s < 0:
        return {"verdict": "oscillatory", "p": None, "ext": None,
                "gci": 0.5 * (max(f) - min(f)) / abs(f1)}
    p = 1.0
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        pn = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(pn - p) < 1e-10:
            break
        p = pn
    ext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
    gci = 1.25 * abs((f1 - f2) / f1) / (r21 ** p - 1)
    return {"verdict": "monotonic", "p": p, "ext": ext, "gci": gci, "r21": r21, "r32": r32}


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    base = Path(args[0] if args else "output/tjunction_gci")
    base.mkdir(parents=True, exist_ok=True)
    if "--report-only" not in argv:
        for h in SIZES:
            run(base, h)
    vs = [values(base / f"run_{h}") for h in SIZES]
    if any(v is None for v in vs):
        print("runs missing")
        return 1
    n = [v["cells"] for v in reversed(vs)]    # fine, medium, coarse
    lines = ["# Grid convergence of the T-junction mixer, Re 200 / 100", "",
             "Laminar, steady, polyhedra converted from Gmsh tetrahedra; GCI_fine after Celik et al. (2008) "
             "with r = (N_coarse/N_fine)^(1/3). Pending: comparison with published measurements.", "",
             "| cells | p inlet_main Pa | p inlet_branch Pa | theta(+0.035) − theta(−0.035) | iterations |",
             "|---|---|---|---|---|"]
    for v in vs:
        lines.append(f"| {v['cells']} | {v['p_main']:.5f} | {v['p_branch']:.5f} | {v['dtheta']:.5f} | "
                     f"{v['iterations']}{'' if v['converged'] else '*'} |")
    lines.append("")
    res = {"values": vs}
    for key, label in [("p_main", "p inlet_main"), ("p_branch", "p inlet_branch"), ("dtheta", "scalar difference")]:
        g = celik3([vs[2][key], vs[1][key], vs[0][key]], n)
        res[key] = g
        p = "–" if g["p"] is None else f"{g['p']:.2f}"
        ext = "–" if g["ext"] is None else f"{g['ext']:.5f}"
        lines.append(f"- {label}: p = {p} ({g['verdict']}), GCI_fine {100 * g['gci']:.2f} %, extrapolated {ext}")
    text = "\n".join(lines) + "\n"
    (base / "report.md").write_text(text)
    (base / "results.json").write_text(json.dumps(res, indent=1))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
