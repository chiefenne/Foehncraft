"""Grid convergence of the asymmetric plane diffuser (D45, stage 8):
separation and reattachment on the lower wall and the pressure recovery
C_p along the upper wall, k-omega SST low-Re on three nested structured
meshes (examples/meshes/diffuser.py, levels 0 to 2: 6 960, 27 840 and
111 360 cells), observed order, extrapolated value and GCI_fine after
Celik et al. (2008), ratio 2. No reference value is used (D201); the
comparison with the published measurements of Buice & Eaton is pending.

    python3 examples/studies/diffuser_gci.py [OUT] [--levels 0,1,2] [--jobs 3] [--report-only] [--sharp]

--sharp runs diffuser_sst_sharp.toml on the meshes with sharp corners
(diffuser.py --sharp: 8 400, 33 600 and 134 400 cells).

The separated region is the longest stretch of the lower wall with
tau_x < 0. C_p(x) = (p(x) - p_ref) / (1/2 rho U_c^2) on the upper wall with p_ref the
wall pressure at x = -10 H, taken at x = 21 H (end of the diffuser) and
x = 45 H (in the outlet duct). The levels run at the same time (--jobs).
Needs the `foehn` executable (FOEHN=target/release/foehn by default).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plate_channel_gci import celik, gci_txt, pct  # noqa: E402, F401

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "examples/cases/flow/diffuser_sst.toml"
SHARP = "--sharp" in sys.argv
if SHARP:
    CASE = ROOT / "examples/cases/flow/diffuser_sst_sharp.toml"
FOEHN = os.environ.get("FOEHN", str(ROOT / "target/release/foehn"))
MESHER = ROOT / "examples/meshes/diffuser.py"


def run(base: Path, level: int) -> None:
    d = base / f"run_{level}"
    if (d / "flow.json").exists():
        return
    d.mkdir(parents=True, exist_ok=True)
    su2 = base / f"diffuser_{level}.su2"
    subprocess.run([sys.executable, str(MESHER), str(su2), "--level", str(level)] + (["--sharp"] if SHARP else []),
                   check=True)
    cmd = [FOEHN, "flow", str(CASE), "--out", str(d), "--threads", "1",
           "--set", f'mesh.file="{su2}"']
    with open(base / f"run_{level}.log", "w") as log:
        r = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    print(f"level {level}: exit {r.returncode}", flush=True)


def wall(d: Path, name: str = "upper"):
    rows = [x for x in (d / f"wall_{name}.csv").read_text().splitlines() if x and not x.startswith("#")]
    head = rows[0].split(",")
    return sorted((dict(zip(head, map(float, r.split(",")))) for r in rows[1:]), key=lambda r: r["x"])


def bubble(d: Path) -> tuple[float, float]:
    """Start and end of the main separated region of the lower wall: the
    longest stretch with tau_x < 0 (a sharp corner leaves the shear near zero
    behind the corner, where small bubbles come and go with the mesh)."""
    w = wall(d, "lower")
    cuts, sign = [], w[0]["tau_x"] < 0.0
    for a, b in zip(w, w[1:]):
        if (b["tau_x"] < 0.0) != sign:
            t = a["tau_x"] / (a["tau_x"] - b["tau_x"])
            cuts.append((a["x"] + t * (b["x"] - a["x"]), not sign))
            sign = not sign
    runs = [(x0, x1) for (x0, neg), (x1, _) in zip(cuts, cuts[1:]) if neg]
    return max(runs, key=lambda r: r[1] - r[0])


def at(pts, x: float) -> float:
    for a, b in zip(pts, pts[1:]):
        if a["x"] <= x <= b["x"]:
            t = (x - a["x"]) / (b["x"] - a["x"])
            return a["pressure"] + t * (b["pressure"] - a["pressure"])
    raise ValueError(x)


def values(d: Path) -> dict | None:
    if not (d / "flow.json").exists():
        return None
    s = json.loads((d / "flow.json").read_text())["solves"][-1]
    x_sep, x_reat = bubble(d)
    w = wall(d)
    p0 = at(w, -10.0)
    return {"cells": s["cells"], "iterations": s["iterations"], "converged": s.get("converged"),
            "x_sep": x_sep, "x_reat": x_reat,
            "cp_21": 2 * (at(w, 21.0) - p0), "cp_45": 2 * (at(w, 45.0) - p0)}


def main(argv: list[str]) -> int:
    levels = [0, 1, 2]
    jobs = 3
    skip = set()
    for flag in ("--levels", "--jobs"):
        if flag in argv:
            i = argv.index(flag)
            skip |= {i, i + 1}
            if flag == "--levels":
                levels = [int(x) for x in argv[i + 1].split(",")]
            else:
                jobs = int(argv[i + 1])
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and i not in skip]
    base = Path(args[0] if args else "output/diffuser_gci")
    base.mkdir(parents=True, exist_ok=True)
    if "--report-only" not in argv:
        with ThreadPoolExecutor(jobs) as ex:
            list(ex.map(lambda lv: run(base, lv), levels))
    vs = [values(base / f"run_{lv}") for lv in levels]
    if any(v is None for v in vs):
        print("runs missing")
        return 1
    lines = [f"# Grid convergence of the asymmetric diffuser{' with sharp corners' if SHARP else ''}, Re_H 20 000", "",
             "k-omega SST low-Re; GCI_fine after Celik et al. (2008), ratio 2. Pending: comparison with published measurements.", "",
             "| cells | x_sep/H | x_reat/H | C_p(21 H) | C_p(45 H) | iterations |", "|---|---|---|---|---|---|"]
    for v in vs:
        lines.append(f"| {v['cells']} | {v['x_sep']:.3f} | {v['x_reat']:.2f} | {v['cp_21']:.4f} | {v['cp_45']:.4f} | "
                     f"{v['iterations']}{'' if v['converged'] else '*'} |")
    lines.append("")
    res = {"values": vs}
    if len(vs) >= 3:
        for key, label in [("x_sep", "separation x/H"), ("x_reat", "reattachment x/H"),
                           ("cp_21", "C_p at 21 H"), ("cp_45", "C_p at 45 H")]:
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
