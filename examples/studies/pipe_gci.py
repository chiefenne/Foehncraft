"""Grid convergence of the turbulent pipe at Re_b = 5300 (D45, stage 8):
friction velocity u_tau/U_b of k-omega SST 2003 (low-Re) on three
nested meshes, observed order p, extrapolated value and GCI_fine
(Celik et al. 2008), and the finest value against the DNS of
El Khoury et al. (2013), Table 2: Re_tau 181, u_tau/U_b = 0.0683
(C_f = 2 (u_tau/U_b)^2). A published DNS scalar only (D201).

    python3 examples/studies/pipe_gci.py [OUT] [--levels 0,1,2] [--report-only]

The base case is examples/cases/flow/pipe_turbulent_sst.toml (48 cells
across the radius, graded). The solution does not depend on x, so the
mesh is refined across the radius only; its nodes are those of level L of
the base mesh (same nested grading as the channel study). The runs
converge to the scaled residual 1e-8 (the monitored u_tau moves 3e-5
relative between 1e-8 and 4 000 iterations).
Writes OUT/report.md (also printed) and OUT/results.json. Needs the `foehn`
executable (FOEHN=target/release/foehn by default).
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plate_channel_gci import celik, gci_txt, nested_ratio, pct  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "examples/cases/flow/pipe_turbulent_sst.toml"
FOEHN = os.environ.get("FOEHN", str(ROOT / "target/release/foehn"))
R = 0.5
UTAU_DNS = 0.0683   # u_tau/U_b, El Khoury et al. (2013), Table 2, Re_b 5300
RE_TAU_DNS = 181


def mesh_sets(level: int) -> list[str]:
    gen = tomllib.loads(CASE.read_text())["mesh"]["generate"]
    nx, ny, nz = gen["cells"]
    gy = gen["grading"]["y"]
    ratio = nested_ratio(gy["ratio"], ny, False, level)
    return [f"mesh.generate.cells=[{nx},{ny << level},{nz}]",
            f"mesh.generate.grading.y.ratio={ratio!r}"]


def run(base: Path, level: int) -> None:
    d = base / f"pipe_sst_{level}"
    if (d / "flow.json").exists():
        return
    d.mkdir(parents=True, exist_ok=True)
    cmd = [FOEHN, "flow", str(CASE), "--out", str(d), "--threads", "1"]
    # 1e-8 on every mesh: u_tau of the finest mesh at 1e-8 and after 4 000
    # iterations (1e-9 not reached, residuals still falling steadily) differ
    # by 3e-5 relative, far below the GCI
    cmd += ["--set", "solver.tolerance=1e-8"]
    for s in mesh_sets(level):
        cmd += ["--set", s]
    with open(d.with_suffix(".log"), "w") as log:
        r = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    print(f"level {level}: exit {r.returncode}", flush=True)


def values(d: Path) -> dict | None:
    if not (d / "flow.json").exists():
        return None
    s = json.loads((d / "flow.json").read_text())["solves"][-1]
    f = s["periodic"][0]["driving_force"]       # N/m3, rho = 1
    ut = math.sqrt(f * R / 2.0)                 # tau_w = F R / 2
    return {"cells": s["cells"], "iterations": s["iterations"],
            "converged": s.get("converged"), "utau": ut, "cf": 2 * ut * ut}


def main(argv: list[str]) -> int:
    levels = [0, 1, 2]
    if "--levels" in argv:
        levels = [int(x) for x in argv[argv.index("--levels") + 1].split(",")]
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or argv[i - 1] != "--levels")]
    base = Path(args[0] if args else "output/pipe_gci")
    if "--report-only" not in argv:
        for lv in levels:
            run(base, lv)
    vs = [values(base / f"pipe_sst_{lv}") for lv in levels]
    if any(v is None for v in vs):
        print("runs missing")
        return 1
    g = celik(vs[-1]["utau"], vs[-2]["utau"], vs[-3]["utau"])
    gc = celik(vs[-1]["cf"], vs[-2]["cf"], vs[-3]["cf"])
    lines = ["# Grid convergence of the turbulent pipe, Re_b = 5300", "",
             f"k-omega SST 2003 low-Re; levels {', '.join(map(str, levels))}; GCI_fine after Celik et al. (2008), ratio 2 across the radius.",
             f"Reference: DNS of El Khoury et al. (2013) Table 2, Re_tau {RE_TAU_DNS}, u_tau/U_b = {UTAU_DNS} (C_f = {2 * UTAU_DNS**2:.5f}).", "",
             "| cells | u_tau/U_b | C_f·10³ | Re_tau | iterations, converged |", "|---|---|---|---|---|"]
    for v in vs:
        lines.append(f"| {v['cells']} | {v['utau']:.5f} | {1e3 * v['cf']:.4f} | {v['utau'] * R * 5300:.1f} | {v['iterations']}{'' if v['converged'] else '*'} |")
    lines += ["", f"u_tau/U_b: p = {g['p'] if g['p'] is None else round(g['p'], 2)} ({g['verdict']}), GCI_fine {gci_txt(g.get('gci'))}, "
              f"extrapolated {g['ext'] if g['ext'] is None else round(g['ext'], 5)}; finest against DNS {pct(vs[-1]['utau'] / UTAU_DNS - 1)}.",
              f"C_f: p = {gc['p'] if gc['p'] is None else round(gc['p'], 2)}, GCI_fine {gci_txt(gc.get('gci'))}; finest against DNS {pct(vs[-1]['cf'] / (2 * UTAU_DNS**2) - 1)}.", ""]
    text = "\n".join(lines)
    (base / "report.md").write_text(text)
    (base / "results.json").write_text(json.dumps({"values": vs, "gci_utau": g, "gci_cf": gc}, indent=1))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
