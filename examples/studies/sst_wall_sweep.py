"""The first-cell y⁺ sweep of the SST wall treatments on the flat plate
(stage 6 specification 9.3; verification manual, k-ω SST part 6b and 6g).

    python3 examples/studies/sst_wall_sweep.py OUT [--jobs N] [--treatment automatic]
        [--yplus 1.7,4.6,9.3,27,92] [--level 1] [--iterations 10000] [--set KEY=VALUE ...]

Runs `tmr_plate_sst.toml` on mesh level `--level` (64 × 48 cells for 1)
with the wall grading set so that the first cell at x = 0.97 has the
given y⁺ (from u_τ of the reference C_f), and reports C_f(0.97) against the
reference in OUT/report.md (by default our own `low-re` result on the finest
plate mesh, level 3, first-cell y⁺ 0.09: 2.7005e-3, the plate GCI study of
D201's rework; no number of another solver), with the y⁺ the run reached. Finished
runs (a flow.json in their folder) are kept; delete a folder to rerun it.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FOEHN = str(ROOT / "target/release/foehn")
CASE = ROOT / "examples/cases/flow/tmr_plate_sst.toml"
CF_REF = 2.7005e-3     # own low-re result, level 3 (plate_channel_gci.py)
NU = 2.0e-7
HEIGHT = 1.0


def ratio_for(y_plus: float, cells: int) -> float:
    """The grading ratio (last/first cell height) of `cells` cells over the
    height 1 m whose first cell centre lies at y⁺ (u_τ from the reference C_f)."""
    u_tau = math.sqrt(CF_REF / 2.0)
    h1 = 2.0 * y_plus * NU / u_tau

    def first(q):
        return HEIGHT * (q - 1.0) / (q**cells - 1.0)

    lo, hi = 1.0 + 1e-9, 2.0
    while first(hi) > h1:
        hi *= 1.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if first(mid) > h1:
            lo = mid
        else:
            hi = mid
    q = 0.5 * (lo + hi)
    return q ** (cells - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--treatment", default="automatic")
    ap.add_argument("--yplus", default="1.7,4.6,9.3,27,92")
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--iterations", type=int, default=10000)
    ap.add_argument("--set", action="append", default=[])
    a = ap.parse_args()
    base = Path(a.out)
    base.mkdir(parents=True, exist_ok=True)
    cells = 24 * 2**a.level
    targets = [float(v) for v in a.yplus.split(",")]

    def run(yp):
        d = base / f"yplus_{yp:g}".replace(".", "p")
        if not (d / "flow.json").exists():
            r = ratio_for(yp, cells)
            cmd = [FOEHN, "flow", str(CASE), "--out", str(d), "--threads", "1",
                   "--set", f"turbulence.wall_treatment=\"{a.treatment}\"",
                   "--set", f"mesh.generate.grading.y.ratio={r}",
                   "--set", f"verification.levels=[{a.level}]",
                   "--set", f"solver.max_iterations={a.iterations}",
                   "--set", "verification.require_convergence=false"]
            for s in a.set:
                cmd += ["--set", s]
            d.mkdir(parents=True, exist_ok=True)
            with open(d.with_suffix(".log"), "w") as log:
                subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        f = d / "flow.json"
        if not f.exists():
            return yp, None
        s = json.loads(f.read_text())["solves"][-1]
        p = next(w for w in s["wall_probes"] if w["name"] == "x097")
        return yp, {"cf": p["skin_friction"], "y_plus": p["y_plus"],
                    "iterations": s["iterations"], "converged": s.get("converged")}

    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(run, targets))
    lines = [f"# SST `{a.treatment}` on the plate, level {a.level}", "",
             "| target y⁺ | y⁺ reached | C_f ·10³ | vs reference | iterations, converged |",
             "|---|---|---|---|---|"]
    for yp, v in res:
        if v is None:
            lines.append(f"| {yp:g} | failed | | | |")
            continue
        lines.append(f"| {yp:g} | {v['y_plus']:.2f} | {1e3 * v['cf']:.4f} | "
                     f"{100 * (v['cf'] / CF_REF - 1):+.1f} % | {v['iterations']}, {v['converged']} |")
    (base / "report.md").write_text("\n".join(lines) + "\n")
    (base / "results.json").write_text(json.dumps(dict((str(k), v) for k, v in res), indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
