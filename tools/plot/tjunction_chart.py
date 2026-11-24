"""Figure of the T-junction mixer (gallery 196): speed and passive scalar on
the mid plane z = 0, from sample lines of a run of
examples/cases/flow/tjunction_mixing.toml on a mesh of the study.

    python3 -m tools.plot.tjunction_chart OUT.png MESH.msh [WORKDIR]

The script runs the case with 100 vertical sample lines in the main pipe
and 15 in the branch (`output.lines`), reads their values and draws them
with matplotlib (pip install matplotlib numpy). Needs the `foehn`
executable (FOEHN=target/release/foehn by default).
"""

import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from tools.plot.colormaps import RAINBOW  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
FOEHN = os.environ.get("FOEHN", str(ROOT / "target/release/foehn"))


def sample_lines():
    out = []
    for i, x in enumerate(np.linspace(0.01, 0.99, 100)):
        out.append((f"m{i}", [float(x), -0.0495, 0.0], [float(x), 0.0495, 0.0], 60))
    for i, x in enumerate(np.linspace(0.277, 0.323, 15)):
        out.append((f"b{i}", [float(x), 0.05, 0.0], [float(x), 0.295, 0.0], 60))
    return out


def main(out: str, mesh: str, work: str = "/tmp/tjunction_chart") -> None:
    w = Path(work)
    w.mkdir(parents=True, exist_ok=True)
    text = (ROOT / "examples/cases/flow/tjunction_mixing.toml").read_text()
    text = text.replace('file = "../../meshes/tjunction_gci_coarse.msh"', f'file = "{Path(mesh).resolve()}"')
    lines = ", ".join(f'{{ name = "{n}", from = {a}, to = {b}, points = {p} }}' for n, a, b, p in sample_lines())
    text = text[: text.index("lines = [")] + f"lines = [ {lines} ]\n"
    (w / "case.toml").write_text(text)
    subprocess.run([FOEHN, "flow", str(w / "case.toml"), "--out", str(w / "run")], check=True,
                   stdout=subprocess.DEVNULL)
    xs, ys, sp, th = [], [], [], []
    for n, *_ in sample_lines():
        rows = [l for l in (w / "run" / f"line_{n}.csv").read_text().splitlines() if l and not l.startswith("#")]
        head = rows[0].split(",")
        for r in rows[1:]:
            d = dict(zip(head, map(float, r.split(","))))
            xs.append(d["x"]); ys.append(d["y"])
            sp.append((d["u"] ** 2 + d["v"] ** 2 + d["w"] ** 2) ** 0.5); th.append(d["theta"])
    xs, ys = np.array(xs), np.array(ys)
    tri = mtri.Triangulation(xs, ys)
    cx = xs[tri.triangles].mean(axis=1)
    cy = ys[tri.triangles].mean(axis=1)
    # triangles that span the empty corners between main pipe and branch
    tri.set_mask((cy > 0.0505) & (np.abs(cx - 0.3) > 0.0262)
                 | ((cy > 0.0505) & (np.ptp(xs[tri.triangles], axis=1) > 0.06))
                 | ((cy > 0.0) & (np.ptp(ys[tri.triangles], axis=1) > 0.12)))
    fig, axs = plt.subplots(2, 1, figsize=(9, 5.4), dpi=200, sharex=True)
    for ax, val, lab, cmap in [(axs[0], sp, "speed |u| (m/s)", RAINBOW), (axs[1], th, "passive scalar θ", RAINBOW)]:
        t = ax.tricontourf(tri, val, levels=40, cmap=cmap)
        fig.colorbar(t, ax=ax, label=lab, pad=0.02)
        ax.set_aspect("equal")
        ax.set_ylabel("y (m)")
    axs[1].set_xlabel("x (m)")
    fig.tight_layout()
    fig.savefig(out)


if __name__ == "__main__":
    main(*sys.argv[1:])
