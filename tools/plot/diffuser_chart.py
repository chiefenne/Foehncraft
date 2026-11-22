"""Figure of the asymmetric diffuser (gallery 197): speed and turbulent
viscosity ratio from sample lines of a run of
examples/cases/flow/diffuser_sst.toml, and the wall pressure coefficient.

    python3 -m tools.plot.diffuser_chart OUT.png [WORKDIR]

Runs the case with vertical sample lines (`output.lines`) between the
walls at 160 stations from x = -5 H to x = 50 H and draws the speed and
the turbulent viscosity ratio with matplotlib (pip install matplotlib
numpy). Needs the `foehn` executable (FOEHN=target/release/foehn).
"""

import math
import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
FOEHN = os.environ.get("FOEHN", str(ROOT / "target/release/foehn"))
TAN = math.tan(math.radians(10.0))
L_DIFF = 3.7 / TAN


def ylow(x):
    return -min(max(x, 0.0), L_DIFF) * TAN


def stations():
    return [float(x) for x in np.concatenate([np.linspace(-5, 0, 20), np.linspace(0.2, 21, 60), np.linspace(21.5, 50, 80)])]


def read(path):
    rows = [x for x in path.read_text().splitlines() if x and not x.startswith("#")]
    head = rows[0].split(",")
    return [dict(zip(head, map(float, r.split(",")))) for r in rows[1:]]


def main(out: str, work: str = "/tmp/diffuser_chart") -> None:
    w = Path(work)
    w.mkdir(parents=True, exist_ok=True)
    text = (ROOT / "examples/cases/flow/diffuser_sst.toml").read_text()
    text = text.replace('file = "../../meshes/diffuser_coarse.su2"', f'file = "{ROOT / "examples/meshes/diffuser_coarse.su2"}"')
    text = text.replace('table = "diffuser_inlet_profile.csv"', f'table = "{ROOT / "examples/cases/flow/diffuser_inlet_profile.csv"}"')
    lines = ", ".join(
        f'{{ name = "s{i}", from = [{x}, {ylow(x) + 0.004}, 0.5], to = [{x}, 0.996, 0.5], points = 70 }}'
        for i, x in enumerate(stations()))
    text = text.replace("wall_csv = true\n", f"wall_csv = true\nlines = [ {lines} ]\n", 1)
    (w / "case.toml").write_text(text)
    subprocess.run([FOEHN, "flow", str(w / "case.toml"), "--out", str(w / "run")], check=True, stdout=subprocess.DEVNULL)
    xs, ys, sp, nt = [], [], [], []
    for i, _ in enumerate(stations()):
        for d in read(w / "run" / f"line_s{i}.csv"):
            xs.append(d["x"]); ys.append(d["y"])
            sp.append((d["u"] ** 2 + d["v"] ** 2) ** 0.5); nt.append(d["nu_t"] / 5e-5)
    xs, ys = np.array(xs), np.array(ys)
    tri = mtri.Triangulation(xs, ys)
    cx, cy = xs[tri.triangles].mean(axis=1), ys[tri.triangles].mean(axis=1)
    tri.set_mask(cy < np.array([ylow(x) for x in cx]) - 1e-9)
    fig, axs = plt.subplots(3, 1, figsize=(10, 7.6), dpi=200, gridspec_kw={"height_ratios": [1, 1, 0.8]})
    for ax, v, lab, cmap in [(axs[0], sp, "speed |u| / U_c", "turbo"), (axs[1], nt, "turbulent viscosity ratio ν_t / ν", "magma")]:
        t = ax.tricontourf(tri, v, levels=40, cmap=cmap)
        fig.colorbar(t, ax=ax, label=lab, pad=0.01)
        ax.set_aspect("equal"); ax.set_ylabel("y / H")
    for name, c in [("upper", "tab:blue"), ("lower", "tab:red")]:
        rows = sorted(read(w / "run" / f"wall_{name}.csv"), key=lambda r: r["x"])
        x = np.array([r["x"] for r in rows]); p = np.array([r["pressure"] for r in rows])
        p0 = np.interp(-10.0, x, p)
        axs[2].plot(x, 2 * (p - p0), color=c, label=f"{name} wall")
    axs[2].set_xlim(-5, 50); axs[2].set_xlabel("x / H"); axs[2].set_ylabel("C_p"); axs[2].grid(alpha=0.3); axs[2].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out)


if __name__ == "__main__":
    main(*sys.argv[1:])
