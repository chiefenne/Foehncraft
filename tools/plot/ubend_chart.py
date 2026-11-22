"""Figure of the U-duct (gallery 204): speed and turbulent viscosity ratio
from sample lines of a run of examples/cases/flow/ubend_sst.toml (or
ubend_sst_sharp.toml with --sharp), and the wall pressure coefficient
along both walls.

    python3 -m tools.plot.ubend_chart OUT.png [WORKDIR] [--sharp]

Runs the case on the committed coarse mesh with sample lines across the
duct (`output.lines`): straight lines in the legs, radial lines in the
rounded bend (sharp: vertical lines through the turn gap), and draws the
fields with matplotlib (pip install matplotlib numpy). Needs the `foehn`
executable (FOEHN=target/release/foehn).
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
NU = 1e-5
E = 0.002   # sample lines stop this far from the walls


def lines(sharp: bool):
    """(from, to) of each sample line."""
    out = []
    for x in np.concatenate([np.linspace(-6.0, -1.0, 12), np.linspace(-0.9, -0.02, 12)]):
        out.append(((x, -1.5 + E), (x, -0.5 - E)))
    if sharp:
        for x in np.linspace(0.005, 0.995, 40):
            out.append(((x, -1.5 + E), (x, 1.5 - E)))
    else:
        for p in np.linspace(-0.5 * math.pi, 0.5 * math.pi, 61):
            c, s = math.cos(p), math.sin(p)
            out.append(((0.5 * c + E * c, 0.5 * s + E * s), (1.5 * c - E * c, 1.5 * s - E * s)))
    for x in np.concatenate([np.linspace(-0.02, -2.0, 30), np.linspace(-2.2, -12.0, 30)]):
        out.append(((x, 0.5 + E), (x, 1.5 - E)))
    return out


def read(path):
    rows = [x for x in path.read_text().splitlines() if x and not x.startswith("#")]
    head = rows[0].split(",")
    return [dict(zip(head, map(float, r.split(",")))) for r in rows[1:]]


def outside(x, y, sharp: bool):
    """True where a point lies in the divider or beyond the outer wall."""
    if x < 0.0:
        return abs(y) < 0.5
    if sharp:
        return False
    r = math.hypot(x, y)
    return r < 0.5 or r > 1.5


def run(work: Path, sharp: bool) -> Path:
    name = "ubend_sst_sharp" if sharp else "ubend_sst"
    mesh = "ubend_sharp_coarse.su2" if sharp else "ubend_coarse.su2"
    text = (ROOT / f"examples/cases/flow/{name}.toml").read_text()
    text = text.replace(f'file = "../../meshes/{mesh}"', f'file = "{ROOT / "examples/meshes" / mesh}"')
    text = text.replace('table = "ubend_inlet_profile.csv"', f'table = "{ROOT / "examples/cases/flow/ubend_inlet_profile.csv"}"')
    ls = ", ".join(f'{{ name = "s{i}", from = [{a[0]}, {a[1]}, 0.5], to = [{b[0]}, {b[1]}, 0.5], points = 90 }}'
                   for i, (a, b) in enumerate(lines(sharp)))
    text = text.replace("wall_csv = true\n", f"wall_csv = true\nlines = [ {ls} ]\n", 1)
    work.mkdir(parents=True, exist_ok=True)
    (work / "case.toml").write_text(text)
    if not (work / "run" / "flow.json").exists():
        subprocess.run([FOEHN, "flow", str(work / "case.toml"), "--out", str(work / "run")], check=True,
                       stdout=subprocess.DEVNULL)
    return work / "run"


def fields(d: Path, sharp: bool):
    xs, ys, sp, nt = [], [], [], []
    for i, _ in enumerate(lines(sharp)):
        for r in read(d / f"line_s{i}.csv"):
            xs.append(r["x"]); ys.append(r["y"])
            sp.append(math.hypot(r["u"], r["v"])); nt.append(r["nu_t"] / NU)
    xs, ys = np.array(xs), np.array(ys)
    tri = mtri.Triangulation(xs, ys)
    cx, cy = xs[tri.triangles].mean(axis=1), ys[tri.triangles].mean(axis=1)
    big = np.max(np.hypot(np.diff(xs[tri.triangles], axis=1, append=xs[tri.triangles][:, :1]),
                          np.diff(ys[tri.triangles], axis=1, append=ys[tri.triangles][:, :1])), axis=1)
    tri.set_mask(np.array([outside(a, b, sharp) for a, b in zip(cx, cy)]) | (big > 0.6))
    return tri, np.array(sp), np.array(nt)


def main(argv) -> None:
    sharp = "--sharp" in argv
    args = [a for a in argv if not a.startswith("--")]
    out = args[0]
    work = Path(args[1] if len(args) > 1 else f"/tmp/ubend_chart{'_sharp' if sharp else ''}")
    d = run(work, sharp)
    sys.path.insert(0, str(ROOT / "examples/studies"))
    import ubend_gci as g  # noqa: E402  (stations along the walls)
    g.SHARP = sharp
    tri, sp, nt = fields(d, sharp)
    fig = plt.figure(figsize=(10, 8.2), dpi=200)
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.9])
    for k, (v, lab, cmap) in enumerate([(sp, "speed |u| / U_b", "turbo"), (nt, "turbulent viscosity ratio ν_t / ν", "magma")]):
        ax = fig.add_subplot(gs[k])
        t = ax.tricontourf(tri, v, levels=48, cmap=cmap)
        fig.colorbar(t, ax=ax, label=lab, pad=0.01)
        ax.set_aspect("equal"); ax.set_xlim(-8, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_ylabel("y / H")
        ax.set_facecolor("#e8e8e8")
    ax = fig.add_subplot(gs[2])
    outer, inner = g.along(d, "outer"), g.along(d, "inner")
    p0 = g.p_at(outer, g.station(-2.0, -1.5))
    for w, name, c in [(inner, "inner wall", "tab:red"), (outer, "outer wall", "tab:blue")]:
        s = np.array([r[0] for r in w]); p = np.array([r[2] for r in w])
        ax.plot(s, 2 * (p - p0), color=c, label=name)
    if not sharp:
        ax.axvspan(0.0, math.pi, color="0.92", zorder=0)
        ax.text(0.5 * math.pi, 0.95, "bend (θ in radians)", ha="center", va="top", transform=ax.get_xaxis_transform())
        ax.set_xlabel("position along the wall: x in the inlet leg, θ in the bend, π + distance behind the bend (H)")
    else:
        ax.set_xlabel("position along the wall (H): x in the inlet leg, then along the divider end or the end wall")
    ax.set_ylabel("C_p"); ax.set_xlim(-4, 16); ax.grid(alpha=0.3); ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out)


if __name__ == "__main__":
    main(sys.argv[1:])
