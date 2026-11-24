"""The consistent QCR stress divergence (D236 step 2), gallery 194: the
error of the explicit QCR2000 stress force against the manufactured
solution of `turbulence::mms`, per cell on tetrahedra with 15 % perturbed
vertices, before (face-neighbour least-squares gradient, linear face
interpolation) and after (quadratic least-squares gradient on the flagged
cells, skewness-corrected face stress), and the convergence of the L2
error on perturbed hexahedra and tetrahedra. The cells shown are one
quarter of the test window (x, y < 0.5).

From the per-cell dumps and the logs of two test runs, one on a build
before the change (`old`) and one after (`new`):

    FOEHN_MMS_QCR_DUMP=output/qcr_new cargo test --release -p foehncraft-flow \\
        --lib qcr_stress_on_perturbed -- --nocapture 2> output/qcr_new.log
    pip install matplotlib numpy
    python3 -m tools.plot.qcr_consistency docs/figures/gallery/qcr_consistency.png \\
        output/qcr_old output/qcr_new

Each folder holds `k-omega-sst-qcr2000_<mesh>_0.15_<level>.csv` and the
log `<folder>.log` with the lines `perturbed <mesh>: [(h, error), ...]`.
"""
import ast
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.tri import Triangulation
from tools.plot.colormaps import RAINBOW


def cells(folder, mesh, level):
    d = np.genfromtxt(folder / f"k-omega-sst-qcr2000_{mesh}_0.15_{level}.csv",
                      delimiter=",", names=True)
    big = max(np.abs(d["fx"]).max(), np.abs(d["fy"]).max())
    e = np.hypot(d["ex"], d["ey"]) / big
    keep = (d["x"] < 0.5) & (d["y"] < 0.5)
    return d["x"][keep], d["y"][keep], e[keep]


def errors(folder):
    text = Path(str(folder) + ".log").read_text()
    out = {}
    for mesh, pairs in re.findall(r"perturbed (\w+): (\[.*?\])", text):
        out[mesh] = np.array(ast.literal_eval(pairs))
    return out


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, old, new = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    fig = plt.figure(figsize=(10, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.8])
    norm = LogNorm(1e-3, 3e-1)
    for k, (folder, title) in enumerate([(old, "before: linear face stress"),
                                         (new, "after: quadratic fit, skew-corrected")]):
        ax = fig.add_subplot(gs[0, k])
        x, y, e = cells(folder, "tet", 2)
        tri = Triangulation(x, y)
        t = tri.triangles
        edge = np.max([np.hypot(x[t[:, i]] - x[t[:, j]], y[t[:, i]] - y[t[:, j]])
                       for i, j in [(0, 1), (1, 2), (2, 0)]], axis=0)
        tri.set_mask(edge > 4.0 * np.median(edge))
        pc = ax.tripcolor(tri, np.maximum(e, 1e-4), shading="gouraud",
                          cmap=RAINBOW, norm=norm)
        ax.set_xlim(0.14, 0.43)
        ax.set_ylim(0.14, 0.43)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        fin = errors(folder)["tet"][-1, 1] * 100
        ax.set_title(f"{title}\nL2 error of the window {fin:.1f} %", fontsize=11)
    fig.colorbar(pc, ax=fig.axes[:2], shrink=0.9, pad=0.02,
                 label="|f_h − f| / max |f|")
    ax = fig.add_subplot(gs[1, :])
    styles = {"hex": ("o", "tab:blue"), "tet": ("^", "tab:orange")}
    for folder, ls, tag in [(old, "--", "before"), (new, "-", "after")]:
        for mesh, e in errors(folder).items():
            m, c = styles[mesh]
            ax.loglog(e[:, 0], e[:, 1] * 100, ls, marker=m, color=c,
                      label=f"perturbed {'hexahedra' if mesh == 'hex' else 'tetrahedra'}, {tag}")
    h = np.array([0.0025, 0.014])
    ax.loglog(h, 0.3 * (h / h[0]), ":", color="0.4", label="order 1")
    ax.loglog(h, 0.12 * (h / h[0]) ** 2, "-.", color="0.4", label="order 2")
    ax.set_xlabel("cell size h")
    ax.set_ylabel("L2 error of the QCR force (% of max |f|)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9, ncol=2, loc="lower right")
    fig.savefig(out, dpi=100)


if __name__ == "__main__":
    main()
