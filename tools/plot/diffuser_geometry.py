"""Geometry, boundary conditions and mesh of the asymmetric diffuser
(gallery 199) as one figure.

    python3 -m tools.plot.diffuser_geometry OUT.png

Reads examples/meshes/diffuser_coarse.su2 (matplotlib, numpy).
"""

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def read_su2(path):
    lines = path.read_text().splitlines()
    ne = int(lines[1].split("=")[1])
    quads = [list(map(int, l.split()[1:5])) for l in lines[2:2 + ne]]
    npn = int(lines[2 + ne].split("=")[1])
    pts = np.array([list(map(float, l.split()[:2])) for l in lines[3 + ne:3 + ne + npn]])
    return pts, quads


def main(out):
    pts, quads = read_su2(ROOT / "examples/meshes/diffuser_coarse.su2")
    t = math.tan(math.radians(10))
    ld = 3.7 / t
    fig, axs = plt.subplots(3, 1, figsize=(11, 10), dpi=170, gridspec_kw={"height_ratios": [1.5, 1, 1]})
    ax = axs[0]
    poly = [(-20, 0), (0, 0), (ld, -3.7), (ld + 50, -3.7), (ld + 50, 1), (-20, 1)]
    ax.fill(*zip(*poly), color="#dbe9f6", ec="k", lw=1.5)
    ax.annotate("", (-20, -0.9), (0, -0.9), arrowprops=dict(arrowstyle="<->"))
    ax.text(-10, -1.7, "inlet duct 20 H", ha="center")
    ax.annotate("", (0, -4.4), (ld, -4.4), arrowprops=dict(arrowstyle="<->"))
    ax.text(ld / 2, -5.5, f"diffuser {ld:.2f} H (10°)", ha="center")
    ax.annotate("", (ld, -4.4), (ld + 50, -4.4), arrowprops=dict(arrowstyle="<->"))
    ax.text(ld + 25, -5.5, "outlet duct 50 H", ha="center")
    ax.annotate("", (-21.5, 0), (-21.5, 1), arrowprops=dict(arrowstyle="<->"))
    ax.text(-22, 0.5, "H = 1", ha="right", va="center")
    ax.annotate("", (ld + 52, -3.7), (ld + 52, 1), arrowprops=dict(arrowstyle="<->"))
    ax.text(ld + 53, -1.35, "4.7 H", va="center")
    ax.text(-19.5, 1.3, "inlet: u = 1 m/s, Tu 5 %, ν_t/ν 10", va="bottom", fontsize=8)
    ax.text(ld + 50, 1.3, "outlet: static p = 0", ha="right", va="bottom", fontsize=8)
    ax.text(30, 1.3, "upper wall (no slip, y = 1)", ha="center", va="bottom", fontsize=8)
    ax.text(18, -2.2, "lower wall (no slip)", ha="center", fontsize=8, rotation=-10)
    ax.set_xlim(-26, ld + 62); ax.set_ylim(-6, 3.4); ax.set_aspect("auto")
    ax.set_title("Geometry (units of H; Re_H = U_b H/ν = 20 000, ν = 5e-5 m²/s; z: symmetry planes)")
    for a, (x0, x1, y0, y1, ttl) in zip(axs[1:], [(-3, 12, -2.6, 1.4, "Mesh at the start of the diffuser"),
                                              (-0.6, 0.6, -0.12, 1.05, "Inlet-duct cross-section and wall clustering at x = 0 (zoom)")]):
        segs = []
        for q in quads:
            p = pts[q]
            if p[:, 0].max() < x0 or p[:, 0].min() > x1:
                continue
            segs.append(np.vstack([p, p[0]]))
        a.add_collection(LineCollection(segs, colors="k", linewidths=0.3))
        a.set_xlim(x0, x1); a.set_ylim(y0, y1); a.set_title(ttl, fontsize=9); a.set_aspect("auto")
    fig.tight_layout()
    fig.savefig(out)


if __name__ == "__main__":
    main(sys.argv[1])
