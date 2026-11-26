"""A sphere cut out of box meshes, L13 stage 0 (D225 item 5): the volume
fraction and face apertures of the cut cells on a slice, and the
convergence of the exact cut against a linear level-set cut (gallery
`docs/gallery/immersed.md`).

    pip install matplotlib numpy
    cargo run --release -p foehncraft-mesh --example cut_sphere -- output/cut
    python3 -m tools.plot.cut_sphere output/cut docs/figures/gallery/immersed_cut_sphere.png
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Circle

from tools.plot.colormaps import RAINBOW

CENTRE, R = (0.51, 0.47), 0.3


def main(run, png):
    s = np.loadtxt(os.path.join(run, "slice.csv"), delimiter=",", skiprows=1)
    w = np.loadtxt(os.path.join(run, "walls.csv"), delimiter=",", skiprows=1)
    conv = json.load(open(os.path.join(run, "convergence.json")))
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9), constrained_layout=True)

    ax = axes[0]
    x0, y0, h, frac, ap = s.T
    quads = [[(a, b), (a + c, b), (a + c, b + c), (a, b + c)] for a, b, c in zip(x0, y0, h)]
    pc = PolyCollection(quads, array=frac, cmap=RAINBOW, edgecolors="#00000030", linewidths=0.4)
    pc.set_clim(0, 1)
    ax.add_collection(pc)
    ax.add_patch(Circle(CENTRE, np.sqrt(R * R - (0.5 - 0.5 / 24 - 0.5) ** 2), fill=False, ec="w", lw=1.2, ls="--"))
    # wall pieces summed per cell: area vector at the area-weighted centroid
    cell = np.floor(w[:, :2] * 24).astype(int)
    keys, inv = np.unique(cell, axis=0, return_inverse=True)
    inv = inv.ravel()
    area = np.zeros((len(keys), 3))
    first = np.zeros((len(keys), 2))
    mag = np.linalg.norm(w[:, 3:6], axis=1)
    np.add.at(area, inv, w[:, 3:6])
    np.add.at(first, inv, w[:, :2] * mag[:, None])
    tot = np.zeros(len(keys))
    np.add.at(tot, inv, mag)
    pts = first / tot[:, None]
    a = area[:, :2]
    tip = pts + 0.03 * a / np.maximum(np.linalg.norm(area, axis=1, keepdims=True), 1e-30)
    ax.add_collection(LineCollection(np.stack([pts, tip], 1), colors="k", linewidths=0.9))
    ax.plot(pts[:, 0], pts[:, 1], "k.", ms=2)
    ax.set_xlim(0.1, 0.92)
    ax.set_ylim(0.06, 0.88)
    ax.set_aspect("equal")
    fig.colorbar(pc, ax=ax, shrink=0.85, label="volume fraction inside the body")
    ax.set_title("cut cells on 24³ hexahedra (slice), wall\narea vector per cell, exact circle dashed")

    ax = axes[1]
    m = (ap > 0) & (ap < 1)
    quads_m = [q for q, k in zip(quads, m) if k]
    pc2 = PolyCollection(quads_m, array=ap[m], cmap="inferno", edgecolors="#00000030", linewidths=0.4)
    pc2.set_clim(0, 1)
    ax.add_collection(pc2)
    ax.add_patch(Circle(CENTRE, np.sqrt(R * R - (0.5 - 0.5 / 24 - 0.5) ** 2), fill=False, ec="#5fd0ff", lw=1.2, ls="--"))
    ax.set_xlim(0.1, 0.92)
    ax.set_ylim(0.06, 0.88)
    ax.set_aspect("equal")
    ax.set_facecolor("#20242c")
    fig.colorbar(pc2, ax=ax, shrink=0.85, label="aperture of the low-x face (cut faces only)")
    ax.set_title("face apertures: part of each face inside,\nthe same number for both cells of a face")

    ax = axes[2]
    rows = conv["rows"]
    n = np.array([r["n"] for r in rows], float)
    for key, label, colour, mk in [("max_linear", "linear level set, max", "#d0402b", "o"),
                                   ("mean_linear", "linear level set, mean", "#d0402b", "s"),
                                   ("max_exact", "exact cut, max", "#2b7bd0", "o"),
                                   ("mean_exact", "exact cut, mean", "#2b7bd0", "s")]:
        y = np.array([r[key] for r in rows])
        ax.loglog(n, y, marker=mk, color=colour, ls="-" if "max" in key else "--", label=label)
    ax.text(9, 2.2e-4, "exact cut: the facet error\nof the STL, ∝ 1/h", color="#2b7bd0", fontsize=8)
    ax.set_xlabel("cells per side")
    ax.set_ylabel("volume fraction error of the cut cells")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_title(f"against the analytic sphere (icosphere of\n{conv['triangles']} triangles for the exact cut)")
    fig.savefig(png, dpi=110)
    print("wrote", png)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
