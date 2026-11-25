"""Pictures of the curved reed stop as a solid (gallery entry of MS5 solids,
L53).

structure_reed_stop.png: the stop of examples/cases/structure/reed_stop.toml
(twenty-node hexahedra from a Gmsh file). Left: bent by 0.3 MPa on its
underside, the displacement drawn 15 times larger over the grey unloaded
stop, coloured by the von Mises stress (Rainbow Desaturated). Right: the
first two modes, coloured by the size of the displacement.

    pip install matplotlib numpy
    foehn call structure.run --input examples/cases/tools/structure_reed_stop.json
    python3 -m tools.plot.structure_solid output/reed_stop docs/figures/gallery/structure_reed_stop.png
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from tools.plot.colormaps import colormap

RAINBOW = colormap()

# faces of the quadratic cells (VTK order): corners, then edge nodes
HEX_FACES = [
    [0, 3, 2, 1, 11, 10, 9, 8],
    [4, 5, 6, 7, 12, 13, 14, 15],
    [0, 1, 5, 4, 8, 17, 12, 16],
    [2, 3, 7, 6, 10, 19, 14, 18],
    [0, 4, 7, 3, 16, 15, 19, 11],
    [1, 2, 6, 5, 9, 18, 13, 17],
]
TET_FACES = [[0, 2, 1, 6, 5, 4], [0, 1, 3, 4, 8, 7], [1, 2, 3, 5, 9, 8], [0, 3, 2, 7, 9, 6]]


def read_vtk(path):
    """Points, cells, cell types and point data of a legacy VTK
    unstructured grid (the structure.vtk of a case with solids)."""
    with open(path) as f:
        lines = f.read().split("\n")
    pts, cells, types, vec, sca = None, [], [], {}, {}
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if t and t[0] == "DATASET" and t[1] != "UNSTRUCTURED_GRID":
            sys.exit(f"{path} has no solids (DATASET {t[1]})")
        if t and t[0] == "POINTS":
            n = int(t[1])
            pts = np.array([[float(v) for v in lines[i + 1 + k].split()] for k in range(n)])
            i += n
        elif t and t[0] == "CELLS":
            m = int(t[1])
            cells = [[int(v) for v in lines[i + 1 + k].split()[1:]] for k in range(m)]
            i += m
        elif t and t[0] == "CELL_TYPES":
            m = int(t[1])
            types = [int(lines[i + 1 + k]) for k in range(m)]
            i += m
        elif t and t[0] == "VECTORS":
            n = len(pts)
            vec[t[1]] = np.array([[float(v) for v in lines[i + 1 + k].split()] for k in range(n)])
            i += n
        elif t and t[0] == "SCALARS":
            n = len(pts)
            sca[t[1]] = np.array([float(lines[i + 2 + k]) for k in range(n)])
            i += n + 1
        i += 1
    return pts, cells, types, vec, sca


def boundary(cells, types):
    """The faces on the outside of the solid cells, as node lists in the
    order of their outline (corner, edge node, corner, …)."""
    seen = {}
    for c, t in zip(cells, types):
        table = HEX_FACES if t == 25 else TET_FACES if t == 24 else []
        for f in table:
            nodes = [c[a] for a in f]
            k = len(nodes) // 2
            key = tuple(sorted(nodes[:k]))
            outline = [x for pair in zip(nodes[:k], nodes[k:]) for x in pair]
            seen.setdefault(key, []).append(outline)
    return [v[0] for v in seen.values() if len(v) == 1]


def need(have, names, path):
    missing = [n for n in names if n not in have]
    if missing:
        sys.exit(f"{path} lacks {missing}; it has: {', '.join(sorted(have))}")


def draw(ax, pts, faces, u, scale, value, norm):
    mm = 1e3
    ax.add_collection3d(
        Poly3DCollection([pts[f] * mm for f in faces], facecolor="0.8", edgecolor="none", alpha=0.25)
    )
    moved = pts + scale * u
    val = np.array([value[f].mean() for f in faces])
    ax.add_collection3d(
        Poly3DCollection(
            [moved[f] * mm for f in faces],
            facecolor=RAINBOW(norm(val)),
            edgecolor=(0, 0, 0, 0.18),
            linewidth=0.25,
        )
    )
    allp = np.vstack([pts * mm, moved * mm])
    lo, hi = allp.min(axis=0), allp.max(axis=0)
    mid, half = (lo + hi) / 2, (hi - lo).max() / 2
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-62)
    ax.set_xlabel("x (mm)", fontsize=8, labelpad=-4)
    ax.set_ylabel("y (mm)", fontsize=8, labelpad=-4)
    ax.set_zlabel("z (mm)", fontsize=8, labelpad=-4)
    ax.tick_params(labelsize=7, pad=-2)
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)
    ax.set_facecolor((0, 0, 0, 0))


def main(out_dir, png):
    path = os.path.join(out_dir, "structure.vtk")
    pts, cells, types, vec, sca = read_vtk(path)
    need({**vec, **sca}, ["displacement", "von_mises", "mode_1", "mode_2"], path)
    with open(os.path.join(out_dir, "result.json")) as f:
        res = json.load(f)
    freqs = [m["frequency"] for m in res["modal"]["modes"]]
    faces = boundary(cells, types)
    fig = plt.figure(figsize=(11.0, 6.2), dpi=100)
    # static: von Mises in MPa
    ax = fig.add_axes([0.0, 0.1, 0.56, 0.78], projection="3d")
    vm = sca["von_mises"] / 1e6
    norm = Normalize(0.0, vm.max())
    draw(ax, pts, faces, vec["displacement"], 15.0, vm, norm)
    ax.set_title(
        "0.3 MPa on the underside: von Mises stress (displacement × 15)",
        fontsize=10,
    )
    cax = fig.add_axes([0.1, 0.07, 0.36, 0.022])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=RAINBOW), cax=cax, orientation="horizontal")
    cb.set_label("von Mises stress (MPa)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    tip = res["static"]["points"]["tip"]["displacement"]
    ax.text2D(
        0.03,
        0.93,
        f"tip lift {tip[2] * 1e3:.3f} mm, largest stress {vm.max():.0f} MPa at the root",
        transform=ax.transAxes,
        fontsize=8,
    )
    # modes
    for k, name in enumerate(["mode_1", "mode_2"]):
        ax = fig.add_axes([0.58, 0.47 - 0.45 * k, 0.4, 0.4], projection="3d")
        u = vec[name]
        mag = np.linalg.norm(u, axis=1)
        draw(ax, pts, faces, u, 4e-3, mag, Normalize(0.0, mag.max()))
        ax.set_title(f"mode {k + 1}: {freqs[k]:.0f} Hz (colour: size of the displacement)", fontsize=9, pad=2)
    fig.suptitle(
        "A curved reed stop as a solid: 160 twenty-node hexahedra read from Gmsh",
        fontsize=12,
        y=0.985,
    )
    os.makedirs(os.path.dirname(png) or ".", exist_ok=True)
    fig.savefig(png, dpi=100)
    print(f"wrote {png}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
