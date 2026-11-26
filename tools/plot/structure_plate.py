"""Pictures of the reed petal as a shell (gallery entry of MS5 shells, L53).

structure_plate_modes.png: the first four modes of the 30 mm by 12 mm reed on clamp
springs (examples/cases/structure/reed_plate.toml, MITC4), each drawn as
the moved shell surface over the grey flat reed. Colour: the size of the
displacement (scaled mode shape), Rainbow Desaturated.

    pip install matplotlib numpy
    foehn call structure.run --input examples/cases/tools/structure_reed_plate.json
    python3 -m tools.plot.structure_plate modes output/reed_plate docs/figures/gallery/structure_plate_modes.png
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


def read_vtk(path):
    """Points, polygons and the point vectors of a legacy VTK polydata file."""
    with open(path) as f:
        tok = f.read().split("\n")
    pts, polys, vec = None, [], {}
    i = 0
    while i < len(tok):
        t = tok[i].split()
        if t and t[0] == "POINTS":
            n = int(t[1])
            pts = np.array([[float(v) for v in tok[i + 1 + k].split()] for k in range(n)])
            i += n
        elif t and t[0] in ("POLYGONS", "LINES"):
            m = int(t[1])
            cells = [[int(v) for v in tok[i + 1 + k].split()[1:]] for k in range(m)]
            if t[0] == "POLYGONS":
                polys = cells
            i += m
        elif t and t[0] == "VECTORS":
            n = len(pts)
            vec[t[1]] = np.array([[float(v) for v in tok[i + 1 + k].split()] for k in range(n)])
            i += n
        i += 1
    return pts, polys, vec


def need(vec, names, path):
    missing = [n for n in names if n not in vec]
    if missing:
        sys.exit(f"{path} lacks {missing}; it has: {', '.join(sorted(vec))}")


def surface(ax, pts, polys, u, scale, norm, label=None):
    """The flat reed in grey and the moved one coloured by |u|."""
    mm = 1e3
    flat = [pts[p] * mm for p in polys]
    ax.add_collection3d(
        Poly3DCollection(flat, facecolor="0.85", edgecolor="none", alpha=0.35)
    )
    moved = pts + scale * u
    faces = [moved[p] * mm for p in polys]
    mag = np.linalg.norm(u, axis=1)
    val = np.array([mag[p].mean() for p in polys])
    pc = Poly3DCollection(
        faces, facecolor=RAINBOW(norm(val)), edgecolor=(0, 0, 0, 0.15), linewidth=0.2
    )
    ax.add_collection3d(pc)
    allp = np.vstack([pts * mm, moved * mm])
    lo, hi = allp.min(axis=0), allp.max(axis=0)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ext = (hi - lo) / max(hi - lo)
    ax.set_box_aspect((ext[0], ext[1], max(ext[2], 0.3)))
    ax.locator_params(axis="z", nbins=4)
    ax.locator_params(axis="y", nbins=4)
    ax.set_xlabel("x (mm)", labelpad=-4)
    ax.set_ylabel("y (mm)", labelpad=-4)
    ax.set_zlabel("z (mm)", labelpad=-4)
    ax.tick_params(labelsize=7, pad=-2)
    if label:
        ax.set_title(label, fontsize=10, pad=0)


ORDINAL = ["first", "second", "third", "fourth", "fifth"]


def kind(pts, u):
    """Bending or twist, from the lift of the two long edges, and its order
    from the sign changes of the lift along the reed."""
    y = pts[:, 1]
    lo, hi = np.isclose(y, y.min()), np.isclose(y, y.max())
    a = u[lo][np.argsort(pts[lo, 0]), 2]
    b = u[hi][np.argsort(pts[hi, 0]), 2]
    twist = abs(a[-1] + b[-1]) < 0.5 * abs(a[-1] - b[-1])
    line = (b - a) if twist else (a + b)
    big = np.abs(line).max()
    s = np.sign(line[np.abs(line) > 1e-3 * big])
    order = 1 + int(np.count_nonzero(s[1:] != s[:-1]))
    return f"{ORDINAL[min(order, 5) - 1]} {'twist' if twist else 'bending'}"


def modes(out_dir, fig_path):
    path = os.path.join(out_dir, "structure.vtk")
    pts, polys, vec = read_vtk(path)
    want = [f"mode_{k}" for k in (1, 2, 3, 4)]
    need(vec, want, path)
    with open(os.path.join(out_dir, "result.json")) as f:
        res = json.load(f)
    freqs = [m["frequency"] for m in res["modal"]["modes"]]
    size = max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))
    fig = plt.figure(figsize=(10.0, 7.0))
    fig.subplots_adjust(left=0.0, right=0.88, bottom=0.03, top=0.9, wspace=0.05, hspace=0.12)
    norm = Normalize(0.0, 1.0)
    for k, name in enumerate(want):
        ax = fig.add_subplot(2, 2, k + 1, projection="3d")
        ax.view_init(elev=24, azim=-58)
        u = vec[name]
        surface(
            ax, pts, polys, u, 0.25 * size, norm,
            f"mode {k + 1}, {kind(pts, u)}: {freqs[k]:.0f} Hz",
        )
    sm = plt.cm.ScalarMappable(norm=norm, cmap=RAINBOW)
    cax = fig.add_axes((0.9, 0.25, 0.02, 0.5))
    fig.colorbar(sm, cax=cax, label="|displacement| (scaled mode shape)")
    fig.suptitle(
        "Reed petal 30 × 12 × 0.2 mm on clamp springs, MITC4 shells: the first four modes",
        fontsize=11,
    )
    fig.savefig(fig_path, dpi=110)
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] != "modes":
        sys.exit(__doc__)
    modes(sys.argv[2], sys.argv[3])
