"""Pictures of the structural dynamics and the thermal load (gallery entry
of MS5 part 3, L53).

structure_reed_pulse.png: the reed petal of reed_plate_pulse.toml struck by
a 1 ms pressure pulse: four snapshots of the moving shell coloured by its
lift (a signed field centred on zero: Rainbow Uniform), and the lift of
the tip against time with the snapshots marked.

structure_reed_hot.png: the reed of reed_plate_hot.toml, 80 K warmer and
6 K hotter on top, bent by its thermal strain; colour: the lift (mm),
Rainbow Desaturated, heights stretched.

    pip install matplotlib numpy
    foehn call structure.run --input examples/cases/tools/structure_reed_plate_pulse.json
    foehn call structure.run --input examples/cases/tools/structure_reed_plate_hot.json
    python3 -m tools.plot.structure_dynamics pulse output/reed_plate_pulse docs/figures/gallery/structure_reed_pulse.png
    python3 -m tools.plot.structure_dynamics hot output/reed_plate_hot docs/figures/gallery/structure_reed_hot.png
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

from tools.plot.colormaps import colormap, uniform_colormap
from tools.plot.structure_plate import read_vtk

RAINBOW = colormap()
UNIFORM = uniform_colormap()
MM = 1e3


def need(vec, names, path):
    missing = [n for n in names if n not in vec]
    if missing:
        sys.exit(f"{path} lacks {missing}; it has: {', '.join(sorted(vec))}")


def draw(ax, pts, polys, u, value, cmap, norm, zlim=None, zaspect=0.3):
    flat = [pts[p] * MM for p in polys]
    ax.add_collection3d(Poly3DCollection(flat, facecolor="0.85", edgecolor="none", alpha=0.3))
    moved = pts + u
    faces = [moved[p] * MM for p in polys]
    val = np.array([value[p].mean() for p in polys])
    ax.add_collection3d(
        Poly3DCollection(faces, facecolor=cmap(norm(val)), edgecolor=(0, 0, 0, 0.12), linewidth=0.2)
    )
    allp = np.vstack([pts * MM, moved * MM])
    lo, hi = allp.min(axis=0), allp.max(axis=0)
    if zlim:
        lo[2], hi[2] = zlim
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ext = (hi - lo) / max(hi - lo)
    ax.set_box_aspect((ext[0], ext[1], zaspect))
    ax.locator_params(axis="z", nbins=3)
    ax.locator_params(axis="y", nbins=3)
    ax.tick_params(labelsize=7, pad=-2)


def pulse(out_dir, fig_path):
    path = os.path.join(out_dir, "structure.vtk")
    pts, polys, vec = read_vtk(path)
    with open(os.path.join(out_dir, "result.json")) as f:
        tr = json.load(f)["transient"]
    snaps = tr["snapshots"]
    names = [f"transient_{i + 1}" for i in range(len(snaps))]
    need(vec, names, path)
    big = max(np.abs(vec[n][:, 2]).max() for n in names)
    norm = Normalize(-big * MM, big * MM)
    fig = plt.figure(figsize=(11.0, 7.4))
    zl = (-big * MM * 1.05, big * MM * 1.05)
    for k, (name, t) in enumerate(zip(names, snaps)):
        ax = fig.add_axes((0.0 + 0.226 * k, 0.36, 0.235, 0.58), projection="3d")
        ax.set_facecolor((1, 1, 1, 0))
        ax.view_init(elev=18, azim=-62)
        u = vec[name]
        draw(ax, pts, polys, u, u[:, 2] * MM, UNIFORM, norm, zl, zaspect=0.75)
        ax.set_title(f"t = {t * 1e3:.1f} ms", fontsize=10, pad=-12)
        ax.set_xlabel("x (mm)", fontsize=8, labelpad=-6)
        ax.set_zlabel("lift (mm)", fontsize=8, labelpad=-6)
    cax = fig.add_axes((0.925, 0.46, 0.012, 0.38))
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=UNIFORM), cax=cax, label="lift (mm)")
    ax = fig.add_axes((0.08, 0.08, 0.84, 0.28))
    t = np.array(tr["times"]) * 1e3
    w = np.array([p[2] for p in tr["points"]["tip"]["displacement"]]) * MM
    ax.axvspan(0.0, 1.0, color="#FFE8A0", label="pressure pulse, 3 kPa peak")
    ax.plot(t, w, color="#2A3F9E", lw=1.6, label="tip lift, generalised-α")
    tw = np.interp(np.array(snaps) * 1e3, t, w)
    ax.scatter(np.array(snaps) * 1e3, tw, c=UNIFORM(norm(tw)), s=60, edgecolor="k", zorder=3, label="snapshots")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("tip lift (mm)")
    ax.set_xlim(t[0], t[-1])
    ax.legend(loc="upper right", fontsize=8, ncol=3)
    ax.grid(alpha=0.3)
    fig.suptitle(
        "Reed petal on clamp springs struck by a 1 ms pressure pulse: it rings at 154 Hz, damped 2 %",
        fontsize=11,
    )
    fig.savefig(fig_path, dpi=110)
    print(f"wrote {fig_path}")


def hot(out_dir, fig_path):
    path = os.path.join(out_dir, "structure.vtk")
    pts, polys, vec = read_vtk(path)
    need(vec, ["displacement"], path)
    u = vec["displacement"]
    lift = u[:, 2] * MM
    norm = Normalize(lift.min(), lift.max())
    fig = plt.figure(figsize=(8.6, 5.2))
    ax = fig.add_axes((0.0, 0.02, 0.84, 0.88), projection="3d")
    ax.view_init(elev=22, azim=-55)
    draw(ax, pts, polys, u, lift, RAINBOW, norm)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("lift (mm)")
    cax = fig.add_axes((0.86, 0.22, 0.02, 0.55))
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=RAINBOW), cax=cax, label="lift (mm)")
    fig.suptitle(
        f"Reed 80 K warmer, 6 K hotter on top: bent towards its seat by {-lift.min():.2f} mm (heights stretched)",
        fontsize=11,
    )
    fig.savefig(fig_path, dpi=110)
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("pulse", "hot"):
        sys.exit(__doc__)
    {"pulse": pulse, "hot": hot}[sys.argv[1]](sys.argv[2], sys.argv[3])
