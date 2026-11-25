"""Pictures of the heat sink warmed by its chip as a structure (gallery
entry of MS6, L53): the temperatures of the solid region of the turbulent
conjugate heat transfer case mapped onto the twenty-node hexahedra of the
structure.

structure_heat_sink_hot.png: top left, the mapped temperature on the
structure (its own mesh); top right, the base bonded to the board at its
underside, the displacement drawn 300 times larger over the grey shape,
coloured by the von Mises stress (Rainbow Desaturated); bottom, the
temperature along the top of the base: the flow case's cells (its line
output) and the structural nodes on that line.

    pip install matplotlib numpy
    foehn flow examples/cases/flow/cht_heat_sink_turbulent.toml --out output/cht_heat_sink_turbulent
    foehn call structure.run --input examples/cases/tools/structure_heat_sink_hot.json
    python3 -m tools.plot.structure_heat_sink output/cht_heat_sink_turbulent output/heat_sink_hot docs/figures/gallery/structure_heat_sink_hot.png
"""
import csv
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
from tools.plot.structure_solid import boundary, need, read_vtk

RAINBOW = colormap()
# drawn with the model's y (up the fins) as the vertical axis
SWAP = [0, 2, 1]


def draw(ax, pts, faces, u, scale, value, norm, ghost):
    mm = 1e3
    p = pts[:, SWAP] * mm
    if ghost:
        ax.add_collection3d(
            Poly3DCollection([p[f] for f in faces], facecolor="0.75", edgecolor="none", alpha=0.25)
        )
    moved = (pts + scale * u)[:, SWAP] * mm
    val = np.array([value[f].mean() for f in faces])
    ax.add_collection3d(
        Poly3DCollection(
            [moved[f] for f in faces],
            facecolor=RAINBOW(norm(val)),
            edgecolor=(0, 0, 0, 0.15),
            linewidth=0.2,
        )
    )
    lo, hi = p.min(axis=0), p.max(axis=0)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1] - 6, hi[1] + 6)
    ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1] + 12, hi[2] - lo[2]), zoom=1.1)
    ax.view_init(elev=20, azim=-62)
    ax.set_xlabel("x (mm)", fontsize=8, labelpad=2)
    ax.set_zlabel("y (mm)", fontsize=8, labelpad=-2)
    ax.set_yticks([])
    ax.tick_params(labelsize=7, pad=0)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_alpha(0.0)
    ax.set_facecolor((0, 0, 0, 0))


def colorbar(fig, rect, norm, label):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=RAINBOW), cax=cax, orientation="horizontal")
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)


def main(flow_dir, out_dir, png):
    path = os.path.join(out_dir, "structure.vtk")
    pts, cells, types, vec, sca = read_vtk(path)
    need({**vec, **sca}, ["displacement", "von_mises", "temperature"], path)
    with open(os.path.join(out_dir, "result.json")) as f:
        res = json.load(f)
    faces = boundary(cells, types)
    fig = plt.figure(figsize=(12.0, 7.6), dpi=100)

    ax = fig.add_axes([0.0, 0.36, 0.5, 0.56], projection="3d")
    t = sca["temperature"]
    tn = Normalize(t.min(), t.max())
    draw(ax, pts, faces, vec["displacement"], 0.0, t, tn, False)
    ax.set_title("temperature from the flow case on the structure's nodes", fontsize=10, pad=0)
    colorbar(fig, [0.1, 0.4, 0.3, 0.018], tn, "temperature (K)")

    ax = fig.add_axes([0.5, 0.36, 0.5, 0.56], projection="3d")
    vm = sca["von_mises"] / 1e6
    vn = Normalize(0.0, np.percentile(vm, 99.5))
    draw(ax, pts, faces, vec["displacement"], 300.0, vm, vn, True)
    ax.set_title("bonded at its underside: von Mises stress (displacement × 300)", fontsize=10, pad=0)
    colorbar(fig, [0.6, 0.4, 0.3, 0.018], vn, "von Mises stress (MPa; the colours end at the 99.5 % value)")
    tip = res["static"]["points"]["fin_3_tip"]["displacement"]
    st = res["static"]
    ax.text2D(
        0.02,
        0.9,
        f"fin tip rises {tip[1] * 1e6:.1f} µm; largest stress {st['max_von_mises']['value'] / 1e6:.0f} MPa "
        "at the end of the bond line",
        transform=ax.transAxes,
        fontsize=8,
    )

    # along the top of the base: the flow cells and the structural nodes
    ax = fig.add_axes([0.08, 0.07, 0.86, 0.24])
    with open(os.path.join(flow_dir, "line_base_top.csv")) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    xf = np.array([float(r["x"]) for r in rows]) * 1e3
    tf = np.array([float(r["temperature"]) for r in rows])
    ax.plot(xf, tf, color="#1f77b4", lw=1.6, label="flow case: cells along y = 4.5 mm (line output)")
    on = (np.abs(pts[:, 1] - 0.005) < 1e-9) & (np.abs(pts[:, 2] - 0.005) < 1e-9)
    order = np.argsort(pts[on, 0])
    ax.plot(
        pts[on, 0][order] * 1e3,
        t[on][order],
        "o",
        ms=3.2,
        color="#d62728",
        label="structure: mapped nodes on the top of the base (y = 5 mm)",
    )
    ax.set_xlabel("x (mm)", fontsize=9)
    ax.set_ylabel("temperature (K)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower center")
    m = res["temperature_fields"][0]
    fig.suptitle(
        f"A heat sink warmed by its chip: {m['source_cells']} cells of the aluminium mapped onto "
        f"{m['nodes']} nodes of twenty-node hexahedra",
        fontsize=12,
        y=0.985,
    )
    os.makedirs(os.path.dirname(png) or ".", exist_ok=True)
    fig.savefig(png, dpi=100)
    print(f"wrote {png}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
