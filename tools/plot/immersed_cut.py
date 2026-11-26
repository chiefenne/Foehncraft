"""Cut cells, L13 stage 2 (gallery and recipe 313).

immersed_cut_cylinder.png: the Schäfer & Turek cylinder at Re = 20 on a
plain box mesh with the cylinder cut out of it: the speed, the vorticity
and a close-up of the cut mesh coloured by the pressure, the exact circle
dashed.

immersed_cut_benchmark.png: the mesh benchmark of the cut (a sphere cut out
of hexahedral, tetrahedral and polyhedral meshes): the error of the fluid
volume against the mesh size for linear crossings (marching cubes), exact
crossings and exact crossings with the surface point of dual contouring;
read from the Markdown table of `cargo run --release -p foehncraft-mesh
--example cut_bench`.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_cut OUT_DIR CUT_RUN [BENCH_TABLE.md]

with the output directory of

    target/release/foehn flow examples/cases/flow/immersed_cut_cylinder_re20.toml --out output/cut
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle

from tools.plot.colormaps import RAINBOW, UNIFORM

XMAX = 0.9
CENTRE, R = (0.2, 0.2), 0.05


def read(path):
    """The cells of a 2D run as their polygons on the zmin boundary, with
    the cell data."""
    import h5py
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]["internal"]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
        b = f["VTKHDF"]["zmin"]
        pts = b["Points"][:]
        off = b["Polygons"]["Offsets"][:]
        conn = b["Polygons"]["Connectivity"][:]
        owner = b["CellData"]["OwnerCell"][:]
        wall = f["VTKHDF"]["cylinder"]["Points"][:]
    polys = [pts[conn[off[i]:off[i + 1]], :2] for i in range(len(off) - 1)]
    return polys, owner, data, wall


def circle(ax, colour="w"):
    ax.add_patch(Circle(CENTRE, R, fill=False, ec=colour, lw=1.0, ls="--"))


def panel(fig, ax, polys, values, cmap, label, title, vmin, vmax, edges="face"):
    lw = 0.3 if edges != "face" else 0.0
    pc = PolyCollection(polys, array=values, cmap=cmap, edgecolors=edges, linewidths=lw,
                        antialiased=edges != "face")
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.colorbar(pc, ax=ax, label=label, shrink=0.9, pad=0.01)


def fields(out, run):
    polys, owner, d, _ = read(os.path.join(run, "fields.vtkhdf"))
    speed = np.hypot(d["velocity_x"], d["velocity_y"])[owner]
    wz = d["vorticity_z"][owner]
    p = d["pressure"][owner]
    xmin = np.array([q[:, 0].min() for q in polys])
    keep = xmin < XMAX
    sel = [q for q, k in zip(polys, keep) if k]
    fig = plt.figure(figsize=(12.5, 6.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.9, 1.0])
    a0 = fig.add_subplot(gs[0, 0])
    a1 = fig.add_subplot(gs[1, 0])
    a2 = fig.add_subplot(gs[:, 1])
    panel(fig, a0, sel, speed[keep], RAINBOW, "speed (m/s)",
          "Cut cells, box mesh h = D/20: speed", 0.0, speed[keep].max())
    w = np.percentile(np.abs(wz[keep]), 99.0)
    panel(fig, a1, sel, wz[keep], UNIFORM, "vorticity (1/s)", "Vorticity", -w, w)
    for a in (a0, a1):
        a.set_xlim(0.0, XMAX)
        a.set_ylim(0.0, 0.41)
    cx = np.array([q[:, 0].mean() for q in polys])
    cy = np.array([q[:, 1].mean() for q in polys])
    near = (np.abs(cx - CENTRE[0]) < 0.075) & (np.abs(cy - CENTRE[1]) < 0.075)
    pn = p[near]
    panel(fig, a2, [q for q, k in zip(polys, near) if k], pn, RAINBOW, "pressure (Pa)",
          "The cut mesh near the cylinder: pressure", pn.min(), pn.max(),
          edges=(0.15, 0.15, 0.15, 0.7))
    circle(a2, "k")
    a2.set_xlim(CENTRE[0] - 0.07, CENTRE[0] + 0.07)
    a2.set_ylim(CENTRE[1] - 0.07, CENTRE[1] + 0.07)
    path = os.path.join(out, "immersed_cut_cylinder.png")
    fig.savefig(path, dpi=110)
    print(path)


def bench(out, table):
    rows = []
    for line in open(table):
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(c) < 10 or not c[1].isdigit():
            continue
        rows.append((c[0], int(c[1]), c[3], float(c[4]), float(c[6])))
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    style = {"linear": ("#d1495b", "o"), "exact": ("#2e86ab", "s"),
             "exact+point": ("#3b8b3b", "^")}
    dash = {"hex": "-", "tet": "--", "poly": ":"}
    for ax, k, title in ((axs[0], 3, "fluid volume"), (axs[1], 4, "wall area")):
        for mesh in ("hex", "tet", "poly"):
            for m, (col, mk) in style.items():
                r = [x for x in rows if x[0] == mesh and x[2] == m]
                h = [2.0 / x[1] for x in r]
                ax.loglog(h, [x[k] for x in r], dash[mesh], color=col, marker=mk,
                          label=f"{mesh}, {m}")
        h = np.array([2 / 32, 2 / 8])
        ax.loglog(h, 3.0 * h**2, "k-", lw=0.8, label="slope 2")
        ax.set_xlabel("cell width h (sphere radius 0.55)")
        ax.set_ylabel(f"relative error of the {title}")
        ax.set_title(f"Sphere cut out of a box: {title}", fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
    axs[1].legend(fontsize=7, ncol=2)
    path = os.path.join(out, "immersed_cut_benchmark.png")
    fig.savefig(path, dpi=110)
    print(path)


def main():
    out, run = sys.argv[1:3]
    os.makedirs(out, exist_ok=True)
    fields(out, run)
    if len(sys.argv) > 3:
        bench(out, sys.argv[3])


if __name__ == "__main__":
    main()
