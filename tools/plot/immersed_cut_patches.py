"""Boundary conditions on cut bodies, L13 stage 2e (gallery and recipe
entry `cut-patches`).

immersed_cut_patches.png: a cylinder that blows through a slot in the
Schäfer–Turek channel, an STL body whose two named solids are two boundaries
with their own conditions (the speed, the cut faces of the slot and of the
wall near the slot), and Taylor–Couette flow with both cylinders cut out of a
plain box mesh (the inner one turning with `omega`, the outer one inverted):
the azimuthal velocity of every cell against the exact solution, and the
velocity error and the torque error against the mesh. Read from the output
directories of

    target/release/foehn flow examples/cases/flow/immersed_cut_blown_cylinder.toml --out BLOWN
    target/release/foehn flow examples/cases/flow/immersed_cut_taylor_couette.toml --out TC

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_cut_patches OUT_DIR BLOWN TC
"""
import json
import math
import os
import sys

import numpy as np

from tools.plot.immersed_cut_plate import read_section

#: Taylor–Couette: r1 turning at W, r2 at rest, μ, depth
R1, R2, W, MU, DEPTH = 0.05, 0.10, 10.0, 5e-4, 0.01
#: the blown cylinder and its slot (degrees)
CX, CY, RC = 0.2, 0.2, 0.05
SLOT = (45.0, 75.0)


def exact_ut(r):
    """u_θ = A r + B / r of the gap, outer cylinder at rest."""
    a = -W * R1 * R1 / (R2 * R2 - R1 * R1)
    b = W * R1 * R1 * R2 * R2 / (R2 * R2 - R1 * R1)
    return a * r + b / r


def boundary_faces(path, name):
    """The polygons of boundary `name` of a VTKHDF file (x, y of the faces
    on the lower z plane, as segments)."""
    import h5py
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"][name]
        pts = g["Points"][:]
        off = g["Polygons"]["Offsets"][:]
        conn = g["Polygons"]["Connectivity"][:]
    z0 = pts[:, 2].min()
    segs = []
    for k in range(len(off) - 1):
        q = pts[conn[off[k]:off[k + 1]]]
        q = q[np.abs(q[:, 2] - z0) < 1e-9, :2]
        if len(q) >= 2:
            segs.append(q)
    return segs


def plot(out, blown, tc):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PolyCollection

    from tools.plot.colormaps import RAINBOW

    fig = plt.figure(figsize=(14.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[2.2, 1.0])

    # the blown cylinder: speed in the channel
    ax = fig.add_subplot(gs[0, 0])
    path = os.path.join(blown, "fields.vtkhdf")
    polys, cells, d = read_section(path)
    keep = [i for i, q in enumerate(polys) if q[:, 0].min() < 1.0]
    sel = [polys[i] for i in keep]
    sp = d["speed"][cells[keep]]
    pc = PolyCollection(sel, array=sp, cmap=RAINBOW, edgecolors="none")
    pc.set_clim(0.0, max(0.45, float(sp.max())))
    ax.add_collection(pc)
    a = np.linspace(0.0, 2.0 * math.pi, 300)
    ax.fill(CX + RC * np.cos(a), CY + RC * np.sin(a), color="0.35", zorder=2)
    s = np.radians(np.linspace(*SLOT, 40))
    ax.plot(CX + RC * np.cos(s), CY + RC * np.sin(s), "-", color="#ff2d55", lw=3.0, zorder=3,
            label="patch `slot`: velocity inlet, 0.4 m/s")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.41)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("a cylinder that blows through a slot (Re 20): one STL body, two named\n"
                 "solids, two boundaries with their own conditions, on a plain box mesh",
                 fontsize=11)
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(pc, ax=ax, shrink=0.75, pad=0.01, label="speed (m/s)")

    # the slot: the cut faces of each patch
    ax = fig.add_subplot(gs[0, 1])
    near = [i for i, q in enumerate(polys)
            if math.hypot(q[:, 0].mean() - CX, q[:, 1].mean() - CY) < RC + 0.03]
    zc = PolyCollection([polys[i] for i in near], array=d["speed"][cells[near]], cmap=RAINBOW,
                        edgecolors=(0, 0, 0, 0.5), linewidths=0.5)
    zc.set_clim(pc.get_clim())
    ax.add_collection(zc)
    for name, col, lab in (("cylinder-wall", "0.15", "cylinder-wall (wall)"),
                           ("cylinder-slot", "#ff2d55", "cylinder-slot (inlet)")):
        segs = boundary_faces(path, name)
        ax.add_collection(LineCollection(segs, colors=col, linewidths=3.0, label=lab))
    ax.set_xlim(CX + 0.0, CX + 0.075)
    ax.set_ylim(CY + 0.0, CY + 0.075)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_title("the slot: every cut face takes the solid\nof the closest triangle", fontsize=11)
    ax.legend(loc="upper right", fontsize=8)

    # Taylor–Couette: u_θ of every cell, finest level
    ax = fig.add_subplot(gs[1, 0])
    polys, cells, d = read_section(os.path.join(tc, "fields.vtkhdf"))
    c = np.array([q.mean(axis=0) for q in polys])
    r = np.hypot(c[:, 0], c[:, 1])
    ux, uy = d["velocity_x"][cells], d["velocity_y"][cells]
    ut = (-c[:, 1] * ux + c[:, 0] * uy) / r
    rr = np.linspace(R1, R2, 200)
    ax.plot(rr, exact_ut(rr) / (W * R1), "k-", lw=1.2, label="exact, A r + B / r")
    ax.plot(r, ut / (W * R1), ".", color="#d1495b", ms=2.0, label="every cell, 176 cells across the box")
    ax.set_xlabel("r (m)")
    ax.set_ylabel("u_θ / (ω r₁)")
    ax.set_title("Taylor–Couette, both cylinders cut out of a box mesh: the inner one\n"
                 "turns (`omega`), the outer one is inverted (fluid inside)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # errors against the mesh
    ax = fig.add_subplot(gs[1, 1])
    res = json.load(open(os.path.join(tc, "result.json")))
    exact_t = 4 * math.pi * MU * W * R1 ** 2 * R2 ** 2 * DEPTH / (R2 ** 2 - R1 ** 2)
    n = np.array([44 * 2 ** i for i in range(len(res["solves"]))])
    l2 = np.array([s["errors"]["l2"] for s in res["solves"]])
    ti, to = [], []
    for s in res["solves"]:
        t = {b["name"]: b.get("torque") for b in s["boundaries"]}
        ti.append(abs(-t["inner"] / exact_t - 1.0))
        to.append(abs(t["outer"] / exact_t - 1.0))
    o = res["orders"]
    ax.loglog(n, l2, "o-", color="#30638e", label=f"velocity L2 (relative): order {o[0]:.2f}, {o[1]:.2f}")
    ax.loglog(n, ti, "s-", color="#d1495b", label="torque, inner cylinder")
    ax.loglog(n, to, "^-", color="#edae49", label="torque, outer cylinder")
    ax.loglog(n, l2[0] * (n / n[0]) ** -2.0, "k--", lw=0.9, label="second order")
    ax.set_xlabel("cells across the box")
    ax.set_ylabel("relative error")
    ax.set_title("errors against the exact solution", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    path = os.path.join(out, "immersed_cut_patches.png")
    fig.savefig(path, dpi=100)
    print(path)
    print("cells", list(n), "L2", list(l2), "orders", o)
    print("torque error inner", ti, "outer", to, "exact", exact_t)


def main(argv):
    if len(argv) < 4:
        sys.exit(__doc__)
    import matplotlib

    matplotlib.use("Agg")
    plot(argv[1], argv[2], argv[3])


if __name__ == "__main__":
    main(sys.argv)
