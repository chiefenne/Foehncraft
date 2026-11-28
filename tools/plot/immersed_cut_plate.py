"""Wall functions on a cut body, L13 stage 2c (gallery and recipe entry
`cut-plate`).

immersed_cut_plate.png: a flat plate tilted by 2° against a hexahedral
mesh, cut out of the mesh as a half-space (`shape = "plane"`), turbulent
(k-ε, scalable wall functions, Re_L = 6·10⁶): the speed near the plate
with the cut cells outlined, the skin friction along the plate against
White's correlation and against the same plate fitted to the mesh, and the
y⁺ of every cut wall face (the y⁺ map of the body). Read from the output
directories of

    target/release/foehn flow examples/cases/flow/immersed_cut_plate_keps.toml --out CUT
    target/release/foehn flow examples/cases/flow/immersed_cut_plate_fitted.toml --out FIT

and, for the comparison without the sampled wall function, of the cut case
with `[immersed.cut] wall_sample = false` (`--set immersed.cut.wall_sample=false`):

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_cut_plate OUT_DIR CUT FIT [NOSAMPLE]
"""
import math
import os
import sys

import numpy as np

NU = 2.0e-7
ANGLE = math.radians(2.0)
#: The plate is the plane through (0, Y0) tilted by ANGLE; s is the distance
#: along it from the inlet (x = 0).
Y0 = 0.02


def read_wall(path):
    """The columns of a wall CSV of Foehncraft (comment lines start with #)."""
    rows = [line for line in open(path) if not line.startswith("#")]
    head = rows[0].strip().split(",")
    data = np.array([[float(v) for v in r.split(",")] for r in rows[1:]])
    return {k: data[:, i] for i, k in enumerate(head)}


def white(s):
    """White's turbulent flat-plate local skin friction, fully turbulent from
    the leading edge: c_f = 0.455 / ln²(0.06 Re_x) (White, Viscous Fluid
    Flow, 3rd ed., 2006, eq. 6-78)."""
    return 0.455 / np.log(0.06 * s / NU) ** 2


def wall_cf(path, tilted):
    """s along the plate, c_f = 2 |τ_w| / (ρ U²) of the faces, the y⁺ of the
    faces, all sorted by s (one value per face column of the 2D mesh)."""
    d = read_wall(path)
    s = d["x"] / math.cos(ANGLE) if tilted else d["x"]
    cf = 2.0 * np.hypot(d["tau_x"], d["tau_y"])
    # the mesh is one cell thick: keep one face per column (z near ½)
    o = np.argsort(s)
    return s[o], cf[o], d["y_plus"][o]


def read_section(path):
    """The cells of a 2D run (one cell layer in z) as polygons on its lower
    z face, the cell of each polygon and the cell data: a polyhedron gives
    its faces that lie on that plane (several for a merged cell)."""
    import h5py
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]["internal"]
        pts = g["Points"][:]
        types = g["Types"][:]
        off = g["Offsets"][:]
        conn = g["Connectivity"][:]
        foff = g["FaceOffsets"][:]
        fconn = g["FaceConnectivity"][:]
        poff = g["PolyhedronOffsets"][:]
        p2f = g["PolyhedronToFaces"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    z0 = pts[:, 2].min()
    polys, cells = [], []
    for c in range(len(types)):
        if types[c] == 42:
            rings = [fconn[foff[k]:foff[k + 1]] for k in p2f[poff[c]:poff[c + 1]]]
            # a merged cell keeps the faces of its pieces: draw each one
            for r in rings:
                if np.all(np.abs(pts[r, 2] - z0) < 1e-9):
                    polys.append(pts[r, :2])
                    cells.append(c)
        else:
            # an extruded cell: its nodes on z = z0, ordered around their centre
            q = pts[conn[off[c]:off[c + 1]]]
            q = q[np.abs(q[:, 2] - z0) < 1e-9, :2]
            m = q.mean(axis=0)
            polys.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
            cells.append(c)
    return polys, np.array(cells), data


def plot(out, cut, fit, nosample=None):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    from tools.plot.colormaps import RAINBOW

    fig = plt.figure(figsize=(15.0, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0])

    # the speed near the plate, the cut cells outlined, the exact plane dashed
    ax = fig.add_subplot(gs[0])
    polys, cells, d = read_section(os.path.join(cut, "fields.vtkhdf"))
    x0, x1, y0, y1 = 0.90, 1.02, 0.050, 0.078
    keep = [i for i, q in enumerate(polys)
            if q[:, 0].max() > x0 and q[:, 0].min() < x1 and q[:, 1].min() < y1]
    sel = [polys[i] for i in keep]
    speed = d["speed"][cells[keep]]
    pc = PolyCollection(sel, array=speed, cmap=RAINBOW, edgecolors=(0, 0, 0, 0.35),
                        linewidths=0.4)
    pc.set_clim(0.0, 1.0)
    ax.add_collection(pc)
    xs = np.array([x0, x1])
    ax.plot(xs, Y0 + xs * math.tan(ANGLE), "w--", lw=1.2, label="exact plate")
    ax.fill_between(xs, y0, Y0 + xs * math.tan(ANGLE), color="0.35", zorder=0)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m), stretched")
    ax.set_title("speed near the trailing part of the plate: the cells cut\n"
                 "along the 2° plate (merged where < 0.3 of a cell)", fontsize=11)
    ax.legend(loc="upper left", fontsize=8)
    fig.colorbar(pc, ax=ax, shrink=0.9, pad=0.01, label="speed / U∞")

    # the skin friction along the plate
    ax = fig.add_subplot(gs[1])
    s = np.linspace(0.05, 1.2, 200)
    ax.plot(s, 1e3 * white(s), "k-", lw=1.6, label="White (2006), turbulent from the edge")
    fs, fcf, _ = wall_cf(os.path.join(fit, "wall_ymin.csv"), False)
    ax.plot(fs, 1e3 * fcf, "-", color="#30638e", lw=1.4, label="body-fitted plate, same mesh")
    if nosample:
        ns, ncf, _ = wall_cf(os.path.join(nosample, "wall_plate.csv"), True)
        ax.plot(ns, 1e3 * ncf, ".", color="#edae49", ms=3, label="cut, wall cell itself")
    cs, ccf, cyp = wall_cf(os.path.join(cut, "wall_plate.csv"), True)
    ax.plot(cs, 1e3 * ccf, ".", color="#d1495b", ms=3, label="cut, sampled one layer out")
    ax.set_xlim(0.0, 1.2)
    ax.set_ylim(1.0, 6.0)
    ax.set_xlabel("s along the plate (m)")
    ax.set_ylabel("c_f × 10³")
    ax.set_title("skin friction of every wall face\n(Re_s = s U∞ / ν up to 6·10⁶)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # the y+ map of the body
    ax = fig.add_subplot(gs[2])
    sc = ax.scatter(cs, cyp, c=cyp, cmap=RAINBOW, s=7, vmin=0, vmax=150)
    _, _, fyp = wall_cf(os.path.join(fit, "wall_ymin.csv"), False)
    ax.plot(fs, fyp, "-", color="#30638e", lw=1.2, label="body-fitted plate")
    ax.axhline(30.0, color="k", lw=0.8, ls=":", label="y⁺ = 30")
    ax.set_xlim(0.0, 1.2)
    ax.set_ylim(0, 160)
    ax.set_xlabel("s along the plate (m)")
    ax.set_ylabel("y⁺ of the first cell (exact distance to the plate)")
    ax.set_title(f"y⁺ map of the cut plate: {len(cs)} wall faces,\n"
                 f"min {cyp.min():.0f}, mean {cyp.mean():.0f}, max {cyp.max():.0f} (LIM011)",
                 fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.colorbar(sc, ax=ax, shrink=0.9, pad=0.01, label="y⁺")
    path = os.path.join(out, "immersed_cut_plate.png")
    fig.savefig(path, dpi=100)
    print(path)
    for at in (0.3, 0.5, 0.7, 0.97, 1.15):
        m = np.abs(cs - at) < 0.015
        mf = np.abs(fs - at) < 0.015
        w = white(at)
        print(f"s {at:4.2f}: White {1e3 * w:.3f}, fitted {100 * (fcf[mf].mean() / w - 1):+5.1f} %,"
              f" cut {100 * (ccf[m].mean() / w - 1):+5.1f} %"
              f" ({100 * (ccf[m].min() / w - 1):+5.1f} to {100 * (ccf[m].max() / w - 1):+5.1f} %)")


def main(argv):
    if len(argv) < 4:
        sys.exit(__doc__)
    import matplotlib

    matplotlib.use("Agg")
    plot(argv[1], argv[2], argv[3], argv[4] if len(argv) > 4 else None)


if __name__ == "__main__":
    main(sys.argv)
