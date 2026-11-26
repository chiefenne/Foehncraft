"""STL bodies cut out of a mesh, L13 stage 2 (gallery and recipe 322).

immersed_stl_benchmark.png: the cut of three triangle surfaces (an
icosphere, a box turned against the mesh, a torus) out of hexahedral,
tetrahedral and polyhedral meshes of 8³ to 32³ cells, against the exact
moments of the same triangles: the error of the fluid volume and of the
wall area for marching cubes, exact crossings, the closest-point wall
point, the dual-contouring wall point and the face points; read from the
JSON of `cargo run --release -p foehncraft-mesh --example cut_stl_bench`.

immersed_stl_box_walls.png: the wall faces of the turned box cut out of a
hexahedral mesh, with the closest point (the edges and corners rounded off)
and with the dual-contouring wall point and the face points (kept sharp),
coloured by how far the face centroids lie from the surface; read from the
JSON of `cargo run --release -p foehncraft-mesh --example cut_stl_walls`.

immersed_stl_sphere_array.png: Stokes flow through a simple cubic array of
STL spheres cut out of the mesh: the pressure on the cut wall faces and the
speed on a slice behind the sphere; from the output directory of

    target/release/foehn flow examples/cases/flow/immersed_cut_stl_sphere_array.toml --out RUN

The STL file of the case is written by the `stl` command (an icosphere of
four subdivisions, 5120 triangles, binary):

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_stl stl examples/cases/flow/immersed_sphere.stl
    python3 -m tools.plot.immersed_stl bench OUT_DIR BENCH.json
    python3 -m tools.plot.immersed_stl walls OUT_DIR WALLS.json
    python3 -m tools.plot.immersed_stl fields OUT_DIR RUN [SWEEP_RUN ...]
"""
import json
import os
import struct
import sys

import numpy as np

CENTRE, RADIUS, LEVELS = (0.5, 0.5, 0.5), 0.2, 4


def icosphere(c, r, levels):
    """The triangles of an icosphere of radius `r` about `c`, outward (the
    same as `tests/surface.rs`)."""
    t = (1.0 + 5.0 ** 0.5) / 2.0
    v = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t],
         [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]]
    v = [list(np.array(p, float) / np.linalg.norm(p)) for p in v]
    f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9],
         [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2],
         [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10],
         [8, 6, 7], [9, 8, 1]]
    for _ in range(levels):
        mid = {}

        def m(a, b):
            k = (min(a, b), max(a, b))
            if k not in mid:
                p = (np.array(v[a]) + np.array(v[b])) / 2.0
                v.append(list(p / np.linalg.norm(p)))
                mid[k] = len(v) - 1
            return mid[k]

        g = []
        for a, b, cc in f:
            ab, bc, ca = m(a, b), m(b, cc), m(cc, a)
            g += [[a, ab, ca], [b, bc, ab], [cc, ca, bc], [ab, bc, ca]]
        f = g
    v = np.array(c) + r * np.array(v)
    return v[np.array(f)]


def write_stl(path, tris):
    """Binary STL: 80-byte header, count, then normal, corners and a zero
    attribute per triangle."""
    with open(path, "wb") as out:
        out.write(b"Foehncraft icosphere".ljust(80, b" "))
        out.write(struct.pack("<I", len(tris)))
        for t in tris:
            n = np.cross(t[1] - t[0], t[2] - t[0])
            n = n / np.linalg.norm(n)
            out.write(struct.pack("<12fH", *n, *t[0], *t[1], *t[2], 0))


METHODS = {
    "MC": ("#d1495b", "o", "marching cubes"),
    "exact": ("#edae49", "s", "exact crossings"),
    "closest": ("#66a182", "v", "+ closest point"),
    "DC": ("#00798c", "^", "+ dual contouring"),
    "DC+face": ("#30638e", "D", "+ face points"),
}
DASH = {"hex": "-", "tet": "--", "poly": ":"}
BODIES = {"sphere": "icosphere", "box": "turned box (sharp)", "torus": "torus"}


def bench(out, path):
    """Volume and wall-area errors against h, one column per body."""
    import matplotlib.pyplot as plt

    rows = json.load(open(path))
    fig, axs = plt.subplots(2, 3, figsize=(13.0, 7.6), constrained_layout=True)
    for j, body in enumerate(BODIES):
        for i, (key, what) in enumerate((("volume", "fluid volume"), ("area", "wall area"))):
            ax = axs[i, j]
            for mesh in ("hex", "tet", "poly"):
                for m, (col, mk, label) in METHODS.items():
                    r = sorted((x for x in rows if x["body"] == body and x["mesh"] == mesh
                                and x["method"] == m), key=lambda x: x["n"])
                    if not r:
                        continue
                    ax.loglog([2.0 / x["n"] for x in r], [x[key] for x in r], DASH[mesh],
                              color=col, marker=mk, ms=4, lw=1.2,
                              label=f"{label}" if mesh == "hex" else None)
            h = np.array([2 / 32, 2 / 8])
            ref = min(x[key] for x in rows if x["body"] == body and x["n"] == 32)
            ax.loglog(h, ref * (h / h[0]) ** 2, "k-", lw=0.8, label="slope 2" if j == 0 else None)
            ax.loglog(h, 3 * ref * (h / h[0]), "k:", lw=0.8, label="slope 1" if j == 0 else None)
            ax.set_title(f"{BODIES[body]}: {what}", fontsize=11)
            ax.set_xlabel("cell width h")
            if j == 0:
                ax.set_ylabel(f"relative error of the {what}")
            ax.grid(True, which="both", alpha=0.3)
    axs[0, 0].legend(fontsize=8, loc="lower right")
    from matplotlib.lines import Line2D
    axs[1, 0].legend([Line2D([], [], color="k", ls=d) for d in DASH.values()],
                     [f"{k} mesh" for k in DASH], fontsize=8, loc="lower right")
    fig.suptitle("Triangle surfaces cut out of box meshes of 8³ to 32³ cells, against the "
                 "exact moments of the same triangles", fontsize=12)
    path = os.path.join(out, "immersed_stl_benchmark.png")
    fig.savefig(path, dpi=100)
    print(path)


def walls(out, path):
    """The section z = z0 of the turned box cut out of the hexahedral mesh:
    the cut cells and the wall faces with the closest point against dual
    contouring with face points, the walls coloured by the distance of
    their face centroids to the box, the exact section dashed."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    from tools.plot.colormaps import RAINBOW

    d = json.load(open(path))
    exact = np.array(d["exact"])
    lo, hi = exact.reshape(-1, 2).min(0), exact.reshape(-1, 2).max(0)
    pad = 0.12
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 6.0), constrained_layout=True)
    norm = Normalize(0.0, 0.25)
    titles = {"closest": "closest point: the corners are cut off",
              "sharp": "dual contouring and face points: kept sharp"}
    area = 8 * (0.55 * 0.4 + 0.55 * 0.3 + 0.4 * 0.3)
    lc = None

    def draw(ax, name):
        ax.add_collection(LineCollection(np.array(d[name]["cells"]), colors="#7d8b99",
                                         linewidths=0.8))
        ax.add_collection(LineCollection(exact, colors="k", linewidths=1.0, linestyles="--"))
        walls = LineCollection(np.array(d[name]["walls"]), cmap=RAINBOW, norm=norm,
                               linewidths=2.6)
        walls.set_array(np.array(d[name]["distance"]))
        ax.add_collection(walls)
        ax.set_aspect("equal")
        return walls

    for ax, name in zip(axs, ("closest", "sharp")):
        lc = draw(ax, name)
        ax.set_xlim(lo[0] - pad, hi[0] + pad)
        ax.set_ylim(lo[1] - pad, hi[1] + pad)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_title(f"{titles[name]}\nwall area {100 * (d[name]['area'] / area - 1):+.1f} % "
                     "of the exact (whole box)", fontsize=11)
        # the lowest corner, enlarged
        k = np.argmin(exact.reshape(-1, 2)[:, 1])
        c = exact.reshape(-1, 2)[k]
        ins = ax.inset_axes([0.6, 0.04, 0.36, 0.36])
        draw(ins, name)
        ins.set_xlim(c[0] - 0.14, c[0] + 0.14)
        ins.set_ylim(c[1] - 0.06, c[1] + 0.22)
        ins.set_xticks([])
        ins.set_yticks([])
        ax.indicate_inset_zoom(ins, edgecolor="#444444")
    fig.colorbar(lc, ax=axs, shrink=0.8, pad=0.01, extend="max",
                 label="distance of the wall face centroid to the box (cells)")
    fig.suptitle(f"Turned box cut out of a {round(2 / d['h'])}³ hexahedral mesh: cells and "
                 f"wall faces in the plane z = {d['z0']:g} m (box dashed)", fontsize=12)
    path = os.path.join(out, "immersed_stl_box_walls.png")
    fig.savefig(path, dpi=100)
    print(path)

# VTK cell types with fixed faces: tetra, hexahedron, wedge, pyramid
FACES = {
    10: [[0, 1, 3], [1, 2, 3], [2, 0, 3], [0, 2, 1]],
    12: [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
    13: [[0, 1, 2], [3, 5, 4], [0, 3, 4, 1], [1, 4, 5, 2], [2, 5, 3, 0]],
    14: [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
}


def read_run(path):
    """The cells of a 3D run as rings of point ids per face, the points,
    the cell data and the wall faces of `sphere` with their owner cell."""
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
        w = f["VTKHDF"]["sphere"]
        wpts = w["Points"][:]
        woff = w["Polygons"]["Offsets"][:]
        wconn = w["Polygons"]["Connectivity"][:]
        owner = w["CellData"]["OwnerCell"][:]
    wall = [wpts[wconn[woff[i]:woff[i + 1]]] for i in range(len(woff) - 1)]

    def faces(c):
        if types[c] == 42:
            return [fconn[foff[k]:foff[k + 1]] for k in p2f[poff[c]:poff[c + 1]]]
        nodes = conn[off[c]:off[c + 1]]
        return [nodes[r] for r in FACES[int(types[c])]]

    return pts, faces, data, wall, owner


def section(pts, faces, n_cells, axis, x0):
    """The polygons of the cells in the plane x[axis] = x0, with their
    cells: the crossings of the face edges, chained around each cell."""
    lo = np.zeros(n_cells)
    polys, cells = [], []
    other = [k for k in range(3) if k != axis]
    for c in range(n_cells):
        segs = []
        for ring in faces(c):
            z = pts[ring, axis] - x0
            if z.min() >= 0 or z.max() < 0:
                continue
            xs = []
            for i in range(len(ring)):
                a, b = ring[i], ring[(i + 1) % len(ring)]
                da, db = pts[a, axis] - x0, pts[b, axis] - x0
                if (da < 0) != (db < 0):
                    t = da / (da - db)
                    xs.append(((min(a, b), max(a, b)),
                               pts[a] + t * (pts[b] - pts[a])))
            if len(xs) >= 2:
                d = xs[-1][1] - xs[0][1]
                xs.sort(key=lambda e: e[1] @ d)
                segs += [(xs[i], xs[i + 1]) for i in range(0, len(xs) - 1, 2)]
        if len(segs) < 3:
            continue
        nb, at = {}, {}
        for (ka, pa), (kb, pb) in segs:
            nb.setdefault(ka, []).append(kb)
            nb.setdefault(kb, []).append(ka)
            at[ka], at[kb] = pa, pb
        start = segs[0][0][0]
        loop, prev, cur = [start], None, start
        while True:
            nxt = [k for k in nb[cur] if k != prev]
            if not nxt or nxt[0] == start:
                break
            prev, cur = cur, nxt[0]
            loop.append(cur)
            if len(loop) > len(nb):
                break
        polys.append(np.array([at[k][other] for k in loop]))
        cells.append(c)
    del lo
    return polys, np.array(cells)


def drag(run):
    """The mesh, its cells per side and the drag coefficient
    K = F/(6π μ a U) with F = G L³ = 1 N of a run of the sphere array case."""
    import tomllib
    r = json.load(open(os.path.join(run, "result.json")))
    u = r["solves"][0]["periodic"][0]["mass_flow"]
    with open(os.path.join(run, "case.resolved.toml"), "rb") as f:
        g = tomllib.load(f)["mesh"]["generate"]
    return g["type"], g["cells"][0], 1.0 / (6.0 * np.pi * RADIUS * u)


def fields(out, run, sweep=()):
    """Stokes flow through the sphere array: the pressure on the cut wall
    faces, the speed in the plane y = ½ and the drag coefficient against
    the cell width."""
    import matplotlib.pyplot as plt
    from matplotlib import colormaps
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    from matplotlib.patches import Circle
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    from tools.plot.colormaps import RAINBOW, UNIFORM

    pts, faces, d, wall, owner = read_run(os.path.join(run, "fields.vtkhdf"))
    n_cells = len(d["pressure"])
    fig = plt.figure(figsize=(15.0, 5.0), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.0])

    # the cut wall, coloured by the pressure, lit
    ax = fig.add_subplot(gs[0], projection="3d")
    p = d["pressure"][owner]
    m = np.abs(p).max()
    norm = Normalize(-m, m)
    col = colormaps[UNIFORM](norm(p))
    light = np.array([-0.4, -0.8, 0.6])
    light /= np.linalg.norm(light)
    for i, f in enumerate(wall):
        n = sum(np.cross(f[a], f[(a + 1) % len(f)]) for a in range(len(f)))
        n = n / max(np.linalg.norm(n), 1e-300)
        col[i, :3] *= 0.55 + 0.45 * max(n @ light, 0.0)
    ax.add_collection3d(Poly3DCollection(wall, facecolors=col, edgecolors=(0, 0, 0, 0.25),
                                         linewidths=0.15))
    c, r = np.array(CENTRE), RADIUS * 1.05
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-62)
    ax.set_axis_off()
    ax.set_title(f"pressure on the {len(wall)} cut wall faces\n(flow along +x)", fontsize=11)
    fig.colorbar(ScalarMappable(norm, colormaps[UNIFORM]), ax=ax, shrink=0.7, pad=0.0,
                 label="pressure p − p̄ (Pa)")

    # the speed in the plane y = ½, the cut cells outlined
    ax = fig.add_subplot(gs[1])
    polys, cells = section(pts, faces, n_cells, 1, CENTRE[1] + 1e-9)
    speed = d["speed"][cells]
    pc = PolyCollection(polys, array=speed, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    near = [q for q in polys
            if np.min(np.hypot(q[:, 0] - CENTRE[0], q[:, 1] - CENTRE[2])) < RADIUS + 0.06]
    ax.add_collection(PolyCollection(near, facecolors="none", edgecolors=(0, 0, 0, 0.5),
                                     linewidths=0.4))
    ax.add_patch(Circle((CENTRE[0], CENTRE[2]), RADIUS, fill=False, ec="w", lw=1.0, ls="--"))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(CENTRE[2], 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title("speed in the plane y = 0.5 m (upper half), cells near\n"
                 "the sphere outlined, the exact sphere dashed", fontsize=11)
    fig.colorbar(pc, ax=ax, shrink=0.9, pad=0.01, label="speed (m/s)")

    # the drag coefficient against the cell width
    ax = fig.add_subplot(gs[2])
    k_ref = 2.152
    runs = [drag(x) for x in sweep] if sweep else [drag(run)]
    style = {"hex": ("#30638e", "o"), "tet": ("#d1495b", "s"), "poly": ("#66a182", "^")}

    def series(h, e, colr, mk, ls, label):
        e = np.array(e)
        ax.loglog(h, np.abs(e), ls, color=colr, marker=mk, ms=7, mew=1.5, label=label)
        for hi, ei in zip(h, e):
            if ei > 0:
                ax.loglog([hi], [ei], mk, color=colr, ms=7, mfc="white", mew=1.5)

    for mesh, (colr, mk) in style.items():
        rs = sorted((x for x in runs if x[0] == mesh), key=lambda x: x[1])
        if rs:
            series([1.0 / x[1] for x in rs], [100 * (x[2] / k_ref - 1) for x in rs], colr, mk,
                   "-", f"cut, STL sphere, {mesh}")
    pen = [(16, 6.2), (24, 5.5), (32, 3.9)]
    series([1 / n for n, _ in pen], [e for _, e in pen], "#edae49", "D", "--",
           "penalised, hex (gallery 291)")
    h = np.array([1 / 48, 1 / 16])
    ax.loglog(h, 0.9 * (h / h[0]) ** 2, "k-", lw=0.8, label="slope 2")
    ax.loglog(h, 3.9 * (h / (1 / 32)), "k:", lw=0.8, label="slope 1")
    ax.set_xlabel("cell width h = L / n (m)")
    ax.set_ylabel("|K / 2.152 − 1| (%): filled below, open above")
    ax.set_title("drag coefficient against the series\n(Hasimoto 1959, Sangani & Acrivos 1982)",
                 fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    path = os.path.join(out, "immersed_stl_sphere_array.png")
    fig.savefig(path, dpi=100)
    print(path)


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__)
    cmd = argv[1]
    if cmd == "stl":
        tris = icosphere(CENTRE, RADIUS, LEVELS)
        write_stl(argv[2], tris)
        print(f"{argv[2]}: {len(tris)} triangles")
        return
    import matplotlib

    matplotlib.use("Agg")
    if cmd == "bench":
        bench(argv[2], argv[3])
    elif cmd == "walls":
        walls(argv[2], argv[3])
    elif cmd == "fields":
        fields(argv[2], argv[3], argv[4:])
    else:
        sys.exit(f"unknown command `{cmd}`: stl, bench, walls or fields")


if __name__ == "__main__":
    main(sys.argv)
