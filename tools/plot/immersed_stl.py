"""STL bodies cut out of a mesh, L13 stage 2 (gallery and recipe 320).

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
    python3 -m tools.plot.immersed_stl fields OUT_DIR RUN
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
        fields(argv[2], argv[3])
    else:
        sys.exit(f"unknown command `{cmd}`: stl, bench, walls or fields")


if __name__ == "__main__":
    main(sys.argv)
