"""Field images of an axisymmetric flow in its meridional plane: the cells
of a mesh of revolution (`[mesh.revolved]`) averaged around the axis, drawn
in the (r, z) half-plane and mirrored to -r. Gallery 103 (D186) was drawn
with it.

    pip install matplotlib h5py numpy [pillow]
    python3 -m tools.plot.meridional_fields OUT.png FIELDS [options]

FIELDS is a `fields.vtkhdf` written by `foehn flow` with `output.vtkhdf =
true`. The cells of one meridional cell (one ring of sectors) are found by
the (r, z) of their centroid; each ring is drawn as the polygon of its first
cell's nodes in (r, z), coloured by the mean over the ring. One panel per
quantity:

    swirl       u_theta, the circumferential velocity
    meridional  |(u_r, u_z)|, with the streamlines of (u_r, u_z)
    axial       u_z
    radial      u_r
    pressure    p, from its minimum
    speed       |u|
    NAME        any cell dataset of the file

Options: --quantities swirl,meridional (default), --interface r0,z0,r1,z1;...
(draws the meridian of the sliding interface, pieces separated by `;`),
--axis-z z0,z1 (draws the axis), --title, --labels (panel titles, separated by `|`), --limits lo:hi,... (colour limits per panel), --mesh (cell
edges), --no-mirror, --dpi, --palette (a 256-colour PNG, needs pillow).
"""
import argparse

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.tri import LinearTriInterpolator, Triangulation  # noqa: E402
from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.vtkhdf import internal  # noqa: E402

LABELS = {
    "swirl": ("u_θ (m/s)", RAINBOW),
    "meridional": ("|(u_r, u_z)| (m/s)", RAINBOW),
    "axial": ("u_z (m/s)", RAINBOW),
    "radial": ("u_r (m/s)", "RdBu_r"),
    "pressure": ("p − p_min (Pa)", "cividis"),
    "speed": ("|u| (m/s)", RAINBOW),
}
SIGNED = ("radial",)


def read(path):
    """Rings of cells: (r, z) polygons, the cell data averaged over each ring,
    and the components of the velocity in (r, θ, z)."""
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    n = len(off) - 1
    cen = np.empty((n, 3))
    polys = []
    for c in range(n):
        ids = np.unique(conn[off[c]:off[c + 1]])
        q = pts[ids]
        cen[c] = q.mean(axis=0)
        polys.append(np.column_stack([np.hypot(q[:, 0], q[:, 1]), q[:, 2]]))
    r = np.hypot(cen[:, 0], cen[:, 1])
    th = np.arctan2(cen[:, 1], cen[:, 0])
    if "velocity_x" in data:
        u, v, w = data["velocity_x"], data["velocity_y"], data["velocity_z"]
        data["radial"] = u * np.cos(th) + v * np.sin(th)
        data["swirl"] = -u * np.sin(th) + v * np.cos(th)
        data["axial"] = w
        data["meridional"] = np.hypot(data["radial"], w)
        data["speed"] = np.sqrt(u * u + v * v + w * w)
    # the rings: cells with one (r, z) of the centroid
    scale = max(r.max(), np.ptp(cen[:, 2]))
    key = np.round(np.column_stack([r, cen[:, 2]]) / (1e-7 * scale)).astype(np.int64)
    _, ring, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    count = np.bincount(inv)
    mean = {k: np.bincount(inv, weights=v) / count for k, v in data.items()
            if v.ndim == 1 and len(v) == n}
    shapes = []
    for c in ring:
        q = polys[c]
        # the polygon in (r, z): the hull of the node positions
        q = np.unique(np.round(q / (1e-9 * scale)) * (1e-9 * scale), axis=0)
        m = q.mean(axis=0)
        q = q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))]
        shapes.append(hull(q))
    rz = np.column_stack([r[ring], cen[ring, 2]])
    return shapes, rz, mean


def cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def hull(q):
    """The convex hull of points in the plane (monotone chain)."""
    p = sorted(map(tuple, q))
    if len(p) <= 2:
        return np.array(p)

    def half(seq):
        out = []
        for x in seq:
            while len(out) >= 2 and cross2(np.subtract(out[-1], out[-2]),
                                           np.subtract(x, out[-2])) <= 0:
                out.pop()
            out.append(x)
        return out

    lo, hi = half(p), half(reversed(p))
    return np.array(lo[:-1] + hi[:-1])


def value(name, mean):
    if name == "pressure":
        return mean["pressure"] - mean["pressure"].min()
    if name in mean:
        return mean[name]
    raise SystemExit(f"unknown quantity {name} (datasets: {', '.join(sorted(mean))})")


def streamlines(ax, rz, mean, window, mirror):
    """Streamlines of the ring-mean (u_r, u_z) on a grid (masked across solids)."""
    x, y = rz[:, 0], rz[:, 1]
    ur, uz = mean["radial"], mean["axial"]
    if mirror:
        x, y = np.concatenate([x, -x]), np.concatenate([y, y])
        ur, uz = np.concatenate([ur, -ur]), np.concatenate([uz, uz])
    tri = Triangulation(x, y)
    tx, ty = x[tri.triangles], y[tri.triangles]
    edge = np.max(np.hypot(tx - np.roll(tx, 1, axis=1), ty - np.roll(ty, 1, axis=1)), axis=1)
    tri.set_mask(edge > 3.0 * np.median(edge))
    x0, x1, y0, y1 = window
    gx, gy = np.meshgrid(np.linspace(x0, x1, 300), np.linspace(y0, y1, 300))
    gu = LinearTriInterpolator(tri, ur)(gx, gy)
    gv = LinearTriInterpolator(tri, uz)(gx, gy)
    ax.streamplot(gx, gy, gu, gv, density=1.4, color=(1, 1, 1, 0.6), linewidth=0.6,
                  arrowsize=0.6, broken_streamlines=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("fields")
    ap.add_argument("--quantities", default="swirl,meridional")
    ap.add_argument("--interface", default=None,
                    help="meridian of the interface: r0,z0,r1,z1;r0,z0,r1,z1;...")
    ap.add_argument("--axis-z", default=None, help="z range of the axis line: z0,z1")
    ap.add_argument("--title", default="")
    ap.add_argument("--labels", default=None)
    ap.add_argument("--limits", default=None)
    ap.add_argument("--mesh", action="store_true")
    ap.add_argument("--no-mirror", action="store_true")
    ap.add_argument("--dpi", type=int, default=170)
    ap.add_argument("--palette", action="store_true")
    a = ap.parse_args()

    shapes, rz, mean = read(a.fields)
    names = a.quantities.split(",")
    mirror = not a.no_mirror
    allp = np.concatenate(shapes)
    rmax = allp[:, 0].max()
    window = [-rmax if mirror else allp[:, 0].min(), rmax, allp[:, 1].min(), allp[:, 1].max()]
    titles = a.labels.split("|") if a.labels else names
    limits = a.limits.split(",") if a.limits else [""] * len(names)
    w = (window[1] - window[0]) / (window[3] - window[2])
    fig, axes = plt.subplots(1, len(names), figsize=(4.2 * w * len(names) + 1.2, 4.6),
                             squeeze=False)
    for ax, name, title, lim in zip(axes[0], names, titles, limits):
        val = value(name, mean)
        label, cmap = LABELS.get(name, (name, RAINBOW))
        polys = list(shapes)
        vals = list(val)
        if mirror:
            polys += [np.column_stack([-p[:, 0], p[:, 1]]) for p in shapes]
            vals += vals
        vals = np.array(vals)
        lo, hi = (None, None)
        if lim.strip(":"):
            lo, hi = (float(x) for x in lim.split(":"))
        elif name in SIGNED:
            m = np.abs(vals).max()
            lo, hi = -m, m
        pc = PolyCollection(polys, array=vals, cmap=cmap,
                            edgecolors=(0, 0, 0, 0.25) if a.mesh else "face",
                            linewidths=0.2 if a.mesh else 0.3)
        pc.set_clim(lo, hi)
        ax.add_collection(pc)
        if name == "meridional":
            streamlines(ax, rz, mean, window, mirror)
        if a.interface:
            for seg in a.interface.split(";"):
                r0, z0, r1, z1 = (float(x) for x in seg.split(","))
                ax.plot([r0, r1], [z0, z1], color="white", lw=1.4, ls="--")
                if mirror:
                    ax.plot([-r0, -r1], [z0, z1], color="white", lw=1.4, ls="--")
        if a.axis_z:
            z0, z1 = (float(x) for x in a.axis_z.split(","))
            ax.plot([0, 0], [z0, z1], color=(1, 1, 1, 0.5), lw=0.6, ls="-.")
        ax.set_xlim(window[0], window[1])
        ax.set_ylim(window[2], window[3])
        ax.set_aspect("equal")
        ax.set_facecolor("#3a3a3a")
        ax.set_xlabel("r (m)")
        ax.set_ylabel("z (m)")
        ax.set_title(title, fontsize=10)
        fig.colorbar(pc, ax=ax, shrink=0.85, label=label)
    if a.title:
        fig.suptitle(a.title, fontsize=11)
    fig.tight_layout()
    if a.palette:
        import io

        from PIL import Image

        buf = io.BytesIO()
        fig.savefig(buf, dpi=a.dpi, format="png")
        buf.seek(0)
        Image.open(buf).convert("RGB").quantize(colors=256, dither=1).save(a.out, optimize=True)
    else:
        fig.savefig(a.out, dpi=a.dpi)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
