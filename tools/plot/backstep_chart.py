"""Backward-facing step of Driver & Seegmiller (1985), gallery 74: the
streamwise velocity with streamlines and the turbulent kinetic energy of
several turbulence models on the same mesh, one row per model, with the
reattachment point from the wall shear and the experiment's 6.26 H.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.backstep_chart OUT_PREFIX LABEL=DIR [LABEL=DIR ...]

DIR is the output directory of a run of backstep_zeta_f0.toml (or of the
same case with another model) with `output.vtkhdf = true`; the script
reads `fields.vtkhdf` and `wall_ymin.csv`. It writes OUT_PREFIX_velocity.png
and OUT_PREFIX_k.png and prints x_r/H per model, e.g.

    python3 -m tools.plot.backstep_chart docs/figures/gallery/backstep \\
        "k-ε (scalable)=output/keps" "k-ε-ζ-f₀=output/zf0" "k-ω SST=output/sst"
"""
import csv
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.vtkhdf import internal  # noqa: E402

X_EXP = 6.26           # Driver & Seegmiller (1985), x_r/H
WINDOW = (-3.0, 14.0, 0.0, 3.0)
U_REF = 1.0


def fields(d):
    """Front-face quads (2D mesh, one hexahedron thick) and cell data."""
    with h5py.File(os.path.join(d, "fields.vtkhdf"), "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        assert np.all(np.diff(off) == 8), "hexahedra only"
        cells = conn.reshape(-1, 8)
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    xyz = pts[cells]                       # (n, 8, 3)
    zmin = xyz[:, :, 2].min(axis=1, keepdims=True)
    # the four nodes on the lower z face, ordered by angle about their centre
    quads = []
    for c, zm in zip(xyz, zmin):
        q = c[np.abs(c[:, 2] - zm[0]) < 1e-12][:, :2]
        m = q.mean(axis=0)
        q = q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))]
        quads.append(q)
    quads = np.array(quads)
    return quads, quads.mean(axis=1), data


def reattachment(d):
    with open(os.path.join(d, "wall_ymin.csv")) as f:
        rows = [r for r in f if not r.startswith("#")]
    w = sorted((float(r["x"]), float(r["tau_x"])) for r in csv.DictReader(rows))
    xr = None
    for (x0, t0), (x1, t1) in zip(w, w[1:]):
        if x0 > 0 and t0 < 0 <= t1:
            xr = x0 + (x1 - x0) * (-t0) / (t1 - t0)
    return xr


def solid(x, y):
    return (x < 0) & (y < 1)


def grid_velocity(cen, data):
    """u, v on a uniform grid (linear on the triangulated cell centres)."""
    tri = mtri.Triangulation(cen[:, 0], cen[:, 1])
    tc = cen[tri.triangles].mean(axis=1)
    # no triangle across the step corner or with very long edges
    e = cen[tri.triangles]
    edge = np.max(np.linalg.norm(e - np.roll(e, 1, axis=1), axis=2), axis=1)
    tri.set_mask(solid(tc[:, 0], tc[:, 1]) | (edge > 2.0))
    gx = np.linspace(WINDOW[0], WINDOW[1], 700)
    gy = np.linspace(WINDOW[2], WINDOW[3], 160)
    X, Y = np.meshgrid(gx, gy)
    u = mtri.LinearTriInterpolator(tri, data["velocity_x"])(X, Y)
    v = mtri.LinearTriInterpolator(tri, data["velocity_y"])(X, Y)
    m = solid(X, Y) | np.ma.getmaskarray(u)
    return gx, gy, np.ma.array(u.filled(0), mask=m), np.ma.array(v.filled(0), mask=m)


def panel(ax, quads, values, cmap, vmin, vmax):
    keep = (quads[:, :, 0].max(axis=1) > WINDOW[0]) & (quads[:, :, 0].min(axis=1) < WINDOW[1]) \
        & (quads[:, :, 1].min(axis=1) < WINDOW[3])
    pc = PolyCollection(quads[keep], array=values[keep], cmap=cmap, edgecolors="face",
                        linewidths=0.2)
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    ax.fill([WINDOW[0], 0, 0, WINDOW[0]], [0, 0, 1, 1], color="0.55", zorder=3)
    ax.set_xlim(WINDOW[0], WINDOW[1])
    ax.set_ylim(WINDOW[2], WINDOW[3])
    ax.set_aspect("equal")
    ax.set_ylabel("y/H")
    return pc


def marks(ax, xr):
    ax.plot([X_EXP], [0.0], marker="^", color="k", ms=7, clip_on=False, zorder=5)
    if xr:
        ax.plot([xr], [0.0], marker="^", color="#d62728", ms=7, clip_on=False, zorder=5)


def main():
    prefix, runs = sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]]
    loaded = []
    for label, d in runs:
        quads, cen, data = fields(d)
        xr = reattachment(d)
        loaded.append((label, quads, cen, data, xr))
        if xr:
            print(f"{label:22s} x_r/H {xr:.2f} ({100 * (xr / X_EXP - 1):+.1f} % from {X_EXP})")
        else:
            print(f"{label:22s} no reattachment")
    n = len(loaded)

    # velocity with streamlines
    fig, axes = plt.subplots(n, 1, figsize=(11, 2.05 * n + 0.5), sharex=True,
                             constrained_layout=True)
    for ax, (label, quads, cen, data, xr) in zip(np.atleast_1d(axes), loaded):
        pc = panel(ax, quads, data["velocity_x"] / U_REF, "RdBu_r", -0.3, 1.1)
        gx, gy, u, v = grid_velocity(cen, data)
        seeds = [(WINDOW[0] + 0.05, y) for y in np.linspace(1.08, 2.9, 9)]
        xb = 0.5 * (xr or X_EXP)
        seeds += [(xb, y) for y in np.linspace(0.06, 0.75, 6)] + [(0.25, 0.08), (0.15, 0.25)]
        ax.streamplot(gx, gy, u, v, start_points=seeds, density=8, color="k",
                      linewidth=0.5, arrowstyle="-", broken_streamlines=False,
                      integration_direction="both", maxlength=6.0)
        ax.contour(gx, gy, u, levels=[0.0], colors="w", linewidths=0.9)
        marks(ax, xr)
        xt = f"x_r/H = {xr:.2f}" if xr else "no reattachment"
        ax.set_title(f"{label}:  {xt}  (experiment {X_EXP})", loc="left", fontsize=10)
    np.atleast_1d(axes)[-1].set_xlabel("x/H")
    fig.colorbar(pc, ax=axes, shrink=0.8, pad=0.01, label="U/U_ref")
    fig.savefig(prefix + "_velocity.png", dpi=170)
    plt.close(fig)

    # turbulent kinetic energy
    kmax = max(np.percentile(d["k"], 99.9) for _, _, _, d, _ in loaded) / U_REF**2
    fig, axes = plt.subplots(n, 1, figsize=(11, 2.05 * n + 0.5), sharex=True,
                             constrained_layout=True)
    for ax, (label, quads, cen, data, xr) in zip(np.atleast_1d(axes), loaded):
        pc = panel(ax, quads, data["k"] / U_REF**2, RAINBOW, 0.0, kmax)
        gx, gy, u, v = grid_velocity(cen, data)
        ax.contour(gx, gy, u, levels=[0.0], colors="c", linewidths=0.9)
        marks(ax, xr)
        ax.set_title(f"{label}", loc="left", fontsize=10)
    np.atleast_1d(axes)[-1].set_xlabel("x/H")
    fig.colorbar(pc, ax=axes, shrink=0.8, pad=0.01, label="k/U_ref²")
    fig.savefig(prefix + "_k.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
