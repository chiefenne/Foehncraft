"""Radiation R1a (gallery entry of the PR): the box furnace with seven
cooled tubes under the roof (D255), on a prescribed flame. Temperature and
radiative source on the mid-depth plane, the incident flux on the side wall
and the net flux into each tube; and the error of Sobol against
pseudo-random emission (verification case 6).

    pip install gmsh matplotlib numpy
    python3 examples/meshes/furnace.py --size 0.15 --out output/furnace.msh
    cargo run --release -p foehncraft-radiation --example furnace -- \\
        output/furnace.msh output/furnace 4000000
    cargo run --release -p foehncraft-radiation --example sobol_study > output/sobol.csv
    python3 -m tools.plot.radiation_furnace docs/figures/gallery output/furnace output/sobol.csv

The figures are radiation_furnace_plane.png, radiation_furnace_walls.png
and radiation_sobol.png.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle

from tools.plot.colormaps import RAINBOW

LENGTH, DEPTH, HEIGHT = 4.0, 2.0, 2.5
TUBE_R, TUBE_Z = 0.08, 2.1
TUBE_X = [0.9 + 0.4 * i for i in range(7)]
BURNER = (1.0, 0.9, 0.25)
OUTLET = (0.6, 1.4, 1.75, 2.35)


def read(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return rows


def column(rows, key, cast=float):
    return np.array([cast(r[key]) for r in rows])


def plane_triangulation(x, z):
    """Delaunay triangles of the points (x, z) without those across a tube
    or longer than the cells."""
    tri = mtri.Triangulation(x, z)
    t = tri.triangles
    cx, cz = x[t].mean(axis=1), z[t].mean(axis=1)
    edge = np.max(
        [np.hypot(x[t[:, i]] - x[t[:, (i + 1) % 3]], z[t[:, i]] - z[t[:, (i + 1) % 3]]) for i in range(3)],
        axis=0,
    )
    mask = edge > 0.35
    for tx in TUBE_X:
        mask |= np.hypot(cx - tx, cz - TUBE_Z) < TUBE_R * 1.05
    tri.set_mask(mask)
    return tri


def tubes(ax):
    for tx in TUBE_X:
        ax.add_patch(Circle((tx, TUBE_Z), TUBE_R, facecolor="0.85", edgecolor="0.25", lw=0.8, zorder=5))


def frame(ax, title):
    ax.set_xlim(0, LENGTH)
    ax.set_ylim(0, HEIGHT)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title(title, fontsize=11)


def plane(out, cells):
    """Temperature and radiative source on the plane y = 1 m (cells within
    a band about it)."""
    x, y, z = (column(cells, k) for k in "xyz")
    band = np.abs(y - DEPTH / 2) < 0.12
    tri = plane_triangulation(x[band], z[band])
    t = column(cells, "temperature")[band] if "temperature" in cells[0] else None
    src = column(cells, "source")[band] * 1e-3
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 9.6), constrained_layout=True)
    ax = axes[0]
    # the band of cell centres stops half a cell short of the walls, where
    # the gas is at its lowest temperature
    ax.set_facecolor(plt.get_cmap(RAINBOW)(0.0))
    cs = ax.tricontourf(tri, t, levels=48, cmap=RAINBOW)
    fig.colorbar(cs, ax=ax, label="temperature (K)", shrink=0.9, ticks=np.arange(1300, 2000, 100))
    tubes(ax)
    _yb, zb, rb = BURNER
    ax.plot([0, 0], [zb - rb, zb + rb], color="k", lw=4, solid_capstyle="butt")
    ax.annotate("burner", (0.03, zb), xytext=(0.25, 0.25), fontsize=9,
                arrowprops={"arrowstyle": "-", "lw": 0.6})
    ax.plot([LENGTH, LENGTH], [OUTLET[2], OUTLET[3]], color="w", lw=4, solid_capstyle="butt")
    ax.annotate("outlet", (LENGTH - 0.03, 0.5 * (OUTLET[2] + OUTLET[3])), xytext=(3.25, 1.45),
                fontsize=9, arrowprops={"arrowstyle": "-", "lw": 0.6})
    frame(ax, "Prescribed flame: temperature on the mid-depth plane")
    ax = axes[1]
    lim = np.percentile(np.abs(src), 99.5)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=max(lim * 0.25, src.max(initial=0.0)))
    cs = ax.tricontourf(tri, np.clip(src, -lim, lim), levels=48, cmap="RdBu_r", norm=norm)
    fig.colorbar(cs, ax=ax, label="radiative source, absorbed − emitted (kW/m³)", shrink=0.9)
    tubes(ax)
    frame(ax, "Monte Carlo radiation: net source of the gas")
    path = os.path.join(out, "radiation_furnace_plane.png")
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(path)


def walls(out, faces):
    """Incident flux on the side wall y = 0 and the net flux into each tube."""
    name = np.array([r["boundary"] for r in faces])
    x, y, z = (column(faces, k) for k in "xyz")
    inc = column(faces, "incident") * 1e-3
    net = column(faces, "net") * 1e-3
    area = column(faces, "area")
    side = (name == "walls") & (y < 1e-6)
    fig, axes = plt.subplots(
        2, 1, figsize=(8.2, 8.0), constrained_layout=True, gridspec_kw={"height_ratios": [1.35, 1]}
    )
    ax = axes[0]
    tri = plane_triangulation(x[side], z[side])
    cs = ax.tricontourf(tri, inc[side], levels=48, cmap=RAINBOW)
    fig.colorbar(cs, ax=ax, label="incident radiative flux (kW/m²)", shrink=0.9)
    tubes(ax)
    frame(ax, "Incident flux on the side wall y = 0")
    ax = axes[1]
    on = name == "tubes"
    tube = np.array([np.argmin([abs(xi - tx) for tx in TUBE_X]) for xi in x[on]])
    q = np.array([np.sum(net[on][tube == i] * area[on][tube == i]) / np.sum(area[on][tube == i]) for i in range(7)])
    # around the tube: net flux against the angle from the top (towards the roof)
    angle = np.degrees(np.arctan2(x[on] - np.array(TUBE_X)[tube], z[on] - TUBE_Z)) % 360.0
    for i in (0, 3, 6):
        sel = tube == i
        order = np.argsort(angle[sel])
        a, v = angle[sel][order], net[on][sel][order]
        bins = np.linspace(0, 360, 25)
        k = np.digitize(a, bins) - 1
        mean = np.array([v[k == j].mean() if np.any(k == j) else np.nan for j in range(24)])
        ax.plot(0.5 * (bins[1:] + bins[:-1]), mean, marker="o", ms=3, lw=1.2,
                label=f"tube {i + 1} (x = {TUBE_X[i]:.1f} m), mean {q[i]:.0f} kW/m²")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 45))
    ax.set_xlabel("angle around the tube from the top, towards the outlet (degrees)")
    ax.set_ylabel("net flux into the tube (kW/m²)")
    ax.set_title("Cooled tubes: the side facing the flame takes most", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    path = os.path.join(out, "radiation_furnace_walls.png")
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(path)
    for i in range(7):
        print(f"  tube {i + 1}: mean net flux {q[i]:.1f} kW/m²")


def sobol(out, path):
    rows = read(path)
    n = column(rows, "bundles")
    fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
    styles = {"slab": ("C0", "slab wall flux (case 1, τ = 1)"), "plates": ("C3", "plate exchange (case 3)")}
    for case, (c, label) in styles.items():
        s = column(rows, f"{case}_sobol")
        r = column(rows, f"{case}_random")
        ax.loglog(n, r, "--o", color=c, ms=4, label=f"{label}, pseudo-random")
        ax.loglog(n, s, "-s", color=c, ms=4, label=f"{label}, scrambled Sobol")
    ref = column(rows, "slab_random")[0] * (n / n[0]) ** -0.5
    ax.loglog(n, ref, ":", color="0.5", lw=1, label=r"$N^{-1/2}$")
    ax.set_xlabel("bundles N")
    ax.set_ylabel("root-mean-square relative error")
    ax.set_title("Emission points: scrambled Sobol against pseudo-random", fontsize=11)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    p = os.path.join(out, "radiation_sobol.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    print(p)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    cells = read(os.path.join(run, "cells.csv"))
    plane(out, cells)
    walls(out, read(os.path.join(run, "walls.csv")))
    if len(sys.argv) > 3:
        sobol(out, sys.argv[3])


if __name__ == "__main__":
    main()
