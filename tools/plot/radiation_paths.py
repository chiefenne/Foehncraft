"""Radiation R1a (gallery entry of the PR): the box furnace of D255 with the
paths of chosen bundles, traced by the same call that gives the wall
fluxes (``Domain::trace``). Every point of a path is where the call's own
transport left a cell, reflected at a wall or ended the bundle; the line
fades with the energy the bundle still carries, as the grey gas absorbs it.

    pip install gmsh matplotlib numpy
    python3 examples/meshes/furnace.py --size 0.15 --out output/furnace.msh
    cargo run --release -p foehncraft-radiation --example furnace -- \\
        output/furnace.msh output/furnace 16000000 60
    python3 -m tools.plot.radiation_paths docs/figures/gallery output/furnace

The figures are radiation_paths_3d.png (the paths in the chamber over the
net flux into the far walls), radiation_paths_plane.png (the paths seen
from the side over the gas temperature) and radiation_paths_walls.png (the
net flux into all six walls).
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects
from matplotlib.colors import Normalize
from matplotlib.patches import Circle, Rectangle
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from tools.plot.colormaps import RAINBOW, SECOND
from tools.plot.radiation_furnace import (
    BURNER,
    DEPTH,
    HEIGHT,
    LENGTH,
    OUTLET,
    TUBE_R,
    TUBE_X,
    TUBE_Z,
    WALLS,
    column,
    read,
    smooth_grid,
)

#: Colours of the paths: bundles emitted by the gas and by a wall.
GAS, WALL = "#2ee6ff", "#9dff3a"
#: The faintest a path is drawn, whatever energy it still carries.
FAINTEST = 0.1


def paths(rows):
    """The traced paths as (source, points (n, 3), energy fraction (n,), kinds)."""
    out = {}
    for r in rows:
        p = out.setdefault(int(r["path"]), {"source": r["source"], "x": [], "f": [], "at": []})
        p["x"].append([float(r["x"]), float(r["y"]), float(r["z"])])
        p["f"].append(float(r["fraction"]))
        p["at"].append(r["at"])
    return [(p["source"], np.array(p["x"]), np.array(p["f"]), p["at"]) for _, p in sorted(out.items())]


def segments(traced, dims, faintest=FAINTEST):
    """Line segments of all paths in the coordinates ``dims`` with a colour
    whose opacity is the mean energy fraction along the segment."""
    seg, rgba = [], []
    for source, x, f, _at in traced:
        c = matplotlib.colors.to_rgb(GAS if source == "gas" else WALL)
        for i in range(len(x) - 1):
            seg.append(x[i : i + 2][:, dims])
            rgba.append((*c, max(faintest, 0.5 * (f[i] + f[i + 1]))))
    return np.array(seg), np.array(rgba)


def faces(rows):
    name = np.array([r["boundary"] for r in rows])
    x, y, z = (column(rows, k) for k in "xyz")
    return name, x, y, z, column(rows, "net") * 1e-3, column(rows, "area")


def walls_of(f):
    """The six walls of the box as (title, mask, u, v, extent of u and v,
    labels): the burner belongs to the front wall, the opening is left out."""
    name, x, y, z, _net, _area = f
    eps = 1e-6
    w = np.isin(name, WALLS)
    front = (name == "front") | (name == "burner")
    return {
        "floor": ("floor z = 0", w & (z < eps), x, y, (LENGTH, DEPTH), ("x (m)", "y (m)")),
        "roof": ("roof z = 2.5 m", w & (z > HEIGHT - eps), x, y, (LENGTH, DEPTH), ("x (m)", "y (m)")),
        "front": ("front x = 0, burner", front, y, z, (DEPTH, HEIGHT), ("y (m)", "z (m)")),
        "side0": ("side wall y = 0", w & (y < eps), x, z, (LENGTH, HEIGHT), ("x (m)", "z (m)")),
        "side2": ("side wall y = 2 m", w & (y > DEPTH - eps), x, z, (LENGTH, HEIGHT), ("x (m)", "z (m)")),
        "back": ("back x = 4 m, opening", w & (x > LENGTH - eps), y, z, (DEPTH, HEIGHT), ("y (m)", "z (m)")),
    }


def on_grid(u, v, val, area, extent):
    """The net flux averaged over patches of about 0.1 m of the wall
    (area-weighted Gaussian, see radiation_furnace.smooth_grid)."""
    return smooth_grid(u, v, val, area, extent)


def flux_scale(f):
    """One colour scale for all walls, from the averaged flux that is drawn
    (the single faces scatter far wider than the average)."""
    _name, _x, _y, _z, net, area = f
    g = np.concatenate([
        on_grid(u[m], v[m], net[m], area[m], extent)[2].ravel()
        for _t, m, u, v, extent, _lab in walls_of(f).values()
    ])
    g = g[np.isfinite(g)]
    return Normalize(vmin=min(0.0, np.percentile(g, 0.5)), vmax=np.percentile(g, 99.9))


def cylinder(ax, x0, z0, r, color):
    a, y = np.meshgrid(np.linspace(0, 2 * np.pi, 49), np.linspace(0, DEPTH, 2))
    ax.plot_surface(x0 + r * np.sin(a), y, z0 + r * np.cos(a), color=color, shade=True, linewidth=0,
                    edgecolor="none", antialiased=False, rstride=1, cstride=1)


def box_edges(ax):
    L, D, H = LENGTH, DEPTH, HEIGHT
    for a, b in [
        ((0, 0, 0), (L, 0, 0)), ((0, D, 0), (L, D, 0)), ((0, 0, H), (L, 0, H)), ((0, D, H), (L, D, H)),
        ((0, 0, 0), (0, D, 0)), ((L, 0, 0), (L, D, 0)), ((0, 0, H), (0, D, H)), ((L, 0, H), (L, D, H)),
        ((0, 0, 0), (0, 0, H)), ((L, 0, 0), (L, 0, H)), ((0, D, 0), (0, D, H)), ((L, D, 0), (L, D, H)),
    ]:
        ax.plot(*zip(a, b), color="0.35", lw=0.7)


def chamber(out, f, traced, bundles):
    """The paths in the chamber, seen from above the near side wall, over
    the net flux into the three far walls (floor, front wall with the
    burner, side wall y = 2 m)."""
    norm = flux_scale(f)
    cmap = plt.get_cmap(SECOND)
    walls = walls_of(f)
    _name, _x, _y, _z, net, area = f
    fig = plt.figure(figsize=(10.0, 6.4))
    ax = fig.add_axes([-0.03, -0.05, 0.9, 0.98], projection="3d", computed_zorder=False)
    ax.set_box_aspect((LENGTH, DEPTH, HEIGHT), zoom=1.17)
    ax.view_init(elev=18, azim=-66)
    for key, place in (("floor", lambda gu, gv: (gu, gv, 0 * gu)),
                       ("side2", lambda gu, gv: (gu, DEPTH + 0 * gu, gv)),
                       ("front", lambda gu, gv: (0 * gu, gu, gv))):
        _t, m, u, v, extent, _lab = walls[key]
        gu, gv, g = on_grid(u[m], v[m], net[m], area[m], extent)
        ax.plot_surface(*place(gu, gv), facecolors=cmap(norm(g)), rstride=1, cstride=1,
                        linewidth=0, antialiased=False, shade=False)
    for tx in TUBE_X:
        cylinder(ax, tx, TUBE_Z, TUBE_R, "0.78")
    y0, y1, z0, z1 = OUTLET
    ax.plot([LENGTH] * 5, [y0, y1, y1, y0, y0], [z0, z0, z1, z1, z0], color="0.2", lw=1.2)
    box_edges(ax)
    seg, rgba = segments(traced, [0, 1, 2])
    ax.add_collection3d(Line3DCollection(seg, colors=rgba, linewidths=1.3))
    for source, x, _f, _at in traced:
        ax.scatter(*x[0], s=9, color=GAS if source == "gas" else WALL, edgecolor="k", linewidth=0.3, depthshade=False)
    ax.set_xlim(0, LENGTH)
    ax.set_ylim(0, DEPTH)
    ax.set_zlim(0, HEIGHT)
    ax.set_axis_off()
    halo = [patheffects.withStroke(linewidth=2.2, foreground="k")]
    for x, y, z, label in ((0.0, BURNER[0], BURNER[1] + BURNER[2] + 0.12, "burner"),
                           (LENGTH, 0.5 * (y0 + y1), z1 + 0.2, "opening"),
                           (TUBE_X[3], 0.0, TUBE_Z - 0.32, "cooled tubes")):
        ax.text(x, y, z, label, fontsize=9, color="w", ha="center", path_effects=halo, zorder=10)
    gas = sum(1 for s, *_ in traced if s == "gas")
    fig.text(0.04, 0.945, f"{len(traced)} of the {bundles} bundles of one call, traced: "
             f"{gas} from the hot gas (cyan), {len(traced) - gas} from the walls (green)", fontsize=11)
    fig.text(0.04, 0.905, "a path fades as the grey gas absorbs its energy; the far walls show the net flux into the wall",
             fontsize=9, color="0.3")
    cax = fig.add_axes([0.87, 0.2, 0.022, 0.55])
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, label="net radiative flux into the wall (kW/m²)")
    path = os.path.join(out, "radiation_paths_3d.png")
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(path)


def side(out, cells, traced):
    """The paths seen from the side (x, z) over the gas temperature on the
    mid-depth plane."""
    x, y, z = (column(cells, k) for k in "xyz")
    thin = np.abs(y - DEPTH / 2) < 0.12
    gx, gz, t = smooth_grid(x[thin], z[thin], column(cells, "temperature")[thin],
                            column(cells, "volume")[thin], (LENGTH, HEIGHT), sigma=0.06)
    fig, ax = plt.subplots(figsize=(8.4, 5.6), constrained_layout=True)
    ax.set_facecolor(plt.get_cmap(RAINBOW)(0.0))
    cs = ax.contourf(gx, gz, t, levels=64, cmap=RAINBOW)
    fig.colorbar(cs, ax=ax, label="gas temperature on the mid-depth plane (K)", shrink=0.9,
                 ticks=np.arange(1300, 2000, 100))
    seg, rgba = segments(traced, [0, 2], faintest=0.06)
    rgba[:, :3] = 1.0
    ax.add_collection(matplotlib.collections.LineCollection(seg, colors=rgba, linewidths=0.9))
    for _source, x0, _f, at in traced:
        ax.plot(*x0[0, [0, 2]], "o", ms=2.6, color="w", mec="k", mew=0.3)
        if at[-1] == "escaped":
            ax.plot(*x0[-1, [0, 2]], ">", ms=5, color="w", mec="k", mew=0.4)
    for tx in TUBE_X:
        ax.add_patch(Circle((tx, TUBE_Z), TUBE_R, facecolor="0.85", edgecolor="0.25", lw=0.8, zorder=5))
    _yb, zb, rb = BURNER
    ax.plot([0, 0], [zb - rb, zb + rb], color="k", lw=4, solid_capstyle="butt")
    ax.plot([LENGTH, LENGTH], [OUTLET[2], OUTLET[3]], color="w", lw=4, solid_capstyle="butt")
    ax.set_xlim(0, LENGTH)
    ax.set_ylim(0, HEIGHT)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title("The traced bundles seen from the side, over the prescribed flame", fontsize=11)
    path = os.path.join(out, "radiation_paths_plane.png")
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(path)


def six_walls(out, f):
    """The net flux into each wall on one scale; the walls are laid out in
    two rows (floor, roof, front; the side walls and the back)."""
    norm = flux_scale(f)
    cmap = plt.get_cmap(SECOND)
    walls = walls_of(f)
    _name, _x, _y, _z, net, area = f
    fig, axes = plt.subplots(
        2, 3, figsize=(12.0, 5.8), constrained_layout=True,
        gridspec_kw={"width_ratios": [LENGTH, LENGTH, DEPTH], "height_ratios": [DEPTH, HEIGHT]},
    )
    order = [["floor", "roof", "front"], ["side0", "side2", "back"]]
    for i, row in enumerate(order):
        for j, key in enumerate(row):
            ax = axes[i][j]
            title, m, u, v, extent, (lu, lv) = walls[key]
            gu, gv, g = on_grid(u[m], v[m], net[m], area[m], extent)
            ax.contourf(gu, gv, np.clip(g, norm.vmin, norm.vmax), levels=64, cmap=cmap, norm=norm)
            if key in ("side0", "side2"):
                for tx in TUBE_X:
                    ax.add_patch(Circle((tx, TUBE_Z), TUBE_R, facecolor="0.85", edgecolor="0.25", lw=0.6))
            if key == "front":
                yb, zb, rb = BURNER
                ax.add_patch(Circle((yb, zb), rb, fill=False, edgecolor="w", lw=0.9, ls="--"))
            if key == "back":
                y0, y1, z0, z1 = OUTLET
                ax.add_patch(Rectangle((y0, z0), y1 - y0, z1 - z0, facecolor="0.9", edgecolor="0.25", lw=0.6))
            mean = np.average(net[m], weights=area[m])
            ax.set_title(f"{title}: {mean:.0f} kW/m²", fontsize=10)
            ax.set_xlim(0, extent[0])
            ax.set_ylim(0, extent[1])
            ax.set_aspect("equal")
            ax.set_xlabel(lu, fontsize=8)
            ax.set_ylabel(lv, fontsize=8)
            ax.tick_params(labelsize=7)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, shrink=0.8,
                 label="net radiative flux into the wall (kW/m²)")
    path = os.path.join(out, "radiation_paths_walls.png")
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(path)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    bundles = sys.argv[3] if len(sys.argv) > 3 else "16 million"
    os.makedirs(out, exist_ok=True)
    f = faces(read(os.path.join(run, "walls.csv")))
    traced = paths(read(os.path.join(run, "paths.csv")))
    chamber(out, f, traced, bundles)
    side(out, read(os.path.join(run, "cells.csv")), traced)
    six_walls(out, f)


if __name__ == "__main__":
    main()
