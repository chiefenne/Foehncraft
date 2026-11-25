"""Topology search of a pipe bend, L13 stage 1 (gallery and recipe 289).

immersed_topology_bend.png: the design found by the density-field topology
search (Borrvall & Petersson 2003 bend, Stokes flow): the solid fraction
1 − γ̃ with the extracted surface γ̃ = ½ (topology.stl, marching
tetrahedra), the speed of the flow in the final design, and the objective
and grey level over the design updates.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_topology OUT_PNG RUN_DIR

with the output directory of

    target/release/foehn flow examples/cases/flow/immersed_topology_bend.toml --out RUN_DIR
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from tools.plot.colormaps import RAINBOW, UNIFORM
from tools.plot.turbulence_fields import read


def need(data, name, path):
    if name not in data:
        sys.exit(f"{path}: no field `{name}`; the file has: {', '.join(sorted(data))}")
    return data[name]


def stl_section(path):
    """The edges of topology.stl in the plane of its smallest z."""
    pts = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if s and s[0] == "vertex":
                pts.append([float(v) for v in s[1:]])
    tri = np.array(pts).reshape(-1, 3, 3)
    z0 = tri[:, :, 2].min()
    segs = []
    for t in tri:
        on = np.abs(t[:, 2] - z0) < 1e-9
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if on[a] and on[b]:
                segs.append([t[a, :2], t[b, :2]])
    return np.array(segs), len(tri)


def history(path):
    rows = [ln for ln in open(path) if ln[0].isdigit()]
    return np.array([[float(v) for v in ln.split(",")] for ln in rows])


def main(argv):
    if len(argv) != 2:
        sys.exit(__doc__)
    out, run = argv
    quads, _, data = read(os.path.join(run, "fields.vtkhdf"))
    solid = need(data, "solid_fraction", run)
    speed = need(data, "speed", run)
    segs, ntri = stl_section(os.path.join(run, "topology.stl"))
    h = history(os.path.join(run, "topology.csv"))

    fig = plt.figure(figsize=(15, 4.8), dpi=150)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05], wspace=0.42)
    for k, (vals, cmap, label, title) in enumerate([
        (solid, UNIFORM, r"solid fraction $1-\tilde\gamma$", r"Design: solid fraction and surface $\tilde\gamma = 1/2$"),
        (speed, RAINBOW, "|u| (m/s)", "Speed in the final design"),
    ]):
        ax = fig.add_subplot(gs[0, k])
        pc = PolyCollection(quads, array=vals, cmap=cmap, edgecolors="face", linewidths=0)
        if k == 0:
            pc.set_clim(0, 1)
        ax.add_collection(pc)
        ax.add_collection(LineCollection(segs, colors="w" if k == 1 else "k", linewidths=1.1))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_title(title, fontsize=10)
        cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(label)
        # inlet and outlet
        ax.annotate("in", (0.0, 0.8), (0.07, 0.95), color="k" if k == 0 else "w",
                    arrowprops=dict(arrowstyle="->", color="k" if k == 0 else "w"), fontsize=9)
        ax.annotate("out", (0.8, 0.0), (0.88, 0.12), color="k" if k == 0 else "w",
                    arrowprops=dict(arrowstyle="->", color="k" if k == 0 else "w"), fontsize=9)
    ax = fig.add_subplot(gs[0, 2])
    ax.semilogy(h[:, 0], h[:, 2], "o-", ms=3, color="#1f5fa8", label="objective Φ (W)")
    ax.set_xlabel("design update")
    ax.set_ylabel("Φ (W)")
    ax.grid(alpha=0.3, which="both")
    ax2 = ax.twinx()
    ax2.plot(h[:, 0], h[:, 5], "s-", ms=3, color="#d9541e", label="grey level")
    ax2.set_ylabel(r"grey level $4\tilde\gamma(1-\tilde\gamma)$")
    ax2.set_ylim(0, 1)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines], loc="upper right", fontsize=8)
    ax.set_title(f"{len(h)} design updates, fluid fraction {h[-1, 3]:.3f}", fontsize=10)
    fig.suptitle("Topology search of a pipe bend (Stokes flow, Brinkman density field, "
                 f"{len(quads)} cells; topology.stl {ntri} triangles)", fontsize=11)
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.13, top=0.85)
    fig.savefig(out, facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
