"""The deforming mesh on ranks and from a restart file (stage MS0 part 2,
gallery and recipe 307).

motion_channel_ranks.png: the breathing channel of
motion_channel_ranks.toml. Left: the speed at a largest swing of the
mesh on the moved mesh, with the cuts between the three ranks as white
lines. Middle: the dye at the end of the run. Right: the axial velocity at
the probe in the middle against time, on one rank, on three ranks and
continued from restart_000075.h5 on two ranks: the three curves lie on top
of each other because the runs are bitwise equal.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.mesh_motion_ranks OUT_DIR ONE_RANK THREE_RANKS CONTINUED

with the output directories of

    F=examples/cases/flow/motion_channel_ranks.toml
    target/release/foehn flow $F --ranks 1 --out output/one --set 'output.fields_times=[1.0625]'
    target/release/foehn flow $F --out output/three --set 'output.fields_times=[1.0625]'
    target/release/foehn flow $F --ranks 2 --restart output/three/restart_000075.h5 --out output/cont
"""
import csv
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read


def probes(run):
    with open(os.path.join(run, "probes.csv")) as f:
        rows = list(csv.reader(line for line in f if not line.startswith("#")))
    head, body = rows[0], rows[1:]
    return {k: np.array([float(r[i]) for r in body]) for i, k in enumerate(head)}


def ranks(run):
    return read(os.path.join(run, "ranks.vtkhdf"))[2]["rank"]


def cuts(quads, rank):
    """The edges between cells of different ranks."""
    owner = {}
    segs = []
    for c, q in enumerate(quads):
        for i in range(4):
            a, b = q[i], q[(i + 1) % 4]
            key = tuple(sorted([tuple(np.round(a, 9)), tuple(np.round(b, 9))]))
            if key in owner and rank[owner[key]] != rank[c]:
                segs.append([a, b])
            owner[key] = c
    return segs


def field(ax, quads, values, label, cmap, rank=None, vmin=None, vmax=None):
    pc = PolyCollection(quads, array=values, cmap=cmap, edgecolors="k", linewidths=0.15)
    if vmin is not None:
        pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    if rank is not None:
        ax.add_collection(LineCollection(cuts(quads, rank), colors="w", linewidths=2.0))
    ax.set_xlim(-0.02, 2.02)
    ax.set_ylim(-0.02, 0.52)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    plt.colorbar(pc, ax=ax, label=label, orientation="horizontal", pad=0.22, shrink=0.7, aspect=30)


def main():
    out, one, three, cont = sys.argv[1:5]
    fig = plt.figure(figsize=(16.0, 3.7), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.25, 1.0])
    # the field file of step 107 (t = 1.07 s, fields_times = [1.0625])
    swing = min(
        glob.glob(os.path.join(three, "fields_0*.vtkhdf")),
        key=lambda f: abs(int(os.path.basename(f)[7:13]) - 107),
    )
    quads, _, d = read(swing)
    ax = fig.add_subplot(gs[0])
    field(ax, quads, d["speed"], "speed (m/s)", RAINBOW, rank=ranks(three))
    ax.set_title("speed at a largest swing (t = 1.07 s); white: cuts between the 3 ranks", fontsize=10)
    quads, _, d = read(os.path.join(three, "fields.vtkhdf"))
    ax = fig.add_subplot(gs[1])
    field(ax, quads, d["dye"], "dye (-)", "viridis", vmin=0.0, vmax=1.0)
    ax.set_title("dye at the end (t = 1.5 s), 3 ranks", fontsize=10)
    ax = fig.add_subplot(gs[2])
    for run, lab, kw in [
        (one, "1 rank", dict(color="#1f5fbf", lw=3.0)),
        (three, "3 ranks", dict(color="#f08c00", lw=1.6)),
        (cont, "from restart_000075.h5 on 2 ranks", dict(color="k", lw=0, marker="o", ms=3)),
    ]:
        p = probes(run)
        ax.plot(p["time"], p["mid_u"], label=lab, **kw)
    ax.axvline(0.75, color="0.6", ls="--", lw=1)
    ax.text(0.76, ax.get_ylim()[0], " restart", color="0.4", va="bottom")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("u at the probe (1, 0.25) m (m/s)")
    ax.set_title("bitwise equal: one rank, three ranks, continued run", fontsize=10)
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)
    path = os.path.join(out, "motion_channel_ranks.png")
    fig.savefig(path, dpi=110)
    print(path)


if __name__ == "__main__":
    main()
