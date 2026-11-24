"""A water spray that drives the air (Lagrangian sprays, two-way coupling SP2; gallery and recipe 280).

spray_two_way.png: left, the air speed with arrows of the air velocity and the
droplets coloured by their diameter at the end of the run; right, over time,
the force of the spray on the air, the parcel passes of each step and the
coupling lag of its last pass.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_two_way.toml --out output/spray_two_way
    python3 -m tools.plot.spray_two_way OUT_DIR output/spray_two_way
"""
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.spray_fields import parcels
from tools.plot.turbulence_fields import read


def field(d, name):
    if name in d:
        return d[name]
    have = ", ".join(sorted(d))
    raise SystemExit(f"fields.vtkhdf has no cell field `{name}`; it has {have}")


def plot(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    vel = np.column_stack([field(d, "velocity_x"), field(d, "velocity_y")])
    speed = np.hypot(vel[:, 0], vel[:, 1])
    last = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))[-1]
    x, _, diam, attrs = parcels(last)
    csv = np.genfromtxt(os.path.join(run, "spray.csv"), delimiter=",", names=True, comments=None, skip_header=2)
    for c in ("passes", "coupling_lag", "gas_force_y"):
        if c not in csv.dtype.names:
            raise SystemExit(f"spray.csv has no column `{c}`; it has {', '.join(csv.dtype.names)}")
    fig = plt.figure(figsize=(12.0, 6.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, width_ratios=[1.05, 1.0])
    ax = fig.add_subplot(gs[:, 0])
    pc = PolyCollection(quads, array=speed, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="air speed (m/s)", shrink=0.8, pad=0.01)
    centre = quads.mean(axis=1)
    # every third cell in each direction: arrows of the air velocity
    nx = len(np.unique(np.round(centre[:, 0], 9)))
    ix = np.arange(len(centre))
    pick = ((ix % nx) % 3 == 1) & ((ix // nx) % 3 == 1)
    ax.quiver(
        centre[pick, 0], centre[pick, 1], vel[pick, 0], vel[pick, 1],
        color="white", scale=25.0, width=0.0025, alpha=0.85, zorder=2,
    )
    sc = ax.scatter(
        x[:, 0], x[:, 1], c=diam * 1e6, s=6, cmap="magma", vmin=50, vmax=400,
        edgecolors="none", zorder=3,
    )
    fig.colorbar(sc, ax=ax, label="droplet diameter (µm)", shrink=0.8, pad=0.01)
    ax.plot([0.55, 0.65], [1.5, 1.5], color="tab:red", lw=5, solid_capstyle="butt", zorder=4)
    ax.set_xlim(0.0, 1.2)
    ax.set_ylim(0.0, 1.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Spray pulling the air at t = {attrs['time']:.2f} s ({len(x)} parcels)", fontsize=11)
    t = csv["time"]
    a = fig.add_subplot(gs[0, 1])
    a.plot(t, -csv["gas_force_y"], color="tab:blue", lw=1.8)
    a.set_ylabel("downward force\non the air (N)")
    a.set_title("Momentum the spray hands to the air (ledger pair, D244)", fontsize=10)
    a.grid(alpha=0.3)
    b = fig.add_subplot(gs[1, 1], sharex=a)
    b.step(t, csv["passes"], where="post", color="tab:green", lw=1.5)
    b.set_ylabel("parcel passes\nof the step")
    b.set_ylim(0.5, 3.5)
    b.set_yticks([1, 2, 3])
    b.set_title("Passes inside the time step (converged, at most 3; D252)", fontsize=10)
    b.grid(alpha=0.3)
    c = fig.add_subplot(gs[2, 1], sharex=a)
    c.semilogy(t, np.maximum(csv["coupling_lag"], 1e-12), color="tab:purple", lw=1.2)
    c.axhline(1e-3, color="k", ls="--", lw=1.0, label="tolerance")
    c.set_ylabel("coupling lag\nof the last pass")
    c.set_xlabel("time (s)")
    c.legend(fontsize=9, loc="upper right")
    c.grid(alpha=0.3, which="both")
    fig.savefig(os.path.join(out, "spray_two_way.png"), dpi=100)
    plt.close(fig)


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(out, run)


if __name__ == "__main__":
    main()
