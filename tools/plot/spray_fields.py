"""Water droplets in an air duct (Lagrangian sprays, SP1; gallery and recipe of the spray layer).

spray_duct.png: the air speed of the transient SST run with the parcels of
the spray coloured by their speed, and below it the mass book of the
liquid (released, in the flow, stuck on the lower wall, left through the
outlet) over time.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_duct.toml --out output/spray_duct
    python3 -m tools.plot.spray_fields OUT_DIR output/spray_duct
"""
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.turbulence_fields import read


def parcels(path):
    import h5py

    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        return g["Points"][:], g["PointData/velocity"][:], g["PointData/diameter"][:], dict(f.attrs)


def duct(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    last = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))[-1]
    x, v, diam, attrs = parcels(last)
    speed = np.hypot(v[:, 0], v[:, 1])
    book = np.genfromtxt(os.path.join(run, "spray.csv"), delimiter=",", names=True, comments=None, skip_header=2)
    fig, axs = plt.subplots(
        2, 1, figsize=(11.0, 6.4), constrained_layout=True, gridspec_kw={"height_ratios": [1.0, 0.62]}
    )
    ax = axs[0]
    pc = PolyCollection(quads, array=d["speed"], cmap="Blues", edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="air speed (m/s)", shrink=0.9, pad=0.01)
    sc = ax.scatter(
        x[:, 0], x[:, 1], c=speed, s=14, cmap="turbo", edgecolors="k", linewidths=0.25, zorder=3
    )
    cb = fig.colorbar(sc, ax=ax, label="droplet speed (m/s)", shrink=0.9, pad=0.01)
    cb.ax.tick_params(labelsize=8)
    ax.plot([0.06, 0.14], [0.4, 0.4], color="tab:red", lw=5, solid_capstyle="butt", zorder=4)
    ax.annotate("nozzle", (0.1, 0.4), (0.1, 0.45), color="tab:red", ha="center", fontsize=10)
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 0.4)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    t = attrs["time"]
    ax.set_title(
        f"{len(x)} parcels of 150 µm water droplets at t = {t:.2f} s in air flowing at 4 m/s", fontsize=11
    )
    ax = axs[1]
    ax.plot(book["time"], book["injected_kg"] * 1e3, color="k", lw=2, label="released")
    ax.plot(book["time"], book["in_flow_kg"] * 1e3, color="tab:blue", lw=2, label="in the flow")
    ax.plot(book["time"], book["deposited_kg"] * 1e3, color="tab:orange", lw=2, label="stuck on the lower wall")
    ax.plot(book["time"], book["escaped_kg"] * 1e3, color="tab:green", lw=2, label="left through the outlet")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("liquid (g)")
    ax.set_title("Mass book of the liquid: released = in the flow + stuck + left + lost", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out, "spray_duct.png"), dpi=100)
    plt.close(fig)


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    duct(out, run)


if __name__ == "__main__":
    main()
