"""Hot air cooled by a water spray (Lagrangian sprays, two-way coupling with
sensible heat SP2b; gallery and recipe of the SP2b pull request).

spray_cooling.png: top, the air temperature in the duct with the droplets
coloured by their temperature at the end of the run; bottom left, the
mass-flow mean air temperature along the duct; bottom right, over time, the
heat the air gives the spray and the mean droplet temperature.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_cooling.toml --out output/spray_cooling
    python3 -m tools.plot.spray_cooling OUT_DIR output/spray_cooling
"""
import glob
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read


def field(d, name):
    if name in d:
        return d[name]
    have = ", ".join(sorted(d))
    raise SystemExit(f"fields.vtkhdf has no cell field `{name}`; it has {have}")


def parcels(path):
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        pd = g["PointData"]
        if "temperature" not in pd:
            have = ", ".join(sorted(pd))
            raise SystemExit(f"{path} has no point field `temperature`; it has {have}")
        return g["Points"][:], pd["temperature"][:], dict(f.attrs)


def plot(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    temp = field(d, "temperature")
    ux = field(d, "velocity_x")
    last = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))[-1]
    x, tp, attrs = parcels(last)
    csv = np.genfromtxt(os.path.join(run, "spray.csv"), delimiter=",", names=True, comments=None, skip_header=2)
    for c in ("gas_heat_w", "t_mean_k"):
        if c not in csv.dtype.names:
            raise SystemExit(f"spray.csv has no column `{c}`; it has {', '.join(csv.dtype.names)}")
    fig = plt.figure(figsize=(12.0, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0])
    ax = fig.add_subplot(gs[0, :])
    lo, hi = temp.min(), temp.max()
    pc = PolyCollection(quads, array=temp - 273.15, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    pc.set_clim(lo - 273.15, hi - 273.15)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="air temperature (°C)", shrink=1.0, pad=0.01)
    sc = ax.scatter(
        x[:, 0], x[:, 1], c=tp - 273.15, s=9, cmap="Greys", vmin=20, vmax=hi - 273.15,
        edgecolors="white", linewidths=0.3, zorder=3,
    )
    fig.colorbar(sc, ax=ax, label="droplet temperature (°C)", shrink=1.0, pad=0.01)
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 0.4)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Air at 300 °C cooled by a water spray, t = {attrs['time']:.2f} s ({len(x)} parcels)", fontsize=11)
    # mass-flow mean temperature per column of cells
    centre = quads.mean(axis=1)
    xs = np.round(centre[:, 0], 9)
    cols = np.unique(xs)
    tm = np.array([np.sum(ux[xs == c] * temp[xs == c]) / np.sum(ux[xs == c]) for c in cols])
    a = fig.add_subplot(gs[1, 0])
    a.plot(cols, tm - 273.15, color="tab:red", lw=2.0)
    a.axvline(0.15, color="k", ls=":", lw=1.0, label="nozzle")
    a.set_xlabel("x (m)")
    a.set_ylabel("mean air temperature (°C)")
    a.set_title("Mass-flow mean air temperature along the duct", fontsize=10)
    a.legend(fontsize=9)
    a.grid(alpha=0.3)
    b = fig.add_subplot(gs[1, 1])
    t = csv["time"]
    b.plot(t, csv["gas_heat_w"] / -1000.0, color="tab:blue", lw=1.8, label="heat from the air (kW)")
    b.set_xlabel("time (s)")
    b.set_ylabel("heat from the air (kW)", color="tab:blue")
    b.grid(alpha=0.3)
    b2 = b.twinx()
    b2.plot(t, csv["t_mean_k"] - 273.15, color="tab:orange", lw=1.5)
    b2.set_ylabel("mean droplet temperature (°C)", color="tab:orange")
    b.set_title("Heat the spray takes from the air (ledger pair)", fontsize=10)
    fig.savefig(os.path.join(out, "spray_cooling.png"), dpi=100)
    plt.close(fig)


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(out, run)


if __name__ == "__main__":
    main()
