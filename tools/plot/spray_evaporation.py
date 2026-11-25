"""Water droplets evaporating in hot air (Lagrangian sprays, evaporation
SP4 and SP5; gallery and recipe of the SP5 pull request).

spray_evaporation.png: top, the water vapour mass fraction in the duct with
the droplets as dots sized and coloured by their diameter at the end of the
run; bottom left, the mass-flow mean air temperature and vapour mass
fraction along the duct; bottom right, over time, the vapour the air
receives from the spray per second (ledger pair) against the 2 g/s the
nozzle sprays, and the Sauter mean diameter of the droplets in the flow.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_evaporation.toml --out output/spray_evaporation
    python3 -m tools.plot.spray_evaporation OUT_DIR output/spray_evaporation
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
        if "diameter" not in pd:
            have = ", ".join(sorted(pd))
            raise SystemExit(f"{path} has no point field `diameter`; it has {have}")
        return g["Points"][:], pd["diameter"][:], dict(f.attrs)


def column(csv, c):
    if c not in csv.dtype.names:
        raise SystemExit(f"spray.csv has no column `{c}`; it has {', '.join(csv.dtype.names)}")
    return csv[c]


def plot(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    y = field(d, "Y_H2O")
    temp = field(d, "temperature")
    ux = field(d, "velocity_x")
    last = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))[-1]
    x, dp, attrs = parcels(last)
    csv = np.genfromtxt(os.path.join(run, "spray.csv"), delimiter=",", names=True, comments=None, skip_header=2)
    fig = plt.figure(figsize=(12.0, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.0])
    ax = fig.add_subplot(gs[0, :])
    pc = PolyCollection(quads, array=y * 100.0, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="water vapour mass fraction (%)", shrink=1.0, pad=0.01)
    sc = ax.scatter(
        x[:, 0], x[:, 1], c=dp * 1e6, s=4.0 + (dp * 1e6) ** 2 / 40.0, cmap="Greys", vmin=0.0, vmax=80.0,
        edgecolors="white", linewidths=0.3, zorder=3,
    )
    fig.colorbar(sc, ax=ax, label="droplet diameter (µm)", shrink=1.0, pad=0.01)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.2)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(
        f"Water droplets evaporating in air at 300 °C, t = {attrs['time']:.2f} s ({len(x)} parcels)",
        fontsize=11,
    )
    # mass-flow means per column of cells
    centre = quads.mean(axis=1)
    xs = np.round(centre[:, 0], 9)
    cols = np.unique(xs)
    tm = np.array([np.sum(ux[xs == c] * temp[xs == c]) / np.sum(ux[xs == c]) for c in cols])
    ym = np.array([np.sum(ux[xs == c] * y[xs == c]) / np.sum(ux[xs == c]) for c in cols])
    a = fig.add_subplot(gs[1, 0])
    a.plot(cols, tm - 273.15, color="tab:red", lw=2.0)
    a.axvline(0.1, color="k", ls=":", lw=1.0)
    a.text(0.105, tm.min() - 273.15, "nozzle", fontsize=9, va="bottom")
    a.set_xlabel("x (m)")
    a.set_ylabel("mean air temperature (°C)", color="tab:red")
    a.grid(alpha=0.3)
    a2 = a.twinx()
    a2.plot(cols, ym * 100.0, color="tab:blue", lw=2.0)
    a2.set_ylabel("mean vapour mass fraction (%)", color="tab:blue")
    a.set_title("Mass-flow means along the duct", fontsize=10)
    b = fig.add_subplot(gs[1, 1])
    t = column(csv, "time")
    b.plot(t, column(csv, "gas_vapour_kg_s") * 1000.0, color="tab:blue", lw=1.8, label="vapour to the air")
    b.axhline(2.0, color="k", ls="--", lw=1.0, label="water sprayed")
    b.set_xlabel("time (s)")
    b.set_ylabel("mass flow (g/s)")
    b.set_ylim(0.0, 2.3)
    b.grid(alpha=0.3)
    b.legend(fontsize=9, loc="lower right")
    b2 = b.twinx()
    b2.plot(t, column(csv, "d32_m") * 1e6, color="tab:gray", lw=1.5)
    b2.set_ylabel("Sauter mean diameter in the flow (µm)", color="tab:gray")
    b.set_title("Vapour the spray gives the air (ledger pair)", fontsize=10)
    fig.savefig(os.path.join(out, "spray_evaporation.png"), dpi=100)
    plt.close(fig)


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(out, run)


if __name__ == "__main__":
    main()
