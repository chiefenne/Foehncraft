"""Water sprayed into a chamber of hot, still air, with the liquid length and
the vapour penetration of the ECN modelling standard ([spray.ecn], Spray A
stage SA3, D315; gallery and recipe of the SA3 pull request).

spray_chamber_ecn.png: left, the water vapour mass fraction in the chamber
with the droplets as dots at the end of the run, the two liquid lengths and
the vapour penetration of that step as lines; middle, the projected
liquid volume along the axis (the largest across, the liquid volume per
area seen along z, in bins of [spray.ecn] bin, recomputed here from the
parcel file for the picture) with the thresholds and the liquid lengths
the solver wrote; right, the liquid lengths and the vapour
penetration over time from spray.csv.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_chamber_ecn.toml --out output/spray_chamber_ecn
    python3 -m tools.plot.spray_chamber_ecn OUT_DIR output/spray_chamber_ecn
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

# [spray.ecn] and the nozzle of examples/cases/flow/spray_chamber_ecn*.toml
PLV = [3e-8, 3e-7]
BIN = 2.5e-3
NOZZLE = np.array([0.0005, 0.08])


def field(d, name):
    if name in d:
        return d[name]
    have = ", ".join(sorted(d))
    raise SystemExit(f"fields.vtkhdf has no cell field `{name}`; it has {have}")


def parcels(path):
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        pd = g["PointData"]
        for c in ["diameter", "droplets"]:
            if c not in pd:
                have = ", ".join(sorted(pd))
                raise SystemExit(f"{path} has no point field `{c}`; it has {have}")
        return g["Points"][:], pd["diameter"][:], pd["droplets"][:], dict(f.attrs)


def column(csv, c):
    if c not in csv.dtype.names:
        raise SystemExit(f"spray.csv has no column `{c}`; it has {', '.join(csv.dtype.names)}")
    return csv[c]


def profile(x, d, n):
    """The largest PLV across per axial bin, as the solver bins it."""
    s = x[:, 0] - NOZZLE[0]
    t = x[:, 1] - NOZZLE[1]
    v = n * np.pi * d**3 / 6.0
    keep = s >= 0.0
    i = np.floor(s[keep] / BIN).astype(int)
    j = np.floor(t[keep] / BIN).astype(int)
    if len(i) == 0:
        return np.zeros(0)
    j0 = j.min()
    grid = np.zeros((i.max() + 1, j.max() - j0 + 1))
    np.add.at(grid, (i, j - j0), v[keep])
    return grid.max(axis=1) / BIN**2


def plot(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    y = field(d, "Y_H2O")
    files = glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf"))
    if not files:
        raise SystemExit(f"{run}/spray has no parcel files: run the case first (see the docstring)")
    x, dp, n, attrs = parcels(max(files))
    csv = np.genfromtxt(os.path.join(run, "spray.csv"), delimiter=",", names=True, comments=None, skip_header=2)
    t = column(csv, "time")
    # numpy drops the dots of the column names
    l1 = column(csv, "nozzleliquid_length_1_m")
    l2 = column(csv, "nozzleliquid_length_2_m")
    vp = column(csv, "nozzlevapour_penetration_m")
    k = int(np.argmin(np.abs(t - attrs["time"])))

    fig = plt.figure(figsize=(16.0, 4.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    pc = PolyCollection(quads, array=y * 100.0, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="water vapour mass fraction (%)", shrink=0.8, pad=0.01)
    ax.scatter(x[:, 0], x[:, 1], s=1.5, color="white", alpha=0.7, lw=0, zorder=3)
    for val, ls, lab in [(l1[k], "--", f"liquid length, PLV {PLV[0]:g} m"), (l2[k], "-", f"liquid length, PLV {PLV[1]:g} m")]:
        ax.axvline(NOZZLE[0] + val, color="k", ls=ls, lw=1.4, label=f"{lab}: {val * 1000:.0f} mm")
    ax.axvline(NOZZLE[0] + vp[k], color="crimson", ls="-.", lw=1.8, label=f"vapour penetration (0.1 %): {vp[k] * 1000:.0f} mm")
    ax.set_xlim(0.0, 0.25)
    ax.set_ylim(0.0, 0.16)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.85)
    ax.set_title(
        f"Water sprayed into still air at 300 °C, t = {attrs['time'] * 1000:.0f} ms ({len(x)} parcels)",
        fontsize=10,
    )

    a = fig.add_subplot(gs[0, 1])
    p = profile(x, dp, n)
    sc = (np.arange(len(p)) + 0.5) * BIN
    a.semilogy(sc * 1000.0, np.where(p > 0, p, np.nan), color="tab:blue", lw=1.6, drawstyle="steps-mid",
               label="largest PLV across (from the parcels)")
    for th, val, ls in [(PLV[0], l1[k], "--"), (PLV[1], l2[k], "-")]:
        a.axhline(th, color="0.4", ls=":", lw=1.0)
        a.axvline(val * 1000.0, color="k", ls=ls, lw=1.2)
    xmax = max(sc[-1] * 1000.0 + 10.0, 50.0) if len(sc) else 100.0
    for th in PLV:
        a.text(xmax * 0.98, th * 1.12, f"threshold {th:g} m", fontsize=8, color="0.3", ha="right")
    a.set_xlim(0.0, xmax)
    a.set_xlabel("distance from the nozzle along the axis (mm)")
    a.set_ylabel("projected liquid volume (m³/m²)")
    a.set_title("Projected liquid volume and the liquid lengths (lines)", fontsize=10)
    a.grid(alpha=0.3, which="both")
    a.legend(fontsize=8, loc="lower left")

    b = fig.add_subplot(gs[0, 2])
    b.plot(t * 1000.0, l1 * 1000.0, color="k", ls="--", lw=1.3, label=f"liquid length, PLV {PLV[0]:g} m")
    b.plot(t * 1000.0, l2 * 1000.0, color="k", lw=1.6, label=f"liquid length, PLV {PLV[1]:g} m")
    b.plot(t * 1000.0, vp * 1000.0, color="crimson", ls="-.", lw=1.8, label="vapour penetration (0.1 %)")
    b.set_xlabel("time after the start of injection (ms)")
    b.set_ylabel("distance from the nozzle (mm)")
    b.set_title("The liquid stops, the vapour goes on", fontsize=10)
    b.grid(alpha=0.3)
    b.legend(fontsize=8, loc="upper left")
    fig.savefig(os.path.join(out, "spray_chamber_ecn.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: python3 -m tools.plot.spray_chamber_ecn OUT_DIR RUN_DIR")
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(out, run)


if __name__ == "__main__":
    main()
