"""Liquid and gas through datums of a sprayed duct (sprays SP1b; gallery and recipe of the datums).

spray_datums.png: top, the air speed with the parcels coloured by their
diameter and the datums of the case drawn in; bottom left, the cumulative
liquid through the datums over time (datums.csv); bottom right, the mass
distribution over the droplet diameter of the liquid that left through the
outlet and of the liquid that landed on the floor (spray_hist_<datum>.csv),
with the input Rosin-Rammler distribution of the nozzle.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_datums.toml --out output/spray_datums
    python3 -m tools.plot.spray_datums OUT_DIR output/spray_datums
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
from tools.plot.turbulence_fields import read


def table(path):
    rows = [r for r in open(path, encoding="utf-8") if not r.startswith("#")]
    head = rows[0].strip().split(",")
    data = np.array([[float(x) if x != "inf" else np.inf for x in r.strip().split(",")] for r in rows[1:]])
    cols = {k: data[:, i] for i, k in enumerate(head)}
    return cols


def need(cols, key, path):
    if key not in cols:
        raise SystemExit(f"{path}: no column {key}; it has {', '.join(cols)}")
    return cols[key]


def parcels(path):
    import h5py

    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        return g["Points"][:], g["PointData/diameter"][:], dict(f.attrs)


def figure(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    last = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))[-1]
    x, diam, attrs = parcels(last)
    dpath = os.path.join(run, "datums.csv")
    dc = table(dpath)
    fig = plt.figure(figsize=(11.0, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85])
    ax = fig.add_subplot(gs[0, :])
    pc = PolyCollection(quads, array=d["speed"], cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="air speed (m/s)", shrink=0.9, pad=0.01)
    order = np.argsort(diam)
    sc = ax.scatter(
        x[order, 0],
        x[order, 1],
        c=diam[order] * 1e6,
        s=4 + 0.12 * diam[order] * 1e6,
        cmap="inferno",
        edgecolors="white",
        linewidths=0.3,
        zorder=3,
    )
    fig.colorbar(sc, ax=ax, label="droplet diameter (µm)", shrink=0.9, pad=0.01)
    style = dict(color="white", lw=2.2, zorder=4)
    ax.plot([0.6, 0.6], [0.0, 0.4], ls="--", **style)
    ax.plot([1.2, 1.2], [0.0, 0.4], **style)
    ax.plot([0.9, 0.9], [0.1, 0.3], ls=":", **style)
    ax.plot([1.6, 1.6], [0.0, 0.4], ls="-.", **style)
    ax.add_patch(plt.Rectangle((0.0, 0.2), 0.3, 0.2, fill=False, ec="white", lw=1.6, ls="--", zorder=4))
    for xl, label in [(0.6, "x06 cut"), (0.9, "disc"), (1.2, "x12 faces"), (1.6, "x16 STL")]:
        ax.text(xl + 0.015, 0.36, label, color="white", fontsize=9, zorder=5)
    ax.text(0.31, 0.21, "nozzle_box", color="white", fontsize=9, zorder=5)
    ax.plot([0.06, 0.14], [0.4, 0.4], color="tab:red", lw=5, solid_capstyle="butt", zorder=4)
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 0.4)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(
        f"{len(x)} parcels at t = {attrs['time']:.2f} s: water, Rosin-Rammler sizes (X = 120 µm, n = 2.5), "
        "air at 4 m/s; datums in white",
        fontsize=11,
    )
    ax = fig.add_subplot(gs[1, 0])
    t = need(dc, "time", dpath)
    for name, label, c in [
        ("x06", "x06 (plane cut)", "tab:blue"),
        ("x12", "x12 (plane on faces)", "tab:green"),
        ("x16", "x16 (STL surface)", "tab:purple"),
        ("outlet", "outlet (boundary)", "k"),
        ("floor", "floor (lower wall)", "tab:orange"),
    ]:
        net = need(dc, f"{name}.liquid_kg_pos", dpath) - need(dc, f"{name}.liquid_kg_neg", dpath)
        ax.plot(t, net * 1e3, lw=2, color=c, label=label)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("liquid through the datum (g)")
    ax.set_title("Cumulative liquid through each datum", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[1, 1])
    for name, label, c in [("outlet", "left through the outlet", "k"), ("floor", "stuck on the floor", "tab:orange")]:
        hp = os.path.join(run, f"spray_hist_{name}.csv")
        h = table(hp)
        b = need(h, "bin", hp)
        inner = b >= 0
        inner &= b < b.max()
        lo, hi, m = h["d_lo"][inner], h["d_hi"][inner], h["mass"][inner]
        keys = np.unique(lo)
        mass = np.array([m[lo == k].sum() for k in keys])
        his = np.array([hi[lo == k][0] for k in keys])
        frac = mass / mass.sum() / np.log(his / keys)
        ax.step(keys * 1e6, frac, where="post", lw=2, color=c, label=label)
    dd = np.geomspace(10e-6, 400e-6, 400)
    X, n = 120e-6, 2.5
    pdf = n * (dd / X) ** n * np.exp(-((dd / X) ** n))
    pdf /= np.trapezoid(pdf, np.log(dd))
    ax.plot(dd * 1e6, pdf, color="tab:red", lw=1.5, ls="--", label="nozzle (input)")
    ax.set_xscale("log")
    ax.set_xlim(8, 500)
    ax.set_xlabel("droplet diameter (µm)")
    ax.set_ylabel("mass per ln d (normalised)")
    ax.set_title("Where the liquid went, by droplet size", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3, which="both")
    fig.savefig(os.path.join(out, "spray_datums.png"), dpi=100)
    plt.close(fig)


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    figure(out, run)


if __name__ == "__main__":
    main()
