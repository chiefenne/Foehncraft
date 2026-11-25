"""Droplets spread by the turbulence of a duct (Lagrangian sprays, eddy-lifetime
dispersion SP2c; gallery and recipe of the SP2c pull request).

spray_dispersion.png: top, the turbulent kinetic energy k of the air in the
duct with the droplets at the end of the run, coloured by their diameter;
bottom left, the sideways spread (standard deviation of y about the nozzle
axis) along the duct for three size classes, pooled over the parcel files of
the second half of the run; bottom right, where the droplets of each class
cross the plane x = 2.5 m.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_dispersion.toml --out output/spray_dispersion
    python3 -m tools.plot.spray_dispersion OUT_DIR output/spray_dispersion
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

AXIS_Y = 0.15
CLASSES = [
    (20e-6, 60e-6, "20–60 µm", "#2a9df4"),
    (60e-6, 150e-6, "60–150 µm", "#f4a62a"),
    (150e-6, 300e-6, "150–300 µm", "#d6336c"),
]


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


def plot(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    k = field(d, "k")
    files = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))
    if not files:
        raise SystemExit(f"no parcel files in {run}/spray")
    x, dia, attrs = parcels(files[-1])
    # pool the second half of the run (the plume is steady there)
    pool_x, pool_d = [], []
    for f in files[len(files) // 2:]:
        xi, di, _ = parcels(f)
        pool_x.append(xi)
        pool_d.append(di)
    px = np.concatenate(pool_x)
    pd = np.concatenate(pool_d)

    fig = plt.figure(figsize=(12.0, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.62, 1.0])
    ax = fig.add_subplot(gs[0, :])
    pc = PolyCollection(quads, array=k, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="k of the air (m²/s²)", shrink=1.0, pad=0.01)
    sc = ax.scatter(
        x[:, 0], x[:, 1], c=dia * 1e6, s=7, cmap="Greys", vmin=0, vmax=300,
        edgecolors="white", linewidths=0.25, zorder=3,
    )
    fig.colorbar(sc, ax=ax, label="droplet diameter (µm)", shrink=1.0, pad=0.01)
    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(0.0, 0.3)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(
        f"Droplets shot along the axis, spread by the eddies, t = {attrs['time']:.2f} s ({len(x)} parcels)",
        fontsize=11,
    )

    edges = np.linspace(0.2, 3.0, 29)
    centres = 0.5 * (edges[1:] + edges[:-1])
    a = fig.add_subplot(gs[1, 0])
    for lo, hi, label, colour in CLASSES:
        sel = (pd >= lo) & (pd < hi)
        sig = []
        for e0, e1 in zip(edges[:-1], edges[1:]):
            s = sel & (px[:, 0] >= e0) & (px[:, 0] < e1)
            sig.append(np.sqrt(np.mean((px[s, 1] - AXIS_Y) ** 2)) * 100 if s.sum() > 20 else np.nan)
        a.plot(centres, sig, "-o", ms=3, color=colour, label=label)
    a.set_xlabel("x (m)")
    a.set_ylabel("sideways spread σ_y (cm)")
    a.set_title("Spread about the nozzle axis by droplet size", fontsize=11)
    a.grid(alpha=0.3)
    a.legend(loc="upper left")

    b = fig.add_subplot(gs[1, 1])
    bins = np.linspace(0.0, 0.3, 13)
    band = (px[:, 0] > 2.3) & (px[:, 0] < 2.7)
    for lo, hi, label, colour in CLASSES:
        s = band & (pd >= lo) & (pd < hi)
        if s.sum() > 0:
            b.hist(px[s, 1], bins=bins, density=True, histtype="step", lw=2, color=colour, label=label)
    b.axvline(AXIS_Y, color="0.3", ls=":", lw=1)
    b.set_xlabel("y (m)")
    b.set_ylabel("probability density (1/m)")
    b.set_title("Where the droplets are at x = 2.3 to 2.7 m", fontsize=11)
    b.grid(alpha=0.3)
    b.legend(loc="upper right")
    fig.savefig(out, dpi=100)
    print(f"wrote {out}")


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(os.path.join(out, "spray_dispersion.png"), run)


if __name__ == "__main__":
    main()
