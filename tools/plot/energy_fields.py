"""The energy equation of stage 7a1 (gallery entry of the PR): the heated
lid-driven cavity (temperature with the streamlines), the conduction box
with k(T) (temperature and conductivity), and the Couette flow heated by
its own friction (temperature across the gap against the exact parabola
on three meshes).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.energy_fields OUT_DIR CAVITY_DIR MMS_DIR COUETTE_DIR

e.g.

    target/release/foehn flow examples/cases/flow/cavity_heated.toml --out output/cavity_heated
    target/release/foehn flow examples/cases/flow/conduction_mms.toml --out output/conduction \\
        --set 'verification.levels=[2]' --set output.vtkhdf=true
    for l in 0 1 2; do
      target/release/foehn flow examples/cases/flow/couette_viscous_heating.toml \\
          --out output/couette_$l --set "verification.levels=[$l]" \\
          --set 'output.lines=[{name="across",from=[1.0,0.0,0.05],to=[1.0,1.0,0.05],points=64}]'
    done
    python3 -m tools.plot.energy_fields docs/figures/gallery output/cavity_heated \\
        output/conduction output/couette

The figures are energy_cavity.png and energy_checks.png.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.turbulence_fields import read
from tools.plot.colormaps import RAINBOW


def grid(quads, values):
    """Cell values of a structured 2-D hex mesh as an (ny, nx) array, with
    the cell-centre coordinates."""
    c = quads.mean(axis=1)
    xs = np.unique(np.round(c[:, 0], 12))
    ys = np.unique(np.round(c[:, 1], 12))
    ix = np.searchsorted(xs, np.round(c[:, 0], 12))
    iy = np.searchsorted(ys, np.round(c[:, 1], 12))
    out = np.full((len(ys), len(xs)), np.nan)
    out[iy, ix] = values
    return xs, ys, out


def field(ax, quads, values, cmap, label, fig):
    pc = PolyCollection(quads, array=values, cmap=cmap, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max())
    ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
    ax.set_aspect("equal")
    fig.colorbar(pc, ax=ax, label=label, shrink=0.85, pad=0.02)
    return pc


def cavity(out, path):
    quads, _, d = read(os.path.join(path, "fields.vtkhdf"))
    t = d["temperature"]
    fig, ax = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
    field(ax, quads, t, "inferno", "temperature (K)", fig)
    xs, ys, ux = grid(quads, d["velocity_x"])
    _, _, uy = grid(quads, d["velocity_y"])
    ax.streamplot(xs, ys, ux, uy, color="w", linewidth=0.6, density=1.3, arrowsize=0.7)
    ax.set_title("Lid at 310 K moving right, bottom at 300 K, Re = 100, Pr = 0.71")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.savefig(os.path.join(out, "energy_cavity.png"), dpi=130)


def line(path):
    with open(path) as f:
        rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")]
    head = rows[0]
    data = np.array([[float(v) for v in r] for r in rows[1:]])
    return {k: data[:, i] for i, k in enumerate(head)}


def checks(out, mms, couette):
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), constrained_layout=True)
    quads, _, d = read(os.path.join(mms, "fields.vtkhdf"))
    t = d["temperature"]
    field(axes[0], quads, t, RAINBOW, "temperature (K)", fig)
    axes[0].set_title("Conduction with k(T) = 1 + 0.01 (T − 300)")
    k = 1.0 + 0.01 * (t - 300.0)
    field(axes[1], quads, k, RAINBOW, "conductivity k(T) (W/(m K))", fig)
    axes[1].set_title("k follows T in every outer iteration")
    ax = axes[2]
    y = np.linspace(0.0, 1.0, 201)
    ax.plot(300.0 + y * (1.0 - y), y, "k-", lw=1.2, label="exact, Br = 1")
    for lvl, m in zip((0, 1, 2), ("o", "s", "^")):
        p = os.path.join(f"{couette}_{lvl}", "line_across.csv")
        if not os.path.exists(p):
            continue
        r = line(p)
        cells = 8 * 2**lvl
        yc = (np.arange(cells) + 0.5) / cells
        tc = np.interp(yc, r["y"], r["temperature"])
        ax.plot(tc, yc, m, ms=4, mfc="none", label=f"{cells} cells across")
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("y/h")
    ax.set_title("Couette flow heated by its own friction")
    ax.legend(loc="lower right", fontsize=8)
    for a in axes[:2]:
        a.set_xlabel("x (m)")
        a.set_ylabel("y (m)")
    fig.savefig(os.path.join(out, "energy_checks.png"), dpi=120)


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    out, cav, mms, couette = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    cavity(out, cav)
    checks(out, mms, couette)


if __name__ == "__main__":
    main()
