"""Immersed bodies by Brinkman penalisation, L13 stage 1 (gallery and recipe 270).

immersed_cylinder_re20.png: the Schäfer & Turek cylinder at Re = 20 on a
plain box mesh with the cylinder as a penalised body, beside the same case
on the body-fitted mesh: the speed of both, the vorticity of the penalised
run, and a close-up of the solid fraction with the mesh and the exact
circle.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_fields OUT_DIR PENALISED_RUN FITTED_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/immersed_cylinder_re20.toml --out output/imm
    target/release/foehn flow examples/cases/flow/cylinder_re20.toml --out output/fitted
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle

from tools.plot.colormaps import RAINBOW, UNIFORM
from tools.plot.turbulence_fields import read

XMAX = 0.9  # crop of the channel (m)
CENTRE, R = (0.2, 0.2), 0.05


def crop(quads, *values):
    keep = quads[:, :, 0].min(axis=1) < XMAX
    return (quads[keep],) + tuple(v[keep] for v in values)


def circle(ax, colour="w"):
    ax.add_patch(Circle(CENTRE, R, fill=False, ec=colour, lw=1.0, ls="--"))


def panel(fig, ax, quads, values, cmap, label, title, vmin=None, vmax=None, edges="face"):
    lw = 0.15 if edges != "face" else 0.0
    pc = PolyCollection(quads, array=values, cmap=cmap, edgecolors=edges, linewidths=lw,
                        antialiased=edges != "face")
    if vmin is not None:
        pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.colorbar(pc, ax=ax, label=label, shrink=0.9, pad=0.01)
    return pc


def main():
    out, pen, fit = sys.argv[1:4]
    qp, _, dp = read(os.path.join(pen, "fields.vtkhdf"))
    qf, _, df = read(os.path.join(fit, "fields.vtkhdf"))
    sp = np.hypot(dp["velocity_x"], dp["velocity_y"])
    sf = np.hypot(df["velocity_x"], df["velocity_y"])
    qp1, sp1, wz, chi = crop(qp, sp, dp["vorticity_z"], dp["solid_fraction"])
    qf1, sf1 = crop(qf, sf)
    top = max(sp1.max(), sf1.max())

    fig, axs = plt.subplots(2, 2, figsize=(12.5, 5.4), constrained_layout=True)
    panel(fig, axs[0, 0], qp1, sp1, RAINBOW, "speed (m/s)",
          "Penalised, box mesh of 36 080 cells: speed", 0.0, top)
    circle(axs[0, 0])
    panel(fig, axs[0, 1], qf1, sf1, RAINBOW, "speed (m/s)",
          "Body-fitted, 24 576 cells: speed", 0.0, top)
    w = np.percentile(np.abs(wz), 99.0)
    panel(fig, axs[1, 0], qp1, wz, UNIFORM, "vorticity (1/s)",
          "Penalised: vorticity", -w, w)
    circle(axs[1, 0], "k")
    # close-up of the solid fraction with the cell edges
    near = (np.abs(qp[:, :, 0].mean(axis=1) - CENTRE[0]) < 0.085) & (
        np.abs(qp[:, :, 1].mean(axis=1) - CENTRE[1]) < 0.085)
    ax = axs[1, 1]
    panel(fig, ax, qp[near], dp["solid_fraction"][near], RAINBOW, "solid fraction χ",
          "Solid fraction χ per cell, exact circle dashed", 0.0, 1.0,
          edges=(0.3, 0.3, 0.3, 0.6))
    circle(ax)
    ax.set_xlim(CENTRE[0] - 0.08, CENTRE[0] + 0.08)
    ax.set_ylim(CENTRE[1] - 0.08, CENTRE[1] + 0.08)
    for a in axs.flat[:3]:
        a.set_xlim(0.0, XMAX)
        a.set_ylim(0.0, 0.41)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "immersed_cylinder_re20.png")
    fig.savefig(path, dpi=110)
    print(path)


if __name__ == "__main__":
    main()
