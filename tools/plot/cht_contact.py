"""Contact resistance and orthotropic solids of stage 7d2 (gallery entry 160).

cht_circuit_board.png: a 1 W chip on a thermal pad on a circuit board in an
air channel, the temperature with the board's copper planes (k = 25 along,
0.35 through the board) and with a plain laminate (k = 0.35), and the
temperature along the top of the board for both and without the pad.
cht_orthotropic_patch.png: the patch test on a perturbed mesh (an
orthotropic laminate with axes turned by 30°, a fluid at rest, a block
behind a contact resistance), the temperature with its isotherms and the
error against the exact piecewise-linear solution.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_contact OUT_DIR BOARD_RUN PLAIN_RUN NOPAD_RUN PATCH_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/cht_circuit_board.toml --out output/cht_board
    target/release/foehn flow examples/cases/flow/cht_circuit_board.toml --out output/cht_board_plain \\
        --set 'materials.board-laminate.conductivity=0.35'
    target/release/foehn flow examples/cases/flow/cht_circuit_board.toml --out output/cht_board_nopad \\
        --set 'interfaces.pad.contact_resistance=0.0'
    target/release/foehn flow examples/cases/flow/cht_orthotropic_patch.toml --out output/cht_patch
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle

from tools.plot.energy_fields import line
from tools.plot.turbulence_fields import read
from tools.plot.colormaps import RAINBOW


def temperature_summary(run):
    with open(os.path.join(run, "result.json")) as f:
        d = json.load(f)
    return next(s for s in d["solves"][0]["scalars"] if s["name"] == "temperature")


def outline(ax):
    ax.add_patch(Rectangle((0.0, 0.0), 0.12, 0.0016, fill=False, ec="w", lw=0.9))
    ax.add_patch(Rectangle((0.04, 0.0016), 0.02, 0.0016, fill=False, ec="w", lw=0.9))


def board(out, run, plain, nopad):
    fig = plt.figure(figsize=(15.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.1])
    pc = None
    for i, (r, title) in enumerate(
        [
            (run, "Board with copper planes: k = 25 W/(m K) along it, 0.35 through it"),
            (plain, "Plain laminate: k = 0.35 W/(m K) in every direction"),
        ]
    ):
        quads, _, d = read(os.path.join(r, "fields.vtkhdf"))
        ax = fig.add_subplot(gs[i])
        pc = PolyCollection(
            quads * 1000.0, array=d["temperature"], cmap="inferno", edgecolors="face",
            linewidths=0.2,
        )
        pc.set_clim(300.0, 425.0)
        ax.add_collection(pc)
        ax.add_patch(Rectangle((0.0, 0.0), 120.0, 1.6, fill=False, ec="c", lw=0.8))
        ax.add_patch(Rectangle((40.0, 1.6), 20.0, 1.6, fill=False, ec="c", lw=0.8))
        ax.set_xlim(0.0, 120.0)
        ax.set_ylim(0.0, 6.0)
        ax.set_aspect(3.0)
        s = temperature_summary(r)
        ax.set_title(f"{title} (y stretched 3 ×); chip at {s['max']:.1f} K", fontsize=11)
        ax.set_ylabel("y (mm)")
        ax.tick_params(labelbottom=False)
        ax.annotate("air, 1 m/s, 300 K →", (3.0, 4.6), color="w", fontsize=10)
    fig.colorbar(pc, ax=fig.axes, label="temperature (K)", shrink=0.9, pad=0.01)
    ax = fig.add_subplot(gs[2])
    for r, label, style in [
        (run, "copper planes, thermal pad R″ = 2·10⁻⁴ m²K/W", ("-", "#d95f02")),
        (nopad, "copper planes, no pad", ("--", "#7570b3")),
        (plain, "plain laminate, thermal pad", ("-", "#1b9e77")),
    ]:
        top = line(os.path.join(r, "line_board_top.csv"))
        chip = temperature_summary(r)["max"]
        ax.plot(top["x"] * 1000.0, top["temperature"], style[0], color=style[1], lw=1.8,
                label=f"{label} (chip {chip:.1f} K)")
    ax.axvspan(40.0, 60.0, color="0.85", zorder=0, label="chip")
    ax.set_xlim(0.0, 120.0)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("top of the board (K)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.savefig(os.path.join(out, "cht_circuit_board.png"), dpi=110)
    plt.close(fig)


def patch(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    t = d["temperature"]
    c = quads.mean(axis=1)
    x, y = c[:, 0], c[:, 1]
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.0), constrained_layout=True)
    ax = axes[0]
    pc = PolyCollection(quads, array=t, cmap=RAINBOW, edgecolors="k", linewidths=0.3)
    ax.add_collection(pc)
    ax.tricontour(x, y, t, levels=np.arange(300.0, 375.0, 4.0), colors="w", linewidths=0.8)
    for xi in (0.4, 0.6):
        ax.axvline(xi, color="c", lw=1.2)
    ax.text(0.02, 0.96, "laminate\nK turned 30°", color="w", va="top", fontsize=10)
    ax.text(0.42, 0.96, "fluid\nat rest", color="w", va="top", fontsize=10)
    ax.text(0.62, 0.96, "block behind\nR″ = 0.01 m²K/W", color="k", va="top", fontsize=10)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Temperature and isotherms on a mesh perturbed by 25 %")
    fig.colorbar(pc, ax=ax, label="temperature (K)", shrink=0.85)
    ax = axes[1]
    err = d["temperature_error"]
    lim = max(np.abs(err).max(), 1e-14)
    pe = PolyCollection(quads, array=err, cmap="RdBu_r", edgecolors="face", linewidths=0.2)
    pe.set_clim(-lim, lim)
    ax.add_collection(pe)
    for xi in (0.4, 0.6):
        ax.axvline(xi, color="k", lw=1.0)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_title(f"Error against the exact solution: at most {lim:.1e} K")
    fig.colorbar(pe, ax=ax, label="T − T_exact (K)", shrink=0.85)
    fig.savefig(os.path.join(out, "cht_orthotropic_patch.png"), dpi=100)
    plt.close(fig)


def main():
    out, run, plain, nopad, patch_run = sys.argv[1:6]
    board(out, run, plain, nopad)
    patch(out, patch_run)


if __name__ == "__main__":
    main()
