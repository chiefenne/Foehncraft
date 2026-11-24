"""Porous regions of the porous plan, stage PM1 (D270; gallery and recipe of
the porous regions).

porous_honeycomb_channel.png: an air channel with a ceramic honeycomb block,
its channels along the flow, beside a block of the same resistance in every
direction: the speed with the streamlines and the pressure on the axis.
porous_block_convergence.png: the pressure along a duct with a porous block
and the error of the pressure drop against the cell size.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.porous_fields OUT_DIR ALONG_RUN ISO_RUN BLOCK_RUN_50 BLOCK_RUN_100 BLOCK_RUN_200

with the output directories of

    target/release/foehn flow examples/cases/flow/porous_honeycomb_channel.toml --out output/along
    target/release/foehn flow examples/cases/flow/porous_honeycomb_channel.toml --set regions.brick.porous.honeycomb.blocking=1.0 --out output/iso
    target/release/foehn flow examples/cases/flow/porous_darcy_block.toml --set 'mesh.generate.cells=[50, 4, 1]' --out output/block50
    target/release/foehn flow examples/cases/flow/porous_darcy_block.toml --out output/block100
    target/release/foehn flow examples/cases/flow/porous_darcy_block.toml --set 'mesh.generate.cells=[200, 4, 1]' --out output/block200
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

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read

# the brick of porous_honeycomb_channel.toml (cm)
BRICK = (22.5, 4.0, 15.0, 12.0)
MU, RHO = 1.8e-5, 1.2


def grid(quads, values, nx, ny):
    c = quads.mean(axis=1)
    order = np.lexsort((c[:, 0], c[:, 1]))
    c = c[order]
    return c[:, 0].reshape(ny, nx), c[:, 1].reshape(ny, nx), values[order].reshape(ny, nx)


def brick(ax, across):
    x, y, w, h = BRICK
    ax.add_patch(Rectangle((x, y), w, h, fill=False, ec="k", lw=1.2, ls="--"))


def speed_panel(fig, ax, run, title):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    nx, ny = 240, 80
    x, y, ux = grid(quads, d["velocity_x"], nx, ny)
    _, _, uy = grid(quads, d["velocity_y"], nx, ny)
    q = quads * 100.0
    sp = np.hypot(d["velocity_x"], d["velocity_y"]) * 100.0
    pc = PolyCollection(q, array=sp, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    pc.set_clim(0.0, 16.0)
    ax.add_collection(pc)
    ax.streamplot(x[0] * 100.0, y[:, 0] * 100.0, ux, uy, color="w", density=1.5,
                  linewidth=0.7, arrowsize=0.8)
    fig.colorbar(pc, ax=ax, label="superficial speed (cm/s)", shrink=0.95, pad=0.01)
    brick(ax, False)
    ax.set_xlim(0.0, 60.0)
    ax.set_ylim(0.0, 20.0)
    ax.set_aspect("equal")
    ax.set_ylabel("y (cm)")
    ax.set_title(title, fontsize=11)
    return x[0], y[:, 0], grid(quads, d["pressure"], nx, ny)[2], ux


def channels(out, along, iso):
    fig, axs = plt.subplots(3, 1, figsize=(11.0, 11.0), constrained_layout=True,
                            gridspec_kw={"height_ratios": [1, 1, 0.8]})
    xa, ya, pa, ua = speed_panel(fig, axs[0], along,
                                 "Honeycomb: channels along x, 10⁴ times the resistance across them")
    xi, yi, pi, ui = speed_panel(fig, axs[1], iso,
                                 "The same resistance in every direction (isotropic porous block)")
    j = len(ya) // 2
    axs[2].plot(xa * 100.0, pa[j], color="tab:blue", label="honeycomb, pressure on the axis")
    axs[2].plot(xi * 100.0, pi[j], color="tab:orange", label="isotropic, pressure on the axis")
    axs[2].axvspan(22.5, 37.5, color="0.88", zorder=0)
    axs[2].set_xlim(0.0, 60.0)
    axs[2].set_xlabel("x (cm)")
    axs[2].set_ylabel("pressure (Pa)")
    k = np.argmin(np.abs(xa - 0.45))
    axs[2].set_title(
        f"Pressure along the axis; the block is shaded. Speed 7.5 cm behind the block on the axis: "
        f"{ua[j, k] * 100:.1f} cm/s (honeycomb), {ui[j, k] * 100:.1f} cm/s (isotropic)",
        fontsize=10,
    )
    axs[2].legend(fontsize=9)
    fig.suptitle(
        "Air at 10 cm/s through a channel 20 cm high with a ceramic honeycomb block (dashed, "
        "6.9 mm channels, ε = 0.26)",
        fontsize=12,
    )
    fig.savefig(os.path.join(out, "porous_honeycomb_channel.png"), dpi=110)
    plt.close(fig)


def block(out, runs):
    exact = (MU * 1e6 + 0.5 * RHO * 10.0) * 0.2
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.4), constrained_layout=True)
    hs, errs = [], []
    for run, n, c in zip(runs, [50, 100, 200], ["tab:blue", "tab:orange", "tab:green"]):
        quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
        xc = quads.mean(axis=1)[:, 0]
        sel = quads.mean(axis=1)[:, 1] < 0.1 / 4
        o = np.argsort(xc[sel])
        axs[0].plot(xc[sel][o], d["pressure"][sel][o], ".-", color=c, ms=3, lw=0.8, label=f"{n} cells")
        r = json.load(open(os.path.join(run, "result.json")))["solves"][0]
        p_in = next(b["mean_pressure"] for b in r["boundaries"] if b["name"] == "xmin")
        hs.append(1.0 / n)
        errs.append(abs(p_in - exact))
    xs = np.array([0.0, 0.4, 0.6, 1.0])
    axs[0].plot(xs, [exact, exact, 0.0, 0.0], "k--", lw=1.0, label="exact")
    axs[0].axvspan(0.4, 0.6, color="0.85", zorder=0)
    axs[0].set_xlabel("x (m)")
    axs[0].set_ylabel("pressure (Pa)")
    axs[0].set_title("Pressure along the duct; the block is shaded", fontsize=11)
    axs[0].legend(fontsize=9)
    hs, errs = np.array(hs), np.array(errs)
    axs[1].loglog(hs, errs / exact, "o-", color="tab:red", label="error of the drop")
    axs[1].loglog(hs, errs[0] / exact * (hs / hs[0]) ** 2, "k:", label="second order")
    axs[1].set_xlabel("cell size (m)")
    axs[1].set_ylabel("relative error of Δp")
    axs[1].set_title(f"Δp = {exact:.1f} Pa exact; {errs[1] / exact * 100:.2f} % off on 100 cells", fontsize=11)
    axs[1].legend(fontsize=9)
    fig.suptitle("A porous block in a duct: Darcy–Forchheimer drop against the exact value", fontsize=12)
    fig.savefig(os.path.join(out, "porous_block_convergence.png"), dpi=110)
    plt.close(fig)


def main():
    out, along_run, iso_run, *block_runs = sys.argv[1:]
    channels(out, along_run, iso_run)
    block(out, block_runs)


if __name__ == "__main__":
    main()
