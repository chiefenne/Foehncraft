"""Turbulent conjugate heat transfer of stage 7d5 (gallery entry of the PR).

cht_heat_sink.png: a finned aluminium heat sink in a turbulent air channel,
the temperature with the streamlines (top) and the speed (bottom).
cht_turbulent_channel.png: the turbulent channel over a conducting wall
(check V7.36), the temperature field (left, height stretched) and the
profile across the wall and the fluid for k-ε and SST (right).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_turbulent OUT_DIR SINK_RUN CHANNEL_KEPS_RUN CHANNEL_SST_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/cht_heat_sink_turbulent.toml --out output/cht_sink
    target/release/foehn flow examples/cases/flow/cht_turbulent_channel.toml --out output/cht_channel_keps
    target/release/foehn flow examples/cases/flow/cht_turbulent_channel.toml --out output/cht_channel_sst \\
        --set turbulence.model=k-omega-sst --set turbulence.wall_treatment=automatic
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.energy_fields import grid, line
from tools.plot.turbulence_fields import read
from matplotlib.collections import PolyCollection


def panel(ax, quads, values, cmap, label, fig, xlim, ylim, vmin=None, vmax=None):
    c = quads.mean(axis=1)
    keep = (c[:, 0] > xlim[0] - 2e-3) & (c[:, 0] < xlim[1] + 2e-3)
    pc = PolyCollection(
        quads[keep], array=values[keep], cmap=cmap, edgecolors="face", linewidths=0.2
    )
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fig.colorbar(pc, ax=ax, label=label, shrink=0.9, pad=0.015)
    return pc


def outline(ax, xs, ys, solid, color):
    ax.contour(xs, ys, solid, levels=[0.5], colors=color, linewidths=1.0)


def sink(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    xs, ys, solid = grid(quads, (d["zone"] > 0).astype(float))
    _, _, ux = grid(quads, d["velocity_x"])
    _, _, uy = grid(quads, d["velocity_y"])
    xlim, ylim = (0.05, 0.26), (0.0, 0.05)
    fig, axes = plt.subplots(2, 1, figsize=(12.0, 6.6), constrained_layout=True)
    ax = axes[0]
    panel(ax, quads, d["temperature"], "inferno", "temperature (K)", fig, xlim, ylim)
    sel = (xs > xlim[0] - 0.01) & (xs < xlim[1] + 0.01)
    ax.streamplot(
        xs[sel], ys, ux[:, sel], uy[:, sel], color="w", linewidth=0.45, density=1.6,
        arrowsize=0.5, broken_streamlines=True,
    )
    outline(ax, xs, ys, solid, "c")
    ax.set_aspect("equal")
    ax.set_title(
        "Finned aluminium heat sink in a turbulent air channel (5 m/s, k-ω SST): "
        "2000 W/m² from the chip below the base"
    )
    ax.set_ylabel("y (m)")
    ax = axes[1]
    panel(ax, quads, d["speed"], "viridis", "speed (m/s)", fig, xlim, ylim)
    ax.streamplot(
        xs[sel], ys, ux[:, sel], uy[:, sel], color="w", linewidth=0.35, density=1.6,
        arrowsize=0.5,
    )
    outline(ax, xs, ys, solid, "#ff7f0e")
    ax.set_aspect("equal")
    ax.set_title("Speed: the flow skims over the fins, slow vortices sit between them")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.savefig(os.path.join(out, "cht_heat_sink.png"), dpi=110)


def channel(out, keps, sst):
    quads, _, d = read(os.path.join(keps, "fields.vtkhdf"))
    fig, axes = plt.subplots(
        1, 2, figsize=(14.5, 4.6), constrained_layout=True, gridspec_kw={"width_ratios": [1.7, 1]}
    )
    ax = axes[0]
    panel(ax, quads, d["temperature"], "inferno", "temperature (K)", fig, (0.0, 20.0), (0.0, 0.45))
    ax.axhline(0.1, color="c", lw=1.0, ls="--")
    ax.axvline(15.0, color="w", lw=0.8, ls=":")
    ax.text(19.8, 0.05, "solid wall, t/k_s = 0.5 m²K/W", color="k", ha="right", va="center",
            fontsize=9)
    ax.text(0.3, 0.43, "fluid, Re = 2·10⁵, Pr = 0.71, 300 K in (the channel goes on to y = 1.1 m)",
            color="w", va="top", fontsize=9)
    ax.set_title("Turbulent channel over a conducting wall at 350 K (k-ε, lowest 0.45 m, stretched)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax = axes[1]
    for run, c, lab in ((keps, "#1b9e77", "k-ε, scalable"), (sst, "#d95f02", "k-ω SST, automatic")):
        p = os.path.join(run, "line_across.csv")
        if not os.path.exists(p):
            continue
        y, t = cells(line(p))
        # no line across the interface: the jump there is the wall function
        s = y < 0.1
        ax.plot(t[s], y[s], "o-", color=c, lw=1.4, ms=4, label=lab)
        ax.plot(t[~s], y[~s], "o-", color=c, lw=1.4, ms=4)
    ax.axhline(0.1, color="0.4", lw=0.8, ls="--")
    ax.text(301, 0.11, "interface", color="0.3", fontsize=9, va="bottom")
    ax.set_ylim(0.0, 0.4)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("y (m)")
    ax.set_title("Across the wall and the fluid at x = 15 m")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out, "cht_turbulent_channel.png"), dpi=110)


def cells(v):
    """The cell values along a line (one per run of equal values) at the
    middle of each run: the cell centres of a line across a hex mesh."""
    y, t = v["y"], v["temperature"]
    out_y, out_t, start = [], [], 0
    for i in range(1, len(t) + 1):
        if i == len(t) or t[i] != t[start]:
            out_y.append(0.5 * (y[start] + y[i - 1]))
            out_t.append(t[start])
            start = i
    return np.array(out_y), np.array(out_t)


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    out, sink_run, keps, sst = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    sink(out, sink_run)
    channel(out, keps, sst)


if __name__ == "__main__":
    main()
