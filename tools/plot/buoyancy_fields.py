"""Gravity and Boussinesq buoyancy of stage 7b1 (gallery and recipe 186).

buoyancy_enclosure_plume.png: a warm strip on the floor of a closed air
enclosure, the temperature with the streamlines of the plume and the speed.
buoyancy_box_at_rest.png: a closed box at rest with a stable stratification
along an oblique gravity on a perturbed mesh, the static pressure and the
speed, which stays at round-off.
buoyancy_pressure_kinds.png: a tall channel with gravity, the reduced
pressure p' the solver works with beside the static pressure p, and both
along the channel axis for an outlet given as p' = 0 and as p = -rho g y.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.buoyancy_fields OUT_DIR PLUME_RUN BOX_RUN REDUCED_RUN STATIC_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/buoyancy_enclosure_plume.toml --out output/plume
    target/release/foehn flow examples/cases/flow/buoyancy_box_at_rest.toml --out output/box
    target/release/foehn flow examples/cases/flow/buoyancy_channel_reduced.toml --out output/chan_reduced
    target/release/foehn flow examples/cases/flow/buoyancy_channel_static.toml --out output/chan_static
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.turbulence_fields import read


def field(run):
    return read(os.path.join(run, "fields.vtkhdf"))


def grid(quads, values, nx, ny):
    """Cell values of a structured hex mesh (x fastest) as (ny, nx)."""
    c = quads.mean(axis=1)
    return c[:, 0].reshape(ny, nx), c[:, 1].reshape(ny, nx), values.reshape(ny, nx)


def plume(out, run):
    quads, _, d = field(run)
    nx, ny = 100, 60
    x, y, ux = grid(quads, d["velocity_x"], nx, ny)
    _, _, uy = grid(quads, d["velocity_y"], nx, ny)
    speed = np.hypot(ux, uy)
    fig, axs = plt.subplots(1, 2, figsize=(14.0, 4.6), constrained_layout=True)
    q = quads * 100.0
    pc = PolyCollection(q, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.2)
    axs[0].add_collection(pc)
    axs[0].streamplot(
        x[0] * 100.0, y[:, 0] * 100.0, ux, uy, color="w", density=1.4, linewidth=0.7, arrowsize=0.8,
    )
    fig.colorbar(pc, ax=axs[0], label="temperature (K)", shrink=0.95, pad=0.01)
    axs[0].set_title(f"Temperature and streamlines; warmest air {d['temperature'].max():.1f} K", fontsize=11)
    pc2 = PolyCollection(q, array=d["speed"] * 100.0, cmap="viridis", edgecolors="face", linewidths=0.2)
    axs[1].add_collection(pc2)
    step = 5
    axs[1].quiver(
        x[::step, ::step] * 100.0, y[::step, ::step] * 100.0, ux[::step, ::step], uy[::step, ::step],
        color="w", scale=0.35, width=0.0025,
    )
    fig.colorbar(pc2, ax=axs[1], label="speed (cm/s)", shrink=0.95, pad=0.01)
    axs[1].set_title(f"Speed and velocity; largest {speed.max() * 100.0:.2f} cm/s", fontsize=11)
    for ax in axs:
        ax.set_xlim(0.0, 10.0)
        ax.set_ylim(0.0, 6.0)
        ax.set_aspect("equal")
        ax.set_xlabel("x (cm)")
        ax.plot([2.9, 4.1], [0.05, 0.05], color="tab:red", lw=4, solid_capstyle="butt")
    axs[0].set_ylabel("y (cm)")
    fig.suptitle(
        "Air in a closed enclosure 10 cm × 6 cm over a floor strip up to 5 K warmer (red), "
        "ceiling 300 K, Ra ≈ 1.1·10⁵",
        fontsize=12,
    )
    fig.savefig(os.path.join(out, "buoyancy_enclosure_plume.png"), dpi=100)
    plt.close(fig)


def box(out, run):
    quads, _, d = field(run)
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 5.0), constrained_layout=True)
    pc = PolyCollection(quads, array=d["static_pressure"], cmap="plasma", edgecolors="k", linewidths=0.25)
    axs[0].add_collection(pc)
    fig.colorbar(pc, ax=axs[0], label="static pressure p (Pa)", shrink=0.95, pad=0.01)
    axs[0].annotate(
        "", xy=(0.62, 0.38), xytext=(0.80, 0.62),
        arrowprops=dict(arrowstyle="-|>", color="w", lw=2.5, mutation_scale=22),
    )
    axs[0].annotate("g", (0.80, 0.50), color="w", fontsize=15, fontweight="bold")
    axs[0].set_title("Static pressure: layered across gravity (37° from the y axis)", fontsize=11)
    s = np.log10(np.maximum(d["speed"], 1e-20))
    pc2 = PolyCollection(quads, array=s, cmap="turbo", edgecolors="k", linewidths=0.25)
    pc2.set_clim(-17.0, -14.0)
    axs[1].add_collection(pc2)
    fig.colorbar(pc2, ax=axs[1], label="log₁₀ speed (m/s)", shrink=0.95, pad=0.01)
    axs[1].set_title(f"Speed after 1000 iterations: at most {d['speed'].max():.1e} m/s", fontsize=11)
    for ax in axs:
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
    axs[0].set_ylabel("y (m)")
    fig.suptitle(
        "A closed box with a stable stratification (14 K, warm above) on a mesh perturbed by 20 %",
        fontsize=12,
    )
    fig.savefig(os.path.join(out, "buoyancy_box_at_rest.png"), dpi=100)
    plt.close(fig)


def kinds(out, reduced, static):
    qa, _, a = field(reduced)
    qb, _, b = field(static)
    fig = plt.figure(figsize=(13.0, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.4, 1.4])
    panels = [
        (a["pressure"], "reduced pressure p′\n(what the solver works with)", "viridis"),
        (a["static_pressure"] / 1000.0, "static pressure p\n= p′ − ρ g y", "magma"),
    ]
    units = ["p′ (Pa)", "p (kPa)"]
    for i, (v, title, cmap) in enumerate(panels):
        ax = fig.add_subplot(gs[i])
        pc = PolyCollection(qa * [10.0, 1.0], array=v, cmap=cmap, edgecolors="face", linewidths=0.2)
        ax.add_collection(pc)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 2.0)
        ax.set_aspect(1.0)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x (dm, 10 × stretched)")
        if i == 0:
            ax.set_ylabel("y (m)")
        else:
            ax.tick_params(labelleft=False)
        fig.colorbar(pc, ax=ax, label=units[i], shrink=0.9, pad=0.02, location="bottom")
    lines = [
        (a, qa, "outlet p′ = 0 (\"reduced\", default)", "-"),
        (b, qb, "outlet p = −19 620 Pa (\"static\")", "--"),
    ]
    for j, (key, scale, label) in enumerate(
        [("pressure", 1.0, "p′ on the axis (Pa)"), ("static_pressure", 1e-3, "p on the axis (kPa)")]
    ):
        ax = fig.add_subplot(gs[2 + j])
        for d, q, lab, ls in lines:
            c = q.mean(axis=1)
            m = np.abs(c[:, 0] - 0.045) < 1e-9
            o = np.argsort(c[m, 1])
            ax.plot(d[key][m][o] * scale, c[m, 1][o], ls, lw=2.2, label=lab)
        ax.set_xlabel(label)
        ax.set_ylim(0.0, 2.0)
        ax.grid(alpha=0.3)
        ax.tick_params(labelleft=False)
        if j == 0:
            ax.legend(fontsize=9, loc="upper right")
            ax.set_title("p′: friction only", fontsize=11)
        else:
            ax.set_title("p: friction + hydrostatics", fontsize=11)
    fig.suptitle(
        "Water-like liquid rising at 0.05 m/s in a channel 0.1 m × 2 m, gravity downwards",
        fontsize=12,
    )
    fig.savefig(os.path.join(out, "buoyancy_pressure_kinds.png"), dpi=100)
    plt.close(fig)


def main():
    out, plume_run, box_run, reduced, static = sys.argv[1:6]
    os.makedirs(out, exist_ok=True)
    plume(out, plume_run)
    box(out, box_run)
    kinds(out, reduced, static)


if __name__ == "__main__":
    main()
