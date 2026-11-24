"""Low-Mach variable density of stage 7c (gallery and recipe 195).

vd_cavity_large_dt.png: air in a closed square cavity with walls at 960 K
and 240 K, the temperature with the streamlines, the density with the
velocity, and the vertical velocity across the middle of the cavity.
vd_checks.png: the heated duct along its axis against the exact outlet
values, the closed box heated through one wall, and p₀(t) of both closed
boxes against the exact line of the box heated uniformly.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.variable_density OUT_DIR CAVITY_RUN DUCT_RUN BOX_RUN SIDE_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/vd_cavity_large_dt.toml --out output/cavity \\
        --set 'mesh.generate.cells=[128,128,1]' --set 'verification.levels=[0]'
    target/release/foehn flow examples/cases/flow/vd_heated_duct.toml --out output/duct
    target/release/foehn flow examples/cases/flow/vd_closed_box_heating.toml --out output/box
    target/release/foehn flow examples/cases/flow/vd_closed_box_side_heated.toml --out output/side
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

R_GAS = 8.314462618 / 0.0289647
CP = 1005.0
P_OP = 101_325.0


def field(run):
    return read(os.path.join(run, "fields.vtkhdf"))


def grid(quads, values, nx, ny):
    """Cell values of a structured hex mesh (x fastest) as (ny, nx)."""
    c = quads.mean(axis=1)
    return c[:, 0].reshape(ny, nx), c[:, 1].reshape(ny, nx), values.reshape(ny, nx)


def uniform(x, y, v, n):
    """A tensor-product field resampled on an n × n uniform grid (for the streamlines)."""
    xs = np.linspace(x[0].min(), x[0].max(), n)
    ys = np.linspace(y[:, 0].min(), y[:, 0].max(), n)
    rows = np.array([np.interp(xs, x[0], r) for r in v])
    cols = np.array([np.interp(ys, y[:, 0], c) for c in rows.T]).T
    return xs, ys, cols


def history(run):
    with open(os.path.join(run, "history.csv")) as f:
        rows = [r for r in csv.reader(line for line in f if not line.startswith("#"))]
    head, body = rows[0], rows[1:]
    return {k: np.array([float(r[i]) for r in body]) for i, k in enumerate(head)}


def cavity(out, run):
    quads, _, d = field(run)
    n = round(float(np.sqrt(len(d["temperature"]))))
    x, y, ux = grid(quads, d["velocity_x"], n, n)
    _, _, uy = grid(quads, d["velocity_y"], n, n)
    _, _, t = grid(quads, d["temperature"], n, n)
    q = quads * 100.0
    fig = plt.figure(figsize=(15.0, 4.9), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05])
    ax = fig.add_subplot(gs[0])
    pc = PolyCollection(q, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.1)
    ax.add_collection(pc)
    xs, ys, uxs = uniform(x, y, ux, 160)
    _, _, uys = uniform(x, y, uy, 160)
    ax.streamplot(xs * 100.0, ys * 100.0, uxs, uys, color="w", density=1.3, linewidth=0.7, arrowsize=0.8)
    fig.colorbar(pc, ax=ax, label="temperature (K)", shrink=0.95, pad=0.01)
    ax.set_title("Temperature and streamlines", fontsize=11)
    ax2 = fig.add_subplot(gs[1])
    pc2 = PolyCollection(q, array=d["density"], cmap=RAINBOW, edgecolors="face", linewidths=0.1)
    ax2.add_collection(pc2)
    step = 8
    ax2.quiver(
        x[::step, ::step] * 100.0, y[::step, ::step] * 100.0, ux[::step, ::step], uy[::step, ::step],
        color="w", scale=2.2, width=0.004,
    )
    fig.colorbar(pc2, ax=ax2, label="density ρ(T, p₀) (kg/m³)", shrink=0.95, pad=0.01)
    ax2.set_title(
        f"Density {d['density'].min():.2f} to {d['density'].max():.2f} kg/m³ and velocity", fontsize=11
    )
    for a in (ax, ax2):
        a.set_xlim(0.0, 3.2)
        a.set_ylim(0.0, 3.2)
        a.set_aspect("equal")
        a.set_xlabel("x (cm)")
        a.plot([0.02, 0.02], [0.0, 3.2], color="tab:red", lw=5, solid_capstyle="butt")
        a.plot([3.18, 3.18], [0.0, 3.2], color="tab:blue", lw=5, solid_capstyle="butt")
    ax.set_ylabel("y (cm)")
    # the vertical velocity across the middle: the two boundary layers
    ax3 = fig.add_subplot(gs[2])
    j = n // 2
    mid = 0.5 * (uy[j - 1] + uy[j]) if n % 2 == 0 else uy[j]
    tmid = 0.5 * (t[j - 1] + t[j]) if n % 2 == 0 else t[j]
    ax3.plot(x[0] * 100.0, mid * 100.0, color="tab:purple", lw=2.0, label="vertical velocity")
    ax3.axhline(0.0, color="0.6", lw=0.8)
    ax3.set_xlabel("x (cm)")
    ax3.set_ylabel("vertical velocity at mid-height (cm/s)", color="tab:purple")
    ax3.set_xlim(0.0, 3.2)
    tw = ax3.twinx()
    tw.plot(x[0] * 100.0, tmid, color="tab:orange", lw=2.0, ls="--", label="temperature")
    tw.set_ylabel("temperature at mid-height (K)", color="tab:orange")
    ax3.set_title("Across the middle: a broad hot layer, a thin cold one", fontsize=11)
    fig.suptitle(
        "Air in a closed cavity of 3.2 cm, hot wall 960 K (red), cold wall 240 K (blue): "
        "variable density, p₀ = 94 975 Pa, 128 × 128 cells",
        fontsize=12,
    )
    fig.savefig(os.path.join(out, "vd_cavity_large_dt.png"), dpi=100)
    plt.close(fig)


def checks(out, duct, box, side):
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)
    # the heated duct along its axis
    quads, _, d = field(duct)
    x, _, t = grid(quads, d["temperature"], 100, 4)
    _, _, rho = grid(quads, d["density"], 100, 4)
    _, _, u = grid(quads, d["velocity_x"], 100, 4)
    ax = axs[0]
    ax.axvspan(0.4, 0.6, color="tab:red", alpha=0.12, label="heat source")
    ax.plot(x[0], t.mean(axis=0) / 300.0, color="tab:orange", lw=2.2, label="T / T_in")
    ax.plot(x[0], u.mean(axis=0) / 1.0, color="tab:purple", lw=2.2, ls="--", label="u / u_in")
    ax.plot(x[0], rho.mean(axis=0) / (P_OP / (R_GAS * 300.0)), color="tab:green", lw=2.2, label="ρ / ρ_in")
    t_out = 300.0 + 5.9131e5 * 0.2 / (P_OP / (R_GAS * 300.0) * CP)
    for v, c in [(t_out / 300.0, "tab:orange"), (300.0 / t_out, "tab:green")]:
        ax.axhline(v, color=c, lw=0.8, ls=":")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("ratio to the inlet value")
    ax.set_title("Heated duct: T, u up by 4/3, ρ down by 3/4 (dotted: exact)", fontsize=11)
    ax.legend(fontsize=9, loc="center left")
    # the closed box heated through one wall
    quads, _, d = field(side)
    x, y, ux = grid(quads, d["velocity_x"], 32, 32)
    _, _, uy = grid(quads, d["velocity_y"], 32, 32)
    ax = axs[1]
    pc = PolyCollection(quads * 100.0, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.3)
    ax.add_collection(pc)
    s = 3
    ax.quiver(
        x[::s, ::s] * 100.0, y[::s, ::s] * 100.0, ux[::s, ::s], uy[::s, ::s],
        color="w", scale=2.0e-3, width=0.004,
    )
    fig.colorbar(pc, ax=ax, label="temperature (K)", shrink=0.95, pad=0.01)
    ax.plot([0.01, 0.01], [0.0, 2.0], color="tab:red", lw=5, solid_capstyle="butt")
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 2.0)
    ax.set_aspect("equal")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_title(
        f"Closed box, left wall 400 K, at 2 s: flow towards the wall, {d['speed'].max() * 1e3:.2f} mm/s",
        fontsize=11,
    )
    # p₀(t) of both boxes
    ax = axs[2]
    hb, hs = history(box), history(side)
    rho0 = P_OP / (R_GAS * 300.0)
    tt = np.linspace(0.0, 5.0, 50)
    ax.plot(tt, P_OP * (1.0 + 1.0e4 * tt / (rho0 * (CP - R_GAS)) / 300.0) / 1000.0, color="0.4", lw=1.2,
            label="exact, uniform heating (c_v)")
    ax.plot(np.r_[0.0, hb["time"]], np.r_[P_OP, hb["p0"]] / 1000.0, "o", color="tab:blue", ms=5,
            label="box heated uniformly, 10⁴ W/m³")
    ax.plot(np.r_[0.0, hs["time"]], np.r_[P_OP, hs["p0"]] / 1000.0, color="tab:red", lw=2.2,
            label="box heated through the left wall")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("thermodynamic pressure p₀ (kPa)")
    ax.set_title("p₀ of the closed boxes; the mass stays to round-off", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    fig.savefig(os.path.join(out, "vd_checks.png"), dpi=100)
    plt.close(fig)


def main():
    out, cav, duct, box, side = sys.argv[1:6]
    os.makedirs(out, exist_ok=True)
    cavity(out, cav)
    checks(out, duct, box, side)


if __name__ == "__main__":
    main()
