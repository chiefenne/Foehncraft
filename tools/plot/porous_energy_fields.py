"""Two temperatures in porous regions, stage PM2 of the porous plan (D271;
gallery and recipe of the porous regions).

porous_schumann.png: gas and matrix temperature along Schumann's bed at
three times against the exact solution (Schumann 1929).
porous_storage_charge.png: a honeycomb store charged with hot air: the
matrix temperature at three times, the gas temperature with the
streamlines, and the temperatures on the axis.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.porous_energy_fields OUT_DIR SCHUMANN_RUN STORAGE_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/porous_schumann.toml --out output/schumann
    target/release/foehn flow examples/cases/flow/porous_storage_charge.toml --out output/storage
"""
import glob
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

# porous_schumann.toml
RHO_C_U, C_S, HA, EPS, U = 1.2 * 1005.0 * 0.1, 0.6 * 2300.0 * 900.0, 603.0, 0.4, 0.1


def need(d, name, path):
    if name not in d:
        sys.exit(f"{path}: no field `{name}`; the file has: {', '.join(sorted(d))}")
    return d[name]


def i0e(z):
    """e^-z I0(z) by its series (z up to about 40)."""
    q = 0.25 * z * z
    term = np.ones_like(z)
    s = np.ones_like(z)
    for k in range(1, 400):
        term = term * q / (k * k)
        s = s + term
    return s * np.exp(-z)


def schumann_j(x, y, n=2000):
    if x <= 0.0:
        return 1.0
    s = np.linspace(0.0, x, n + 1)
    z = 2.0 * np.sqrt(y * s)
    f = np.exp(-(y + s) + z) * i0e(z)
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return 1.0 - (x / n / 3.0) * np.dot(w, f)


def schumann(xi, eta):
    if eta < 0.0:
        return 0.0, 0.0
    return schumann_j(xi, eta), 1.0 - schumann_j(eta, xi)


def fields_in_order(run):
    """The fields files of a run in time order: fields_<step> then the final."""
    steps = sorted(glob.glob(os.path.join(run, "fields_*.vtkhdf")),
                   key=lambda p: int(os.path.basename(p)[7:-7]))
    return steps + [os.path.join(run, "fields.vtkhdf")]


def schumann_figure(out, run):
    fig, ax = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    xs = np.linspace(0.0, 1.0, 201)
    colors = ["#2b6cb0", "#d69e2e", "#c53030"]
    for path, t, c in zip(fields_in_order(run), [1000.0, 2000.0, 3000.0], colors):
        quads, _, d = read(path)
        x = quads.mean(axis=1)[:, 0]
        o = np.argsort(x)
        tg = need(d, "temperature", path)[o]
        ts = need(d, "matrix_temperature", path)[o]
        ex = np.array([schumann(HA * v / RHO_C_U, HA * (t - EPS * v / U) / C_S) for v in xs])
        ax.plot(xs, 300.0 + 100.0 * ex[:, 0], "-", color=c, lw=1.2)
        ax.plot(xs, 300.0 + 100.0 * ex[:, 1], "--", color=c, lw=1.2)
        ax.plot(x[o][::3], tg[::3], "o", color=c, ms=4, mfc="none", label=f"gas, {t:.0f} s")
        ax.plot(x[o][::3], ts[::3], "s", color=c, ms=3.5, label=f"matrix, {t:.0f} s")
    ax.plot([], [], "k-", lw=1.2, label="exact gas (Schumann 1929)")
    ax.plot([], [], "k--", lw=1.2, label="exact matrix")
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("temperature (K)")
    ax.legend(fontsize=8.5, ncol=2)
    ax.set_title("Air at 400 K into a ceramic bed at 300 K, 0.1 m/s, no conduction: "
                 "every third of 100 cells against the exact solution", fontsize=10.5)
    fig.savefig(os.path.join(out, "porous_schumann.png"), dpi=110)
    plt.close(fig)


def field_panel(ax, quads, values, title, lo, hi, stack_only=False):
    if stack_only:
        # the matrix is only in the stack; the plenums are grey
        y = quads.mean(axis=1)[:, 1]
        values = np.where((y > 0.15) & (y < 0.85), values, np.nan)
    cmap = matplotlib.colormaps[RAINBOW].copy()
    cmap.set_bad("0.85")
    pc = PolyCollection(quads * 100.0, array=values, cmap=cmap, edgecolors="face",
                        linewidths=0.2)
    pc.set_clim(lo, hi)
    ax.add_collection(pc)
    ax.add_patch(Rectangle((0.0, 15.0), 40.0, 70.0, fill=False, ec="k", lw=1.0, ls="--"))
    ax.set_xlim(0.0, 40.0)
    ax.set_ylim(0.0, 100.0)
    ax.set_aspect("equal")
    ax.set_xticks([0, 20, 40])
    ax.set_title(title, fontsize=10)
    return pc


def storage_figure(out, run):
    paths = fields_in_order(run)
    times = [1500.0, 3000.0, 6000.0]
    fig = plt.figure(figsize=(13.0, 6.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 1, 0.08, 2.1])
    lo, hi = 293.0, 573.0
    pc = None
    for i, (path, t) in enumerate(zip(paths, times)):
        quads, _, d = read(path)
        ax = fig.add_subplot(gs[0, i])
        pc = field_panel(ax, quads, need(d, "matrix_temperature", path), f"matrix, {t:.0f} s",
                         lo, hi, stack_only=True)
        if i == 0:
            ax.set_ylabel("y (cm)")
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("x (cm)")
    # the gas at the middle time, with its streamlines
    quads, _, d = read(paths[1])
    ax = fig.add_subplot(gs[0, 3])
    field_panel(ax, quads, need(d, "temperature", paths[1]), f"gas, {times[1]:.0f} s", lo, hi)
    c = quads.mean(axis=1)
    order = np.lexsort((c[:, 0], c[:, 1]))
    nx, ny = 40, 100
    gx = c[order, 0].reshape(ny, nx)[0] * 100.0
    gy = c[order, 1].reshape(ny, nx)[:, 0] * 100.0
    ux = need(d, "velocity_x", paths[1])[order].reshape(ny, nx)
    uy = need(d, "velocity_y", paths[1])[order].reshape(ny, nx)
    ax.streamplot(gx, gy, ux, uy, color="w", density=1.0, linewidth=0.6, arrowsize=0.7)
    ax.set_yticklabels([])
    ax.set_xlabel("x (cm)")
    fig.colorbar(pc, cax=fig.add_subplot(gs[0, 4]), label="temperature (K)")
    # the axis at the three times
    ax = fig.add_subplot(gs[0, 5])
    colors = ["#2b6cb0", "#d69e2e", "#c53030"]
    for path, t, col in zip(paths, times, colors):
        quads, _, d = read(path)
        c = quads.mean(axis=1)
        sel = np.abs(c[:, 0] - 0.205) < 0.006
        o = np.argsort(-c[sel, 1])
        y = c[sel, 1][o] * 100.0
        ax.plot(need(d, "temperature", path)[sel][o], y, "-", color=col, lw=1.3,
                label=f"gas, {t:.0f} s")
        tm = np.where((y > 15.0) & (y < 85.0), need(d, "matrix_temperature", path)[sel][o],
                      np.nan)
        ax.plot(tm, y, "--", color=col, lw=1.3,
                label=f"matrix, {t:.0f} s")
    ax.axhspan(15.0, 85.0, color="0.92", zorder=0)
    ax.set_ylim(0.0, 100.0)
    ax.set_xlim(lo - 10.0, hi + 10.0)
    ax.set_xlabel("temperature (K) at x = 20 cm")
    ax.set_ylabel("y (cm)")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("on the axis; the stack is shaded", fontsize=10)
    fig.suptitle("Charging a honeycomb store: air at 573 K from the top right into ceramic "
                 "bricks at 293 K (6.9 mm channels, ε = 0.26)", fontsize=12)
    fig.savefig(os.path.join(out, "porous_storage_charge.png"), dpi=110)
    plt.close(fig)


def main():
    out, schumann_run, storage_run = sys.argv[1:4]
    schumann_figure(out, schumann_run)
    storage_figure(out, storage_run)


if __name__ == "__main__":
    main()
