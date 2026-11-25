"""Species transport of stage SF1 (gallery and recipe 278).

species_coaxial_jet.png: the laminar methane jet in a coflow of air on its
axisymmetric wedge, mirrored about the axis: the methane mass fraction,
the density of the mixture and the speed, and the mass fractions and the
density along the axis.
species_binary_diffusion.png: the decay of a cosine mode of CH4 in N2
against the exact solution.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.species_fields OUT_DIR JET_RUN DIFFUSION_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/species_coaxial_jet.toml --out output/jet
    target/release/foehn flow examples/cases/flow/species_binary_diffusion.toml --out output/binary
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from tools.plot.colormaps import RAINBOW


def cells(run):
    """Cell centroids (x, r) and the cell data of a wedge run."""
    with h5py.File(os.path.join(run, "fields.vtkhdf"), "r") as f:
        g = f["VTKHDF/internal"]
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    cen = np.array([pts[conn[a:b]].mean(axis=0) for a, b in zip(off[:-1], off[1:])])
    r = np.hypot(cen[:, 1], cen[:, 2])
    return cen[:, 0], r, data


def line(run, name):
    path = os.path.join(run, f"line_{name}.csv")
    rows = [l for l in open(path) if not l.startswith("#")]
    head = rows[0].strip().split(",")
    vals = np.array([[float(v) for v in l.split(",")] for l in rows[1:]])
    return {h: vals[:, i] for i, h in enumerate(head)}


def jet(out, run):
    x, r, d = cells(run)
    # mirrored about the axis, in mm
    X = np.concatenate([x, x]) * 1000.0
    R = np.concatenate([r, -r]) * 1000.0
    tri = mtri.Triangulation(X, R)
    # drop the long triangles across the axis gap and outside
    xy = np.c_[X, R][tri.triangles]
    span = np.ptp(xy, axis=1).max(axis=1)
    tri.set_mask(span > 12.0)
    fig = plt.figure(figsize=(12.0, 9.6), constrained_layout=True)
    gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1.25])
    panels = [
        ("Y_CH4", 1.0, "viridis", "methane mass fraction Y_CH4 (-)"),
        ("density", 1.0, "magma", "density of the mixture (kg/m³)"),
        ("speed", 1.0, RAINBOW, "speed (m/s)"),
    ]
    for i, (k, s, cmap, label) in enumerate(panels):
        ax = fig.add_subplot(gs[i])
        v = np.concatenate([d[k], d[k]]) * s
        tc = ax.tripcolor(tri, v, shading="gouraud", cmap=cmap)
        fig.colorbar(tc, ax=ax, label=label, shrink=0.9, pad=0.01)
        ax.set_xlim(0.0, X.max())
        ax.set_ylim(-30.0, 30.0)
        ax.set_aspect("equal")
        ax.set_ylabel("r (mm)")
        ax.plot([0, 0], [-5, 5], color="tab:red", lw=5, solid_capstyle="butt")
        if i == 0:
            ax.set_title(
                "Methane jet (red: nozzle, r = 5 mm, 1 m/s) in a coflow of air at 0.2 m/s, "
                "300 K, 1 atm; mirrored about the axis",
                fontsize=11,
            )
    ax = fig.add_subplot(gs[3])
    a = line(run, "axis")
    s_mm = a["x"] * 1000.0
    y_n2 = 1.0 - a["Y_CH4"] - a["Y_O2"]
    ax.plot(s_mm, a["Y_CH4"], color="tab:green", lw=2, label="Y_CH4")
    ax.plot(s_mm, a["Y_O2"], color="tab:blue", lw=2, label="Y_O2")
    ax.plot(s_mm, y_n2, color="tab:gray", lw=2, label="Y_N2 = 1 − Y_CH4 − Y_O2")
    ax.set_xlabel("distance from the nozzle on the axis (mm)")
    ax.set_ylabel("mass fraction (-)")
    ax.set_xlim(0.0, 300.0)
    ax.set_ylim(0.0, 1.02)
    ax.grid(alpha=0.3)
    # the density from the line's mass fractions: rho = p M / (R T)
    m = 1.0 / (a["Y_CH4"] / 0.01604246 + a["Y_O2"] / 0.0319988 + y_n2 / 0.0280134)
    rho = 101325.0 * m / (8.314462618 * a["temperature"])
    ax2 = ax.twinx()
    ax2.plot(s_mm, rho, color="tab:red", lw=1.6, ls="--", label="ρ = p₀ M/(R T)")
    ax2.set_ylabel("density (kg/m³)", color="tab:red")
    ax2.set_ylim(0.6, 1.2)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="center left", fontsize=9)
    ax.set_title("Along the axis: air mixes into the jet and the density rises from that of methane (0.65 kg/m³)", fontsize=11)
    path = os.path.join(out, "species_coaxial_jet.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)


def binary(out, run):
    x, r, d = cells(run)
    row = r < 0.05  # one row of the two cell rows of the slab
    x = x[row]
    d = {k: v[row] for k, v in d.items()}
    order = np.argsort(x)
    xs = x[order]
    y = d["Y_CH4"][order]
    t = 50.0
    exact = 0.5 + 0.4 * np.cos(np.pi * xs) * np.exp(-1e-3 * np.pi**2 * t)
    fig, ax = plt.subplots(figsize=(8.0, 4.4), constrained_layout=True)
    ax.plot(xs, 0.5 + 0.4 * np.cos(np.pi * xs), color="0.6", lw=1.2, ls=":", label="start, t = 0")
    ax.plot(xs, exact, color="k", lw=1.4, label="exact at t = 50 s")
    ax.plot(xs[::3], y[::3], "o", ms=4, color="tab:green", label="Y_CH4, 100 cells, BDF2")
    ax.plot(xs[::3], 1.0 - y[::3], "s", ms=4, color="tab:gray", label="Y_N2 = 1 − Y_CH4 (closure)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("mass fraction (-)")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3)
    err = np.sqrt(np.mean((y - exact) ** 2))
    ax.set_title(f"CH4 diffusing in N2 at rest, D = 1e-3 m²/s; RMS error {err:.1e}", fontsize=11)
    ax.legend(fontsize=9, loc="lower left")
    path = os.path.join(out, "species_binary_diffusion.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)


def main():
    out, jet_run, bin_run = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    jet(out, jet_run)
    binary(out, bin_run)


if __name__ == "__main__":
    main()
