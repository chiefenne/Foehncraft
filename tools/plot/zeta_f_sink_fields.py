"""The sink of k near walls of the ζ-f models (PR #141), gallery 74: k on
the TMR flat plate with the papers' sink −ε and with the wall sink
max(k/T_K, 2νk/d²) (k-ε-ζ-f), and k at the stagnation point of the
impinging jet with the sink k/T_K of 6d2 and the wall sink (k-ε-ζ-f₀),
two panels each on one colour scale, from the `fields.vtkhdf` of the runs.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.zeta_f_sink_fields plate OUT.png EPSILON_RUN WALL_RUN
    python3 -m tools.plot.zeta_f_sink_fields jet OUT.png KOLMOGOROV_RUN WALL_RUN

with the runs (`--set output.vtkhdf=true` each):

    foehn flow examples/cases/flow/tmr_plate_zeta_f.toml --out output/pl_eps \\
        --set turbulence.k_sink='"epsilon"' --set solver.max_iterations=3000 \\
        --set verification.require_convergence=false
    foehn flow examples/cases/flow/tmr_plate_zeta_f.toml --out output/pl_wall
    foehn flow examples/cases/flow/impinging_jet.toml --out output/jet_kol \\
        --set turbulence.k_sink='"kolmogorov"' \\
        --set turbulence.constants.c_epsilon1_zeta=0.012 \\
        --set turbulence.constants.c_epsilon1_power=0.5
    foehn flow examples/cases/flow/impinging_jet.toml --out output/jet_wall
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.vtkhdf import internal  # noqa: E402


def centroids(path):
    """Cell centroids (mean of the cell's points) and the cell data."""
    import h5py
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    c = np.add.reduceat(pts[conn], off[:-1], axis=0) / np.diff(off)[:, None]
    return c, data


def plate(out, runs):
    window = (-0.05, 1.0, 0.0, 0.02)
    fig, axes = plt.subplots(len(runs), 1, figsize=(10, 2.3 * len(runs)), sharex=True,
                             constrained_layout=True)
    norm = LogNorm(1e-9, 3e-3)
    levels = np.logspace(-9, np.log10(3e-3), 40)
    for ax, (label, d) in zip(axes, runs):
        c, data = centroids(d + "/fields.vtkhdf")
        x, y, k = c[:, 0], c[:, 1], data["k"]
        keep = (x > window[0] - 0.05) & (x < window[1] + 0.05) & (y < 2 * window[3])
        tri = mtri.Triangulation(x[keep], y[keep] / window[3])
        cs = ax.tricontourf(tri, np.clip(k[keep], 1.01e-9, 2.9e-3), levels=levels, norm=norm,
                            cmap="inferno", extend="both")
        ax.set_xlim(window[0], window[1])
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.01", "0.02"])
        ax.set_ylabel("y / m")
        ax.axvline(0.0, color="w", lw=0.6, ls=":")
        ax.set_title(label, fontsize=10, loc="left")
    axes[-1].set_xlabel("x / m  (the plate starts at x = 0)")
    cb = fig.colorbar(cs, ax=axes, shrink=0.9, pad=0.01, label="k / U∞²")
    cb.set_ticks([1e-9, 1e-7, 1e-5, 1e-3])
    fig.suptitle("TMR flat plate, k-ε-ζ-f: the k of the boundary layer with two sinks of k",
                 fontsize=11)
    fig.savefig(out, dpi=110)
    print(out)


def jet(out, runs):
    window = (0.0, 1.5, 0.0, 0.012)   # r/D, x/D (the plate at x = 0, the pipe exit at 2)
    fig, axes = plt.subplots(len(runs), 1, figsize=(10, 2.5 * len(runs)), sharex=True,
                             constrained_layout=True)
    loaded = []
    for label, d in runs:
        c, data = centroids(d + "/fields.vtkhdf")
        r = np.hypot(c[:, 1], c[:, 2])
        keep = (r < window[1] + 0.2) & (c[:, 0] < 2 * window[3])
        # the axis: mirror the cells next to it so the triangulation reaches r = 0
        near = keep & (r < 0.1)
        r = np.concatenate([r[keep], -r[near]])
        x = np.concatenate([c[keep, 0], c[near, 0]])
        k = np.concatenate([data["k"][keep], data["k"][near]])
        loaded.append((label, r, x / window[3], k))
    norm = LogNorm(1e-5, 2e-2)
    levels = np.logspace(-5, np.log10(2e-2), 41)
    for ax, (label, r, x, k) in zip(axes, loaded):
        tri = mtri.Triangulation(r, x)
        cs = ax.tricontourf(tri, np.clip(k, 1.01e-5, 1.99e-2), levels=levels, norm=norm,
                            cmap=RAINBOW, extend="both")
        ax.set_xlim(window[0], window[1])
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.006", "0.012"])
        ax.set_ylabel("x / D")
        ax.set_title(label, fontsize=10, loc="left")
    axes[-1].set_xlabel("r / D  (the jet axis at r = 0, the plate at x = 0)")
    cb = fig.colorbar(cs, ax=axes, shrink=0.9, pad=0.01, label="k / U_b²")
    cb.set_ticks([1e-5, 1e-4, 1e-3, 1e-2])
    fig.suptitle("Impinging jet at H/D = 2, k-ε-ζ-f₀: k in the wall layer of the plate with two sinks of k",
                 fontsize=11)
    fig.savefig(out, dpi=110)
    print(out)


def main():
    if len(sys.argv) != 5 or sys.argv[1] not in ("plate", "jet"):
        sys.exit(__doc__)
    kind, out, a, b = sys.argv[1:]
    if kind == "plate":
        plate(out, [("sink −ε (the papers; k held at its floor near the leading edge)", a),
                    ("sink max(k/T_K, 2νk/d²) (the default since #141)", b)])
    else:
        jet(out, [("sink k/T_K (6d2): Nu₀ 160", a),
                  ("sink max(k/T_K, 2νk/d²) (since #141): Nu₀ 196", b)])


if __name__ == "__main__":
    main()
