"""The steady residual of the adjoint (L17 stage A0), gallery 271: left
the speed of the converged lid-driven cavity on perturbed polygons with
streamlines, middle the momentum residual |R_u| per cell of the start
fields after one outer iteration, right the three norms of the one
residual R(q) along the SIMPLEC run and the coupled run of the same case.
From the CSV files of the example ``adjoint_a0``:

    cargo run --release -p foehncraft-flow --example adjoint_a0 -- \\
        examples/cases/flow/adjoint_residual_cavity.toml output/adjoint_a0
    pip install matplotlib numpy
    python3 -m tools.plot.adjoint_a0 docs/figures/gallery/adjoint_a0_residual.png output/adjoint_a0
"""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from tools.plot.colormaps import RAINBOW


def read(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    out, root = Path(sys.argv[1]), Path(sys.argv[2])
    cells = read(root / "cells.csv")
    hs = read(root / "history_simplec.csv")
    hc = read(root / "history_coupled.csv")
    tri = mtri.Triangulation(cells["x"], cells["y"])
    speed = np.hypot(cells["ux"], cells["uy"])

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), gridspec_kw={"width_ratios": [1, 1, 1.25]})
    ax = axes[0]
    t = ax.tricontourf(tri, speed, levels=40, cmap=RAINBOW)
    gx, gy = np.meshgrid(np.linspace(0.01, 0.99, 80), np.linspace(0.01, 0.99, 80))
    interp = lambda v: mtri.LinearTriInterpolator(tri, v)(gx, gy)
    ax.streamplot(gx, gy, interp(cells["ux"]), interp(cells["uy"]), color="white",
                  linewidth=0.6, density=1.1, arrowsize=0.6)
    fig.colorbar(t, ax=ax, shrink=0.85, label="speed |u| (m/s)")
    ax.set_title("converged speed, Re 100, 576 polygons", fontsize=10)

    ax = axes[1]
    r = np.log10(np.maximum(cells["r_start"], 1e-30))
    t = ax.tricontourf(tri, r, levels=40, cmap=RAINBOW)
    fig.colorbar(t, ax=ax, shrink=0.85, label="log10 |R_u| per cell (N)")
    ax.set_title("|R_u| per cell after one iteration", fontsize=10)
    for a in axes[:2]:
        a.set_aspect("equal")
        a.set_xlim(0, 1)
        a.set_ylim(0, 1)
        a.set_xlabel("x (m)")
        a.set_ylabel("y (m)")

    ax = axes[2]
    colours = {"momentum": "#1f77b4", "continuity": "#2ca02c", "mass_flux": "#d62728"}
    labels = {"momentum": "R_u", "continuity": "R_c", "mass_flux": "R_ṁ"}
    for key, c in colours.items():
        ax.semilogy(hs["iteration"], hs[key], color=c, lw=1.4, label=f"{labels[key]}, SIMPLEC")
        ax.semilogy(hc["iteration"], hc[key], color=c, lw=1.4, ls="--", label=f"{labels[key]}, coupled")
    ax.set_xlabel("outer iteration")
    ax.set_ylabel("relative norm of R(q)")
    ax.set_title("one residual R(q), two algorithms", fontsize=10)
    ax.grid(True, which="major", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
