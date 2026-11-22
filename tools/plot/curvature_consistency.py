"""The consistent curvature correction (D236 step 3), gallery 197: the
error of k against the manufactured solution of the curvature correction
of SST (`turbulence::mms`), per cell on tetrahedra with 15 % perturbed
vertices, before (DS_ij/Dt from the least-squares gradient of the
least-squares velocity gradient) and after (S_ij from a second-order
velocity gradient, its gradient from a quadratic least-squares fit of
S_ij on the flagged cells), and the
convergence of the L2 errors of k and ω on regular and perturbed
tetrahedra.

From the per-cell dumps and the logs of two test runs, one on a build
before the change (`old`) and one after (`new`):

    FOEHN_MMS_CURV_DUMP=output/curv_new cargo test --release -p foehncraft-flow \\
        --lib curvature_correction_on_tetrahedra -- --nocapture 2> output/curv_new.log
    pip install matplotlib numpy
    python3 -m tools.plot.curvature_consistency docs/figures/gallery/curvature_consistency.png \\
        output/curv_old output/curv_new

Each folder holds `curvature_tet_<perturb>_<level>.csv` and the log
`<folder>.log` with the lines `curvature tet perturb <p>: h <h> L2 [k, ω]`.
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.tri import Triangulation


def cells(folder, perturb, level):
    d = np.genfromtxt(folder / f"curvature_tet_{perturb}_{level}.csv",
                      delimiter=",", names=True)
    return d["x"], d["y"], np.abs(d["ek"])


def errors(folder):
    text = Path(str(folder) + ".log").read_text()
    out = {}
    for p, h, k, w in re.findall(
            r"curvature tet perturb ([\d.]+): h ([\d.e+-]+) L2 \[([\d.e+-]+), ([\d.e+-]+)\]", text):
        out.setdefault(p, []).append((float(h), float(k), float(w)))
    return {p: np.array(v) for p, v in out.items()}


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, old, new = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    fig = plt.figure(figsize=(10, 9.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.8])
    norm = LogNorm(1e-5, 3e-2)
    for k, (folder, title) in enumerate([(old, "before: gradient of the gradient"),
                                         (new, "after: quadratic fit of S_ij")]):
        ax = fig.add_subplot(gs[0, k])
        x, y, e = cells(folder, "0.15", 2)
        pc = ax.tripcolor(Triangulation(x, y), np.maximum(e, 1e-6), shading="gouraud",
                          cmap="inferno", norm=norm)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        fin = errors(folder)["0.15"][-1, 1] * 100
        ax.set_title(f"{title}\nL2 error of k {fin:.3f} %", fontsize=11)
    fig.colorbar(pc, ax=fig.axes[:2], shrink=0.9, pad=0.02,
                 label="|k_h − k| / max k, perturbed tetrahedra, finest level")
    ax = fig.add_subplot(gs[1, :])
    for folder, ls, tag in [(old, "--", "before"), (new, "-", "after")]:
        for p, e in errors(folder).items():
            mesh = "perturbed tetrahedra" if float(p) > 0 else "tetrahedra"
            c = "tab:red" if float(p) > 0 else "tab:purple"
            ax.loglog(e[:, 0], e[:, 1] * 100, ls, marker="o", color=c, label=f"k, {mesh}, {tag}")
            ax.loglog(e[:, 0], e[:, 2] * 100, ls, marker="s", color=c, alpha=0.5,
                      label=f"ω, {mesh}, {tag}")
    h = np.array([0.0024, 0.0102])
    ax.loglog(h, 0.015 * (h / h[0]) ** 2, "-.", color="0.4", label="order 2")
    ax.set_xlabel("cell size h")
    ax.set_ylabel("L2 error (% of the largest value)")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_ylim(0.003, 2.5)
    ax.legend(fontsize=8, ncol=3, loc="lower right")
    fig.savefig(out, dpi=100)


if __name__ == "__main__":
    main()
