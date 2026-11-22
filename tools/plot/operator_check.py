"""The operator exactness check (D235), gallery 188: the error of the
velocity gradient for a quadratic field (relative to h |H|) and of its
gradient, the velocity Hessian of the curvature correction (relative to
|H|), per cell of the cavity of `cavity_sst_operator_check.toml` on
regular hexahedra and on perturbed tetrahedra (triangular prisms). Cells
above `solver.operator_tolerance` (0.01) are the flagged ones. From the
`operators.csv` and `fields.vtkhdf` of two runs:

    F=examples/cases/flow/cavity_sst_operator_check.toml
    target/release/foehn flow $F --out output/ops_tet --set output.vtkhdf=true
    target/release/foehn flow $F --out output/ops_hex --set output.vtkhdf=true \\
        --set mesh.generate.type=\\"hex\\" --set mesh.generate.perturb=0
    pip install matplotlib h5py numpy
    python3 -m tools.plot.operator_check docs/figures/gallery/operator_check.png \\
        output/ops_hex output/ops_tet
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm


def polygons(path):
    """The bottom face of every cell of a one-layer 2D mesh, in cell order."""
    import h5py
    with h5py.File(path, "r") as f:
        g = f["VTKHDF/internal"]
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
    out = []
    for a, b in zip(off[:-1], off[1:]):
        c = pts[conn[a:b]]
        zm = c[:, 2].min()
        q = c[np.abs(c[:, 2] - zm) < 1e-12][:, :2]
        m = q.mean(axis=0)
        out.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
    return out


def errors(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["gradient"], d["hessian"]


def panel(ax, polys, v, title):
    pc = PolyCollection(polys, array=np.maximum(v, 1e-4), cmap="magma",
                        norm=LogNorm(1e-3, 2.0), edgecolors="0.35", linewidths=0.2)
    ax.add_collection(pc)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    n = int((v > 0.01).sum())
    ax.set_title(f"{title}: {n} of {len(v)} cells above 0.01", fontsize=11)
    return pc


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, hexd, tetd = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    ph, pt = polygons(hexd / "fields.vtkhdf"), polygons(tetd / "fields.vtkhdf")
    gh, hh = errors(hexd / "operators.csv")
    gt, ht = errors(tetd / "operators.csv")
    fig, axes = plt.subplots(3, 1, figsize=(9, 11.5), constrained_layout=True)
    panel(axes[0], ph, hh, "regular hexahedra: Hessian error")
    panel(axes[1], pt, gt, "perturbed tetrahedra: gradient error")
    pc = panel(axes[2], pt, ht, "perturbed tetrahedra: Hessian error")
    fig.colorbar(pc, ax=axes, shrink=0.6, pad=0.02,
                 label="error for a quadratic field (gradient / h|H|, Hessian / |H|)")
    fig.savefig(out, dpi=100)


if __name__ == "__main__":
    main()
