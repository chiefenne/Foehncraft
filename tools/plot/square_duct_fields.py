"""The square duct with and without the quadratic constitutive relation of
SST (`square_duct_qcr.toml`, Re_b 3 500, gallery 114): the streamwise
velocity in the cross-section with the streamlines of the secondary flow,
linear SST on the left, QCR2000 on the right. The quarter duct of the
case is mirrored about its two symmetry planes to the whole duct.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.square_duct_fields OUT.png QCR_DIR SST_DIR

e.g.

    target/release/foehn flow examples/cases/flow/square_duct_qcr.toml --out output/duct_qcr
    target/release/foehn flow examples/cases/flow/square_duct_qcr.toml --out output/duct_sst \\
        --set turbulence.qcr=\\"none\\"
    python3 -m tools.plot.square_duct_fields docs/figures/gallery/square_duct_qcr.png \\
        output/duct_qcr output/duct_sst
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import h5py
import matplotlib.pyplot as plt
import numpy as np
from tools.plot.colormaps import RAINBOW
from tools.plot.vtkhdf import internal


def section(d):
    """u, v, w on the cross-section of the first cell layer, mirrored to
    the whole duct (y, z from −1 to 1, the walls at ±1)."""
    with h5py.File(os.path.join(d, "fields.vtkhdf"), "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        cells = g["Connectivity"][:].reshape(-1, 8)
        u, v, w = (g["CellData"][f"velocity_{c}"][:] for c in "xyz")
    c = pts[cells].mean(axis=1)
    sel = c[:, 0] < c[:, 0].min() + 1e-9
    y, z = c[sel, 1], c[sel, 2]
    ys, zs = np.unique(np.round(y, 12)), np.unique(np.round(z, 12))
    iy = np.searchsorted(ys, np.round(y, 12))
    iz = np.searchsorted(zs, np.round(z, 12))
    grid = np.zeros((3, len(zs), len(ys)))
    for k, a in enumerate((u[sel], v[sel], w[sel])):
        grid[k, iz, iy] = a
    # the case has the walls at 0 and the symmetry planes at 1: Y = y − 1
    Y = np.concatenate([ys - 1.0, 1.0 - ys[::-1]])
    Z = np.concatenate([zs - 1.0, 1.0 - zs[::-1]])
    U = np.block([[grid[0], grid[0][:, ::-1]], [grid[0][::-1], grid[0][::-1, ::-1]]])
    V = np.block([[grid[1], -grid[1][:, ::-1]], [grid[1][::-1], -grid[1][::-1, ::-1]]])
    W = np.block([[grid[2], grid[2][:, ::-1]], [-grid[2][::-1], -grid[2][::-1, ::-1]]])
    return Y, Z, U, V, W


def even(a, y, z, yy, zz):
    """`a` on the even grid (yy, zz) that streamplot needs: interpolated
    along y, then along z."""
    rows = np.array([np.interp(yy, y, r) for r in a])
    return np.array([np.interp(zz, z, c) for c in rows.T]).T


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, on, off = sys.argv[1:]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), constrained_layout=True)
    for ax, d, title in ((axes[0], off, "k-ω SST"), (axes[1], on, "k-ω SST with QCR2000")):
        Y, Z, U, V, W = section(d)
        pc = ax.pcolormesh(Y, Z, U, cmap=RAINBOW, vmin=0.0, vmax=1.4, shading="nearest")
        s = np.hypot(V, W)
        smax = s.max()
        yy = np.linspace(Y[0], Y[-1], len(Y))
        zz = np.linspace(Z[0], Z[-1], len(Z))
        if smax > 1e-6:
            Vi, Wi = even(V, Y, Z, yy, zz), even(W, Y, Z, yy, zz)
            ax.streamplot(yy, zz, Vi, Wi, color=np.hypot(Vi, Wi) / smax, cmap="Reds",
                          density=1.6, linewidth=0.8, arrowsize=0.7)
        ax.set_title(f"{title}: largest secondary velocity {smax * 100:.2f} % of U_b",
                     fontsize=10)
        ax.set_aspect("equal")
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xlabel("y / h")
    axes[0].set_ylabel("z / h")
    fig.colorbar(pc, ax=axes, label="u / U_b", shrink=0.9, pad=0.01)
    fig.suptitle("Square duct, Re_b = 2hU_b/ν = 3 500: streamwise velocity and secondary flow",
                 fontsize=11)
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == "__main__":
    main()
