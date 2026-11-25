"""A dirty STL body on a slice, L13 stage 0 (D225): winding number, signed
distance and the error of the fast winding number (gallery
`docs/gallery/immersed.md`).

    pip install matplotlib numpy
    cargo run --release -p foehncraft-mesh --example surface_slice -- output/surface
    python3 -m tools.plot.surface_slice output/surface docs/figures/gallery/immersed_surface_winding.png
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from tools.plot.colormaps import RAINBOW, UNIFORM


def main(run, png):
    d = np.loadtxt(os.path.join(run, "slice.csv"), delimiter=",", skiprows=1)
    tr = np.loadtxt(os.path.join(run, "trace.csv"), delimiter=",", skiprows=1)
    summary = json.load(open(os.path.join(run, "summary.json")))
    xs, zs = np.unique(d[:, 0]), np.unique(d[:, 1])
    shape = (len(zs), len(xs))
    w = d[:, 2].reshape(shape)
    we = d[:, 3].reshape(shape)
    phi = d[:, 4].reshape(shape)
    ext = (xs[0], xs[-1], zs[0], zs[-1])
    segs = tr[:, :4].reshape(-1, 2, 2)
    flipped = tr[:, 4] > 0

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.7), constrained_layout=True)

    def surface(ax, colour):
        ax.add_collection(LineCollection(segs[~flipped], colors=colour, linewidths=1.1))
        ax.add_collection(LineCollection(segs[flipped], colors="#ff2fd0", linewidths=2.0))
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("z (m)")

    ax = axes[0]
    im = ax.imshow(w, origin="lower", extent=ext, cmap=RAINBOW, vmin=-0.1, vmax=2.1)
    ax.contour(xs, zs, w, levels=[0.5], colors="w", linewidths=1.0, linestyles="--")
    surface(ax, "k")
    fig.colorbar(im, ax=ax, shrink=0.85, label="generalised winding number")
    ax.set_title("winding number: 1 inside, 2 in the overlap,\nhole closed smoothly (dashed: w = ½)")

    ax = axes[1]
    lim = np.abs(phi).max()
    im = ax.imshow(phi, origin="lower", extent=ext, cmap=UNIFORM, vmin=-lim, vmax=lim)
    ax.contour(xs, zs, phi, levels=np.arange(-1.0, 1.01, 0.1), colors="k", linewidths=0.4)
    surface(ax, "w")
    fig.colorbar(im, ax=ax, shrink=0.85, label="signed distance φ (m), < 0 inside")
    ax.set_title("signed distance: exact closest point,\nsign from the winding number")

    ax = axes[2]
    err = np.log10(np.maximum(np.abs(w - we), 1e-12))
    im = ax.imshow(err, origin="lower", extent=ext, cmap="inferno", vmin=-10, vmax=-2)
    surface(ax, "#8fd8ff")
    fig.colorbar(im, ax=ax, shrink=0.85, label="log10 |fast − exact|")
    ax.set_title(
        f"fast (β = 2) vs exact winding number, {summary['triangles']}\n"
        f"triangles: max {summary['fast_vs_exact_max']:.1e}, "
        f"{summary['exact_winding_s'] / summary['fast_winding_s']:.0f}× faster"
    )
    fig.savefig(png, dpi=110)
    print("wrote", png)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
