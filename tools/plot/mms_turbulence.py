"""Manufactured solutions of the k-ε-ζ-f₀ model, the curvature correction
and the QCR (gallery 184): the rotation function f_r1 of the prescribed
flow of the curvature test, and the error of each check against the cell
size. The error values are those printed by

    cargo test --release -p foehncraft-flow --lib turbulence::mms -- --nocapture

    pip install matplotlib numpy
    python3 -m tools.plot.mms_turbulence docs/figures/gallery/mms_turbulence.png
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PI, U0, BETA = np.pi, 0.25, 1.0
H = 1e-4


def vel(x, y):
    sx, cx, sy, cy = np.sin(PI * x), np.cos(PI * x), np.sin(PI * y), np.cos(PI * y)
    s2x, c2x, s2y, c2y = (np.sin(2 * PI * x), np.cos(2 * PI * x),
                          np.sin(2 * PI * y), np.cos(2 * PI * y))
    return (U0 * (sx * cy + BETA * s2x * cy + BETA * sx * c2y),
            -U0 * (cx * sy + 2 * BETA * c2x * sy + 0.5 * BETA * cx * s2y))


def tensors(x, y):
    g = np.zeros((2, 2) + np.shape(x))
    for j, (dx, dy) in enumerate(((H, 0.0), (0.0, H))):
        up = vel(x + dx, y + dy)
        um = vel(x - dx, y - dy)
        for i in range(2):
            g[i, j] = (up[i] - um[i]) / (2 * H)
    return 0.5 * (g + g.swapaxes(0, 1)), 0.5 * (g - g.swapaxes(0, 1))


def fr1(x, y, omega=10.0):
    x = np.asarray(x) + 1e-9
    y = np.asarray(y) + 1e-9
    """f_r1 of Smirnov & Menter (2009), c_r1 = 1, c_r2 = 2, c_r3 = 1."""
    s, w = tensors(x, y)
    sm = np.sqrt(2 * (s * s).sum((0, 1)))
    om = np.sqrt(2 * (w * w).sum((0, 1)))
    rs = sm / om
    h = 1e-3
    ds = np.zeros_like(s)
    u, v = vel(x, y)
    ds = (u * (tensors(x + h, y)[0] - tensors(x - h, y)[0])
          + v * (tensors(x, y + h)[0] - tensors(x, y - h)[0])) / (2 * h)
    num = (np.einsum("ik...,jk...->ij...", w, s) * ds).sum((0, 1))
    d2 = np.maximum(sm ** 2, 0.09 * omega ** 2)
    rt = 2 * num / (om * d2 ** 1.5)
    f = 4 * rs / (1 + rs) * (1 - np.arctan(2 * rt)) - 1
    return np.clip(f, 0.0, 1.25), rt


# (h, L2 error) of the three levels, from the test output
DATA = {
    "ζ-f₀, k": {
        "hex": [(1.3975e-2, 3.721e-4), (6.9877e-3, 9.584e-5), (3.4939e-3, 2.342e-5)],
        "poly": [(1.3153e-2, 4.833e-4), (6.7760e-3, 1.243e-4), (3.4401e-3, 3.014e-5)],
        "tet": [(9.8821e-3, 2.244e-4), (4.9411e-3, 5.095e-5), (2.4705e-3, 1.164e-5)],
    },
    "curvature correction, k": {
        "hex": [(1.3975e-2, 1.260e-2), (6.9877e-3, 3.435e-3), (3.4939e-3, 9.081e-4)],
        "poly": [(1.3153e-2, 1.448e-2), (6.7760e-3, 3.893e-3), (3.4401e-3, 1.010e-3)],
        "tet": [(9.8821e-3, 1.362e-2), (4.9411e-3, 1.343e-2), (2.4705e-3, 1.296e-2)],
    },
    "QCR2000, force": {
        "hex": [(1.3975e-2, 3.355e-2), (6.9877e-3, 7.354e-3), (3.4939e-3, 1.767e-3)],
        "poly": [(1.3153e-2, 6.953e-2), (6.7760e-3, 1.868e-2), (3.4401e-3, 4.546e-3)],
        "tet": [(9.8821e-3, 2.532e-1), (4.9411e-3, 1.170e-1), (2.4705e-3, 5.706e-2)],
        "hex, perturbed 15 %": [(1.3975e-2, 8.494e-2), (6.9877e-3, 5.746e-2), (3.4939e-3, 5.501e-2)],
    },
}
COLOURS = {"hex": "#1f77b4", "poly": "#2ca02c", "tet": "#d62728",
           "hex, perturbed 15 %": "#9467bd"}


def main(out):
    fig = plt.figure(figsize=(10, 8), dpi=100, constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax = fig.add_subplot(gs[0, 0])
    x, y = np.meshgrid(np.linspace(0, 1, 400), np.linspace(0, 1, 400))
    f, _ = fr1(x, y)
    im = ax.pcolormesh(x, y, f, cmap="viridis", shading="auto", vmin=0, vmax=1.25)
    u, v = vel(x, y)
    ax.streamplot(x[0], y[:, 0], u, v, color="white", linewidth=0.8, density=1.4,
                  arrowsize=0.8)
    ax.set_aspect("equal")
    ax.set_title("Flow of the curvature test, f_r1")
    fig.colorbar(im, ax=ax, shrink=0.8, label="f_r1")
    for k, (title, series) in enumerate(DATA.items()):
        a = fig.add_subplot(gs[(k + 1) // 2, (k + 1) % 2])
        for name, pts in series.items():
            h, e = zip(*pts)
            p = np.log(e[1] / e[2]) / np.log(h[1] / h[2])
            a.loglog(h, e, "o-", color=COLOURS[name], label=f"{name} (p {p:.2f})")
        h0 = series["hex"][2][0]
        e0 = series["hex"][2][1]
        a.loglog([h0, 4 * h0], [e0, 16 * e0], "k:", lw=0.8, label="slope 2")
        a.set_title(title)
        a.set_xlabel("cell size h")
        a.set_ylabel("relative L2 error")
        a.grid(True, which="both", alpha=0.3)
        a.legend(fontsize=7)
    fig.savefig(out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "docs/figures/gallery/mms_turbulence.png")
