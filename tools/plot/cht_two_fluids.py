"""Two fluids in one conjugate case (D194, stage MF-A, gallery entry 206).

cht_counterflow_two_fluids.png: the water / oil counterflow exchanger
across a steel plate: the temperature (top, the height stretched), the
velocity along the channels (middle, the oil flows against x), and the
bulk temperatures of the two fluids along the exchanger against the
straight lines of ε–NTU with the fully developed Nu = 5.385 (bottom).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_two_fluids OUT_DIR RUN

with the output directory of

    target/release/foehn flow examples/cases/flow/cht_counterflow_two_fluids.toml --out output/two_fluids
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.turbulence_fields import read
from tools.plot.colormaps import UNIFORM

L, H, T, STRETCH = 0.5, 0.002, 0.0005, 25.0
TW_IN, TO_IN = 340.0, 300.0


def eps_ntu():
    kw, ko, kp = 0.6094956547006251, 0.13, 15.0
    hw, ho = 5.385 * kw / (2 * H), 5.385 * ko / (2 * H)
    u = 1.0 / (1.0 / hw + T / kp + 1.0 / ho)
    c = 996.5570590375 * 0.0025 * H * 4180.413119939875
    ntu = u * L / c
    return ntu / (1.0 + ntu)


def bulk(quads, d, lo, hi):
    """Velocity-weighted mean temperature per column of cells in y ∈ (lo, hi)."""
    xc = quads[:, :, 0].mean(axis=1)
    yc = quads[:, :, 1].mean(axis=1)
    sel = (yc > lo) & (yc < hi)
    xs = np.unique(np.round(xc[sel], 9))
    out = []
    for x in xs:
        m = sel & (np.abs(xc - x) < 1e-9)
        u = np.abs(d["velocity_x"][m])
        out.append((u * d["temperature"][m]).sum() / u.sum())
    return xs, np.array(out)


def field(ax, quads, values, cmap, clim, label):
    q = [[(x, STRETCH * y) for x, y in c] for c in quads]
    pc = PolyCollection(q, array=values, cmap=cmap, edgecolors="face", linewidths=0.2)
    pc.set_clim(*clim)
    ax.add_collection(pc)
    ax.set_xlim(0.0, L)
    ax.set_ylim(0.0, STRETCH * (2 * H + T))
    for y in (H, H + T):
        ax.axhline(STRETCH * y, color="w", lw=0.7, ls="--")
    ax.set_yticks([STRETCH * v for v in (0.0, H, H + T, 2 * H + T)])
    ax.set_yticklabels(["0", "2", "2.5", "4.5"])
    ax.set_ylabel(f"y (mm, ×{STRETCH:g})")
    plt.colorbar(pc, ax=ax, label=label, pad=0.01)


def main():
    out, run = sys.argv[1], sys.argv[2]
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    fig, axs = plt.subplots(3, 1, figsize=(11.0, 9.0), constrained_layout=True,
                            height_ratios=[1.0, 1.0, 1.25])
    field(axs[0], quads, d["temperature"], "inferno", (300.0, 340.0), "temperature (K)")
    axs[0].text(0.01, STRETCH * 0.5 * H, "water →", color="w", va="center")
    axs[0].text(0.01, STRETCH * (H + 0.5 * T), "steel", color="w", va="center", fontsize=8)
    axs[0].text(0.49, STRETCH * (1.5 * H + T), "← oil", color="w", va="center", ha="right")
    field(axs[1], quads, 1000.0 * d["velocity_x"], UNIFORM, (-9.5, 9.5), "velocity x (mm/s)")
    axs[1].set_xlabel("")
    xw, tw = bulk(quads, d, 0.0, H)
    xo, to = bulk(quads, d, H + T, 2 * H + T)
    e = eps_ntu()
    x = np.linspace(0.0, L, 2)
    ax = axs[2]
    ax.plot(xw, tw, "o", ms=2.5, color="#d95f02", label="water, Foehncraft")
    ax.plot(x, TW_IN - 40.0 * e * x / L, "-", color="#d95f02", lw=1.0, alpha=0.7,
            label=f"water, ε–NTU (ε = {e:.4f})")
    ax.plot(xo, to, "o", ms=2.5, color="#1b9e77", label="oil, Foehncraft")
    ax.plot(x, TO_IN + 40.0 * e * (1.0 - x / L), "-", color="#1b9e77", lw=1.0, alpha=0.7,
            label="oil, ε–NTU")
    ax.set_xlim(0.0, L)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("bulk temperature (K)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.savefig(os.path.join(out, "cht_counterflow_two_fluids.png"), dpi=90)
    eps = (TW_IN - tw[-1]) / 40.0
    print(f"ε last column {eps:.4f}, ε–NTU {e:.4f}")


if __name__ == "__main__":
    main()
