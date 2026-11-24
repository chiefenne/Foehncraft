"""Absolute pressure values and an operating pressure per fluid part (D239,
gallery and recipe 220).

pressure_values_two_parts.png: the water / oil counterflow exchanger of
`cht_two_fluids_absolute.toml` with its outlets at 1 bar (water) and
1.5 bar (oil): the absolute pressure the case file speaks (top, the height
stretched), the gauge pressure p - p_op of each part the solver works with
(middle), and the pressure along the middle of each channel against the
outlet values (bottom).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.pressure_values OUT_DIR RUN

with the output directory of

    target/release/foehn flow examples/cases/flow/cht_two_fluids_absolute.toml --out output/two_fluids_absolute
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.turbulence_fields import read
from tools.plot.colormaps import RAINBOW

L, H, T, STRETCH = 0.5, 0.002, 0.0005, 25.0


def panel(ax, quads, values, cmap, clim, label, plate=True):
    q = [[(x, STRETCH * y) for x, y in c] for c in quads]
    pc = PolyCollection(q, array=values, cmap=cmap, edgecolors="face", linewidths=0.2)
    pc.set_clim(*clim)
    ax.add_collection(pc)
    ax.set_xlim(0.0, L)
    ax.set_ylim(0.0, STRETCH * (2 * H + T))
    if plate:
        ax.axhspan(STRETCH * H, STRETCH * (H + T), color="0.55")
        ax.text(0.25, STRETCH * (H + 0.5 * T), "steel plate", color="w", va="center",
                ha="center", fontsize=8)
    ax.set_yticks([STRETCH * v for v in (0.0, H, H + T, 2 * H + T)])
    ax.set_yticklabels(["0", "2", "2.5", "4.5"])
    ax.set_ylabel(f"y (mm, ×{STRETCH:g})")
    plt.colorbar(pc, ax=ax, label=label, pad=0.01)


def main():
    out, run = sys.argv[1], sys.argv[2]
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    frame = json.load(open(os.path.join(run, "flow.json")))["pressure_values"]
    op = {p["regions"][0]: p["operating_pressure"] for p in frame["parts"]}
    yc = quads[:, :, 1].mean(axis=1)
    xc = quads[:, :, 0].mean(axis=1)
    water, oil = yc < H, yc > H + T
    fluid = water | oil
    p = d["pressure"]
    gauge = np.where(water, p - op["water"], p - op["oil"])
    fig, axs = plt.subplots(3, 1, figsize=(11.0, 9.0), constrained_layout=True,
                            height_ratios=[1.0, 1.0, 1.2])
    panel(axs[0], quads[fluid], p[fluid] / 1e5, RAINBOW, (0.95, 1.55),
          "absolute pressure (bar)")
    axs[0].text(0.01, STRETCH * 0.5 * H, "water →  outlet 1.0 bar at x = 0.5 m", color="w",
                va="center")
    axs[0].text(0.49, STRETCH * (1.5 * H + T), "oil ←  outlet 1.5 bar at x = 0", color="k",
                va="center", ha="right")
    axs[0].set_title("the case file: absolute values (pressure_values = \"absolute\")")
    g = gauge[fluid]
    lim = np.abs(g).max()
    panel(axs[1], quads[fluid], g, "coolwarm", (-lim, lim), "p − p_op of the part (Pa)")
    axs[1].set_title("the solver: gauge pressure of each part, p_op from its own outlet "
                     f"(water {op['water']:.0f} Pa, oil {op['oil']:.0f} Pa)")
    twin = axs[2].twinx()
    lines = []
    for ax, name, sel, colour in ((axs[2], "water", water, "tab:blue"),
                                  (twin, "oil", oil, "tab:orange")):
        ys = np.unique(np.round(yc[sel], 9))
        ymid = ys[np.argmin(np.abs(ys - ys.mean()))]
        mid = sel & (np.abs(yc - ymid) < 1e-9)
        order = np.argsort(xc[mid])
        lines += ax.plot(xc[mid][order], p[mid][order] - op[name], color=colour, lw=2,
                         label=f"{name}: p − {op[name] / 1e5:.1f} bar along the middle "
                               "of its channel")
        ax.set_ylim(bottom=0.0)
    twin.set_ylabel("oil: p − p_op (Pa)", color="tab:orange")
    axs[2].set_xlim(0.0, L)
    axs[2].set_xlabel("x (m)")
    axs[2].set_ylabel("water: p − p_op (Pa)", color="tab:blue")
    axs[2].legend(lines, [ln.get_label() for ln in lines], loc="upper center")
    axs[2].set_title("each part ends at its own outlet pressure (gauge 0); "
                     "the same bits as the gauge case")
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "pressure_values_two_parts.png"), dpi=110)


if __name__ == "__main__":
    main()
