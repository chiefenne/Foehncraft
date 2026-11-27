"""ECN Spray G, condition G1, grid orientation test (D337; gallery and recipe
of its pull request): the same case run with the cell lines of the grid
along the plume and at 45° to it.

spray_g_grid.png: left, the fuel vapour mass fraction on the plane of the
plume at SNAPSHOT ms (the sector's plane `zmin`) for the two runs, with
the cell edges drawn, the plume leaving the tip at 33° from the injector
axis (horizontal); right, the liquid penetration (projected liquid volume
above 0.2e-3 mm³/mm²) and the vapour penetration (vapour mass fraction
above 0.1 %) of the runs against the measurements of the Engine Combustion
Network (see tools/plot/spray_g_g1.py for their sources). A third run (the
case spray_g_g1.toml, whose lines run along the injector axis) is drawn
thin when given. The script prints the penetrations at the last time and
the Sauter mean diameter 15 mm from the tip of each run.

    target/release/foehn flow examples/cases/flow/spray_g_g1_grid_aligned.toml --out output/spray_g_grid_aligned
    target/release/foehn flow examples/cases/flow/spray_g_g1_grid_turned.toml --out output/spray_g_grid_turned
    python3 -m tools.plot.spray_g_grid OUT_DIR ALIGNED_RUN TURNED_RUN [AXIS_RUN]
"""
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.spray_g_g1 import (  # noqa: E402
    LIQUID_LOW,
    PDI_D32,
    VAPOUR_PEN,
    cells,
    d32_at_15mm,
    nearest,
    penetrations,
    plume_plane,
    step_of,
    times,
    vapour_field,
)

SNAPSHOT = 0.6  # ms


def side_view(ax, run, title):
    tmap = times(run)
    f = nearest(sorted(glob.glob(os.path.join(run, "fields_*.vtkhdf")), key=step_of), SNAPSHOT, tmap)
    _, xyz, d = cells(f)
    quads, keep = plume_plane(xyz)
    y = vapour_field(d, f)[keep]
    pc = PolyCollection(quads, array=np.maximum(y, 1e-5), cmap=RAINBOW, norm=LogNorm(1e-3, 0.3),
                        edgecolors="0.6", linewidths=0.1)
    ax.add_collection(pc)
    ax.set_xlim(0.0, 45.0)
    ax.set_ylim(0.0, 36.0)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("distance from the tip along the axis (mm)")
    return pc


def main():
    args = sys.argv[1:]
    if len(args) not in (3, 4):
        sys.exit(__doc__)
    out, aligned, turned = args[:3]
    axis = args[3] if len(args) == 4 else None
    os.makedirs(out, exist_ok=True)
    runs = [("grid along the plume (33°)", aligned, "tab:red", "tab:blue", 2.0),
            ("grid at 45° to the plume", turned, "tab:purple", "tab:green", 2.0)]
    if axis:
        runs.append(("grid along the injector axis (spray_g_g1)", axis, "0.35", "0.35", 1.0))
    pens = {name: penetrations(run) for name, run, *_ in runs}

    fig = plt.figure(figsize=(12.0, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    pc = side_view(ax1, aligned, f"vapour at {SNAPSHOT} ms, cell lines along the plume (33°) and across it")
    side_view(ax2, turned, f"vapour at {SNAPSHOT} ms, cell lines at 45° to the plume (-12° and 78°)")
    ax1.set_ylabel("distance from the axis (mm)")
    ax2.set_ylabel("distance from the axis (mm)")
    fig.colorbar(pc, ax=[ax1, ax2], label="IC8H18 vapour mass fraction, plane of the plume", shrink=0.8,
                 pad=0.01, aspect=40, location="bottom")

    a = fig.add_subplot(gs[0, 1])
    m = np.array(LIQUID_LOW)
    a.plot(m[:, 0], m[:, 1], "o", color="k", ms=4, label="measured")
    for name, _, c1, _, lw in runs:
        liq = pens[name][0]
        a.plot(liq[:, 0], liq[:, 1], color=c1, lw=lw, label=name)
    a.axvline(0.78, color="0.6", lw=0.8, ls=":")
    a.set_ylabel("liquid penetration (mm)")
    a.set_xlim(0.0, 1.25)
    a.set_ylim(bottom=0.0)
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="lower right")
    a.set_title("Liquid penetration (projected volume above 0.2e-3 mm³/mm²)", fontsize=10)

    b = fig.add_subplot(gs[1, 1])
    m = np.array(VAPOUR_PEN)
    b.plot(m[:, 0], m[:, 1], "o", color="k", ms=4, label="measured (schlieren)")
    for name, _, _, c2, lw in runs:
        vap = pens[name][1]
        b.plot(vap[:, 0], vap[:, 1], color=c2, lw=lw, label=name)
    b.axvline(0.78, color="0.6", lw=0.8, ls=":")
    b.set_xlabel("time after the start of injection (ms)")
    b.set_ylabel("vapour penetration (mm)")
    b.set_xlim(0.0, 1.25)
    b.set_ylim(bottom=0.0)
    b.grid(alpha=0.3)
    b.legend(fontsize=8, loc="lower right")
    b.set_title("Vapour penetration (mass fraction above 0.1 %)", fontsize=10)

    path = os.path.join(out, "spray_g_grid.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)
    for name, run, *_ in runs:
        liq, vap = pens[name]
        print(f"{name}: liquid {liq[-1, 1]:.2f} mm at {liq[-1, 0]:.2f} ms, vapour {vap[-1, 1]:.2f} mm "
              f"at {vap[-1, 0]:.2f} ms, D32 at 15 mm {d32_at_15mm(run):.1f} µm (measured {PDI_D32:.0f} µm)")
        for t in (0.4, 0.8, 1.2):
            i = np.argmin(np.abs(liq[:, 0] - t))
            j = np.argmin(np.abs(vap[:, 0] - t))
            print(f"    {t:.1f} ms: liquid {liq[i, 1]:.2f} mm, vapour {vap[j, 1]:.2f} mm")


if __name__ == "__main__":
    main()
