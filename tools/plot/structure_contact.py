"""Pictures of the reed valve with its seat and stop (gallery entry of MS5
contact, L53).

structure_reed_valve.png: the reed of reed_valve_stop.toml lifted off its
seat by a 3 ms pressure pulse. Top: four snapshots of the shell coloured by
its lift (Rainbow Desaturated) under the stop, drawn as the plane the
mid-surface may reach (the stop's face less half the plate thickness).
Middle: the lift of the tip with the stop and the seat. Bottom: the
contact forces of the seat and of the stop at every step.

    pip install matplotlib numpy
    foehn call structure.run --input examples/cases/tools/structure_reed_valve_stop.json
    python3 -m tools.plot.structure_contact output/reed_valve_stop docs/figures/gallery/structure_reed_valve.png
"""

import csv
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from tools.plot.colormaps import colormap
from tools.plot.structure_plate import read_vtk

RAINBOW = colormap()
MM = 1e3


def need(names, have, path):
    missing = [n for n in names if n not in have]
    if missing:
        sys.exit(f"{path} lacks {missing}; it has: {', '.join(sorted(have))}")


def stop_level(x, case):
    """Height (m) the mid-surface reaches under the stop at x."""
    p, n, off = case["point"], case["normal"], case["offset"]
    # (x − p)·n = off on the plane through p with normal n (y = 0)
    return p[2] + (off - (x - p[0]) * n[0]) / n[2]


def main(out_dir, fig_path):
    path = os.path.join(out_dir, "structure.vtk")
    pts, polys, vec = read_vtk(path)
    with open(os.path.join(out_dir, "result.json")) as f:
        res = json.load(f)
    tr = res["transient"]
    snaps = tr["snapshots"]
    names = [f"transient_{i + 1}" for i in range(len(snaps))]
    need(names, vec, path)
    with open(os.path.join(out_dir, "transient.csv")) as f:
        rows = list(csv.DictReader(f))
    need(["time", "tip_uz", "contact_seat", "contact_stop"], rows[0], "transient.csv")
    t = np.array([float(r["time"]) for r in rows]) * 1e3
    tip = np.array([float(r["tip_uz"]) for r in rows]) * MM
    seat = np.array([float(r["contact_seat"]) for r in rows])
    stop = np.array([float(r["contact_stop"]) for r in rows])
    # the stop of the case file: through (0, 0, 0.3 mm), rising 1 mm over the
    # reed, its face 0.1 mm from the mid-surface
    s = 1.0 / 30.0
    plane = {"point": [0.0, 0.0, 0.3e-3], "normal": [s / np.hypot(1, s), 0.0, -1 / np.hypot(1, s)], "offset": 0.1e-3}
    xs = np.array([pts[:, 0].min(), pts[:, 0].max()])
    ys = np.array([pts[:, 1].min(), pts[:, 1].max()])
    top = max(vec[n][:, 2].max() for n in names) * MM
    norm = Normalize(0.0, top)
    fig = plt.figure(figsize=(11.0, 8.0))
    for k, (name, ts) in enumerate(zip(names, snaps)):
        ax = fig.add_axes((0.0 + 0.226 * k, 0.47, 0.235, 0.47), projection="3d")
        ax.set_facecolor((1, 1, 1, 0))
        ax.view_init(elev=20, azim=-62)
        flat = [pts[p] * MM for p in polys]
        ax.add_collection3d(Poly3DCollection(flat, facecolor="0.55", edgecolor="none", alpha=0.35))
        u = vec[name]
        moved = (pts + u) * MM
        val = np.array([u[p, 2].mean() for p in polys]) * MM
        ax.add_collection3d(
            Poly3DCollection([moved[p] for p in polys], facecolor=RAINBOW(norm(val)), edgecolor=(0, 0, 0, 0.12), linewidth=0.2)
        )
        X, Y = np.meshgrid(xs, ys)
        Z = stop_level(X, plane)
        ax.plot_surface(X * MM, Y * MM, Z * MM, color="#8A8F98", alpha=0.25, edgecolor="#555", linewidth=0.5)
        ax.set_xlim(xs[0] * MM, xs[1] * MM)
        ax.set_ylim(ys[0] * MM, ys[1] * MM)
        ax.set_zlim(-0.05, 1.35)
        ax.set_box_aspect((1.0, 0.4, 0.62))
        ax.locator_params(axis="z", nbins=3)
        ax.locator_params(axis="y", nbins=3)
        ax.tick_params(labelsize=7, pad=-2)
        ax.set_title(f"t = {ts * 1e3:.2f} ms", fontsize=10, pad=-10)
        ax.set_xlabel("x (mm)", fontsize=8, labelpad=-6)
        ax.set_zlabel("lift (mm)", fontsize=8, labelpad=-6)
    cax = fig.add_axes((0.93, 0.56, 0.012, 0.3))
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=RAINBOW), cax=cax, label="lift (mm)")
    ax = fig.add_axes((0.08, 0.30, 0.84, 0.2))
    ax.axvspan(0.0, 3.5, color="#FFE8A0", label="pressure pulse, 3 kPa")
    ax.axhline(stop_level(xs[1], plane) * MM, color="#555", ls="--", lw=1.0, label="stop at the tip")
    ax.axhline(0.0, color="0.3", lw=1.0, label="seat")
    ax.plot(t, tip, color="#2A3F9E", lw=1.5, label="tip lift")
    tw = np.interp(np.array(snaps) * 1e3, t, tip)
    ax.scatter(np.array(snaps) * 1e3, tw, c=RAINBOW(norm(tw)), s=55, edgecolor="k", zorder=3, label="snapshots")
    ax.set_ylabel("tip lift (mm)")
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(-0.1, 1.45)
    ax.tick_params(labelbottom=False)
    ax.legend(loc="upper right", fontsize=8, ncol=5)
    ax.grid(alpha=0.3)
    ax = fig.add_axes((0.08, 0.07, 0.84, 0.2))
    ax.vlines(t, 0, stop, color="#D7301F", lw=0.9, label="stop")
    ax.vlines(t, 0, seat, color="#1A9850", lw=0.9, label="seat")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_ylim(0, 300)
    ax.set_xlim(t[0], t[-1])
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("contact force (N)")
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    pen = tr["contact"]["penetration"]
    fig.suptitle(
        f"Reed valve: lifted against its stop, slammed back onto its seat (largest penetration {pen:.0e} m)",
        fontsize=11,
    )
    fig.savefig(fig_path, dpi=100)
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
