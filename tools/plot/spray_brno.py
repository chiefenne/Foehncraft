"""The Brno air-blast water spray (D322; gallery and recipe of the samples
pull request).

spray_brno.png: below the axis, the gas speed on the symmetry plane x = 0
of the quarter model; above it, the paths of the droplet streams laid into
the plane of their own distance from the axis, coloured by the droplet
diameter.
spray_brno_compare.png: the chart of examples/studies/brno_spray/compare.py.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_brno.toml --out output/spray_brno
    python3 -m tools.plot.spray_brno OUT_DIR output/spray_brno
"""
import os
import subprocess
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import internal

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def plane_x0(path, name):
    """The cell faces on x = 0 as (z, y) quads and the cell values of `name`."""
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        cd = g["CellData"]
        if name not in cd:
            raise SystemExit(f"fields.vtkhdf has no cell field `{name}`; it has {', '.join(sorted(cd))}")
        val = cd[name][:]
    cells = conn.reshape(-1, 8)
    xyz = pts[cells]
    on = np.abs(xyz[:, :, 0].min(axis=1)) < 1e-12
    quads, v = [], []
    for c, a in zip(xyz[on], val[on]):
        q = c[np.abs(c[:, 0]) < 1e-12][:, [2, 1]]
        m = q.mean(axis=0)
        quads.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
        v.append(a)
    return np.array(quads), np.array(v)


def streams(path):
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        x = g["Points"][:]
        pd = g["PointData"]
        for k in ("diameter", "stream_id"):
            if k not in pd:
                raise SystemExit(f"{path} has no point field `{k}`; it has {', '.join(sorted(pd))}")
        return x, pd["diameter"][:], pd["stream_id"][:]


STEP = 8


def plot(out, run):
    quads, w = plane_x0(os.path.join(run, "fields.vtkhdf"), "velocity_z")
    q = quads.copy()
    q[:, :, 1] *= -1.0
    q *= 1e3
    x, d, sid = streams(os.path.join(run, "spray", "streams.vtkhdf"))
    fig, ax = plt.subplots(figsize=(12.0, 5.6), constrained_layout=True)
    ax.add_patch(plt.Rectangle((20.0, 0.0), 60.0, 16.0, color="#14141e", zorder=0))
    pc = PolyCollection(q, array=w, cmap=RAINBOW, edgecolors="face", linewidths=0.2, clim=(0.0, 110.0))
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="gas speed along the axis (m/s)", pad=0.01)
    starts = np.flatnonzero(np.r_[True, sid[1:] != sid[:-1]])
    ends = np.r_[starts[1:], len(sid)]
    segs, dd = [], []
    for a, b in list(zip(starts, ends))[::STEP]:
        r = np.hypot(x[a:b, 0], x[a:b, 1]) * 1e3
        p = np.column_stack([x[a:b, 2] * 1e3, r])
        segs.extend(np.stack([p[:-1], p[1:]], axis=1))
        dd.extend(0.5 * (d[a : b - 1] + d[a + 1 : b]) * 1e6)
    lc = LineCollection(segs, array=np.array(dd), cmap="plasma", clim=(5.0, 45.0), linewidths=0.6, alpha=0.8)
    ax.add_collection(lc)
    fig.colorbar(lc, ax=ax, label="droplet diameter (µm)", pad=0.01)
    for z in (40.0, 60.0):
        ax.axvline(z, color="white", lw=1.0, ls="--")
        ax.text(z + 0.5, -15.0, f"PDA plane z = {z:g} mm", color="white", fontsize=9, va="bottom")
    ax.axhline(0.0, color="white", lw=0.8)
    ax.text(21.0, 15.0, "droplet streams", color="white", fontsize=10, va="top")
    ax.text(21.0, -9.0, "gas", color="white", fontsize=10, va="top")
    ax.set_xlim(20.0, 80.0)
    ax.set_ylim(-16.0, 16.0)
    ax.set_aspect("equal")
    ax.set_xlabel("distance from the atomizer z (mm)")
    ax.set_ylabel("distance from the axis (mm)")
    ax.set_title(
        "Brno air-blast water spray from its measured plane z = 20 mm: gas speed and droplet streams "
        f"({len(starts)} streams, every {STEP}th drawn)",
        fontsize=11,
    )
    fig.savefig(os.path.join(out, "spray_brno.png"), dpi=100)
    plt.close(fig)
    subprocess.run(
        [sys.executable, os.path.join(ROOT, "examples", "studies", "brno_spray", "compare.py"), run,
         "--png", os.path.join(out, "spray_brno_compare.png")],
        check=True,
    )


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(out, run)


if __name__ == "__main__":
    main()
