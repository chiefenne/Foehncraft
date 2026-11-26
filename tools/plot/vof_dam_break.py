"""The collapse of a water column in air (volume of fluid, L3 stage MP1;
gallery and recipe of the MP1 pull request).

vof_dam_break.png: left, the volume fraction of the water at four times of
the run; right, the surge front Z = z/a against T = t sqrt(2 g / a) of one
or more runs (meshes) with the measurements of Martin & Moyce (1952).

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/vof_dam_break.toml --out output/vof_dam_break
    python3 -m tools.plot.vof_dam_break OUT_DIR output/vof_dam_break [FINER_RUN ...]
"""
import glob
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read

# Martin, J. C. & Moyce, W. J. (1952), Part IV. An experimental study of the
# collapse of liquid columns on a rigid horizontal plane, Phil. Trans. R.
# Soc. A 244, 312-324, table 2 (n² = 2, a = 2¼ in): the surge front Z = z/a
# at T = n t sqrt(g/a). Read from the JSTOR text of the table (to verify
# against the scan, docs/verification/DATA_SOURCES.md section 3h).
MEASURED_T = [1.19, 1.56, 1.87, 2.26, 2.57, 2.88, 3.25, 3.59, 3.89, 4.20, 4.54, 4.85]
MEASURED_Z = [1.44, 1.89, 2.33, 2.78, 3.22, 3.67, 4.11, 4.56, 5.00, 5.44, 5.89, 6.33]

A = 0.05715  # the column width (m)
G = 9.81


def field(d, name):
    if name in d:
        return d[name]
    have = ", ".join(sorted(d))
    raise SystemExit(f"the fields file has no cell field `{name}`; it has {have}")


def vof_csv(run):
    path = os.path.join(run, "vof.csv")
    if not os.path.exists(path):
        raise SystemExit(f"{run} has no vof.csv (not a volume-of-fluid run?)")
    csv = np.genfromtxt(path, delimiter=",", names=True, comments=None, skip_header=2)
    for c in ("time", "x_max"):
        if c not in csv.dtype.names:
            raise SystemExit(f"vof.csv has no column `{c}`; it has {', '.join(csv.dtype.names)}")
    return csv


def front(run):
    """T and Z of a run: the right end of the cells with alpha > 0.5."""
    csv = vof_csv(run)
    quads, _, _ = read(os.path.join(run, "fields.vtkhdf"))
    h = np.ptp(quads[0][:, 0])
    t = csv["time"] * np.sqrt(2.0 * G / A)
    return t, (csv["x_max"] + 0.5 * h) / A, h


def plot(out, runs):
    run = runs[0]
    snaps = sorted(glob.glob(os.path.join(run, "fields_*.vtkhdf")))
    if len(snaps) < 4:
        raise SystemExit(f"{run} has {len(snaps)} field snapshots; 4 or more are needed ([output] fields_every)")
    pick = [snaps[i] for i in np.linspace(1, len(snaps) - 1, 4).round().astype(int)]
    fig = plt.figure(figsize=(13.0, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, width_ratios=[1.35, 1.0])
    pc = None
    for i, path in enumerate(pick):
        quads, _, d = read(path)
        with h5py.File(path, "r") as f:
            time = float(f.attrs["time"]) if "time" in f.attrs else float("nan")
        ax = fig.add_subplot(gs[i, 0])
        pc = PolyCollection(
            quads, array=field(d, "volume_fraction"), cmap=RAINBOW, edgecolors="face", linewidths=0.1
        )
        pc.set_clim(0.0, 1.0)
        ax.add_collection(pc)
        ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max())
        ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
        ax.set_aspect("equal")
        ax.set_ylabel("y (m)")
        ax.set_title(f"t = {time:.3f} s, T = {time * np.sqrt(2 * G / A):.2f}", fontsize=10)
        if i == len(pick) - 1:
            ax.set_xlabel("x (m)")
        else:
            ax.set_xticklabels([])
    fig.colorbar(pc, ax=fig.axes, location="left", label="volume fraction of water", shrink=0.6, pad=0.01)
    b = fig.add_subplot(gs[:, 1])
    b.plot(MEASURED_T, MEASURED_Z, "ko", ms=6, mfc="white", label="Martin & Moyce 1952 (measured)")
    colours = ["tab:blue", "tab:red", "tab:green"]
    for k, r in enumerate(runs):
        t, z, h = front(r)
        b.plot(t, z, color=colours[k % 3], lw=2.0, label=f"THINC, h = a/{A / h:.0f}")
    b.set_xlabel("T = t √(2g/a)")
    b.set_ylabel("surge front Z = z/a")
    b.set_xlim(0.0, 5.5)
    b.set_ylim(1.0, 8.0)
    b.grid(alpha=0.3)
    b.legend(fontsize=9, loc="upper left")
    b.set_title("Surge front of a water column 2a high (n² = 2)", fontsize=11)
    fig.suptitle("Collapse of a water column in air (volume of fluid, density ratio 828)", fontsize=12)
    fig.savefig(os.path.join(out, "vof_dam_break.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    out, runs = sys.argv[1], sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    plot(out, runs)


if __name__ == "__main__":
    main()
