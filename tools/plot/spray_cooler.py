"""Draws the picture of the spray-cooler chat (docs/agents.md, gallery 312).

Top: the air temperature in the duct with the droplet streams of the last
pass, coloured by diameter. Bottom left: the vapour mass fraction. Bottom
right: what the assistant read from each run, the mass-flow mean outlet
temperature against the water flow, with the 150 °C target.

    python3 tools/plot/spray_cooler.py RUNS_DIR OUT.png

RUNS_DIR holds the run folders w3.0, w3.5 and w3.55 of the session
(output/agent/spray_cooler after a replay); the outlet temperatures are
those of examples/cases/tools/mcp/spray_cooler/results.csv.
"""

import csv
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.plot.colormaps import RAINBOW  # noqa: E402


def cells(path):
    """Cell centres (x, y) and the cell fields of a 2D hex run."""
    with h5py.File(path, "r") as f:
        g = f["VTKHDF/internal"]
        pts = g["Points"][:]
        con = g["Connectivity"][:].reshape(-1, 8)
        c = pts[con].mean(axis=1)
        data = {k: g["CellData"][k][:] for k in ("temperature", "Y_H2O")}
    return c[:, 0], c[:, 1], data


def grid(x, y, v):
    xs, ys = np.unique(np.round(x, 9)), np.unique(np.round(y, 9))
    ix = np.searchsorted(xs, np.round(x, 9))
    iy = np.searchsorted(ys, np.round(y, 9))
    a = np.full((len(ys), len(xs)), np.nan)
    a[iy, ix] = v
    return xs, ys, a


def main(runs, out):
    run = Path(runs) / "w3.55"
    x, y, d = cells(run / "fields.vtkhdf")
    with h5py.File(run / "spray" / "streams.vtkhdf", "r") as f:
        sp = f["VTKHDF/Points"][:]
        dia = f["VTKHDF/PointData/diameter"][:]
    pick = np.arange(0, len(sp), 7)

    fig = plt.figure(figsize=(12, 8.2), dpi=150)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.25], hspace=0.55, wspace=0.25)

    ax = fig.add_subplot(gs[0, :])
    xs, ys, t = grid(x, y, d["temperature"] - 273.15)
    m = ax.pcolormesh(xs, ys, t, cmap=RAINBOW, shading="nearest", vmin=140, vmax=300)
    ax.scatter(sp[pick, 0], sp[pick, 1], color="white", s=0.4, linewidths=0, alpha=0.6)
    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, 0.2)
    ax.set_title("Air at 300 °C, 5 m/s, cooled by 3.55 g/s of water from three nozzles (white: droplet streams)")
    ax.set_ylabel("y (m)")
    fig.colorbar(m, ax=ax, label="temperature (°C)", pad=0.01)

    ax = fig.add_subplot(gs[1, :])
    xs, ys, w = grid(x, y, d["Y_H2O"])
    m = ax.pcolormesh(xs, ys, w, cmap="inferno", shading="nearest")
    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, 0.2)
    ax.set_title("Water vapour mass fraction (both maps: height stretched)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.colorbar(m, ax=ax, label="Y_H2O (-)", pad=0.01)

    ax = fig.add_subplot(gs[2, 0])
    last = np.isclose(x, x.max())
    order = np.argsort(y[last])
    ax.plot(d["temperature"][last][order] - 273.15, y[last][order], color="#1f6fb4")
    ax.axvline(150, color="#c0392b", ls="--", lw=1)
    ax.set_xlabel("temperature at the outlet (°C)")
    ax.set_ylabel("y (m)")
    ax.set_title("Across the outlet: cool core, hot layers at the walls")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[2, 1])
    rows = [r for r in csv.DictReader(open(ROOT / "examples/cases/tools/mcp/spray_cooler/results.csv")) if not r["run"].endswith("-b")]
    q = [float(r["water_g_s"]) for r in rows]
    tc = [float(r["outlet_mass_mean_C"]) for r in rows]
    ax.plot(q, tc, "o-", color="#1f6fb4")
    for a, b, off in zip(q, tc, [(8, -12), (8, 6), (8, -12)]):
        ax.annotate(f"{a:.2f} g/s: {b:.1f} °C", (a, b), textcoords="offset points", xytext=off, fontsize=8, ha="left")
    ax.set_xlim(2.95, 3.75)
    ax.axhline(150, color="#c0392b", ls="--", lw=1)
    ax.text(2.97, 150.6, "target 150 °C", color="#c0392b", fontsize=8)
    ax.set_xlabel("water flow (g/s per 0.1 m of duct depth)")
    ax.set_ylabel("outlet, mass-flow mean (°C)")
    ax.set_ylim(145, 174)
    ax.set_title("What the assistant read from each run")
    ax.grid(alpha=0.3)

    fig.savefig(out, bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
