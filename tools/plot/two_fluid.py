"""Bubbly pipe flow of the two-fluid model (L3 stage MP3; gallery and recipe
of the MP3 pull request).

pipe: two_fluid_pipe.png, the gas fraction α_g along the pipe for j_G =
0.1 and 0.4 m/s (axis to wall, mirrored to show the whole pipe, the
radius stretched), and the radial profiles at L/D = 35 against the
measurements of Wang et al. (1987).

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/bubbly_pipe_wang_jg01.toml --out output/bubbly_pipe/jg01
    target/release/foehn flow examples/cases/flow/bubbly_pipe_wang_jg04.toml --out output/bubbly_pipe/jg04
    python3 -m tools.plot.two_fluid pipe OUT_DIR output/bubbly_pipe
"""
import csv
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.colormaps import RAINBOW

RADIUS = 0.028575
DIAMETER = 2.0 * RADIUS
DATA = os.path.join(
    os.path.dirname(__file__), "..", "..", "examples", "studies", "bubbly_pipe", "data", "wang1987_void.csv"
)
RUNS = (("jg01", 0.1), ("jg04", 0.4))


def cells(path):
    """Cell centres (x, r) and the cell fields, in the original cell order."""
    with h5py.File(path, "r") as f:
        g = f["VTKHDF/internal"] if "internal" in f["VTKHDF"] else f["VTKHDF"]
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        offs = g["Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    # hexahedra and the prisms at the axis of the wedge
    n = np.diff(offs)
    c = np.add.reduceat(pts[conn], offs[:-1], axis=0) / n[:, None]
    order = np.argsort(data["OriginalCellIds"]) if "OriginalCellIds" in data else np.arange(len(c))
    r = np.hypot(c[:, 1], c[:, 2])
    return c[order, 0], r[order], {k: v[order] for k, v in data.items()}


def need(d, name, path):
    if name not in d:
        have = ", ".join(sorted(d))
        raise SystemExit(f"{path} has no cell field `{name}`; it has {have}")
    return d[name]


def line(path):
    with open(path) as f:
        rows = [r for r in f if not r.startswith("#")]
    rd = list(csv.DictReader(rows))
    if not rd or "alpha_gas" not in rd[0]:
        have = ", ".join(rd[0].keys()) if rd else "nothing"
        raise SystemExit(f"{path} has no column `alpha_gas`; it has {have}")
    return np.array([float(r["y"]) for r in rd]), np.array([float(r["alpha_gas"]) for r in rd])


def measured():
    out = {}
    with open(DATA) as f:
        rows = [r for r in f if not r.startswith("#")]
    for r in csv.DictReader(rows):
        out.setdefault(float(r["jg"]), []).append((float(r["r_over_R"]), float(r["alpha"])))
    return out


def pipe(out_dir, base):
    data = measured()
    fig = plt.figure(figsize=(14.0, 6.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.3, 1.0], hspace=0.35, wspace=0.18)
    colours = {0.1: "#1f6fb4", 0.4: "#d1495b"}
    axp = fig.add_subplot(gs[:, 1])
    for row, (key, jg) in enumerate(RUNS):
        run = os.path.join(base, key)
        fpath = os.path.join(run, "fields.vtkhdf")
        x, r, d = cells(fpath)
        a = need(d, "volume_fraction", fpath)
        xs, rs = np.unique(np.round(x, 9)), np.unique(np.round(r, 9))
        nx, nr = len(xs), len(rs)
        A = a.reshape(nr, nx)
        ax = fig.add_subplot(gs[row, 0])
        # mirrored about the axis: the whole pipe, the radius stretched
        rr = np.concatenate([-rs[::-1], rs]) / RADIUS
        AA = np.concatenate([A[::-1], A])
        vmax = max(0.05, float(np.quantile(a, 0.995)))
        m = ax.pcolormesh(xs / DIAMETER, rr, AA, cmap=RAINBOW, vmin=0.0, vmax=vmax, shading="nearest")
        ax.axvline(35.0, color="white", lw=1.0, ls="--")
        ax.set_ylabel("r/R")
        ax.set_title(f"gas fraction, $j_G$ = {jg} m/s, $j_L$ = 1.08 m/s (radius stretched)", fontsize=10)
        if row == 1:
            ax.set_xlabel("x/D (flow upwards, gas enters uniform at x = 0)")
        fig.colorbar(m, ax=ax, pad=0.01, label="α_g")
        lp = os.path.join(run, "line_ld35.csv")
        y, al = line(lp)
        axp.plot(y / RADIUS, al, color=colours[jg], lw=2.2, label=f"Foehncraft, $j_G$ = {jg} m/s")
        px, py = zip(*data[jg])
        axp.plot(px, py, "o", color=colours[jg], mfc="white", mew=1.6, label=f"Wang et al. 1987, $j_G$ = {jg} m/s")
    axp.set_xlim(0.0, 1.0)
    axp.set_ylim(bottom=0.0)
    axp.set_xlabel("r/R")
    axp.set_ylabel("void fraction α_g at L/D = 35")
    axp.grid(alpha=0.3)
    axp.legend(fontsize=8, frameon=False, loc="upper left")
    out = os.path.join(out_dir, "two_fluid_pipe.png")
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"wrote {out}")


def main():
    if len(sys.argv) < 4 or sys.argv[1] != "pipe":
        raise SystemExit(__doc__)
    os.makedirs(sys.argv[2], exist_ok=True)
    pipe(sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
