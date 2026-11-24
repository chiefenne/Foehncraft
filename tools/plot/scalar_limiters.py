"""Bounded convection of the scalars and the energy equation (D237, gallery
entry 211).

scalar_limiters.png: the oblique step of `scalar_oblique_step_mlp.toml`
on perturbed polyhedra with the face-neighbour limiters (Venkatakrishnan,
Barth–Jespersen) and the MLP vertex limiter with Nishikawa's R3 (top), the
profile across the step at x = 0.8 (bottom left) and the largest excursion
beyond the bounds in the three cases of the entry (bottom right; the
T-junction and the heat exchanger numbers are those of the gallery text,
from the runs listed there).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.scalar_limiters OUT_DIR VENK BJ MLP MLP0

with the output directories of

    C=examples/cases/flow/scalar_oblique_step_mlp.toml
    target/release/foehn flow $C --out output/limiters/mlp
    target/release/foehn flow $C --out output/limiters/mlp0 --set solver.scalar.limiter_epsilon=0.0
    target/release/foehn flow $C --out output/limiters/venk --set 'solver.scalar.limiter="venkatakrishnan"'
    target/release/foehn flow $C --out output/limiters/bj --set 'solver.scalar.limiter="barth-jespersen"'
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from tools.plot.colormaps import RAINBOW
from tools.plot.vtkhdf import internal

X_CUT = 0.8
SLOPE = 0.6

# largest excursion beyond the bounds, % of the range of the field, as run
# (gallery text 211): T-junction coarse mesh (theta in [0, 1]) and the
# exchanger air inlet (573.15 K, range 280 K), restarts of 150 iterations
EXCURSION = {
    "T-junction": {"Venkatakrishnan": 0.66, "Barth–Jespersen": 0.03, "MLP + R3": 0.19,
                   "MLP + R3, ε = 0": None},
    "heat exchanger": {"Venkatakrishnan": 4.50, "Barth–Jespersen": 15.75, "MLP + R3": 0.21,
                       "MLP + R3, ε = 0": 0.09},
}


def polygons(path):
    """The bottom polygon of every cell of a 2D case one cell thick and its
    cell data."""
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    polys = []
    for a, b in zip(off[:-1], off[1:]):
        p = pts[conn[a:b]]
        q = p[np.abs(p[:, 2] - p[:, 2].min()) < 1e-12][:, :2]
        m = q.mean(axis=0)
        polys.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
    return polys, data


def extremes(run):
    with open(os.path.join(run, "flow.json")) as f:
        s = json.load(f)["solves"][0]
    c = s["scalars"][0]
    return c["min"], c["max"], s["iterations"]


def field(ax, polys, c, title):
    pc = PolyCollection(polys, array=c, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    pc.set_clim(0.0, 1.0)
    ax.add_collection(pc)
    # cells more than 0.1 % beyond the bounds, outlined
    over = [p for p, v in zip(polys, c) if v > 1.0 + 1e-3]
    under = [p for p, v in zip(polys, c) if v < -1e-3]
    ax.add_collection(PolyCollection(over, facecolors="none", edgecolors="#ff2d55", linewidths=1.4))
    ax.add_collection(PolyCollection(under, facecolors="none", edgecolors="#00e5ff", linewidths=1.4))
    ax.plot([0.0, 1.0], [0.0, SLOPE], color="w", lw=0.8, ls="--")
    ax.axvline(X_CUT, color="w", lw=0.6, ls=":")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    return pc


def profile(polys, c):
    """Cells crossed by the line x = X_CUT, ordered by y."""
    ys, vs = [], []
    for p, v in zip(polys, c):
        if p[:, 0].min() <= X_CUT <= p[:, 0].max():
            ys.append(p[:, 1].mean())
            vs.append(v)
    o = np.argsort(ys)
    return np.array(ys)[o], np.array(vs)[o]


def main():
    out, venk, bj, mlp, mlp0 = sys.argv[1:6]
    runs = [("Venkatakrishnan", venk, "#d95f02"), ("Barth–Jespersen", bj, "#7570b3"),
            ("MLP + R3", mlp, "#1b9e77"), ("MLP + R3, ε = 0", mlp0, "#e7298a")]
    fields = {name: polygons(os.path.join(run, "fields.vtkhdf")) for name, run, _ in runs}
    fig = plt.figure(figsize=(12.0, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.9])
    for k, (name, run, _) in enumerate(runs[:3]):
        polys, d = fields[name]
        lo, hi, it = extremes(run)
        sub = {"Venkatakrishnan": "faces", "Barth–Jespersen": "faces", "MLP + R3": "vertices, ε = 0.01 Δc"}[name]
        pc = field(fig.add_subplot(gs[0, k]), polys, d["c"],
                   f"{name} ({sub})\nc in [{lo:.4f}, {hi:.4f}], {it} iterations")
    fig.colorbar(pc, ax=fig.axes[:3], label="c", shrink=0.8, pad=0.01)
    ax = fig.add_subplot(gs[1, 0:2])
    y_step = SLOPE * X_CUT
    ax.step([0.0, y_step, 1.0], [0.0, 1.0, 1.0], where="post", color="k", lw=1.0, label="exact")
    for name, _, col in runs:
        y, v = profile(*[fields[name][0], fields[name][1]["c"]])
        ax.plot(y, v, "o-", ms=3, lw=1.0, color=col, label=name)
    ax.axhline(1.0, color="0.6", lw=0.5)
    ax.axhline(0.0, color="0.6", lw=0.5)
    ax.set_xlim(y_step - 0.25, y_step + 0.25)
    ax.set_xlabel(f"y at x = {X_CUT}")
    ax.set_ylabel("c")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    ax = fig.add_subplot(gs[1, 2])
    cases = ["oblique step"] + list(EXCURSION)
    names = [n for n, _, _ in runs]
    obl = {}
    for name, run, _ in runs:
        lo, hi, _ = extremes(run)
        obl[name] = 100.0 * max(hi - 1.0, -lo, 0.0)
    table = {"oblique step": obl, **EXCURSION}
    w = 0.2
    for k, (name, _, col) in enumerate(runs):
        xs = np.arange(len(cases)) + (k - 1.5) * w
        vals = [table[c][name] for c in cases]
        ax.bar([x for x, v in zip(xs, vals) if v is not None],
               [v if v >= 1e-3 else 0.0 for v in vals if v is not None], w, color=col, label=name)
        for x, v in zip(xs, vals):
            if v is None:
                ax.text(x, 0.03, "stalls", rotation=90, ha="center", va="bottom", fontsize=7)
            elif v < 1e-3:
                ax.text(x, 0.025, "0", ha="center", va="bottom", fontsize=7)
    ax.set_yscale("log")
    ax.set_ylim(0.01, 40.0)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels(cases, fontsize=9)
    ax.set_ylabel("largest excursion beyond the bounds (% of range)")
    ax.grid(alpha=0.3, axis="y")
    fig.savefig(os.path.join(out, "scalar_limiters.png"), dpi=100)
    for name in names:
        print(f"{name}: oblique step {table['oblique step'][name]:.4f} %")


if __name__ == "__main__":
    main()
