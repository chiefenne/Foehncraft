"""Water and steam tube circuits, stage SF7b (gallery and recipe of the
circuits).

circuits_evaporator.png: a once-through evaporator tube at 7 MPa: the
tube drawn as a band coloured by the void fraction, the qualities and the
void along the tube, the bulk, saturation and wall temperatures with the
heat transfer regimes, and the CHF ratio.
circuits_drum_loop.png: the natural circulation loop of a drum boiler:
the loop drawn with its tubes coloured by the void fraction, the void and
the quality along a riser, and the pressure around the loop.

    pip install matplotlib numpy
    python3 -m tools.plot.circuits OUT_DIR EVAPORATOR_RUN DRUM_RUN

with the output directories of

    echo '{"case": "examples/cases/circuits/evaporator_tube.toml", "out": "output/evaporator"}' > e.json
    target/release/foehn call circuits.run --input e.json
    echo '{"case": "examples/cases/circuits/drum_natural_circulation.toml", "out": "output/drum"}' > d.json
    target/release/foehn call circuits.run --input d.json
"""
import csv
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from tools.plot.colormaps import RAINBOW

COLUMNS = ["tube", "s_m", "p_Pa", "T_K", "x_eq", "x_flow", "void", "T_wall_K", "regime", "chf_ratio"]


def read(run):
    path = os.path.join(run, "profiles.csv")
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        sys.exit(f"{path}: no segments")
    missing = [c for c in COLUMNS if c not in rows[0]]
    if missing:
        sys.exit(f"{path}: no column(s) {', '.join(missing)}; the file has: {', '.join(rows[0])}")
    tubes = {}
    for r in rows:
        t = tubes.setdefault(r["tube"], {k: [] for k in rows[0] if k != "tube"})
        for k, v in r.items():
            if k == "tube":
                continue
            t[k].append(v if k == "regime" else float(v))
    for t in tubes.values():
        for k in t:
            if k != "regime":
                t[k] = np.array(t[k])
    with open(os.path.join(run, "result.json")) as f:
        result = json.load(f)
    return tubes, result


def tube(tubes, name, run):
    if name not in tubes:
        sys.exit(f"{run}: no tube `{name}`; the run has: {', '.join(tubes)}")
    return tubes[name]


def band(ax, x0, y0, x1, y1, values, width, vmin=0.0, vmax=1.0):
    """A tube from (x0, y0) to (x1, y1) as a thick line coloured by values."""
    n = len(values)
    xs = np.linspace(x0, x1, n + 1)
    ys = np.linspace(y0, y1, n + 1)
    seg = np.stack([np.column_stack([xs[:-1], ys[:-1]]), np.column_stack([xs[1:], ys[1:]])], axis=1)
    lc = LineCollection(seg, cmap=RAINBOW, linewidths=width, capstyle="projecting")
    lc.set_array(np.asarray(values))
    lc.set_clim(vmin, vmax)
    ax.add_collection(lc)
    return lc


REGIME_COLOURS = {
    "liquid": "#d8e6f3",
    "subcooled boiling": "#cfe8cf",
    "saturated boiling": "#f6e7c1",
    "post-CHF": "#f4c7c3",
    "superheated vapour": "#e6d3ef",
}


def evaporator(out, run):
    tubes, _ = read(run)
    t = tube(tubes, "evaporator", run)
    s = t["s_m"]
    fig = plt.figure(figsize=(12, 6.2))
    gs = fig.add_gridspec(3, 3, width_ratios=[0.35, 1, 1], hspace=0.45, wspace=0.35)
    # the tube, vertical, coloured by the void
    ax = fig.add_subplot(gs[:, 0])
    lc = band(ax, 0, 0, 0, s[-1] + s[0], t["void"], 26)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-0.3, s[-1] + 0.5)
    ax.set_xticks([])
    ax.set_ylabel("height (m)")
    ax.set_title("void fraction", fontsize=10)
    ax.annotate("water in\n530 K", (0, 0), (0.15, -0.25), fontsize=8, ha="left")
    ax.annotate("wet steam out", (0, s[-1]), (0.15, s[-1] + 0.2), fontsize=8, ha="left")
    fig.colorbar(lc, ax=ax, fraction=0.5, pad=0.02, shrink=0.8)
    # regimes as background
    def shade(ax):
        start = 0
        for i in range(1, len(s) + 1):
            if i == len(s) or t["regime"][i] != t["regime"][start]:
                lo = s[start] - (s[1] - s[0]) / 2
                hi = s[i - 1] + (s[1] - s[0]) / 2
                ax.axvspan(lo, hi, color=REGIME_COLOURS.get(t["regime"][start], "white"), lw=0)
                start = i
    ax = fig.add_subplot(gs[0, 1:])
    shade(ax)
    ax.plot(s, t["x_eq"], "k--", label="equilibrium quality")
    ax.plot(s, t["x_flow"], color="#1f5fa8", lw=2, label="flowing quality (Levy)")
    ax.plot(s, t["void"], color="#c03030", lw=2, label="void fraction (drift flux)")
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_ylabel("quality, void")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_xlim(0, s[-1] + s[0])
    ax = fig.add_subplot(gs[1, 1:])
    shade(ax)
    tsat = np.array([np.nan if r == "superheated vapour" else 0 for r in t["regime"]])
    ax.plot(s, t["T_K"], color="#1f5fa8", lw=2, label="bulk")
    ax.plot(s, t["T_wall_K"], color="#c03030", lw=2, label="wall (Gnielinski, Chen)")
    ax.set_ylabel("temperature (K)")
    ax.set_xlim(0, s[-1] + s[0])
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for r, c in REGIME_COLOURS.items() if r in t["regime"]]
    labels = [r for r in REGIME_COLOURS if r in t["regime"]]
    l1 = ax.legend(fontsize=8, loc="lower right")
    ax.add_artist(l1)
    ax.legend(handles, labels, fontsize=8, loc="upper left", ncol=3)
    del tsat
    ax = fig.add_subplot(gs[2, 1:])
    ax.plot(s, t["chf_ratio"], color="#7a3fa0", lw=2)
    ax.axhline(1.3, color="grey", ls="--", lw=1)
    ax.text(s[-1], 1.35, "margin 1.3", ha="right", fontsize=8, color="grey")
    ax.set_ylabel("CHF ratio (Bowring)")
    ax.set_xlabel("distance from the inlet (m)")
    ax.set_xlim(0, s[-1] + s[0])
    ax.set_yscale("log")
    fig.suptitle("Once-through evaporator tube, 7 MPa, 400 kW/m², G = 1114 kg/(m² s)", fontsize=12)
    fig.savefig(os.path.join(out, "circuits_evaporator.png"), dpi=110, bbox_inches="tight")


def drum(out, run):
    tubes, result = read(run)
    r = tube(tubes, "risers", run)
    dc = tube(tubes, "downcomer", run)
    rel = tube(tubes, "relief", run)
    fig = plt.figure(figsize=(12, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.3], hspace=0.4, wspace=0.25)
    ax = fig.add_subplot(gs[:, 0])
    # downcomer from the drum (x = 0, z = 22) down to the lower header (z = 0),
    # risers from the lower header up to the upper header (z = 20), relief to the drum
    lc = band(ax, 0.0, 22.0, 0.0, 0.0, dc["void"], 14)
    ax.plot([0.0, 4.0], [0.0, 0.0], color="#555", lw=6, solid_capstyle="round")
    for k in range(5):
        xk = 1.6 + 0.55 * k
        band(ax, xk, 0.0, xk, 20.0, r["void"], 5)
    ax.plot([1.4, 4.0], [20.0, 20.0], color="#555", lw=6, solid_capstyle="round")
    band(ax, 4.0, 20.0, 0.6, 22.0, rel["void"], 10)
    ax.add_patch(plt.Rectangle((-1.0, 21.2), 2.0, 1.6, color="#888", zorder=3))
    ax.text(0, 22.0, "drum", ha="center", va="center", color="white", zorder=4, fontsize=9)
    ax.text(-0.5, 11, "downcomer", rotation=90, ha="right", va="center", fontsize=9)
    ax.text(4.3, 10, "30 risers,\n150 kW/m²", ha="left", va="center", fontsize=9)
    ax.set_xlim(-2, 7.5)
    ax.set_ylim(-1.5, 24)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_ylabel("height (m)")
    fig.colorbar(lc, ax=ax, label="void fraction", fraction=0.05, pad=0.02)
    drum_node = next(n for n in result["nodes"] if n["kind"] == "drum")["drum"]
    ax.set_title(
        f"steam {drum_node['steam_flow']['value']:.2f} kg/s, circulation ratio {drum_node['circulation_ratio']:.1f}",
        fontsize=10,
    )
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(r["s_m"], r["x_eq"], "k--", label="equilibrium quality")
    ax.plot(r["s_m"], r["x_flow"], color="#1f5fa8", lw=2, label="flowing quality")
    ax.plot(r["s_m"], r["void"], color="#c03030", lw=2, label="void fraction")
    ax.set_ylabel("quality, void (riser)")
    ax.legend(fontsize=8, loc="center right")
    ax.set_xlim(0, 20)
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(r["s_m"], r["T_wall_K"], color="#c03030", lw=2, label="wall")
    ax.plot(r["s_m"], r["T_K"], color="#1f5fa8", lw=2, label="bulk")
    ax.set_xlabel("distance along a riser (m)")
    ax.set_ylabel("temperature (K)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 20)
    fig.suptitle("Natural circulation of a drum boiler at 7 MPa", fontsize=12)
    fig.savefig(os.path.join(out, "circuits_drum_loop.png"), dpi=110, bbox_inches="tight")


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, evap, drum_run = sys.argv[1:]
    os.makedirs(out, exist_ok=True)
    evaporator(out, evap)
    drum(out, drum_run)


if __name__ == "__main__":
    main()
