"""The agent-driven heat sink study (gallery 181) as a figure: the
temperature fields of the fin pitches the assistant ran, on one colour scale
(top), and the mean temperature of the base top and the pressure drop of the
channel over the fin pitch (bottom), read the way the assistant read them
(`line_base_top.csv`, `boundaries.csv`).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.heat_sink_study OUT.png STUDY_DIR

STUDY_DIR holds the runs `pitch_<mm>` of the recorded study
(docs/agents.md, "A recorded study"): output/agent/study.
"""
import csv
import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.cht_turbulent import outline, panel
from tools.plot.energy_fields import grid
from tools.plot.turbulence_fields import read


def rows(path):
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    return list(csv.DictReader(lines))


def numbers(run):
    t = [float(r["temperature"]) for r in rows(os.path.join(run, "line_base_top.csv"))]
    b = {r["boundary"]: float(r["mean_pressure"]) for r in rows(os.path.join(run, "boundaries.csv"))}
    return float(np.mean(t)), float(np.max(t)), b["xmin"] - b["xmax"]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out, study = sys.argv[1:3]
    runs = sorted(glob.glob(os.path.join(study, "pitch_*")),
                  key=lambda r: int(re.search(r"pitch_(\d+)", r).group(1)))
    runs = [r for r in runs if os.path.exists(os.path.join(r, "fields.vtkhdf"))
            and re.fullmatch(r"pitch_\d+", os.path.basename(r))]
    pitch = [int(re.search(r"(\d+)", os.path.basename(r)).group(1)) for r in runs]
    nums = [numbers(r) for r in runs]
    fig = plt.figure(figsize=(12.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05])
    xlim, ylim = (0.07, 0.19), (0.0, 0.04)
    lo = min(float(np.min(read(os.path.join(r, "fields.vtkhdf"))[2]["temperature"])) for r in runs)
    hi = max(float(np.max(read(os.path.join(r, "fields.vtkhdf"))[2]["temperature"])) for r in runs)
    for k, (r, p, n) in enumerate(zip(runs, pitch, nums)):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        quads, _, d = read(os.path.join(r, "fields.vtkhdf"))
        xs, ys, solid = grid(quads, (d["zone"] > 0).astype(float))
        panel(ax, quads, d["temperature"], "inferno", "temperature (K)", fig, xlim, ylim, lo, hi)
        outline(ax, xs, ys, solid, "c")
        ax.set_aspect("equal")
        ax.set_title(f"pitch {p} mm: base top {n[0]:.1f} K, pressure drop {n[2]:.1f} Pa", fontsize=10)
        ax.set_ylabel("y (m)")
        if k // 2 == 1:
            ax.set_xlabel("x (m)")
    ax = fig.add_subplot(gs[2, :])
    ax.plot(pitch, [n[0] for n in nums], "o-", color="#d1495b", lw=1.8, label="base top, mean (K)")
    ax.plot(pitch, [n[1] for n in nums], "s--", color="#d1495b", lw=1.0, alpha=0.7, label="base top, maximum (K)")
    ax.set_xlabel("fin pitch (mm); three fins 3 mm thick and 20 mm tall, 5 m/s, k-ω SST, 2000 W/m² chip")
    ax.set_ylabel("temperature (K)", color="#d1495b")
    ax2 = ax.twinx()
    ax2.plot(pitch, [n[2] for n in nums], "^-", color="#1f6fb4", lw=1.8, label="pressure drop (Pa)")
    ax2.set_ylabel("pressure drop inlet to outlet (Pa)", color="#1f6fb4")
    ax.grid(alpha=0.3)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=9)
    ax.set_title("What the assistant read: one 2D mesh of 15 000 cells per pitch, first-cell y⁺ 5 to 18 "
                 "(LIM002), a sweep for illustration, not a design", fontsize=10)
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")
    for p, n in zip(pitch, nums):
        print(f"pitch {p} mm: mean {n[0]:.2f} K, max {n[1]:.2f} K, dp {n[2]:.2f} Pa")


if __name__ == "__main__":
    main()
