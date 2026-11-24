"""Start values per region and per fluid part (D238, gallery and recipe 221).

start_values_two_parts.png: the water / oil counterflow exchanger of
`cht_two_fluids_start.toml` at its start (before the first iteration):
the temperature of every region from its own table (the water and the oil
at their inlet temperatures, the steel plate between them), the pressure of
each fluid part from its own boundaries (the water at its outlet, the oil
on the Laplace solution between its total-pressure inlet and its outlet),
the oil pressure along its channel at the start and converged, and the
residual history against a uniform start.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.start_values OUT_DIR START CONVERGED UNIFORM

with the output directories of the case run with `max_iterations = 0`
(START), as it is (CONVERGED), and its gauge twin with `[initial]
pressure = 0` and no region tables (UNIFORM).
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.pressure_values import H, L, STRETCH, T, panel
from tools.plot.turbulence_fields import read
from tools.plot.colormaps import RAINBOW


def residuals(run):
    rows = [l for l in open(os.path.join(run, "residuals.csv")) if not l.startswith("#")]
    head = rows[0].strip().split(",")
    data = np.array([[float(v) for v in r.split(",")] for r in rows[1:]])
    return data[:, 0], data[:, head.index("continuity")], data[:, head.index("energy")]


def main():
    out, start, conv, uniform = sys.argv[1:5]
    quads, _, d0 = read(os.path.join(start, "fields.vtkhdf"))
    _, _, d1 = read(os.path.join(conv, "fields.vtkhdf"))
    frame = json.load(open(os.path.join(start, "flow.json")))["pressure_values"]
    op = {p["regions"][0]: p["operating_pressure"] for p in frame["parts"]}
    yc = quads[:, :, 1].mean(axis=1)
    xc = quads[:, :, 0].mean(axis=1)
    water, oil = yc < H, yc > H + T
    fluid = water | oil
    fig = plt.figure(figsize=(11.0, 9.5), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.1])
    ax0 = fig.add_subplot(gs[0, :])
    panel(ax0, quads, d0["temperature"], "coolwarm", (300.0, 340.0), "temperature (K)",
          plate=False)
    for y, text in ((0.5 * H, "water 340 K"), (H + 0.5 * T, "steel plate 320 K"),
                    (1.5 * H + T, "oil 300 K")):
        ax0.text(0.25, STRETCH * y, text, color="k", va="center", ha="center", fontsize=8)
    ax0.set_title("start temperature per region: water and oil at their own inlets "
                  "(\"inlet\"), the plate at 320 K")
    ax1 = fig.add_subplot(gs[1, :])
    g0 = np.where(water, d0["pressure"] - op["water"], d0["pressure"] - op["oil"])
    panel(ax1, quads[fluid], g0[fluid], RAINBOW, (0.0, 46.0), "p − p_op of the part (Pa)")
    ax1.set_title("start pressure per fluid part: the water at its outlet (1 bar), the oil "
                  "on a Laplace solve between 1.50046 and 1.5 bar")
    ax2 = fig.add_subplot(gs[2, 0])
    ys = np.unique(np.round(yc[oil], 9))
    ymid = ys[np.argmin(np.abs(ys - ys.mean()))]
    mid = oil & (np.abs(yc - ymid) < 1e-9)
    o = np.argsort(xc[mid])
    ax2.plot(xc[mid][o], (d0["pressure"][mid][o] - op["oil"]), lw=2, label="start (Laplace)")
    ax2.plot(xc[mid][o], (d1["pressure"][mid][o] - op["oil"]), "--", lw=2, label="converged")
    ax2.set_xlim(0.0, L)
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("oil: p − 1.5 bar (Pa)")
    ax2.set_title("oil pressure along the middle of its channel")
    ax2.legend()
    ax3 = fig.add_subplot(gs[2, 1])
    for run, label, colour in ((conv, "start per region and part", "tab:blue"),
                               (uniform, "uniform start", "tab:orange")):
        it, c, e = residuals(run)
        ax3.semilogy(it, e, lw=1.8, color=colour, label=f"{label}: {int(it[-1])} iterations")
        ax3.semilogy(it, c, lw=0.8, ls=":", color=colour)
    ax3.axhline(1e-9, color="0.4", lw=0.8)
    ax3.set_xlabel("iteration")
    ax3.set_ylabel("energy (solid), continuity (dotted)")
    ax3.set_title("scaled residuals to the criterion 1e-9")
    ax3.legend()
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "start_values_two_parts.png"), dpi=110)


if __name__ == "__main__":
    main()
