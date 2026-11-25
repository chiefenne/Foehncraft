"""The deforming mesh of stage MS0 (gallery and recipe 296).

motion_piston.png: the evenly compressed cylinder (L22). Left: p₀ against
the volume over one crank period for the adiabatic walls, on the isentrope
p₀ V^γ = const, and for the walls held at 300 K, whose loop falls below it.
Middle: the temperature of the cooled case at the top dead centre on the
compressed mesh, with velocity arrows. Right: the error of p₀ at the top
dead centre and after the whole period against the time step, falling as
Δt² (BDF2).
motion_free_stream.png: the uniform flow on the swinging mesh of
motion_free_stream.toml at the largest swing, coloured by |u − U|/U, and
that deviation over time (round-off).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.mesh_motion OUT_DIR ADIABATIC_RUN COOLED_RUN SWEEP_DIR FREE_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/motion_piston_adiabatic.toml --out output/adiabatic
    target/release/foehn flow examples/cases/flow/motion_piston_cooled_walls.toml --out output/cooled
    for dt in 0.02 0.01 0.005 0.0025; do for end in 0.5 1.0; do
      target/release/foehn flow examples/cases/flow/motion_piston_adiabatic.toml \\
        --out output/sweep/dt${dt}_end${end} --set time.time_step=$dt --set time.end_time=$end
    done; done
    target/release/foehn flow examples/cases/flow/motion_free_stream.toml --out output/free \\
        --set 'output.fields_times=[0.125]' --set output.history_every=1
"""
import csv
import glob
import json
import os
import sys

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.turbulence_fields import read
from tools.plot.colormaps import RAINBOW

R_GAS = 8.314462618 / 0.0289647
CP = 1005.0
GAMMA = CP / (CP - R_GAS)
P0 = 101_325.0
V0 = 0.02 * 0.04 * 0.002


def history(run):
    with open(os.path.join(run, "history.csv")) as f:
        rows = [r for r in csv.reader(line for line in f if not line.startswith("#"))]
    head, body = rows[0], rows[1:]
    return {k: np.array([float(r[i]) for r in body]) for i, k in enumerate(head)}


def result(run):
    with open(os.path.join(run, "result.json")) as f:
        return json.load(f)["solves"][0]


def piston(out, adiabatic, cooled, sweep):
    fig = plt.figure(figsize=(15.0, 4.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.0])
    # p–V
    ax = fig.add_subplot(gs[0])
    v = np.linspace(V0 / 4, V0, 200)
    ax.plot(v / V0, P0 * (V0 / v) ** GAMMA / 1e3, color="0.6", lw=4, label="isentrope p₀V^γ = const")
    for run, lab, col in [(adiabatic, "adiabatic walls", "#1f5fbf"), (cooled, "walls at 300 K", "#d9480f")]:
        h = history(run)
        vv = np.concatenate([[V0], h["volume"]])
        pp = np.concatenate([[P0], h["p0"]])
        ax.plot(vv / V0, pp / 1e3, "-", color=col, lw=1.4, marker="o", ms=2.2, label=lab)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([0.25, 0.5, 1.0], ["1/4", "1/2", "1"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks([100, 200, 400, 700], ["100", "200", "400", "700"])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("V / V₀")
    ax.set_ylabel("p₀ (kPa)")
    ax.set_title("p₀ over one crank period")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    # the temperature at the top dead centre
    ax = fig.add_subplot(gs[1])
    quads, _, d = read(os.path.join(cooled, "fields_000050.vtkhdf"))
    q = quads * 100.0
    t = d["temperature"]
    pc = PolyCollection(q, array=t, cmap=RAINBOW, edgecolors="none", antialiased=False)
    ax.add_collection(pc)
    c = quads.mean(axis=1) * 100.0
    ux, uy = d["velocity_x"], d["velocity_y"]
    # every fourth column and every eighth row of the 40 × 80 cells
    k = np.arange(len(t))
    sel = (k % 40 % 4 == 2) & (k // 40 % 8 == 4)
    ax.quiver(c[sel, 0], c[sel, 1], ux[sel], uy[sel], color="w", scale=0.4, width=0.005)
    ax.set_xlim(0, 2.0)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_title("walls at 300 K, top dead centre (t = 0.5 s)")
    fig.colorbar(pc, ax=ax, label="T (K)", shrink=0.6)
    # the error against Δt
    ax = fig.add_subplot(gs[2])
    for end, lab, col in [("0.5", "top dead centre", "#1f5fbf"), ("1.0", "after the period", "#2b8a3e")]:
        pts = []
        for run in glob.glob(os.path.join(sweep, f"dt*_end{end}")):
            dt = float(os.path.basename(run)[2:].split("_")[0])
            s = result(run)
            vol = s["transient"]["mesh_motion"]["volume"]
            pe = P0 * (V0 / vol) ** GAMMA
            pts.append((dt, abs(s["density"]["p0"] / pe - 1.0)))
        pts.sort()
        dts, err = np.array(pts).T
        ax.loglog(dts, err, "o-", color=col, label=f"|p₀/p₀,exact − 1|, {lab}")
    ref = np.array([0.0025, 0.02])
    ax.loglog(ref, 2.7e-3 * (ref / 0.02) ** 2, "k--", lw=1, label="Δt²")
    ax.set_xlabel("Δt (s)")
    ax.set_title("adiabatic walls: second order in time")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.savefig(os.path.join(out, "motion_piston.png"), dpi=120)


def free_stream(out, run):
    fig = plt.figure(figsize=(13.0, 4.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0])
    ax = fig.add_subplot(gs[0])
    files = sorted(glob.glob(os.path.join(run, "fields_0*.vtkhdf")))
    quads, _, d = read(files[0])
    dev = np.hypot(d["velocity_x"] - 1.0, d["velocity_y"]) + 1e-18
    pc = PolyCollection(quads, array=np.log10(dev), cmap=RAINBOW, edgecolors="k", linewidths=0.4)
    ax.add_collection(pc)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 0.52)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("uniform flow 1 m/s on the swinging mesh, largest swing (t = 0.125 s)")
    fig.colorbar(pc, ax=ax, label="log₁₀ |u − U| / U")
    ax = fig.add_subplot(gs[1])
    s = result(run)
    h = history(run)
    ax.semilogy(h["time"], np.maximum(h["momentum_x"], 1e-18), label="momentum residual (x)")
    ax.semilogy(h["time"], np.maximum(h["continuity"], 1e-18), label="continuity residual")
    ax.axhline(abs(s["transient"]["mesh_motion"]["closure_max"]), color="k", ls="--", lw=1,
               label="space conservation, worst cell")
    ax.set_ylim(1e-18, 1e-12)
    ax.set_xlabel("t (s)")
    ax.set_title("residuals stay at round-off")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.savefig(os.path.join(out, "motion_free_stream.png"), dpi=110)


def main():
    out, adiabatic, cooled, sweep, free = sys.argv[1:6]
    piston(out, adiabatic, cooled, sweep)
    free_stream(out, free)


if __name__ == "__main__":
    main()
