"""Thermal acceleration of solids in transient CHT (D213, gallery entry 180).

cht_heated_wall_accelerated.png: the heated wall of entry 159 with
`inertia_scale = "auto"`: the temperature at the end of the accelerated
run (top) and the heated face over time against the unaccelerated run
(bottom).
cht_heated_wall_periodic.png: the same wall with a periodic inlet
temperature: the temperature at the end (top) and the cycle mean of the
heated face over time for the true wall, the restore with a ramp and the
restore as a jump (bottom).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_acceleration OUT_DIR ACC REF PER PER_JUMP PER_REF

with the output directories of

    target/release/foehn flow examples/cases/flow/cht_heated_wall_accelerated.toml --out output/acc
    target/release/foehn flow examples/cases/flow/cht_heated_wall_transient.toml --out output/ref
    target/release/foehn flow examples/cases/flow/cht_heated_wall_periodic_accelerated.toml --out output/per
    target/release/foehn flow examples/cases/flow/cht_heated_wall_periodic_accelerated.toml --out output/per_jump \\
        --set energy.acceleration.restore_ramp=0
    target/release/foehn flow examples/cases/flow/cht_heated_wall_periodic_accelerated.toml --out output/per_ref \\
        --set energy.acceleration.inertia_scale=1.0 --set time.end_time=4000.0
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.energy_fields import line
from tools.plot.turbulence_fields import read


def acceleration(run):
    with open(os.path.join(run, "result.json")) as f:
        return json.load(f)["solves"][0]["transient"]["acceleration"]


def field(fig, ax, run, title, lim):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    pc = PolyCollection(
        quads, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.2
    )
    pc.set_clim(*lim)
    ax.add_collection(pc)
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.0, 1.2)
    ax.set_aspect("equal")
    ax.axhline(0.2, color="c", lw=0.9, ls="--")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(title, fontsize=10)
    fig.colorbar(pc, ax=ax, label="temperature (K)", shrink=0.9, pad=0.01)


def marks(ax, a, ymin, ymax):
    for k, (key, label, c) in enumerate((
        ("accelerated_from", "f measured", "#1b9e77"),
        ("settled_at", "settled", "#7570b3"),
        ("restored_at", "true ρc_p from", "#d95f02"),
    )):
        t = a.get(key)
        if t is None:
            continue
        ax.axvline(t, color=c, lw=0.9, ls="--")
        y = ymin + (0.04 + 0.12 * k) * (ymax - ymin)
        ax.text(t, y, f" {label} {t:.0f} s", color=c, fontsize=8)


def steady(out, acc, ref):
    fig = plt.figure(figsize=(11.0, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.15])
    a = acceleration(acc)
    w = a["regions"][0]
    field(
        fig, fig.add_subplot(gs[0]), acc,
        "Accelerated wall at t = 600 s: 1000 W/m² into the wall from below, the channel flow "
        "above it", (300.0, 324.0),
    )
    ax = fig.add_subplot(gs[1])
    pr, pa = line(os.path.join(ref, "probes.csv")), line(os.path.join(acc, "probes.csv"))
    ax.plot(pr["time"], pr["heated_face_T"], "-", color="0.3", lw=1.8, label="true wall")
    ax.plot(pa["time"], pa["heated_face_T"], "-", color="#d95f02", lw=1.8,
            label=f"accelerated: f = {w['factor']:.3f} on ρc_p (auto, τ = {w['time_constant']:.0f} s)")
    marks(ax, a, 300.0, 325.0)
    ax.set_xlim(0.0, 3000.0)
    ax.set_ylim(300.0, 325.0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("heated face temperature (K)")
    end = pr["heated_face_T"][-1]
    ax.set_title(
        f"The same end state, five times sooner: {pa['heated_face_T'][-1]:.3f} K at 600 s against "
        f"{end:.3f} K at 3000 s", fontsize=10,
    )
    ax.legend(loc="center right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out, "cht_heated_wall_accelerated.png"), dpi=90)


def cycle_mean(p, period=50.0):
    t, x = p["time"], p["heated_face_T"]
    n = int(round(t[-1] / period))
    tc, m = [], []
    for k in range(1, n + 1):
        sel = (t > (k - 1) * period + 1e-9) & (t <= k * period + 1e-9)
        tc.append(k * period)
        m.append(x[sel].mean())
    return np.array(tc), np.array(m)


def periodic(out, per, jump, ref):
    fig = plt.figure(figsize=(11.0, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.15])
    a = acceleration(per)
    field(
        fig, fig.add_subplot(gs[0]), per,
        "Periodic inlet temperature 300 ± 50 K over 50 s, at t = 1500 s", (270.0, 335.0),
    )
    ax = fig.add_subplot(gs[1])
    tr, mr = cycle_mean(line(os.path.join(ref, "probes.csv")))
    tp, mp = cycle_mean(line(os.path.join(per, "probes.csv")))
    tj, mj = cycle_mean(line(os.path.join(jump, "probes.csv")))
    final = mr[-1]
    ax.plot(tr, mr - final, "-", color="0.3", lw=1.8, label="true wall")
    ax.plot(tj, mj - final, "-", color="#7570b3", lw=1.8, label="accelerated, restored as a jump")
    ax.plot(tp, mp - final, "-", color="#d95f02", lw=1.8,
            label="accelerated, restored over 3 periods (default)")
    lo = np.floor(min((mj - final)[tj > 500].min(), -1.0))
    marks(ax, a, lo, 0.5)
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xlim(0.0, 4000.0)
    ax.set_ylim(lo, 0.5)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("cycle mean of the heated face\nminus the true periodic state (K)")
    ax.set_title(
        f"At 1500 s: ramp {mp[-1] - final:+.4f} K, jump {mj[-1] - final:+.3f} K "
        "(the jump leaves part of the swing in the deep wall)", fontsize=10,
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out, "cht_heated_wall_periodic.png"), dpi=90)


def main():
    out, acc, ref, per, jump, pref = sys.argv[1:7]
    steady(out, acc, ref)
    periodic(out, per, jump, pref)


if __name__ == "__main__":
    main()
