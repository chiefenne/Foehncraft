"""Cycle-averaged coupling of a periodic flow with its solids (D204, gallery
entry 190).

cht_pulsating_channel_cycle.png: the pulsating channel over a heated wall:
the temperature at the end of the loop (top) and the cycle mean of the
heated face over time for the loop (mean-flow solid phase, the default),
the loop with the first-cell Robin solid phase and the fully coupled run
(bottom).
cht_cycle_convergence.png: the change of the wall temperature and the
interface heat flow per update of the two solid phases.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_cycle OUT_DIR LOOP ROBIN REF

with the output directories of

    target/release/foehn flow examples/cases/flow/cht_pulsating_channel_cycle.toml --out output/cyc
    target/release/foehn flow examples/cases/flow/cht_pulsating_channel_cycle.toml --out output/cyc_robin \\
        --set energy.cycle.solid_model=robin --set energy.cycle.max_updates=12
    target/release/foehn flow examples/cases/flow/cht_pulsating_channel_cycle.toml --out output/cyc_ref \\
        --set energy.cycle.enabled=false --set time.end_time=4000.0
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.energy_fields import line
from tools.plot.turbulence_fields import read

PERIOD = 50.0
STRETCH = 2.5


def cycle_mean(p, name="heated_face_T"):
    t, x = p["time"], p[name]
    n = int(round(t[-1] / PERIOD))
    tc, m = [], []
    for k in range(1, n + 1):
        sel = (t > (k - 1) * PERIOD + 1e-9) & (t <= k * PERIOD + 1e-9)
        tc.append(k * PERIOD)
        m.append(x[sel].mean())
    return np.array(tc), np.array(m)


def updates(run):
    rows = [r for r in open(os.path.join(run, "cycle_coupling.csv")) if not r.startswith("#")]
    head = rows[0].strip().split(",")
    data = np.array([[float(v) if v else np.nan for v in r.strip().split(",")] for r in rows[1:]])
    return {h: data[:, i] for i, h in enumerate(head)}


def field(fig, ax, run, title, lim):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    quads = [[(x, STRETCH * y) for x, y in q] for q in quads]
    pc = PolyCollection(
        quads, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.2
    )
    pc.set_clim(*lim)
    ax.add_collection(pc)
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.0, STRETCH * 0.4)
    ax.set_aspect("equal")
    ax.axhline(STRETCH * 0.2, color="c", lw=0.9, ls="--")
    ticks = np.arange(0.0, 0.41, 0.1)
    ax.set_yticks(STRETCH * ticks)
    ax.set_yticklabels([f"{v:.1f}" for v in ticks])
    ax.set_xlabel("x (m)")
    ax.set_ylabel(f"y (m, ×{STRETCH:g})")
    ax.set_title(title, fontsize=10)
    fig.colorbar(pc, ax=ax, label="temperature (K)", shrink=0.9, pad=0.01)


def history(out, loop, robin, ref):
    fig = plt.figure(figsize=(11.0, 8.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[0.75, 1.15])
    field(
        fig, fig.add_subplot(gs[0]), loop,
        "End of the loop at t = 500 s: the wall (below the dashed line) at its cycle-mean "
        "state, 1000 W/m² from below, the pulsating flow above it", (290.0, 345.0),
    )
    ax = fig.add_subplot(gs[1])
    tr, mr = cycle_mean(line(os.path.join(ref, "probes.csv")))
    tl, ml = cycle_mean(line(os.path.join(loop, "probes.csv")))
    tb, mb = cycle_mean(line(os.path.join(robin, "probes.csv")))
    final = mr[-1]
    ax.plot(tr, mr - final, "-", color="0.3", lw=1.8, label="fully coupled")
    ax.plot(tb, mb - final, "-", color="#7570b3", lw=1.8,
            label="loop, first-cell Robin solid phase (12 updates)")
    ax.plot(tl, ml - final, "-", color="#d95f02", lw=2.2,
            label="loop, mean-flow solid phase (default, 4 updates)")
    for t in updates(loop)["time"]:
        ax.axvline(t, color="#d95f02", lw=0.6, ls=":")
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xlim(0.0, 2000.0)
    ax.set_ylim(-35.0, 3.0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("cycle mean of the heated face\nminus the periodic state (K)")
    ax.set_title(
        f"The loop is at the periodic state at 500 s ({ml[-1] - final:+.3f} K); the fully coupled "
        "run is within 0.1 K at 1350 s, the Robin solid phase at 1450 s", fontsize=10,
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out, "cht_pulsating_channel_cycle.png"), dpi=90)


def convergence(out, loop, robin, ref):
    tr, qr = cycle_mean(line(os.path.join(ref, "heat.csv")), "interface-wall_heat")
    q_ref = qr[-1]
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    for run, c, label in (
        (robin, "#7570b3", "first-cell Robin"),
        (loop, "#d95f02", "mean-flow (default)"),
    ):
        u = updates(run)
        axs[0].semilogy(u["update"], u["wall_change_max"], "o-", color=c, lw=1.6, label=label)
        axs[1].semilogy(
            u["update"], np.abs(u["interface-wall_heat"] - q_ref) / q_ref, "o-", color=c,
            lw=1.6, label=label,
        )
    axs[0].axhline(0.01, color="0.4", lw=0.9, ls="--")
    axs[0].text(12.2, 0.012, "tolerance 0.01 K", ha="right", fontsize=8, color="0.3")
    axs[0].set_ylabel("largest change of the wall temperature (K)")
    axs[1].set_ylabel("|cycle-mean heat flow − fully coupled| / fully coupled")
    for ax, t in zip(axs, ("Wall temperature per update", "Interface heat flow per update")):
        ax.set_xlabel("update")
        ax.set_xticks(range(1, 13))
        ax.set_title(t, fontsize=10)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=9)
    fig.savefig(os.path.join(out, "cht_cycle_convergence.png"), dpi=90)


def main():
    out, loop, robin, ref = sys.argv[1:5]
    history(out, loop, robin, ref)
    convergence(out, loop, robin, ref)


if __name__ == "__main__":
    main()
