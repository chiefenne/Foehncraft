"""The frozen-flow march (D118) and the cycle-averaged loop with a transient
solid phase (D204), gallery entry 208.

cht_heated_wall_frozen.png: the heated wall under a fluid whose viscosity
falls with the temperature: the temperature at the end of the warm-up
(top), the heated face and the fluid near the wall over time for the
frozen-flow march and the fully coupled run, the coupled phases of the
march shaded (bottom left), and the difference of the march against the
fully coupled run (bottom right).
cht_frozen_steps.png: the time step over time and the steps taken, march
against fully coupled.
cht_pulsating_channel_warmup.png: the pulsating channel with a period of
5 s warmed up by the loop with a transient solid phase: the temperature at
the end (top), the wall at the end of every march against the cycle means
of the fully coupled run (bottom left), and their difference (bottom
right).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_frozen OUT_DIR FROZEN REF WARMUP WARMUP_REF

with the output directories of

    target/release/foehn flow examples/cases/flow/cht_heated_wall_frozen.toml --out output/frz
    target/release/foehn flow examples/cases/flow/cht_heated_wall_frozen.toml --out output/frz_ref \\
        --set energy.frozen_flow.enabled=false
    target/release/foehn flow examples/cases/flow/cht_pulsating_channel_warmup.toml --out output/wu
    target/release/foehn flow examples/cases/flow/cht_pulsating_channel_warmup.toml --out output/wu_ref \\
        --set energy.cycle.enabled=false
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.cht_cycle import updates
from tools.plot.energy_fields import line
from tools.plot.turbulence_fields import read

PERIOD = 5.0
ORANGE, GREY, GREEN, PURPLE = "#d95f02", "0.3", "#1b9e77", "#7570b3"


def frozen_summary(run):
    r = json.load(open(os.path.join(run, "result.json")))
    return r["solves"][0]["transient"]["frozen_flow"]


def field(fig, ax, run, stretch, lim, ymax, title, wall):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    quads = [[(x, stretch * y) for x, y in q] for q in quads]
    pc = PolyCollection(
        quads, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.2
    )
    pc.set_clim(*lim)
    ax.add_collection(pc)
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.0, stretch * ymax)
    ax.set_aspect("equal")
    ax.axhline(stretch * wall, color="c", lw=0.9, ls="--")
    ticks = np.arange(0.0, ymax + 1e-9, 0.2 if ymax > 0.5 else 0.1)
    ax.set_yticks(stretch * ticks)
    ax.set_yticklabels([f"{v:.1f}" for v in ticks])
    ax.set_xlabel("x (m)")
    ax.set_ylabel(f"y (m, ×{stretch:g})" if stretch != 1.0 else "y (m)")
    ax.set_title(title, fontsize=10)
    fig.colorbar(pc, ax=ax, label="temperature (K)", shrink=0.9, pad=0.01)


def interp(t, tr, vr):
    return np.interp(t, tr, vr)


def heated_wall(out, frz, ref):
    pf, pr = line(os.path.join(frz, "probes.csv")), line(os.path.join(ref, "probes.csv"))
    s = frozen_summary(frz)
    fig = plt.figure(figsize=(11.0, 8.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.0], width_ratios=[1.2, 1.0])
    ax = fig.add_subplot(gs[0, :])
    field(
        fig, ax, frz, 1.0, (300.0, 322.0), 1.2,
        "The wall and the channel at t = 3000 s, after the frozen-flow march (the wall below "
        "the dashed line, heated from below with 1000 W/m²)", 0.2,
    )
    ax = fig.add_subplot(gs[1, 0])
    for p in s["phases"]:
        if p["kind"] == "coupled":
            ax.axvspan(p["start"], max(p["end"], p["start"] + 8.0), color=ORANGE, alpha=0.18, lw=0)
    for name, label, c in (
        ("heated_face_T", "heated face", ORANGE),
        ("near_wall_T", "fluid 4 cm above the wall", PURPLE),
    ):
        ax.plot(pr["time"], pr[name], "-", color=GREY, lw=1.6)
        ax.plot(pf["time"], pf[name], "o", color=c, ms=2.6, mfc="none", label=f"{label}, march")
    ax.plot([], [], "-", color=GREY, lw=1.6, label="fully coupled (6000 steps)")
    ax.axvspan(0, 0, color=ORANGE, alpha=0.18, label=f"coupled steps of the march ({s['refreshes']} refreshes)")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("temperature at x = 3 m (K)")
    ax.set_title(
        f"The march: {s['frozen_steps']} steps of 10 s with the flow held, "
        f"{s['coupled_steps']} coupled steps of 0.5 s", fontsize=10,
    )
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_xlim(0.0, 3000.0)
    ax = fig.add_subplot(gs[1, 1])
    worst = 0.0
    for name, label, c in (
        ("heated_face_T", "heated face", ORANGE),
        ("interface_T", "interface", GREEN),
        ("near_wall_T", "fluid near the wall", PURPLE),
        ("channel_T", "channel middle", GREY),
    ):
        d = pf[name] - interp(pf["time"], pr["time"], pr[name])
        worst = max(worst, np.abs(d).max())
        ax.plot(pf["time"], d, "-", color=c, lw=1.5, label=label)
    for p in s["phases"][1:]:
        if p["kind"] == "coupled":
            ax.axvline(p["start"], color=ORANGE, lw=0.6, ls=":")
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("march − fully coupled (K)")
    ax.set_title(f"Largest difference {worst:.3f} K over a rise of 21 K", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_xlim(0.0, 3000.0)
    fig.savefig(os.path.join(out, "cht_heated_wall_frozen.png"), dpi=90)


def steps(out, frz, ref):
    hf, hr = line(os.path.join(frz, "history.csv")), line(os.path.join(ref, "history.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.0), constrained_layout=True)
    ax = axs[0]
    ax.semilogy(hr["time"], hr["dt"], "-", color=GREY, lw=1.6, label="fully coupled")
    ax.semilogy(hf["time"], hf["dt"], "-", color=ORANGE, lw=1.4, drawstyle="steps-pre", label="frozen-flow march")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("time step (s)")
    ax.set_title("The flow's step while coupled, the march's with the flow held", fontsize=10)
    ax.legend(fontsize=9, loc="center right")
    ax.grid(alpha=0.3, which="both")
    ax.set_xlim(0.0, 3000.0)
    ax = axs[1]
    ax.plot(hr["time"], hr["step"], "-", color=GREY, lw=1.8, label=f"fully coupled: {int(hr['step'][-1])} steps")
    ax.plot(hf["time"], hf["step"], "-", color=ORANGE, lw=1.8, label=f"frozen-flow march: {int(hf['step'][-1])} steps")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("steps taken")
    ax.set_title("Steps to the end of the warm-up", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_xlim(0.0, 3000.0)
    fig.savefig(os.path.join(out, "cht_frozen_steps.png"), dpi=90)


def march_points(run, names):
    """The probes at the end of every march: the row of the update step,
    placed at time + march_jump."""
    p = line(os.path.join(run, "probes.csv"))
    u = updates(run)
    t, v = [], {n: [] for n in names}
    for tu, j in zip(u["time"], u["march_jump"]):
        k = int(np.argmin(np.abs(p["time"] - tu)))
        t.append(tu + (0.0 if np.isnan(j) else j))
        for n in names:
            v[n].append(p[n][k])
    return np.array(t), {n: np.array(x) for n, x in v.items()}


def cycle_means(p, name, at):
    t, x = p["time"], p[name]
    return np.array([x[(t > a - PERIOD + 1e-9) & (t <= a + 1e-9)].mean() for a in at])


def warmup(out, run, ref):
    pr = line(os.path.join(ref, "probes.csv"))
    names = ("heated_face_T", "wall_depth_T", "wall_upstream_T", "wall_middle_T", "wall_downstream_T")
    tm, vm = march_points(run, names)
    fig = plt.figure(figsize=(11.0, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0], width_ratios=[1.2, 1.0])
    ax = fig.add_subplot(gs[0, :])
    field(
        fig, ax, run, 2.5, (290.0, 345.0), 0.4,
        "End of the warm-up at t = 1500 s: the wall (below the dashed line) after 15 updates, the "
        "pulsating flow (period 5 s) above it", 0.2,
    )
    ax = fig.add_subplot(gs[1, 0])
    ends = np.arange(PERIOD, pr["time"][-1] + 1e-9, PERIOD)
    for name, label, c in (
        ("heated_face_T", "heated face (x = 3 m)", ORANGE),
        ("wall_depth_T", "mid-depth of the wall (x = 2 m)", GREEN),
        ("wall_upstream_T", "below the surface, x = 0.5 m", PURPLE),
    ):
        ax.plot(ends, cycle_means(pr, name, ends), "-", color=GREY, lw=1.5)
        ax.plot(tm, vm[name], "o", color=c, ms=5.0, label=f"{label}, loop")
    ax.plot([], [], "-", color=GREY, lw=1.5, label="fully coupled, cycle means")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("temperature (K)")
    ax.set_title("The wall at the end of every march against the fully coupled warm-up", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_xlim(0.0, 1500.0)
    ax = fig.add_subplot(gs[1, 1])
    worst = 0.0
    for name, label, c in (
        ("heated_face_T", "heated face", ORANGE),
        ("wall_depth_T", "mid-depth", GREEN),
        ("wall_upstream_T", "surface, x = 0.5 m", PURPLE),
        ("wall_middle_T", "surface, x = 2 m", "#e7298a"),
        ("wall_downstream_T", "surface, x = 3.5 m", "#66a61e"),
    ):
        d = vm[name] - cycle_means(pr, name, tm)
        worst = max(worst, np.abs(d).max())
        ax.plot(tm, d, "o-", color=c, ms=3.5, lw=1.2, label=label)
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("loop − fully coupled cycle mean (K)")
    ax.set_title(f"Largest difference {worst:.2f} K over a rise of 35 K", fontsize=10)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    ax.set_xlim(0.0, 1500.0)
    fig.savefig(os.path.join(out, "cht_pulsating_channel_warmup.png"), dpi=90)


def main():
    out, frz, ref, wu, wu_ref = sys.argv[1:6]
    heated_wall(out, frz, ref)
    steps(out, frz, ref)
    warmup(out, wu, wu_ref)


if __name__ == "__main__":
    main()
