"""The cycle-averaged loop with the skin swing and fully coupled cycles at
the end (D204, gallery entry 202).

cht_pulsating_channel_skin.png: the thin-skin wall under the pulsating
channel: the temperature at the end of the run (top), the cycle mean of a
probe in the skin minus the periodic state over time for the loop and the
fully coupled run (bottom left) and the temperature over the last cycle at
the three skin probes (bottom right).
cht_skin_swing.png: the swing of the surface temperature along the wall,
estimated by the loop from the effusivity and from the fully coupled cycles
(left), and the change of the wall temperature and of the swing per update
(right).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.cht_skin OUT_DIR LOOP REF

with the output directories of

    target/release/foehn flow examples/cases/flow/cht_pulsating_channel_skin.toml --out output/skin
    target/release/foehn flow examples/cases/flow/cht_pulsating_channel_skin.toml --out output/skin_ref \\
        --set energy.cycle.enabled=false --set time.end_time=10000.0 --set output.fields_every=4000
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.cht_cycle import cycle_mean, updates
from tools.plot.energy_fields import line
from tools.plot.turbulence_fields import read

PERIOD = 50.0
STRETCH = 2.5
SKIN = (
    ("skin_upstream_T", "x = 0.5 m", "#1b9e77"),
    ("skin_middle_T", "x = 2.0 m", "#d95f02"),
    ("skin_downstream_T", "x = 3.5 m", "#7570b3"),
)


def faces(run):
    rows = [r for r in open(os.path.join(run, "cycle_faces.csv")) if not r.startswith("#")]
    head = rows[0].strip().split(",")
    data = np.array(
        [[float(v) if v else np.nan for v in r.strip().split(",")[2:]] for r in rows[1:]]
    )
    d = {h: data[:, i] for i, h in enumerate(head[2:])}
    order = np.argsort(d["x"])
    return {h: v[order] for h, v in d.items()}


def last_cycle(p, name):
    t = p["time"]
    sel = t > t[-1] - PERIOD + 1e-9
    return (t[sel] - (t[-1] - PERIOD)) / PERIOD, p[name][sel]


def history(out, loop, ref):
    fig = plt.figure(figsize=(11.0, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0], width_ratios=[1.25, 1.0])
    ax = fig.add_subplot(gs[0, :])
    quads, _, d = read(os.path.join(loop, "fields.vtkhdf"))
    quads = [[(x, STRETCH * y) for x, y in q] for q in quads]
    pc = PolyCollection(
        quads, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.2
    )
    pc.set_clim(290.0, 350.0)
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
    ax.set_title(
        "End of the run at t = 1250 s, after the loop and three fully coupled cycles: the wall "
        "(below the dashed line), its thin skin swinging with the flow", fontsize=10,
    )
    fig.colorbar(pc, ax=ax, label="temperature (K)", shrink=0.9, pad=0.01)

    ax = fig.add_subplot(gs[1, 0])
    pr = line(os.path.join(ref, "probes.csv"))
    pl = line(os.path.join(loop, "probes.csv"))
    tr, mr = cycle_mean(pr, "skin_middle_T")
    tl, ml = cycle_mean(pl, "skin_middle_T")
    final = mr[-1]
    ax.plot(tr, mr - final, "-", color="0.3", lw=1.8, label="fully coupled")
    ax.plot(tl, ml - final, "-", color="#d95f02", lw=2.2, label="loop with the skin swing")
    u = updates(loop)
    for t in u["time"]:
        ax.axvline(t, color="#d95f02", lw=0.5, ls=":")
    ax.axvspan(u["time"][-1], tl[-1], color="#d95f02", alpha=0.12, lw=0)
    ax.text(u["time"][-1] + 60.0, -1.2, "fully coupled\ncycles", fontsize=8, color="#a04000")
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xlim(0.0, 7000.0)
    ax.set_ylim(-8.0, 1.0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("cycle mean at x = 2 m, 2.5 mm below the surface,\nminus the periodic state (K)")
    ax.set_title(
        f"The loop ends at 1250 s ({ml[-1] - final:+.3f} K); fully coupled within 0.1 K at "
        "4400 s", fontsize=10,
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    for name, label, c in SKIN:
        x, y = last_cycle(pr, name)
        ax.plot(x, y, "-", color=c, lw=1.8, label=f"{label}, fully coupled")
        x, y = last_cycle(pl, name)
        ax.plot(x, y, "o", color=c, ms=4.0, mfc="none", label=f"{label}, loop")
    ax.set_xlabel("phase of the last cycle (t/T)")
    ax.set_ylabel("temperature 2.5 mm below the surface (K)")
    ax.set_title("The swing over the last cycle", fontsize=10)
    ax.legend(fontsize=7.5, ncol=1, loc="center right")
    ax.grid(alpha=0.3)
    ax.set_xlim(0.0, 1.8)
    ax.set_xticks(np.arange(0.0, 1.01, 0.2))
    fig.savefig(os.path.join(out, "cht_pulsating_channel_skin.png"), dpi=90)


def swing(out, loop):
    f = faces(loop)
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    ax = axs[0]
    res = 0.5 * (f["resolved_t_wall_max"] - f["resolved_t_wall_min"])
    ax.plot(f["x"], res, "-", color="0.25", lw=2.0, label="fully coupled cycles at the end")
    ax.plot(f["x"], f["t_wall_amplitude"], "--", color="#d95f02", lw=1.8,
            label="loop: from the effusivity √(k ρ c_p)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("swing of the surface temperature, half its range (K)")
    ax.set_title("Swing of the wall surface along the channel", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(bottom=0.0)
    ax = axs[1]
    u = updates(loop)
    ax.semilogy(u["update"], u["wall_change_max"], "o-", color="#d95f02", lw=1.6,
                label="mean wall temperature")
    ax.semilogy(u["update"], u["swing_change"], "s-", color="#1b9e77", lw=1.6,
                label="harmonics of the swing")
    ax.axhline(0.01, color="0.4", lw=0.9, ls="--")
    ax.text(u["update"][-1], 0.012, "tolerance 0.01 K", ha="right", fontsize=8, color="0.3")
    ax.set_xlabel("update")
    ax.set_ylabel("largest change per update (K)")
    ax.set_xticks(u["update"])
    ax.set_title("Convergence of the loop", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")
    fig.savefig(os.path.join(out, "cht_skin_swing.png"), dpi=90)


def main():
    out, loop, ref = sys.argv[1:4]
    history(out, loop, ref)
    swing(out, loop)


if __name__ == "__main__":
    main()
