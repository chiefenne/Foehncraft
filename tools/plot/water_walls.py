"""Water and steam circuits behind the walls, stage SF7c (gallery and
recipe of the coupled walls).

circuit_evaporator_pass.png: the convective evaporator pass of
examples/cases/flow/circuit_evaporator_pass.toml: the gas temperature in
the pass with the two riser banks beside it coloured by their void
fraction (left), and along the height the heat flux into the water, the
temperatures of the gas-side wall, the inner wall and the water, and the
void fraction of each bank (right).

    pip install matplotlib h5py numpy
    python3 -m tools.plot.water_walls OUT_DIR PASS_RUN

with the output directory of

    target/release/foehn flow examples/cases/flow/circuit_evaporator_pass.toml --out output/evaporator_pass
"""
import csv
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read

STRETCH = 8.0  # the width of the pass drawn 8 times its scale
BANK = 0.025  # drawn width of a riser bank (m, true scale)


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(l for l in f if not l.startswith("#")))


def walls(run):
    out = {}
    for r in rows(os.path.join(run, "circuit", "walls.csv")):
        out.setdefault(r["boundary"], []).append(
            (float(r["y_m"]), float(r["q_W_m2"]), float(r["T_wall_K"]), float(r["T_water_K"]))
        )
    return {b: np.array(sorted(v)) for b, v in out.items()}


def profiles(run):
    out = {}
    for r in rows(os.path.join(run, "circuit", "profiles.csv")):
        out.setdefault(r["tube"], []).append(
            (
                float(r["s_m"]),
                float(r["void"]),
                float(r["x_eq"]),
                float(r["T_K"]),
                float(r["T_wall_K"]),
            )
        )
    return {t: np.array(sorted(v)) for t, v in out.items()}


def main(out, run):
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    res = json.load(open(os.path.join(run, "result.json")))
    c = res["solves"][0]["circuits"]
    w = walls(run)
    p = profiles(run)
    fig = plt.figure(figsize=(13.0, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.25, 1.0, 1.0, 1.0])
    # left: the pass, width stretched, with the riser banks beside it
    ax = fig.add_subplot(gs[0])
    q = quads.copy()
    q[:, :, 0] = q[:, :, 0] * STRETCH
    pc = PolyCollection(q, array=d["temperature"], cmap="inferno", edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    cb = fig.colorbar(pc, ax=ax, location="bottom", shrink=0.9, pad=0.01)
    cb.set_label("gas temperature (K)")
    gap = 0.012 * STRETCH
    width = BANK * STRETCH
    cmap = matplotlib.colormaps[RAINBOW]
    norm = matplotlib.colors.Normalize(0.0, max(p[t][:, 1].max() for t in ("risers_left", "risers_right")))
    for tube, x0 in (("risers_left", -gap - width), ("risers_right", 0.2 * STRETCH + gap)):
        s = p[tube][:, 0]
        ds = 6.0 / len(s)
        for si, void in zip(s, p[tube][:, 1]):
            ax.add_patch(
                plt.Rectangle((x0, si - ds / 2), width, ds, color=cmap(norm(void)), linewidth=0)
            )
        ax.add_patch(plt.Rectangle((x0, 0.0), width, 6.0, fill=False, edgecolor="0.3", linewidth=0.8))
    sm = plt.cm.ScalarMappable(norm=norm, cmap=RAINBOW)
    cb2 = fig.colorbar(sm, ax=ax, location="bottom", shrink=0.9, pad=0.01)
    cb2.set_label("void fraction of the risers")
    ax.set_xlim(-gap - width - 0.05, 0.2 * STRETCH + gap + width + 0.05)
    ax.set_ylim(0.0, 6.0)
    ax.set_xticks([0.0, 0.1 * STRETCH, 0.2 * STRETCH], ["0", "0.1", "0.2"])
    ax.set_xlabel(f"x (m), width drawn {STRETCH:.0f}× its scale")
    ax.set_ylabel("height (m)")
    ax.set_title("Flue gas rises between\ntwo water walls", fontsize=10)
    colors = {"xmin": "#d62728", "xmax": "#1f77b4"}
    names = {"xmin": "left wall", "xmax": "right wall"}
    tubes = {"xmin": "risers_left", "xmax": "risers_right"}
    # heat flux into the water
    ax = fig.add_subplot(gs[1])
    for b in ("xmin", "xmax"):
        ax.plot(w[b][:, 1] / 1e3, w[b][:, 0], color=colors[b], label=names[b])
    ax.set_xlabel("heat flux into the water, per gas-side area (kW/m²)")
    ax.set_ylim(0.0, 6.0)
    ax.set_title("Wall heat flux", fontsize=10)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)
    # temperatures
    ax = fig.add_subplot(gs[2])
    for b in ("xmin", "xmax"):
        ax.plot(w[b][:, 2], w[b][:, 0], color=colors[b], label=f"{names[b]}, gas side")
        pr = p[tubes[b]]
        ax.plot(pr[:, 4], pr[:, 0], color=colors[b], ls="--", label=f"{names[b]}, inner wall")
    pr = p["risers_left"]
    ax.plot(pr[:, 3], pr[:, 0], color="k", lw=1.2, label="water (both banks)")
    ax.set_xlabel("temperature (K)")
    ax.set_ylim(0.0, 6.0)
    ax.set_title("Gas-side wall, inner wall, water", fontsize=10)
    ax.legend(loc="center right", fontsize=7.5)
    ax.grid(alpha=0.3)
    # void
    ax = fig.add_subplot(gs[3])
    for b in ("xmin", "xmax"):
        pr = p[tubes[b]]
        ax.plot(pr[:, 1], pr[:, 0], color=colors[b], label=f"{tubes[b]}: void")
        ax.plot(pr[:, 2] * 10.0, pr[:, 0], color=colors[b], ls=":", label=f"{tubes[b]}: 10 x_eq")
    ax.axvline(0.0, color="0.5", lw=0.6)
    ax.set_xlabel("void fraction, 10 × equilibrium quality")
    ax.set_ylim(0.0, 6.0)
    ax.set_title("Steam in the risers", fontsize=10)
    ax.legend(loc="upper left", fontsize=7.5)
    ax.grid(alpha=0.3)
    drum = c["drums"][0]
    tl, tr = c["tubes"]
    fig.suptitle(
        f"Convective evaporator pass on a drum loop at 7 MPa (SF7c): steam {drum['steam_flow']:.3f} kg/s, "
        f"circulation ratio {drum['circulation_ratio']:.0f}\nleft bank {tl['heat'] / 1e3:.0f} kW and "
        f"{tl['flow']:.1f} kg/s, right bank {tr['heat'] / 1e3:.0f} kW and {tr['flow']:.1f} kg/s; "
        f"{c['solves']} circuit solves, coupling residual {c['residual']:.1e}",
        fontsize=10.5,
    )
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "circuit_evaporator_pass.png")
    fig.savefig(path, dpi=110)
    print(path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
