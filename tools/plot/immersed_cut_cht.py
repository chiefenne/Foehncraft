"""Conjugate heat transfer through cut cells, L13 stage 2d (gallery and
recipe entry `cut-cht`).

immersed_cut_cht.png: a heated aluminium pin cut out of the plain box mesh
of the Schäfer–Turek channel and kept as a solid region (the temperature
of fluid and pin, the cut cells of both sides outlined), the temperature
of the concentric-cylinder case against the exact solution, and the L2
error of that case on 20, 40 and 80 cells across for the conductivity
ratios 9, 90 and 900. Read from the output directories of

    target/release/foehn flow examples/cases/flow/immersed_cut_pin_cht.toml --out PIN
    for k in 9 90 900; do
      target/release/foehn flow examples/cases/flow/immersed_cut_cylinder_cht.toml --out CYL$k \\
        --set "verification.levels=[0,1,2]" --set constants.ks=$k \\
        --set materials.core-solid.conductivity=$k
    done

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_cut_cht OUT_DIR PIN CYL9 CYL90 CYL900
"""
import json
import math
import os
import sys

import numpy as np

from tools.plot.immersed_cut_plate import read_section

#: the concentric-cylinder case (examples/cases/flow/immersed_cut_cylinder_cht.toml)
R, Q, KF, TR = 0.25, 1000.0, 1.0, 400.0
#: the pin of the channel case
PIN = (0.2, 0.2, 0.05)


def exact(r, ks):
    """The concentric-cylinder temperature: parabolic in the solid, the
    logarithm of r in the fluid."""
    return np.where(r < R, TR + Q * (R * R - r * r) / (4.0 * ks),
                    TR - Q * R * R / (2.0 * KF) * np.log(np.maximum(r, 1e-12) / R))


def centroids(polys):
    return np.array([q.mean(axis=0) for q in polys])


def sweep(path):
    """The cells across and the L2 temperature error of each level of a sweep."""
    r = json.load(open(os.path.join(path, "result.json")))
    l2 = [s["scalars"][0]["l2"] for s in r["solves"]]
    n = [20 * 2 ** i for i in range(len(l2))]
    return np.array(n), np.array(l2), r["scalar_orders"][0][1]


def plot(out, pin, cyl):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    from tools.plot.colormaps import RAINBOW

    fig = plt.figure(figsize=(14.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[2.2, 1.0])

    # the heated pin: temperature of fluid and solid, the cut cells outlined
    ax = fig.add_subplot(gs[0, 0])
    polys, cells, d = read_section(os.path.join(pin, "fields.vtkhdf"))
    x0, x1, y0, y1 = 0.0, 1.0, 0.0, 0.41
    keep = [i for i, q in enumerate(polys) if q[:, 0].min() < x1]
    t = d["temperature"][cells[keep]]
    zone = d["zone"][cells[keep]]
    sel = [polys[i] for i in keep]
    cx, cy, rp = PIN
    near = np.array([math.hypot(q[:, 0].mean() - cx, q[:, 1].mean() - cy) < rp + 0.03
                     for q in sel])
    pc = PolyCollection(sel, array=t, cmap=RAINBOW, edgecolors="none")
    pc.set_clim(300.0, t.max())
    ax.add_collection(pc)
    a = np.linspace(0.0, 2.0 * math.pi, 200)
    ax.plot(cx + rp * np.cos(a), cy + rp * np.sin(a), "w--", lw=0.9, label="exact pin")
    # a zoom on the front of the pin: every cell outlined, cut and solid
    zx0, zx1, zy0, zy1 = cx - rp - 0.004, cx - rp + 0.026, cy + 0.012, cy + 0.042
    ax.add_patch(plt.Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fill=False, ec="k", lw=0.8))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ns = len(np.unique(cells[keep][zone == 1]))
    ax.set_title(f"heated aluminium pin in the Schäfer–Turek channel (Re 20, Pr 1):\n"
                 f"one plain box mesh, the pin cut out and kept as {ns} solid cells",
                 fontsize=11)
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(pc, ax=ax, shrink=0.75, pad=0.01, label="temperature (K)")

    ax = fig.add_subplot(gs[0, 1])
    ring = [q for q, k in zip(sel, near) if k]
    zi = PolyCollection(ring, array=t[near], cmap=RAINBOW, edgecolors=(0, 0, 0, 0.9),
                        linewidths=0.8)
    zi.set_clim(300.0, t.max())
    ax.add_collection(zi)
    ax.plot(cx + rp * np.cos(a), cy + rp * np.sin(a), "w:", lw=0.8, alpha=0.8)
    ax.set_xlim(zx0, zx1)
    ax.set_ylim(zy0, zy1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_title("the front of the pin (small box, left): fluid and\npin cells share each cut face (exact circle dotted)",
                 fontsize=11)

    # the concentric cylinders: every cell against the exact solution
    ax = fig.add_subplot(gs[1, 0])
    rr = np.linspace(0.0, 0.7, 400)
    colours = {9: "#30638e", 90: "#d1495b", 900: "#edae49"}
    for k, path in cyl.items():
        polys, cells, d = read_section(os.path.join(path, "fields.vtkhdf"))
        c = centroids(polys)
        r = np.hypot(c[:, 0] - 0.5, c[:, 1] - 0.5)
        tt = d["temperature"][cells]
        ax.plot(rr, exact(rr, k), "-", color="k", lw=0.9)
        ax.plot(r, tt, ".", color=colours[k], ms=2.5, label=f"k_s / k_f = {k}")
    ax.axvline(R, color="0.4", lw=0.8, ls=":")
    ax.text(R + 0.01, 371, "interface", fontsize=8, color="0.3")
    ax.set_xlim(0.0, 0.7)
    ax.set_xlabel("r (m)")
    ax.set_ylabel("temperature (K)")
    ax.set_title("concentric cylinders, 80 cells across:\n"
                 "every cell against the exact solution (lines)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

    # the order of the error
    ax = fig.add_subplot(gs[1, 1])
    for k, path in cyl.items():
        n, l2, orders = sweep(path)
        ax.loglog(n, l2, "o-", color=colours[k],
                  label=f"k ratio {k}: order {orders[0]:.2f}, {orders[1]:.2f}")
        print(f"k ratio {k}: L2 {', '.join(f'{e:.4e}' for e in l2)}, orders {orders}")
    n = np.array([20.0, 80.0])
    ax.loglog(n, 0.12 * (n / 20.0) ** -2, "k--", lw=0.9, label="second order")
    ax.set_xlabel("cells across the box")
    ax.set_ylabel("L2 error of the temperature (K)")
    ax.set_title("concentric cylinders cut both ways:\nthe error falls with h²", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    path = os.path.join(out, "immersed_cut_cht.png")
    fig.savefig(path, dpi=100)
    print(path)


def main(argv):
    if len(argv) < 6:
        sys.exit(__doc__)
    import matplotlib

    matplotlib.use("Agg")
    plot(argv[1], argv[2], {9: argv[3], 90: argv[4], 900: argv[5]})


if __name__ == "__main__":
    main(sys.argv)
