"""The wall heat quantities of `[output] derived` on the laminar isothermal
flat plate (case plate_heated_laminar, gallery 346): the temperature of the
thermal boundary layer, the local Nusselt number Nu_x = q_w x/(k ΔT)
against the similarity solution of Pohlhausen (1921) and his correlation,
and the heat transfer coefficient of the wall cell h = q_w/(T_w − T_P)
beside q_w/ΔT.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.wall_heat OUT.png RUN_DIR

with the output directory of

    target/release/foehn flow examples/cases/flow/plate_heated_laminar.toml --out output/plate_heat
"""
import sys

import matplotlib

matplotlib.use("Agg")
import h5py  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.turbulence_fields import read  # noqa: E402
from tools.verification.pohlhausen import solve  # noqa: E402

PR = 0.71
NU = 1e-4            # m²/s
U = 1.0              # m/s
K = 0.14084507042253522
DT = 1.0             # K, T_w − T∞
STRETCH = 6.0        # across the layer in the field panel


def plate(path):
    """x of the face centroids of the boundary block `plate` and its datasets."""
    with h5py.File(path, "r") as f:
        g = f["VTKHDF/plate"]
        pts = g["Points"][:]
        conn = g["Polygons/Connectivity"][:]
        off = g["Polygons/Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    x = np.array([pts[conn[a:b]].mean(axis=0)[0] for a, b in zip(off[:-1], off[1:])])
    order = np.argsort(x)
    return x[order], {k: v[order] for k, v in data.items()}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    path = f"{run}/fields.vtkhdf"
    quads, _, cells = read(path)
    x, w = plate(path)
    _, theta0, _ = solve(PR)

    fig = plt.figure(figsize=(10, 3.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1])

    # the thermal layer, stretched across
    ax = fig.add_subplot(gs[0])
    q = quads.copy()
    q[:, :, 1] *= STRETCH
    keep = (q[:, :, 1].mean(axis=1) < 0.12 * STRETCH) & (q[:, :, 0].mean(axis=1) > -0.1)
    t = cells["temperature"] - 300.0
    pc = PolyCollection(q[keep], array=t[keep], cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    pc.set_clim(0.0, DT)
    ax.add_collection(pc)
    ax.set_xlim(-0.1, 1.0)
    ax.set_ylim(0.0, 0.12 * STRETCH)
    ax.set_aspect("equal")
    ax.plot([0, 1], [0, 0], color="k", lw=3, solid_capstyle="butt")
    ax.set_yticks(np.array([0, 0.04, 0.08, 0.12]) * STRETCH, ["0", "0.04", "0.08", "0.12"])
    ax.set_xlabel("x (m)")
    ax.set_ylabel(f"y (m, {STRETCH:g}× stretched)")
    ax.set_title("T − T∞ (K) over the heated plate", fontsize=9)
    fig.colorbar(pc, ax=ax, orientation="horizontal", shrink=0.8, pad=0.02)

    # the local Nusselt number from q_w
    ax = fig.add_subplot(gs[1])
    re = U * x / NU
    nu = w["wall_heat_flux"] * x / (K * DT)
    xs = np.linspace(1e-3, 1.0, 400)
    ax.plot(xs, theta0 * np.sqrt(U * xs / NU), "k--", lw=1.2,
            label=f"similarity {theta0:.4f} Re_x^½")
    ax.plot(xs, 0.332 * PR ** (1 / 3) * np.sqrt(U * xs / NU), ":", color="0.45", lw=1.2,
            label="0.332 Re_x^½ Pr^⅓")
    ax.plot(x, nu, color="#d62728", lw=1.6, label="FoehnCraft q_w x/(k ΔT)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 32)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Nu_x")
    ax.set_title("local Nusselt number", fontsize=9)
    ax.legend(fontsize=7, loc="lower right", frameon=False)
    inset = ax.inset_axes([0.14, 0.62, 0.42, 0.3])
    m = (x > 0.02) & (x < 0.98)
    inset.plot(x[m], 100 * (nu[m] / (theta0 * np.sqrt(re[m])) - 1), color="#d62728", lw=1)
    inset.axhline(0, color="k", lw=0.6)
    inset.set_ylim(-2, 2)
    inset.set_xlim(0, 1)
    inset.tick_params(labelsize=6)
    inset.set_ylabel("%", fontsize=6)

    # the coefficient of the wall cell against q_w / ΔT
    ax = fig.add_subplot(gs[2])
    ax.plot(x, w["heat_transfer_coefficient"], color="#1f77b4", lw=1.6,
            label="heat_transfer_coefficient\nq_w/(T_w − T_P) = k/d_n")
    ax.plot(x, w["wall_heat_flux"] / DT, color="#d62728", lw=1.6, label="q_w/(T_w − T∞)")
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("W/(m² K)")
    ax.set_title("heat transfer coefficients", fontsize=9)
    ax.set_ylim(2, 2e4)
    ax.legend(fontsize=7, loc="upper right", frameon=False)

    fig.suptitle("Laminar flat plate at Re_L = 10⁴, Pr 0.71, wall 1 K above the free stream",
                 fontsize=10)
    fig.savefig(out, dpi=100)
    print(out)


if __name__ == "__main__":
    main()
