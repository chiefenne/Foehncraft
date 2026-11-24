"""Laminar natural convection benchmarks of stage 7b2 (gallery and recipe
of the natural convection PR).

natural_convection_cavity.png: the differentially heated square cavity at
Ra = 1e3, 1e4, 1e5, 1e6: temperature with streamlines, and the mean Nusselt
number of the hot wall extrapolated from three meshes against Ra.
natural_convection_plate.png: the isothermal vertical plate (Ostrach): the
temperature and the velocity near the plate, the local Nusselt number
along the plate and the profiles across the layer in similarity form,
against the similarity solution of tools/verification/ostrach.py.
mixed_convection_channel.png: fully developed mixed convection between a
hot and a cold wall (Aung & Worku), net flow up and down, against the
exact profiles.

    pip install matplotlib h5py numpy scipy
    python3 -m tools.plot.natural_convection OUT_DIR CAV1E3 CAV1E4 CAV1E5 CAV1E6 PLATE UP DOWN

with the output directories of

    C=examples/cases/flow/natural_convection_cavity_ra1e5.toml
    target/release/foehn flow $C --out output/cav1e3 --set 'gravity.vector=[0,-0.71,0]'
    target/release/foehn flow $C --out output/cav1e4 --set 'gravity.vector=[0,-7.1,0]'
    target/release/foehn flow $C --out output/cav1e5
    target/release/foehn flow $C --out output/cav1e6 --set 'gravity.vector=[0,-710,0]'
    target/release/foehn flow examples/cases/flow/natural_convection_plate.toml --out output/plate
    M=examples/cases/flow/mixed_convection_channel.toml
    target/release/foehn flow $M --out output/up
    target/release/foehn flow $M --out output/down \\
        --set 'boundary.ymin.mass_flow=-2.5e-4' --set 'analytic.mean_velocity=-5e-3'
"""
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

# mean Nusselt numbers of other codes, context only (D201): de Vahl Davis
# (1983) for Ra 1e3 to 1e6, Le Quéré (1991) at 1e6
DE_VAHL_DAVIS = {1e3: 1.118, 1e4: 2.243, 1e5: 4.519, 1e6: 8.800}
LE_QUERE_1E6 = 8.825
DEPTH = 0.01          # cavity depth (m); k = 1, ΔT = 1: Nu = Q/DEPTH


def result(run):
    with open(os.path.join(run, "result.json")) as f:
        return json.load(f)


def estimate(run, quantity):
    """(value of the finest mesh, extrapolated value, GCI fraction)."""
    for e in result(run).get("uncertainty", []):
        if e["quantity"] == quantity:
            return e["values"][0], e["extrapolated"], e["gci"]
    raise KeyError(quantity)


def tensor(quads, values):
    """Cell values of a structured 2D hex mesh as (y, x, grid)."""
    c = quads.mean(axis=1)
    xs = np.unique(np.round(c[:, 0], 12))
    ys = np.unique(np.round(c[:, 1], 12))
    ix = np.searchsorted(xs, np.round(c[:, 0], 12))
    iy = np.searchsorted(ys, np.round(c[:, 1], 12))
    g = np.full((ys.size, xs.size), np.nan)
    g[iy, ix] = values
    return xs, ys, g


def uniform(xs, ys, g, n=200, lo=(0.0, 0.0), hi=(1.0, 1.0)):
    """Bilinear resampling of a tensor grid onto a uniform one."""
    xu = np.linspace(lo[0], hi[0], n)
    yu = np.linspace(lo[1], hi[1], n)
    a = np.array([np.interp(xu, xs, row) for row in g])
    return xu, yu, np.array([np.interp(yu, ys, col) for col in a.T]).T


def cavity(out, runs):
    ras = [1e3, 1e4, 1e5, 1e6]
    fig = plt.figure(figsize=(16.0, 4.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 1, 1.25])
    nus = []
    for i, (ra, run) in enumerate(zip(ras, runs)):
        quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
        ax = fig.add_subplot(gs[i])
        pc = PolyCollection(quads, array=d["temperature"] - 300.0, cmap=RAINBOW,
                            edgecolors="face", linewidths=0.1)
        pc.set_clim(-0.5, 0.5)
        ax.add_collection(pc)
        xs, ys, ux = tensor(quads, d["velocity_x"])
        _, _, uy = tensor(quads, d["velocity_y"])
        xu, yu, uxu = uniform(xs, ys, ux)
        _, _, uyu = uniform(xs, ys, uy)
        ax.streamplot(xu, yu, uxu, uyu, color="w", density=1.1, linewidth=0.6, arrowsize=0.6)
        nu, nu_x, gci = estimate(run, "scalars.temperature.xmin.heat_flow")
        nus.append((ra, nu_x / DEPTH, gci))
        ax.set_title(f"Ra = 10^{int(np.log10(ra))}:  Nu = {nu_x / DEPTH:.3f}", fontsize=11)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_aspect("equal")
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.set_xlabel("x/L")
        if i == 0:
            ax.set_ylabel("y/L")
        else:
            ax.tick_params(labelleft=False)
    fig.colorbar(pc, ax=fig.axes[:4], label="(T − T_mean)/ΔT", shrink=0.85, pad=0.01)
    ax = fig.add_subplot(gs[4])
    r = np.array([n[0] for n in nus])
    v = np.array([n[1] for n in nus])
    ax.loglog(r, v, "o-", color="tab:red", lw=2, ms=7, label="Foehncraft, extrapolated (3 meshes)")
    ax.loglog(list(DE_VAHL_DAVIS), list(DE_VAHL_DAVIS.values()), "k+", ms=12, mew=1.6,
              label="de Vahl Davis 1983 (context)")
    for ra, nu, gci in nus:
        ax.annotate(f"GCI {100 * gci:.3f} %", (ra, nu), textcoords="offset points",
                    xytext=(-8, 9), ha="right", fontsize=8)
    ax.set_xlabel("Rayleigh number")
    ax.set_ylabel("mean Nusselt number of the hot wall")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("Mean Nu against Ra", fontsize=11)
    fig.suptitle("Air (Pr 0.71) in a square cavity, left wall hot, right wall cold, top and "
                 "bottom adiabatic: temperature and streamlines on 128² cells", fontsize=12)
    fig.savefig(os.path.join(out, "natural_convection_cavity.png"), dpi=100)
    plt.close(fig)
    return nus


def plate(out, run):
    from tools.verification.ostrach import profile, solve

    nu_ = 1e-4
    gbdt = 1.0
    k = 1e-4 * 1000.0 / 0.71
    _f2, h1 = solve(0.71)
    quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
    fig = plt.figure(figsize=(15.0, 5.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[0.55, 0.55, 1.2, 1.2])
    for i, (key, scale, lab, cmap) in enumerate([
        ("temperature", 1.0, "T − T∞ (K)", "inferno"),
        ("speed", 100.0, "speed (cm/s)", RAINBOW),
    ]):
        ax = fig.add_subplot(gs[i])
        val = d[key] - 300.0 if key == "temperature" else d[key] * scale
        pc = PolyCollection(quads * [5.0, 1.0], array=val, cmap=cmap, edgecolors="face",
                            linewidths=0.05)
        ax.add_collection(pc)
        ax.set_xlim(0.0, 0.5)
        ax.set_ylim(-0.2, 1.0)
        ax.set_aspect("equal")
        ax.plot([0.0, 0.0], [0.0, 1.0], color="tab:red", lw=5, solid_capstyle="butt")
        ax.set_xticks([0.0, 0.25, 0.5])
        ax.set_xticklabels(["0", "5", "10"])
        ax.set_xlabel("x (cm, 5 × stretched)")
        if i == 0:
            ax.set_ylabel("height (m)")
        else:
            ax.tick_params(labelleft=False)
        ax.set_title("Temperature" if i == 0 else "Speed", fontsize=11)
        fig.colorbar(pc, ax=ax, label=lab, shrink=0.9, pad=0.02, location="bottom")
    # local Nusselt number along the plate
    w = np.genfromtxt(os.path.join(run, "wall_plate.csv"), delimiter=",", names=True,
                      comments="#")
    o = np.argsort(w["y"])
    x = w["y"][o]
    q = w["temperature_flux"][o]
    nux = q * x / k
    gr = gbdt * x**3 / nu_**2
    ref = -h1 * (gr / 4.0) ** 0.25
    ax = fig.add_subplot(gs[2])
    ax.plot(x, nux, color="tab:red", lw=2.4, label="Foehncraft, finest mesh")
    ax.plot(x, ref, "k--", lw=1.4, label="similarity solution, −H′(0) = %.5f" % -h1)
    for name, xp in [("x025", 0.25), ("x050", 0.5), ("x075", 0.75)]:
        _q, qx, gci = estimate(run, f"wall_probes.{name}.temperature_flux")
        n = qx * xp / k
        r = -h1 * (gbdt * xp**3 / nu_**2 / 4.0) ** 0.25
        ax.plot(xp, n, "o", color="tab:blue", ms=7)
        ax.annotate(f"{100 * (n / r - 1):+.2f} %\nGCI {100 * gci:.3f} %", (xp, n),
                    textcoords="offset points", xytext=(6, -26), fontsize=8)
    ax.plot([], [], "o", color="tab:blue", label="extrapolated (3 meshes) at the wall probes")
    ax.set_xlabel("distance from the leading edge x (m)")
    ax.set_ylabel("local Nusselt number q x/(k ΔT)")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 35.0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("Local Nusselt number along the plate", fontsize=11)
    # profiles in similarity form at x = 0.5
    eta, fp, h = profile(0.71)
    xs, ys, vy = tensor(quads, d["velocity_y"])
    _, _, tt = tensor(quads, d["temperature"])
    j = np.argmin(np.abs(ys - 0.5))
    y0 = ys[j]
    scale = (gbdt * y0**3 / nu_**2 / 4.0) ** 0.25 / y0
    e = xs * scale
    uref = 2.0 * nu_ / y0 * np.sqrt(gbdt * y0**3 / nu_**2)
    ax = fig.add_subplot(gs[3])
    ax.plot(e, vy[j] / uref, "o", color="tab:red", ms=3.5, label="F′ = v/(2 ν Gr_x^½/x)")
    ax.plot(e, tt[j] - 300.0, "s", color="tab:blue", ms=3.5, label="H = (T − T∞)/ΔT")
    ax.plot(eta, fp, "k-", lw=1.2, label="similarity solution")
    ax.plot(eta, h, "k-", lw=1.2)
    ax.set_xlim(0.0, 8.0)
    ax.set_xlabel("η = (y/x) (Gr_x/4)^¼")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title(f"Across the layer at x = {y0:.3f} m (cell centres)", fontsize=11)
    fig.suptitle("Free convection along a vertical plate 1 K above the ambient, Pr 0.71, "
                 "Gr_L = 10⁸ (Ostrach 1953)", fontsize=12)
    fig.savefig(os.path.join(out, "natural_convection_plate.png"), dpi=100)
    plt.close(fig)


def mixed(out, up, down):
    gbdt, nu_, um = 0.01, 0.01, 5e-3
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True,
                            gridspec_kw={"width_ratios": [0.7, 1, 1]})
    quads, _, d = read(os.path.join(up, "fields.vtkhdf"))
    pc = PolyCollection(quads, array=d["velocity_y"] * 1000.0, cmap="RdBu_r",
                        edgecolors="k", linewidths=0.15)
    lim = np.abs(d["velocity_y"]).max() * 1000.0
    pc.set_clim(-lim, lim)
    axs[0].add_collection(pc)
    axs[0].set_xlim(0.0, 1.0)
    axs[0].set_ylim(0.0, 0.25)
    axs[0].set_aspect("equal")
    axs[0].set_xlabel("x/D")
    axs[0].set_ylabel("y (periodic)")
    axs[0].set_title("Vertical velocity, net flow up (32 × 8, perturbed)", fontsize=11)
    fig.colorbar(pc, ax=axs[0], label="v (mm/s)", location="bottom", shrink=0.9)
    xe = np.linspace(0.0, 1.0, 201)
    for ax, run, sign, title in [(axs[1], up, 1.0, "net flow up: reversal at the cold wall"),
                                 (axs[2], down, -1.0, "net flow down: reversal at the hot wall")]:
        quads, _, d = read(os.path.join(run, "fields.vtkhdf"))
        c = quads.mean(axis=1)
        xi = 2.0 * xe - 1.0
        v = gbdt / (48.0 * nu_) * (xi**3 - xi) + 6.0 * sign * um * xe * (1.0 - xe)
        ax.plot(xe, v * 1000.0, "k-", lw=1.5, label="exact (Poiseuille + buoyant cubic)")
        ax.plot(c[:, 0], d["velocity_y"] * 1000.0, "o", color="tab:red", ms=3, alpha=0.6,
                label="Foehncraft, every cell")
        ax.axhline(sign * um * 1000.0, color="tab:blue", ls=":", label="mean velocity U_m")
        ax.axhline(0.0, color="0.5", lw=0.8)
        ax.set_xlabel("x/D (hot wall at 0, cold wall at 1)")
        ax.set_ylabel("v (mm/s)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        orders = result(run).get("orders")
        extra = f"; observed order {orders[-1]:.2f}" if orders else ""
        ax.set_title(title + extra, fontsize=11)
    fig.suptitle("Fully developed mixed convection between a hot and a cold vertical wall, "
                 "Gr/Re = 200 (Aung & Worku 1986)", fontsize=12)
    fig.savefig(os.path.join(out, "mixed_convection_channel.png"), dpi=100)
    plt.close(fig)


def main():
    out = sys.argv[1]
    cav = sys.argv[2:6]
    plate_run, up, down = sys.argv[6:9]
    os.makedirs(out, exist_ok=True)
    for ra, nu, gci in cavity(out, cav):
        ref = DE_VAHL_DAVIS[ra]
        print(f"Ra {ra:.0e}: Nu {nu:.4f} (GCI {100 * gci:.3f} %), de Vahl Davis {ref} "
              f"({100 * (nu / ref - 1):+.2f} %)")
    plate(out, plate_run)
    mixed(out, up, down)


if __name__ == "__main__":
    main()
