"""Evaporation and boiling by the volume of fluid (L3 stage MP2; gallery and
recipe of the MP2 pull request).

checks: vof_phase_change_checks.png, the vapour layer of the Stefan and the
sucking problem and the radius of Scriven's bubble against the exact
solutions (Welch & Wilson 2000; Scriven 1959), one panel each.

film: vof_film_boiling.png, the temperature with the interface (α = 0.5)
of a film boiling run at four times, and the Nusselt number of the wall
(heat.csv) over time against the correlation of Berenson (1961).

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/vof_stefan.toml --out output/vof_stefan
    target/release/foehn flow examples/cases/flow/vof_sucking.toml --out output/vof_sucking
    target/release/foehn flow examples/cases/flow/vof_scriven.toml --out output/vof_scriven
    target/release/foehn flow examples/cases/flow/vof_film_boiling.toml --out output/vof_film_boiling
    python3 -m tools.plot.vof_phase_change checks OUT_DIR output/vof_stefan output/vof_sucking output/vof_scriven
    python3 -m tools.plot.vof_phase_change film OUT_DIR output/vof_film_boiling
"""
import glob
import math
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read

# the exact constants of the example cases (their headers derive them)
STEFAN = dict(length=1.0e-3, area=4.0e-10, lam=0.06774198562683174, alpha=1.978878679494709e-5)
SUCKING = dict(length=1.2e-3, area=9.0e-12, beta=3.3893598713587387, alpha=1.6758045748921316e-7)
SCRIVEN = dict(r0=2.0e-6, r1=5.0e-4, angle=2.0, beta=4.0669305415576495, alpha=1.6758045748921316e-7)
# film boiling: λ₀ = √(σ/(g Δρ)), the vapour conductivity, the wall superheat
FILM = dict(l0=0.007230164703039606, k_v=1.0, dT=5.0, rho_v=5.0, rho_l=200.0, mu_v=0.005, h_lv=1.0e4, g=9.81)


def vof_csv(run, need=("time", "volume")):
    path = os.path.join(run, "vof.csv")
    if not os.path.exists(path):
        raise SystemExit(f"{run} has no vof.csv (not a volume-of-fluid run?)")
    csv = np.genfromtxt(path, delimiter=",", names=True, comments=None, skip_header=2)
    for c in need:
        if c not in csv.dtype.names:
            raise SystemExit(f"vof.csv has no column `{c}`; it has {', '.join(csv.dtype.names)}")
    return csv


def field(d, name):
    if name in d:
        return d[name]
    have = ", ".join(sorted(d))
    raise SystemExit(f"the fields file has no cell field `{name}`; it has {have}")


def layer(run, c, exact):
    """t, computed and exact thickness of a planar vapour layer."""
    csv = vof_csv(run)
    t = csv["time"]
    return t, c["length"] - csv["volume"] / c["area"], exact(t)


def bubble(run):
    """t, computed and exact radius of the bubble on the sphere sector."""
    c = SCRIVEN
    csv = vof_csv(run)
    t2 = math.tan(math.radians(c["angle"]))
    om = 4.0 * t2 * t2 / (1.0 + 2.0 * t2 * t2) ** 1.5
    total = om * (c["r1"] ** 3 - c["r0"] ** 3) / 3.0
    vapour = total - csv["volume"] + om * c["r0"] ** 3 / 3.0
    t = csv["time"]
    return t, np.cbrt(3.0 * vapour / om), 2.0 * c["beta"] * np.sqrt(c["alpha"] * t)


def checks(out, runs):
    if len(runs) != 3:
        raise SystemExit("checks needs three runs: Stefan, sucking, Scriven")
    s, k = STEFAN, SUCKING
    panels = [
        ("Stefan problem (wall 10 K above T_sat)", "vapour layer δ (mm)",
         layer(runs[0], s, lambda t: 2.0 * s["lam"] * np.sqrt(s["alpha"] * t))),
        ("Sucking problem (liquid 2 K superheated)", "vapour layer δ (mm)",
         layer(runs[1], k, lambda t: 2.0 * k["beta"] * np.sqrt(k["alpha"] * t))),
        ("Scriven's bubble (liquid 1.25 K superheated)", "bubble radius R (mm)", bubble(runs[2])),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 6.8), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[2.2, 1.0]))
    for j, (title, ylabel, (t, comp, ex)) in enumerate(panels):
        a, b = axes[0, j], axes[1, j]
        a.plot(t * 1e3, ex * 1e3, color="black", lw=2.0, label="exact")
        step = max(1, len(t) // 25)
        a.plot(t[::step] * 1e3, comp[::step] * 1e3, "o", color="tab:red", ms=5, mfc="white", mew=1.6,
               label="Foehncraft, THINC + interface energy")
        a.set_title(title, fontsize=11)
        a.set_ylabel(ylabel)
        a.grid(alpha=0.3)
        a.legend(fontsize=9, loc="upper left")
        dev = 100.0 * (comp / ex - 1.0)
        b.plot(t * 1e3, dev, color="tab:blue", lw=1.6)
        b.axhline(0.0, color="black", lw=0.8)
        b.set_xlabel("t (ms)")
        b.set_ylabel("deviation (%)")
        b.set_ylim(-5.0, 5.0)
        b.grid(alpha=0.3)
    fig.suptitle("Evaporation of water at 1 atm against exact solutions (volume of fluid, density ratio 1604)",
                 fontsize=12)
    fig.savefig(os.path.join(out, "vof_phase_change_checks.png"), dpi=100)
    plt.close(fig)


def berenson():
    f = FILM
    g = f["rho_v"] * (f["rho_l"] - f["rho_v"]) * f["g"] * f["h_lv"] * f["l0"] ** 3 / (f["k_v"] * f["mu_v"] * f["dT"])
    return 0.425 * g ** 0.25


def wall_nusselt(run):
    """t and Nu = q λ₀/(λ_v ΔT) of the heated wall from heat.csv: the
    conduction into the domain is that of the wall (the liquid at the
    outlet is at T_sat), over the wall area."""
    path = os.path.join(run, "heat.csv")
    if not os.path.exists(path):
        raise SystemExit(f"{run} has no heat.csv (a run without [energy]?)")
    csv = np.genfromtxt(path, delimiter=",", names=True, comments=None, skip_header=2)
    names = csv.dtype.names
    for c in ("time", "conduction"):
        if c not in names:
            raise SystemExit(f"heat.csv has no column `{c}`; it has {', '.join(names)}")
    quads, _, _ = read(os.path.join(run, "fields.vtkhdf"))
    with h5py.File(os.path.join(run, "fields.vtkhdf"), "r") as h:
        from tools.plot.turbulence_fields import internal
        z = internal(h)["Points"][:, 2]
    area = np.ptp(quads[:, :, 0]) * np.ptp(z)
    f = FILM
    return csv["time"], csv["conduction"] / area * f["l0"] / (f["k_v"] * f["dT"])


def film(out, run):
    snaps = sorted(glob.glob(os.path.join(run, "fields_*.vtkhdf")))
    if len(snaps) < 4:
        raise SystemExit(f"{run} has {len(snaps)} field snapshots; 4 or more are needed ([output] fields_every)")
    pick = [snaps[i] for i in np.linspace(len(snaps) // 4, len(snaps) - 1, 4).round().astype(int)]
    fig = plt.figure(figsize=(14.0, 8.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 4, height_ratios=[3.0, 1.0])
    pc = None
    for i, path in enumerate(pick):
        quads, _, d = read(path)
        with h5py.File(path, "r") as h:
            time = float(h.attrs["time"]) if "time" in h.attrs else float("nan")
        ax = fig.add_subplot(gs[0, i])
        pc = PolyCollection(quads, array=field(d, "temperature"), cmap=RAINBOW, edgecolors="face",
                            linewidths=0.1)
        pc.set_clim(500.0, 505.0)
        ax.add_collection(pc)
        cx = quads[:, :, 0].mean(axis=1)
        cy = quads[:, :, 1].mean(axis=1)
        xs, ys = np.unique(cx.round(9)), np.unique(cy.round(9))
        if len(xs) * len(ys) == len(cx):
            a = np.full((len(ys), len(xs)), np.nan)
            a[np.searchsorted(ys, cy.round(9)), np.searchsorted(xs, cx.round(9))] = field(d, "volume_fraction")
            ax.contour(xs, ys, a, levels=[0.5], colors="white", linewidths=1.6)
        ax.set_xlim(quads[:, :, 0].min(), quads[:, :, 0].max())
        ax.set_ylim(quads[:, :, 1].min(), quads[:, :, 1].max())
        ax.set_aspect("equal")
        ax.set_title(f"t = {time:.2f} s", fontsize=11)
        ax.set_xlabel("x (m)")
        if i == 0:
            ax.set_ylabel("y (m)")
        else:
            ax.set_yticklabels([])
    fig.colorbar(pc, ax=fig.axes[:4], location="right", label="temperature (K), interface white",
                 shrink=0.8, pad=0.01)
    b = fig.add_subplot(gs[1, :])
    t, nu = wall_nusselt(run)
    b.plot(t, nu, color="tab:red", lw=1.4, label="Foehncraft, wall mean")
    late = t > 0.25 * t[-1]
    if late.any():
        b.axhline(nu[late].mean(), color="tab:red", ls="--", lw=1.2,
                  label=f"time mean after {0.25 * t[-1]:.2f} s: {nu[late].mean():.2f}")
    nb = berenson()
    b.axhline(nb, color="black", lw=1.6, label=f"Berenson 1961: {nb:.2f}")
    b.set_xlabel("t (s)")
    b.set_ylabel("Nu = q λ₀/(λ_v ΔT)")
    b.set_ylim(0.0, max(4.0 * nb, 1.2 * np.percentile(nu, 95)))
    b.grid(alpha=0.3)
    b.legend(fontsize=9, loc="upper right")
    fig.suptitle("Film boiling on a horizontal wall 5 K above T_sat (volume of fluid, near-critical fluid)",
                 fontsize=12)
    fig.savefig(os.path.join(out, "vof_film_boiling.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) < 4 or sys.argv[1] not in ("checks", "film"):
        raise SystemExit(__doc__)
    what, out, runs = sys.argv[1], sys.argv[2], sys.argv[3:]
    os.makedirs(out, exist_ok=True)
    if what == "checks":
        checks(out, runs)
    else:
        film(out, runs[0])


if __name__ == "__main__":
    main()
