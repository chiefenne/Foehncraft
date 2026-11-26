"""The drops take turbulence energy from the gas (D330; gallery and recipe
of its pull request): ECN Spray G, condition G1, on the k-ω SST model, run
without and with `[spray.coupling] turbulence = "dispersion-work"`.

spray_turbulence_two_way.png: left, a side view at `SNAPSHOT` ms with the
injector at the left and its axis horizontal: the turbulence kinetic energy
k of the gas on the plane of plumes 1 and 5 (the sector's plane `zmin`),
above the axis the run without the term, below it (mirrored) the run with
it, both on one logarithmic scale; right, the turbulence
kinetic energy of the gas of the whole injector over time, and the liquid
and vapour penetration of both runs against the measurements of the ECN
(the points of `tools/plot/spray_g_g1.py`).

The script prints the numbers of the gallery text: the turbulence kinetic
energy of both runs, the work the eddies did on the drops (`result.json`,
`gas_turbulence_work`, times the 16 sectors of the injector) and the
penetrations at 0.4, 0.6, 0.8 and 1.1 ms.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1
    target/release/foehn flow examples/cases/flow/spray_g_g1.toml \\
        --set spray.coupling.turbulence=dispersion-work --out output/spray_g_g1_work
    python3 -m tools.plot.spray_turbulence_two_way OUT_DIR output/spray_g_g1 output/spray_g_g1_work
"""
import glob
import json
import os
import sys
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.spray_g_g1 import (  # noqa: E402
    LIQUID_LOW,
    VAPOUR_PEN,
    cells,
    nearest,
    penetrations,
    plume_plane,
    step_of,
    times,
)
from tools.plot.spray_g_metrics import SECTOR  # noqa: E402

SNAPSHOT = 0.6  # ms, the side view
SECTORS = 360.0 / SECTOR  # sectors of the whole injector
REPORT = (0.4, 0.6, 0.8, 1.1)  # ms, the penetrations printed


def field_files(run):
    files = sorted(glob.glob(os.path.join(run, "fields_*.vtkhdf")), key=step_of)
    if not files:
        raise SystemExit(f"{run} has no fields_*.vtkhdf")
    return files


def cell_field(d, name, path):
    if name not in d:
        raise SystemExit(f"{path} has no cell field `{name}`; it has {', '.join(sorted(d))}")
    return d[name]


def volumes(xyz):
    """Volumes of the hexahedra (nodes in VTK order) from six tetrahedra."""
    v = np.zeros(len(xyz))
    for a, b, c in ((1, 2, 6), (1, 6, 5), (4, 5, 6), (4, 6, 7), (3, 7, 6), (3, 6, 2)):
        p0 = xyz[:, 0]
        v += np.einsum("ij,ij->i", xyz[:, a] - p0, np.cross(xyz[:, b] - p0, xyz[:, c] - p0)) / 6.0
    return np.abs(v)


def density(run):
    """The gas density of the run (the case keeps it constant)."""
    with open(os.path.join(run, "case.resolved.toml"), "rb") as f:
        rho = tomllib.load(f).get("fluid", {}).get("density")
    if not isinstance(rho, (int, float)):
        raise SystemExit(f"{run}/case.resolved.toml has no constant [fluid] density")
    return float(rho)


def tke_history(run):
    """(t ms, turbulence kinetic energy of the gas of the whole injector, mJ)."""
    tmap = times(run)
    rho = density(run)
    out = []
    vol = None
    for p in field_files(run):
        _, xyz, d = cells(p)
        if vol is None:
            vol = volumes(xyz)
        out.append((tmap[step_of(p)] * 1e3, rho * float((cell_field(d, "k", p) * vol).sum()) * SECTORS * 1e3))
    return np.array(out)


def work(run):
    """The work of the eddies on the gas of the whole injector (mJ, negative)
    from result.json, 0 without the term."""
    with open(os.path.join(run, "result.json")) as f:
        r = json.load(f)
    c = (r["solves"][-1].get("spray") or {}).get("coupling") or {}
    return c.get("gas_turbulence_work", 0.0) * SECTORS * 1e3


def at(series, t):
    return float(np.interp(t, series[:, 0], series[:, 1]))


def plot(out, run0, run1):
    k_lo, k_hi = 1e-2, 1e3
    fig = plt.figure(figsize=(12.0, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.0])
    ax = fig.add_subplot(gs[:, 0])
    pc = None
    for run, sign in ((run0, 1.0), (run1, -1.0)):
        ffile = nearest(field_files(run), SNAPSHOT, times(run))
        _, xyz, d = cells(ffile)
        quads, keep = plume_plane(xyz)
        quads = quads * np.array([1.0, sign])
        pc = PolyCollection(quads, array=np.clip(cell_field(d, "k", ffile)[keep], k_lo, k_hi), cmap=RAINBOW,
                            norm=LogNorm(k_lo, k_hi), edgecolors="face", linewidths=0.2)
        ax.add_collection(pc)
    t_snap = times(run0)[step_of(nearest(field_files(run0), SNAPSHOT, times(run0)))] * 1e3
    fig.colorbar(pc, ax=ax, label="turbulence kinetic energy k of the gas (m²/s²), plane of plumes 1 and 5",
                 shrink=0.85, pad=0.01, location="bottom", aspect=40)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.text(2.0, 36.0, "without the term", fontsize=10, va="top")
    ax.text(2.0, -36.0, "with dispersion-work: the drops take energy from k", fontsize=10, va="bottom")
    ax.set_xlim(0.0, 50.0)
    ax.set_ylim(-40.0, 40.0)
    ax.set_aspect("equal")
    ax.set_xlabel("distance from the injector tip along its axis (mm)")
    ax.set_ylabel("distance from the injector axis (mm)")
    ax.set_title(f"ECN Spray G (G1), k-ω SST, t = {t_snap:.2f} ms", fontsize=10.5)

    e0, e1 = tke_history(run0), tke_history(run1)
    a = fig.add_subplot(gs[0, 1])
    a.plot(e0[:, 0], e0[:, 1], color="tab:purple", lw=1.4, ls="--", label="without the term")
    a.plot(e1[:, 0], e1[:, 1], color="tab:purple", lw=2.0, label="with the term")
    a.axvline(0.78, color="0.6", lw=0.8, ls=":")
    a.set_ylabel("k of the gas (mJ)")
    a.set_xlim(0.0, 1.25)
    a.set_ylim(bottom=0.0)
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="upper left")
    a.set_title("Turbulence kinetic energy of the gas, whole injector", fontsize=10)

    (l0, v0), (l1, v1) = penetrations(run0), penetrations(run1)
    b = fig.add_subplot(gs[1, 1])
    m = np.array(LIQUID_LOW)
    b.plot(m[:, 0], m[:, 1], "o", color="k", ms=4, label="measured, 0.2e-3 mm³/mm²")
    b.plot(l0[:, 0], l0[:, 1], color="tab:red", lw=1.4, ls="--", label="without the term")
    b.plot(l1[:, 0], l1[:, 1], color="tab:red", lw=2.0, label="with the term")
    b.axvline(0.78, color="0.6", lw=0.8, ls=":")
    b.set_ylabel("liquid penetration (mm)")
    b.set_xlim(0.0, 1.25)
    b.set_ylim(bottom=0.0)
    b.grid(alpha=0.3)
    b.legend(fontsize=8, loc="lower right")
    b.set_title("Liquid along the injector axis (end of injection dotted)", fontsize=10)

    c = fig.add_subplot(gs[2, 1])
    m = np.array(VAPOUR_PEN)
    c.plot(m[:, 0], m[:, 1], "o", color="k", ms=4, label="measured (schlieren)")
    c.plot(v0[:, 0], v0[:, 1], color="tab:blue", lw=1.4, ls="--", label="without the term")
    c.plot(v1[:, 0], v1[:, 1], color="tab:blue", lw=2.0, label="with the term")
    c.axvline(0.78, color="0.6", lw=0.8, ls=":")
    c.set_xlabel("time after the start of injection (ms)")
    c.set_ylabel("vapour penetration (mm)")
    c.set_xlim(0.0, 1.25)
    c.set_ylim(bottom=0.0)
    c.grid(alpha=0.3)
    c.legend(fontsize=8, loc="lower right")
    c.set_title("Vapour (mass fraction 0.1 %)", fontsize=10)

    path = os.path.join(out, "spray_turbulence_two_way.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)
    for t in REPORT:
        print(f"{t:4.2f} ms  k {at(e0, t):7.2f} -> {at(e1, t):7.2f} mJ  liquid {at(l0, t):5.1f} -> "
              f"{at(l1, t):5.1f} mm (measured {at(np.array(LIQUID_LOW), t):5.1f})  vapour {at(v0, t):5.1f} -> "
              f"{at(v1, t):5.1f} mm (measured {at(np.array(VAPOUR_PEN), t):5.1f})")
    print(f"work of the eddies on the drops over the run: {work(run0):.2f} -> {work(run1):.2f} mJ")


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    plot(out, sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
