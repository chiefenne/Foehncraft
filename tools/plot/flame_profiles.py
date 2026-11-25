"""The laminar premixed methane-air flame of stage SF2b (gallery and recipe
of the chemistry in the flow cells).

flame_methane_air_structure.png: the GRI-Mech 3.0 flame at the end of its
run: the temperature and the OH mass fraction along the tube as colour
strips, then the temperature and the main species, and the radicals and
intermediates on a log axis, against the distance from the flame (where
the temperature is 1500 K).
flame_methane_air_speed.png: the burning velocity from the fuel
consumption over time for GRI-Mech 3.0 and Jones & Lindstedt, with the
measured value of stoichiometric methane and air at 298 K and 1 atm.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.flame_profiles OUT_DIR GRI_RUN JL4_RUN

with the output directories of

    target/release/foehn flow examples/cases/flow/flame_methane_air_gri30.toml --out output/gri
    target/release/foehn flow examples/cases/flow/flame_methane_air_jl4.toml --out output/jl4
"""
import os
import sys
from itertools import pairwise

import matplotlib

matplotlib.use("Agg")
import h5py
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.colormaps import RAINBOW

# measured burning velocity of stoichiometric CH4-air, 298 K, 1 atm: the
# value and standard deviation of the experimental target compiled for
# FFCM-1 (Smith, Tao, Wang 2016; heat-flux and counterflow measurements)
MEASURED = (0.360, 0.017)


def fields(run):
    """Cell centres along x and the cell data, sorted by x."""
    with h5py.File(os.path.join(run, "fields.vtkhdf"), "r") as f:
        g = f["VTKHDF/internal"]
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    x = np.array([pts[conn[a:b], 0].mean() for a, b in pairwise(off)])
    o = np.argsort(x)
    return x[o], {k: v[o] for k, v in data.items() if v.ndim == 1}


def speed(run, y_fuel=0.05519):
    """Time and burning velocity S_c = -m_F/(rho_u Y_F A) of every step."""
    with open(os.path.join(run, "chemistry.csv")) as f:
        rows = [l for l in f if not l.startswith("#")]
    head = rows[0].strip().split(",")
    v = np.array([[float(s) for s in l.split(",")] for l in rows[1:]])
    t = v[:, head.index("time")]
    m = -v[:, head.index("production_CH4")]
    # inlet: 300 K, 1 atm, the fresh mixture (N2 the rest); tube section 1e-4 m square
    w = 1.0 / (y_fuel / 16.043e-3 + 0.2200 / 31.998e-3 + (1 - y_fuel - 0.22) / 28.014e-3)
    rho = 101325.0 * w / (8.314462618 * 300.0)
    return t, m / (rho * y_fuel * 1e-8)


def structure(out, run):
    x, d = fields(run)
    t = d["temperature"]
    x0 = np.interp(1500.0, t, x)
    s = (x - x0) * 1e3
    fig = plt.figure(figsize=(12.0, 6.6))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.32, 0.32, 0.12, 3.0], hspace=0.08, wspace=0.42)
    for row, (key, label, scale) in enumerate(
        [("temperature", "T (K)", 1.0), ("Y_OH", "Y_OH (‰)", 1e3)]
    ):
        a = fig.add_subplot(gs[row, :])
        v = d[key] * scale
        a.pcolormesh(s, [0, 1], np.vstack([v, v]), cmap=RAINBOW, shading="nearest")
        a.set_xlim(-2.0, 3.0)
        a.set_yticks([])
        a.set_xticks([])
        a.text(
            -1.97, 0.5, f"{label}: {max(v.min(), 0.0):.0f} to {v.max():.4g}", va="center",
            color="white", fontsize=9,
        )
    a = fig.add_subplot(gs[3, 0])
    b = a.twinx()
    for k, c in [("CH4", "#1f4e9c"), ("O2", "#2e8b57"), ("H2O", "#d9661f"),
                 ("CO2", "#8b1a1a"), ("CO", "#9b59b6"), ("H2", "#17a2b8")]:
        a.plot(s, d[f"Y_{k}"], color=c, lw=1.8, label=k)
    b.plot(s, t, color="black", lw=2.2, ls="--", label="T")
    a.set_xlim(-2.0, 3.0)
    a.set_xlabel("distance from the 1500 K point (mm)")
    a.set_ylabel("mass fraction")
    b.set_ylabel("temperature (K)")
    a.legend(loc="center left", fontsize=8, frameon=False)
    b.legend(loc="center right", fontsize=8, frameon=False)
    a.set_title("GRI-Mech 3.0: temperature and main species", fontsize=10)
    a = fig.add_subplot(gs[3, 1])
    for k, c in [("OH", "#d62728"), ("H", "#1f77b4"), ("O", "#2ca02c"),
                 ("HO2", "#ff7f0e"), ("CH2O", "#9467bd"), ("CH3", "#8c564b"),
                 ("HCO", "#e377c2"), ("NO", "#17becf")]:
        a.semilogy(s, np.maximum(d[f"Y_{k}"], 1e-12), color=c, lw=1.6, label=k)
    a.set_xlim(-2.0, 3.0)
    a.set_ylim(1e-8, 2e-2)
    a.set_xlabel("distance from the 1500 K point (mm)")
    a.set_ylabel("mass fraction")
    a.legend(loc="lower right", fontsize=8, ncol=2, frameon=False)
    a.set_title("radicals, intermediates and NO", fontsize=10)
    p = os.path.join(out, "flame_methane_air_structure.png")
    fig.savefig(p, dpi=100, bbox_inches="tight")
    print(p)


def history(out, gri, jl4):
    fig, a = plt.subplots(figsize=(10.0, 4.6))
    m, sd = MEASURED
    a.axhspan(m - sd, m + sd, color="#999999", alpha=0.3, lw=0)
    a.axhline(m, color="#555555", lw=1.2, ls="--", label=f"measured {m:.2f} ± {sd:.3f} m/s")
    for run, label, c in [(gri, "GRI-Mech 3.0", "#1f4e9c"), (jl4, "Jones & Lindstedt 4-step", "#d9661f")]:
        t, v = speed(run)
        a.plot(t * 1e3, v, color=c, lw=1.8, label=label)
    a.set_xlabel("time (ms)")
    a.set_ylabel("burning velocity from fuel consumption (m/s)")
    a.set_ylim(0.0, 0.55)
    a.legend(loc="lower right", fontsize=9, frameon=False)
    p = os.path.join(out, "flame_methane_air_speed.png")
    fig.savefig(p, dpi=100, bbox_inches="tight")
    print(p)


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    out, gri, jl4 = argv
    os.makedirs(out, exist_ok=True)
    structure(out, gri)
    history(out, gri, jl4)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
