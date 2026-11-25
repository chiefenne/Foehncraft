"""Water and steam states of IAPWS-IF97, stage SF7a (gallery and recipe
of the materials).

steam_ph_chart.png: the pressure-enthalpy chart of water and steam from
`foehn steam`, coloured by the density, with the saturation dome, the
isotherms and the states of the two steam cases.
steam_superheater_channel.png: steam at 4 MPa superheated between two
plates at 760 K: temperature with streamlines, density, and the IF97
properties along the temperature at 4 MPa.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.steam_states OUT_DIR SUPERHEATER_RUN [FOEHN]

with the output directory of

    target/release/foehn flow examples/cases/flow/steam_superheater_channel.toml --out output/superheater

FOEHN defaults to target/release/foehn.
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read

FOEHN = "target/release/foehn"


def steam(**kw):
    """One `foehn steam` query: its JSON (`state` or `saturation`), or None
    outside IF97."""
    args = [FOEHN, "--json", "steam"]
    for k, v in kw.items():
        args += [f"--{k}", repr(float(v))]
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    return json.loads(r.stdout)


def many(queries):
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as ex:
        return list(ex.map(lambda q: steam(**q), queries))


def need(d, name, path):
    if name not in d:
        sys.exit(f"{path}: no field `{name}`; the file has: {', '.join(sorted(d))}")
    return d[name]


def ph_chart(out):
    hs = np.linspace(50e3, 4200e3, 150)
    ps = np.logspace(np.log10(1e3), np.log10(100e6), 110)
    res = many([{"pressure": p, "enthalpy": h} for p in ps for h in hs])
    rho = np.array([r["state"]["density"] if r else np.nan for r in res]).reshape(len(ps), len(hs))
    fig, ax = plt.subplots(figsize=(9.6, 6.2), constrained_layout=True)
    m = ax.pcolormesh(hs / 1e3, ps / 1e6, rho, cmap=RAINBOW, norm=LogNorm(5e-3, 1100.0),
                      shading="auto")
    fig.colorbar(m, ax=ax, label="density (kg/m³)")
    # the saturation dome
    pd = np.logspace(np.log10(1e3), np.log10(22.0e6), 80)
    sat = many([{"pressure": p} for p in pd])
    hl = [s["saturation"]["liquid"]["enthalpy"] / 1e3 for s in sat]
    hv = [s["saturation"]["vapour"]["enthalpy"] / 1e3 for s in sat]
    ax.plot(hl, pd / 1e6, "k-", lw=1.6)
    ax.plot(hv, pd / 1e6, "k-", lw=1.6)
    ax.plot([hl[-1], hv[-1]], [pd[-1] / 1e6] * 2, "k-", lw=1.6)
    for x in (0.25, 0.5, 0.75):
        q = many([{"pressure": p, "quality": x} for p in pd])
        ax.plot([s["state"]["enthalpy"] / 1e3 for s in q], pd / 1e6, "k:", lw=0.8)
    # isotherms
    pi = np.logspace(np.log10(1e3), np.log10(100e6), 90)
    for t in (300.0, 400.0, 500.0, 600.0, 647.096, 700.0, 800.0, 1000.0):
        st = many([{"pressure": p, "temperature": t} for p in pi])
        h = np.array([s["state"]["enthalpy"] / 1e3 if s else np.nan for s in st])
        ax.plot(h, pi / 1e6, "w-", lw=0.9)
        k = np.nanargmin(np.abs(pi - 0.3e6)) if t > 450.0 else np.nanargmin(np.abs(pi - 40e6))
        if np.isfinite(h[k]):
            ax.text(h[k], pi[k] / 1e6, "T_c" if t == 647.096 else f"{t:.0f} K", color="w", fontsize=8, ha="center",
                    va="bottom", rotation=80)
    # the states of the cases at 4 MPa
    for t, lab, off, ha in ((540.0, "channel inlet", (-6, -14), "right"),
                            (600.0, "duct inlet", (-4, 7), "right"),
                            (750.0, "duct outlet", (5, 7), "left")):
        s = steam(pressure=4e6, temperature=t)["state"]
        ax.plot(s["enthalpy"] / 1e3, 4.0, "o", ms=6, mfc="w", mec="k")
        ax.annotate(lab, (s["enthalpy"] / 1e3, 4.0), xytext=off, textcoords="offset points",
                    fontsize=8.5, color="w", ha=ha)
    ax.set_yscale("log")
    ax.set_xlim(hs[0] / 1e3, hs[-1] / 1e3)
    ax.set_ylim(ps[0] / 1e6, ps[-1] / 1e6)
    ax.set_xlabel("specific enthalpy (kJ/kg)")
    ax.set_ylabel("pressure (MPa)")
    ax.set_title("Water and steam of IAPWS-IF97 from `foehn steam`: density over (h, p), the "
                 "saturation dome with x = 0.25, 0.5, 0.75, isotherms", fontsize=10)
    fig.savefig(os.path.join(out, "steam_ph_chart.png"), dpi=110)
    plt.close(fig)


def panel(ax, quads, values, label, fig):
    pc = PolyCollection(quads * 1e3, array=values, cmap=RAINBOW, edgecolors="face",
                        linewidths=0.2)
    ax.add_collection(pc)
    ax.set_xlim(0.0, 60.0)
    ax.set_ylim(0.0, 10.0)
    ax.set_aspect("equal")
    fig.colorbar(pc, ax=ax, label=label, fraction=0.025, pad=0.01)
    return pc


def superheater(out, run):
    path = os.path.join(run, "fields.vtkhdf")
    quads, _, d = read(path)
    fig = plt.figure(figsize=(12.0, 5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[2.3, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, quads, need(d, "temperature", path), "temperature (K)", fig)
    c = quads.mean(axis=1)
    order = np.lexsort((c[:, 0], c[:, 1]))
    nx, ny = 240, 40
    gx = c[order, 0].reshape(ny, nx)[0] * 1e3
    gy = c[order, 1].reshape(ny, nx)[:, 0] * 1e3
    ux = need(d, "velocity_x", path)[order].reshape(ny, nx)
    uy = need(d, "velocity_y", path)[order].reshape(ny, nx)
    ax.streamplot(gx, gy, ux, uy, color="w", density=(2.4, 0.5), linewidth=0.6, arrowsize=0.6)
    ax.set_ylabel("y (mm)")
    ax.set_title("steam at 540 K and 6 mm/s between plates at 760 K, 4 MPa", fontsize=10)
    ax = fig.add_subplot(gs[1, 0])
    panel(ax, quads, need(d, "density", path), "density (kg/m³)", fig)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("density from IF97: the heated steam expands and speeds up", fontsize=10)
    # the IF97 properties at 4 MPa
    ax = fig.add_subplot(gs[:, 1])
    ts = np.linspace(524.0, 800.0, 70)
    st = [r["state"] for r in many([{"pressure": 4e6, "temperature": t} for t in ts])]
    s0 = st[0]
    for key, lab, col in (("density", "ρ", "#2b6cb0"), ("specific_heat", "c_p", "#c53030"),
                          ("viscosity", "μ", "#d69e2e"), ("conductivity", "k", "#2f855a"),
                          ("prandtl", "Pr", "#6b46c1")):
        ax.plot(ts, [s[key] / s0[key] for s in st], "-", color=col, lw=1.6, label=lab)
    ax.axvline(523.5075, color="0.4", ls=":", lw=1.0)
    ax.text(527.0, 0.62, "T_s(4 MPa)", fontsize=8, color="0.3", rotation=90)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("value / saturated vapour value")
    ax.set_title("IF97 steam at 4 MPa, as the case sees it", fontsize=10)
    ax.legend(fontsize=9)
    fig.savefig(os.path.join(out, "steam_superheater_channel.png"), dpi=110)
    plt.close(fig)


def main():
    global FOEHN
    out, run = sys.argv[1:3]
    if len(sys.argv) > 3:
        FOEHN = sys.argv[3]
    ph_chart(out)
    superheater(out, run)


if __name__ == "__main__":
    main()
