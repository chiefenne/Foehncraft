"""The surface of an n-dodecane droplet in nitrogen at high pressure (ECN
Spray A, stage SA1; gallery and recipe of the SA1 pull request).

spray_a_surface.png: left, the fuel vapour mole fraction at the droplet
surface over Raoult's law against the droplet temperature at 1, 3, 6 and
8 MPa (Peng-Robinson, D311), each curve ending at the critical mixing
temperature or the boiling point; middle, the latent heat over its value at
saturation (D312) and the nitrogen dissolved in the surface liquid at the
same pressures; right, the bubble pressure over the dissolved nitrogen on
the seven isotherms of the measured data of Garcia-Cordova et al. (2011)
with 21 of their points.

    pip install matplotlib numpy
    cargo run --release -p foehncraft-particles --example spray_a_surface -- output/spray_a_surface
    python3 -m tools.plot.spray_a_surface output/spray_a_surface docs/figures/gallery
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.colormaps import RAINBOW


def table(path, need):
    if not os.path.exists(path):
        raise SystemExit(f"{path} is missing: run the example first (see the docstring)")
    d = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    for c in need:
        if c not in d.dtype.names:
            raise SystemExit(f"{path} has no column `{c}`; it has {', '.join(d.dtype.names)}")
    return d


def plot(run, out):
    s = table(os.path.join(run, "surface.csv"), ["p_MPa", "T", "y_fuel", "y_raoult", "x_N2", "latent_ratio", "end"])
    b = table(os.path.join(run, "bubble.csv"), ["T", "x_N2", "p_kPa"])
    m = table(os.path.join(run, "measured.csv"), ["T", "p_kPa", "x_N2"])
    cmap = plt.get_cmap(RAINBOW)
    pressures = sorted(set(s["p_MPa"]))
    pcol = {p: cmap(0.08 + 0.84 * i / max(len(pressures) - 1, 1)) for i, p in enumerate(pressures)}

    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.9), constrained_layout=True)
    a, c, e = axs
    for p in pressures:
        k = s["p_MPa"] == p
        t = s["T"][k]
        f = s["y_fuel"][k] / s["y_raoult"][k]
        a.plot(t, f, color=pcol[p], lw=2.2, label=f"{p:g} MPa")
        end = s["end"][k][-1]
        a.plot(t[-1], f[-1], "o", color=pcol[p], ms=6, mec="k", mew=0.6)
        a.annotate(
            f"{t[-1]:.1f} K ({end})", (t[-1], f[-1]), textcoords="offset points",
            xytext=(-6, 6), ha="right", fontsize=8, color="k",
        )
        c.plot(t, s["latent_ratio"][k], color=pcol[p], lw=2.2, label=f"L/L_sat, {p:g} MPa")
        c.plot(t, s["x_N2"][k] * 100.0 / 40.0, color=pcol[p], lw=1.4, ls="--")
    a.axhline(1.0, color="k", lw=1.0, ls=":")
    a.text(310, 1.02, "Raoult's law", fontsize=9)
    a.set_xlim(300, 670)
    a.set_ylim(0.8, 3.2)
    a.set_xlabel("droplet temperature (K)")
    a.set_ylabel("surface fuel vapour fraction / Raoult")
    a.set_title("n-dodecane in N$_2$: fuel at the droplet surface", fontsize=10)
    a.grid(alpha=0.3)
    a.legend(fontsize=9, loc="upper left")

    c.set_xlim(300, 670)
    c.set_ylim(0.0, 1.05)
    c.set_xlabel("droplet temperature (K)")
    c.set_ylabel("latent heat / latent heat at saturation (solid)")
    c2 = c.twinx()
    c2.set_ylim(0.0, 1.05 * 40.0)
    c2.set_ylabel("dissolved N$_2$ in the surface liquid, mol % (dashed)")
    c.set_title("Latent heat falls to zero at the critical mixing state", fontsize=10)
    c.grid(alpha=0.3)

    temps = sorted(set(b["T"]))
    tcol = {t: cmap(0.05 + 0.9 * i / max(len(temps) - 1, 1)) for i, t in enumerate(temps)}
    for t in temps:
        k = b["T"] == t
        e.plot(b["x_N2"][k], b["p_kPa"][k] / 1e3, color=tcol[t], lw=1.8, label=f"{t:g} K")
        km = np.isclose(m["T"], t)
        e.plot(m["x_N2"][km], m["p_kPa"][km] / 1e3, "o", color=tcol[t], ms=6, mec="k", mew=0.6)
    e.set_xlim(0.0, 0.32)
    e.set_ylim(0.0, 20.0)
    e.set_xlabel("dissolved N$_2$ mole fraction")
    e.set_ylabel("bubble pressure (MPa)")
    e.set_title("Peng-Robinson (lines) against measured bubble points (dots)", fontsize=10)
    e.grid(alpha=0.3)
    e.legend(fontsize=8, loc="upper left", ncol=2)
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "spray_a_surface.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: python3 -m tools.plot.spray_a_surface RUN_DIR OUT_DIR")
    plot(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
