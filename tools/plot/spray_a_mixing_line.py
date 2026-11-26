"""The adiabatic mixing line of n-dodecane and nitrogen at Spray A
conditions and the liquid length of Siebers (1999) from it (ECN Spray A,
stage SA1 part 2, checks VA.4 and VA.5; gallery and recipe of the pull
request).

spray_a_mixing_line.png: left, the mixing temperature over the mixture
fraction of 363 K n-dodecane in 900 K nitrogen at 6 MPa, frozen (one phase
of the Peng-Robinson equation), in equilibrium (Peng-Robinson flash) and
ideal (ideal gas, Raoult's law), the two-phase region shaded and the
saturation mixture fraction of each line marked; middle, the liquid mass
fraction of the mixture in the two-phase region; right, the liquid length
of Siebers' scaling law with the saturation mixture fraction of each model
against the measured Spray A liquid lengths at the four conditions.

    pip install matplotlib numpy
    cargo run --release -p foehncraft-particles --example spray_a_mixing_line -- output/spray_a_mixing_line
    python3 -m tools.plot.spray_a_mixing_line output/spray_a_mixing_line docs/figures/gallery
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
    ln = table(os.path.join(run, "line.csv"), ["model", "Z", "T", "liquid", "rho", "split"])
    ll = table(
        os.path.join(run, "liquid_length.csv"),
        ["model", "T_gas", "rho_gas", "Z_sat", "L_mm", "theta_deg", "measured_mm"],
    )
    cmap = plt.get_cmap(RAINBOW)
    col = {"frozen": cmap(0.12), "equilibrium": cmap(0.88), "ideal": cmap(0.62)}
    label = {
        "frozen": "frozen (one phase, Peng-Robinson)",
        "equilibrium": "equilibrium (Peng-Robinson flash)",
        "ideal": "ideal gas, Raoult's law",
    }

    fig, axs = plt.subplots(1, 3, figsize=(16.0, 4.9), constrained_layout=True)
    a, b, c = axs
    zsat = {}
    for m in ["frozen", "ideal", "equilibrium"]:
        k = ln["model"] == m
        z, t, s = ln["Z"][k], ln["T"][k], ln["split"][k] == 1
        a.plot(z, t, color=col[m], lw=2.4 if m == "equilibrium" else 1.8, label=label[m],
               ls="--" if m == "frozen" else "-")
        zsat[m] = z[np.argmax(s)] if s.any() else None
    for m in ["equilibrium", "ideal"]:
        k = (ll["model"] == m) & (ll["T_gas"] == 900) & (ll["rho_gas"] == 22.8) & (ll["theta_deg"] == 21.5)
        zs = float(ll["Z_sat"][k][0])
        kk = ln["model"] == m
        ts = np.interp(zs, ln["Z"][kk], ln["T"][kk])
        a.plot(zs, ts, "o", color=col[m], ms=8, mec="k", mew=0.8)
        right = m == "equilibrium"
        a.annotate(f"Z$_{{sat}}$ = {zs:.3f} ({m})", (zs, ts), textcoords="offset points",
                   xytext=(10, 6) if right else (-10, -4), ha="left" if right else "right", fontsize=9)
    ke = (ln["model"] == "equilibrium") & (ln["split"] == 1)
    a.axvspan(ln["Z"][ke].min(), ln["Z"][ke].max(), color="0.85", alpha=0.6, lw=0)
    a.text(0.62, 860, "two phases\n(equilibrium)", fontsize=9, ha="center", va="top")
    a.set_xlim(0, 1)
    a.set_ylim(340, 920)
    a.set_xlabel("mixture fraction Z (fuel mass fraction from the injector)")
    a.set_ylabel("mixing temperature (K)")
    a.set_title("363 K n-dodecane into 900 K N$_2$ at 6 MPa", fontsize=10)
    a.grid(alpha=0.3)
    a.legend(fontsize=8.5, loc="lower left")

    for m in ["ideal", "equilibrium"]:
        k = (ln["model"] == m) & (ln["split"] == 1)
        b.plot(ln["Z"][k], ln["liquid"][k], color=col[m], lw=2.4 if m == "equilibrium" else 1.8,
               label=label[m])
    b.plot([0, 1], [0, 1], color="k", lw=0.8, ls=":")
    b.text(0.36, 0.43, "Z: all the fuel liquid", rotation=33, fontsize=8)
    b.set_xlim(0.25, 1.0)
    b.set_ylim(0, 1.0)
    b.set_xlabel("mixture fraction Z")
    b.set_ylabel("mass fraction of the liquid phase")
    b.set_title("Liquid phase: the high-pressure liquid takes up N$_2$", fontsize=10)
    b.grid(alpha=0.3)
    b.legend(fontsize=8.5, loc="upper left")

    conds = []
    for t, r in zip(ll["T_gas"], ll["rho_gas"]):
        if (t, r) not in conds:
            conds.append((t, r))
    xs = np.arange(len(conds))
    for j, m in enumerate(["ideal", "equilibrium"]):
        for th, mk in [(21.5, "o"), (18.7, "^")]:
            y = []
            for t, r in conds:
                k = (ll["model"] == m) & (ll["T_gas"] == t) & (ll["rho_gas"] == r) & (ll["theta_deg"] == th)
                y.append(float(ll["L_mm"][k][0]))
            c.plot(xs + (j - 0.5) * 0.18, y, mk, color=col[m], ms=8, mec="k", mew=0.6,
                   label=f"Siebers, {m}, {th:g}°", ls="none")
    meas = []
    for t, r in conds:
        k = (ll["T_gas"] == t) & (ll["rho_gas"] == r)
        meas.append(float(ll["measured_mm"][k][0]))
    c.bar(xs, meas, width=0.62, color="0.88", edgecolor="0.4", label="measured (DBI, ECN)", zorder=0)
    c.set_xticks(xs, [f"{t:g} K\n{r:g} kg/m$^3$" for t, r in conds], fontsize=9)
    c.set_ylabel("liquid length (mm)")
    c.set_ylim(0, 24)
    c.set_title("Siebers (1999) liquid length with our Z$_{sat}$ against Spray A", fontsize=10)
    c.grid(alpha=0.3, axis="y")
    c.legend(fontsize=8, loc="upper left")
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "spray_a_mixing_line.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: python3 -m tools.plot.spray_a_mixing_line RUN_DIR OUT_DIR")
    plot(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
