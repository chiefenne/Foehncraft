"""Gas-phase chemistry of stage SF2 (gallery and recipe of the reactor).

chemistry_reactor_traces.png: stoichiometric methane and air at 1 atm from
1400 K in the constant-pressure adiabatic reactor with GRI-Mech 3.0 and the
two global mechanisms: the temperature against the equilibrium of each
mechanism, and the main species and radicals of GRI-Mech 3.0.
chemistry_ignition_map.png: the temperature of the GRI-Mech 3.0 reactor
over time and start temperature (1100 to 1900 K), with the ignition delay.

    pip install matplotlib numpy
    python3 -m tools.plot.reactor_traces OUT_DIR [FOEHN]

FOEHN defaults to target/release/foehn; the script runs the three reactor
cases of examples/cases/chemistry and the start-temperature sweep itself.
"""
import csv
import json
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tools.plot.colormaps import RAINBOW

CASES = "examples/cases/chemistry"
MECHS = [
    ("gri30", "GRI-Mech 3.0 (53 species, 325 reactions)", "#1f4e9c"),
    ("jl4", "Jones & Lindstedt 1988 (4 steps)", "#d9661f"),
    ("wd2", "Westbrook & Dryer 1981 (2 steps)", "#2e8b57"),
]


def run(foehn, case, csv_path):
    out = subprocess.run(
        [foehn, "--json", "chemistry", "reactor", case, "--csv", csv_path],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return json.loads(out)


def read(path):
    with open(path) as f:
        r = csv.reader(f)
        head = next(r)
        v = np.array([[float(x) for x in row] for row in r])
    return head, v


def traces(out, foehn):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.0, 4.8))
    for key, label, colour in MECHS:
        case = os.path.join(CASES, f"reactor_methane_air_{key}.toml")
        p = os.path.join(out, f"reactor_{key}.csv")
        res = run(foehn, case, p)
        head, v = read(p)
        t, T = v[:, 0], v[:, 1]
        a.semilogx(t[1:] * 1e3, T[1:], color=colour, lw=2.0, label=label)
        a.axhline(res["equilibrium_temperature"], color=colour, lw=1.0, ls="--")
        if key == "gri30":
            names = ["CH4", "O2", "H2O", "CO2", "CO", "H2", "OH", "H", "O", "NO"]
            cmap = plt.get_cmap(RAINBOW)
            for i, n in enumerate(names):
                x = v[:, head.index("X_" + n)]
                b.loglog(t[1:] * 1e3, np.maximum(x[1:], 1e-12), lw=1.8,
                         color=cmap(i / (len(names) - 1)), label=n)
    a.set_xlabel("time (ms)")
    a.set_ylabel("temperature (K)")
    a.set_xlim(1e-4, 100)
    a.set_title("CH₄–air, φ = 1, 1 atm, from 1400 K\n(dashed: equilibrium of each mechanism)")
    a.legend(loc="upper left", fontsize=8.5, frameon=False)
    a.grid(True, which="both", alpha=0.25)
    b.set_xlabel("time (ms)")
    b.set_ylabel("mole fraction")
    b.set_xlim(1e-4, 100)
    b.set_ylim(1e-8, 1)
    b.set_title("GRI-Mech 3.0: species and radicals")
    b.legend(loc="lower left", ncol=2, fontsize=8.5, frameon=False)
    b.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "chemistry_reactor_traces.png"), dpi=100)
    plt.close(fig)


def ignition_map(out, foehn):
    base = open(os.path.join(CASES, "reactor_methane_air_gri30.toml")).read()
    base = base.split("[reactor.check]")[0]
    temps = np.arange(1100.0, 1901.0, 10.0)
    grid_t = np.logspace(-6, -1, 400)
    rows, delays = [], []
    for t0 in temps:
        case = os.path.join(out, f"sweep_{int(t0)}.toml")
        with open(case, "w") as f:
            f.write(base.replace("temperature = 1400.0", f"temperature = {t0}"))
        p = os.path.join(out, f"sweep_{int(t0)}.csv")
        res = run(foehn, case, p)
        _, v = read(p)
        rows.append(np.interp(grid_t, v[1:, 0], v[1:, 1]))
        delays.append(res["ignition_delay"])
    T = np.array(rows)
    fig, ax = plt.subplots(figsize=(10.0, 5.6))
    m = ax.pcolormesh(grid_t * 1e3, temps, T, cmap=RAINBOW, shading="nearest", rasterized=True)
    d = np.array([x if x is not None else np.nan for x in delays])
    ax.plot(d * 1e3, temps, color="white", lw=2.2, ls="--", label="ignition delay (largest dT/dt)")
    ax.set_xscale("log")
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("start temperature (K)")
    ax.set_title("GRI-Mech 3.0, CH₄–air φ = 1, 1 atm: temperature of the reactor")
    ax.legend(loc="upper right", frameon=False, labelcolor="white")
    c = fig.colorbar(m, ax=ax)
    c.set_label("temperature (K)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "chemistry_ignition_map.png"), dpi=100)
    plt.close(fig)
    for t0, x in zip(temps, delays):
        print(f"T0 = {t0:.0f} K: ignition delay {x}")


def main(argv):
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    out = argv[0]
    foehn = argv[1] if len(argv) > 1 else "target/release/foehn"
    os.makedirs(out, exist_ok=True)
    traces(out, foehn)
    ignition_map(out, foehn)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
