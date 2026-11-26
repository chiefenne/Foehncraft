"""The real-fluid gas film of an n-dodecane droplet in nitrogen at high
pressure (ECN Spray A, stage SA2, D313; gallery and recipe of the pull
request).

spray_a_film.png: left, the dense-fluid factors eta/eta0 and lambda/lambda0
of Chung et al. (1988) for nitrogen at 6 and 15 MPa with the Peng-Robinson
density (lines) against the NIST reference ratios to 0.1 MPa (points);
middle, the film of a droplet in 900 K nitrogen by the 1/3 rule, the
real-fluid properties over the ideal-gas ones against the droplet
temperature at 6 and 15 MPa; right, a 10 um droplet of 363 K n-dodecane
in still 900 K nitrogen at 6 MPa with either film: squared diameter and
temperature.

    pip install matplotlib numpy
    cargo run --release -p foehncraft-particles --example spray_a_film -- output/spray_a_film
    python3 -m tools.plot.spray_a_film output/spray_a_film docs/figures/gallery
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
    n2 = table(os.path.join(run, "nitrogen.csv"), ["p_MPa", "T", "eta_factor", "lambda_factor"])
    nist = table(os.path.join(run, "nist.csv"), ["p_MPa", "T", "eta_ratio", "lambda_ratio"])
    film = table(os.path.join(run, "film.csv"), ["p_MPa", "T_drop", "rho", "eta", "lambda", "cp", "D"])
    drop = table(os.path.join(run, "droplet.csv"), ["film", "t_ms", "d_um", "T"])
    cmap = plt.get_cmap(RAINBOW)

    fig, axs = plt.subplots(1, 3, figsize=(16.0, 4.9), constrained_layout=True)
    a, b, c = axs
    for p, ls in [(6, "-"), (15, "--")]:
        k = n2["p_MPa"] == p
        kn = nist["p_MPa"] == p
        a.plot(n2["T"][k], n2["eta_factor"][k], ls, color=cmap(0.12), lw=2.0, label=f"viscosity, Chung, {p} MPa")
        a.plot(n2["T"][k], n2["lambda_factor"][k], ls, color=cmap(0.85), lw=2.0,
               label=f"conductivity, Chung, {p} MPa")
        fill = p == 6
        a.plot(nist["T"][kn], nist["eta_ratio"][kn], "o", color=cmap(0.12), ms=7, mec="k", mew=0.7,
               mfc=cmap(0.12) if fill else "white", label=f"viscosity, NIST, {p} MPa")
        a.plot(nist["T"][kn], nist["lambda_ratio"][kn], "s", color=cmap(0.85), ms=7, mec="k", mew=0.7,
               mfc=cmap(0.85) if fill else "white", label=f"conductivity, NIST, {p} MPa")
    a.set_xlabel("temperature (K)")
    a.set_ylabel("property at p over property at 0.1 MPa")
    a.set_title("Nitrogen: dense-fluid factors of Chung et al. (1988) and NIST", fontsize=10)
    a.set_xlim(280, 1220)
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="upper right", ncol=2)

    names = [("rho", "density", 0.05), ("eta", "viscosity", 0.25), ("lambda", "conductivity", 0.45),
             ("cp", "c$_p$", 0.65), ("D", "diffusivity", 0.9)]
    for p, ls in [(6, "-"), (15, "--")]:
        k = film["p_MPa"] == p
        for col, lab, f in names:
            b.plot(film["T_drop"][k], film[col][k], ls, color=cmap(f), lw=2.0,
                   label=lab if p == 6 else None)
        tend = float(film["T_end"][k][0]) if "T_end" in film.dtype.names else None
        if tend:
            b.axvline(tend, color="0.5", lw=0.8, ls=ls)
            b.text(tend + (2 if p == 6 else -3), 0.905, f"critical\nmixing,\n{p} MPa", fontsize=7.5,
                   ha="left" if p == 6 else "right", va="bottom", color="0.35")
    b.axhline(1.0, color="k", lw=0.8, ls=":")
    b.set_xlabel("droplet temperature (K)")
    b.set_ylabel("real-fluid film over ideal-gas film")
    b.set_title("Film by the 1/3 rule in 900 K N$_2$ (solid 6 MPa, dashed 15 MPa)", fontsize=10)
    b.set_ylim(0.9, 1.16)
    b.set_xlim(365, 668)
    b.grid(alpha=0.3)
    b.legend(fontsize=8.5, loc="lower left", ncol=2)

    c2 = c.twinx()
    for name, f, ls in [("real-fluid", 0.12, "-"), ("ideal-gas", 0.62, "--")]:
        k = drop["film"] == name
        t, d, tp = drop["t_ms"][k], drop["d_um"][k], drop["T"][k]
        c.plot(t, (d / 10.0) ** 2, ls, color=cmap(f), lw=2.2, label=f"(d/d$_0$)$^2$, {name} film")
        c2.plot(t, tp, ls, color=cmap(0.88 if name == "real-fluid" else 0.75), lw=1.6,
                label=f"temperature, {name} film")
    c.set_xlabel("time (ms)")
    c.set_ylabel("(d/d$_0$)$^2$")
    c2.set_ylabel("droplet temperature (K)")
    c.set_ylim(0, 1.05)
    c.set_xlim(-0.01, 0.84)
    c2.set_ylim(350, 660)
    c.set_title("10 µm n-dodecane droplet, 363 K, in still 900 K N$_2$ at 6 MPa", fontsize=10)
    c.grid(alpha=0.3)
    h1, l1 = c.get_legend_handles_labels()
    h2, l2 = c2.get_legend_handles_labels()
    c.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right", bbox_to_anchor=(0.995, 0.55))
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "spray_a_film.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: python3 -m tools.plot.spray_a_film RUN_DIR OUT_DIR")
    plot(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
