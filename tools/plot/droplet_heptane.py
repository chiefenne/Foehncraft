"""Single n-heptane droplets in still nitrogen at 0.1 MPa (D320, D321;
gallery and recipe of the liquid library pull request).

droplet_heptane.png: left, d^2/d0^2 against t/d0^2 at five gas temperatures
with the 1/3 rule film (solid) and the far-field properties (dashed);
middle, the droplet temperature of the same runs (the wet bulb); right, the
evaporation rate K against the gas temperature for both film states, the
microgravity fit of Chauveau et al. (2019) and the fibre-free extrapolation
of Chauveau et al. (2008).

    pip install matplotlib numpy
    cargo run --release -p foehncraft-flow --example droplet_heptane -- output/droplet_heptane
    python3 -m tools.plot.droplet_heptane output/droplet_heptane docs/figures/gallery
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
    h = table(os.path.join(run, "history.csv"), ["T_gas", "film", "t_over_d0sq_s_per_mm2", "d2_over_d02", "T_drop"])
    r = table(os.path.join(run, "rates.csv"), ["T_gas", "K_one_third", "K_far_field", "K_microgravity", "K_fibre_free"])
    cmap = plt.get_cmap(RAINBOW)
    temps = sorted(set(h["T_gas"]))
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.9), constrained_layout=True)
    for i, tg in enumerate(temps):
        c = cmap(0.08 + 0.84 * i / max(len(temps) - 1, 1))
        for film, ls in [("one-third", "-"), ("far-field", "--")]:
            m = (h["T_gas"] == tg) & (h["film"] == film)
            lab = f"{tg:.0f} K" if film == "one-third" else None
            axs[0].plot(h["t_over_d0sq_s_per_mm2"][m], h["d2_over_d02"][m], ls, color=c, lw=1.8, label=lab)
            axs[1].plot(h["t_over_d0sq_s_per_mm2"][m], h["T_drop"][m], ls, color=c, lw=1.8, label=lab)
    axs[0].set_xlabel("t / d0² (s/mm²)")
    axs[0].set_ylabel("d² / d0²")
    axs[0].set_ylim(0, 1.05)
    axs[0].set_xlim(left=0)
    axs[0].legend(title="N2, 0.1 MPa (solid 1/3 rule, dashed far field)", fontsize=8, title_fontsize=8)
    axs[1].set_xlabel("t / d0² (s/mm²)")
    axs[1].set_ylabel("droplet temperature (K)")
    axs[1].set_xlim(left=0)
    t = r["T_gas"]
    axs[2].plot(t, r["K_one_third"], "o-", color=cmap(0.15), lw=2, label="1/3 rule film (default)")
    axs[2].plot(t, r["K_far_field"], "s--", color=cmap(0.85), lw=2, label="far-field properties")
    axs[2].plot(t, r["K_fibre_free"], "-", color="0.25", lw=2, label="measured, fibre-free extrapolation (2008)")
    axs[2].plot(t, r["K_microgravity"], ":", color="0.25", lw=2.2, label="measured, microgravity fit (2019)")
    axs[2].set_xlabel("gas temperature (K)")
    axs[2].set_ylabel("K (mm²/s)")
    axs[2].set_ylim(bottom=0)
    axs[2].legend(fontsize=8)
    for a in axs:
        a.grid(alpha=0.3)
    fig.suptitle("n-heptane droplets, d0 = 0.37 mm, in still nitrogen at 0.1 MPa", fontsize=11)
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "droplet_heptane.png"), dpi=100)


if __name__ == "__main__":
    plot(sys.argv[1] if len(sys.argv) > 1 else "output/droplet_heptane",
         sys.argv[2] if len(sys.argv) > 2 else "docs/figures/gallery")
