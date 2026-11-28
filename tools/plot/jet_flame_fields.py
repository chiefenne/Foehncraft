"""Turbulent jet flame of stage SF3 (gallery and recipe 351).

jet_flame_edc_fields.png: the piloted methane jet flame with the geometry
and inlets of Sandia Flame D, steady RANS with the eddy dissipation
concept, mirrored about the axis: temperature, CO mass fraction and heat
release rate.
jet_flame_edc_profiles.png: temperature and mass fractions along the axis
and the radial temperature at x/D = 15, 30 and 45, with the eddy
dissipation model (EDM_RUN, optional) dashed.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.jet_flame_fields OUT_DIR RUN [EDM_RUN]

with the output directories of

    target/release/foehn flow examples/cases/flow/jet_flame_edc_jl4.toml --out output/jet_flame
    target/release/foehn flow examples/cases/flow/jet_flame_edm_wd2.toml --out output/jet_flame_edm
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from tools.plot.colormaps import RAINBOW
from tools.plot.species_fields import cells

D = 0.0072  # main jet diameter (m)


def fields(out, run):
    x, r, d = cells(run)
    keep = x < 0.45
    x, r = x[keep], r[keep]
    d = {k: v[keep] for k, v in d.items() if v.ndim == 1}
    # mirrored about the axis, in x/D and r/D
    X = np.concatenate([x, x]) / D
    R = np.concatenate([r, -r]) / D
    tri = mtri.Triangulation(X, R)
    xy = np.c_[X, R][tri.triangles]
    span = np.ptp(xy, axis=1).max(axis=1)
    tri.set_mask(span > 3.0)
    fig = plt.figure(figsize=(12.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(3, 1)
    panels = [
        ("temperature", 1.0, RAINBOW, "temperature (K)"),
        ("Y_CO", 1.0, "inferno", "CO mass fraction (-)"),
        ("heat_release_rate", 1e-6, "magma", "heat release rate (MW/m³, log)"),
    ]
    for i, (k, s, cmap, label) in enumerate(panels):
        ax = fig.add_subplot(gs[i])
        v = np.concatenate([d[k], d[k]]) * s
        if k == "Y_CO":
            v = np.clip(v, 0.0, None)
        norm = None
        if k == "heat_release_rate":
            # the lips of the pilot burn 100 times faster than the jet
            v = np.clip(v, 0.1, None)
            norm = matplotlib.colors.LogNorm(vmin=0.1, vmax=v.max())
        tc = ax.tripcolor(tri, v, shading="gouraud", cmap=cmap, norm=norm)
        fig.colorbar(tc, ax=ax, label=label, shrink=0.9, pad=0.01)
        ax.set_xlim(0.0, 60.0)
        ax.set_ylim(-6.0, 6.0)
        ax.set_aspect("equal")
        ax.set_ylabel("r/D")
        if i == 2:
            ax.set_xlabel("x/D (D = 7.2 mm)")
        if i == 0:
            ax.set_title(
                "Piloted CH4/air jet flame (Sandia Flame D inlets, Re 22 400), steady RANS, "
                "k-epsilon, eddy dissipation concept with JL4; mirrored about the axis",
                fontsize=11,
            )
    path = os.path.join(out, "jet_flame_edc_fields.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)


def axis_cells(run):
    """The cells next to the axis, by x: x and the cell data."""
    x, r, d = cells(run)
    near = r < r.min() * 1.5
    o = np.argsort(x[near])
    return x[near][o], {k: v[near][o] for k, v in d.items() if v.ndim == 1}


def column(run, xd):
    """The column of cells nearest x/D = xd, by r: r and the temperature."""
    x, r, d = cells(run)
    xs = np.unique(x)
    s = np.abs(x - xs[np.argmin(np.abs(xs - xd * D))]) < 1e-9
    o = np.argsort(r[s])
    return r[s][o], d["temperature"][s][o]


def profiles(out, run, edm=None):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    x, a = axis_cells(run)
    s = x / D
    ax.plot(s, a["temperature"], color="tab:red", lw=2, label="T, EDC (JL4)")
    if edm:
        xe, e = axis_cells(edm)
        ax.plot(xe / D, e["temperature"], color="tab:red", lw=1.4, ls=":", label="T, EDM (WD2)")
    ax.set_xlabel("x/D on the axis")
    ax.set_ylabel("temperature (K)", color="tab:red")
    ax.set_xlim(0.0, 69.0)
    ax.set_ylim(250.0, 2400.0)
    ax.grid(alpha=0.3)
    ax2 = ax.twinx()
    for k, c in [("Y_CH4", "tab:green"), ("Y_O2", "tab:blue"), ("Y_CO", "tab:orange"),
                 ("Y_CO2", "tab:gray"), ("Y_H2O", "tab:cyan")]:
        ax2.plot(s, a[k], color=c, lw=1.5, ls="--", label=k)
    ax2.set_ylabel("mass fraction, EDC (-)")
    ax2.set_ylim(0.0, 0.36)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left", ncol=2)
    i = int(np.argmax(a["temperature"]))
    ax.set_title(
        f"Along the axis: EDC peak {a['temperature'][i]:.0f} K at x/D = {s[i]:.0f}", fontsize=11
    )
    for xd, c in [(15, "tab:purple"), (30, "tab:olive"), (45, "tab:brown")]:
        r, t = column(run, xd)
        bx.plot(r / D, t, color=c, lw=2, label=f"x/D = {xd}, EDC")
        if edm:
            r, t = column(edm, xd)
            bx.plot(r / D, t, color=c, lw=1.4, ls=":", label=f"x/D = {xd}, EDM")
    bx.set_xlabel("r/D")
    bx.set_ylabel("temperature (K)")
    bx.set_xlim(0.0, 5.0)
    bx.set_ylim(250.0, 2400.0)
    bx.grid(alpha=0.3)
    bx.legend(fontsize=8, ncol=2 if edm else 1)
    bx.set_title("Radial temperature", fontsize=11)
    path = os.path.join(out, "jet_flame_edc_profiles.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)


def main():
    out, run = sys.argv[1:3]
    edm = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(out, exist_ok=True)
    fields(out, run)
    profiles(out, run, edm)


if __name__ == "__main__":
    main()
