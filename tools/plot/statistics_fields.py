"""Statistics of a transient run (D71) behind a cylinder, gallery 110: the
time mean of u, the RMS of v and the resolved shear stress ⟨u'v'⟩ of the
vortex street at Re = 100, one panel each, from `statistics.vtkhdf`; with
the phase averages of the lift coefficient and of v at the wake probe from
`statistics_phase.csv` when given.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.statistics_fields OUT.png STATISTICS.vtkhdf [PHASE.csv]

e.g. for the cylinder of gallery 110 (recipe 110):

    target/release/foehn flow examples/cases/flow/cylinder_re100_statistics.toml
    python3 -m tools.plot.statistics_fields docs/figures/gallery/transient_statistics.png \\
        output/cylinder_re100_statistics/statistics.vtkhdf \\
        output/cylinder_re100_statistics/statistics_phase.csv
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from tools.plot.turbulence_fields import read  # noqa: E402
from tools.plot.colormaps import RAINBOW  # noqa: E402

U_MEAN = 1.0       # mean inlet velocity (m/s)
D = 0.1            # cylinder diameter (m)
CENTRE = (0.2, 0.2)
WINDOW = (0.0, 1.2, 0.0, 0.41)


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    out, path = sys.argv[1], sys.argv[2]
    phase = sys.argv[3] if len(sys.argv) == 4 else None
    quads, _, data = read(path)
    panels = [
        ("mean_velocity_x", "⟨u⟩ / Ū", "RdBu_r", 1.0 / U_MEAN, "signed"),
        ("rms_velocity_y", "v_rms / Ū", RAINBOW, 1.0 / U_MEAN, "pos"),
        ("covariance_velocity_xy", "⟨u'v'⟩ / Ū²", "PuOr_r", 1.0 / U_MEAN**2, "signed"),
    ]
    rows = len(panels) + (1 if phase else 0)
    fig, axes = plt.subplots(rows, 1, figsize=(8.4, 2.45 * rows), constrained_layout=True)
    for ax, (name, label, cmap, scale, kind) in zip(axes, panels):
        v = data[name] * scale
        x = (quads[:, :, 0] - CENTRE[0]) / D
        y = (quads[:, :, 1] - CENTRE[1]) / D
        q = np.stack([x, y], axis=-1)
        keep = (quads[:, :, 0].min(axis=1) < WINDOW[1])
        if kind == "signed":
            m = np.percentile(np.abs(v[keep]), 99.5)
            lo, hi = -m, m
        else:
            lo, hi = 0.0, np.percentile(v[keep], 99.5)
        pc = PolyCollection(q[keep], array=v[keep], cmap=cmap, edgecolors="face", linewidths=0.2)
        if name == "mean_velocity_x":
            lo, hi = np.percentile(v[keep], 0.5), np.percentile(v[keep], 99.5)
            pc.set_norm(TwoSlopeNorm(0.0, lo, hi))
        else:
            pc.set_clim(lo, hi)
        ax.add_collection(pc)
        ax.add_patch(plt.Circle((0, 0), 0.5, color="0.55", zorder=3))
        ax.set_xlim((WINDOW[0] - CENTRE[0]) / D, (WINDOW[1] - CENTRE[0]) / D)
        ax.set_ylim((WINDOW[2] - CENTRE[1]) / D, (WINDOW[3] - CENTRE[1]) / D)
        ax.set_aspect("equal")
        ax.set_ylabel("y/D")
        fig.colorbar(pc, ax=ax, label=label, shrink=0.9, pad=0.01)
    axes[len(panels) - 1].set_xlabel("x/D")
    if phase:
        ax = axes[-1]
        names, rows_ = read_phase(phase)
        ph = rows_[:, 0]
        for n, label, colour in (("coefficients.cl", "C_L", "tab:blue"),
                                 ("probes.wake.v", "v / Ū at the wake probe", "tab:green")):
            if n in names:
                c = rows_[:, names.index(n)]
                ax.plot(np.r_[ph, ph + 1], np.r_[c, c], "o-", ms=3, color=colour, label=label)
        ax.set_xlabel("phase (periods)")
        ax.set_ylabel("phase average")
        ax.set_xlim(0, 2)
        ax.grid(alpha=0.3)
        ax.legend(loc="center right", fontsize=8)
    fig.suptitle("Cylinder at Re = 100, t = 5 to 12 s: statistics of the vortex street (D71)", fontsize=11)
    fig.savefig(out, dpi=110)
    print(out)


def read_phase(path):
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    names = lines[0].strip().split(",")
    rows = np.array([[float(x) for x in ln.split(",")] for ln in lines[1:] if ln.strip()])
    return names, rows


if __name__ == "__main__":
    main()
