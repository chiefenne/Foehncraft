"""Derived fields of the turbulence models (D72) behind the backward-facing
step, gallery 108: the Reynolds shear stress u'v' of the Boussinesq
estimate, the eddy-viscosity ratio and the turbulence time scale, one panel
each, from a file written by `foehn post derive`.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.turbulence_fields OUT.png DERIVED.vtkhdf

e.g. for the k-ε-ζ-f₀ step of gallery 74:

    target/release/foehn flow examples/cases/flow/backstep_zeta_f0.toml \\
        --set run.restart_final=true
    target/release/foehn post derive examples/cases/flow/backstep_zeta_f0.toml \\
        --fields turbulence-scales,reynolds-stress,turbulence-wall
    python3 -m tools.plot.turbulence_fields docs/figures/gallery/turbulence_derived.png \\
        output/backstep_zeta_f0/derived.vtkhdf
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from tools.plot.vtkhdf import internal  # noqa: E402

from tools.plot import backstep_chart as bc  # noqa: E402
from tools.plot.colormaps import UNIFORM  # noqa: E402

U_REF = 1.05   # u at x/H = -4, y/H = 5 (the reference of the step's C_p)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out, path = sys.argv[1], sys.argv[2]
    quads, _, data = read(path)
    panels = [
        ("reynolds_stress_xy", "u'v' / U_ref²  (−2 ν_t S_xy)", UNIFORM, lambda v: v / U_REF**2, True),
        ("viscosity_ratio", "ν_t / ν", "magma", lambda v: v, False),
        ("time_scale", "T U_ref / H  (time scale of the model)", "viridis", lambda v: v * U_REF, False),
    ]
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 2.4 * len(panels)), sharex=True,
                             constrained_layout=True)
    for ax, (name, label, cmap, f, signed) in zip(axes, panels):
        v = f(data[name])
        keep = (quads[:, :, 0].mean(axis=1) > bc.WINDOW[0]) & (quads[:, :, 0].mean(axis=1) < bc.WINDOW[1]) \
            & (quads[:, :, 1].mean(axis=1) < bc.WINDOW[3])
        if signed:
            m = np.percentile(np.abs(v[keep]), 99.5)
            lo, hi = -m, m
        else:
            lo, hi = 0.0, np.percentile(v[keep], 99.5)
        pc = bc.panel(ax, quads, v, cmap, lo, hi)
        fig.colorbar(pc, ax=ax, label=label, shrink=0.9, pad=0.01)
    axes[-1].set_xlabel("x/H")
    fig.suptitle("Backward-facing step, k-ε-ζ-f₀: derived fields of the turbulence model", fontsize=11)
    fig.savefig(out, dpi=110)
    print(out)


def read(path):
    """As backstep_chart.fields, from any file."""
    import h5py
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        assert np.all(np.diff(off) == 8), "hexahedra only"
        cells = conn.reshape(-1, 8)
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    xyz = pts[cells]
    zmin = xyz[:, :, 2].min(axis=1, keepdims=True)
    quads = []
    for c, zm in zip(xyz, zmin):
        q = c[np.abs(c[:, 2] - zm[0]) < 1e-12][:, :2]
        m = q.mean(axis=0)
        quads.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
    return np.array(quads), None, data


if __name__ == "__main__":
    main()
