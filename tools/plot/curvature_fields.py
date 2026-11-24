"""The curvature correction of SST (D170) behind the backward-facing step,
gallery 113: the rotation function f_r1 that multiplies the production of
k and ω, the eddy-viscosity ratio with the correction, and its change
against plain SST, one panel each, with the reattachment points of both
runs. From two files written by `foehn post derive`.

    pip install matplotlib h5py numpy
    python3 -m tools.plot.curvature_fields OUT.png ON.vtkhdf OFF.vtkhdf [XR_ON XR_OFF]

e.g. (level 1 of `backstep_validation.toml`, 44 160 cells):

    S="--set turbulence.model=\\"k-omega-sst\\" --set turbulence.wall_treatment=\\"automatic\\"
       --set turbulence.ambient={k=1.5e-4,omega=0.54} --set verification.levels=[1]"
    target/release/foehn flow examples/cases/flow/backstep_validation.toml --out output/cc_on $S \\
        --set turbulence.curvature_correction=true --set run.restart_final=true
    target/release/foehn flow examples/cases/flow/backstep_validation.toml --out output/cc_off $S \\
        --set run.restart_final=true
    target/release/foehn post derive examples/cases/flow/backstep_validation.toml --out output/cc_on $S \\
        --set turbulence.curvature_correction=true --fields curvature-correction,turbulence-scales
    target/release/foehn post derive examples/cases/flow/backstep_validation.toml --out output/cc_off $S \\
        --fields turbulence-scales
    python3 -m tools.plot.curvature_fields docs/figures/gallery/curvature_correction.png \\
        output/cc_on/derived.vtkhdf output/cc_off/derived.vtkhdf 6.869 6.743
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from tools.plot import backstep_chart as bc
from tools.plot.turbulence_fields import read
from tools.plot.colormaps import RAINBOW


def main():
    if len(sys.argv) not in (4, 6):
        sys.exit(__doc__)
    out, on, off = sys.argv[1:4]
    xr = [float(v) for v in sys.argv[4:6]] if len(sys.argv) == 6 else None
    quads, _, d_on = read(on)
    _, _, d_off = read(off)
    fr1 = d_on["rotation_function"]
    r_on, r_off = d_on["viscosity_ratio"], d_off["viscosity_ratio"]
    centre = quads.mean(axis=1)
    keep = (centre[:, 0] > bc.WINDOW[0]) & (centre[:, 0] < bc.WINDOW[1]) & (centre[:, 1] < bc.WINDOW[3])
    rel = (r_on - r_off) / np.maximum(r_off, 1e-30)
    fig, axes = plt.subplots(3, 1, figsize=(10, 7.2), sharex=True, constrained_layout=True)
    pc = bc.panel(axes[0], quads, fr1, "RdBu_r", 0.0, 1.25)
    pc.set_norm(TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=1.25))
    fig.colorbar(pc, ax=axes[0], label="f_r1  (production factor)", shrink=0.9, pad=0.01)
    hi = np.percentile(r_on[keep], 99.5)
    pc = bc.panel(axes[1], quads, r_on, RAINBOW, 0.0, hi)
    fig.colorbar(pc, ax=axes[1], label="ν_t / ν  with the correction", shrink=0.9, pad=0.01)
    m = min(np.percentile(np.abs(rel[keep]), 99.5), 0.5)
    pc = bc.panel(axes[2], quads, rel, "PuOr_r", -m, m)
    fig.colorbar(pc, ax=axes[2], label="change of ν_t against plain SST", shrink=0.9, pad=0.01)
    if xr:
        for ax in axes:
            ax.plot([xr[1]], [0.04], marker="v", color="0.15", ms=7, zorder=5)
            ax.plot([xr[0]], [0.04], marker="v", color="#d62728", ms=7, zorder=5)
        axes[0].text(xr[0] + 0.15, 0.25, f"x_r/H {xr[0]:.2f} (corrected)", color="#d62728",
                     fontsize=8, zorder=6)
        axes[0].text(xr[1] - 0.15, 0.25, f"{xr[1]:.2f} (SST)", color="0.15", fontsize=8,
                     ha="right", zorder=6)
    axes[-1].set_xlabel("x/H")
    fig.suptitle("Backward-facing step, k-ω SST with the curvature correction of "
                 "Smirnov & Menter (2009)", fontsize=11)
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == "__main__":
    main()
