"""Near-body refinement with hanging-node faces, L13 stage 2f (gallery and
recipe entry `refine`).

immersed_refine.png: the Schäfer–Turek cylinder on a D/10 box mesh refined
twice near the cylinder (the speed, the cells outlined), the temperature
error of the manufactured conduction solution on a perturbed mesh refined
near a circle (a zoom on its upper left), its convergence against the same case without refinement, and
the lift of the refined cylinder against the uniform cut meshes of gallery
313. Read from the output directories of

    target/release/foehn flow examples/cases/flow/immersed_refine_cylinder_re20.toml --out CYL
    target/release/foehn flow examples/cases/flow/immersed_refine_mms.toml --out MMS \\
        --set 'verification.levels=[1]'

    pip install matplotlib h5py numpy
    python3 -m tools.plot.immersed_refine OUT_DIR CYL MMS
"""
import os
import sys

import numpy as np

from tools.plot.immersed_cut_plate import read_section

# L2 temperature errors of immersed_refine_mms.toml on 16, 32, 64 input
# cells across, refined (the case file) and without [immersed.refine]
MMS_H = np.array([1 / 16, 1 / 32, 1 / 64])
MMS_REFINED = np.array([5.5409e-2, 1.4578e-2, 3.7617e-3])
MMS_UNIFORM = np.array([1.0813e-1, 2.4101e-2, 5.8615e-3])
# lift of the cut cylinder: uniform box meshes D/10, D/20, D/40 (gallery
# 313) and the refined D/10 mesh of the case file (cells after the cut)
UNIFORM_CELLS = np.array([8944, 35780, 143088])
UNIFORM_CL = np.array([0.0078, 0.0103, 0.0106])
UNIFORM_CD = np.array([5.554, 5.588, 5.595])


def field(path, name):
    import h5py
    with h5py.File(path, "r") as f:
        return f["VTKHDF"]["internal"]["CellData"][name][:]


def benchmark(out):
    """C_D, C_L, Δp and the cells of a run, from its result.json."""
    import json
    r = json.load(open(os.path.join(out, "result.json")))
    s = r["solves"][0]
    v = [x["value"] for x in s["benchmark"]["values"]]
    return v, s["cells"]


def plot(out, cyl, mms):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    from tools.plot.colormaps import RAINBOW

    fig = plt.figure(figsize=(14.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0])

    # the cylinder: speed with the cells outlined
    ax = fig.add_subplot(gs[0, 0])
    polys, cells, d = read_section(os.path.join(cyl, "fields.vtkhdf"))
    x0, x1, y0, y1 = 0.02, 0.62, 0.0, 0.41
    keep = [i for i, q in enumerate(polys) if q[:, 0].min() < x1 and q[:, 0].max() > x0]
    pc = PolyCollection([polys[i] for i in keep], array=d["speed"][cells[keep]], cmap=RAINBOW,
                        edgecolors=(0, 0, 0, 0.45), linewidths=0.3)
    pc.set_clim(0.0, 0.42)
    ax.add_collection(pc)
    t = np.linspace(0, 2 * np.pi, 200)
    ax.fill(0.2 + 0.05 * np.cos(t), 0.2 + 0.05 * np.sin(t), color="0.35")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    v, n = benchmark(cyl)
    ax.set_title(f"cylinder at Re 20 on a D/10 box mesh refined twice near it: {n} cells\n"
                 f"C_D {v[0]:.3f}, C_L {v[1]:.5f}, Δp {v[2]:.4f} Pa", fontsize=11)
    fig.colorbar(pc, ax=ax, shrink=0.85, pad=0.01, label="speed (m/s)")

    # the manufactured solution: the error field across the levels
    ax = fig.add_subplot(gs[0, 1])
    polys, cells, d = read_section(os.path.join(mms, "fields.vtkhdf"))
    err = d["temperature_error"][cells]
    lim = np.abs(err).max()
    keep = [i for i, q in enumerate(polys) if q[:, 0].min() < 0.6 and q[:, 1].max() > 0.5]
    polys, err = [polys[i] for i in keep], err[keep]
    pc = PolyCollection(polys, array=err, cmap="RdBu_r", edgecolors=(0, 0, 0, 0.35), linewidths=0.25)
    pc.set_clim(-lim, lim)
    ax.add_collection(pc)
    ax.plot(0.45 + 0.28 * np.cos(t), 0.55 + 0.28 * np.sin(t), "k--", lw=1.0, label="cut surface")
    ax.set_xlim(0.04, 0.6)
    ax.set_ylim(0.5, 0.98)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("conduction MMS on a perturbed 32 × 32 mesh refined near a circle\n"
                 "(upper left part): T − T_exact, no jump at the faces between levels", fontsize=11)
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(pc, ax=ax, shrink=0.85, pad=0.01, label="T − T_exact (K)")

    # convergence of the manufactured solution
    ax = fig.add_subplot(gs[1, 1])
    ax.loglog(MMS_H, MMS_UNIFORM, "o-", color="#30638e", label="without refinement")
    ax.loglog(MMS_H, MMS_REFINED, "s-", color="#d1495b", label="refined twice near the circle")
    ax.loglog(MMS_H, MMS_REFINED[0] * (MMS_H / MMS_H[0]) ** 2, "k--", lw=0.9, label="second order")
    o = np.log2(MMS_REFINED[:-1] / MMS_REFINED[1:])
    ax.set_xlabel("width of the input cells (m)")
    ax.set_ylabel("L2 temperature error (K)")
    ax.set_title(f"manufactured solution: order {o[0]:.2f}, {o[1]:.2f} with hanging faces", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # the lift against the uniform meshes
    ax = fig.add_subplot(gs[1, 0])
    ax.axhspan(0.0104, 0.0110, color="0.85", label="reference range (Schäfer & Turek 1996)")
    ax.semilogx(UNIFORM_CELLS, UNIFORM_CL, "o-", color="#30638e", label="uniform box, cut (D/10, D/20, D/40)")
    ax.semilogx([n], [v[1]], "s", color="#d1495b", ms=9, label="D/10 refined to D/40 near the cylinder")
    for c, l, lab in zip(UNIFORM_CELLS, UNIFORM_CL, ["D/10", "D/20", "D/40"]):
        ax.annotate(lab, (c, l), textcoords="offset points", xytext=(6, -12), fontsize=8)
    ax.set_xlabel("cells")
    ax.set_ylabel("lift coefficient C_L")
    ax.set_ylim(0.0070, 0.0115)
    ax.set_title("lift of the cut cylinder: the refined mesh matches D/40\n"
                 "with a seventh of its cells", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    path = os.path.join(out, "immersed_refine.png")
    fig.savefig(path, dpi=100)
    print(path)


def main(argv):
    if len(argv) < 4:
        sys.exit(__doc__)
    import matplotlib

    matplotlib.use("Agg")
    plot(argv[1], argv[2], argv[3])


if __name__ == "__main__":
    main(sys.argv)
