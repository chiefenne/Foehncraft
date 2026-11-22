"""Skin friction along the inner wall of the U-duct (gallery 204): the
grid study of the rounded bend, the curvature correction of SST and the
sharp-edged variant.

    python3 -m tools.plot.ubend_friction OUT.png GCI_DIR CC_DIR SHARP_DIR

GCI_DIR holds run_0 .. run_2 of examples/studies/ubend_gci.py, CC_DIR run_1
of the same with `--set turbulence.curvature_correction=true`, SHARP_DIR
run_1 of `ubend_gci.py --sharp`. The abscissa is the distance along the
inner wall from the point where it leaves the turn (negative upstream:
along the inner arc of radius 0.5 H, or along the divider end).
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "examples/studies"))
import ubend_gci as g  # noqa: E402


def curve(d: Path, sharp: bool):
    g.SHARP = sharp
    w = g.along(d, "inner")
    exit_ = g.station(0.0, 0.5)
    xs, cf = [], []
    for st, tau, *_ in w:
        if st < exit_:
            s = -(exit_ - st) * (1.0 if sharp else 0.5)   # radians times the inner radius
        else:
            s = st - exit_
        xs.append(s)
        cf.append(2.0 * tau)                              # rho = 1, U_b = 1
    return xs, cf


def main(out, gci, cc, sharp) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.2), dpi=200)
    for lv, c, lab in [(0, "#9ecae1", "8 832 cells"), (1, "#4292c6", "35 328 cells"), (2, "#08306b", "141 312 cells")]:
        x, y = curve(Path(gci) / f"run_{lv}", False)
        ax.plot(x, y, color=c, lw=1.6, label=f"rounded, {lab}")
    x, y = curve(Path(cc) / "run_1", False)
    ax.plot(x, y, color="tab:orange", lw=1.6, ls="--", label="rounded, curvature correction, 35 328 cells")
    x, y = curve(Path(sharp) / "run_1", True)
    ax.plot(x, y, color="tab:red", lw=1.4, label="sharp edges, 50 400 cells")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.axvline(0.0, color="0.6", lw=0.8, ls=":")
    ax.text(-0.05, 0.97, "turn", ha="right", va="top", transform=ax.get_xaxis_transform(), color="0.4")
    ax.text(0.05, 0.97, "outlet leg", ha="left", va="top", transform=ax.get_xaxis_transform(), color="0.4")
    ax.set_xlim(-1.0, 5.0)
    ax.set_ylim(-0.008, 0.012)
    ax.set_xlabel("distance along the inner wall from the end of the turn (H)")
    ax.set_ylabel("C_f = 2 τ_w / ρ U_b² (along the flow)")
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out)


if __name__ == "__main__":
    main(*sys.argv[1:5])
