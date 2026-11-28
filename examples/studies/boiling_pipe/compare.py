"""Subcooled boiling in a heated tube (stage MP4, D359): the cross-section
averaged void fraction of the two-fluid model against the measurements of
Bartolomei et al. (1982), as void fraction over the equilibrium quality.

    foehn flow examples/cases/flow/boiling_pipe_bartolomei_b7.toml --out OUT/b7
    (likewise b5, b9, b11, b12, b13)
    python3 examples/studies/boiling_pipe/compare.py OUT [--chart FILE] [--fields NAME]

Reads OUT/<case>/fields.vtkhdf (or the file named by --fields) of every
case present and the measured points of
examples/studies/boiling_pipe/data/bartolomei1982_void.csv. The averaged
void of an axial station is Σ α r / Σ r over its cells (the volumes of a
wedge grow with r); the equilibrium quality at x follows from the heat
balance, x_e = (c_p (T_in − T_sat) + 4 q x/(G D))/h_fg, with the mean
c_p of the case file. Prints per case the mean and largest deviation from
the measured points (the run interpolated at their x_e) and the quality
where the averaged void first passes 0.02 in the run and in the data.
`--chart FILE` writes a figure (SVG or PNG by the suffix).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "bartolomei1982_void.csv"
D = 0.012

# p (Pa), G, q, T_in, T_sat, h_fg (IAPWS-IF97), mean c_p of the case file
CASES = {
    "b5": (6.81e6, 998.0, 0.44e6, 521.0, 557.124, 1.517425e6, 5046.2),
    "b7": (6.84e6, 961.0, 1.13e6, 466.0, 557.420, 1.515479e6, 4773.3),
    "b9": (7.01e6, 996.0, 1.98e6, 434.0, 559.077, 1.504487e6, 4676.7),
    "b11": (14.74e6, 1847.0, 0.77e6, 598.0, 613.910, 1.018263e6, 7217.5),
    "b12": (14.75e6, 2123.0, 1.13e6, 583.0, 613.964, 1.017591e6, 6673.9),
    "b13": (14.70e6, 2014.0, 1.72e6, 545.0, 613.693, 1.020945e6, 5917.0),
}


def read_data() -> dict[str, list[tuple[float, float]]]:
    out: dict[str, list[tuple[float, float]]] = {}
    with DATA.open() as f:
        rows = [r for r in f if not r.startswith("#")]
    for r in csv.DictReader(rows):
        key = r["case"].lower().replace("-", "")
        out.setdefault(key, []).append((float(r["x_e"]), float(r["void"])))
    return out


def averaged(path: Path):
    """(x, averaged α) per axial station of a field file."""
    import h5py
    import numpy as np

    with h5py.File(path, "r") as f:
        g = f["VTKHDF/internal"] if "internal" in f["VTKHDF"] else f["VTKHDF"]
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        offs = g["Offsets"][:]
        names = list(g["CellData"])
        if "volume_fraction" not in names:
            sys.exit(f"{path}: no cell field volume_fraction; it has {', '.join(sorted(names))}")
        a = g["CellData"]["volume_fraction"][:]
    n = np.diff(offs)
    c = np.add.reduceat(pts[conn], offs[:-1], axis=0) / n[:, None]
    x = np.round(c[:, 0], 9)
    r = np.hypot(c[:, 1], c[:, 2])
    xs = np.unique(x)
    out = []
    for xv in xs:
        m = x == xv
        out.append((float(xv), float((a[m] * r[m]).sum() / r[m].sum())))
    return out


def quality(case: str, x: float) -> float:
    _, g, q, tin, ts, hfg, cp = CASES[case]
    return (cp * (tin - ts) + 4.0 * q * x / (g * D)) / hfg


def interp(xs: list[float], ys: list[float], x: float) -> float | None:
    if x < xs[0] or x > xs[-1]:
        return None
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i]) if xs[i + 1] > xs[i] else 0.0
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


def onset(xs: list[float], ys: list[float], level: float = 0.02) -> float | None:
    for i in range(len(xs) - 1):
        if ys[i] < level <= ys[i + 1]:
            return xs[i] + (level - ys[i]) / (ys[i + 1] - ys[i]) * (xs[i + 1] - xs[i])
    return None


def main() -> None:
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    out = Path(args[0])
    chart = Path(args[args.index("--chart") + 1]) if "--chart" in args else None
    fname = args[args.index("--fields") + 1] if "--fields" in args else "fields.vtkhdf"
    data = read_data()
    runs = {}
    print("case  points  mean|dev|  max|dev|  onset x_e run  onset x_e data")
    for case in CASES:
        path = out / case / fname
        if not path.exists():
            continue
        st = averaged(path)
        xe = [quality(case, x) for x, _ in st]
        al = [a for _, a in st]
        runs[case] = (xe, al)
        pts = data.get(case, [])
        dev = []
        for q, a in pts:
            m = interp(xe, al, q)
            if m is not None:
                dev.append(abs(m - a))
        dpts = sorted(pts)
        o_run = onset(xe, al)
        o_dat = onset([p[0] for p in dpts], [p[1] for p in dpts])
        fmt = lambda v: f"{v:+.3f}" if v is not None else "  -   "
        mean = sum(dev) / len(dev) if dev else float("nan")
        print(
            f"{case:5} {len(dev):6}  {mean:9.3f}  {max(dev) if dev else float('nan'):8.3f}"
            f"  {fmt(o_run):>13}  {fmt(o_dat):>14}"
        )
    if chart and runs:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
        colours = ["#d62728", "#1f77b4", "#2ca02c"]
        for ax, group in zip(axes, (("b5", "b7", "b9"), ("b11", "b12", "b13"))):
            for col, case in zip(colours, group):
                p = CASES[case]
                label = f"{case.upper()}: {p[0] / 1e6:.2f} MPa, q {p[2] / 1e6:.2f} MW/m²"
                if case in runs:
                    ax.plot(*runs[case], color=col, lw=1.8, label=label)
                pts = data.get(case, [])
                ax.plot([q for q, _ in pts], [a for _, a in pts], "o", mfc="none", color=col, ms=4)
            ax.set_xlabel("equilibrium quality $x_e$")
            ax.set_xlim(-0.2, 0.1)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8, loc="upper left")
        axes[0].set_ylabel("cross-section averaged void fraction")
        axes[0].set_ylim(-0.02, 0.65)
        fig.suptitle("Two-fluid wall boiling (lines) against Bartolomei et al. 1982 (circles)", fontsize=10)
        fig.tight_layout()
        fig.savefig(chart, dpi=150)
        print(f"wrote {chart}")


if __name__ == "__main__":
    main()
