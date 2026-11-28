"""Bubbly up-flow in a vertical pipe (stage MP3, D356): the radial void
profiles of the two-fluid model at L/D = 35 against the measurements of
Wang et al. (1987), j_L = 1.08 m/s, j_G = 0.1 and 0.4 m/s.

    foehn flow examples/cases/flow/bubbly_pipe_wang_jg01.toml --out OUT/jg01
    foehn flow examples/cases/flow/bubbly_pipe_wang_jg04.toml --out OUT/jg04
    python3 examples/studies/bubbly_pipe/compare.py OUT [--chart FILE]

Reads OUT/jg01/line_ld35.csv and OUT/jg04/line_ld35.csv (the line from the
axis to the wall at x = 35 D) and the measured points of
examples/studies/bubbly_pipe/data/wang1987_void.csv. Prints per case the
area-averaged void fraction of the run, the core value (r/R < 0.6) and
the wall peak (position and value) of the run and of the measurement.
`--chart FILE` writes a figure (SVG or PNG by the suffix) of both
profiles with the measured points.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "wang1987_void.csv"
CASES = {"jg01": 0.1, "jg04": 0.4}
RADIUS = 0.028575


def read_line(path: Path) -> dict[str, list[float]]:
    """The columns of a `foehncraft` line file."""
    rows = [r for r in path.read_text().splitlines() if r and not r.startswith("#")]
    head = rows[0].split(",")
    cols: dict[str, list[float]] = {h: [] for h in head}
    for r in rows[1:]:
        for h, v in zip(head, r.split(",")):
            cols[h].append(float(v))
    if "alpha_gas" not in cols:
        avail = ", ".join(head)
        sys.exit(f"{path}: no column alpha_gas (a two-fluid run?); available: {avail}")
    return cols


def read_data() -> dict[float, list[tuple[float, float]]]:
    out: dict[float, list[tuple[float, float]]] = {}
    with DATA.open() as f:
        rows = [r for r in f if not r.startswith("#")]
    for r in csv.DictReader(rows):
        out.setdefault(float(r["jg"]), []).append((float(r["r_over_R"]), float(r["alpha"])))
    return out


def area_mean(r: list[float], a: list[float]) -> float:
    """∫ α 2 r dr / R² by the trapezoidal rule."""
    num = sum(0.5 * (a[i] * r[i] + a[i + 1] * r[i + 1]) * (r[i + 1] - r[i]) for i in range(len(r) - 1))
    return 2.0 * num / r[-1] ** 2


def summary(points: list[tuple[float, float]]) -> tuple[float, float, float]:
    core = [a for x, a in points if x < 0.6]
    x, a = max(points, key=lambda p: p[1])
    return sum(core) / len(core), x, a


def main() -> None:
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    out = Path(args[0])
    chart = Path(args[args.index("--chart") + 1]) if "--chart" in args else None
    data = read_data()
    runs = {}
    for key, jg in CASES.items():
        path = out / key / "line_ld35.csv"
        if not path.exists():
            print(f"{key}: {path} missing, skipped")
            continue
        c = read_line(path)
        rr = [y / RADIUS for y in c["y"]]
        runs[key] = (rr, c["alpha_gas"])
        sim = list(zip(rr, c["alpha_gas"]))
        sc, sx, sa = summary(sim)
        mc, mx, ma = summary(data[jg])
        print(
            f"j_G = {jg} m/s: area mean {area_mean(rr, c['alpha_gas']):.4f}; "
            f"core {sc:.4f} (measured {mc:.4f}); "
            f"wall peak {sa:.3f} at r/R {sx:.3f} (measured {ma:.3f} at {mx:.3f})"
        )
    if chart and runs:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        colours = {"jg01": "#1f6fb4", "jg04": "#d1495b"}
        for key, (rr, a) in runs.items():
            jg = CASES[key]
            ax.plot(rr, a, color=colours[key], lw=2, label=f"Foehncraft, $j_G$ = {jg} m/s")
            x, y = zip(*data[jg])
            ax.plot(x, y, "o", color=colours[key], mfc="white", label=f"Wang et al. (1987), $j_G$ = {jg} m/s")
        ax.set_xlabel("r/R")
        ax.set_ylabel("void fraction α")
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, frameon=False)
        ax.set_title("Bubbly up-flow, D = 57 mm, $j_L$ = 1.08 m/s, L/D = 35", fontsize=10)
        fig.tight_layout()
        fig.savefig(chart, dpi=150)
        print(f"wrote {chart}")


if __name__ == "__main__":
    main()
