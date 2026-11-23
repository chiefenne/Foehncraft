"""Exact reductions over ranks (MPI step W2), gallery 219: for 2 to 16
ranks, how far the old sum (pairwise per rank, then in rank order) lands
from the exact one-rank value, in units in the last place, for the domain
volume, its first moment and the boundary area of
`mesh/exact_reductions.toml`; the exact sums over the ranks are at 0 on
every rank count. From the `mesh_report.json` of one run per rank count:

    F=examples/cases/mesh/exact_reductions.toml
    for n in $(seq 2 16); do
        target/release/foehn mesh $F --parts $n --out output/exact_reductions/p$n
    done
    pip install matplotlib
    python3 -m tools.plot.exact_reductions docs/figures/gallery/exact_reductions_ulps.png \\
        output/exact_reductions
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    out, root = Path(sys.argv[1]), Path(sys.argv[2])
    rows = {}
    for rep in sorted(root.glob("p*/mesh_report.json")):
        r = json.loads(rep.read_text())
        parts = r["partition"]["parts"]
        rows[parts] = {x["name"]: x for x in r["reductions"]}
    if not rows:
        print(f"no p*/mesh_report.json under {root}", file=sys.stderr)
        return 1
    parts = sorted(rows)
    names = list(rows[parts[0]])
    fig, ax = plt.subplots(figsize=(8.0, 4.48), dpi=125)
    colours = ["#1f77b4", "#d62728", "#2ca02c"]
    for name, col in zip(names, colours):
        ulps = [rows[p][name]["pairwise_ulps"] for p in parts]
        ax.plot(parts, ulps, "o-", color=col, lw=1.6, ms=5,
                label=f"old: {name}")
    exact_bad = [p for p in parts for n in names if not rows[p][n]["bitwise_equal"]]
    ax.plot(parts, [0] * len(parts), "s-", color="black", lw=2.2,
            ms=6, label="exact sums (W2), all three", zorder=5)
    ax.set_xlabel("ranks")
    ax.set_ylabel("distance from the one-rank value [ulp]")
    ax.set_xticks(parts)
    ax.set_ylim(-0.4, 3.2)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=9, loc="upper center", ncol=2)
    ax.set_title("Sums over ranks: the old rank-order sum moves, the exact sum does not"
                 + ("" if not exact_bad else "  (EXACT DIFFERS!)"), fontsize=11)
    fig.tight_layout()
    fig.savefig(out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
