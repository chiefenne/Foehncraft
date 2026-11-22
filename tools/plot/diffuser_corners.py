"""Figure of the diffuser with rounded and with sharp corners (gallery 199):
skin friction C_f = tau_x / (1/2 rho U_c^2) along the lower wall for
three runs on level 1, and the start of the diffuser enlarged.

    python3 -m tools.plot.diffuser_corners OUT.png ROUNDED SHARP SHARP_PLAIN

ROUNDED is a run of examples/cases/flow/diffuser_sst.toml, SHARP one of
diffuser_sst_sharp.toml (cells shortened towards the corners) and
SHARP_PLAIN one of diffuser_sst.toml on a sharp-cornered mesh with the
plain 0.25 H cells (docs/recipes/turbulence.md, recipe 199). Each holds the
wall_lower.csv of the run (matplotlib, numpy).
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

RUNS = [
    ("rounded corners (R = 4.3 H)", "#30638e", "-"),
    ("sharp corners, short cells at the corners", "#d1495b", "-"),
    ("sharp corners, plain 0.25 H cells", "#edae49", "--"),
]


def wall(d: Path):
    rows = [x for x in (d / "wall_lower.csv").read_text().splitlines() if x and not x.startswith("#")]
    head = rows[0].split(",")
    a = np.array([[float(v) for v in r.split(",")] for r in rows[1:]])
    a = a[np.argsort(a[:, head.index("x")])]
    return a[:, head.index("x")], 2.0 * a[:, head.index("tau_x")]


def main(out: str, runs: list[str]) -> None:
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), dpi=180, gridspec_kw={"width_ratios": [2.2, 1]})
    for (label, colour, style), d in zip(RUNS, runs):
        x, cf = wall(Path(d))
        for ax in axs:
            ax.plot(x, 1e3 * cf, style, color=colour, lw=1.6, label=label)
    for ax, (x0, x1), ttl in zip(axs, [(-5, 40), (-1, 4)], ["lower wall", "start of the diffuser (enlarged)"]):
        ax.axhline(0.0, color="k", lw=0.6)
        ax.set_xlim(x0, x1)
        ax.set_xlabel("x / H")
        ax.set_title(ttl, fontsize=10)
        ax.grid(alpha=0.3)
    axs[0].axvspan(0.0, 20.98, color="#dddddd", alpha=0.4, lw=0)
    axs[0].text(10.5, 3.2, "diffuser", ha="center", fontsize=9)
    axs[0].set_ylabel("C_f × 1000")
    axs[0].set_ylim(-1.0, 3.8)
    axs[1].set_ylim(-0.6, 13.0)
    axs[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Diffuser, level 1: skin friction of the lower wall with rounded and with sharp corners", fontsize=11)
    fig.tight_layout()
    fig.savefig(out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:5])
