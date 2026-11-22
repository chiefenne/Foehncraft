"""Figure of the transient diffuser runs (gallery 199): the velocity u at
four probes over time for two transient runs on level 1 that start from
steady restart files, one from the converged steady state and one from a
state where SIMPLEC at relaxation 0.5 had stalled after 1 700 iterations.

    python3 -m tools.plot.diffuser_transient OUT.png CONVERGED_RUN STALLED_RUN

Each RUN directory holds the probes.csv of a transient run of
examples/cases/flow/diffuser_sst.toml with the probes `sep` (10 H, -1.5 H),
`shear` (15 H, -1 H), `reat` (28 H, -3.5 H) and `core` (15 H, 0.5 H) and
`--restart` from the steady restart file (docs/recipes/turbulence.md,
recipe 199).
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PROBES = [
    ("sep", "separated region (10 H, −1.5 H)", "#d1495b"),
    ("shear", "shear layer (15 H, −1 H)", "#edae49"),
    ("reat", "near reattachment (28 H, −3.5 H)", "#00798c"),
    ("core", "jet core (15 H, 0.5 H)", "#30638e"),
]


def read(path: Path):
    rows = [x for x in path.read_text().splitlines() if x and not x.startswith("#")]
    head = rows[0].split(",")
    data = np.array([[float(v) for v in r.split(",")] for r in rows[1:]])
    return {k: data[:, i] for i, k in enumerate(head)}


def main(out: str, converged: str, stalled: str) -> None:
    a, b = read(Path(converged) / "probes.csv"), read(Path(stalled) / "probes.csv")
    fig, axs = plt.subplots(2, 2, figsize=(10, 6.4), sharex=True)
    for ax, (name, label, colour) in zip(axs.flat, PROBES):
        k = f"{name}_u"
        ax.plot(b["time"], b[k], color=colour, lw=2.2, label="from the stalled state (relaxation 0.5)")
        ax.plot(a["time"], a[k], color="k", lw=1.4, ls="--", label="from the converged state (relaxation 0.8)")
        ax.set_title(label, fontsize=10)
        ax.set_ylabel("u / U_c")
        ax.grid(alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel("t U_c / H")
    axs[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle("Diffuser, level 1 (27 840 cells): transient runs started from steady restart files", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=200)


if __name__ == "__main__":
    main(*sys.argv[1:4])
