"""Figure of the turbulent pipe (gallery 193): the axial velocity and the
eddy viscosity on a cross section (the wedge profile revolved) and the mean
velocity in wall units against the log law of El Khoury et al. (2013)
(kappa = 0.41, B = 5.2) and the DNS centreline value.

    python3 -m tools.plot.pipe_chart OUT.png RUN_DIR

RUN_DIR holds the line_profile.csv and flow.json of the run of
examples/cases/flow/pipe_turbulent_sst.toml.
"""

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from tools.plot.colormaps import RAINBOW  # noqa: E402

R, NU = 0.5, 1.0 / 5300.0


def read(path: Path):
    rows = [l for l in path.read_text().splitlines() if l and not l.startswith("#")]
    head = rows[0].split(",")
    data = np.array([[float(x) for x in r.split(",")] for r in rows[1:]])
    return {k: data[:, i] for i, k in enumerate(head)}


def main(out: str, run: str) -> None:
    d = Path(run)
    p = read(d / "line_profile.csv")
    f = json.loads((d / "flow.json").read_text())["solves"][-1]["periodic"][0]["driving_force"]
    ut = math.sqrt(f * R / 2.0)
    r, u, nut = p["y"], p["u"], p["nu_t"]
    th = np.linspace(0, 2 * np.pi, 361)
    rr, tt = np.meshgrid(r, th)
    fig = plt.figure(figsize=(11, 3.8), dpi=200)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.35], wspace=0.7)
    for i, (val, lab, cmap) in enumerate([(u, "u / U_b", RAINBOW), (nut / NU, "$\\nu_t$ / $\\nu$", RAINBOW)]):
        ax = fig.add_subplot(gs[i])
        v = np.tile(val, (len(th), 1))
        m = ax.pcolormesh(rr * np.cos(tt), rr * np.sin(tt), v[:-1, :-1], cmap=cmap, shading="flat")
        ax.set_aspect("equal")
        ax.axis("off")
        fig.colorbar(m, ax=ax, shrink=0.75, label=lab)
    ax = fig.add_subplot(gs[2])
    yp = (R - r) * ut / NU
    ok = yp > 0.1
    ax.semilogx(yp[ok], u[ok] / ut, lw=2, label="k-ω SST, 384 cells")
    ypl = np.logspace(0, math.log10(R * ut / NU), 100)
    ax.semilogx(ypl, np.where(ypl < 5, ypl, np.nan), "k:", lw=1, label="u⁺ = y⁺")
    ax.semilogx(ypl[ypl > 30], (np.log(ypl) / 0.41 + 5.2)[ypl > 30], "k--", lw=1, label="log law (κ 0.41, B 5.2)")
    ax.plot([R * ut / NU], [19.13], "o", color="crimson", label="DNS U_cl⁺ 19.13")
    ax.set_xlabel("y⁺")
    ax.set_ylabel("u⁺")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3)
    fig.savefig(out, bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
