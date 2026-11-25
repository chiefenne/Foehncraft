"""Hot air cooled by a water spray in the steady spray mode (D260, SP2s;
gallery and recipe of the SP2s pull request).

spray_cooling_steady.png: top, the steady air temperature in the duct with
the stream paths coloured by the droplet temperature; bottom left, the
mass-flow mean air temperature along the duct, steady against the transient
run of spray_cooling.toml at its last step; bottom right, over the stream
passes, the source residual and the heat the air gives the spray, against
the time mean of the transient run.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_cooling_steady.toml --out output/spray_cooling_steady
    target/release/foehn flow examples/cases/flow/spray_cooling.toml --out output/spray_cooling
    python3 -m tools.plot.spray_cooling_steady OUT_DIR output/spray_cooling_steady [output/spray_cooling]

The transient run is optional; without it the comparison curves are left out.
The time mean of the transient heat is taken after 0.3 s.
"""
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from tools.plot.colormaps import RAINBOW
from tools.plot.turbulence_fields import read


def field(d, name, path):
    if name in d:
        return d[name]
    raise SystemExit(f"{path} has no cell field `{name}`; it has {', '.join(sorted(d))}")


def streams(path):
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        pd = g["PointData"]
        for n in ("stream_id", "temperature"):
            if n not in pd:
                raise SystemExit(f"{path} has no point field `{n}`; it has {', '.join(sorted(pd))}")
        return g["Points"][:], pd["stream_id"][:], pd["temperature"][:]


def csv(path, cols):
    d = np.genfromtxt(path, delimiter=",", names=True, comments=None, skip_header=2)
    d = np.atleast_1d(d)
    for c in cols:
        if c not in d.dtype.names:
            raise SystemExit(f"{path} has no column `{c}`; it has {', '.join(d.dtype.names)}")
    return d


def mixed(path):
    quads, _, d = read(path)
    temp = field(d, "temperature", path)
    ux = field(d, "velocity_x", path)
    xs = np.round(quads.mean(axis=1)[:, 0], 9)
    cols = np.unique(xs)
    tm = np.array([np.sum(ux[xs == c] * temp[xs == c]) / np.sum(ux[xs == c]) for c in cols])
    return quads, temp, cols, tm


def plot(out, run, transient):
    quads, temp, cols, tm = mixed(os.path.join(run, "fields.vtkhdf"))
    x, sid, tp = streams(os.path.join(run, "spray", "streams.vtkhdf"))
    sp = csv(os.path.join(run, "spray.csv"), ("pass", "source_residual", "gas_heat_w"))
    fig = plt.figure(figsize=(12.0, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0])
    ax = fig.add_subplot(gs[0, :])
    lo, hi = temp.min(), temp.max()
    pc = PolyCollection(quads, array=temp - 273.15, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    pc.set_clim(lo - 273.15, hi - 273.15)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="air temperature (°C)", shrink=1.0, pad=0.01)
    # every fourth stream as a line, coloured by the droplet temperature
    segs, vals = [], []
    ids = np.unique(sid)
    for i in ids[::4]:
        m = sid == i
        p = x[m][:, :2]
        t = tp[m] - 273.15
        segs.extend(np.stack([p[:-1], p[1:]], axis=1))
        vals.extend(0.5 * (t[:-1] + t[1:]))
    lc = LineCollection(segs, array=np.array(vals), cmap="Greys", linewidths=0.6, zorder=3)
    lc.set_clim(20, hi - 273.15)
    ax.add_collection(lc)
    fig.colorbar(lc, ax=ax, label="droplet temperature (°C)", shrink=1.0, pad=0.01)
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 0.4)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(
        f"Steady spray mode: air at 300 °C cooled by water streams ({len(ids)} streams, every fourth drawn)",
        fontsize=11,
    )
    a = fig.add_subplot(gs[1, 0])
    a.plot(cols, tm - 273.15, color="tab:red", lw=2.0, label="steady mode")
    b = fig.add_subplot(gs[1, 1])
    b.semilogy(sp["pass"], sp["source_residual"], color="tab:purple", lw=1.6, marker="o", ms=2.5)
    b.set_xlabel("stream pass")
    b.set_ylabel("source residual", color="tab:purple")
    b.grid(alpha=0.3, which="both")
    b2 = b.twinx()
    b2.plot(sp["pass"], sp["gas_heat_w"] / -1000.0, color="tab:blue", lw=1.8, label="steady mode")
    b2.set_ylabel("heat from the air (kW)", color="tab:blue")
    if transient:
        _, _, tc, ttm = mixed(os.path.join(transient, "fields.vtkhdf"))
        a.plot(tc, ttm - 273.15, color="tab:gray", lw=1.4, ls="--", label="transient run, last step")
        tr = csv(os.path.join(transient, "spray.csv"), ("time", "gas_heat_w"))
        late = tr["time"] > 0.3
        q = tr["gas_heat_w"][late] / -1000.0
        b2.axhline(q.mean(), color="tab:gray", ls="--", lw=1.2, label=f"transient mean after 0.3 s ({q.mean():.2f} kW)")
        b2.legend(fontsize=8, loc="center right")
    a.axvline(0.15, color="k", ls=":", lw=1.0, label="nozzle")
    a.set_xlabel("x (m)")
    a.set_ylabel("mean air temperature (°C)")
    a.set_title("Mass-flow mean air temperature along the duct", fontsize=10)
    a.legend(fontsize=9)
    a.grid(alpha=0.3)
    b.set_title("Convergence of the stream passes", fontsize=10)
    fig.savefig(os.path.join(out, "spray_cooling_steady.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    transient = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(out, exist_ok=True)
    plot(out, run, transient)


if __name__ == "__main__":
    main()
