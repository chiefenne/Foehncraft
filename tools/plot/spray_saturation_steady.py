"""More water than the air can take, in the steady spray mode (V28.11;
SP5 with D260; gallery and recipe of the SP5 pull request).

spray_saturation_steady.png: top, the steady water vapour mass fraction in
the duct with the stream paths coloured by the droplet diameter; bottom
left, the mass-flow mean air temperature and vapour mass fraction along the
duct against the exact adiabatic saturation of the inlet air; bottom right,
over the stream passes, the source residual and the vapour the air receives
against the exact one.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_saturation_steady.toml --out output/spray_saturation_steady
    python3 -m tools.plot.spray_saturation_steady OUT_DIR output/spray_saturation_steady

The exact state is computed here from the inlet (400 K, Y_H2O = 0.001,
17.6 g/s of air, 1010 J/(kg K)) with the saturation pressure of the code
(IAPWS-95 auxiliary equation) and the latent heat of water with a constant
vapour specific heat: c_p dT = -L(T) dY / (1 - Y) up to Y = Y_s(T).
"""
import math
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

A = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
M_GAS = 1.0 / (0.2328 / 0.031998 + 0.7662 / 0.028014)
M_VAP = 0.018015268


def y_sat(t):
    th = 1.0 - t / 647.096
    s = math.sqrt(th)
    p = 22.064e6 * math.exp(
        647.096 / t * (A[0] * th + A[1] * th * s + A[2] * th**3 + A[3] * th**3 * s + A[4] * th**4 + A[5] * th**7.5)
    )
    x = p / 101325.0
    return x * M_VAP / (x * M_VAP + (1.0 - x) * M_GAS)


def exact(t=400.0, y=0.001, cp=1010.0, air=0.0176):
    latent = lambda tt: 2.4417e6 + (1864.0 - 4181.0) * (tt - 298.15)
    dy = 1e-6
    while y < y_sat(t):
        t -= latent(t) * dy / (1.0 - y) / cp
        y += dy
    dry = air * (1.0 - 0.001)
    return t, y, dry / (1.0 - y) - dry / (1.0 - 0.001)


def field(d, name, path):
    if name in d:
        return d[name]
    raise SystemExit(f"{path} has no cell field `{name}`; it has {', '.join(sorted(d))}")


def streams(path):
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        pd = g["PointData"]
        for n in ("stream_id", "diameter"):
            if n not in pd:
                raise SystemExit(f"{path} has no point field `{n}`; it has {', '.join(sorted(pd))}")
        return g["Points"][:], pd["stream_id"][:], pd["diameter"][:]


def csv(path, cols):
    d = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True, comments=None, skip_header=2))
    for c in cols:
        if c not in d.dtype.names:
            raise SystemExit(f"{path} has no column `{c}`; it has {', '.join(d.dtype.names)}")
    return d


def plot(out, run):
    path = os.path.join(run, "fields.vtkhdf")
    quads, _, d = read(path)
    y = field(d, "Y_H2O", path)
    temp = field(d, "temperature", path)
    ux = field(d, "velocity_x", path)
    xs = np.round(quads.mean(axis=1)[:, 0], 9)
    cols = np.unique(xs)
    tm = np.array([np.sum(ux[xs == c] * temp[xs == c]) / np.sum(ux[xs == c]) for c in cols])
    ym = np.array([np.sum(ux[xs == c] * y[xs == c]) / np.sum(ux[xs == c]) for c in cols])
    x, sid, dp = streams(os.path.join(run, "spray", "streams.vtkhdf"))
    sp = csv(os.path.join(run, "spray.csv"), ("pass", "source_residual", "gas_vapour_kg_s"))
    t_ex, y_ex, m_ex = exact()
    fig = plt.figure(figsize=(12.0, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.0])
    ax = fig.add_subplot(gs[0, :])
    pc = PolyCollection(quads, array=y * 100.0, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="water vapour mass fraction (%)", shrink=1.0, pad=0.01)
    segs, vals = [], []
    ids = np.unique(sid)
    for i in ids[::5]:
        m = sid == i
        p = x[m][:, :2]
        dd = dp[m] * 1e6
        segs.extend(np.stack([p[:-1], p[1:]], axis=1))
        vals.extend(0.5 * (dd[:-1] + dd[1:]))
    lc = LineCollection(segs, array=np.array(vals), cmap="Greys", linewidths=0.6, zorder=3)
    lc.set_clim(0.0, 25.0)
    ax.add_collection(lc)
    fig.colorbar(lc, ax=ax, label="droplet diameter (µm)", shrink=1.0, pad=0.01)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(
        f"Steady spray mode: 2 g/s of water in 17.6 g/s of air at 400 K ({len(ids)} streams, every fifth drawn)",
        fontsize=11,
    )
    a = fig.add_subplot(gs[1, 0])
    a.plot(cols, tm, color="tab:red", lw=2.0, label="air temperature")
    a.axhline(t_ex, color="tab:red", ls="--", lw=1.0, label=f"adiabatic saturation {t_ex:.1f} K")
    a.set_ylim(300.0, 405.0)
    a.set_xlabel("x (m)")
    a.set_ylabel("mean air temperature (K)", color="tab:red")
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="center right")
    a2 = a.twinx()
    a2.plot(cols, ym * 100.0, color="tab:blue", lw=2.0)
    a2.axhline(y_ex * 100.0, color="tab:blue", ls="--", lw=1.0)
    a2.set_ylim(0.0, 4.5)
    a2.set_ylabel(f"mean vapour mass fraction (%), exact {y_ex * 100:.2f}", color="tab:blue")
    a.set_title("Mass-flow means along the duct against the exact end state", fontsize=10)
    b = fig.add_subplot(gs[1, 1])
    b.semilogy(sp["pass"], sp["source_residual"], color="tab:purple", lw=1.6, marker="o", ms=2.5)
    b.set_xlabel("stream pass")
    b.set_ylabel("source residual", color="tab:purple")
    b.grid(alpha=0.3, which="both")
    b2 = b.twinx()
    b2.plot(sp["pass"], sp["gas_vapour_kg_s"] * 1000.0, color="tab:blue", lw=1.8)
    b2.axhline(m_ex * 1000.0, color="tab:blue", ls="--", lw=1.0)
    b2.set_ylabel(f"vapour to the air (g/s), exact {m_ex * 1000:.3f}", color="tab:blue")
    b.set_title("The stream passes: the residual falls to 0.1 and stalls", fontsize=10)
    fig.savefig(os.path.join(out, "spray_saturation_steady.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(out, run)


if __name__ == "__main__":
    main()
