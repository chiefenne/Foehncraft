"""A diesel-like spray breaking up in dense nitrogen (Lagrangian sprays,
breakup SP6, D300 to D308; gallery and recipe of the SP6 pull request).

spray_breakup.png: top, a side view of the spray at the end of the run, the
gas speed on the symmetry plane z = 0 (mirrored about the spray axis to
show the whole jet) with the parcels as dots coloured by the diameter of
their drops; bottom left, the liquid penetration (97 % of the liquid, and
the farthest parcel) over time against the correlation of Hiroyasu and
Arai (1990); bottom right, the Sauter mean diameter along the spray at the
end of the run. A second run (the same case with `kh-rt`) is drawn dashed
for comparison when given.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_breakup.toml --out output/spray_breakup
    target/release/foehn flow examples/cases/flow/spray_breakup.toml --out output/spray_breakup_khrt \\
        --set spray.breakup.model=kh-rt
    python3 -m tools.plot.spray_breakup OUT_DIR output/spray_breakup [output/spray_breakup_khrt]
"""
import glob
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.vtkhdf import internal  # noqa: E402

# the case (examples/cases/flow/spray_breakup.toml and its injector file)
RHO_L = 830.0       # liquid density (kg/m³)
RHO_G = 15.0        # nitrogen density (kg/m³)
DP = 60e6           # injection pressure above the gas (Pa)
D0 = 1.5e-4         # hole diameter (m)


def hiroyasu_arai(t):
    """Spray tip penetration of Hiroyasu and Arai (1990), SAE 900475 (m)."""
    tb = 28.65 * RHO_L * D0 / np.sqrt(RHO_G * DP)
    early = 0.39 * np.sqrt(2.0 * DP / RHO_L) * t
    late = 2.95 * (DP / RHO_G) ** 0.25 * np.sqrt(D0 * t)
    return np.where(t < tb, early, late)


def plane(path):
    """The cells on the symmetry plane z = 0: their x-y quads and data."""
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        if not np.all(np.diff(off) == 8):
            raise SystemExit(f"{path}: hexahedra only")
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    xyz = pts[conn.reshape(-1, 8)]
    keep = np.abs(xyz[:, :, 2].min(axis=1)) < 1e-12
    quads = []
    for c in xyz[keep]:
        q = c[np.abs(c[:, 2]) < 1e-12][:, :2]
        m = q.mean(axis=0)
        quads.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
    return np.array(quads), {k: v[keep] for k, v in data.items()}


def field(d, name):
    if name in d:
        return d[name]
    raise SystemExit(f"fields.vtkhdf has no cell field `{name}`; it has {', '.join(sorted(d))}")


def parcels(run):
    files = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))
    if not files:
        raise SystemExit(f"{run}/spray has no parcels_*.vtkhdf")
    with h5py.File(files[-1], "r") as f:
        g = f["VTKHDF"]
        pd = g["PointData"]
        for name in ("diameter", "droplets"):
            if name not in pd:
                raise SystemExit(f"{files[-1]} has no point field `{name}`; it has {', '.join(sorted(pd))}")
        return g["Points"][:], pd["diameter"][:], pd["droplets"][:], float(f.attrs["time"])


def csv(run):
    c = np.genfromtxt(os.path.join(run, "spray.csv"), delimiter=",", names=True, comments=None, skip_header=2)
    for name in ("time", "holepenetration_m", "holefarthest_m"):
        if name not in c.dtype.names:
            raise SystemExit(f"spray.csv has no column `{name}`; it has {', '.join(c.dtype.names)}")
    return c


def smd_along(x, d, n, edges):
    """Sauter mean diameter of the drops in slabs of x (µm)."""
    out = np.full(len(edges) - 1, np.nan)
    for j in range(len(edges) - 1):
        s = (x >= edges[j]) & (x < edges[j + 1])
        if s.sum() >= 5:
            out[j] = (n[s] * d[s] ** 3).sum() / (n[s] * d[s] ** 2).sum() * 1e6
    return out


def plot(out, run, other):
    quads, d = plane(os.path.join(run, "fields.vtkhdf"))
    speed = np.sqrt(field(d, "velocity_x") ** 2 + field(d, "velocity_y") ** 2 + field(d, "velocity_z") ** 2)
    x, dp, n, t_end = parcels(run)
    c = csv(run)
    fig = plt.figure(figsize=(12.0, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.0])
    ax = fig.add_subplot(gs[0, :])
    mirrored = quads * np.array([1.0, -1.0])
    both = np.concatenate([quads, mirrored]) * 1e3
    vals = np.concatenate([speed, speed])
    pc = PolyCollection(both, array=vals, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
    pc.set_clim(0.0, np.percentile(speed, 99.8))
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label="gas speed (m/s)", shrink=1.0, pad=0.01)
    # the parcels seen from the side, both halves of the jet
    order = np.argsort(dp)[::-1]
    xs, ys, ds = x[order, 0] * 1e3, x[order, 1] * 1e3, dp[order] * 1e6
    sign = np.where(np.arange(len(xs)) % 2 == 0, 1.0, -1.0)
    sc = ax.scatter(
        xs, sign * ys, c=ds, s=0.5 + ds / 8.0, cmap="Greys", norm=LogNorm(vmin=1.0, vmax=150.0),
        edgecolors="none", alpha=0.9, zorder=3,
    )
    fig.colorbar(sc, ax=ax, label="drop diameter (µm)", shrink=1.0, pad=0.01)
    reach = max(float(np.nanmax(c["holefarthest_m"])) * 1e3 * 1.15, 20.0)
    ax.set_xlim(0.0, min(reach, 100.0))
    ax.set_ylim(-12.0, 12.0)
    ax.set_aspect("equal")
    ax.set_xlabel("distance from the nozzle (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(
        f"KH-ACT breakup of a 60 MPa fuel spray in nitrogen at 15 kg/m³, t = {t_end * 1e3:.2f} ms "
        f"({len(x)} parcels)",
        fontsize=11,
    )
    a = fig.add_subplot(gs[1, 0])
    tt = np.linspace(0.0, float(c["time"][-1]), 400)
    a.plot(tt * 1e3, hiroyasu_arai(tt) * 1e3, color="k", ls="--", lw=1.2, label="Hiroyasu and Arai (1990)")
    a.plot(c["time"] * 1e3, c["holepenetration_m"] * 1e3, color="tab:red", lw=2.0, label="KH-ACT, 97 % of the liquid")
    a.plot(c["time"] * 1e3, c["holefarthest_m"] * 1e3, color="tab:red", lw=1.0, alpha=0.6, label="KH-ACT, farthest parcel")
    if other:
        o = csv(other)
        a.plot(o["time"] * 1e3, o["holepenetration_m"] * 1e3, color="tab:blue", lw=1.6, ls="--",
               label="KH-RT, 97 % of the liquid")
    a.set_xlabel("time after the start of injection (ms)")
    a.set_ylabel("liquid penetration (mm)")
    a.grid(alpha=0.3)
    a.legend(fontsize=8.5, loc="upper left")
    a.set_title("Penetration of the liquid", fontsize=10)
    b = fig.add_subplot(gs[1, 1])
    edges = np.arange(0.0, min(reach, 100.0) + 2.0, 2.0) * 1e-3
    mid = 0.5 * (edges[1:] + edges[:-1]) * 1e3
    b.plot(mid, smd_along(x[:, 0], dp, n, edges), color="tab:red", lw=2.0, marker="o", ms=3, label="KH-ACT")
    if other:
        xo, do, no, _ = parcels(other)
        b.plot(mid, smd_along(xo[:, 0], do, no, edges), color="tab:blue", lw=1.6, ls="--", marker="s", ms=3,
               label="KH-RT")
    b.set_ylim(bottom=0.0)
    b.set_xlabel("distance from the nozzle (mm)")
    b.set_ylabel("Sauter mean diameter (µm)")
    b.grid(alpha=0.3)
    b.legend(fontsize=9)
    b.set_title(f"Drop size along the spray at t = {t_end * 1e3:.2f} ms (2 mm slabs)", fontsize=10)
    path = os.path.join(out, "spray_breakup.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    other = sys.argv[3] if len(sys.argv) == 4 else None
    os.makedirs(out, exist_ok=True)
    plot(out, run, other)


if __name__ == "__main__":
    main()
