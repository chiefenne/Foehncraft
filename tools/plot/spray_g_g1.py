"""ECN Spray G, condition G1, on a 22.5° sector (D309; gallery and recipe of
its pull request): the liquid and vapour penetration along the injector
axis against the measurements of the Engine Combustion Network.

spray_g_g1.png: left, a side view at `SNAPSHOT` ms with the injector at the
left and its axis horizontal: above the axis the fuel vapour mass fraction
on the plane of plumes 1 and 5 (the sector's plane `zmin`, plume 5 by the
symmetry of the injector), below it the projected liquid volume of all 8
plumes seen perpendicular to that plane, as the extinction imaging of the
ECN sees the spray, with its contour at 0.2e-3 mm³/mm²; right, the liquid
penetration (projected liquid volume above 0.2e-3 and 2e-3 mm³/mm²) and the
vapour penetration (vapour mass fraction above 0.1 %) over time against the
measured ones.

The measured penetrations are quoted from the Sandia Spray G data of the
ECN (https://ecn.sandia.gov/data/sandia-spray-g-data/, iso-octane, G1:
files liquid_data_2e-4_thr.mat, liquid_data_2e-3_thr.mat, column 16, and
vapor_data.mat, column 12; high-speed extinction and schlieren imaging by
J. Hwang and L. M. Pickett, cited as Hwang et al., Fuel 280, 2020, 118359),
at the frames nearest to the listed times; their uncertainties are below
0.3 mm before 1 ms.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1
    python3 -m tools.plot.spray_g_g1 OUT_DIR output/spray_g_g1
"""
import glob
import os
import re
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.spray_g_metrics import (  # noqa: E402
    PIXEL,
    SECTOR,
    VAPOUR_THRESHOLD,
    liquid_penetration,
    plv_image,
    vapour_penetration,
)
from tools.plot.vtkhdf import internal  # noqa: E402

VAPOUR = "IC8H18"
SNAPSHOT = 0.6  # ms, the side view

# measured (ECN, see the module text): (time after the start of injection
# in ms, penetration along the injector axis in mm)
LIQUID_LOW = [  # projected liquid volume above 0.2e-3 mm³/mm²
    (0.045, 3.21), (0.104, 8.63), (0.149, 12.00), (0.298, 20.44), (0.402, 23.78),
    (0.506, 26.92), (0.595, 29.32), (0.699, 31.68), (0.804, 33.82), (0.997, 36.68),
    (1.101, 37.92), (1.205, 38.89),
]
LIQUID_HIGH = [  # above 2e-3 mm³/mm²; the dense liquid is gone after 1 ms
    (0.045, 3.02), (0.104, 8.45), (0.149, 11.81), (0.298, 19.06), (0.402, 22.11),
    (0.506, 24.59), (0.595, 25.76), (0.699, 26.38), (0.804, 26.96), (0.997, 19.45),
]
VAPOUR_PEN = [  # schlieren boundary
    (0.045, 3.87), (0.104, 8.80), (0.149, 12.21), (0.298, 20.49), (0.402, 23.88),
    (0.506, 26.95), (0.595, 29.41), (0.699, 32.02), (0.804, 34.38), (0.997, 38.31),
    (1.101, 40.22), (1.205, 42.01),
]


def step_of(path):
    return int(re.search(r"_(\d+)\.vtkhdf$", path).group(1))


def times(run):
    """Step -> time (s) from history.csv."""
    out = {}
    with open(os.path.join(run, "history.csv")) as f:
        for line in f:
            cols = line.split(",")
            if cols[0].strip().isdigit():
                out[int(cols[0])] = float(cols[1])
    return out


def parcel_files(run):
    files = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")), key=step_of)
    if not files:
        raise SystemExit(f"{run}/spray has no parcels_*.vtkhdf")
    return files


def parcels(path):
    with h5py.File(path, "r") as f:
        pd = f["VTKHDF"]["PointData"]
        for name in ("diameter", "droplets"):
            if name not in pd:
                raise SystemExit(f"{path} has no point field `{name}`; it has {', '.join(sorted(pd))}")
        return f["VTKHDF"]["Points"][:], pd["diameter"][:], pd["droplets"][:], float(f.attrs["time"])


def cells(path):
    """Cell centres, nodes of the hexahedra and cell data of a field file."""
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        if not np.all(np.diff(off) == 8):
            raise SystemExit(f"{path}: hexahedra only")
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    xyz = pts[conn.reshape(-1, 8)]
    return xyz.mean(axis=1), xyz, data


def vapour_field(d, path):
    for name in (VAPOUR, f"Y_{VAPOUR}", f"mass_fraction_{VAPOUR}"):
        if name in d:
            return d[name]
    raise SystemExit(f"{path} has no cell field of `{VAPOUR}`; it has {', '.join(sorted(d))}")


def plume_plane(xyz):
    """The cells with a face on the plume plane φ = −SECTOR/2 (`zmin`): their
    (x, r) quads in mm and a mask of those cells."""
    phi0 = -0.5 * np.radians(SECTOR)
    phi = np.arctan2(xyz[:, :, 2], xyz[:, :, 1])
    on = np.abs(phi - phi0) < 1e-9
    keep = on.sum(axis=1) == 4
    quads = []
    for c, o in zip(xyz[keep], on[keep]):
        q = np.column_stack([c[o, 0], np.hypot(c[o, 1], c[o, 2])]) * 1e3
        m = q.mean(axis=0)
        quads.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
    return np.array(quads), keep


def penetrations(run):
    """(t ms, liquid 0.2e-3, liquid 2e-3) from the parcel files and (t ms,
    vapour) from the field files, in mm."""
    liq = []
    for p in parcel_files(run):
        x, d, n, t = parcels(p)
        liq.append((t * 1e3, liquid_penetration(x, d, n, 0.2e-3), liquid_penetration(x, d, n, 2e-3)))
    tmap = times(run)
    vap = [(0.0, 0.0)]
    for p in sorted(glob.glob(os.path.join(run, "fields_*.vtkhdf")), key=step_of):
        c, _, d = cells(p)
        vap.append((tmap[step_of(p)] * 1e3, vapour_penetration(c, vapour_field(d, p), VAPOUR_THRESHOLD)))
    return np.array([(0.0, 0.0, 0.0)] + liq), np.array(vap)


def nearest(files, t_ms, tmap=None):
    best, gap = None, np.inf
    for p in files:
        t = tmap[step_of(p)] * 1e3 if tmap is not None else parcels(p)[3] * 1e3
        if abs(t - t_ms) < gap:
            best, gap = p, abs(t - t_ms)
    return best


def plot(out, run):
    liq, vap = penetrations(run)
    tmap = times(run)
    ffile = nearest(sorted(glob.glob(os.path.join(run, "fields_*.vtkhdf")), key=step_of), SNAPSHOT, tmap)
    pfile = nearest(parcel_files(run), SNAPSHOT)
    _, xyz, d = cells(ffile)
    quads, keep = plume_plane(xyz)
    y = vapour_field(d, ffile)[keep]
    x, dp, n, t_p = parcels(pfile)
    img = plv_image(x, dp, n)

    fig = plt.figure(figsize=(12.0, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0])
    ax = fig.add_subplot(gs[:, 0])
    pc = PolyCollection(quads, array=np.maximum(y, 1e-5), cmap=RAINBOW, norm=LogNorm(1e-3, 0.3),
                        edgecolors="face", linewidths=0.2)
    ax.add_collection(pc)
    nx, ny = img.shape
    ext = [0.0, nx * PIXEL, -0.5 * ny * PIXEL, 0.5 * ny * PIXEL]
    lower = np.where(np.arange(ny) * PIXEL - 0.5 * ny * PIXEL < 0.0, 1.0, np.nan)
    shown = np.ma.masked_invalid(np.where(img > 1e-5, img, np.nan) * lower[None, :])
    im = ax.imshow(shown.T, origin="lower", extent=ext, cmap="bone_r", norm=LogNorm(1e-4, 5e-2),
                   interpolation="nearest", aspect="auto")
    xs = (np.arange(nx) + 0.5) * PIXEL
    ys = (np.arange(ny) + 0.5) * PIXEL - 0.5 * ny * PIXEL
    ax.contour(xs, ys, np.where(ys[None, :] < 0.0, img, 0.0).T, levels=[0.2e-3], colors="tab:red",
               linewidths=1.0)
    fig.colorbar(im, ax=ax, label="projected liquid volume of all 8 plumes, mm³/mm² (below)", shrink=0.85,
                 pad=0.01, location="bottom", aspect=40)
    fig.colorbar(pc, ax=ax, label=f"{VAPOUR} vapour mass fraction, plane of plumes 1 and 5 (above)",
                 shrink=0.85, pad=0.01, location="bottom", aspect=40)
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_xlim(0.0, 50.0)
    ax.set_ylim(-40.0, 40.0)
    ax.set_aspect("equal")
    ax.set_xlabel("distance from the injector tip along its axis (mm)")
    ax.set_ylabel("distance from the injector axis (mm)")
    ax.set_title(f"ECN Spray G (G1), t = {t_p * 1e3:.2f} ms: vapour above, liquid below", fontsize=10.5)

    a = fig.add_subplot(gs[0, 1])
    m = np.array(LIQUID_LOW)
    a.plot(m[:, 0], m[:, 1], "o", color="k", ms=4, label="measured, 0.2e-3 mm³/mm²")
    m = np.array(LIQUID_HIGH)
    a.plot(m[:, 0], m[:, 1], "s", color="0.45", mfc="none", ms=4, label="measured, 2e-3 mm³/mm²")
    a.plot(liq[:, 0], liq[:, 1], color="tab:red", lw=2.0, label="Foehncraft, 0.2e-3")
    a.plot(liq[:, 0], liq[:, 2], color="tab:red", lw=1.4, ls="--", label="Foehncraft, 2e-3")
    a.axvline(0.78, color="0.6", lw=0.8, ls=":")
    a.set_ylabel("liquid penetration (mm)")
    a.set_xlim(0.0, 1.25)
    a.set_ylim(bottom=0.0)
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="lower right")
    a.set_title("Liquid: projected liquid volume above the threshold", fontsize=10)

    b = fig.add_subplot(gs[1, 1])
    m = np.array(VAPOUR_PEN)
    b.plot(m[:, 0], m[:, 1], "o", color="k", ms=4, label="measured (schlieren)")
    b.plot(vap[:, 0], vap[:, 1], color="tab:blue", lw=2.0, label="Foehncraft, vapour mass fraction 0.1 %")
    b.axvline(0.78, color="0.6", lw=0.8, ls=":")
    b.set_xlabel("time after the start of injection (ms)")
    b.set_ylabel("vapour penetration (mm)")
    b.set_xlim(0.0, 1.25)
    b.set_ylim(bottom=0.0)
    b.grid(alpha=0.3)
    b.legend(fontsize=8, loc="lower right")
    b.set_title("Vapour, along the injector axis (end of injection dotted)", fontsize=10)

    path = os.path.join(out, "spray_g_g1.png")
    fig.savefig(path, dpi=100)
    plt.close(fig)
    print(path)
    for t, lo, hi in liq:
        print(f"liquid {t:6.3f} ms  {lo:6.2f} mm (0.2e-3)  {hi:6.2f} mm (2e-3)")
    for t, v in vap:
        print(f"vapour {t:6.3f} ms  {v:6.2f} mm")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    plot(out, run)


if __name__ == "__main__":
    main()
