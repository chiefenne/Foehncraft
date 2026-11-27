"""Spray G with its device drawn (injector devices D343 to D345; gallery
and recipe of their pull request).

spray_g_device.png: the 8 plumes of the Spray G G1 sector run (every parcel
mirrored about the plume plane and turned by k·45°, as in
`spray_g_metrics.full_injector`, kept in the case frame), coloured by the
parcel speed, under the standard gasoline tip placed by the device of
`spray_g_device_injectors.toml` (tool `spray.place`, the same numbers as a
run). Left an oblique view of the whole spray, right the first 5 mm at the
tip with the hole axes.

standard_injector_tips.png: the two standard tips of assets/injectors/
(gasoline 8 holes, diesel 7 holes), oblique, with their default hole axes.

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_g_g1.toml --set time.end_time=6e-4 --out output/spray_g_g1
    python3 -m tools.plot.spray_g_device OUT_DIR output/spray_g_g1
"""
import glob
import json
import math
import os
import re
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CASES = os.path.join(ROOT, "examples", "cases", "flow")
FOEHN = os.path.join(ROOT, "target", "release", "foehn")
SECTOR = 22.5  # degrees: the sector of the G1 case, plume 1 on φ = −SECTOR/2
PLUMES = 8


def call(tool, payload):
    r = subprocess.run(
        [FOEHN, "call", tool, "--json", "--quiet", "--input", "-"],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=True,
    )
    env = json.loads(r.stdout)
    return env.get("output", env.get("result", env))


def read_vtu(path):
    """Points, cells (lists of point indices), types and cell fields of an
    ASCII .vtu as written by the solver."""
    text = open(path).read()

    def array(name):
        m = re.search(rf'Name="{name}"[^>]*>(.*?)</DataArray>', text, re.S)
        if m is None:
            have = sorted(set(re.findall(r'Name="([^"]+)"', text)))
            raise SystemExit(f"{path}: no array `{name}`; it has {', '.join(have)}")
        return np.array(m.group(1).split(), dtype=float)

    pts = np.array(
        re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S).group(1).split(), dtype=float
    ).reshape(-1, 3)
    conn = array("connectivity").astype(int)
    off = array("offsets").astype(int)
    cells = np.split(conn, off[:-1])
    fields = {k: array(k).astype(int) for k in ("device", "part", "hole", "injector")}
    return pts, cells, array("types").astype(int), fields


def placed(injectors):
    """The placed device of an injector file: triangles (n, 3, 3) and the
    hole lines [(exit, end)], from `spray.place`."""
    out = tempfile.mkdtemp(prefix="devices_")
    call("spray.place", {"file": injectors, "out": out})
    pts, cells, types, _ = read_vtu(os.path.join(out, "spray", "devices.vtu"))
    tris = np.array([pts[c] for c, t in zip(cells, types) if t == 5])
    lines = [pts[c] for c, t in zip(cells, types) if t == 3]
    return tris, lines


def parcel_file(run):
    files = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))
    if not files:
        raise SystemExit(f"{run}/spray has no parcels_*.vtkhdf")
    return files[-1]


def parcels(path):
    import h5py

    with h5py.File(path, "r") as f:
        pd = f["VTKHDF"]["PointData"]
        if "velocity" not in pd:
            raise SystemExit(f"{path} has no point field `velocity`; it has {', '.join(sorted(pd))}")
        return f["VTKHDF"]["Points"][:], pd["velocity"][:], float(f.attrs["time"])


def full_injector(x, v):
    """The 8 plumes from the sector, in the case frame (the axis along +x,
    the azimuth φ = atan2(z, y)): each parcel and its mirror about the plume
    plane, turned by k·45°."""
    phi0 = -0.5 * math.radians(SECTOR)
    r = np.hypot(x[:, 1], x[:, 2])
    phi = np.arctan2(x[:, 2], x[:, 1])
    speed = np.linalg.norm(v, axis=1)
    xs, ss = [], []
    for sign in (1.0, -1.0):
        for k in range(PLUMES):
            a = phi0 + sign * (phi - phi0) + k * 2.0 * math.pi / PLUMES
            xs.append(np.column_stack([x[:, 0], r * np.cos(a), r * np.sin(a)]))
            ss.append(speed)
    return np.concatenate(xs), np.concatenate(ss)


def shade(tris, light=(0.4, -0.5, 0.75), base=(0.62, 0.64, 0.68)):
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)
    lt = np.array(light) / np.linalg.norm(light)
    i = 0.35 + 0.65 * np.abs(n @ lt)
    return np.clip(np.outer(i, base), 0, 1)


def draw_tip(ax, tris, lines, view, line_len=None):
    """Tip and hole axes in a 3D axis, coordinates mapped by `view`."""
    t = view(tris.reshape(-1, 3)).reshape(-1, 3, 3)
    pc = Poly3DCollection(t, facecolors=shade(t), edgecolors="none", linewidths=0)
    ax.add_collection3d(pc)
    for p in lines:
        a, b = p[0], p[-1]
        if line_len is not None:
            b = a + line_len * (b - a) / np.linalg.norm(b - a)
        q = view(np.array([a, b]))
        ax.plot(q[:, 0], q[:, 1], q[:, 2], color="k", lw=1.4)


def spray_figure(out, run):
    tris, lines = placed(os.path.join(CASES, "spray_g_device_injectors.toml"))
    x, v, t = parcels(parcel_file(run))
    X, S = full_injector(x, v)
    # the injector points down: draw (y, z, -x) in mm
    def view(p):
        p = np.asarray(p) * 1e3
        return np.column_stack([p[:, 1], p[:, 2], -p[:, 0]])

    fig = plt.figure(figsize=(13.0, 6.4), constrained_layout=True)
    norm = plt.Normalize(0, np.percentile(S, 99))
    cmap = plt.get_cmap("plasma")
    for col, (lim, title, glyph) in enumerate([
        (26.0, f"Spray G, 8 plumes at t = {t * 1e3:.2f} ms under the standard gasoline tip", 6.0),
        (5.0, "The first 5 mm: tip, hole axes and the young plumes", 2.0),
    ]):
        ax = fig.add_subplot(1, 2, col + 1, projection="3d")
        P = view(X)
        keep = (np.abs(P[:, 0]) < lim) & (np.abs(P[:, 1]) < lim) & (P[:, 2] > -lim)
        if col == 1:
            # thin out the dense core for the close-up
            keep &= np.arange(len(P)) % 2 == 0
        order = np.argsort(P[keep, 0] - P[keep, 1])
        ax.scatter(*P[keep][order].T, s=0.6 if col == 0 else 1.5, c=cmap(norm(S[keep][order])),
                   depthshade=False, linewidths=0)
        draw_tip(ax, tris, lines, view, line_len=glyph * 1e-3)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim * (1.6 if col == 0 else 1.2), 0.25 * lim * (0.6 if col == 0 else 1.6))
        ax.set_box_aspect((1, 1, 0.92))
        ax.view_init(elev=22, azim=-63)
        ax.set_xlabel("y (mm)")
        ax.set_ylabel("z (mm)")
        ax.set_zlabel("−x, along the injector axis (mm)")
        ax.set_title(title, fontsize=11)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=fig.axes, shrink=0.6, label="parcel speed (m/s)")
    fig.savefig(os.path.join(out, "spray_g_device.png"), dpi=100)
    plt.close(fig)


def tips_figure(out):
    from tools.data.injector_standard import STANDARD, build, holes

    fig = plt.figure(figsize=(12.0, 5.6), constrained_layout=True)
    titles = {
        "standard-gasoline": "standard-gasoline: 8 holes of 0.165 mm at 37°",
        "standard-diesel": "standard-diesel: 7 holes of 0.11 mm at 77.5°",
    }
    for col, (name, s) in enumerate(STANDARD.items()):
        ax = fig.add_subplot(1, 2, col + 1, projection="3d")
        tris = np.array(build(s))
        r = s["exit"] + 1.3  # a close-up of the holes
        keep = (np.hypot(tris[:, :, 0], tris[:, :, 1]).max(axis=1) < r) & (tris[:, :, 2].min(axis=1) > -r)
        lines = [np.array([p, p + n]) for p, n in holes(s)]
        draw_tip(ax, tris[keep], lines, lambda p: np.asarray(p), line_len=0.9)
        ax.set_xlim(-r, r)
        ax.set_ylim(-r, r)
        ax.set_zlim(-r, r)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=38, azim=-60)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_zlabel("z, into the chamber (mm)")
        ax.set_title(titles[name], fontsize=11)
    fig.savefig(os.path.join(out, "standard_injector_tips.png"), dpi=100)
    plt.close(fig)


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    spray_figure(out, run)
    tips_figure(out)


if __name__ == "__main__":
    main()
