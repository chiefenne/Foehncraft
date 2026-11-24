"""Field images of a 2D flow with a rotor (one cell layer along z): the
moving reference frame and the sliding mesh of gallery 94 to 96 and 98.

    pip install matplotlib h5py numpy [pillow]
    python3 -m tools.plot.rotor_fields OUT.png FIELDS [FIELDS ...] [options]

FIELDS is a `fields.vtkhdf` (or `fields_<step>.vtkhdf`) written by
`foehn flow` with `output.vtkhdf = true`; add `output.derived =
["vorticity"]` for the vorticity panel. Each cell is drawn as its face on
the lowest z plane (hexahedra, prisms and the polyhedra of an imprinted
interface alike), coloured by the cell value. One panel per quantity, or
with several FIELDS files (a time series) one panel per file with the
first quantity:

    speed       |u| (absolute velocity), with streamlines
    relative    |u - Omega x r| in the rotor (r < --rotor-radius; the
                stator is left grey): the velocity a blade sees, with
                streamlines
    vorticity   the z component of the vorticity
    pressure    p, from its minimum
    NAME        any cell dataset of the file, e.g. the fields of rotating
                machinery of `output.derived` (D183):
                relative_velocity_magnitude (with the streamlines of the
                relative velocity, w in the turning zones and c elsewhere),
                relative_flow_angle, flow_angle, rothalpy,
                relative_total_pressure, velocity_circumferential, ...;
                NAME@rotor draws it in the rotor only (r < --rotor-radius);
                k and nu_t of a turbulent run
    viscosity_ratio  nu_t / nu of a turbulent run (needs --nu)

Options: --quantities speed,vorticity (default), --omega (rad/s, for
`relative`), --nu (m²/s, for `viscosity_ratio`), --rotor-radius (m, draws the interface circle and bounds
`relative`), --outer-radius (m, masks the streamlines outside the vessel),
--mesh (draw the cell edges), --window x0,x1,y0,y1 (zoom), --title,
--labels (panel titles, comma separated), --limits lo:hi,... (colour limits
per panel), --no-streamlines, --dpi, --palette (a 256-colour PNG with
dithering, about a third of the size; needs pillow; for the docs).

Gallery 94 to 96 and 98 to 100 were drawn with it; examples/studies/mixer.py
has the commands of 98 to 100.
"""
import argparse
import math

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.tri import LinearTriInterpolator, Triangulation  # noqa: E402
from tools.plot.colormaps import RAINBOW, UNIFORM  # noqa: E402

LABELS = {
    "speed": ("|u| (m/s)", RAINBOW),
    "relative": ("|u − Ω×r| (m/s)", RAINBOW),
    "vorticity": ("ω_z (1/s)", UNIFORM),
    "pressure": ("p − p_min (Pa)", RAINBOW),
    "relative_velocity_magnitude": ("|w| (m/s)", RAINBOW),
    "frame_velocity_magnitude": ("|u| = Ω r (m/s)", RAINBOW),
    "relative_flow_angle": ("β (°)", UNIFORM),
    "flow_angle": ("α (°)", UNIFORM),
    "rothalpy": ("I − I_min (m²/s²)", RAINBOW),
    "relative_total_pressure": ("p + ½ρ|w|² − min (Pa)", RAINBOW),
    "velocity_circumferential": ("c_u (m/s)", UNIFORM),
    "relative_velocity_circumferential": ("w_u (m/s)", UNIFORM),
    "velocity_radial": ("c_r (m/s)", UNIFORM),
    "k": ("k (m²/s²)", "inferno"),
    "nu_t": ("ν_t (m²/s)", RAINBOW),
    "viscosity_ratio": ("ν_t / ν", RAINBOW),
}
SIGNED = ("vorticity", "velocity_circumferential", "relative_velocity_circumferential",
          "velocity_radial")
FROM_MIN = ("pressure", "rothalpy", "relative_total_pressure")


def read(path):
    """The lower-z faces of the cells (polygons), their centroids and the
    cell data."""
    with h5py.File(path, "r") as f:
        g = f["VTKHDF/internal"]
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    z0 = pts[:, 2].min()
    tol = 1e-9 + 1e-6 * (pts[:, 2].max() - z0)
    polys = []
    for c in range(len(off) - 1):
        ids = np.unique(conn[off[c]:off[c + 1]])
        q = pts[ids]
        q = q[np.abs(q[:, 2] - z0) < tol][:, :2]
        m = q.mean(axis=0)
        q = q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))]
        polys.append(q)
    cen = np.array([p.mean(axis=0) for p in polys])
    return polys, cen, data


def quantity(name, cen, data, omega, r_rotor, nu=None):
    name, _, region = name.partition("@")
    if region:
        if region != "rotor":
            raise SystemExit(f"unknown region {region} (only @rotor)")
        inside = np.hypot(cen[:, 0], cen[:, 1]) < r_rotor
        val = np.ma.masked_where(~inside, quantity(name, cen, data, omega, r_rotor, nu))
        return val - val.min() if name in FROM_MIN else val
    if name == "viscosity_ratio":
        if not nu or "nu_t" not in data:
            raise SystemExit("viscosity_ratio needs --nu and a turbulent run (nu_t)")
        return data["nu_t"] / nu
    u, v = data["velocity_x"], data["velocity_y"]
    if name == "speed":
        return np.hypot(u, v)
    if name == "relative":
        inside = np.hypot(cen[:, 0], cen[:, 1]) < r_rotor
        val = np.hypot(u + omega * cen[:, 1], v - omega * cen[:, 0])
        return np.ma.masked_where(~inside, val)
    if name == "vorticity":
        return data["vorticity_z"]
    if name in FROM_MIN and name in data:
        return data[name] - data[name].min()
    if name in data:
        return data[name]
    raise SystemExit(f"unknown quantity {name} (datasets: {', '.join(sorted(data))})")


def streamlines(ax, cen, data, omega, r_rotor, relative, window, outline):
    """Streamlines of the cell-centre velocity, interpolated on a grid (the
    points outside the fluid masked)."""
    u, v = data["velocity_x"].copy(), data["velocity_y"].copy()
    if relative == "file":
        u, v = data["relative_velocity_x"].copy(), data["relative_velocity_y"].copy()
        relative = False
    if relative:
        inside = np.hypot(cen[:, 0], cen[:, 1]) < r_rotor
        u = np.where(inside, u + omega * cen[:, 1], u)
        v = np.where(inside, v - omega * cen[:, 0], v)
    keep = np.ones(len(cen), bool)
    if relative:
        keep = np.hypot(cen[:, 0], cen[:, 1]) < r_rotor
    c, u, v = cen[keep], u[keep], v[keep]
    tri = Triangulation(c[:, 0], c[:, 1])
    # mask triangles that bridge a solid (longer edges than the cells around)
    x, y = c[tri.triangles, 0], c[tri.triangles, 1]
    edge = np.max(np.hypot(x - np.roll(x, 1, axis=1), y - np.roll(y, 1, axis=1)), axis=1)
    tri.set_mask(edge > 3.0 * np.median(edge))
    x0, x1, y0, y1 = window
    n = 300
    gx, gy = np.meshgrid(np.linspace(x0, x1, n), np.linspace(y0, y1, n))
    gu = LinearTriInterpolator(tri, u)(gx, gy)
    gv = LinearTriInterpolator(tri, v)(gx, gy)
    r = np.hypot(gx, gy)
    out = (r > outline) if outline is not None else np.zeros_like(r, bool)
    if relative:
        out |= r > r_rotor
    gu = np.ma.masked_where(out, gu)
    gv = np.ma.masked_where(out, gv)
    ax.streamplot(gx, gy, gu, gv, density=1.6, color=(1, 1, 1, 0.55), linewidth=0.5,
                  arrowsize=0.5, broken_streamlines=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("fields", nargs="+")
    ap.add_argument("--quantities", default="speed,vorticity")
    ap.add_argument("--omega", type=float, default=0.0)
    ap.add_argument("--nu", type=float, default=None,
                    help="kinematic viscosity (m²/s), for viscosity_ratio")
    ap.add_argument("--rotor-radius", type=float, default=0.0)
    ap.add_argument("--outer-radius", type=float, default=None,
                    help="radius of the vessel (masks the streamlines outside)")
    ap.add_argument("--mesh", action="store_true")
    ap.add_argument("--window", default=None)
    ap.add_argument("--title", default="")
    ap.add_argument("--labels", default=None, help="panel titles, comma separated")
    ap.add_argument("--limits", default=None,
                    help="colour limits per panel, e.g. 0:0.5,-40:40 (empty: automatic)")
    ap.add_argument("--no-streamlines", action="store_true")
    ap.add_argument("--dpi", type=int, default=170)
    ap.add_argument("--palette", action="store_true",
                    help="write a 256-colour PNG (dithered, needs pillow)")
    a = ap.parse_args()

    names = a.quantities.split(",")
    if len(a.fields) > 1:
        frames = [read(f) for f in a.fields]
        names = [names[0]] * len(frames)
    else:
        frames = [read(a.fields[0])] * len(names)
    allp = np.concatenate(frames[0][0])
    if a.window:
        window = [float(v) for v in a.window.split(",")]
    else:
        window = [allp[:, 0].min(), allp[:, 0].max(), allp[:, 1].min(), allp[:, 1].max()]
    titles = a.labels.split(",") if a.labels else None
    limits = a.limits.split(",") if a.limits else [""] * len(names)
    if len(a.fields) > 1 and len(limits) == 1:
        limits = limits * len(names)
    fig, axes = plt.subplots(1, len(names), figsize=(5.2 * len(names), 5.0),
                             squeeze=False, constrained_layout=True)
    for k, (ax, name) in enumerate(zip(axes[0], names)):
        polys, cen, data = frames[k]
        val = quantity(name, cen, data, a.omega, a.rotor_radius, a.nu)
        base = name.partition("@")[0]
        label, cmap = LABELS.get(base, (base, RAINBOW))
        if limits[k]:
            lo, hi = (float(v) for v in limits[k].split(":"))
        elif base in SIGNED:
            hi = np.percentile(np.abs(val), 98)
            lo = -hi
        elif base.endswith("_angle"):
            lo, hi = 0.0, 180.0
        else:
            lo, hi = 0.0, np.percentile(np.ma.compressed(val), 99.5)
        pc = PolyCollection(polys, array=val, cmap=cmap, edgecolors="face", linewidths=0.2)
        pc.set_clim(lo, hi)
        if a.mesh:
            pc.set_edgecolor((0, 0, 0, 0.35))
            pc.set_linewidth(0.25)
        ax.add_collection(pc)
        if name in ("speed", "relative", "relative_velocity_magnitude") and not a.no_streamlines:
            rel = "file" if name == "relative_velocity_magnitude" else name == "relative"
            streamlines(ax, cen, data, a.omega, a.rotor_radius, rel, window, a.outer_radius)
        if a.rotor_radius > 0:
            t = np.linspace(0, 2 * math.pi, 400)
            ax.plot(a.rotor_radius * np.cos(t), a.rotor_radius * np.sin(t), color="w",
                    lw=0.6, ls=(0, (4, 3)), alpha=0.7)
        ax.set_facecolor("#3a3a3a")
        ax.set_xlim(window[0], window[1])
        ax.set_ylim(window[2], window[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles[k] if titles else label, fontsize=11)
        fig.colorbar(pc, ax=ax, shrink=0.8, label=label)
    if a.title:
        fig.suptitle(a.title, fontsize=12)
    fig.savefig(a.out, dpi=a.dpi)
    if a.palette:
        from PIL import Image
        im = Image.open(a.out).convert("RGB")
        pal = im.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
        im.quantize(palette=pal, dither=Image.Dither.FLOYDSTEINBERG).save(a.out, optimize=True)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
