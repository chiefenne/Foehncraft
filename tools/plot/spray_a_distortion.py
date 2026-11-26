"""ECN Spray A with the Corrected Distortion model (Nguyen, Dahms, Pickett
and Tagliante 2022; stage SA4, D316; gallery and recipe of the SA4 pull
request).

spray_a_distortion.png: left, the n-dodecane vapour mass fraction next to
the symmetry plane z = 0 at the end of the run with the droplets (at their
distance from the axis) as dots, the standard dynamic drag above the axis
and the Corrected Distortion model below it, with the liquid length of each
run and the measured one; middle, the liquid length (PLV 0.2·10⁻³
mm³/mm²) and the vapour penetration (0.1 %) of both runs over time against
the measured quasi-steady liquid lengths; right, the drag and the Sherwood
number of a distorted drop over those of a sphere (the example
`spray_a_distortion`), against the dynamic drag of the standard model.

    pip install matplotlib h5py numpy
    cargo run --release -p foehncraft-particles --example spray_a_distortion -- output/spray_a_distortion
    target/release/foehn flow examples/cases/flow/spray_a_distortion.toml --set spray.breakup.drag=dynamic --out output/spray_a_dynamic
    target/release/foehn flow examples/cases/flow/spray_a_distortion.toml --out output/spray_a_corrected
    python3 -m tools.plot.spray_a_distortion OUT_DIR output/spray_a_distortion output/spray_a_dynamic output/spray_a_corrected
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

from tools.plot.colormaps import RAINBOW  # noqa: E402
from tools.plot.vtkhdf import internal  # noqa: E402

# the nozzle of examples/cases/flow/spray_a_distortion_injectors.toml
NOZZLE_X = 0.0002
# measured quasi-steady liquid lengths of Spray A at 900 K, 22.8 kg/m³
# (docs/requirements/step28-spray-a.md, section 2): Mie 3 % (bkldaAL4) and
# the PLV threshold 0.2e-3 mm³/mm² (nozzle 210370)
MEASURED = [(10.4, "Mie scattering, 3 %"), (11.2, "PLV 0.2·10⁻³ mm³/mm²")]
FUEL = "Y_NC12H26"


def plane(path):
    """The cells of the first layer above z = zmin as quads in (x, y), and
    their cell fields."""
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        conn = g["Connectivity"][:]
        off = g["Offsets"][:]
        if not np.all(np.diff(off) == 8):
            raise SystemExit(f"{path}: hexahedra only")
        cells = conn.reshape(-1, 8)
        data = {k: g["CellData"][k][:] for k in g["CellData"]}
    if FUEL not in data:
        raise SystemExit(f"{path} has no cell field `{FUEL}`; it has {', '.join(sorted(data))}")
    xyz = pts[cells]
    z0 = pts[:, 2].min()
    lo = xyz[:, :, 2].min(axis=1)
    keep = np.abs(lo - z0) < 1e-12
    quads = []
    for c in xyz[keep]:
        q = c[np.abs(c[:, 2] - z0) < 1e-12][:, :2]
        m = q.mean(axis=0)
        quads.append(q[np.argsort(np.arctan2(q[:, 1] - m[1], q[:, 0] - m[0]))])
    return np.array(quads), {k: v[keep] for k, v in data.items()}


def parcels(run):
    files = glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf"))
    if not files:
        raise SystemExit(f"{run}/spray has no parcel files: run the case first (see the docstring)")
    with h5py.File(max(files), "r") as f:
        g = f["VTKHDF"]
        pd = g["PointData"]
        for c in ["diameter", "distortion"]:
            if c not in pd:
                raise SystemExit(f"{max(files)} has no point field `{c}`; it has {', '.join(sorted(pd))}")
        return g["Points"][:], pd["diameter"][:], dict(f.attrs)


def history(run):
    csv = np.genfromtxt(os.path.join(run, "spray.csv"), delimiter=",", names=True, comments=None, skip_header=2)
    out = {}
    # numpy drops the dots and dashes of the column names
    for key, col in [("t", "time"), ("ll", "sprayaliquid_length_1_m"), ("vp", "sprayavapour_penetration_m")]:
        if col not in csv.dtype.names:
            raise SystemExit(f"{run}/spray.csv has no column `{col}`; it has {', '.join(csv.dtype.names)}")
        out[key] = csv[col]
    return out


def plot(out, factors, runs):
    fig = plt.figure(figsize=(16.0, 5.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.75, 1.0, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    vmax = 0.0
    hist = {}
    for name, run, sign in runs:
        quads, d = plane(os.path.join(run, "fields.vtkhdf"))
        q = quads.copy()
        q[:, :, 1] *= sign
        y = d[FUEL] * 100.0
        vmax = max(vmax, np.percentile(y, 99.9))
        pc = PolyCollection(q * 1000.0, array=y, cmap=RAINBOW, edgecolors="face", linewidths=0.2)
        ax.add_collection(pc)
        x, dp, attrs = parcels(run)
        r = np.hypot(x[:, 1], x[:, 2])
        ax.scatter(x[:, 0] * 1000.0, sign * r * 1000.0, s=1.2, color="white", alpha=0.55, lw=0, zorder=3)
        h = history(run)
        hist[name] = h
        k = int(np.argmin(np.abs(h["t"] - attrs["time"])))
        ll = (NOZZLE_X + h["ll"][k]) * 1000.0
        ax.plot([ll, ll], [0.0, sign * 12.0], color="k", lw=1.6)
        ax.text(ll + 0.4, sign * 10.8, f"liquid length {h['ll'][k] * 1000:.1f} mm", fontsize=8,
                va="center", color="k")
        ax.text(1.0, sign * 10.8, name, fontsize=9.5, va="center", color="white", weight="bold")
    for c in ax.collections:
        if isinstance(c, PolyCollection):
            c.set_clim(0.0, vmax)
    fig.colorbar(ax.collections[0], ax=ax, label="n-dodecane vapour mass fraction (%)", shrink=0.85, pad=0.01)
    m = (NOZZLE_X * 1000.0 + MEASURED[1][0])
    ax.plot([m, m], [-12.0, 12.0], color="0.85", ls="--", lw=1.3, label=f"measured liquid length {MEASURED[1][0]} mm (PLV)")
    ax.axhline(0.0, color="0.9", lw=0.8)
    ax.set_xlim(0.0, 50.0)
    ax.set_ylim(-12.0, 12.0)
    ax.set_aspect("equal")
    ax.set_xlabel("distance from the nozzle wall (mm)")
    ax.set_ylabel("distance from the axis (mm)")
    ax.legend(fontsize=8, loc="lower right", framealpha=0.85)
    ax.set_title(f"ECN Spray A, 900 K, 6 MPa, t = {attrs['time'] * 1000:.1f} ms: vapour and droplets", fontsize=10)

    b = fig.add_subplot(gs[0, 1])
    colours = {"dynamic drag": "tab:blue", "Corrected Distortion": "crimson"}
    for name, h in hist.items():
        c = colours.get(name, "k")
        b.plot(h["t"] * 1000.0, h["ll"] * 1000.0, color=c, lw=1.7, label=f"liquid length, {name}")
        b.plot(h["t"] * 1000.0, h["vp"] * 1000.0, color=c, lw=1.4, ls="-.", label=f"vapour penetration, {name}")
    for val, lab in MEASURED:
        b.axhline(val, color="0.35", ls=":", lw=1.2)
        b.text(0.02, val + 0.6, f"measured {val} mm ({lab})", fontsize=7.5, color="0.25")
    b.set_xlabel("time after the start of injection (ms)")
    b.set_ylabel("distance from the nozzle (mm)")
    b.set_title("Liquid length (PLV 0.2·10⁻³ mm³/mm²) and vapour penetration (0.1 %)", fontsize=9.5)
    b.grid(alpha=0.3)
    b.legend(fontsize=7.5, loc="upper left")

    c = fig.add_subplot(gs[0, 2])
    f = np.genfromtxt(os.path.join(factors, "factors.csv"), delimiter=",", names=True)
    lam = 10.0
    for re, ls in [(20.0, ":"), (100.0, "-"), (500.0, "--")]:
        s = f[(f["re"] == re) & (f["lambda"] == lam)]
        c.plot(s["y"], s["drag_corrected"], color="tab:blue", ls=ls, lw=1.6,
               label=f"drag, Re {re:.0f}")
        c.plot(s["y"], s["sherwood"], color="crimson", ls=ls, lw=1.6, label=f"Sherwood, Re {re:.0f}")
    s = f[(f["re"] == 100.0) & (f["lambda"] == lam)]
    c.plot(s["y"], s["drag_dynamic"], color="0.3", lw=1.4, ls="-.", label="dynamic drag (standard)")
    c.axhline(1.0, color="0.6", lw=0.8)
    c.axvline(0.0, color="0.6", lw=0.8)
    c.set_xlabel("distortion y (prolate < 0 < flattened)")
    c.set_ylabel("over the solid sphere")
    c.set_title(f"A distorted drop (λ = μ_l/μ_g = {lam:.0f}, Sc = 2)", fontsize=10)
    c.grid(alpha=0.3)
    c.legend(fontsize=7.5, loc="upper left", ncol=1)
    fig.savefig(os.path.join(out, "spray_a_distortion.png"), dpi=100)
    plt.close(fig)


def main():
    if len(sys.argv) != 5:
        raise SystemExit(
            "usage: python3 -m tools.plot.spray_a_distortion OUT_DIR FACTORS_DIR RUN_DYNAMIC RUN_CORRECTED"
        )
    out, factors, dyn, cd = sys.argv[1:]
    os.makedirs(out, exist_ok=True)
    plot(out, factors, [("dynamic drag", dyn, 1.0), ("Corrected Distortion", cd, -1.0)])


if __name__ == "__main__":
    main()
