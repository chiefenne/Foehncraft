"""A five-hole nozzle described by its device (injector devices, D340 to
D342; gallery and recipe of their pull request).

spray_multihole_device.png: left, an oblique view of the box with the
parcels of the last parcel file coloured by the hole they came from, the
device axis (black), the hole axes (one short line per hole from its exit)
and the mounting point; right, the view along the device axis from the
chamber side, which shows the ring of five holes 72° apart and their plumes.
The holes and the frame come from the tool `spray.devices` (the same
numbers as the run).

    pip install matplotlib h5py numpy
    target/release/foehn flow examples/cases/flow/spray_multihole_device.toml --out output/spray_multihole_device
    python3 -m tools.plot.spray_multihole_device OUT_DIR output/spray_multihole_device
"""
import glob
import json
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INJECTORS = os.path.join(ROOT, "examples", "cases", "flow", "spray_multihole_device_injectors.toml")
COLOURS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]


def devices(foehn, path):
    """The output of `spray.devices` for an injector or case file."""
    r = subprocess.run(
        [foehn, "call", "spray.devices", "--json", "--quiet", "--input", "-"],
        input=json.dumps({"file": path}),
        capture_output=True,
        text=True,
        check=True,
    )
    env = json.loads(r.stdout)
    return env.get("output", env.get("result", env))


def parcels(path):
    import h5py

    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        pd = g["PointData"]
        if "injector" not in pd:
            raise SystemExit(f"{path}: no field `injector`; the fields are {sorted(pd.keys())}")
        return g["Points"][:], pd["injector"][:], dict(f.attrs)


def main():
    out, run = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    foehn = os.path.join(ROOT, "target", "release", "foehn")
    dev = devices(foehn, INJECTORS)
    d = dev["devices"][0]
    files = sorted(glob.glob(os.path.join(run, "spray", "parcels_*.vtkhdf")))
    if not files:
        raise SystemExit(f"no parcel files in {run}/spray")
    x, inj, attrs = parcels(files[-1])
    inj = np.asarray(inj).astype(int).ravel()
    mount, e1, e2, axis = (np.array(d[k]) for k in ("mount", "e1", "e2", "axis"))

    fig = plt.figure(figsize=(13.0, 6.0), constrained_layout=True)
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    for k in range(len(d["holes"])):
        m = inj == k
        # the box's y is up: draw (x, z, y)
        ax.scatter(x[m, 0], x[m, 2], x[m, 1], s=3, color=COLOURS[k % 5], depthshade=False, label=f"hole {k + 1}")
    a = mount + 0.03 * axis
    ax.plot([mount[0], a[0]], [mount[2], a[2]], [mount[1], a[1]], color="k", lw=2)
    for h in d["holes"]:
        p, n = np.array(h["exit"]), np.array(h["direction"])
        q = p + 0.012 * n
        ax.plot([p[0], q[0]], [p[2], q[2]], [p[1], q[1]], color="k", lw=1.2)
    ax.scatter([mount[0]], [mount[2]], [mount[1]], color="k", s=30)
    ax.set_xlim(0, 0.12)
    ax.set_ylim(0, 0.1)
    ax.set_zlim(0, 0.1)
    ax.set_box_aspect((1.2, 1.0, 1.0))
    ax.view_init(elev=18, azim=-58)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_zlabel("y (m)")
    ax.legend(fontsize=8, loc="upper left", markerscale=3)
    ax.set_title(f"{len(x)} parcels at t = {attrs.get('time', 0.0):.3f} s, nozzle tilted 20°", fontsize=11)

    ax = fig.add_subplot(1, 2, 2)
    r = x - mount
    u, v = r @ e1, r @ e2
    for k in range(len(d["holes"])):
        m = inj == k
        ax.scatter(u[m] * 1e3, v[m] * 1e3, s=3, color=COLOURS[k % 5])
    for k, h in enumerate(d["holes"]):
        p = np.array(h["exit"]) - mount
        n = np.array(h["direction"])
        pu, pv = p @ e1 * 1e3, p @ e2 * 1e3
        # the hole axis, drawn 8 mm long from the exit
        ax.plot([pu, pu + 8 * (n @ e1) / np.hypot(n @ e1, n @ e2)], [pv, pv + 8 * (n @ e2) / np.hypot(n @ e1, n @ e2)],
                color="k", lw=1.5)
        m = inj == k
        if m.any():
            cu, cv = np.median(u[m]) * 1e3, np.median(v[m]) * 1e3
            far = np.hypot(cu, cv)
            ax.text(cu * (far + 14) / far, cv * (far + 14) / far, f"hole {h['name']}", ha="center", va="center",
                    fontsize=10, color=COLOURS[k % 5], fontweight="bold")
    ax.plot([0], [0], "k+", ms=12)
    ax.set_xlim(-60, 62)
    ax.set_ylim(-52, 62)
    ax.set_aspect("equal")
    ax.set_xlabel("along e1 (mm)")
    ax.set_ylabel("along e2 (mm)")
    ax.set_title("Seen along the nozzle axis: five holes 72° apart, 30° from the axis", fontsize=11)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out, "spray_multihole_device.png"), dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    main()
