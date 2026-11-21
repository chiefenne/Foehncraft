"""Compares a finished run of the Fyfe 2 L tank with the CFD planes of the
paper (Fyfe et al., OPRD 2024, SI / figshare doi:10.6084/m9.figshare.26405800;
their sliding-mesh k-omega SST run, a snapshot at t = 60 s). The plane files
are not in the repository (licence of the dataset; download them yourself).

    pip install h5py numpy scipy matplotlib
    python3 examples/applications/stirred_tank_fyfe/compare.py output/case \
        "YZ plane 200rpm baffled Viscosity 0.0009.csv" --out docs/gallery/images

Their plane x = 0 has z = 0 at the vessel bottom and z pointing down into the
liquid, so z_ours = -z_theirs and w_ours = -w_theirs. Their rotor phase is
unknown, so the plot shows our x = 0 plane (both sides, y > 0 and y < 0), and
the average over the rotor phase (planes through the axis at 8 angles).
Writes `fyfe_speed_planes.png` and `fyfe_profiles.png`.
"""
import argparse
import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

R, H, TIP = 0.0675, 0.150, 0.5 * 0.045 * 20.943951023931955


def paper_plane(path):
    with open(path, encoding="utf-8") as f:
        for _ in range(5):
            f.readline()
        hdr = [h.strip() for h in f.readline().split(",")][:-1]
        d = np.genfromtxt(f, delimiter=",", filling_values=np.nan)
    i = hdr.index("Velocity in Stn Frame u [ m s^-1 ]")
    y, z = d[:, 1], -d[:, 2]
    return y, z, d[:, i], d[:, i + 1], -d[:, i + 2]   # u (tangential at x=0), v (radial), w up


def ours(run):
    with h5py.File(Path(run) / "fields.vtkhdf", "r") as f:
        g = f["VTKHDF/internal"]
        pts, con = g["Points"][:], g["Connectivity"][:].reshape(-1, 4)
        cen = pts[con].mean(axis=1)
        cd = g["CellData"]
        u = np.stack([cd["velocity_x"][:], cd["velocity_y"][:], cd["velocity_z"][:]], 1)
    return cKDTree(cen), u


def sample(tree, u, y, z, phi):
    """Velocity (tangential, radial, axial) in the plane through the axis at angle phi."""
    e_r = np.array([np.cos(phi), np.sin(phi), 0.0])
    e_t = np.array([-np.sin(phi), np.cos(phi), 0.0])
    p = np.stack([y * e_r[0], y * e_r[1], z], -1)
    d, i = tree.query(p)
    ok = d < 0.004
    v = u[i]
    return (np.where(ok, v @ e_t, np.nan), np.where(ok, v @ e_r, np.nan), np.where(ok, v[:, 2], np.nan))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("plane", type=Path)
    ap.add_argument("--out", type=Path, default=Path("."))
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    py, pz, pu, pv, pw = paper_plane(a.plane)
    tree, u = ours(a.run)
    pspeed = np.sqrt(pu**2 + pv**2 + pw**2)

    # --- the plane x = 0: paper and ours side by side (speed / tip speed)
    gy, gz = np.meshgrid(np.linspace(-R, R, 330), np.linspace(0.0, H, 360))
    tu, tr, tw = sample(tree, u, gy.ravel(), gz.ravel(), np.pi / 2)   # plane x = 0: e_r = +y
    ms = np.sqrt(tu**2 + tr**2 + tw**2).reshape(gy.shape)
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 5.4), dpi=200, sharey=True)
    kw = dict(cmap="turbo", vmin=0, vmax=0.5)
    sc = ax[0].tripcolor(py, pz, pspeed / TIP, shading="gouraud", **kw)
    ax[0].set_title("Fyfe et al. 2024, CFD plane (sliding mesh)", fontsize=9)
    ax[1].pcolormesh(gy, gz, ms / TIP, shading="auto", rasterized=True, **kw)
    ax[1].set_title("Foehncraft, frozen rotor", fontsize=9)
    for x in ax:
        x.set_aspect("equal"); x.set_xlabel("y  (m)")
    ax[0].set_ylabel("height above the bottom  (m)")
    fig.suptitle("Speed / tip speed in the plane x = 0, 200 rpm, baffled", fontsize=10)
    fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02, label="speed / tip speed  (-)")
    fig.savefig(a.out / "fyfe_speed_planes.png"); plt.close(fig)

    # --- profiles: axial velocity on horizontal lines, tangential velocity too
    lines = [0.015, 0.060, 0.100]                       # below the blades, above, mid-liquid
    fig, axs = plt.subplots(2, 3, figsize=(11, 6.4), dpi=200, sharex=True)
    ys = np.linspace(-R + 0.001, R - 0.001, 200)
    stats = {}
    for j, zl in enumerate(lines):
        m = np.abs(pz - zl) < 0.002
        order = np.argsort(py[m])
        phis = np.linspace(0, 2 * np.pi, 9)[:-1]
        for row, (comp, lab) in enumerate(((2, "axial velocity w / tip speed"), (0, "tangential velocity / tip speed"))):
            ax = axs[row, j]
            paper = (pw if comp == 2 else pu)[m][order]
            ax.plot(py[m][order], paper / TIP, ".", ms=2.5, color="0.35", label="Fyfe et al. CFD (t = 60 s)")
            one = sample(tree, u, ys, zl + 0 * ys, np.pi / 2)[comp]
            ax.plot(ys, one / TIP, "-", color="C3", lw=1.2, label="Foehncraft, plane x = 0")
            avg = np.nanmean([sample(tree, u, ys, zl + 0 * ys, p)[comp] for p in phis], axis=0)
            ax.plot(ys, avg / TIP, "--", color="C0", lw=1.2, label="Foehncraft, phase average")
            ax.axhline(0, color="k", lw=0.4)
            if row == 0:
                ax.set_title(f"height {1000 * zl:.0f} mm", fontsize=9)
            if j == 0:
                ax.set_ylabel(lab, fontsize=8)
            if row == 1:
                ax.set_xlabel("y  (m)")
            # difference of the phase average with the paper on the common points
            yy = py[m][order]
            ok = np.abs(yy) < 0.06
            ref = np.interp(yy[ok], ys, avg)
            stats[f"{'w' if comp == 2 else 'u_theta'}_z{1000 * zl:.0f}mm"] = float(
                np.sqrt(np.nanmean((ref - paper[ok]) ** 2)) / TIP)
    axs[0, 0].legend(fontsize=7, loc="best")
    fig.suptitle("Fyfe 2 L tank, 200 rpm, baffled: velocity profiles on horizontal lines", fontsize=10)
    fig.tight_layout(); fig.savefig(a.out / "fyfe_profiles.png"); plt.close(fig)
    print(json.dumps(stats, indent=2))
    if a.json:
        a.json.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
