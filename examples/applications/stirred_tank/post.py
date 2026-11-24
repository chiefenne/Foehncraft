"""Engineering numbers of the stirred tank from a finished run: the power
number from the impeller torque, the pumping number from the axial flow
through the discharge plane below the impeller, the share of the liquid
that moves slower than 5 % of the tip speed (the dead zone of a blend), and
the mean y+ on the walls.

    pip install h5py numpy
    python3 examples/applications/stirred_tank/post.py output/case [--rpm 90] [--json checks.json]

The run directory is the `--out` of `foehn flow`: it needs `boundaries.csv`
and `fields.vtkhdf` (output.vtkhdf = true).
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "python"))
from foehncraft.vtkhdf import internal  # noqa: E402  (fields.vtkhdf in one piece)

RHO = 998.2          # water at 20 °C (the material database)
D = 0.4              # impeller diameter (m), as mesh.py
Z_IMP = 1.0 / 3.0    # impeller mid-plane (m)
DZ_PLANE = 0.07      # discharge plane below the mid-plane (m): below the blades
T = 1.0              # tank diameter (m)


def boundaries(path):
    rows = {}
    with open(path, encoding="utf-8") as f:
        lines = [line for line in f if not line.startswith("#")]
    for row in csv.DictReader(lines):
        rows[row["boundary"]] = row
    return rows


def cells(path):
    """Centroids, volumes and the velocity of the cells (tetrahedra)."""
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts = g["Points"][:]
        con = g["Connectivity"][:]
        off = g["Offsets"][:]
        if not np.all(np.diff(off) == 4):
            raise SystemExit("post.py reads tetrahedral meshes only")
        tet = pts[con.reshape(-1, 4)]
        cen = tet.mean(axis=1)
        vol = np.abs(np.einsum("ij,ij->i", tet[:, 1] - tet[:, 0],
                               np.cross(tet[:, 2] - tet[:, 0], tet[:, 3] - tet[:, 0]))) / 6.0
        cd = g["CellData"]
        u = np.stack([cd["velocity_x"][:], cd["velocity_y"][:], cd["velocity_z"][:]], axis=1)
        walls = {}
        for name in f["VTKHDF"]:
            sub = f["VTKHDF"][name]
            if name not in ("internal", "Assembly") and "CellData" in sub and "y_plus" in sub["CellData"]:
                walls[name] = sub["CellData"]["y_plus"][:]
    return cen, vol, u, walls


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("--rpm", type=float, default=90.0)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()

    n = a.rpm / 60.0
    omega = 2.0 * math.pi * n
    tip = math.pi * n * D
    b = boundaries(a.run / "boundaries.csv")
    rotor = [w for w in ("impeller", "shaft") if w in b]
    torque = sum(float(b[w]["torque"]) for w in rotor)
    power = abs(torque) * omega
    n_p = power / (RHO * n**3 * D**5)

    cen, vol, u, walls = cells(a.run / "fields.vtkhdf")
    r = np.hypot(cen[:, 0], cen[:, 1])
    # downward flow through the disc r < D/2 at the discharge plane: the
    # volume integral of -u_z over a thin slab divided by its height
    z0, h = Z_IMP - DZ_PLANE, 0.01
    slab = (np.abs(cen[:, 2] - z0) < 0.5 * h) & (r < 0.5 * D)
    q = float(np.sum(-u[slab, 2] * vol[slab]) / h)
    n_q = q / (n * D**3)
    speed = np.linalg.norm(u, axis=1)
    dead = float(np.sum(vol[speed < 0.05 * tip]) / np.sum(vol))
    result = {
        "rpm": a.rpm,
        "reynolds": RHO * n * D**2 / 1.0016e-3,
        "torque_N_m": torque,
        "power_W": power,
        "power_number": n_p,
        "pumping_number": n_q,
        "flow_number_ratio": n_q / n_p if n_p > 0 else None,
        "dead_volume_fraction": dead,
        "mean_speed_over_tip": float(np.sum(speed * vol) / np.sum(vol) / tip),
        "y_plus_mean": {w: float(np.mean(y)) for w, y in walls.items()},
    }
    for k, v in result.items():
        print(f"{k:24s} {v}")
    if a.json:
        a.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
