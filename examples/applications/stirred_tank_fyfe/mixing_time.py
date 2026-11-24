"""Mixing time of the tracer run of the Fyfe 2 L tank (`tracer.toml`): from
the fields `fields_<step>.vtkhdf` of the run it takes the volume-weighted
mean of the tracer, the largest deviation of any cell from it, and the first
time after which every cell stays within the tolerance (theta_95 for 5 %).
The measured value of the paper (Table 3, 200 rpm baffled) and the
blend-time correlation of Grenville are drawn as reference lines.

    python3 examples/applications/stirred_tank_fyfe/mixing_time.py RUN [--dt 0.05] [--tol 0.05]
        [--scalar tracer] [--out DIR] [--np 1.3]

Writes `mixing_time.json` and `fyfe_mixing_time.png` into --out (default RUN).
Tetrahedral meshes only (the volume of a cell from its four nodes).
"""
import argparse
import json
import re
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "python"))
from foehncraft.vtkhdf import internal  # noqa: E402  (fields.vtkhdf in one piece)

T_TANK, D_IMP, N_REV = 0.135, 0.045, 200.0 / 60.0   # m, m, 1/s
T_MEASURED = 18.7        # s, Fyfe et al. 2024, Table 3: 200 rpm, baffled, whole-reactor plateau time


def grenville(n_p):
    """theta_95 = 5.20 N_P^(-1/3) N^-1 (T/D)^2 (Grenville; source to verify)."""
    return 5.20 * n_p ** (-1.0 / 3.0) / N_REV * (T_TANK / D_IMP) ** 2


def read(path, scalar):
    with h5py.File(path, "r") as f:
        g = internal(f)
        pts, con = g["Points"][:], g["Connectivity"][:].reshape(-1, 4)
        c = g["CellData"][scalar][:]
    return pts, con, c


def volumes(pts, con):
    a, b, c, d = (pts[con[:, i]] for i in range(4))
    return np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("--dt", type=float, default=0.05, help="time step of the run (s)")
    ap.add_argument("--tol", type=float, default=0.05)
    ap.add_argument("--scalar", default="tracer")
    ap.add_argument("--np", type=float, default=1.3, dest="n_p", help="power number for the Grenville line")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    out = a.out or a.run
    files = sorted((int(m.group(1)), p) for p in a.run.glob("fields_*.vtkhdf") if (m := re.search(r"fields_(\d+)", p.name)))
    if not files:
        raise SystemExit(f"no fields_<step>.vtkhdf in {a.run}")
    rows, vol = [], None
    for step, p in files:
        pts, con, c = read(p, a.scalar)
        if vol is None:
            vol = volumes(pts, con)
        mean = float(np.sum(c * vol) / np.sum(vol))
        rows.append((step * a.dt, mean, float(c.min()), float(c.max())))
    t = np.array([r[0] for r in rows])
    final = rows[-1][1]
    dev = np.array([max(abs(r[3] / final - 1.0), abs(r[2] / final - 1.0)) for r in rows])
    ok = dev <= a.tol
    theta = None
    for i in range(len(t)):
        if ok[i:].all():
            theta = float(t[i])
            break
    res = {"tracer_mean": final, "mean_first": rows[0][1], "tolerance": a.tol, "theta": theta,
           "last_time": float(t[-1]), "measured_s": T_MEASURED, "grenville_s": grenville(a.n_p), "n_p": a.n_p,
           "times": t.tolist(), "max_deviation": dev.tolist()}
    (out / "mixing_time.json").write_text(json.dumps(res, indent=1))
    print(f"tracer mean {final:.5g} (first {rows[0][1]:.5g}); theta_{round(100 * (1 - a.tol))} = "
          f"{'not reached' if theta is None else format(theta, '.1f') + ' s'} (last file at {t[-1]:.1f} s); "
          f"measured {T_MEASURED} s, Grenville {grenville(a.n_p):.1f} s")
    fig, ax = plt.subplots(figsize=(6.4, 3.8), dpi=200)
    ax.semilogy(t, np.maximum(dev, 1e-4), color="C0", lw=1.6, label="largest cell deviation from the mean")
    ax.axhline(a.tol, color="0.4", lw=0.8, ls="--", label=f"{100 * a.tol:.0f} % band")
    ax.axvline(T_MEASURED, color="C3", lw=1.0, label="measured plateau time, 200 rpm baffled (Fyfe et al.)")
    ax.axvline(grenville(a.n_p), color="C2", lw=1.0, ls=":", label=f"Grenville, N_P = {a.n_p}")
    if theta is not None:
        ax.plot([theta], [a.tol], "o", color="C0")
        ax.annotate(f"θ95 = {theta:.1f} s", (theta, a.tol), xytext=(6, 8), textcoords="offset points", fontsize=8)
    ax.set_xlabel("time after the release (s)")
    ax.set_ylabel("max |c / c_mean − 1|")
    ax.set_title("Tracer mixing in the 2 L tank, 200 rpm, baffled (frozen rotor)", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "fyfe_mixing_time.png")


if __name__ == "__main__":
    main()
