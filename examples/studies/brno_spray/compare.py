#!/usr/bin/env python3
"""The Brno spray case against its measurements (D322).

    python3 examples/studies/brno_spray/compare.py OUT [--png FILE]

OUT is the output folder of `foehn flow examples/cases/flow/spray_brno.toml
--out OUT`. The streams of the final pass (spray/streams.vtkhdf) are cut by
the planes z = 40 and 60 mm; each crossing counts with the droplets per second
of its stream, so the statistics are flux-weighted like those of the phase
Doppler probe, which counts the droplets that cross it. The measured points of
both traverses are pooled by their distance from the axis.

Prints, per plane and size class, the mean axial velocity on rings of ± 1 mm
against the measured one, d32 and the liquid volume flux (each over its plane mean), and writes
OUT/brno_compare.csv (and the chart with --png).
"""

import argparse
import csv
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PLANES = [0.040, 0.060]
CLASSES = [(5.0, 10.0), (10.0, 20.0), (20.0, 30.0), (30.0, 50.0)]
RADII = np.arange(0.0, 14.5, 2.0)  # mm


def crossings(path, zc):
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        x = g["Points"][:]
        pd = g["PointData"]
        v = pd["velocity"][:]
        d = pd["diameter"][:]
        rate = pd["droplet_rate"][:]
        sid = pd["stream_id"][:]
    out = []
    starts = np.flatnonzero(np.r_[True, sid[1:] != sid[:-1]])
    ends = np.r_[starts[1:], len(sid)]
    for a, b in zip(starts, ends):
        z = x[a:b, 2]
        k = np.flatnonzero((z[:-1] < zc) & (z[1:] >= zc))
        if len(k) == 0:
            continue
        i = a + k[0]
        f = (zc - x[i, 2]) / (x[i + 1, 2] - x[i, 2])
        p = x[i] + f * (x[i + 1] - x[i])
        out.append(
            (
                np.hypot(p[0], p[1]) * 1e3,
                (d[i] + f * (d[i + 1] - d[i])) * 1e6,
                v[i, 2] + f * (v[i + 1, 2] - v[i, 2]),
                rate[i],
            )
        )
    return np.array(out)


def measured():
    rows = list(csv.DictReader(l for l in open(os.path.join(HERE, "data", "brno_profiles.csv")) if not l.startswith("#")))
    for r in rows:
        for k in r:
            r[k] = float(r[k])
        r["r_mm"] = np.hypot(r["x_mm"], r["y_mm"])
    return rows


def pooled(rows, z, r, lo, hi):
    sel = [m for m in rows if m["z_mm"] == z and abs(m["r_mm"] - r) < 0.01 and m["d_lo_um"] == lo and m["d_hi_um"] == min(hi, 999)]
    if not sel:
        return None
    n = np.array([m["n"] for m in sel])
    return dict(
        u=(n * [m["u_mean"] for m in sel]).sum() / n.sum(),
        d32=np.mean([m["d32_um"] for m in sel]),
        q=np.mean([m["volume_rate"] for m in sel]),
    )


HALF = 1.0  # mm: the rings of the model, r ± 1 mm (the points lie 2 mm apart)


def shape_norm(points):
    pts = [(r, p["q"]) for r, p in points if p is not None]
    r = np.array([p[0] for p in pts])
    q = np.array([p[1] for p in pts])
    return np.trapezoid(q * 2 * np.pi * r, r) / (np.pi * r[-1] ** 2)


def ring(c, r, lo, hi):
    inner, outer = max(r - HALF, 0.0), r + HALF
    s = (c[:, 0] >= inner) & (c[:, 0] < outer) & (c[:, 1] >= lo) & (c[:, 1] < hi)
    if s.sum() < 5:
        return None
    w = c[s, 3]
    d = c[s, 1]
    area = np.pi * (outer**2 - inner**2)
    return dict(
        u=(w * c[s, 2]).sum() / w.sum(),
        d32=(w * d**3).sum() / (w * d**2).sum(),
        q=(w * d**3).sum() / area,
        n=int(s.sum()),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--png")
    a = ap.parse_args()
    meas = measured()
    path = os.path.join(a.out, "spray", "streams.vtkhdf")
    table = []
    cuts = {zc: crossings(path, zc) for zc in PLANES}
    for zc in PLANES:
        z = round(zc * 1e3)
        c = cuts[zc]
        # the flux of the model and of the probe, each over its mean over
        # the plane (∫ q 2πr dr / π R², trapezoids on the points): the shape;
        # the probe area of the PDA is not known, and its flux is not
        # conserved from plane to plane
        q0_sim = shape_norm([(r, ring(c, r, 0.0, 1e9)) for r in RADII])
        q0_meas = shape_norm([(r, pooled(meas, z, r, 0.0, 1e9)) for r in RADII])
        print(f"z = {z} mm: {len(c)} streams cross")
        for r in RADII:
            for lo, hi in CLASSES + [(0.0, 1e9)]:
                m = pooled(meas, z, r, lo, hi)
                s = ring(c, r, lo, hi)
                if m is None or s is None:
                    continue
                row = dict(z_mm=z, r_mm=r, d_lo=lo, d_hi=min(hi, 999), u_meas=m["u"], u_sim=s["u"],
                           d32_meas=m["d32"], d32_sim=s["d32"], q_meas=m["q"] / q0_meas, q_sim=s["q"] / q0_sim, n_sim=s["n"])
                table.append(row)
                if hi > 999:
                    print(f"  r {r:4.1f} mm all: u {s['u']:6.2f} / {m['u']:6.2f} m/s, d32 {s['d32']:5.2f} / {m['d32']:5.2f} um, "
                          f"flux {row['q_sim']:5.3f} / {row['q_meas']:5.3f}")
                else:
                    print(f"  r {r:4.1f} mm {lo:g}-{hi:g} um: u {s['u']:6.2f} / {m['u']:6.2f} m/s ({s['n']} crossings)")
    with open(os.path.join(a.out, "brno_compare.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    if a.png:
        chart(table, a.png)


def chart(table, png):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colours = ["#2b83ba", "#1a9641", "#fdae61", "#d7191c"]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.2), dpi=100)
    for i, z in enumerate([40, 60]):
        a = ax[i, 0]
        for (lo, hi), col in zip(CLASSES, colours):
            rows = [t for t in table if t["z_mm"] == z and t["d_lo"] == lo and t["d_hi"] == hi]
            r = [t["r_mm"] for t in rows]
            a.plot(r, [t["u_sim"] for t in rows], "-", color=col, lw=2, label=f"{lo:g}-{hi:g} µm model")
            a.plot(r, [t["u_meas"] for t in rows], "o", color=col, ms=6, mfc="none", mew=1.6, label="measured")
        a.set_ylabel(f"z = {z} mm\naxial velocity (m/s)")
        a.grid(alpha=0.3)
        rows = [t for t in table if t["z_mm"] == z and t["d_hi"] == 999]
        r = [t["r_mm"] for t in rows]
        ax[i, 1].plot(r, [t["d32_sim"] for t in rows], "-", color="#5e3c99", lw=2, label="model")
        ax[i, 1].plot(r, [t["d32_meas"] for t in rows], "s", color="#5e3c99", mfc="none", mew=1.6, label="measured")
        ax[i, 1].set_ylabel("d32 (µm)")
        ax[i, 1].grid(alpha=0.3)
        ax[i, 2].plot(r, [t["q_sim"] for t in rows], "-", color="#e66101", lw=2, label="model")
        ax[i, 2].plot(r, [t["q_meas"] for t in rows], "^", color="#e66101", mfc="none", mew=1.6, label="measured")
        ax[i, 2].set_ylabel("liquid flux / its plane mean")
        ax[i, 2].grid(alpha=0.3)
    for a in ax[1]:
        a.set_xlabel("distance from the axis (mm)")
    ax[0, 0].legend(fontsize=8, ncol=2)
    ax[0, 1].legend(fontsize=8)
    ax[0, 2].legend(fontsize=8)
    ax[0, 0].set_title("droplet velocity by size class")
    ax[0, 1].set_title("Sauter mean diameter")
    ax[0, 2].set_title("liquid volume flux")
    fig.suptitle("Brno air-blast water spray: model from the measured plane z = 20 mm against the PDA data", fontsize=12)
    fig.tight_layout()
    fig.savefig(png)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
