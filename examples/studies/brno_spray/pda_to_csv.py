#!/usr/bin/env python3
"""The PDA data of the Brno air-blast water spray as the files of the
Brno dispersion case (D322).

Source: E. Rácz, M. Malý, J. Jedelský et al., "Spatio-temporal analysis of
sprays by using Phase Doppler Anemometry data", data set, Brno University of
Technology, Zenodo record 17935932 (licence CC BY 4.0). Operating point
air-blast_1.28kgh_W_1.2bar_100C_1: water, 1.28 kg/h, 1.2 bar air
over-pressure, liquid at 100 °C at the atomizer inlet. Each raw file holds the
droplets that crossed the probe at one point (arrival time, axial velocity
LDA1, the velocity LDA4 along y, diameter).

    python3 examples/studies/brno_spray/pda_to_csv.py DIR

DIR holds the 90 .txt files of the operating point. The script writes

* examples/cases/flow/spray_brno_samples.csv: droplets of the plane
  z = 20 mm for `[[spray.injector]] samples` (D322), drawn by mass;
* examples/studies/brno_spray/data/brno_profiles.csv: the measured
  statistics per point and size class at z = 20, 40 and 60 mm;
* examples/studies/brno_spray/data/brno_gas_20mm.csv: the gas at z = 20 mm,
  taken from the droplets below 5 µm (mean axial velocity and k);
* examples/cases/flow/spray_brno_inlet.csv: the inlet table of the case
  (r, u_n, k, omega) from it, see `inlet` below.

How the plane is made (assumptions, stated in the gallery):

* the spray is taken as round: the points of the y traverse (and the axis)
  give the droplets of each radius, LDA4 times the sign of y is their radial
  velocity; the tangential velocity is left out (the x traverse measures it,
  its mean is below 2 m/s);
* a droplet stands for the liquid volume d³ per sampling time of its point
  and for the ring area of its radius (r ± 0.5 mm); the size dependence of
  the probe area is left out;
* the radius of a drawn droplet is spread uniformly over the area of its
  ring, so that the parcels do not start on circles.
"""

import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CLASSES = [(0.0, 5.0), (5.0, 10.0), (10.0, 20.0), (20.0, 30.0), (30.0, 50.0), (50.0, 1e9)]
N_SAMPLES = 6000


def read(directory):
    pts = []
    for f in sorted(glob.glob(os.path.join(directory, "*.txt"))):
        lines = open(f, encoding="latin-1").read().split("\n")
        try:
            x, y, z = (float(v.replace(" mm", "").replace(",", ".")) for v in lines[3].split(";")[1:4])
        except (ValueError, IndexError):
            continue
        rows = []
        for l in lines[6:]:
            c = l.split("\t")
            if len(c) >= 8:
                try:
                    rows.append([float(v.replace(",", ".")) for v in c[:8]])
                except ValueError:
                    pass
        if not rows:
            continue
        a = np.array(rows)
        ok = a[:, 7] > 0
        a = a[ok]
        pts.append(
            dict(x=x, y=y, z=z, t=a[:, 1].max() / 1e3, u=a[:, 3], w=a[:, 4], d=a[:, 7], name=os.path.basename(f))
        )
    return pts


def profiles(pts, out):
    with open(out, "w") as f:
        f.write("# Brno air-blast water spray, 1.28 kg/h, 1.2 bar, 100 C (Racz et al., Zenodo 17935932, CC BY 4.0)\n")
        f.write("# per point and size class: droplets per second, mean and rms of the axial velocity,\n")
        f.write("# liquid volume per second relative to the axis point of z = 20 mm (d^3 sum / sampling time)\n")
        f.write("# made by examples/studies/brno_spray/pda_to_csv.py\n")
        f.write("z_mm,x_mm,y_mm,d_lo_um,d_hi_um,n,rate_per_s,u_mean,u_rms,v_mean,v_rms,d10_um,d32_um,volume_rate\n")
        ref = None
        for p in sorted(pts, key=lambda p: (p["z"], p["x"], p["y"])):
            if p["z"] == 20 and p["x"] == 0 and p["y"] == 0 and ref is None:
                ref = (p["d"] ** 3).sum() / p["t"]
        for p in sorted(pts, key=lambda p: (p["z"], p["x"], p["y"])):
            for lo, hi in CLASSES + [(0.0, 1e9)]:
                s = (p["d"] >= lo) & (p["d"] < hi)
                n = int(s.sum())
                if n < 20:
                    continue
                d = p["d"][s]
                f.write(
                    f"{p['z']:g},{p['x']:g},{p['y']:g},{lo:g},{min(hi, 999):g},{n},{n / p['t']:.6g},"
                    f"{p['u'][s].mean():.5g},{p['u'][s].std():.5g},{p['w'][s].mean():.5g},{p['w'][s].std():.5g},"
                    f"{d.mean():.5g},{(d**3).sum() / (d**2).sum():.5g},{(d**3).sum() / p['t'] / ref:.6g}\n"
                )


def plane(pts, out_samples, out_gas):
    rng = np.random.default_rng(2026)
    ring = [p for p in pts if p["z"] == 20 and p["x"] == 0]
    by_r = {}
    for p in ring:
        by_r.setdefault(abs(p["y"]), []).append(p)
    r_all, d_all, u_all, v_all, w_all, lo_all, hi_all = [], [], [], [], [], [], []
    gas = []
    for r, ps in sorted(by_r.items()):
        lo, hi = max(r - 0.5, 0.0), r + 0.5
        area = np.pi * (hi * hi - lo * lo)
        for p in ps:
            sign = 1.0 if p["y"] > 0 else (-1.0 if p["y"] < 0 else 1.0)
            w = p["d"] ** 3 / p["t"] * area / len(ps)
            r_all.append(np.full(len(p["d"]), r))
            lo_all.append(np.full(len(p["d"]), lo))
            hi_all.append(np.full(len(p["d"]), hi))
            d_all.append(p["d"])
            u_all.append(p["u"])
            v_all.append(sign * p["w"])
            w_all.append(w)
        small = np.concatenate([p["u"][p["d"] < 5] for p in ps])
        small_w = np.concatenate([p["w"][p["d"] < 5] for p in ps])
        if len(small) >= 20:
            k = 0.5 * (small.var() + 2 * small_w.var())
            gas.append((r, small.mean(), k, len(small)))
    d = np.concatenate(d_all)
    u = np.concatenate(u_all)
    v = np.concatenate(v_all)
    w = np.concatenate(w_all)
    lo = np.concatenate(lo_all)
    hi = np.concatenate(hi_all)
    pick = rng.choice(len(d), size=N_SAMPLES, replace=True, p=w / w.sum())
    rr = np.sqrt(lo[pick] ** 2 + rng.random(N_SAMPLES) * (hi[pick] ** 2 - lo[pick] ** 2))
    with open(out_samples, "w") as f:
        f.write("# The droplets of the plane z = 20 mm of the Brno air-blast water spray (1.28 kg/h,\n")
        f.write("# 1.2 bar, liquid at 100 C): Racz et al., PDA data set, Brno University of Technology,\n")
        f.write("# Zenodo record 17935932, CC BY 4.0. Drawn by liquid volume from the y traverse\n")
        f.write(f"# ({N_SAMPLES} rows of equal mass) by examples/studies/brno_spray/pda_to_csv.py; SI units.\n")
        f.write("r, d, u_axial, u_radial\n")
        for i, j in enumerate(pick):
            f.write(f"{rr[i] * 1e-3:.6e}, {d[j] * 1e-6:.6e}, {u[j]:.5g}, {v[j]:.5g}\n")
    with open(out_gas, "w") as f:
        f.write("# The gas at z = 20 mm: droplets below 5 um of the y traverse (mean axial velocity, k from the\n")
        f.write("# variances of LDA1 and twice LDA4); Racz et al., Zenodo 17935932, CC BY 4.0\n")
        f.write("r_mm,u_gas,k,n\n")
        for r, ug, k, n in gas:
            f.write(f"{r:g},{ug:.5g},{k:.5g},{n}\n")
    inlet(gas, os.path.join(ROOT, "examples", "cases", "flow", "spray_brno_inlet.csv"))
    # rows of equal mass: a row stands for 1/d³ droplets, d32 = N / Σ 1/d
    return len(pick) / (1.0 / d[pick]).sum()


# the eddy viscosity of the round jet, nu_t = 0.028 U_c r_1/2 (the
# self-similar far field, Pope 2000, Turbulent Flows, section 5.1), with the
# measured U_c and half-width at z = 20 mm, for omega = k / nu_t
COFLOW = 1.0  # m/s, outside the measured jet: the gas the jet entrains


def inlet(gas, out):
    r = np.array([g[0] for g in gas]) * 1e-3
    u = np.array([g[1] for g in gas])
    k = np.array([g[2] for g in gas])
    uc = u[0]
    half = np.interp(-0.5 * uc, -u, r)
    nut = 0.028 * uc * half
    # beyond the last measured point: down to the coflow at 12 mm, then flat
    r = np.append(r, [0.012, 0.05])
    u = np.maximum(np.append(u, [COFLOW, COFLOW]), COFLOW)
    k = np.append(k, [1.0, 1.0])
    with open(out, "w") as f:
        f.write("# Inlet of spray_brno.toml at z = 20 mm: the gas velocity and k of the droplets below 5 um\n")
        f.write("# (Racz et al., Zenodo 17935932, CC BY 4.0); omega = k / nu_t with nu_t = 0.028 U_c r_1/2\n")
        f.write(f"# = {nut:.3g} m2/s (U_c = {uc:.4g} m/s, r_1/2 = {half * 1e3:.3g} mm); {COFLOW} m/s of coflow outside\n")
        f.write("# made by examples/studies/brno_spray/pda_to_csv.py\n")
        f.write("r, u_n, k, omega\n")
        for ri, ui, ki in zip(r, u, k):
            f.write(f"{ri:.4g}, {ui:.5g}, {ki:.5g}, {ki / nut:.5g}\n")


def main():
    pts = read(sys.argv[1])
    print(f"{len(pts)} points")
    profiles(pts, os.path.join(HERE, "data", "brno_profiles.csv"))
    d32 = plane(
        pts,
        os.path.join(ROOT, "examples", "cases", "flow", "spray_brno_samples.csv"),
        os.path.join(HERE, "data", "brno_gas_20mm.csv"),
    )
    print(f"plane z = 20 mm: d32 of the drawn droplets {d32:.2f} um")


if __name__ == "__main__":
    main()
