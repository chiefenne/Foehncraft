"""Draw the two standard injector tips of the spray devices (D343, stage ID3
of docs/requirements/step28-injector-devices.md) and write them as binary
STL files in assets/injectors/.

Our own parametric geometry, for the picture: a round body that ends in a
domed tip, with straight holes drilled through the dome. Each hole is a
bore (a tube of the hole diameter and length along the hole axis) whose
mouth is cut out of the dome, so the holes show as dark openings in a
viewer. The surface is not closed at the hole mouths (the bore meets the
dome at a ring of its own): it is a picture, not a wall (D344).

Frame of every file (D343): the tip at the origin, the injector axis along
+z into the chamber, the azimuth zero along +x, millimetres. The hole
pattern of each file is the default pattern of `geometry =
"standard-gasoline"` and `"standard-diesel"`, and the numbers below are the
ones the solver uses for those defaults (crates/foehncraft-particles/src/
device.rs, `standard`); keep both in step.

Published numbers (quoted with their source; everything else is our choice
and says so):

- gasoline: 8 holes of 0.165 mm, 0.17 mm long, drill angle 37° from the
  axis: the ECN Spray G nozzle specification as quoted in D309 (Engine
  Combustion Network, Spray G pages; hole length 0.16 to 0.18 mm). Body
  diameter 6 mm, dome radius 4 mm, exit circle radius 1 mm: our choice (the
  exit radius is the one of `spray_g_g1_injectors.toml`).
- diesel: 7 holes of 0.11 mm, 0.8 mm long, 155° included spray angle
  (77.5° from the axis): our choice, typical of a passenger-car common-rail
  nozzle (no published nozzle is copied). The single axial hole of ECN
  Spray A (0.09 mm, D310) is the pattern `count = 1, polar = 0` on the same
  body. Body diameter 4 mm, a sac dome of radius 1.2 mm with a cone tangent
  to it at 50°, exit circle radius 0.65 mm: our choice.

    python3 tools/data/injector_standard.py          # writes assets/injectors/*.stl
    python3 tools/data/injector_standard.py --check  # fails if the files differ
"""

import math
import os
import struct
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "assets", "injectors")

# name: body diameter, body length, dome radius, holes, polar (deg),
# exit circle radius, hole diameter, hole length (all mm)
STANDARD = {
    "standard-gasoline": dict(body=6.0, length=8.0, dome=4.0, count=8, polar=37.0, exit=1.0, hole=0.165, bore=0.17),
    "standard-diesel": dict(body=4.0, length=8.0, dome=1.2, cone=50.0, count=7, polar=77.5, exit=0.65, hole=0.11, bore=0.8),
}


def dome_z(rho, s):
    """The tip: a spherical cap of radius `dome` through the tip (z = 0);
    with `cone` (degrees, the polar angle of the tangent point on the cap)
    a cone tangent to the cap carries on out to the body, otherwise the cap
    meets the body cylinder. z(ρ) <= 0 below the tip."""
    R = s["dome"]
    rho = np.asarray(rho, dtype=float)
    psi = s.get("cone")
    if psi is None:
        r = np.minimum(rho, R)
        return -(R - np.sqrt(np.maximum(R * R - r * r, 0.0)))
    psi = math.radians(psi)
    rt, zt = R * math.sin(psi), -R * (1 - math.cos(psi))
    cap = -(R - np.sqrt(np.maximum(R * R - np.minimum(rho, rt) ** 2, 0.0)))
    return np.where(rho <= rt, cap, zt - math.tan(psi) * (rho - rt))


def holes(s):
    out = []
    for k in range(s["count"]):
        phi = 2 * math.pi * k / s["count"]
        th = math.radians(s["polar"])
        r = s["exit"]
        exit_pt = np.array([r * math.cos(phi), r * math.sin(phi), float(dome_z(np.array(r), s))])
        n = np.array([math.sin(th) * math.cos(phi), math.sin(th) * math.sin(phi), math.cos(th)])
        out.append((exit_pt, n))
    return out


def tris_of_grid(P):
    """Two triangles per quad of a grid of points P[i, j] (outward when i
    runs up the axis and j counter-clockwise seen from +z)."""
    t = []
    for i in range(P.shape[0] - 1):
        for j in range(P.shape[1] - 1):
            a, b, c, d = P[i, j], P[i, j + 1], P[i + 1, j + 1], P[i + 1, j]
            t.append((a, b, c))
            t.append((a, c, d))
    return t


def tube(p0, n, radius, length, m=16):
    """The bore of a hole: a tube from the exit p0 back along -n."""
    a = np.cross(n, [0.0, 0.0, 1.0] if abs(n[2]) < 0.9 else [1.0, 0.0, 0.0])
    a /= np.linalg.norm(a)
    b = np.cross(n, a)
    ang = np.linspace(0, 2 * math.pi, m + 1)
    ring = np.array([radius * (math.cos(t) * a + math.sin(t) * b) for t in ang])
    P = np.array([p0 - length * n + ring, p0 + ring])
    # inward-facing (the bore is a hole in the body): reverse
    return [(c, b_, a_) for (a_, b_, c) in tris_of_grid(P)]


def build(s):
    rb = s["body"] / 2
    hs = holes(s)
    T = []
    # the body cylinder from z = -length to the dome edge
    zedge = float(dome_z(np.array(rb), s))
    nphi = 72
    phi = np.linspace(0, 2 * math.pi, nphi + 1)
    zs = np.linspace(-s["length"], zedge, 2)
    P = np.array([[[rb * math.cos(p), rb * math.sin(p), z] for p in phi] for z in zs])
    T += tris_of_grid(P)
    # a flat back
    for j in range(nphi):
        T.append((np.array([0, 0, -s["length"]]), P[0, j + 1], P[0, j]))
    # the dome on a fine polar grid, triangles in the hole mouths removed
    nr = 40
    rho = np.linspace(rb, 0.0, nr)
    D = np.array([[[r * math.cos(p), r * math.sin(p), float(dome_z(np.array(r), s))] for p in phi] for r in rho])
    for t in tris_of_grid(D):
        c = sum(t) / 3
        cut = False
        for p0, n in hs:
            d = c - p0
            # distance from the hole axis, near the exit
            dist = np.linalg.norm(d - (d @ n) * n)
            if dist < 0.5 * s["hole"] and abs(d @ n) < 4 * s["hole"]:
                cut = True
                break
        if not cut and np.linalg.norm(np.cross(t[1] - t[0], t[2] - t[0])) > 1e-14:
            T.append(t)
    for p0, n in hs:
        T += tube(p0, n, 0.5 * s["hole"], s["bore"])
    return T


def stl_bytes(T, name):
    head = f"Foehncraft {name}, mm, tip at origin, axis +z, azimuth 0 along +x".encode()[:80]
    out = bytearray(head.ljust(80, b" "))
    out += struct.pack("<I", len(T))
    for a, b, c in T:
        n = np.cross(b - a, c - a)
        ln = np.linalg.norm(n)
        n = n / ln if ln > 0 else n
        out += struct.pack("<12fH", *n, *a, *b, *c, 0)
    return bytes(out)


def main():
    check = "--check" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    bad = []
    for name, s in STANDARD.items():
        data = stl_bytes(build(s), name)
        path = os.path.join(OUT, f"{name}.stl")
        if check:
            if not os.path.exists(path) or open(path, "rb").read() != data:
                bad.append(path)
        else:
            with open(path, "wb") as f:
                f.write(data)
            print(f"{path}: {(len(data) - 84) // 50} triangles")
    if bad:
        print("out of date:", *bad)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
