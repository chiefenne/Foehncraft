"""A vessel with three stirrers in 2D (one cell layer): three rotors with their
own blades, each inside its own interface circle, in a stator that is the
rest of the vessel. Each rotor turns at its own speed about its own axis
(the sliding mesh with several rotors by imprinting, D185).

    pip install gmsh
    python3 examples/meshes/three_mixer.py              # default size
    python3 examples/meshes/three_mixer.py --size 1.2e-3

Writes `three_mixer.msh` next to this script: the vessel is a stadium (two
half circles of radius 0.10 m about x = -0.12 and +0.12 m, joined by
straight walls), the rotors sit at:
- x = -0.12 m (zone `rotor_a`, three straight blades)
- x = +0.12 m (zone `rotor_b`, five straight blades)
- x = 0, y = -0.03 m (zone `rotor_c`, seven circular-arc blades)

Each meets the stator on its circle of radius 0.06 m with its own nodes
(64 rotor and 96 stator intervals). All blades are 5 mm thick and reach
from the hub (radius 10 mm) to 45 mm.

The curved blades of rotor C (`curved_blade`): the camber line is a circle
arc of radius R_CAMBER that leaves the hub radially and bends backward
(the tip trails an anticlockwise rotor by about 21 degrees). The blade is
the band between the two concentric arcs R_CAMBER -/+ T/2, so its thickness
is exactly constant (an offset of a spline through points is not), closed
at the tip by a half circle (two quarter arcs: a circle arc of exactly 180
degrees is not fixed by its start, centre and end) and at the root by a line
well inside the hub. Hub and blades are fused into one tool before the cut
(the roots of seven blades overlap inside the hub). A root that ends on the
hub circle instead leaves two thin cusps per blade where the root line and
the hub circle meet at a small angle: 98 faces above 70 degrees of
non-orthogonality (up to 86) in `foehn mesh`, against 70 at most here.

Boundaries: `blades_a`, `blades_b`, `blades_c` (blades and hubs, turning
with their rotors), `vessel` (the vessel wall), `a_rotor_interface`,
`a_stator_interface`, `b_rotor_interface`, `b_stator_interface`,
`c_rotor_interface`, `c_stator_interface`, `zmin`, `zmax`. Lengths in m.
"""
import argparse
import math
import os

import gmsh

R_VESSEL = 0.10      # radius of the half circles of the vessel
X_ROTOR_A = -0.12    # distance of rotor_a axis from the middle
X_ROTOR_B = 0.12     # distance of rotor_b axis from the middle
Y_ROTOR_C = -0.03    # y position of rotor_c (centered at x=0)
R_IFACE = 0.06       # interface radius
R_HUB = 0.010        # hub radius
R_TIP = 0.045        # blade tip radius
T_BLADE = 0.005      # blade thickness
R_CAMBER = 0.040     # camber radius of the curved blades
DEPTH = 0.005        # thickness of the cell layer
N_ROTOR = 64         # intervals of an interface circle, rotor side
N_STATOR = 96        # stator side
ROTORS = [
    ("a", X_ROTOR_A, 0, 3, False),
    ("b", X_ROTOR_B, 0, 5, False),
    ("c", 0, Y_ROTOR_C, 7, True),
]


def camber_angle(r, rc):
    """The angle phi on the camber circle (centre (R_HUB, -rc), radius rc)
    of the point at the distance r from the axis, on the branch that starts
    at the hub point (R_HUB, 0) (phi = 90 degrees) and runs outward with
    falling phi (r > R_HUB) or inward with growing phi (r < R_HUB)."""
    def dist(phi):
        return math.hypot(R_HUB + rc * math.cos(phi), -rc + rc * math.sin(phi))
    lo, hi = (-0.5 * math.pi, 0.5 * math.pi) if r > R_HUB else (0.5 * math.pi, 1.5 * math.pi)
    # dist rises from the hub point outward (falling phi) and falls inward
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if (dist(mid) > r) == (r > R_HUB):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def curved_blade(rc):
    """One circular-arc blade of constant thickness T_BLADE about the axis
    at the origin, along +x: the camber line is the arc of radius rc that
    leaves the hub radially at (R_HUB, 0) and bends to -y (backward curved
    for an anticlockwise rotor). The blade is the band between the
    concentric arcs rc - T/2 and rc + T/2, closed at the tip by a half
    circle of radius T/2 (no sharp corner, so no sliver cells), and at the
    root by a straight line inside the hub (the cut fuses blade and hub, so
    no tangent contact of blade and hub circle is left)."""
    occ = gmsh.model.occ
    h = 0.5 * T_BLADE
    cx, cy = R_HUB, -rc
    p_tip = camber_angle(R_TIP - h, rc)     # the cap reaches R_TIP
    p_root = camber_angle(0.4 * R_HUB, rc)  # both root corners inside the hub
    centre = occ.addPoint(cx, cy, 0)

    def on(r, phi):
        return occ.addPoint(cx + r * math.cos(phi), cy + r * math.sin(phi), 0)

    outer_root, outer_tip = on(rc + h, p_root), on(rc + h, p_tip)
    inner_root, inner_tip = on(rc - h, p_root), on(rc - h, p_tip)
    tip_centre = on(rc, p_tip)
    # the tip cap is a half circle: two quarter arcs, since a circle arc of
    # exactly 180 degrees is not defined by its start, centre and end
    tx, ty = math.sin(p_tip), -math.cos(p_tip)   # camber tangent, outward
    nose = occ.addPoint(cx + rc * math.cos(p_tip) + h * tx,
                        cy + rc * math.sin(p_tip) + h * ty, 0)
    curves = [
        occ.addCircleArc(outer_root, centre, outer_tip),
        occ.addCircleArc(outer_tip, tip_centre, nose),
        occ.addCircleArc(nose, tip_centre, inner_tip),
        occ.addCircleArc(inner_tip, centre, inner_root),
        occ.addLine(inner_root, outer_root),
    ]
    return occ.addPlaneSurface([occ.addCurveLoop(curves)])


def rotor(x0, y0, blades, curved=False):
    """Create a rotor: straight blades, or circular-arc blades of constant
    thickness (`curved`, camber radius R_CAMBER)."""
    occ = gmsh.model.occ
    disk = (2, occ.addDisk(x0, y0, 0, R_IFACE, R_IFACE))
    tools = [(2, occ.addDisk(x0, y0, 0, R_HUB, R_HUB))]

    for k in range(blades):
        angle = 2.0 * math.pi * k / blades
        if curved:
            s = curved_blade(R_CAMBER)
        else:
            s = occ.addRectangle(0.5 * R_HUB, -0.5 * T_BLADE, 0.0,
                                 R_TIP - 0.5 * R_HUB, T_BLADE)
        occ.rotate([(2, s)], 0, 0, 0, 0, 0, 1, angle)
        occ.translate([(2, s)], x0, y0, 0)
        tools.append((2, s))

    # fuse hub and blades first: the roots of seven curved blades overlap
    # each other inside the hub, and one tool body cuts cleanly
    hub, _ = occ.fuse(tools[:1], tools[1:])
    out, _ = occ.cut([disk], hub)
    return out


def stator():
    """Create the stator (vessel without rotor circles)."""
    occ = gmsh.model.occ
    # Two disks and a connecting rectangle
    a = (2, occ.addDisk(X_ROTOR_A, 0, 0, R_VESSEL, R_VESSEL))
    b = (2, occ.addDisk(X_ROTOR_B, 0, 0, R_VESSEL, R_VESSEL))
    mid = (2, occ.addRectangle(X_ROTOR_A, -R_VESSEL, 0,
                                X_ROTOR_B - X_ROTOR_A, 2 * R_VESSEL))
    vessel, _ = occ.fuse([a], [b, mid])
    # Holes for rotor interfaces
    holes = [(2, occ.addDisk(x0, y0, 0, R_IFACE, R_IFACE))
             for _, x0, y0, _, _ in ROTORS]
    out, _ = occ.cut(vessel, holes)
    return out


def on_circle(dim, tag, x0, y0, r):
    """True for a curve or surface that lies on the circle r about (x0, y0)."""
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    return (abs(xmax - x0 - r) < 1e-6 and abs(xmin - x0 + r) < 1e-6 and
            abs(ymax - y0 - r) < 1e-6 and abs(ymin - y0 + r) < 1e-6)


def which_rotor(dim, tag):
    """The rotor whose interface disk holds the entity, or None."""
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    for name, x0, y0, _, _ in ROTORS:
        if (xmin > x0 - R_IFACE - 1e-6 and xmax < x0 + R_IFACE + 1e-6 and
            ymin > y0 - R_IFACE - 1e-6 and ymax < y0 + R_IFACE + 1e-6):
            return name
    return None


def build(size, path):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("three_mixer")
    occ = gmsh.model.occ

    # Create rotors
    rotor_surfs = {}
    for name, x0, y0, n_blades, curved in ROTORS:
        rotor_surfs[name] = rotor(x0, y0, n_blades, curved)

    stator_surfs = stator()
    occ.synchronize()

    # Set transfinite curves on interface circles
    for name, x0, y0, _, _ in ROTORS:
        for s, n in [(rotor_surfs[name], N_ROTOR), (stator_surfs, N_STATOR)]:
            for c in gmsh.model.getBoundary(s, oriented=False):
                if on_circle(1, c[1], x0, y0, R_IFACE):
                    gmsh.model.mesh.setTransfiniteCurve(c[1], n + 1)

    surfaces = gmsh.model.getEntities(2)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.25 * size)

    # Finer mesh near the blades
    blades = []
    for name in rotor_surfs:
        for c in gmsh.model.getBoundary(rotor_surfs[name], oriented=False):
            is_on_interface = any(
                on_circle(1, c[1], x0, y0, R_IFACE)
                for _, x0, y0, _, _ in ROTORS
            )
            if not is_on_interface:
                blades.append(c[1])

    f = gmsh.model.mesh.field
    if blades:
        d = f.add("Distance")
        f.setNumbers(d, "CurvesList", blades)
        f.setNumber(d, "Sampling", 200)
        t = f.add("Threshold")
        f.setNumber(t, "InField", d)
        f.setNumber(t, "SizeMin", 0.3 * size)
        f.setNumber(t, "SizeMax", size)
        f.setNumber(t, "DistMin", 0.002)
        f.setNumber(t, "DistMax", 0.02)
        f.setAsBackgroundMesh(t)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 8)          # frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)

    # One layer of cells along z
    occ.extrude(surfaces, 0, 0, DEPTH, numElements=[1], recombine=True)
    occ.synchronize()

    vols = {"rotor_a": [], "rotor_b": [], "rotor_c": [], "stator": []}
    for e in gmsh.model.getEntities(3):
        r = which_rotor(*e)
        if r is None:
            vols["stator"].append(e[1])
        else:
            vols[f"rotor_{r}"].append(e[1])

    groups = {
        "blades_a": [], "blades_b": [], "blades_c": [], "vessel": [],
        "zmin": [], "zmax": [],
        "a_rotor_interface": [], "a_stator_interface": [],
        "b_rotor_interface": [], "b_stator_interface": [],
        "c_rotor_interface": [], "c_stator_interface": [],
    }

    seen = {}
    for name, vs in vols.items():
        for v in vs:
            for b in gmsh.model.getBoundary([(3, v)], oriented=False):
                seen.setdefault(b[1], []).append(name)

    for tag, owners in seen.items():
        bb = gmsh.model.getBoundingBox(2, tag)
        owner = owners[0]
        if bb[5] < 1e-6:
            groups["zmin"].append(tag)
        elif bb[2] > DEPTH - 1e-6:
            groups["zmax"].append(tag)
        else:
            # Check if on any interface
            iface_match = None
            for name, x0, y0, _, _ in ROTORS:
                if on_circle(2, tag, x0, y0, R_IFACE):
                    iface_match = name
                    break

            if iface_match:
                side = "stator" if owner == "stator" else "rotor"
                groups[f"{iface_match}_{side}_interface"].append(tag)
            elif owner == "stator":
                groups["vessel"].append(tag)
            else:
                # Blade of which rotor?
                rotor_name = owner.split("_")[1]  # "rotor_a" -> "a"
                groups[f"blades_{rotor_name}"].append(tag)

    for name, tags in groups.items():
        if tags:
            gmsh.model.addPhysicalGroup(2, tags, name=name)

    for name, vs in vols.items():
        gmsh.model.addPhysicalGroup(3, vs, name=name)

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)

    n = sum(len(tags) for vs in vols.values() for v in vs
            for tags in gmsh.model.mesh.getElements(3, v)[1])
    gmsh.finalize()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--size", type=float, default=2.5e-3, help="largest cell size (m)")
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    a = ap.parse_args()
    path = os.path.join(a.out, "three_mixer.msh")
    print(f"{path}: {build(a.size, path)} cells")


if __name__ == "__main__":
    main()
