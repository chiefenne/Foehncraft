"""A shell and tube heat exchanger as one conformal multi-region tet mesh:
hot air on the shell side, water on the tube side, steel tubes, tubesheets,
segmental baffles and housing, every interface shared node for node.

    pip install gmsh
    python3 examples/applications/shell_tube_hx/mesh.py               # level 2, 3 mm, about 1.0 M cells
    python3 examples/applications/shell_tube_hx/mesh.py --level 0     # smoke, about 0.25 M cells
    target/release/foehn mesh examples/applications/shell_tube_hx/mesh.toml

Geometry (mm, the axis is y): `shell_tube_hx.step`, a realistic exchanger
that is not from a real project. Housing: shell of inner radius 45 mm and
wall 6 mm with dished heads (y = -87 to 387 mm), two tube-side nozzles on
top (inner radius 20 mm, at y = -30 and 330 mm) and two shell-side nozzles
on the sides (inner radius 15 mm, x = -100 at y = 28 and x = +100 at
y = 266). Bundle: seven tubes 21 x 3 mm (six on a 30 mm circle and one in
the centre), 300 mm long, flush with two tubesheets 6 mm thick at y = 0 and
294 mm; three segmental baffles 6 mm thick at y = 71, 147 and 223 mm, cut
20 mm off the axis on alternating sides.

Method. The housing is read from the STEP file. Tubes, tubesheets and
baffles are rebuilt as exact OpenCASCADE primitives with the dimensions of
the STEP file (identical to 1e-12 mm): a Boolean of the STEP tubes against
the void inside the housing returns wrong pieces in OpenCASCADE 7.x, the
primitives do not. One fragment of a box, the housing, four nozzle plugs
and all internal parts gives every piece once with shared faces; each piece
takes its zone from the inputs that contain it (housing > tube bore > tube
wall > tubesheet > baffle), the box pieces inside the housing are the
fluids (between the tubesheets: shell side; beyond them: tube side), the
plugs and the box outside are dropped. The flat ends of the plugs close the
nozzles (the four inlet/outlet faces). Tetrahedra (HXT), Netgen optimised,
no prism layers; written in metres.

Cell zones: `shell_side`, `tube_side`, `tubes`, `tubesheets`, `baffles`,
`housing`. Boundaries: `tube_nozzle_front` (y = -30 mm), `tube_nozzle_rear`
(y = 330 mm), `shell_nozzle_front` (x = -100 mm), `shell_nozzle_rear`
(x = +100 mm), `housing_outer` (all outer housing faces, including the
nozzle rims). The interfaces between zones carry no physical group: the
flow solver finds them from the cell zones (recipe 156).
"""
import argparse, itertools, math, os
import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))
TUBES = [(30 * math.cos(math.radians(a)), 30 * math.sin(math.radians(a))) for a in (0, 60, 120, 180, 240, 300)] + [(0.0, 0.0)]
R_OUT, R_IN, LENGTH, R_SHELL, THICK = 10.5, 7.5, 300.0, 45.0, 6.0
BAFFLES = [(71.0, -20.0, +1), (147.0, 20.0, -1), (223.0, -20.0, +1)]  # y start, cut z, +1 keeps z > cut
CAPS = {  # centre of the nozzle end, outward normal, plug radius (between inner and outer radius)
    "tube_nozzle_front": ((0, -30, 100), (0, 0, 1), 23),
    "tube_nozzle_rear": ((0, 330, 100), (0, 0, 1), 23),
    "shell_nozzle_front": ((-100, 28, 23), (-1, 0, 0), 18),
    "shell_nozzle_rear": ((100, 266, -23), (1, 0, 0), 18),
}
ZONES = ("shell_side", "tube_side", "tubes", "tubesheets", "baffles", "housing")


def build(step, size, out):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 2)
    occ = gmsh.model.occ
    solids = [e for e in occ.importShapes(step) if e[0] == 3]
    occ.synchronize()
    housing = [t for d, t in solids if gmsh.model.getEntityName(3, t).endswith("HOUSING")]
    assert len(housing) == 1, "no single HOUSING solid in the STEP file"
    occ.remove([e for e in solids if e[1] != housing[0]], recursive=True)

    role = {housing[0]: "housing"}
    role[occ.addBox(-150, -150, -150, 300, 700, 300)] = "box"
    for p, n, r in CAPS.values():
        role[occ.addCylinder(*p, *(5 * x for x in n), r)] = "plug"
    for x, z in TUBES:
        role[occ.addCylinder(x, 0, z, 0, LENGTH, 0, R_OUT)] = "tube_wall"
        role[occ.addCylinder(x, 0, z, 0, LENGTH, 0, R_IN)] = "bore"
    for y0 in (0.0, LENGTH - THICK):
        role[occ.addCylinder(0, y0, 0, 0, THICK, 0, R_SHELL)] = "tubesheet"
    for y0, zc, side in BAFFLES:
        disc = occ.addCylinder(0, y0, 0, 0, THICK, 0, R_SHELL)
        half = occ.addBox(-60, y0 - 1, zc if side > 0 else -60, 120, THICK + 2, 60 - zc if side > 0 else zc + 60)
        piece, _ = occ.intersect([(3, disc)], [(3, half)])
        role[piece[0][1]] = "baffle"

    inputs = [(3, t) for t in role]
    _, pieces = occ.fragment(inputs, [])
    occ.synchronize()
    owners = {}
    for e, ps in zip(inputs, pieces):
        for _, t in ps:
            owners.setdefault(t, set()).add(role[e[1]])
    drop = [t for t, r in owners.items() if "plug" in r or (r == {"box"} and occ.getBoundingBox(3, t)[0] < -149)]
    occ.remove([(3, t) for t in drop], recursive=True)
    occ.synchronize()
    zone = {z: [] for z in ZONES}
    for t, r in owners.items():
        if t in drop:
            continue
        if "housing" in r:
            z = "housing"
        elif "bore" in r:
            z = "tube_side"
        elif "tube_wall" in r:
            z = "tubes"
        elif "tubesheet" in r:
            z = "tubesheets"
        elif "baffle" in r:
            z = "baffles"
        elif r == {"box"}:
            bb = occ.getBoundingBox(3, t)
            z = "tube_side" if bb[4] < 1e-3 or bb[1] > LENGTH - 1e-3 else "shell_side"
        else:
            raise RuntimeError(f"piece {t} with owners {r}")
        zone[z].append(t)
    print("zone         pieces   volume (L)")
    for z, tags in zone.items():
        print(f"{z:12s} {len(tags):6d}   {sum(occ.getMass(3, t) for t in tags) / 1e6:.6f}")
        gmsh.model.addPhysicalGroup(3, tags, name=z)

    faces = {z: {abs(s) for _, s in gmsh.model.getBoundary([(3, t) for t in v], combined=False, oriented=False)} for z, v in zone.items()}
    print("interface                  area (cm²)")
    for a, b in itertools.combinations(ZONES, 2):
        shared = faces[a] & faces[b]
        if shared:
            print(f"{a:11s} | {b:11s} {sum(occ.getMass(2, s) for s in shared) / 100:9.2f}")
    assert not faces["shell_side"] & faces["tube_side"], "the two fluids touch"

    outer = gmsh.model.getBoundary([(3, t) for v in zone.values() for t in v], combined=True, oriented=False)
    caps = {k: [] for k in CAPS}
    wall = []
    for _, s in outer:
        c = occ.getCenterOfMass(2, s)
        bb = occ.getBoundingBox(2, s)
        ext = max(bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2])
        hit = [k for k, (p, n, r) in CAPS.items() if math.dist(c, p) < 1e-3 and ext < 2 * r - 5]
        (caps[hit[0]] if hit else wall).append(s)
    for k, v in caps.items():
        assert len(v) == 1, (k, v)
        gmsh.model.addPhysicalGroup(2, v, name=k)
    gmsh.model.addPhysicalGroup(2, wall, name="housing_outer")

    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size / 6)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.ScalingFactor", 0.001)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(out)
    total = 0
    for z, tags in zone.items():
        n = sum(len(gmsh.model.mesh.getElements(3, t)[1][0]) for t in tags)
        total += n
        print(f"{z:12s} {n:9d} tetrahedra")
    print(f"total        {total:9d} tetrahedra -> {out}")
    gmsh.finalize()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--level", type=int, default=2, help="0 smoke ... 3 fine; each level divides the cell size by 1.3")
    ap.add_argument("--step", default=os.path.join(HERE, "shell_tube_hx.step"))
    a = ap.parse_args()
    size = 3.0 * 1.3 ** (2 - a.level)  # mm: level 2 = 3 mm
    build(a.step, size, os.path.join(HERE, f"hx_L{a.level}.msh"))
