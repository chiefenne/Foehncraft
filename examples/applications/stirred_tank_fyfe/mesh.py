"""The 2 L lab stirred tank of Fyfe et al. (Org. Process Res. Dev. 2024, 28,
3661-3673, doi:10.1021/acs.oprd.4c00229, CC BY 4.0): a cylindrical vessel of
inner diameter T = 135 mm with a chamfered bottom, liquid height 150 mm, a
four-blade 45 degree pitched-blade turbine of tip diameter 45 mm and a
four-plate baffle insert. Dimensions from the paper's Supporting Information
(vessel, impeller, hub, shaft, clearance, baffle plates); the bottom chamfer
is read from the paper's published CFD plane. The geometry is theirs, not a
Foehncraft design, so the four-fold arrangement is shown with the reference.

    pip install gmsh
    python3 examples/applications/stirred_tank_fyfe/mesh.py                # level 1 (coarse)
    python3 examples/applications/stirred_tank_fyfe/mesh.py --level 2      # medium
    python3 examples/applications/stirred_tank_fyfe/mesh.py --level 0      # CI smoke mesh

Writes `fyfe_L<level>.msh` next to this script. Tetrahedra, no prism layers.
Boundaries: `impeller` (hub and blades inside the rotor zone), `shaft`,
`vessel` (side wall and bottom), `baffles`, `surface` (the liquid level, a
flat symmetry plane). Cell zones: `rotor`, `stator`. Lengths in m, z up from
the vessel bottom (the paper's plane files have z = 0 at the bottom and z
pointing down into the liquid).
"""
import argparse
import math
import os

import gmsh

R = 0.0675                 # vessel radius (T = 135 mm)
H = 0.150                  # liquid height
R_FLAT, CHAMFER_H = 0.035, 0.016   # flat bottom out to R_FLAT, then a cone up to the wall
R_TIP = 0.0225             # impeller tip radius (D = 45 mm)
CHORD = 0.015 / math.sin(math.radians(45.0))    # plate width: 15 mm projected blade height
T_BLADE = 0.0014           # blade thickness
Z_IMP = 0.0325             # blade mid-height: 25 mm clearance to the lower blade edge
R_HUB, Z_HUB0, Z_HUB1 = 0.006, 0.026, 0.038
R_SHAFT = 0.0035
N_BLADES = 4
R_ZONE, Z0_ZONE, Z1_ZONE = 0.027, 0.020, 0.050   # rotor zone cylinder (the paper's inner domain box)
R_COLUMN = 0.012           # the rotor zone goes on around the shaft to the liquid level
BAFFLE_R0, BAFFLE_W, BAFFLE_T = 0.035, 0.012, 0.005
Z_BAFFLE = 0.050           # baffle plates start here and reach the liquid level

LEVEL_SIZES = {  # bulk, rotor zone, blades, baffles, shaft and hub
    0: (0.010, 0.004, 0.0020, 0.003, 0.0025),
}


def sizes(level):
    if level in LEVEL_SIZES:
        return LEVEL_SIZES[level]
    base = (0.0065, 0.0030, 0.0012, 0.0024, 0.0018)
    f = 1.3 ** (level - 1)
    return tuple(s / f for s in base)


def blade(angle):
    """A flat 45 degree plate along +x, turned about z by `angle`: for a
    positive rotation about +z it pushes the liquid down."""
    occ = gmsh.model.occ
    b = occ.addBox(0.0, -0.5 * CHORD, -0.5 * T_BLADE, R_TIP, CHORD, T_BLADE)
    occ.rotate([(3, b)], 0, 0, 0, 1, 0, 0, math.radians(45.0))
    occ.translate([(3, b)], 0, 0, Z_IMP)
    occ.rotate([(3, b)], 0, 0, 0, 0, 0, 1, angle)
    return [(3, b)]


def geometry():
    occ = gmsh.model.occ
    prof = [(0, 0), (R_FLAT, 0), (R, CHAMFER_H), (R, H), (0, H)]
    p = [occ.addPoint(r, 0, z) for r, z in prof]
    loop = occ.addCurveLoop([occ.addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))])
    face = occ.addPlaneSurface([loop])
    tank = [e for e in occ.revolve([(2, face)], 0, 0, 0, 0, 0, 1, 2.0 * math.pi) if e[0] == 3]
    occ.remove([(2, face)], recursive=True)

    hub = occ.addCylinder(0, 0, Z_HUB0, 0, 0, Z_HUB1 - Z_HUB0, R_HUB)
    shaft = occ.addCylinder(0, 0, Z_HUB1, 0, 0, H - Z_HUB1 + 0.02, R_SHAFT)
    blades = []
    for k in range(N_BLADES):
        blades += blade(2.0 * math.pi * k / N_BLADES)
    impeller, _ = occ.fuse([(3, hub)], blades + [(3, shaft)])

    obstacles = []
    for k in range(4):        # on the +x, +y, -x, -y axes (the paper's YZ plane cuts two of them)
        b = occ.addBox(BAFFLE_R0, -0.5 * BAFFLE_T, Z_BAFFLE, BAFFLE_W, BAFFLE_T, H - Z_BAFFLE + 0.02)
        occ.rotate([(3, b)], 0, 0, 0, 0, 0, 1, 0.5 * math.pi * k)
        obstacles.append((3, b))
    fluid, _ = occ.cut(tank, impeller + obstacles)

    # rotor zone: a cylinder of four 90 degree sectors, cut planes midway
    # between the blades, so each lateral face has four sides and takes a
    # structured mesh (see interface())
    corners = [(0, Z0_ZONE), (R_ZONE, Z0_ZONE), (R_ZONE, Z1_ZONE), (R_COLUMN, Z1_ZONE),
               (R_COLUMN, H), (0, H)]
    p = [occ.addPoint(r, 0, z) for r, z in corners]
    loop = occ.addCurveLoop([occ.addLine(p[i], p[(i + 1) % len(p)]) for i in range(len(p))])
    face = occ.addPlaneSurface([loop])
    occ.rotate([(2, face)], 0, 0, 0, 0, 0, 1, math.pi / N_BLADES)
    sectors = []
    for k in range(N_BLADES):
        f = occ.copy([(2, face)])
        occ.rotate(f, 0, 0, 0, 0, 0, 1, 2.0 * math.pi * k / N_BLADES)
        sectors += [e for e in occ.revolve(f, 0, 0, 0, 0, 0, 1, 2.0 * math.pi / N_BLADES)
                    if e[0] == 3]
    occ.remove([(2, face)], recursive=True)
    out, children = occ.fragment(fluid, sectors)
    occ.synchronize()
    fluid_children = set(t for _, t in children[0])
    zone_children = set(t for c in children[1:] for _, t in c)
    rotor = sorted(fluid_children & zone_children)
    stator = sorted(fluid_children - zone_children)
    junk = sorted(zone_children - fluid_children)
    if junk:
        occ.remove([(3, t) for t in junk], recursive=True)
        occ.synchronize()
    return rotor, stator


def classify(rotor, stator):
    owners, count = {}, {}
    for name, vs in (("rotor", rotor), ("stator", stator)):
        for v in vs:
            for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
                owners.setdefault(s, set()).add(name)
                count[s] = count.get(s, 0) + 1
    groups = {n: [] for n in ("impeller", "shaft", "vessel", "baffles", "surface")}
    for s, own in owners.items():
        if count[s] == 2:
            continue
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, s)
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        rc, span = math.hypot(cx, cy), max(x1 - x0, y1 - y0)
        if z0 > H - 1e-6:
            groups["surface"].append(s)
        elif max(abs(x0), abs(x1), abs(y0), abs(y1)) < 1.01 * R_SHAFT and z0 > Z_HUB1 - 1e-6:
            groups["shaft"].append(s)
        elif own == {"rotor"}:
            groups["impeller"].append(s)
        elif span > 0.04:
            groups["vessel"].append(s)
        elif rc > BAFFLE_R0 - 0.01:
            groups["baffles"].append(s)
        else:
            raise SystemExit(f"surface {s} not classified: bbox {x0, y0, z0, x1, y1, z1}")
    for name, tags in groups.items():
        if not tags:
            raise SystemExit(f"no surfaces for {name}")
        gmsh.model.addPhysicalGroup(2, sorted(tags), name=name)
    gmsh.model.addPhysicalGroup(3, rotor, name="rotor")
    gmsh.model.addPhysicalGroup(3, stator, name="stator")
    return groups


def interface(rotor, stator, level):
    """Structured mesh of vertical rectangles on the lateral interface (the
    frame flux through it vanishes), see examples/applications/stirred_tank."""
    h_zone, h_col = sizes(level)[1], sizes(level)[4]
    rot = set(s for v in rotor for _, s in gmsh.model.getBoundary([(3, v)], oriented=False))
    sta = set(s for v in stator for _, s in gmsh.model.getBoundary([(3, v)], oriented=False))
    lateral = []
    for s in sorted(rot & sta):
        bb = gmsh.model.getBoundingBox(2, s)
        z0, z1 = bb[2], bb[5]
        if z1 - z0 < 1e-6:
            continue
        column = z0 > Z1_ZONE - 1e-6
        radius, h = (R_COLUMN, h_col) if column else (R_ZONE, h_zone)
        n_theta = max(3, math.ceil(2.0 * math.pi * radius / N_BLADES / h))
        n_z = max(2, math.ceil((z1 - z0) / h))
        for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
            cb = gmsh.model.getBoundingBox(1, c)
            vertical = cb[5] - cb[2] > 1e-6
            gmsh.model.mesh.setTransfiniteCurve(c, (n_z if vertical else n_theta) + 1)
        gmsh.model.mesh.setTransfiniteSurface(s, "Alternate")
        lateral.append(s)
    if len(lateral) != 2 * N_BLADES:
        raise SystemExit(f"expected {2 * N_BLADES} lateral interface faces, found {lateral}")


def size_fields(groups, level):
    h_bulk, h_zone, h_blade, h_obst, h_shaft = sizes(level)
    f = gmsh.model.mesh.field
    fields = []

    def near(surfs, h, dist_min, dist_max):
        d = f.add("Distance")
        f.setNumbers(d, "SurfacesList", surfs)
        f.setNumber(d, "Sampling", 40)
        t = f.add("Threshold")
        f.setNumber(t, "InField", d)
        f.setNumber(t, "SizeMin", h)
        f.setNumber(t, "SizeMax", h_bulk)
        f.setNumber(t, "DistMin", dist_min)
        f.setNumber(t, "DistMax", dist_max)
        fields.append(t)

    near(groups["impeller"], h_blade, 0.002, 0.015)
    near(groups["baffles"], h_obst, 0.002, 0.02)
    near(groups["shaft"], h_shaft, 0.002, 0.02)
    near(groups["vessel"], 1.5 * h_obst, 0.002, 0.02)
    c = f.add("Cylinder")
    f.setNumber(c, "Radius", R_ZONE + 0.004)
    f.setNumber(c, "XCenter", 0.0)
    f.setNumber(c, "YCenter", 0.0)
    f.setNumber(c, "ZCenter", 0.5 * (Z0_ZONE + Z1_ZONE))
    f.setNumber(c, "XAxis", 0.0)
    f.setNumber(c, "YAxis", 0.0)
    f.setNumber(c, "ZAxis", 0.5 * (Z1_ZONE - Z0_ZONE) + 0.004)
    f.setNumber(c, "VIn", h_zone)
    f.setNumber(c, "VOut", h_bulk)
    fields.append(c)
    m = f.add("Min")
    f.setNumbers(m, "FieldsList", fields)
    f.setAsBackgroundMesh(m)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_bulk)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5 * h_blade)


def build(level, path, threads, verbose=False):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 4 if verbose else 2)
    gmsh.option.setNumber("General.NumThreads", threads)
    gmsh.model.add("stirred_tank_fyfe")
    rotor, stator = geometry()
    groups = classify(rotor, stator)
    size_fields(groups, level)
    interface(rotor, stator, level)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(path)
    n = 0
    for v in rotor + stator:
        types, tags, _ = gmsh.model.mesh.getElements(3, v)
        n += sum(len(t) for t in tags)
    gmsh.finalize()
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--level", type=int, default=1, help="0 smoke, 1 coarse, 2 medium")
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args()
    path = os.path.join(a.out, f"fyfe_L{a.level}.msh")
    n = build(a.level, path, a.threads)
    print(f"{path}: {n} cells")


if __name__ == "__main__":
    main()
