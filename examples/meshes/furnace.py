"""A box furnace with a row of cooled tubes under the roof (radiation R1a,
D255): a tetrahedral Gmsh mesh for the gallery case of the radiation call.

    pip install gmsh
    python3 examples/meshes/furnace.py                 # furnace.msh, default size
    python3 examples/meshes/furnace.py --size 0.08 --out furnace_fine.msh

The chamber is 4 m long (x), 2 m deep (y) and 2.5 m high (z). Seven
horizontal tubes of 0.16 m diameter run through the chamber along y under
the roof (z = 2.1 m, x = 0.9 to 3.3 m). The burner is a disc of 0.5 m
diameter in the front wall (x = 0) at mid-depth, 0.9 m above the floor; the
flue gas leaves through an opening in the back wall (x = 4 m) under the
roof.

Boundaries: `burner` (the disc), `front` (the rest of the front wall),
`walls` (floor, roof, side walls and the rest of the back wall), `tubes`
(the tube surfaces), `outlet` (the opening). Lengths in m. The cells are
`--size` in the chamber and a third of it on the tubes.
"""
import argparse
import os

import gmsh

LENGTH, DEPTH, HEIGHT = 4.0, 2.0, 2.5
TUBE_R = 0.08
TUBE_Z = 2.1
TUBE_X = [0.9 + 0.4 * i for i in range(7)]
BURNER = (DEPTH / 2, 0.9, 0.25)       # y, z of the centre, radius
OUTLET = (0.6, 1.4, 1.75, 2.35)       # y0, y1, z0, z1 in the back wall


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--size", type=float, default=0.12, help="cell size in the chamber (m)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "furnace.msh"))
    a = ap.parse_args()

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("furnace")
    occ = gmsh.model.occ
    box = occ.addBox(0, 0, 0, LENGTH, DEPTH, HEIGHT)
    tubes = [occ.addCylinder(x, -0.1, TUBE_Z, 0, DEPTH + 0.2, 0, TUBE_R) for x in TUBE_X]
    chamber, _ = occ.cut([(3, box)], [(3, t) for t in tubes])
    # the burner disc and the outlet opening as faces of their own
    yb, zb, rb = BURNER
    disc = occ.addDisk(0, yb, zb, rb, rb, zAxis=[1, 0, 0], xAxis=[0, 1, 0])
    y0, y1, z0, z1 = OUTLET
    corners = [(y0, z0), (y1, z0), (y1, z1), (y0, z1)]
    pts = [occ.addPoint(LENGTH, y, z) for y, z in corners]
    lines = [occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
    opening = occ.addPlaneSurface([occ.addCurveLoop(lines)])
    occ.fragment(chamber, [(2, disc), (2, opening)])
    occ.synchronize()

    volumes = gmsh.model.getEntities(3)
    assert len(volumes) == 1, volumes
    eps = 1e-6
    groups = {"burner": [], "front": [], "walls": [], "tubes": [], "outlet": []}
    for dim, tag in gmsh.model.getBoundary(volumes, oriented=False):
        x0, ya, za, x1, yb_, zb_ = gmsh.model.getBoundingBox(dim, tag)
        if x1 < eps:
            r = max(ya - (yb - rb), (yb + rb) - yb_, za - (zb - rb), (zb + rb) - zb_)
            groups["burner" if abs(r) < 1e-3 and yb_ - ya < 2 * rb + 1e-3 else "front"].append(tag)
        elif x0 > LENGTH - eps:
            inside = ya > y0 - eps and yb_ < y1 + eps and za > z0 - eps and zb_ < z1 + eps
            groups["outlet" if inside else "walls"].append(tag)
        elif zb_ - za < 2 * TUBE_R + 1e-3 and za > TUBE_Z - TUBE_R - 1e-3 and x1 - x0 < 2 * TUBE_R + 1e-3:
            groups["tubes"].append(tag)
        else:
            groups["walls"].append(tag)
    for name, tags in groups.items():
        assert tags, name
        gmsh.model.addPhysicalGroup(2, tags, name=name)
    gmsh.model.addPhysicalGroup(3, [volumes[0][1]], name="chamber")

    # finer cells on the tubes and the burner
    fine = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fine, "SurfacesList", groups["tubes"] + groups["burner"])
    ramp = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ramp, "InField", fine)
    gmsh.model.mesh.field.setNumber(ramp, "SizeMin", a.size / 3)
    gmsh.model.mesh.field.setNumber(ramp, "SizeMax", a.size)
    gmsh.model.mesh.field.setNumber(ramp, "DistMin", 0.02)
    gmsh.model.mesh.field.setNumber(ramp, "DistMax", 0.4)
    gmsh.model.mesh.field.setAsBackgroundMesh(ramp)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(a.out)
    n = len(gmsh.model.mesh.getElementsByType(4)[0])
    print(f"{a.out}: {n} tetrahedra; " + ", ".join(f"{k} {len(v)} faces" for k, v in groups.items()))
    gmsh.finalize()


if __name__ == "__main__":
    main()
