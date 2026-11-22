// T-junction mixer of D45 for the grid-convergence study: the geometry of
// tjunction.geo (main pipe diameter 0.1 m, length 1 m; branch diameter
// 0.05 m, length 0.3 m, joining at x = 0.3 m) with a global size h that the
// study sets, tetrahedra as Gmsh makes them (Netgen-optimised tetrahedra
// break the conversion into polyhedra: GEN001 at boundary cells).
//
//   gmsh -3 examples/meshes/tjunction_gci.geo -setnumber h 0.02 -format msh41 -o tj.msh
//
// Boundaries: inlet_main, inlet_branch, outlet, wall; the volume is the
// zone fluid. The size grows from h/2 at the junction edge to h elsewhere.

DefineConstant[ h = {0.02, Min 0.002, Max 0.1, Name "h"} ];
SetFactory("OpenCASCADE");
Mesh.MeshSizeMax = h;
Mesh.MeshSizeMin = h / 4;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MshFileVersion = 4.1;
Mesh.ElementOrder = 1;

Cylinder(1) = {0, 0, 0, 1.0, 0, 0, 0.05};
Cylinder(2) = {0.3, 0, 0, 0, 0.3, 0, 0.025};
BooleanUnion(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

eps = 1e-6;
in_main()   = Surface In BoundingBox{-eps, -1, -1, eps, 1, 1};
outlet()    = Surface In BoundingBox{1 - eps, -1, -1, 1 + eps, 1, 1};
in_branch() = Surface In BoundingBox{0, 0.3 - eps, -1, 1, 0.3 + eps, 1};
all()       = Boundary{ Volume{3}; };
wall()      = all();
wall() -= in_main();
wall() -= outlet();
wall() -= in_branch();

Physical Surface("inlet_main") = in_main();
Physical Surface("inlet_branch") = in_branch();
Physical Surface("outlet") = outlet();
Physical Surface("wall") = wall();
Physical Volume("fluid") = {3};
