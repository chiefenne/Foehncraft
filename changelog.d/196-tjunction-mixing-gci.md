- **T-junction mixing on polyhedra with grid convergence (D45).** New case
  `tjunction_mixing.toml`: laminar T-junction (Re 200 and 100) with a passive
  scalar, Gmsh tetrahedra converted into polyhedra, three meshes of 4 852,
  16 665 and 44 738 cells (`examples/studies/tjunction_gci.py`). The main-inlet
  pressure converges at order 1.77 (GCI 2.4 %); the branch-inlet pressure does
  not (order 0.2, GCI 29 %), and the conversion fails with GEN001 on some
  meshes; both are written as gaps. Validation against measurements is
  pending. Gallery and recipe 196.
