- **Engineering gallery: shell and tube heat exchanger mesh (application
  203).** `examples/applications/shell_tube_hx/mesh.py` turns a STEP
  assembly (housing, seven tubes, two tubesheets, three segmental baffles)
  into one conformal multi-region tet mesh with Gmsh: hot air on the shell
  side, water on the tube side and four metal zones, the four nozzles
  closed as named inlet and outlet faces, every interface shared node for
  node (983 249 cells at level 2). A new page of the engineering gallery,
  heat exchangers, holds the entry with a section and a cut-away of the
  mesh; recipe 203 gives the steps for one's own assembly. The flow run
  waits for two fluids with their own materials in one conjugate case
  (ENE006).
