# Case catalogue

The case files in `examples/cases/` with what each one is, the mesh it needs
and where it is discussed. Pictures are in the [Example gallery](gallery.md),
command sequences in the [Recipes](recipes.md).

How to read the **Mesh** column:

* **synthetic** – the mesh is generated in memory from the `[mesh]` section
  of the case file; no mesh file is needed.
* **generated** – the built-in test mesher makes the mesh from the
  `[mesh.generate]` section (D34).
* **Gmsh** – the mesh is read from a Gmsh MSH 4.1 file (D33); the `.geo`
  source is in `examples/meshes/`, the `.msh` is made from it with Gmsh
  and not stored in the repository.

The catalogue is split by physics; each page has its own table of contents.

| Page | What it covers |
|---|---|
| [Meshes (`examples/cases/mesh/`)](cases/meshes.md) | build or read a mesh, check, renumber and partition it; mesh keys |
| [Scalar cases (`examples/cases/scalar/`)](cases/scalar.md) | steady convection–diffusion of a scalar: manufactured solutions, convection schemes; scalar keys |
| [Poisson studies (`examples/cases/poisson/`)](cases/poisson.md) | the pressure solver on the Poisson equation, with hypre as the reference; Poisson keys |
| [Chemistry cases (`examples/cases/chemistry/`)](cases/chemistry.md) | the 0D reactor and the chemical equilibrium of a mechanism (SF2); reactor keys |
| [Flow cases (`examples/cases/flow/`)](cases/flow.md) | incompressible flow, laminar and with the RANS turbulence models, steady and transient, grouped by physics; what a flow run writes |
| [Flow case-file keys](cases/flow-keys.md) | every key of a flow case, the derived fields (D72), the boundary conditions (D29) and their sensible combinations (D30), expressions, materials and regions (D53, D102, D112), passive scalars, turbulence models (D129–D139) |
| [Python scripts (`examples/python/`, `examples/studies/`)](cases/python.md) | cases built, run and studied with the Python package |
| [Water and steam circuits (`examples/cases/circuits/`)](cases/circuits.md) | one-dimensional water and steam tube circuits of a boiler (SF7b): parallel tubes, a heated riser, an evaporator tube, a drum loop with natural circulation; circuit keys |
| [Structural cases (`examples/cases/structure/`)](cases/structure.md) | the linear structural solver (L53, MS5): beams, shells, quadratic solids from Gmsh, clamp springs, contact with a seat and a stop, transient and thermal, static and modal analyses; structural keys |
| [Tool inputs (`examples/cases/tools/`)](cases/tools.md) | JSON inputs of `foehn call`, one per tool of the registry (D196) |
| [Engineering gallery (`examples/applications/`)](cases/applications.md) | the industrial and laboratory application cases: stirred tanks, the shell and tube heat exchanger |
| [Benchmarks (`examples/cases/bench/`)](cases/benchmarks.md) | gate G1, the kernel benchmark; benchmark keys |

## Sections and keys

Links to the former single page land here; each row points to the new place.

| Section | Page |
|---|---|
| <a id="meshes-examplescasesmesh"></a>[Meshes (`examples/cases/mesh/`)](cases/meshes.md#meshes-examplescasesmesh) | [Meshes (`examples/cases/mesh/`)](cases/meshes.md) |
| <a id="mesh-case-file-keys"></a>[Mesh case-file keys](cases/meshes.md#mesh-case-file-keys) | [Meshes (`examples/cases/mesh/`)](cases/meshes.md) |
| <a id="scalar-cases-examplescasesscalar"></a>[Scalar cases (`examples/cases/scalar/`)](cases/scalar.md#scalar-cases-examplescasesscalar) | [Scalar cases (`examples/cases/scalar/`)](cases/scalar.md) |
| <a id="scalar-case-file-keys"></a>[Scalar case-file keys](cases/scalar.md#scalar-case-file-keys) | [Scalar cases (`examples/cases/scalar/`)](cases/scalar.md) |
| <a id="poisson-studies-examplescasespoisson"></a>[Poisson studies (`examples/cases/poisson/`)](cases/poisson.md#poisson-studies-examplescasespoisson) | [Poisson studies (`examples/cases/poisson/`)](cases/poisson.md) |
| <a id="poisson-case-file-keys"></a>[Poisson case-file keys](cases/poisson.md#poisson-case-file-keys) | [Poisson studies (`examples/cases/poisson/`)](cases/poisson.md) |
| <a id="flow-cases-examplescasesflow"></a>[Flow cases (`examples/cases/flow/`)](cases/flow.md#flow-cases-examplescasesflow) | [Flow cases (`examples/cases/flow/`)](cases/flow.md) |
| <a id="flow-case-file-keys"></a>[Flow case-file keys](cases/flow-keys.md#flow-case-file-keys) | [Flow case-file keys](cases/flow-keys.md) |
| <a id="derived-fields-d72"></a>[Derived fields (D72)](cases/flow-keys.md#derived-fields-d72) | [Flow case-file keys](cases/flow-keys.md) |
| <a id="boundary-conditions-d29"></a>[Boundary conditions (D29)](cases/flow-keys.md#boundary-conditions-d29) | [Flow case-file keys](cases/flow-keys.md) |
| <a id="sensible-combinations-d30"></a>[Sensible combinations (D30)](cases/flow-keys.md#sensible-combinations-d30) | [Flow case-file keys](cases/flow-keys.md) |
| <a id="expressions-materials-and-regions-d53-d102-d112"></a>[Expressions, materials and regions (D53, D102, D112)](cases/flow-keys.md#expressions-materials-and-regions-d53-d102-d112) | [Flow case-file keys](cases/flow-keys.md) |
| <a id="benchmarks-examplescasesbench"></a>[Benchmarks (`examples/cases/bench/`)](cases/benchmarks.md#benchmarks-examplescasesbench) | [Benchmarks (`examples/cases/bench/`)](cases/benchmarks.md) |
| <a id="case-file-keys"></a>[Case-file keys](cases/benchmarks.md#case-file-keys) | [Benchmarks (`examples/cases/bench/`)](cases/benchmarks.md) |

## Adding a feature

Every feature ships with a gallery example, case files and a recipe (D18,
as in PyAeroFlow); see [CONTRIBUTING.md](../CONTRIBUTING.md). Put the case file in the `examples/cases/` folder of its physics and add a row to that page, under the group that fits (a new page, and a row in the first table above, when a new physics begins). Add its keys to the page's key table; the complete list of keys is generated into [keys.md](keys.md).
