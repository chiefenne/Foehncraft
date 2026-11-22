- **Operator exactness check of the mesh (D236).** Before the solve the
  velocity gradient and its gradient (the velocity Hessian of the
  curvature correction) are applied to linear and quadratic test fields on
  the mesh of the run. Cells where they are not exact are flagged:
  `MESH004` (a broken stencil), `MESH005` (with the curvature correction:
  it does not converge there), `MESH006` (with the QCR: the stress
  divergence does not converge there) and the limitation `LIM008`; the 50
  worst cells per operator go to `flow.json` and the `case.check` report,
  every cell to `operators.csv`. New key `solver.operator_tolerance`
  (0.01). No result changes. Gallery and recipe 189.
