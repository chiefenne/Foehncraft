- **Curvature correction consistent on irregular meshes (D236 step 3).**
  DS_ij/Dt of the curvature correction of SST takes S_ij from a velocity
  gradient that is second order everywhere (the quadratic least-squares
  fit of the QCR on the cells where the face-neighbour gradient is first
  order), and on the cells where the least-squares gradient of S_ij is not
  consistent, the gradient of S_ij from a quadratic least-squares fit over
  interior cells. The manufactured solution converges at order 2 on
  regular and 15 % perturbed tetrahedra (error of k 0.54 %, 0.14 %,
  0.036 %) instead of staying at about 1.3 %; a channel in pure shear stays
  unchanged to round-off. The operator check reports
  `quadratic_hessian`; `MESH005` remains only where the fit fails, and a
  run whose cells take the fit writes `operators.csv` and the check in
  `flow.json`. Gallery and recipe 197.
