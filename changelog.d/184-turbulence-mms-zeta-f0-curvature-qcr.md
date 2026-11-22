- **Manufactured solutions of k-ε-ζ-f₀, the curvature correction and the
  QCR (stage 6 specification 9.2).** The harness of PR #148 now covers
  ζ-f₀ (orders 1.85–2.27 on hexahedra, polyhedra, tetrahedra), the
  curvature correction of SST in a second prescribed flow whose strain axes
  turn (1.9–2.1 on hexahedra and polyhedra) and the explicit stress of the
  QCR2000/QCR2013 (order 2 on regular hexahedra and polyhedra). Findings,
  not fixed here: the curvature correction does not converge on tetrahedra
  (DS_ij/Dt from the least-squares gradient of the gradient), and the
  divergence of the explicit QCR stress is first order on tetrahedra and
  does not converge on randomly perturbed meshes. New:
  `Flow::explicit_stress_force`. The grades of the three features list
  the gaps. Gallery and recipe 184.
