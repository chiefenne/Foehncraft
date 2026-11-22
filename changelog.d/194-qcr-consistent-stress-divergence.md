- **QCR stress divergence consistent on irregular meshes (D236 step 2).**
  On the cells where the operator check finds the face-neighbour velocity
  gradient first order, the explicit QCR stress takes a quadratic
  least-squares velocity gradient over the neighbours of the neighbours,
  and the face value of the stress is corrected for skewness with cell
  gradients of the stress. The stress force of the manufactured solution
  now converges on meshes with 15 % perturbed vertices (hexahedra 4.6 %,
  1.5 %, 0.65 %; tetrahedra first order) instead of staying at about 5 %.
  `MESH006` is now reported only for cells where the quadratic fit is rank
  deficient. QCR cases on meshes with flagged cells change slightly (the
  square duct, regular hexahedra with 286 flagged wall cells: largest secondary velocity 1.558 % of U_b, before 1.558 %). Gallery and recipe 194.
