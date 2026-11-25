# Multi-rank MPI: spike, measurements and sized plan

Spike of 2026-11-22 for item 2 of the minimal path to v1.0 (D235, section 7a
of [requirements](../requirements.md)): multi-rank MPI on one machine with
multi-node designed in (O12), the reproducible mode across ranks (D31, D93)
and restart on another rank count. The spike measures and sizes; it builds
no feature. The code of the measurements is in three examples (section 3) and
the decision cards are answered in the pull request thread, not here.

> **2026-11-23 (maintainer decision): multi-rank MPI moved to v1.1.** W1
> (#200) and W2 (#219) are on `main` as groundwork; W3 onward, the scaling
> measurements (W9) and the GPU exact sums are v1.1. Where this note says
> "v1.0" for multi-rank work, read v1.1.

**Result in four sentences.** The parallel layer is a good start (communicator
trait with an MPI backend, METIS partitioning, local meshes with halos, a
tested halo exchange), but none of the numerics runs distributed yet: the
discretisation, the kernels, the linear solvers, the flow engine, output and
restart all assume one partition. Two ranks do **not** give the matrix of one
rank today, not even bit for bit on the geometry, because the local mesh is
re-derived with its own face orientation and face order. With the canonical
construction of section 3 (the fix is small) geometry and matrix rows are
bitwise equal to one rank on 1, 2, 3 and 4 ranks, hexahedra and polyhedra,
with halo values taken from the owner. The full path is about ten pull
requests; the two items with real risk are the distributed multigrid and the
flow-engine integration.

## 1. What exists

All facts are from the repository at `main` after #191.

| Part | State | Where |
| --- | --- | --- |
| Communicator trait, serial and MPI backend (`MpiComm`, feature `mpi`) | Works: rank, size, barrier, gather of one value, paired send and receive; `sum_reproducible` gathers in rank order and sums pairwise | `foehncraft-parallel/src/comm.rs` |
| Partitioning | METIS (one constraint, no weights) and recursive coordinate bisection; deterministic | `foehncraft-mesh/src/partition.rs` |
| Local meshes with halo | `extract`: own cells (inner, then outer, each in reverse Cuthill–McKee order), halo cells grouped by owner rank, one contiguous receive block per neighbour; faces inner, processor, boundary; every entity keeps its global id (`cell_id`, `face_id`, `point_id`); periodic meshes are refused | `foehncraft-mesh/src/decompose.rs` |
| Halo exchange | One message per neighbour carrying all fields, host-staged; GPU-aware transport refused; exchange of geometry (volume, centroid) | `decompose.rs`, example `mpi_halo` in CI with 2 and 4 ranks |
| Reproducible reductions across threads | Block sums plus a fixed tree: bitwise equal for any thread count; **not** equal for any partition | `foehncraft-core/src/repro.rs`, `foehncraft-parallel/src/repro.rs` |
| Discretisation, kernels, backends, multigrid, coupled solver, flow engine | One partition only. `Discretisation::new` stops with an error when the mesh has halo cells; no function of `foehncraft-fv` or `foehncraft-flow` takes a communicator | `foehncraft-fv/src/disc.rs` |
| Restart | The file stores `cell_id` and `face_id`, and reading maps by global id (so a different cell order is already handled); the text still promises bitwise continuation only for "same mesh, same decomposition and device" | `foehncraft-flow/src/restart.rs` |
| Output | VTKHDF writer takes one whole mesh | `foehncraft-mesh/src/vtkhdf.rs` |
| Multigrid | Setup on the host from the assembled matrix, greedy aggregation in local row order (so it gives the same aggregates for any thread count, not for any partition, as D93 (5) already notes); smoothers are pointwise (Jacobi), which D93 requires | `foehncraft-fv/src/amg/` |
| Sparse solves by hypre | Sequential build only | `foehncraft-fv/src/hypre.rs` |
| CI | MPI job with 2 and 4 ranks on Linux and macOS (smoke test, halo exchange) | `.github/workflows/ci.yml` |

Size of the parts that have to change (lines of Rust): mesh 16.5 k, fv 6.8 k,
flow 52.9 k, kernels 10.2 k, tools 4.1 k, CLI 1.9 k. The flow crate has about
98 sum and fold sites, 41 uses of the fixed-order sums and 35 calls of the
backend dot product and min-max reductions that must become rank-aware
reductions; 40 call sites of the backend operators (gradient, limiter, face
coefficients, assembly, matrix-vector product) need a halo exchange before
them or a defined owner.

## 2. What the reproducible mode asks for (D93 constraints met or open)

| D93 constraint | Today | Open |
| --- | --- | --- |
| (1) Global order: every per-cell sum over faces in global face order | Cell faces are in *local* face order (sorted by local owner); the local order differs from the global one | Yes, section 3 |
| (2) No result depends on where the cuts are | Not met (reductions, aggregation, geometry) | Yes |
| (3) Halo cells are copies of the owner's values; processor faces computed once from the same inputs | Halo volumes and centroids are exchanged; face area and face centroid are recomputed locally with a different node order | Yes, section 3 |
| (4) Reductions exact and order-independent | Fixed order per run, not per partition | Yes, section 4 |
| (5) Multigrid setup on the global graph | Local row order | Yes, section 5 |
| (6) Arithmetic: no contracted FMA, one in-house math library | Not examined here | Not part of the multi-rank work |
| (7) Random numbers keyed by global ids | Counter-based generator exists | Check in the audit |
| (8) Iteration control on exact global values | Needs (4) | Yes |
| Restart written on 4 ranks, continued on 2, bitwise the uninterrupted run | Mapping by global id exists; write and read paths are single process | Yes, W6 |

## 3. Feasibility test: do 2 ranks give the matrix of 1 rank?

**Question** (from the coordinator): can 2 ranks on one machine produce
identical matrices to 1 rank for one small case?

**Test.** Example `foehncraft-mesh/examples/mpi_matrix.rs`. Every rank builds
the same perturbed box mesh (12 × 10 × 8, polyhedra, or hexahedra with
`SPIKE_KIND=hex`), partitions it (METIS), extracts its part with two halo
layers and computes the geometry of the local mesh with the halo volumes and
centroids from the exchange. It then builds the rows of a diffusion matrix
(face coefficient |S|²/(d·S), the stage-2 form) for its own cells and compares
every number bit for bit with the rows of the one-rank run, keyed by global
ids. The one-rank run is the same program on the whole mesh.

**Result with the code as it is (`LocalMesh::geometry`).** Not equal, on any
rank count, including one rank (where the local mesh is the global mesh in
another order). Counts of own entities that differ in the last bits, polyhedral
case, 2 ranks (rank 0 of 643 own cells): face area vectors 1522 of 4350, face
centroids 1639, cell volumes 705, cell centroids 839, and **all 643
off-diagonal rows**, because the matrix coefficients inherit the geometry.

**Two causes.**

1. `extract` re-orients every face so that its owner is the lower local cell
   (reversing the node loop after the first node) and the polygon area and
   centroid are summed in that node order: the same face gets different last
   bits on the two sides and against the global mesh. About half of the
   interior faces are flipped.
2. A cell sums its faces (volume, centroid, and in a matrix row the
   diagonal) in local face order, not in global face order.

**Result with the canonical construction** (function `canonical_geometry` of
the example, about 60 lines): every face polygon is evaluated with the node
order and orientation of the global mesh (the local orientation only flips
the sign, exactly), every cell sums its faces in ascending global face id with
the global sense, halo cells receive volume and centroid from their owner by
the exchange. Then face areas, face centroids, cell volumes, cell centroids,
all off-diagonal entries and all diagonals summed in global face order are
**bitwise equal to the one-rank run on 1, 2, 3 and 4 ranks, polyhedra and
hexahedra** (8 runs, all pass). The diagonal summed in *local* face order
still differs on about 40 % of the rows, which is why D93 (1) says global
order. The example runs in the MPI job of CI and fails when any of this stops
being true.

**What it does not show.** It covers geometry and a pure diffusion matrix; the
convection coefficients, the non-orthogonal correction, the gradient (least
squares on the face neighbours) and the limiters need the same treatment and
are in W3. One side effect to expect: canonical geometry changes the one-rank
bits compared with v0.7 in the last place (it is a different, but fixed,
summation order), so any stored hash must be regenerated once; no reference
band is affected.

## 4. Cost of exact reductions (measured)

Example `foehncraft-parallel/examples/exact_sum_spike.rs`: a long accumulator
(the 2098-bit range of binary64 in 68 limbs of 32 bits in 64-bit words,
carries postponed, at most three limbs touched per addition) summed over 20
million values with a wide dynamic range and cancellation, on one core of the
Xeon at 2.3 GHz of this session.

- Partition independence: the exact sum has the same bits for 1, 2, 3, 4 and
  7 parts in any order of merging; the fixed-order pairwise sum of the parts,
  which is what `sum_reproducible` does today, gives different last bits from
  the one-part value for 2, 3, 4 and 7 parts (and from the exact value).
- Cost: 0.82 ns per value for the pairwise sum, 4.9 ns per value for the exact
  sum, a factor of about 6 for a sum over data in cache. A dot product streams
  two vectors (16 bytes per value) and is partly memory-bound, so the factor
  for a whole dot product is smaller; the share of dot products in a
  multigrid-preconditioned iteration is small. The end-to-end cost has to be
  measured in the pull request (D93 asks for it), not guessed.
- Message size: 68 limbs, 544 bytes per reduction, one all-reduce.

## 5. Cost of communication on one machine (measured)

Example `foehncraft-parallel/examples/mpi_latency.rs`, Open MPI 4.1.6 on 4
cores (shared memory), ranks oversubscribed to the cores:

| Ranks | Reproducible sum (gather of one value) | Ring, 100 values | 1 000 | 10 000 | 100 000 |
| --- | --- | --- | --- | --- | --- |
| 2 | 0.6 µs | 1.7 µs | 4.9 µs | 9.5 µs | 74 µs |
| 4 | 1.6 µs | 2.0 µs | 5.2 µs | 8.2 µs | 82 µs |

A halo message of 10 000 values costs about 10 µs and a reduction a few µs.
The share of communication in an iteration is not measured here (it needs the
distributed solver; W5 and W9 report it from the phase timers). What the
numbers say is that on one machine a message costs microseconds, so the work
is restructuring code, not tuning messages. Latencies on a cluster network
were not measured; that part stays "designed for" (O12), and the overlap of
the halo exchange with the inner rows (D38) is in the plan for it.

On one machine ranks do not add speed that threads do not already give: a
rank is a process with its own thread pool. The scaling report therefore has
to compare *ranks × threads* splits of the same cores (for 4 cores: 4×1, 2×2,
1×4), which is also the honest hybrid answer for a multi-node run.

## 6. Sized plan

Sizes are agent working time for code and tests, not elapsed time: CI waits
and long runs are separate lines; "not sized yet" means the spike did not give
a basis. One pull request each unless a line says two.

| ID | Work | Agent time | Run or CI | Risk |
| --- | --- | --- | --- | --- |
| W1 | **Canonical local mesh and geometry** (mesh): faces keep global node order and orientation plus a flip flag; cell faces in global face order; geometry per section 3; tests on 1 to 4 ranks, METIS and RCB, hex, tet, polyhedra, mixed; the spike is the prototype | 30–60 min | 1 CI cycle | low |
| W2 | **Exact reductions** (core, parallel, fv backends): accumulator, `Communicator::sum_exact` (all-reduce of the limbs), replace the dot product and min-max of the CPU backend; audit of the 98 + 41 + 35 sites; measured cost per dot and per iteration in the PR | 2–3.5 h | cost run of minutes | medium (the audit) |
| W3 | **Distributed discretisation and kernels**: `Discretisation` on a local mesh with halo (least-squares weights need layer 1, the gradient and limiter an exchange), processor faces computed from the same inputs on both sides, matrix slots in global face order, kernel ranges for inner, processor and boundary faces; GPU path not sized yet | 1.5–3 h | 1 CI cycle | medium |
| W4 | **Distributed Krylov and multigrid, replicated setup first**: matrix-vector product with the exchange overlapped with the inner rows; vectors own plus halo; the hierarchy is built on every rank from the gathered matrix in global order, which gives the aggregates of the one-rank run by construction (D93 (5)); per-level ownership and halo plans derived from the hierarchy; Jacobi smoothers stay; the reduction multigrid of the coupled momentum sweep (D179) the same way. Memory per rank is the whole matrix for the setup: fine for a few million cells on one machine, not a weak-scaling design. **The scalable setup** (distributed aggregation with priorities hashed from global ids, Bell, Dalton & Olson 2012; Galerkin products in global order) is not sized yet | 2–4 h | 2 CI cycles | high |
| W5 | **Flow engine on ranks** (two pull requests): simplec, coupled, transient, energy, scalars and turbulence with exchanges where fields are read across faces; boundary lists identical on all ranks; wall distance (the exact method needs all wall faces on every rank, or the PDE method); monitors, probes (owner rank), ledger and balances with exact sums; the rule "an observer restores the lagged state" kept | 2–4 h | CI cycles, 1 h of runs | high |
| W6 | **Restart and output**: write gathered in global order (files byte-identical, D84), read by global id (restart on another rank count); parallel output as a second step, not sized yet | 1–2 h | restart 4 → 2 test | medium |
| W7 | **Launch and agent interfaces**: communicator chosen by the run, only rank 0 writes log, journal and run record, rank and thread counts in the provenance (D85), a `ranks` input of the run tool and `foehn run --ranks N` (see card 4), errors with diagnostic codes for what does not run on ranks yet, build and install notes for the `mpi` feature (Open MPI and libclang) | 45–90 min | 1 CI cycle | low |
| W8 | **Verification in CI** (D93): the same cases on 1, 2, 4 ranks and 1, 2, 4 threads compared by hash; on a difference the first field, iteration and global cell id where the runs part; restart written on 4 and continued on 2 | 1–1.5 h | a few minutes per case | medium |
| W9 | **Scaling measurements and report**: strong and weak scaling as far as the hardware allows, ranks × threads splits, the overhead of the mode against the default mode on the cases of D45, T4 with two ranks sharing the GPU only if GPU ranks are in scope; "designed for" label for multi-node | 1 h of writing | 4–8 h of runs, spread over cases | low |

**Status.** W1 done in #200 (2026-11-22): the canonical construction is in the library (`Mesh::face_flip`, cell faces in global order, `Geometry::compute`, `decompose::canonical_check`); the one-rank bits did not change, because the reference is the mesh handed to the decomposition. W2 done in #219 (exact reductions). W3 done in #230 (2026-11-24): flipping a face (W1) kept the geometry exact but not the asymmetric face formulas of the solver (owner weight w against 1 − w), so `extract` now keeps every face in its global orientation, a processor face may be owned by the halo cell and is computed on both ranks from the same inputs; `Discretisation` and `CpuBackend::on_ranks` run the scalar solver on a local mesh (exchange of φ, gradient, ψ and the matrix-product vector, exact dots and extrema), bitwise one rank on 1 to 4 ranks for polyhedra, mixed cells and tetrahedra, linear upwind and bounded central, BiCGStab and CG with Jacobi; `ThreadComm` runs ranks as threads in tests. The matrix-vector product is already distributed, so W4 is the multigrid. W4 done in #231 (2026-11-24): every rank gathers the matrix in global order and builds the whole hierarchy (`Backend::amg_matrix`); `amg::dist::distribute` keeps the own rows per level (a coarse row belongs to the owner of the first fine row of its aggregate) and plans the per-level ghost exchange without communication; the V-cycle and the reduction multigrid (D179) run on the own rows, the coarsest level is gathered and solved on every rank; CG + multigrid and the reduction multigrid are bitwise one rank on 1 to 4 ranks. W5 part 1 done in #234 (2026-11-24): the CPU flow device runs on the local mesh (cell buffers own plus halo, the kind of every buffer agreed by all ranks, a halo exchange after every kernel that writes cells, reductions over own cells and own-owner faces summed exactly), the linear solves use the scalar backend on ranks, and the pressure reference, the mesh skewness and the start fields are global decisions (start fields computed on the global mesh, each rank keeps its cells); steady laminar SIMPLEC is bitwise one rank on 1 to 4 ranks for polyhedra, mixed cells and tetrahedra (threads and `mpi_flow` in CI) and `[run] ranks` of `foehn flow` checks a case. Kernels not yet safe on ranks (turbulence and wall distance, vertex and quadratic gradients, the MLP limiter, transpose products) refuse to run there. W5 part 2 done in #244 (2026-11-24): transient PISO, passive scalars, the energy equation of one fluid, k-ω SST and k-ε (wall functions included) and the coupled solver on ranks, bitwise one rank on 1 to 4 ranks in nine check cases (threads and `mpi_flow`); the turbulence cell kernels and wall functions exchange ν_t, F1 and the wall values, the wall distance is computed on the whole mesh and handed to the ranks (a part may hold no wall; replicated like the multigrid setup), the energy refresh reads halo temperatures, and the coupled solver's pseudo-time length and speed and the multigrid reuse decisions are exact global reductions. Not yet on ranks (refused): ζ-f and elliptic blending, curvature correction and QCR, MLP and quadratic gradients, gravity, turning zones, solid regions, schedules and level sweeps. Monitors, probes and balances need the run loop on ranks and move to W7. W7 part 1 done in #279 (2026-11-25): the run itself on ranks, `[run] rank_count = N`, `foehn flow --ranks N` (card 4 named it `foehn run`, which replays scripts in this CLI), the input `ranks` of `flow.run`, threads of one process or the processes of `mpirun` (a build with the `mpi` feature; only rank 0 prints and journals). Every rank runs the loop of the case; rank 0 reads control files, signals and limits and hands its decision to the others; only rank 0 writes the log, `run.jsonl` and the output. Log columns, probes, boundary summaries, wall values, the ledger, `fields.vtkhdf` and `restart.h5` are computed on the fields gathered into the global order (the W6 gather, without matrices) with the discretisation of the whole mesh, so every output file is byte for byte one rank (2, 3 and 4 ranks in the test, whole `result.json` compared). Steady laminar runs so far; passive scalars, energy and turbulence models (their end summaries read device data of the rank, partition-local wall lists and floor counters) and transient runs are W7 part 2; settings through the control port are refused on ranks. First numbers on 1M cells (example `ranks_scaling`, a steady laminar channel on 200 × 100 × 50 hexahedra, 10 SIMPLEC iterations, setup excluded, this session's 4-core machine, ranks as threads of one process sharing the 4 cores): no ranks 0.84 s per iteration; 1 rank × 4 threads 0.86 s (+2 %); 2 ranks × 2 threads 0.93 s (+10 %); 4 ranks × 1 thread 1.08 s (+28 %); every run bitwise the reference; peak memory 10.4 GB for the whole mesh and all ranks in one process. With the cores fixed this measures the price of the decomposition (halo exchanges, exact reductions over ranks, the replicated multigrid setup at every rebuild), not a speed-up; strong scaling over more cores and nodes is W9. W6 done in #262 (2026-11-24): `Restart::map_fields` drives every array and matrix of a restart through one map (reorder by identity, gather, scatter); `ranks::Layout` gathers the ranks' parts into the order of the global mesh, each face from the rank of its owner cell and each multigrid matrix row by row into the sliced-ELL storage of the global mesh, and the file written from it has the bytes of the one-rank restart (D84); a file in any order is read on any rank count by global identity. Nine check cases written on 4 ranks continue on 2 (and written on 2 continue on 3) bitwise the uninterrupted one-rank run; `[run] ranks` accepts `[run] restart` and checks the gathered file byte for byte. The file is gathered on every rank and written by rank 0 (parallel HDF5 not needed at this size); restarts of conjugate cases and moving meshes stay refused on ranks. The field output of a run on ranks moves with the run itself to W7. W7 part 1 done in #279 (steady laminar runs on ranks, rank 0 writes every output file byte for byte as one rank); W7 part 2 done in #290 (2026-11-25): passive scalars, the energy equation, k-ω SST and k-ε, the rotating frame and transient runs, with the scalar, heat and turbulence summaries gathered into the order of the whole mesh; species, statistics and the conjugate time-saving loops stay refused on ranks. #303 (2026-11-25): species transport, the variable density of stage 7c and `[output.statistics]` on ranks, every output file byte for byte one rank; the MLP vertex limiter runs on ranks (each rank bounds the vertices of its own cells, the ranks sharing a vertex of a processor face combine minimum and maximum), its cell-volume ε scaling sums the domain volume over the ranks; the species, density and history summaries are computed over all ranks. Still refused: a closed volume, gravity, solid and fluid regions (and with them the time-saving loops of conjugate runs and several fluids), turning zones.
W2 done in #219 (2026-11-23): `foehncraft_core::exact` (long accumulator, correctly rounded, a batched path with one integer bin per exponent), `Communicator::sum_exact` and `all_sum_i64` (one all-reduce of 71 integers), exact dot product and order-independent min/max in the CPU backend and the flow CPU device, about 30 sums over cells, faces and patches of the flow crate exact, the example `mpi_reductions` in CI on 1 to 4 ranks. One-rank results change once in the last bits against v0.7 (no reference band affected). Measured cost on this session's Xeon: see the table in the pull request. The GPU keeps block sums until v1.1 (owner decision of 2026-11-23). Left for W5: the hand-written accumulation loops of the flow crate (`s += …` over cells) and the per-rank all-reduce of each global sum; time-series statistics, timings and three-component products need no change.

**Total of W1 to W8:** about 11 to 20 agent hours of active work in nine to
ten pull requests, plus CI waits, plus 4 to 8 hours of runs for W9. The
estimate has its weight in W4 and W5; if W4 needs the scalable setup, add what
section 6 calls not sized yet. Order: W1 → W2 → W3 → W4 → W5 → W6 → W7 → W8 →
W9; W2 and W3 can run in parallel after W1 if two threads are used.

## 7. Flaws of industry behaviour to avoid (for the challenger)

- Partition-dependent smoothers (per-partition Gauss–Seidel, ILU) and sums in
  partition order make a result depend on the rank count; users cannot tell
  a physical difference from a decomposition difference. Foehncraft keeps
  pointwise smoothers and exact reductions (D93), and pays the measured cost.
- Halo values recomputed by the receiving rank instead of copied from the
  owner are a second source of last-bit differences; D93 (3) forbids it.
- The user must not have to choose a partition to get a result: results do
  not depend on the cuts, and the partition is a performance setting.

## 8. Decisions of the maintainer

**Card 1 (2026-11-22): scope.** The maintainer chose "everything in v1.0":
besides fluid-only flows, periodic faces, sliding interfaces, moving
reference frames (MRF) and CHT with several regions also run on ranks in
v1.0. Consequences for the plan, none of which is sized by this spike:

| ID | Work | Basis | Agent time |
| --- | --- | --- | --- |
| W10 | Periodic faces in the decomposition (`decompose` refuses periodic meshes today; the affine map of D109 across a periodic pair that crosses a rank cut; the driving force of the periodic pipe as an exact global reduction) | Code read: one refusal in `decompose.rs`, the rest unknown | not sized yet |
| W11 | MRF zones (a source term in a cell zone, probably partition-local) | Not read in this spike | not sized yet |
| W12 | Sliding interfaces: imprinting and the supermesh of D58 and the rotating zone across rank cuts, torque and flux sums exact | Not read in this spike | not sized yet |
| W13 | CHT views per rank with halo plans and the two-constraint METIS partition (D159, D122), interface faces on rank cuts computed once | Plan exists in the stage 7 page, section 11 | not sized yet |

These four were deferred by the "fluid-only" option to keep the path short;
taking them into v1.0 moves the tag. A first estimate needs a code read of
each, about one hour of agent time per item, which is the proposal for the
next step if the maintainer wants numbers before the order is fixed.

**Card 2: multigrid setup.** Replicated first (W4 as written); the scalable
setup replaces it later without changing any result, because both give the
aggregates of the one-rank run.

**Card 3: GPU.** CPU ranks only in v1.0. GPU runs stay on one rank; GPU-aware
MPI, exact reductions on the GPU and the CPU against GPU hash check with two
ranks sharing the T4 move to v1.1. W9 therefore measures CPU ranks only.

**Card 4: launch.** `foehn run --ranks N` (the `ranks` input of the run tool,
so the CLI, the Python package and the MCP server share it); the run starts
itself under the system MPI launcher. Starting with `mpirun` by hand keeps
working. This is W7.

**Nothing is open from the spike.** What remains is the order of work and the
go to start W1; coding starts only on the maintainer's explicit go.

## References

- A. Bell, S. Dalton, L. N. Olson, *Exposing fine-grained parallelism in
  algebraic multigrid methods*, SIAM J. Sci. Comput. 34, C123–C152 (2012).
- J. Demmel, H. D. Nguyen, *Parallel reproducible summation*, IEEE Trans.
  Comput. 64, 2060–2070 (2015).
- S. Collange, D. Defour, S. Graillat, R. Iakymchuk, *Numerical
  reproducibility for the parallel reduction on multi- and many-core
  architectures*, Parallel Comput. 49, 83–97 (2015).
- U. Kulisch, *An axiomatic approach to rounded computations*, Numer. Math.
  18, 1–17 (1971).
