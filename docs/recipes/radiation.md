# Thermal radiation

Part of the [Recipes](../recipes.md): thermal radiation by Monte Carlo energy bundles (step 29, L15). In R1a the radiation is called from Rust on a prescribed temperature field. The `[radiation]` keys of the case file and the coupling to the energy equation come with R1b.

| # | Recipe |
|---|---|
| 259 | [Compute the radiation of a prescribed temperature field, and check it](#radiation-call) |
| 281 | [Follow chosen bundles of a radiation call](#radiation-paths) |

---

<a id="radiation-call"></a>
## 259. Compute the radiation of a prescribed temperature field, and check it

One call of `foehncraft_radiation` takes a mesh, the temperature (K) and
the absorption coefficient (1/m) per cell, and one wall condition per
boundary of the mesh. The conditions are `Wall::black(T)`,
`Wall::grey(ε, T)`, `Wall::mirror()` for a symmetry plane and
`Wall::open()` for an opening. The call returns the net source per cell
(W/m³) with its statistical error, the incident and net flux per wall face
(W/m²), and the energy account:

```rust
use foehncraft_radiation::{Domain, Settings, Wall};

let walls: Vec<Wall> = mesh.boundaries.iter().map(|b| match b.name.as_str() {
    "outlet" => Wall::open(),
    "tubes" => Wall::grey(0.85, 620.0),
    _ => Wall::grey(0.75, 1150.0),
}).collect();
let s = Settings { bundles: 4_000_000, ..Settings::default() };
let sol = Domain::whole(&mesh, &geo).solve(&temperature, &kappa, &walls, &s)?;
assert!(sol.ledger.closes());
foehncraft_radiation::output::write_vtkhdf("radiation.vtkhdf", &mesh, &sol, &[])?;
```

- **Noise.** The error falls with the number of bundles, faster than
  1/√N thanks to the scrambled Sobol points (gallery 259). `source_error`
  and `wall_net_error` give it per cell and face; `Settings::batches`
  (8 by default) sets how many batches it is estimated from.
- **Repeatable.** The same `Settings::call` gives bitwise the same result
  on any number of threads; another `call` number gives an independent
  estimate.
- **Ranks.** On a rank, use `Domain::rank(&local, &geo, &comm)` with the
  local mesh of `foehncraft_mesh::decompose::extract`. Bundles cross to
  the rank behind a processor face, and the tallies are bitwise those of
  one rank.
- **GPU.** With the CUDA backend registered
  (`foehncraft_cuda::register()`), use `Domain::whole(...).on_gpu(0, ptx)`.
  The CPU = GPU check is `python3 -m tools.bench.radiation_cloud` on a
  machine with an NVIDIA GPU, or the notebook
  `tools/notebooks/radiation_colab.ipynb` on a free Colab T4 (Runtime,
  Change runtime type, T4 GPU, then Run all).

Check the build against the exact solutions and look at the furnace of
gallery 259:

```
cargo test --release -p foehncraft-radiation --test cases -- --nocapture
pip install gmsh
python3 examples/meshes/furnace.py --size 0.15 --out output/furnace.msh
cargo run --release -p foehncraft-radiation --example furnace -- \
    output/furnace.msh output/furnace 4000000
```

The gallery figures of gallery 259 and 281 use KoflCraft's mesh of the same
furnace, with wall layers on every wall and the side walls' layers yielding
to the tube layers (`layers_yield`): 2.11 million cells. The example takes
it as the `foehncraft.vtu` hand-off of `kofl generate-duct` (polyhedral mode):

```
cargo run --release -p foehncraft-radiation --example furnace -- \
    PATH/TO/KOFLCRAFT/OUTPUT/foehncraft.vtu output/furnace_fine 16000000 60
python3 -m tools.plot.radiation_furnace output/figures_fine output/furnace_fine
python3 -m tools.plot.radiation_paths output/figures_fine output/furnace_fine
```

The fine call of 18 million bundles takes 2 h 15 min on 16 threads (8 cores
with Hyper-Threading); the coarse call takes about 5 minutes. The plot
scripts average the Monte Carlo fluxes over about 10 cm of wall and over a
volume band of the chamber: single faces and cells carry the noise of the
bundles they catch, so they are drawn as averages, not one by one.

On MPI ranks (Open MPI or MPICH, with the `mpi` feature):

```
cargo build --release -p foehncraft-radiation --features mpi --example mpi_radiation
RAYON_NUM_THREADS=1 mpirun -n 4 target/release/examples/mpi_radiation
```

[↑ back to top](#thermal-radiation)

---

<a id="radiation-paths"></a>
## 281. Follow chosen bundles of a radiation call

`Domain::trace` follows chosen bundles of the call that `Domain::solve`
makes with the same arguments: the same emitters, weights and Sobol
points, moved by the same kernel. `pick` says how many bundles of each
emitter to follow, counted from its first (zero: none; more than the
emitter has: all). Every path lists the points where the bundle leaves a
cell, is reflected at a wall or ends, with the quanta it still carries,
and the quanta it leaves in each cell and on each wall:

```rust
use foehncraft_radiation::{Domain, Source, Trace};

let domain = Domain::whole(&mesh, &geo);
let sol = domain.solve(&temperature, &kappa, &walls, &s)?;
// the first two bundles of cell 120 and the first bundle of boundary face 7
let traced = domain.trace(&temperature, &kappa, &walls, &s, |e| match e {
    Source::Cell(120) => 2,
    Source::Wall(7) => 1,
    _ => 0,
})?;
for path in &traced.paths {
    for p in &path.points {
        let fraction = p.w as f64 / path.w0 as f64; // energy left
        if p.at == Trace::Escaped { /* left through an opening at p.x */ }
    }
    let watts = path.w0 as f64 * traced.quantum;
}
```

- **Exact share.** The quanta a traced bundle leaves are its share of the
  call's tallies: tracing every bundle (`|_| u32::MAX`) gives the call's
  absorbed power per cell, net flux per face and energy account again, to
  the quantum.
- **Any bundle count.** Bundle k of an emitter starts from the same Sobol
  point whatever `Settings::bundles` is; only its weight changes. The 60
  paths of the furnace (gallery 281) are the same in a call of 200 000
  and of 16 million bundles, so a small call traces the paths of a large
  one.
- **A fair sample.** To draw paths that represent the call, pick emitters
  in proportion to their power, as the example `furnace` does with its
  last argument (`PATHS`, written to `paths.csv`).
- **Whole mesh only.** Tracing runs on the threads of one process; on a
  rank's part of a mesh it returns an error, since a path would stop at
  the processor faces.

```
python3 examples/meshes/furnace.py --size 0.15 --out output/furnace.msh
cargo run --release -p foehncraft-radiation --example furnace -- \
    output/furnace.msh output/furnace 200000 60
python3 -m tools.plot.radiation_paths output/figures output/furnace
```

[↑ back to top](#thermal-radiation)
