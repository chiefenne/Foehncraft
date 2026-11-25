# Thermal radiation

Part of the [Recipes](../recipes.md): thermal radiation by Monte Carlo energy bundles (step 29, L15). In R1a the radiation is called from Rust on a prescribed temperature field. The `[radiation]` keys of the case file and the coupling to the energy equation come with R1b.

| # | Recipe |
|---|---|
| 259 | [Compute the radiation of a prescribed temperature field, and check it](#radiation-call) |

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

On MPI ranks (Open MPI or MPICH, with the `mpi` feature):

```
cargo build --release -p foehncraft-radiation --features mpi --example mpi_radiation
RAYON_NUM_THREADS=1 mpirun -n 4 target/release/examples/mpi_radiation
```

[↑ back to top](#thermal-radiation)
