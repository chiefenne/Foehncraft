//! Thermal radiation by Monte Carlo transport of energy bundles (step 29,
//! R1a: grey gas, grey diffuse walls, symmetry planes and open boundaries).
//!
//! A call emits energy bundles from every cell (power `4κσT⁴V`) and every
//! grey wall face (power `εσT⁴A`), in proportion to their power, and
//! follows each bundle along the cell walk of `foehncraft_kernels::walk`
//! (D257). In each cell the bundle passes it deposits the fraction
//! `1 − exp(−κℓ)` of its weight, `ℓ` the exact length of its path in the
//! cell (track-length estimator); a grey wall absorbs the fraction ε and
//! reflects the rest in a cosine-weighted direction; a symmetry plane
//! reflects it specularly; an open boundary lets it go. Light bundles play
//! Russian roulette.
//!
//! Every weight is a whole number of quanta (the total emitted power over
//! 2⁵²), so the tallies are integer sums that do not depend on the order of
//! the additions, and the energy account closes exactly: emitted = absorbed
//! by the gas + absorbed by the walls + escaped + roulette (D107 (3), (4)).
//! Every random number is keyed by the global identity of the emitter, the
//! bundle and the event (D107 (2), D121), and positions and directions at
//! emission come from Owen-scrambled Sobol sequences (Burley, JCGT 9 (4),
//! 2020). A call therefore gives bitwise the same result on any number of
//! threads and, with the walk handed over at processor faces (D258), on any
//! number of ranks.
//!
//! The bundles move on the CPU threads, or on a GPU (D259,
//! [`Domain::on_gpu`] with the backend that `foehncraft_cuda::register`
//! installs): the same kernel body, integer tallies, so the GPU gives
//! bitwise the tallies of the CPU. [`CpuDevicePath`] runs the data path of
//! the GPU on the CPU threads.
//!
//! [`exact`] has the exact solutions of the verification cases (section 3
//! of `docs/requirements/step29-radiation.md`); [`output`] writes the fields
//! of a call as VTKHDF and CSV; [`Domain::trace`] records the paths of
//! chosen bundles of a call.

pub mod device;
pub mod exact;
pub mod output;
mod solve;
pub mod trace;

pub use device::{CpuDevicePath, DeviceTallies, GpuFactory, Transport, gpu_backend, register_gpu};
pub use foehncraft_kernels::radiation::{Trace, WallKind};
pub use solve::{Domain, Ledger, Points, SIGMA, Settings, Solution, Wall};
pub use trace::{Path, PathPoint, Source, Traced};
