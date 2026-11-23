//! Backends: where the kernel bodies run.
//!
//! A backend owns the device-resident copy of the [`Discretisation`] and a
//! fixed set of named vectors ([`V`]); the solver drives it through the
//! operations of [`Backend`] and never touches the data directly. The CPU
//! backend runs the bodies of `foehncraft_kernels` in (multi-threaded) loops,
//! the CUDA backend (crate `foehncraft-cuda`, D60) launches the same bodies
//! compiled to PTX (`gpu/ptx-kernels`). GPU backends live in their own
//! crates and register with [`register_gpu`]: this crate has no vendor
//! code. Every operation writes each output from one
//! index only, and dot products use fixed blocks, so the results do not
//! depend on the number of threads, and CPU and GPU agree (D31).

pub mod cpu;

use crate::amg::Hierarchy;
use crate::amg::cycle::LevelOps;
use crate::disc::Discretisation;
use std::sync::OnceLock;

/// The named vectors of a backend.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(missing_docs)]
pub enum V {
    /// Cell values φ (also the solution of the linear solve).
    Phi,
    /// Gradient components.
    Gx,
    Gy,
    Gz,
    /// Limiter value per cell.
    Psi,
    /// Face coefficients c_P, c_N, b_f.
    Cp,
    Cn,
    Bf,
    /// Right-hand side and inverse diagonal.
    Rhs,
    Dinv,
    /// BiCGStab work vectors.
    R,
    Rh,
    P,
    Ph,
    Vv,
    S,
    Sh,
    T,
    /// CG work vectors: preconditioned residual, A p, previous residual.
    Z,
    Q,
    Rold,
}

impl V {
    /// Vectors with one value per cell.
    pub const CELL: [V; 18] = [
        V::Phi,
        V::Gx,
        V::Gy,
        V::Gz,
        V::Psi,
        V::Rhs,
        V::Dinv,
        V::R,
        V::Rh,
        V::P,
        V::Ph,
        V::Vv,
        V::S,
        V::Sh,
        V::T,
        V::Z,
        V::Q,
        V::Rold,
    ];
    /// Vectors with one value per face.
    pub const FACE: [V; 3] = [V::Cp, V::Cn, V::Bf];

    /// Index into a backend's vector table.
    pub fn index(self) -> usize {
        self as usize
    }

    /// Number of vectors.
    pub const COUNT: usize = 21;
}

/// Gradient limiter (D25).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Limiter {
    /// No limiter (ψ = 1).
    None,
    /// Venkatakrishnan (1995), ε² = (K h)³.
    Venkatakrishnan {
        /// K.
        k: f64,
    },
    /// Venkatakrishnan–Wang (Wang 2000), ε² = (K (φ_max − φ_min))².
    VenkatakrishnanWang {
        /// K.
        k: f64,
    },
    /// Barth–Jespersen (1989): ψ = min(1, Δ₊/Δ₋).
    BarthJespersen,
}

/// Settings of the face coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceSettings {
    /// 0 first-order upwind, 1 linear upwind.
    pub blend: f64,
    /// Fraction of the bounded central scheme (0: none).
    pub central: f64,
    /// β of the bounding of the central scheme (0: none).
    pub beta: f64,
    /// Apply the non-orthogonal correction.
    pub nonorth: bool,
    /// Limit of the non-orthogonal correction (negative: none).
    pub nonorth_limit: f64,
}

/// Operations of a backend; the multigrid level operations of
/// [`LevelOps`] drive the V-cycle of [`crate::amg::cycle::vcycle`].
pub trait Backend: LevelOps {
    /// Name (`cpu`, `cuda`).
    fn name(&self) -> String;
    /// Device description.
    fn device(&self) -> String;
    /// Number of cells.
    fn n_cells(&self) -> usize;
    /// Copies host values into a vector.
    fn upload(&mut self, v: V, data: &[f64]) -> anyhow::Result<()>;
    /// Copies a vector to the host.
    fn download(&mut self, v: V) -> anyhow::Result<Vec<f64>>;
    /// Waits for all queued work (for phase timers).
    fn sync(&mut self) -> anyhow::Result<()>;

    /// Least-squares gradient of `Phi` into `Gx`, `Gy`, `Gz`.
    fn gradient(&mut self) -> anyhow::Result<()>;
    /// Limiter of `Phi` and its gradient into `Psi`.
    fn limiter(&mut self, l: Limiter) -> anyhow::Result<()>;
    /// Face coefficients into `Cp`, `Cn`, `Bf`.
    fn faces(&mut self, s: FaceSettings) -> anyhow::Result<()>;
    /// Matrix rows, `Rhs` and `Dinv`, with under-relaxation `alpha`.
    fn assemble(&mut self, alpha: f64) -> anyhow::Result<()>;

    /// out = A x.
    fn spmv(&mut self, x: V, out: V) -> anyhow::Result<()>;
    /// out = b − A x.
    fn residual(&mut self, x: V, b: V, out: V) -> anyhow::Result<()>;
    /// (a, b): exact on the CPU (MPI step W2), fixed block order on the GPU.
    fn dot(&mut self, a: V, b: V) -> anyhow::Result<f64>;
    /// Minimum and maximum of a cell vector.
    fn min_max(&mut self, a: V) -> anyhow::Result<(f64, f64)>;
    /// P = R + β (P − ω Vv).
    fn update_p(&mut self, beta: f64, omega: f64) -> anyhow::Result<()>;
    /// out = x − α z.
    fn sub_scaled(&mut self, x: V, z: V, alpha: f64, out: V) -> anyhow::Result<()>;
    /// Phi += α Ph + ω Sh.
    fn update_x(&mut self, alpha: f64, omega: f64) -> anyhow::Result<()>;
    /// out = Dinv x.
    fn precondition(&mut self, x: V, out: V) -> anyhow::Result<()>;
    /// dst = src.
    fn copy(&mut self, src: V, dst: V) -> anyhow::Result<()>;
    /// y = y + α x.
    fn axpy(&mut self, alpha: f64, x: V, y: V) -> anyhow::Result<()>;
    /// y = x + β y.
    fn xpby(&mut self, x: V, beta: f64, y: V) -> anyhow::Result<()>;

    /// The assembled matrix values in the sliced-ELL slots of the
    /// discretisation (for the multigrid setup on the host).
    fn matrix_values(&mut self) -> anyhow::Result<Vec<f64>>;
    /// Replaces the matrix by values in the sliced-ELL slots of the
    /// discretisation (assembled elsewhere, e.g. by the flow solver) and sets
    /// `Dinv` from its diagonal (slot 0).
    fn set_matrix(&mut self, val: &[f64]) -> anyhow::Result<()>;
    /// Stores a multigrid hierarchy on the device, in its precision, in the
    /// active slot (see [`amg_slot`](Self::amg_slot)).
    fn amg_load(&mut self, h: &Hierarchy) -> anyhow::Result<()>;
    /// Makes hierarchy slot `slot` the active one: the V-cycle and
    /// [`amg_load`](Self::amg_load) use it, the others are kept (slot 0:
    /// the pressure; stage 6 PR 6d adds slot 1 for the turbulence
    /// equations, so that their hierarchy does not replace the pressure's).
    fn amg_slot(&mut self, slot: usize) -> anyhow::Result<()> {
        anyhow::ensure!(
            slot == 0,
            "this backend holds one multigrid hierarchy only (slot {slot} asked)"
        );
        Ok(())
    }
    /// The host copy of the discretisation.
    fn discretisation(&self) -> Option<&Discretisation>;
}

/// Builds a GPU backend on device `ordinal` from the discretisation and the
/// kernel module (for CUDA the PTX file of `gpu/ptx-kernels/build.py`).
pub type GpuFactory =
    fn(&Discretisation, usize, &std::path::Path) -> anyhow::Result<Box<dyn Backend>>;

static GPU: OnceLock<(&'static str, GpuFactory)> = OnceLock::new();

/// Registers the GPU backend of this build under `name` (`cuda`); called
/// once by the front end (`foehncraft_cuda::register`). Later calls are
/// ignored.
pub fn register_gpu(name: &'static str, f: GpuFactory) {
    let _ = GPU.set((name, f));
}

/// The name of the registered GPU backend, if any.
pub fn gpu_backend() -> Option<&'static str> {
    GPU.get().map(|g| g.0)
}

/// Builds a backend by name: `cpu` (all threads), `cpu-serial`, or `cuda`
/// (when registered; `ptx` is the module built by `gpu/ptx-kernels/build.py`).
pub fn create(
    name: &str,
    disc: &Discretisation,
    threads: Option<usize>,
    ptx: &std::path::Path,
) -> anyhow::Result<Box<dyn Backend>> {
    match name {
        "cpu" => Ok(Box::new(cpu::CpuBackend::new(disc.clone(), threads)?)),
        "cpu-serial" => Ok(Box::new(cpu::CpuBackend::new(disc.clone(), Some(1))?)),
        "cuda" => match GPU.get() {
            Some((n, f)) if *n == name => f(disc, 0, ptx),
            _ => anyhow::bail!("this build has no GPU backend; build with `--features cuda`"),
        },
        other => anyhow::bail!("unknown backend `{other}` (cpu, cpu-serial, cuda)"),
    }
}
