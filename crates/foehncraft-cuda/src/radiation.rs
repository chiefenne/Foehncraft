//! Radiation (step 29, R1a, D259): the transport of the energy bundles on
//! an NVIDIA GPU, the entry `radiation_transport` of `gpu/ptx-kernels`, one
//! GPU thread per bundle. The walk mesh is uploaded once; each launch
//! uploads the medium and the bundle records, zeroes the integer tallies,
//! runs the kernel body that the CPU runs, and reads the tallies back.
//! The tallies are whole numbers of quanta, so the GPU gives bitwise the
//! tallies of the CPU threads.

use anyhow::{Context, Result, ensure};
use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, DevicePtr, DeviceRepr, LaunchConfig,
    PushKernelArg, ValidAsZeroBits,
};
use cudarc::nvrtc::Ptx;
use foehncraft_kernels::entry::Raw;
use foehncraft_kernels::radiation::{BUNDLE_WORDS, Medium, TransportArgs};
use foehncraft_kernels::walk::WalkMesh;
use foehncraft_radiation::{DeviceTallies, Transport};
use std::path::Path;
use std::sync::Arc;

/// Threads per block (the walk is long and branchy: small blocks).
const THREADS: u32 = 128;

/// The argument block as the kernel parameter (passed by value).
#[repr(transparent)]
struct Arg(TransportArgs);
// SAFETY: TransportArgs is a #[repr(C)] struct of integers and
// floating-point numbers; the PTX entry takes it by value with the same
// layout (one definition in foehncraft-kernels).
unsafe impl DeviceRepr for Arg {}

/// A device copy of a host array and its address for the argument block.
struct Up<T> {
    _buf: CudaSlice<T>,
    raw: Raw<T>,
}

fn up<T: DeviceRepr + ValidAsZeroBits>(s: &Arc<CudaStream>, v: &[T]) -> Result<Up<T>> {
    let buf = if v.is_empty() {
        s.alloc_zeros::<T>(1)?
    } else {
        s.clone_htod(v)?
    };
    let raw = Raw::new(buf.device_ptr(s).0, v.len());
    Ok(Up { _buf: buf, raw })
}

/// Zeroed device array of `n` entries.
fn zeros(s: &Arc<CudaStream>, n: usize) -> Result<Up<u64>> {
    let buf = s.alloc_zeros::<u64>(n.max(1))?;
    let raw = Raw::new(buf.device_ptr(s).0, n);
    Ok(Up { _buf: buf, raw })
}

fn flat(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

/// Reinterprets the address of a flat f64 array as one of 3-vectors.
fn raw3(u: &Up<f64>) -> Raw<[f64; 3]> {
    Raw::new(u.raw.ptr, u.raw.len as usize / 3)
}

/// Reinterprets the address of a flat f64 array as one of 4-vectors.
fn raw4(u: &Up<f64>) -> Raw<[f64; 4]> {
    Raw::new(u.raw.ptr, u.raw.len as usize / 4)
}

/// The radiation transport on one GPU.
pub struct CudaRadiation {
    ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    f: CudaFunction,
    /// The walk mesh on the device (kept alive for the addresses in `mesh`).
    _f64: Vec<Up<f64>>,
    _u32: Vec<Up<u32>>,
    _u8: Vec<Up<u8>>,
    /// The argument block with the walk mesh filled in.
    mesh: TransportArgs,
}

impl CudaRadiation {
    /// Uploads the walk mesh to GPU `ordinal` and loads the PTX module.
    pub fn new(m: &WalkMesh<'_>, ordinal: usize, ptx: &Path) -> Result<Self> {
        ensure!(
            ptx.exists(),
            "PTX module `{}` not found; build it with python3 gpu/ptx-kernels/build.py",
            ptx.display()
        );
        let ctx =
            CudaContext::new(ordinal).context("CUDA device (is the NVIDIA driver installed?)")?;
        let module = ctx
            .load_module(Ptx::from_src(std::fs::read_to_string(ptx)?))
            .context("loading the PTX module")?;
        let f = module
            .load_function("radiation_transport")
            .context("entry point radiation_transport (rebuild the PTX module)")?;
        let s = ctx.default_stream();
        let points = up(&s, &flat(m.points))?;
        let apex = up(&s, &flat(m.apex))?;
        let fan = up(&s, &flat(m.fan))?;
        let face_ptr = up(&s, m.face_ptr)?;
        let face_nodes = up(&s, m.face_nodes)?;
        let owner = up(&s, m.owner)?;
        let neighbour = up(&s, m.neighbour)?;
        let cell_ptr = up(&s, m.cell_ptr)?;
        let cell_faces = up(&s, m.cell_faces)?;
        let plane_flat: Vec<f64> = m.plane.iter().flatten().copied().collect();
        let plane = up(&s, &plane_flat)?;
        let convex = up(&s, m.convex)?;
        let mesh = TransportArgs {
            points: raw3(&points),
            face_ptr: face_ptr.raw,
            face_nodes: face_nodes.raw,
            owner: owner.raw,
            neighbour: neighbour.raw,
            cell_ptr: cell_ptr.raw,
            cell_faces: cell_faces.raw,
            apex: raw3(&apex),
            fan: raw3(&fan),
            plane: raw4(&plane),
            convex: convex.raw,
            n_cells: m.n_cells as u64,
            periodic_start: m.periodic_start as u64,
            ..TransportArgs::default()
        };
        Ok(Self {
            ctx,
            stream: s,
            f,
            _f64: vec![points, apex, fan, plane],
            _u32: vec![face_ptr, face_nodes, owner, neighbour, cell_ptr, cell_faces],
            _u8: vec![convex],
            mesh,
        })
    }

    /// Name and compute capability of the device.
    pub fn device(&self) -> String {
        let name = self.ctx.name().unwrap_or_else(|_| "NVIDIA GPU".into());
        let (ma, mi) = self.ctx.compute_capability().unwrap_or((0, 0));
        format!("{name} (sm_{ma}{mi})")
    }
}

impl Transport for CudaRadiation {
    fn name(&self) -> &str {
        "cuda"
    }

    fn run(
        &mut self,
        m: &WalkMesh<'_>,
        md: &Medium<'_>,
        words: &[u64],
        call: u64,
        batches: usize,
    ) -> Result<DeviceTallies> {
        ensure!(
            m.n_cells as u64 == self.mesh.n_cells && m.points.len() as u64 == self.mesh.points.len,
            "radiation launch on another mesh than the uploaded one"
        );
        ensure!(
            words.len().is_multiple_of(BUNDLE_WORDS),
            "incomplete bundle records"
        );
        let n = words.len() / BUNDLE_WORDS;
        let (nc, nb) = (m.n_cells, md.normal.len());
        let s = self.stream.clone();
        let kappa = up(&s, md.kappa)?;
        let wall: Vec<u8> = md.wall.iter().map(|&k| k as u8).collect();
        let wall = up(&s, &wall)?;
        let emissivity = up(&s, md.emissivity)?;
        let normal = up(&s, &flat(md.normal))?;
        let bundles = up(&s, words)?;
        let cell = zeros(&s, batches * nc)?;
        let incident = zeros(&s, batches * nb)?;
        let absorbed = zeros(&s, batches * nb)?;
        let ledger = zeros(&s, 3)?;
        let a = TransportArgs {
            kappa: kappa.raw,
            wall: wall.raw,
            emissivity: emissivity.raw,
            normal: raw3(&normal),
            first_boundary: md.first_boundary as u64,
            reach: md.reach,
            max_steps: u64::from(md.max_steps),
            bundles: bundles.raw,
            n_bundles: n as u64,
            call,
            n_boundary: nb as u64,
            cell: cell.raw,
            incident: incident.raw,
            absorbed: absorbed.raw,
            ledger: ledger.raw,
            ..self.mesh
        };
        if n > 0 {
            let cfg = LaunchConfig {
                grid_dim: ((n as u32).div_ceil(THREADS), 1, 1),
                block_dim: (THREADS, 1, 1),
                shared_mem_bytes: 0,
            };
            let arg = Arg(a);
            let mut l = s.launch_builder(&self.f);
            l.arg(&arg);
            // SAFETY: radiation_transport takes (TransportArgs); every
            // address in it is a device buffer of this launch or of the
            // uploaded walk mesh, alive until after the synchronisation.
            unsafe { l.launch(cfg) }?;
        }
        s.synchronize()?;
        let back = |u: &Up<u64>, len: usize| -> Result<Vec<i64>> {
            let v = s.clone_dtoh(&u._buf)?;
            Ok(v[..len].iter().map(|&x| x as i64).collect())
        };
        let l = back(&ledger, 3)?;
        Ok(DeviceTallies {
            cell: back(&cell, batches * nc)?,
            incident: back(&incident, batches * nb)?,
            absorbed: back(&absorbed, batches * nb)?,
            ledger: [l[0], l[1], l[2]],
        })
    }
}
