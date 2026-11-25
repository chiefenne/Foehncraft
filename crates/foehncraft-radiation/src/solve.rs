//! One radiation call on a mesh or on one rank's part of it (R1a).

use crate::device::{DeviceTallies, Transport};
use foehncraft_core::exact::ExactSum;
use foehncraft_core::rng::{keyed, stream, uniform_keyed};
use foehncraft_kernels::radiation::{
    BUNDLE_WORDS, Bundle, Fate, Medium, Tally, WallKind, cosine_about, isotropic, owen_sobol4,
    point_in_tet, point_in_triangle, sobol_seed, transport,
};
use foehncraft_kernels::walk::{Walk, WalkMesh, inner};
use foehncraft_mesh::decompose::LocalMesh;
use foehncraft_mesh::geometry::{norm, scale};
use foehncraft_mesh::handover::HandOver;
use foehncraft_mesh::locate::{Locator, step_limit};
use foehncraft_mesh::{Geometry, Mesh};
use foehncraft_parallel::{Communicator, SerialComm};
use rayon::prelude::*;
use std::sync::Mutex;
use std::sync::atomic::{AtomicI64, Ordering};

/// Stefan–Boltzmann constant (W/(m² K⁴), CODATA 2018, exact in the SI).
pub const SIGMA: f64 = 5.670_374_419e-8;

/// Quanta of the total emitted power of a call: 2⁵² (D107 (3)).
const QUANTA: f64 = 4_503_599_627_370_496.0;

/// Numbers after the walk in a hand-over record: w, w0, emitter, index,
/// event, batch.
const PAYLOAD: usize = 6;

/// The radiation condition of one boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wall {
    /// Grey diffuse wall, symmetry plane (mirror) or open boundary.
    pub kind: WallKind,
    /// Emissivity (grey walls; 1 for black).
    pub emissivity: f64,
    /// Temperature (K) of a grey wall; an open boundary radiates nothing
    /// in (surroundings at 0 K).
    pub temperature: f64,
}

impl Wall {
    /// A black wall at temperature `t`.
    pub fn black(t: f64) -> Self {
        Self::grey(1.0, t)
    }

    /// A grey diffuse wall of emissivity `e` at temperature `t`.
    pub fn grey(e: f64, t: f64) -> Self {
        assert!(e > 0.0 && e <= 1.0, "emissivity in (0, 1]");
        Self {
            kind: WallKind::Grey,
            emissivity: e,
            temperature: t,
        }
    }

    /// A symmetry plane.
    pub fn mirror() -> Self {
        Self {
            kind: WallKind::Mirror,
            emissivity: 0.0,
            temperature: 0.0,
        }
    }

    /// An open boundary to cold surroundings.
    pub fn open() -> Self {
        Self {
            kind: WallKind::Open,
            emissivity: 0.0,
            temperature: 0.0,
        }
    }
}

/// Settings of a call.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settings {
    /// Number of bundles over the whole domain (about; each emitter gets
    /// its share of the power, at least one).
    pub bundles: u64,
    /// Batches of the statistical error estimate.
    pub batches: u32,
    /// Number of the call (keys the random numbers, D107 (2)).
    pub call: u64,
    /// Bundles per launch on a device (2²⁰: about 180 MB of records).
    pub launch: u64,
    /// Points of the emission: scrambled Sobol (the default) or
    /// pseudo-random, the reference of the benchmark (case 6).
    pub points: Points,
}

/// The points that place and aim the bundles at emission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Points {
    /// Owen-scrambled Sobol points per emitter (Burley 2020): the bundles of
    /// an emitter are stratified in position and direction.
    #[default]
    Sobol,
    /// Independent keyed pseudo-random numbers (D121), the plain Monte
    /// Carlo estimate.
    Random,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            bundles: 100_000,
            batches: 8,
            call: 0,
            launch: 1 << 20,
            points: Points::Sobol,
        }
    }
}

/// The energy account of a call, in quanta over all ranks: emitted =
/// absorbed by the gas + absorbed by the walls + escaped + roulette + lost,
/// exactly (D107 (4)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Ledger {
    /// Emitted by the gas and the walls.
    pub emitted: i64,
    /// Absorbed by the gas.
    pub absorbed_gas: i64,
    /// Absorbed by the walls.
    pub absorbed_walls: i64,
    /// Left through open boundaries.
    pub escaped: i64,
    /// Removed (positive) or added (negative) by Russian roulette; zero in
    /// expectation.
    pub roulette: i64,
    /// Carried by bundles lost in degenerate geometry (0).
    pub lost: i64,
}

impl Ledger {
    /// Whether the account closes to the quantum.
    pub fn closes(&self) -> bool {
        self.emitted
            == self.absorbed_gas + self.absorbed_walls + self.escaped + self.roulette + self.lost
    }
}

/// The result of a call on the own cells and boundary faces of a rank.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Solution {
    /// Power of one quantum (W).
    pub quantum: f64,
    /// Net radiative source per own cell, absorbed − emitted (W/m³).
    pub source: Vec<f64>,
    /// Standard error of the source (W/m³), from the batch means.
    pub source_error: Vec<f64>,
    /// Absorbed power per own cell (quanta).
    pub absorbed: Vec<i64>,
    /// Emitted power per own cell (quanta).
    pub emitted: Vec<i64>,
    /// Incident flux per boundary face (W/m²), indexed from the first
    /// boundary face; on an open boundary the flux that leaves through it.
    pub wall_incident: Vec<f64>,
    /// Net flux into the wall per boundary face, absorbed − emitted (W/m²).
    pub wall_net: Vec<f64>,
    /// Standard error of the net wall flux (W/m²).
    pub wall_net_error: Vec<f64>,
    /// Absorbed power per boundary face (quanta).
    pub wall_absorbed: Vec<i64>,
    /// Emitted power per boundary face (quanta).
    pub wall_emitted: Vec<i64>,
    /// The energy account over all ranks.
    pub ledger: Ledger,
    /// Bundles emitted over all ranks.
    pub bundles: i64,
    /// Hand-overs between ranks over all ranks.
    pub hand_overs: i64,
}

/// One rank's link to the others.
struct Link<'a> {
    local: &'a LocalMesh,
    comm: &'a dyn Communicator,
}

/// What emits: a cell or a boundary face (local index), its global
/// identity and power.
#[derive(Debug, Clone, Copy)]
struct Emitter {
    id: u64,
    cell: u32,
    /// Boundary face (local index) or `u32::MAX` for a cell.
    face: u32,
    quanta: i64,
    bundles: u32,
}

/// The mesh side of radiation calls: the whole mesh on one rank, or one
/// rank's local mesh.
pub struct Domain<'a> {
    mesh: &'a Mesh,
    geo: &'a Geometry,
    locator: Locator<'a>,
    /// Global position of each own cell.
    position: Vec<u32>,
    n_global: usize,
    link: Option<Link<'a>>,
    shapes: Shapes,
    /// Boundary of each boundary face.
    boundary: Vec<u32>,
    max_steps: u32,
    /// The device that moves the bundles, if not the CPU threads.
    device: Option<Mutex<Box<dyn Transport + 'a>>>,
}

/// Where bundles start: the tetrahedra of the cells and the fan triangles
/// of the boundary faces with their cumulative fractions.
struct Shapes {
    /// Tetrahedra of each own cell: (face, edge) and the cumulative volume
    /// fraction.
    tet_ptr: Vec<u32>,
    tet: Vec<(u32, u32)>,
    tet_cdf: Vec<f64>,
    /// Cumulative area fraction of the fan triangles of each boundary face.
    tri_cdf: Vec<f64>,
    tri_ptr: Vec<u32>,
    /// Unit normal into the domain of each boundary face.
    normal: Vec<[f64; 3]>,
    /// First boundary face.
    first: usize,
    /// Length of a ray segment, longer than any chord of the domain.
    reach: f64,
}

impl<'a> Domain<'a> {
    /// The whole mesh on one rank.
    pub fn whole(mesh: &'a Mesh, geo: &'a Geometry) -> Self {
        Self::build(
            mesh,
            geo,
            None,
            (0..mesh.n_cells as u32).collect(),
            mesh.n_cells,
        )
    }

    /// One rank's part: `local` with its geometry from
    /// [`LocalMesh::geometry`], the communicator of the ranks and the
    /// number of cells of the global mesh.
    pub fn rank(local: &'a LocalMesh, geo: &'a Geometry, comm: &'a dyn Communicator) -> Self {
        let n_global = {
            let mut n = [local.mesh.n_cells as i64];
            comm.all_sum_i64(&mut n);
            n[0] as usize
        };
        let position = local.cell_position[..local.mesh.n_cells].to_vec();
        Self::build(
            &local.mesh,
            geo,
            Some(Link { local, comm }),
            position,
            n_global,
        )
    }

    fn comm(&self) -> &dyn Communicator {
        match &self.link {
            Some(l) => l.comm,
            None => &SerialComm,
        }
    }

    fn build(
        mesh: &'a Mesh,
        geo: &'a Geometry,
        link: Option<Link<'a>>,
        position: Vec<u32>,
        n_global: usize,
    ) -> Self {
        let locator = Locator::new(mesh, geo);
        let wm = locator.walk_mesh();
        // tetrahedra of the cells by volume
        let (mut tet_ptr, mut tet, mut tet_cdf) = (vec![0u32], Vec::new(), Vec::new());
        for c in 0..mesh.n_cells {
            let start = tet.len();
            let mut sum = 0.0;
            for (f, _) in mesh.cell_face_iter(c) {
                for e in 0..mesh.face(f).len() {
                    let t = wm.tet(c, f, e);
                    sum += inner(t[1], t[2], t[3], t[0]).abs();
                    tet.push((f as u32, e as u32));
                    tet_cdf.push(sum);
                }
            }
            for v in &mut tet_cdf[start..] {
                *v /= sum;
            }
            tet_ptr.push(tet.len() as u32);
        }
        // fan triangles of the boundary faces by area, inward normals
        let first = mesh.n_interior_faces();
        let (mut tri_ptr, mut tri_cdf) = (vec![0u32], Vec::new());
        let mut normal = Vec::new();
        let mut boundary = Vec::new();
        for (f, &s) in geo
            .face_area
            .iter()
            .enumerate()
            .take(mesh.n_faces())
            .skip(first)
        {
            let nodes = mesh.face(f);
            let k = nodes.len();
            let m = node_average(mesh, nodes);
            let start = tri_cdf.len();
            let mut sum = 0.0;
            for e in 0..k {
                let p = mesh.points[nodes[e] as usize];
                let q = mesh.points[nodes[(e + 1) % k] as usize];
                sum += 0.5 * norm(cross(sub(p, m), sub(q, m)));
                tri_cdf.push(sum);
            }
            for v in &mut tri_cdf[start..] {
                *v /= sum;
            }
            tri_ptr.push(tri_cdf.len() as u32);
            normal.push(scale(s, -1.0 / norm(s)));
            boundary.push(mesh.boundary_of(f).expect("boundary face") as u32);
        }
        // reach: twice the diagonal of the global box, plus one
        let (lo, hi) = mesh.bounding_box();
        let comm: &dyn Communicator = match &link {
            Some(l) => l.comm,
            None => &SerialComm,
        };
        let mut d2 = 0.0;
        for k in 0..3 {
            let e = comm.max(hi[k]) - comm.min(lo[k]);
            d2 += e * e;
        }
        let reach = 2.0 * d2.sqrt() + 1.0;
        Self {
            mesh,
            geo,
            locator,
            position,
            n_global,
            link,
            shapes: Shapes {
                tet_ptr,
                tet,
                tet_cdf,
                tri_cdf,
                tri_ptr,
                normal,
                first,
                reach,
            },
            boundary,
            max_steps: step_limit(n_global),
            device: None,
        }
    }

    /// Moves the bundles on `device` instead of the CPU threads (D259):
    /// the bundles are emitted on the host and transported by the device
    /// in launches of [`Settings::launch`] bundles. One device per run: not on ranks
    /// (GPU with ranks comes with version 1.1).
    pub fn on_device(mut self, device: Box<dyn Transport + 'a>) -> anyhow::Result<Self> {
        anyhow::ensure!(
            self.link.is_none(),
            "radiation on a device runs on one rank (GPU with ranks comes in version 1.1)"
        );
        self.device = Some(Mutex::new(device));
        Ok(self)
    }

    /// Moves the bundles on the registered GPU backend, device `ordinal`,
    /// with the kernel module `ptx` (`python3 gpu/ptx-kernels/build.py`).
    pub fn on_gpu(self, ordinal: usize, ptx: &std::path::Path) -> anyhow::Result<Self> {
        let gpu = crate::device::create_gpu(&self.locator.walk_mesh(), ordinal, ptx)?;
        self.on_device(gpu)
    }

    /// The name of the device that moves the bundles (`cpu` for the CPU
    /// threads).
    pub fn device_name(&self) -> String {
        match &self.device {
            Some(d) => d
                .lock()
                .map_or_else(|_| "poisoned".into(), |d| d.name().to_string()),
            None => "cpu".into(),
        }
    }

    /// Number of boundary faces of the mesh.
    pub fn n_boundary_faces(&self) -> usize {
        self.shapes.normal.len()
    }

    /// Boundary of each boundary face (index from the first boundary face).
    pub fn boundary_of_faces(&self) -> &[u32] {
        &self.boundary
    }

    /// Solves one call: temperature (K) and absorption coefficient (1/m)
    /// per own cell, the radiation condition per boundary of the mesh.
    pub fn solve(
        &self,
        temperature: &[f64],
        kappa: &[f64],
        walls: &[Wall],
        s: &Settings,
    ) -> anyhow::Result<Solution> {
        let mesh = self.mesh;
        let (nc, nb) = (mesh.n_cells, self.shapes.normal.len());
        assert_eq!(temperature.len(), nc, "one temperature per own cell");
        assert_eq!(kappa.len(), nc, "one absorption coefficient per own cell");
        assert_eq!(walls.len(), mesh.boundaries.len(), "one wall per boundary");
        assert!(s.batches >= 2, "at least two batches");
        let comm = self.comm();
        let first = mesh.n_interior_faces();
        let wall_of = |f: usize| walls[self.boundary[f] as usize];

        // emitted powers and their exact total over the ranks
        let t4 = |t: f64| t * t * t * t;
        let cell_power: Vec<f64> = (0..nc)
            .map(|c| 4.0 * kappa[c] * SIGMA * t4(temperature[c]) * self.geo.cell_volume[c])
            .collect();
        let face_power: Vec<f64> = (0..nb)
            .map(|f| {
                let w = wall_of(f);
                if w.kind == WallKind::Grey {
                    w.emissivity * SIGMA * t4(w.temperature) * norm(self.geo.face_area[first + f])
                } else {
                    0.0
                }
            })
            .collect();
        let mut total = ExactSum::new();
        total.add_slice(&cell_power);
        total.add_slice(&face_power);
        let total = comm.sum_exact(&total);
        if total <= 0.0 {
            return Ok(self.empty(nc, nb));
        }
        let quantum = total / QUANTA;

        // emitters: quanta and bundles from their own power and the total,
        // so the counts do not depend on the partition (D258)
        let n_global = self.n_global as u64;
        let share = |id: u64, p: f64| -> (i64, u32) {
            let quanta = (p / quantum).round() as i64;
            if quanta == 0 {
                return (0, 0);
            }
            let x = s.bundles as f64 * (p / total);
            let mut n = x.floor();
            let u = (keyed([stream::RADIATION, id, (s.call << 32) | 0xffff_ffff, 0xfe]) >> 11)
                as f64
                * (1.0 / (1u64 << 53) as f64);
            if u < x - n {
                n += 1.0;
            }
            let n = (n as i64).clamp(1, quanta.min(u32::MAX as i64));
            (quanta, n as u32)
        };
        let mut emitters: Vec<Emitter> = Vec::new();
        let mut emitted = vec![0i64; nc];
        let mut wall_emitted = vec![0i64; nb];
        for c in 0..nc {
            let id = u64::from(self.position[c]);
            let (quanta, bundles) = share(id, cell_power[c]);
            emitted[c] = quanta;
            if bundles > 0 {
                emitters.push(Emitter {
                    id,
                    cell: c as u32,
                    face: u32::MAX,
                    quanta,
                    bundles,
                });
            }
        }
        for f in 0..nb {
            let id = n_global + mesh.face_id[first + f];
            let (quanta, bundles) = share(id, face_power[f]);
            wall_emitted[f] = quanta;
            if bundles > 0 {
                emitters.push(Emitter {
                    id,
                    cell: mesh.owner[first + f],
                    face: f as u32,
                    quanta,
                    bundles,
                });
            }
        }

        // transport
        let kinds: Vec<WallKind> = (0..nb).map(|f| wall_of(f).kind).collect();
        let emissivity: Vec<f64> = (0..nb).map(|f| wall_of(f).emissivity).collect();
        let medium = Medium {
            kappa,
            wall: &kinds,
            emissivity: &emissivity,
            normal: &self.shapes.normal,
            first_boundary: first,
            reach: self.shapes.reach,
            max_steps: self.max_steps,
        };
        let wm = self.locator.walk_mesh();
        let b = s.batches as usize;
        let tally = Tallies::new(b, nc, nb);
        let shapes = &self.shapes;
        let weight = |e: &Emitter, k: u32| {
            let base = e.quanta / i64::from(e.bundles);
            base + i64::from(i64::from(k) < e.quanta % i64::from(e.bundles))
        };
        let mut queue: Vec<Bundle> = Vec::new();
        if let Some(device) = &self.device {
            // emission on the host, transport on the device, in launches
            let mut device = device
                .lock()
                .map_err(|_| anyhow::anyhow!("radiation device poisoned"))?;
            let mut start = 0;
            while start < emitters.len() {
                let mut end = start;
                let mut n = 0u64;
                while end < emitters.len()
                    && (n == 0 || n + u64::from(emitters[end].bundles) <= s.launch)
                {
                    n += u64::from(emitters[end].bundles);
                    end += 1;
                }
                let words: Vec<u64> = emitters[start..end]
                    .par_iter()
                    .flat_map_iter(|e| {
                        (0..e.bundles)
                            .flat_map(move |k| shapes.emit(e, k, weight(e, k), s, &wm).to_words())
                    })
                    .collect();
                debug_assert_eq!(words.len() as u64, n * BUNDLE_WORDS as u64);
                let t = device.run(&wm, &medium, &words, s.call, b)?;
                tally.add(&t)?;
                start = end;
            }
        } else {
            queue = emitters
                .par_iter()
                .map(|e| {
                    let mut handed = Vec::new();
                    for k in 0..e.bundles {
                        let mut bundle = shapes.emit(e, k, weight(e, k), s, &wm);
                        if transport(&wm, &medium, &mut bundle, s.call, &tally) == Fate::HandOver {
                            handed.push(bundle);
                        }
                    }
                    handed
                })
                .collect::<Vec<_>>()
                .into_iter()
                .flatten()
                .collect();
        }
        let mut hand_overs = 0i64;
        if let Some(link) = &self.link {
            let hand = HandOver::new(link.local, PAYLOAD);
            let mut outbox = hand.outbox();
            loop {
                let mut in_flight = [queue.len() as i64];
                comm.all_sum_i64(&mut in_flight);
                if in_flight[0] == 0 {
                    break;
                }
                hand_overs += queue.len() as i64;
                for x in &queue {
                    let payload = [
                        x.w as f64,
                        x.w0 as f64,
                        x.emitter as f64,
                        f64::from(x.index),
                        f64::from(x.event),
                        f64::from(x.batch),
                    ];
                    hand.pack(&x.walk, &payload, &mut outbox);
                }
                let inbox = hand.exchange(comm, &mut outbox);
                let arrived: Vec<Bundle> = hand
                    .unpack(&inbox)
                    .into_iter()
                    .map(|(walk, p)| Bundle {
                        walk,
                        w: p[0] as i64,
                        w0: p[1] as i64,
                        emitter: p[2] as u64,
                        index: p[3] as u32,
                        event: p[4] as u32,
                        batch: p[5] as u32,
                    })
                    .collect();
                queue = arrived
                    .into_par_iter()
                    .filter_map(|mut x| {
                        (transport(&wm, &medium, &mut x, s.call, &tally) == Fate::HandOver)
                            .then_some(x)
                    })
                    .collect();
            }
        } else {
            assert!(
                queue.is_empty(),
                "a bundle reached a processor face on one rank"
            );
        }

        // results
        let cells = tally.cell_totals();
        let incident = tally.wall_totals(&tally.incident);
        let wall_absorbed = tally.wall_totals(&tally.absorbed);
        let source: Vec<f64> = (0..nc)
            .map(|c| (cells[c] - emitted[c]) as f64 * quantum / self.geo.cell_volume[c])
            .collect();
        let source_error = (0..nc)
            .map(|c| tally.error(&tally.cell, nc, c) * quantum / self.geo.cell_volume[c])
            .collect();
        let area: Vec<f64> = (0..nb)
            .map(|f| norm(self.geo.face_area[first + f]))
            .collect();
        let wall_incident = (0..nb)
            .map(|f| incident[f] as f64 * quantum / area[f])
            .collect();
        let wall_net = (0..nb)
            .map(|f| (wall_absorbed[f] - wall_emitted[f]) as f64 * quantum / area[f])
            .collect();
        let wall_net_error = (0..nb)
            .map(|f| tally.error(&tally.absorbed, nb, f) * quantum / area[f])
            .collect();
        let mut sums = [
            emitted.iter().sum::<i64>() + wall_emitted.iter().sum::<i64>(),
            cells.iter().sum(),
            wall_absorbed.iter().sum(),
            tally.escaped.load(Ordering::Relaxed),
            tally.roulette.load(Ordering::Relaxed),
            tally.lost.load(Ordering::Relaxed),
            emitters.iter().map(|e| i64::from(e.bundles)).sum(),
            hand_overs,
        ];
        comm.all_sum_i64(&mut sums);
        Ok(Solution {
            quantum,
            source,
            source_error,
            absorbed: cells,
            emitted,
            wall_incident,
            wall_net,
            wall_net_error,
            wall_absorbed,
            wall_emitted,
            ledger: Ledger {
                emitted: sums[0],
                absorbed_gas: sums[1],
                absorbed_walls: sums[2],
                escaped: sums[3],
                roulette: sums[4],
                lost: sums[5],
            },
            bundles: sums[6],
            hand_overs: sums[7],
        })
    }

    fn empty(&self, nc: usize, nb: usize) -> Solution {
        let comm = self.comm();
        let mut sums = [0i64; 8];
        comm.all_sum_i64(&mut sums);
        Solution {
            source: vec![0.0; nc],
            source_error: vec![0.0; nc],
            absorbed: vec![0; nc],
            emitted: vec![0; nc],
            wall_incident: vec![0.0; nb],
            wall_net: vec![0.0; nb],
            wall_net_error: vec![0.0; nb],
            wall_absorbed: vec![0; nb],
            wall_emitted: vec![0; nb],
            ..Default::default()
        }
    }
}

impl Shapes {
    /// Bundle `k` of emitter `e` with `w` quanta: position and direction
    /// from two scrambled Sobol sequences of the emitter (D107 (2)), or
    /// from keyed pseudo-random numbers.
    fn emit(&self, e: &Emitter, k: u32, w: i64, s: &Settings, wm: &WalkMesh<'_>) -> Bundle {
        let (u, v) = match s.points {
            Points::Sobol => (
                owen_sobol4(k, sobol_seed(e.id, s.call, 1)),
                owen_sobol4(k, sobol_seed(e.id, s.call, 2)),
            ),
            Points::Random => {
                let r = |d: u64| {
                    uniform_keyed([
                        stream::RADIATION,
                        e.id,
                        (s.call << 32) | u64::from(k),
                        0x100 | d,
                    ])
                };
                ([r(0), r(1), r(2), r(3)], [r(4), r(5), r(6), r(7)])
            }
        };
        let c = e.cell as usize;
        let (x, face, edge, dir) = if e.face == u32::MAX {
            let r = self.tet_ptr[c] as usize..self.tet_ptr[c + 1] as usize;
            let i = r.start + pick(&self.tet_cdf[r.clone()], u[0]);
            let (face, edge) = self.tet[i];
            let t = wm.tet(c, face as usize, edge as usize);
            let x = point_in_tet(&t, [u[1], u[2], u[3]]);
            (x, face, edge, isotropic([v[0], v[1]]))
        } else {
            let f = e.face as usize;
            let face = self.first + f;
            let r = self.tri_ptr[f] as usize..self.tri_ptr[f + 1] as usize;
            let cdf = &self.tri_cdf[r];
            let edge = pick(cdf, u[0]);
            // the rest of the first number within the chosen triangle is
            // the first coordinate in it, so that position and direction
            // are one stratified point of four dimensions
            let lo = if edge == 0 { 0.0 } else { cdf[edge - 1] };
            let a = ((u[0] - lo) / (cdf[edge] - lo)).clamp(0.0, 1.0);
            let t = wm.tet(c, face, edge);
            // the outer triangle (fan point, p, q) of the tetrahedron
            let x = point_in_triangle(t[1], t[2], t[3], [a, u[1]]);
            (
                x,
                face as u32,
                edge as u32,
                cosine_about(self.normal[f], [u[2], u[3]]),
            )
        };
        let target = std::array::from_fn(|i| x[i] + self.reach * dir[i]);
        let batch = (keyed([stream::RADIATION, e.id, (s.call << 32) | u64::from(k), 0xff])
            % u64::from(s.batches)) as u32;
        Bundle {
            walk: Walk::new(x, e.cell, face, edge, target),
            w,
            w0: w,
            emitter: e.id,
            index: k,
            event: 0,
            batch,
        }
    }
}

/// The first index whose cumulative fraction exceeds `u`.
fn pick(cdf: &[f64], u: f64) -> usize {
    cdf.iter().position(|&c| u < c).unwrap_or(cdf.len() - 1)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn node_average(mesh: &Mesh, nodes: &[u32]) -> [f64; 3] {
    foehncraft_mesh::geometry::node_average(&mesh.points, nodes)
}

/// Integer tallies of a call, one array per batch.
struct Tallies {
    batches: usize,
    cell: Vec<AtomicI64>,
    incident: Vec<AtomicI64>,
    absorbed: Vec<AtomicI64>,
    escaped: AtomicI64,
    roulette: AtomicI64,
    lost: AtomicI64,
    nc: usize,
    nb: usize,
}

impl Tallies {
    fn new(batches: usize, nc: usize, nb: usize) -> Self {
        let zeros = |n: usize| (0..n).map(|_| AtomicI64::new(0)).collect();
        Self {
            batches,
            cell: zeros(batches * nc),
            incident: zeros(batches * nb),
            absorbed: zeros(batches * nb),
            escaped: AtomicI64::new(0),
            roulette: AtomicI64::new(0),
            lost: AtomicI64::new(0),
            nc,
            nb,
        }
    }

    /// Adds the tallies of a device.
    fn add(&self, t: &DeviceTallies) -> anyhow::Result<()> {
        anyhow::ensure!(
            t.cell.len() == self.cell.len()
                && t.incident.len() == self.incident.len()
                && t.absorbed.len() == self.absorbed.len(),
            "device tallies of the wrong size"
        );
        for (a, &v) in self.cell.iter().zip(&t.cell) {
            a.fetch_add(v, Ordering::Relaxed);
        }
        for (a, &v) in self.incident.iter().zip(&t.incident) {
            a.fetch_add(v, Ordering::Relaxed);
        }
        for (a, &v) in self.absorbed.iter().zip(&t.absorbed) {
            a.fetch_add(v, Ordering::Relaxed);
        }
        self.escaped.fetch_add(t.ledger[0], Ordering::Relaxed);
        self.roulette.fetch_add(t.ledger[1], Ordering::Relaxed);
        self.lost.fetch_add(t.ledger[2], Ordering::Relaxed);
        Ok(())
    }

    fn cell_totals(&self) -> Vec<i64> {
        (0..self.nc)
            .map(|c| {
                (0..self.batches)
                    .map(|b| self.cell[b * self.nc + c].load(Ordering::Relaxed))
                    .sum()
            })
            .collect()
    }

    fn wall_totals(&self, a: &[AtomicI64]) -> Vec<i64> {
        (0..self.nb)
            .map(|f| {
                (0..self.batches)
                    .map(|b| a[b * self.nb + f].load(Ordering::Relaxed))
                    .sum()
            })
            .collect()
    }

    /// Standard error (quanta) of the total of entry `i` of a batched
    /// array of `n` entries per batch, from the batch means.
    fn error(&self, a: &[AtomicI64], n: usize, i: usize) -> f64 {
        let nb = self.batches as f64;
        let x: Vec<f64> = (0..self.batches)
            .map(|b| a[b * n + i].load(Ordering::Relaxed) as f64)
            .collect();
        let total: f64 = x.iter().sum();
        let var: f64 = x.iter().map(|&v| (nb * v - total).powi(2)).sum::<f64>() / (nb * (nb - 1.0));
        var.sqrt()
    }
}

impl Tally for Tallies {
    fn cell(&self, batch: u32, cell: usize, q: i64) {
        self.cell[batch as usize * self.nc + cell].fetch_add(q, Ordering::Relaxed);
    }

    fn wall(&self, batch: u32, face: usize, incident: i64, absorbed: i64) {
        let i = batch as usize * self.nb + face;
        self.incident[i].fetch_add(incident, Ordering::Relaxed);
        self.absorbed[i].fetch_add(absorbed, Ordering::Relaxed);
    }

    fn escaped(&self, q: i64) {
        self.escaped.fetch_add(q, Ordering::Relaxed);
    }

    fn roulette(&self, q: i64) {
        self.roulette.fetch_add(q, Ordering::Relaxed);
    }

    fn lost(&self, q: i64) {
        self.lost.fetch_add(q, Ordering::Relaxed);
    }
}
