//! The discretisation of a scalar transport equation on one mesh: all
//! arrays the kernel bodies of `foehncraft_kernels::fv` read, in
//! structure-of-arrays form, computed once on the host.
//!
//! - topology (owner, neighbour, cell-to-face lists) copied from the mesh;
//! - geometry: face and cell centroids, the diffusion factor |S|²/(d·S), the
//!   non-orthogonal vector k = S − d |S|²/(d·S) and the interpolation weight
//!   of every face;
//! - least-squares gradient weights per cell-face entry (D24);
//! - boundary coefficients a, b of the face values (see `foehncraft_kernels::fv`);
//! - the sliced-ELL matrix structure (D37): slot 0 of a row is the diagonal,
//!   then the interior faces of the cell in stored order; slices of C = 32
//!   rows, σ = 1 (rows keep the cell order, so the assembly writes rows in
//!   place and the GPU loads of a warp are coalesced).

use crate::problem::{BoundaryKind, Problem};
use anyhow::{Result, bail, ensure};
use foehncraft_kernels::unpack_face;
use foehncraft_mesh::geometry::{dot, norm, scale, sub};
use foehncraft_mesh::{Geometry, Mesh};

/// Slice height of the sliced-ELL matrix: the warp width of NVIDIA GPUs.
pub const CHUNK: usize = 32;

/// Three components in structure-of-arrays form.
#[derive(Debug, Clone, Default)]
pub struct Soa3 {
    /// x components.
    pub x: Vec<f64>,
    /// y components.
    pub y: Vec<f64>,
    /// z components.
    pub z: Vec<f64>,
}

impl Soa3 {
    /// Zero vectors of length `n`.
    pub fn zeros(n: usize) -> Self {
        Self {
            x: vec![0.0; n],
            y: vec![0.0; n],
            z: vec![0.0; n],
        }
    }

    /// From an array of vectors.
    pub fn from_aos(v: &[[f64; 3]]) -> Self {
        Self {
            x: v.iter().map(|a| a[0]).collect(),
            y: v.iter().map(|a| a[1]).collect(),
            z: v.iter().map(|a| a[2]).collect(),
        }
    }

    /// Vector `i`.
    pub fn get(&self, i: usize) -> [f64; 3] {
        [self.x[i], self.y[i], self.z[i]]
    }

    /// Sets vector `i`.
    pub fn set(&mut self, i: usize, v: [f64; 3]) {
        self.x[i] = v[0];
        self.y[i] = v[1];
        self.z[i] = v[2];
    }

    /// Borrowed form for the kernel bodies.
    pub fn view(&self) -> foehncraft_kernels::fv::Vec3<'_, f64> {
        foehncraft_kernels::fv::Vec3 {
            x: &self.x,
            y: &self.y,
            z: &self.z,
        }
    }

    /// Length.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// True when empty.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

/// Sliced-ELL structure of the matrix.
#[derive(Debug, Clone, Default)]
pub struct SellStructure {
    /// Slice height C.
    pub chunk: usize,
    /// Start of each slice in `col`.
    pub slice_ptr: Vec<u32>,
    /// Slots per row in each slice (the longest row of the slice).
    pub slice_len: Vec<u32>,
    /// Column of every slot; padding slots point at the row itself.
    pub col: Vec<u32>,
    /// Number of real (non-padding) entries.
    pub nnz: usize,
}

impl SellStructure {
    /// Number of stored slots including padding.
    pub fn stored(&self) -> usize {
        self.col.len()
    }

    /// Storage position of slot `j` of row `r`.
    pub fn pos(&self, r: usize, j: usize) -> usize {
        let s = r / self.chunk;
        self.slice_ptr[s] as usize + j * self.chunk + r % self.chunk
    }
}

/// All host arrays of the scalar discretisation.
#[derive(Debug, Clone, Default)]
pub struct Discretisation {
    /// Local (own) cells: the rows of the matrix.
    pub n_cells: usize,
    /// Halo cells after the own cells (0 on one partition; MPI step W3):
    /// the cell arrays of the discretisation (centroid, volume) and the
    /// cell vectors of a backend have `n_cells + n_halo` entries.
    pub n_halo: usize,
    /// Faces.
    pub n_faces: usize,
    /// Interior faces.
    pub n_interior: usize,
    /// Owner per face.
    pub owner: Vec<u32>,
    /// Neighbour per interior face.
    pub neighbour: Vec<u32>,
    /// Cell-to-face pointer.
    pub cell_ptr: Vec<u32>,
    /// Packed cell-to-face entries.
    pub cell_faces: Vec<u32>,
    /// Cell centroids.
    pub cell_centroid: Soa3,
    /// Cell volumes.
    pub cell_volume: Vec<f64>,
    /// Face centroids.
    pub face_centroid: Soa3,
    /// Face area vectors.
    pub face_area: Soa3,
    /// |S|²/(d·S) per face.
    pub dgeo: Vec<f64>,
    /// Non-orthogonal vector k per face.
    pub k: Soa3,
    /// Owner weight of the face interpolation.
    pub weight: Vec<f64>,
    /// Least-squares weights per cell-face entry.
    pub lsq: Soa3,
    /// Boundary coefficient a per boundary face.
    pub bc_a: Vec<f64>,
    /// Boundary coefficient b per boundary face.
    pub bc_b: Vec<f64>,
    /// Mass flux per face.
    pub mass_flux: Vec<f64>,
    /// Diffusion coefficient per face.
    pub gamma: Vec<f64>,
    /// Source per cell, integrated over the cell (units of φ times kg/s).
    pub source: Vec<f64>,
    /// Matrix structure.
    pub sell: SellStructure,
    /// Largest non-orthogonality angle in degrees (for the report).
    pub max_nonorth_deg: f64,
    /// First periodic face (`n_interior` without periodic faces; D29).
    pub periodic: usize,
    /// Shift of the neighbour centroid per periodic face
    /// (`foehncraft_mesh::Periodic::shift`).
    pub shift: Soa3,
    /// Partner face per periodic face.
    pub partner: Vec<u32>,
}

/// Inverse of a symmetric 3×3 matrix by the adjugate, `None` when it is
/// singular relative to its size.
fn inverse3(m: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c01 = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    let c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    let det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    let size = (m[0][0] + m[1][1] + m[2][2]).powi(3);
    if det.abs() <= 1e-12 * size || !det.is_finite() {
        return None;
    }
    let inv = 1.0 / det;
    Some([
        [
            c00 * inv,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv,
        ],
        [
            c01 * inv,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv,
        ],
        [
            c02 * inv,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv,
        ],
    ])
}

impl Discretisation {
    /// The centroid of the neighbour of interior face `f` as the face sees
    /// it: moved by the period next to a periodic face (D29).
    pub fn neighbour_centroid(&self, f: usize) -> [f64; 3] {
        let n = self.neighbour[f] as usize;
        let c = &self.cell_centroid;
        let x = [c.x[n], c.y[n], c.z[n]];
        if f >= self.periodic && f < self.n_interior {
            let i = f - self.periodic;
            [
                x[0] + self.shift.x[i],
                x[1] + self.shift.y[i],
                x[2] + self.shift.z[i],
            ]
        } else {
            x
        }
    }

    /// Builds the discretisation of `problem` on `mesh`.
    ///
    /// On a rank's part of a decomposed mesh (MPI step W3) `geo` must hold
    /// the volumes and centroids of the halo cells
    /// ([`foehncraft_mesh::decompose::LocalMesh::geometry`]). Every face,
    /// processor faces included, keeps its global orientation, so the
    /// arrays of a face (diffusion factor, k, weight) are bitwise those of
    /// the one-rank run, and the least-squares weights of an own cell sum
    /// its faces in global face order; the matrix rows are the own cells,
    /// their columns may be halo cells.
    pub fn new(mesh: &Mesh, geo: &Geometry, problem: &Problem) -> Result<Self> {
        let nc = mesh.n_cells;
        let nt = nc + mesh.n_halo;
        ensure!(
            geo.cell_centroid.len() >= nt && geo.cell_volume.len() >= nt,
            "geometry without the {} halo cells of the mesh",
            mesh.n_halo
        );
        ensure!(
            mesh.n_halo == 0 || mesh.periodic.is_none(),
            "periodic faces run on one partition"
        );
        let nf = mesh.n_faces();
        let ni = mesh.n_interior_faces();
        let nb = nf - ni;

        // boundary kinds and values per boundary face
        let mut kind = vec![BoundaryKind::ZeroGradient; nb];
        let mut value = vec![0.0; nb];
        for p in &mesh.boundaries {
            let Some(bc) = problem.boundary.get(&p.name) else {
                bail!(
                    "no boundary condition for boundary `{}` (boundaries: {})",
                    p.name,
                    mesh.boundaries
                        .iter()
                        .map(|p| p.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                );
            };
            for f in p.range() {
                kind[f - ni] = bc.kind;
                let (n, _) = geo.normal_area(f);
                value[f - ni] = bc.value_at(geo.face_centroid[f], n);
            }
        }
        for name in problem.boundary.keys() {
            ensure!(
                mesh.boundary(name).is_some(),
                "boundary condition for `{name}`, but the mesh has no such boundary"
            );
        }

        let xc = &geo.cell_centroid;
        let periodic = mesh.periodic_start();
        let pshift: Vec<[f64; 3]> = mesh
            .periodic
            .as_ref()
            .map(|p| p.shift.clone())
            .unwrap_or_default();
        // the neighbour centroid seen from interior face f (moved next to a
        // periodic face)
        let xn = |f: usize| {
            let x = xc[mesh.neighbour[f] as usize];
            if f >= periodic {
                foehncraft_mesh::geometry::add(x, pshift[f - periodic])
            } else {
                x
            }
        };
        let mut dgeo = vec![0.0; nf];
        let mut k = Soa3::zeros(nf);
        let mut weight = vec![1.0; nf];
        let mut max_nonorth = 0.0f64;
        let (mut pair_nonorth, mut shared) = (0.0f64, false);
        let mut bc_a = vec![0.0; nb];
        let mut bc_b = vec![0.0; nb];
        for f in 0..nf {
            let s = geo.face_area[f];
            let p = mesh.owner[f] as usize;
            let d = if f < ni {
                sub(xn(f), xc[p])
            } else {
                sub(geo.face_centroid[f], xc[p])
            };
            let ds = dot(d, s);
            ensure!(
                ds > 0.0,
                "face {f}: cell centre distance not along the face normal"
            );
            let ss = dot(s, s);
            dgeo[f] = ss / ds;
            k.set(f, sub(s, scale(d, ss / ds)));
            let cosang = (ds / (norm(d) * ss.sqrt())).clamp(-1.0, 1.0);
            max_nonorth = max_nonorth.max(cosang.acos().to_degrees());
            if f < ni {
                let n = scale(s, 1.0 / ss.sqrt());
                let dpf = dot(sub(geo.face_centroid[f], xc[p]), n);
                let dfn = dot(sub(xn(f), geo.face_centroid[f]), n);
                weight[f] = dfn / (dpf + dfn);
            } else {
                let fb = f - ni;
                match kind[fb] {
                    BoundaryKind::FixedValue => {
                        bc_a[fb] = 0.0;
                        bc_b[fb] = value[fb];
                    }
                    BoundaryKind::FixedGradient | BoundaryKind::ZeroGradient => {
                        bc_a[fb] = 1.0;
                        bc_b[fb] = value[fb] * ds / ss.sqrt();
                    }
                }
            }
        }

        // two cells that share several faces (the pieces of a merged cut
        // cell around a bend, L13 D228): the two-point flux of the pair
        // uses the summed area vector, |S|²/(d·S) split over the pieces by
        // area, so that Σ dgeo and Σ k are those of one face with the exact
        // geometry of the pieces (a small piece nearly parallel to d alone
        // would get a large implicit part and an explicit correction ten
        // times its area). Faces of one pair are adjacent in face order
        // (`MeshBuilder` sorts them by owner and neighbour).
        let mut f = 0;
        while f < periodic {
            let (p, n) = (mesh.owner[f], mesh.neighbour[f]);
            let mut e = f + 1;
            while e < periodic && mesh.owner[e] == p && mesh.neighbour[e] == n {
                e += 1;
            }
            if e - f > 1 {
                let d = sub(xn(f), xc[p as usize]);
                let mut s = [0.0; 3];
                let mut area = 0.0;
                for a in &geo.face_area[f..e] {
                    s = foehncraft_mesh::geometry::add(s, *a);
                    area += norm(*a);
                }
                let (ds, ss) = (dot(d, s), dot(s, s));
                let pair = ss / ds;
                // the interpolation weight of the pair, along its normal from
                // the centroid of the pieces projected onto its plane
                // ([`pair_centroid`]): along one piece's own normal it can
                // leave [0, 1] (the piece's plane need not lie between the
                // two centroids)
                let nrm = scale(s, 1.0 / ss.sqrt());
                let xm = pair_centroid(geo, f..e, nrm);
                let dpf = dot(sub(xm, xc[p as usize]), nrm);
                let dfn = dot(sub(xn(f), xm), nrm);
                let w = (dfn / (dpf + dfn)).clamp(0.0, 1.0);
                for g in f..e {
                    dgeo[g] = pair * norm(geo.face_area[g]) / area;
                    k.set(g, sub(geo.face_area[g], scale(d, dgeo[g])));
                    weight[g] = w;
                }
                let cosang = (ds / (norm(d) * ss.sqrt())).clamp(-1.0, 1.0);
                pair_nonorth = pair_nonorth.max(cosang.acos().to_degrees());
                shared = true;
            } else {
                let s = geo.face_area[f];
                let d = sub(xn(f), xc[p as usize]);
                let cosang = (dot(d, s) / (norm(d) * norm(s))).clamp(-1.0, 1.0);
                pair_nonorth = pair_nonorth.max(cosang.acos().to_degrees());
            }
            f = e;
        }
        if shared {
            // the angle the flux sees: that of the pair
            let mut bmax = 0.0f64;
            for f in periodic..nf {
                let s = geo.face_area[f];
                let p = mesh.owner[f] as usize;
                let d = if f < ni {
                    sub(xn(f), xc[p])
                } else {
                    sub(geo.face_centroid[f], xc[p])
                };
                let cosang = (dot(d, s) / (norm(d) * norm(s))).clamp(-1.0, 1.0);
                bmax = bmax.max(cosang.acos().to_degrees());
            }
            max_nonorth = pair_nonorth.max(bmax);
        }

        // least-squares weights
        let lsq = lsq_weights(mesh, geo, &kind)?;

        // physics per face and cell
        let mass_flux: Vec<f64> = (0..nf)
            .map(|f| {
                problem.density * dot(problem.velocity(geo.face_centroid[f]), geo.face_area[f])
            })
            .collect();
        let gamma = vec![problem.diffusivity; nf];
        let source: Vec<f64> = (0..nc)
            .map(|c| problem.source(xc[c]) * geo.cell_volume[c])
            .collect();

        let sell = sell_structure(mesh);
        Ok(Self {
            n_cells: nc,
            n_halo: mesh.n_halo,
            n_faces: nf,
            n_interior: ni,
            owner: mesh.owner.clone(),
            neighbour: mesh.neighbour.clone(),
            cell_ptr: mesh.cell_ptr.clone(),
            cell_faces: mesh.cell_faces.clone(),
            cell_centroid: Soa3::from_aos(&xc[..nt]),
            cell_volume: geo.cell_volume[..nt].to_vec(),
            face_centroid: Soa3::from_aos(&geo.face_centroid),
            face_area: Soa3::from_aos(&geo.face_area),
            dgeo,
            k,
            weight,
            lsq,
            bc_a,
            bc_b,
            mass_flux,
            gamma,
            source,
            sell,
            max_nonorth_deg: max_nonorth,
            periodic,
            shift: Soa3::from_aos(&pshift),
            partner: mesh
                .periodic
                .as_ref()
                .map(|p| p.partner.clone())
                .unwrap_or_default(),
        })
    }

    /// The least-squares gradient weights alone: the `lsq` of
    /// [`Discretisation::new`] with the same arguments, without the rest.
    pub fn least_squares(mesh: &Mesh, geo: &Geometry, problem: &Problem) -> Result<Soa3> {
        let ni = mesh.n_interior_faces();
        let mut kind = vec![BoundaryKind::ZeroGradient; mesh.n_faces() - ni];
        for p in &mesh.boundaries {
            let Some(bc) = problem.boundary.get(&p.name) else {
                bail!("no boundary condition for boundary `{}`", p.name);
            };
            for f in p.range() {
                kind[f - ni] = bc.kind;
            }
        }
        lsq_weights(mesh, geo, &kind)
    }

    /// Borrowed topology for the kernel bodies.
    pub fn topology(&self) -> foehncraft_kernels::fv::Topology<'_> {
        foehncraft_kernels::fv::Topology {
            owner: &self.owner,
            neighbour: &self.neighbour,
            cell_ptr: &self.cell_ptr,
            cell_faces: &self.cell_faces,
            n_interior: self.n_interior,
            periodic: self.periodic,
        }
    }

    /// The assembled matrix as CSR (columns ascending), from the values of
    /// its sliced-ELL slots: slot 0 is the diagonal, then one slot per
    /// interior face of the cell.
    pub fn matrix_csr(&self, val: &[f64]) -> crate::amg::Csr {
        let s = &self.sell;
        let rows = (0..self.n_cells)
            .map(|c| {
                let lo = self.cell_ptr[c] as usize;
                let hi = self.cell_ptr[c + 1] as usize;
                let len = 1 + self.cell_faces[lo..hi]
                    .iter()
                    .filter(|&&e| unpack_face(e).0 < self.n_interior)
                    .count();
                (0..len)
                    .map(|j| {
                        let k = s.pos(c, j);
                        (s.col[k], val[k])
                    })
                    .collect()
            })
            .collect();
        crate::amg::Csr::from_rows(self.n_cells, rows)
    }

    /// Transpose map of the matrix storage (D59, see
    /// `foehncraft_kernels::adjoint`): for every slot the slot of the
    /// mirrored entry, a_PN ↔ a_NP of each interior face; the diagonal and
    /// the padding slots map to themselves.
    pub fn transpose_map(&self) -> Vec<u32> {
        let s = &self.sell;
        let mut map: Vec<u32> = (0..s.stored() as u32).collect();
        // slot of each interior face in the owner row and the neighbour row
        let mut slot = vec![[u32::MAX; 2]; self.n_interior];
        for c in 0..self.n_cells {
            let mut j = 1;
            for &e in &self.cell_faces[self.cell_ptr[c] as usize..self.cell_ptr[c + 1] as usize] {
                let (f, is_nb) = unpack_face(e);
                if f < self.n_interior {
                    slot[f][is_nb as usize] = s.pos(c, j) as u32;
                    j += 1;
                }
            }
        }
        // a periodic face has its slot in the owner row only; its mirror is
        // the slot of the partner face in the partner's owner row
        for f in self.periodic..self.n_interior {
            let q = self.partner[f - self.periodic] as usize;
            slot[f][1] = slot[q][0];
        }
        for [a, b] in slot {
            // on ranks a processor face has a slot in the row of its own
            // cell only (the halo cell has no row): no transposed slot
            if a == u32::MAX || b == u32::MAX {
                debug_assert!(self.n_halo > 0);
                continue;
            }
            map[a as usize] = b;
            map[b as usize] = a;
        }
        map
    }

    /// Borrowed boundary coefficients.
    pub fn boundary(&self) -> foehncraft_kernels::fv::Boundary<'_, f64> {
        foehncraft_kernels::fv::Boundary {
            a: &self.bc_a,
            b: &self.bc_b,
            jump: &[],
        }
    }

    /// ε² of Venkatakrishnan's limiter per cell, (K h)³ with h³ the cell
    /// volume.
    pub fn venkat_eps2(&self, k: f64) -> Vec<f64> {
        let k3 = k * k * k;
        self.cell_volume.iter().map(|v| k3 * v).collect()
    }

    /// Bytes of the device-resident data (for the memory estimate).
    pub fn bytes(&self) -> usize {
        8 * (self.cell_centroid.len() * 4
            + self.face_centroid.len() * 3
            + self.face_area.len() * 3
            + self.n_faces * 6
            + self.lsq.len() * 3
            + self.sell.stored()
            + self.n_cells * 20)
            + 4 * (self.owner.len()
                + self.neighbour.len()
                + self.cell_ptr.len()
                + self.cell_faces.len()
                + self.sell.stored())
    }
}

/// The least-squares gradient weights of `mesh`, one per entry of
/// `cell_faces`, with `kind` the boundary kind of every boundary face
/// (the `lsq` of [`Discretisation::new`]).
pub(crate) fn lsq_weights(mesh: &Mesh, geo: &Geometry, kind: &[BoundaryKind]) -> Result<Soa3> {
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let xc = &geo.cell_centroid;
    let periodic = mesh.periodic_start();
    let pshift: &[[f64; 3]] = mesh.periodic.as_ref().map_or(&[], |p| &p.shift);
    let xn = |f: usize| {
        let x = xc[mesh.neighbour[f] as usize];
        if f >= periodic {
            foehncraft_mesh::geometry::add(x, pshift[f - periodic])
        } else {
            x
        }
    };
    let ne = mesh.cell_faces.len();
    let mut lsq = Soa3::zeros(ne);
    let mut stencil: Vec<([f64; 3], f64)> = Vec::new();
    let mut other: Vec<i64> = Vec::new();
    for c in 0..nc {
        stencil.clear();
        let mut g = [[0.0; 3]; 3];
        let lo = mesh.cell_ptr[c] as usize;
        let hi = mesh.cell_ptr[c + 1] as usize;
        // the cell across each interior face: a neighbour that shares
        // several faces with the cell (merged cut cells, L13 D228) is one
        // point of the stencil, its weight split over those faces
        other.clear();
        other.extend(mesh.cell_faces[lo..hi].iter().map(|&e| {
            let (f, is_nb) = unpack_face(e);
            if f >= ni || f >= periodic {
                -1 - f as i64
            } else if is_nb {
                mesh.owner[f] as i64
            } else {
                mesh.neighbour[f] as i64
            }
        }));
        for (i, &e) in mesh.cell_faces[lo..hi].iter().enumerate() {
            let (f, is_nb) = unpack_face(e);
            let times = other.iter().filter(|&&o| o == other[i]).count() as f64;
            let d = if f < ni {
                if is_nb {
                    sub(xc[mesh.owner[f] as usize], xc[c])
                } else {
                    sub(xn(f), xc[c])
                }
            } else {
                let d = sub(geo.face_centroid[f], xc[c]);
                match kind[f - ni] {
                    BoundaryKind::FixedValue => d,
                    _ => {
                        let s = geo.face_area[f];
                        let n = scale(s, 1.0 / norm(s));
                        scale(n, dot(d, n))
                    }
                }
            };
            let w = 1.0 / (dot(d, d) * times);
            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += w * d[i] * d[j];
                }
            }
            stencil.push((d, w));
        }
        let Some(inv) = inverse3(g) else {
            bail!("cell {c}: least-squares stencil is degenerate");
        };
        for (i, (d, w)) in stencil.iter().enumerate() {
            let mut v = [0.0; 3];
            for (r, vr) in v.iter_mut().enumerate() {
                *vr = w * (inv[r][0] * d[0] + inv[r][1] * d[1] + inv[r][2] * d[2]);
            }
            lsq.set(lo + i, v);
        }
    }
    Ok(lsq)
}

/// The centroid of the faces `faces` of one cell pair projected onto the
/// plane of their summed area vector (unit normal `n`): each piece weighted
/// by its area seen along `n`. Pieces in one plane give their joint
/// centroid; a small piece nearly parallel to `n` hardly counts.
pub fn pair_centroid(geo: &Geometry, faces: std::ops::Range<usize>, n: [f64; 3]) -> [f64; 3] {
    let (mut x, mut w) = ([0.0; 3], 0.0);
    for g in faces.clone() {
        let a = dot(geo.face_area[g], n).max(0.0);
        x = foehncraft_mesh::geometry::add(x, scale(geo.face_centroid[g], a));
        w += a;
    }
    if w > 0.0 {
        scale(x, 1.0 / w)
    } else {
        geo.face_centroid[faces.start]
    }
}

/// Sliced-ELL structure of the matrix of `mesh` (see the module docs).
pub fn sell_structure(mesh: &Mesh) -> SellStructure {
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let row_len = |c: usize| 1 + mesh.cell_face_iter(c).filter(|&(f, _)| f < ni).count();
    let n_slices = nc.div_ceil(CHUNK);
    let mut slice_ptr = Vec::with_capacity(n_slices);
    let mut slice_len = Vec::with_capacity(n_slices);
    let mut total = 0usize;
    for s in 0..n_slices {
        let len = (s * CHUNK..((s + 1) * CHUNK).min(nc))
            .map(row_len)
            .max()
            .unwrap_or(1);
        slice_ptr.push(total as u32);
        slice_len.push(len as u32);
        total += len * CHUNK;
    }
    let mut col = vec![0u32; total];
    // padding (also the lanes past the last row) points at a valid row
    for s in 0..n_slices {
        for j in 0..slice_len[s] as usize {
            for lane in 0..CHUNK {
                let r = (s * CHUNK + lane).min(nc - 1);
                col[slice_ptr[s] as usize + j * CHUNK + lane] = r as u32;
            }
        }
    }
    let mut nnz = 0;
    for c in 0..nc {
        let s = c / CHUNK;
        let pos = |j: usize| slice_ptr[s] as usize + j * CHUNK + c % CHUNK;
        col[pos(0)] = c as u32;
        let mut j = 1;
        for (f, is_nb) in mesh.cell_face_iter(c) {
            if f < ni {
                col[pos(j)] = if is_nb {
                    mesh.owner[f]
                } else {
                    mesh.neighbour[f]
                };
                j += 1;
            }
        }
        nnz += j;
    }
    SellStructure {
        chunk: CHUNK,
        slice_ptr,
        slice_len,
        col,
        nnz,
    }
}
