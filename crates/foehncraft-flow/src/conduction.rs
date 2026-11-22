//! Conduction across region interfaces and in orthotropic solids (stage 7,
//! PR 7d2; D157 §9.1, §9.2, D116, L48).
//!
//! **Interface faces.** Between two regions (fluid–solid, solid–solid) the
//! temperature gradient jumps, and a contact resistance R″_c (m²K/W) or a
//! thin wall t/k may sit in the face. With the normal distances δ_P, δ_N
//! of the two centroids to the face and the face-normal temperatures
//!
//!   T̃_P = T_P + ∇T_P·(r_P − δ_P n̂),   r_P = x_f − x_P   (likewise T̃_N),
//!
//! the heat flux from P to N is the series
//!
//!   q = (T̃_P − T̃_N + c_P R_P + c_N R_N) / (R_P + R″_c + R_N),
//!   R_P = δ_P / k_n,P (or 1/h_w of a turbulent fluid side, §9.4),
//!
//! with k_n = n̂·K·n̂ and c = −(K n̂ − k_n n̂)·∇T the part of the flux of an
//! orthotropic side that its tangential gradient drives (zero for an
//! isotropic material). The face temperatures of the two sides are
//! T_f,a = T̃_P − R_P (q − c_P) and T_f,b = T̃_N + R_N (q − c_N); they
//! differ by q R″_c.
//!
//! **One-sided gradients (§9.2).** The least-squares gradient of a cell
//! next to an interface takes the face temperature of its own side at x_f
//! in place of the cell across, so that every region's gradient sees only
//! its own material. With piecewise-linear T the gradients, the face
//! temperatures and the fluxes are then exact on any mesh and at any
//! conductivity ratio; the pair (gradients, face temperatures) is iterated
//! a few times per sweep from the last values (deferred), and converges
//! with the outer iterations.
//!
//! **Orthotropic conductivity.** The face conductivity of the matrix is
//! k_n of the face normal; the rest of the flux, −(K S − k_n S)·∇T_f
//! (tangential to the face), is a deferred correction with the
//! interpolated gradient, on interior faces of a region and on its
//! boundary faces with a fixed temperature (a convection boundary with the
//! weight (1 − a) of its face value; a flux boundary carries its given
//! flux).
//!
//! All this is a correction of the face constants b_f of the device's face
//! assembly (one number per face: added to P, subtracted from N, so heat
//! flux continuity holds to round-off at every iteration, D117 (1)), and of
//! the gradient of the interface cells, both on the host; the device
//! kernels stay those of every scalar (CPU and GPU alike).

use crate::device::{Buf, Buf3, FlowDevice};
use anyhow::Result;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// `[interfaces.<name>]`: a contact resistance or a thin wall between two
/// regions (stage 7d2, D157 §9.1).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InterfaceSection {
    /// The two regions: solid regions of `[regions]`, a fluid region, or
    /// `fluid` (all fluid cells).
    pub regions: [String; 2],
    /// The thermal contact resistance R″_c (m² K/W, ≥ 0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub contact_resistance: Option<f64>,
    /// A thin wall in the face, as the resistance t/k (no heat capacity;
    /// the full thin wall with layers is D189).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thin_wall: Option<ThinWallSection>,
}

/// `thin_wall = { thickness, conductivity | material }`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ThinWallSection {
    /// t (m, > 0).
    pub thickness: f64,
    /// k (W/(m K), > 0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<f64>,
    /// A material with one conductivity (k at the face temperature).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
}

/// The resistance of an `[interfaces]` entry.
#[derive(Debug, Clone)]
pub(crate) enum Resistance {
    /// R″_c (m² K/W).
    Fixed(f64),
    /// t and k(T) of a thin wall.
    Wall(f64, crate::materials::Curve),
}

impl Resistance {
    /// R″ at the face temperature `t`.
    pub fn at(&self, t: f64) -> f64 {
        match self {
            Resistance::Fixed(r) => *r,
            Resistance::Wall(th, k) => th / k.at(t),
        }
    }
}

/// A pair of regions with interface faces, for the report.
#[derive(Debug, Clone)]
pub(crate) struct Pair {
    /// The name in the report: `interface-<solid>` between the fluid and a
    /// solid, `interface-<a>-<b>` between two solids.
    pub name: String,
}

/// One interface face.
#[derive(Debug, Clone)]
pub(crate) struct IFace {
    /// The face, its owner and neighbour.
    pub f: u32,
    pub o: u32,
    pub n: u32,
    /// |S| and the unit normal.
    pub area: f64,
    pub normal: [f64; 3],
    /// r_P − δ_P n̂ and r_N + δ_N n̂ (the tangential offsets).
    pub tp: [f64; 3],
    pub tn: [f64; 3],
    /// The non-orthogonal vector of the face assembly, k = S − d |S|²/(d·S).
    pub kvec: [f64; 3],
    /// |S|²/(d·S) and the owner weight of the gradient interpolation.
    pub dgeo: f64,
    pub w: f64,
    /// The `[interfaces]` entry (index), if any.
    pub spec: Option<usize>,
    /// The report pair and whether the owner is on the side the heat flow
    /// is counted into (the fluid; of two solids the first in the regions'
    /// order).
    pub pair: usize,
    pub owner_into: bool,
    /// Refreshed per outer iteration: R_P, R″_c, R_N (m² K/W) and the
    /// cross vectors K n̂ − k_n n̂ of the two sides.
    pub r_p: f64,
    pub r_c: f64,
    pub r_n: f64,
    pub a_p: [f64; 3],
    pub a_n: [f64; 3],
}

/// A point of the one-sided least-squares stencil of a cell.
#[derive(Debug, Clone, Copy)]
enum Point {
    /// A cell of the same region.
    Cell(u32),
    /// The face temperature of the own side of an interface face (slot,
    /// whether the cell is the owner).
    Face(u32, bool),
    /// A boundary face (a φ_P + b).
    Boundary(u32),
}

/// A cell with a one-sided gradient: the stencil points and their weight
/// vectors, ∇T = Σ v_j (T_j − T_P).
#[derive(Debug, Clone)]
struct OneCell {
    c: u32,
    pts: Vec<(Point, [f64; 3])>,
}

/// The interface and orthotropic conduction of the energy equation.
#[derive(Debug, Clone)]
pub(crate) struct Conduction {
    /// The interface faces.
    pub faces: Vec<IFace>,
    /// Per interior face: the slot in `faces` (`u32::MAX`: none).
    pub slot: Vec<u32>,
    /// The report pairs.
    pub pairs: Vec<Pair>,
    /// The resistances of the `[interfaces]` entries.
    pub specs: Vec<Resistance>,
    cells: Vec<OneCell>,
    /// Per cell: its entry of `cells` (`u32::MAX`: none).
    cell_slot: Vec<u32>,
    /// The one-sided gradients of the last sweep (per entry of `cells`).
    grad: Vec<[f64; 3]>,
    /// The faces (interior or boundary) inside orthotropic regions with
    /// their cross vectors K S − k_n S, refreshed per outer iteration.
    pub cross: Vec<CrossFace>,
    /// Per interface face of the last correction: q |S| from P to N (W)
    /// and the face temperatures of the owner's and the neighbour's side.
    pub heat: Vec<f64>,
    pub t_o: Vec<f64>,
    pub t_n: Vec<f64>,
    /// Per boundary face: the heat flow into the domain that the cross
    /// terms add (W), of the last correction.
    pub boundary_cross: Vec<f64>,
    /// Host copies of the last sweep (the correction after the face
    /// assembly reads them).
    phi: Vec<f64>,
    g: [Vec<f64>; 3],
    bc: (Vec<f64>, Vec<f64>),
    n_interior: usize,
    /// The Robin conditions of the solid phase of the cycle-averaged loop
    /// (D204), per interface face: `None` keeps the face as it is.
    pub robin: Option<Vec<Option<RobinFace>>>,
}

/// A fluid–solid interface face of the solid phase of the cycle-averaged
/// loop (D204): the heat flow into the solid is Ḡ (T̄_ref − T_s*), T_s* the
/// solid cell's temperature at the face (the offset of the face and the
/// cross term of an orthotropic solid included), Ḡ the cycle mean of the
/// conductance of the fluid side, the contact and the solid half cell in
/// series. The face itself carries no conduction in the matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RobinFace {
    /// Ḡ (W/K).
    pub g: f64,
    /// T̄_ref (K).
    pub t_ref: f64,
    /// The solid cell is the owner of the face.
    pub solid_owner: bool,
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn axpy(a: f64, x: [f64; 3], y: [f64; 3]) -> [f64; 3] {
    [a * x[0] + y[0], a * x[1] + y[1], a * x[2] + y[2]]
}

fn inverse3(m: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    let scale: f64 = m.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
    if !det.is_finite() || det.abs() <= 1e-14 * scale * scale * scale {
        return None;
    }
    let inv = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv,
        ],
    ])
}

impl Conduction {
    /// The interface faces of `mesh` between cells of different regions
    /// (`region` per cell, an index of `names`), with the `[interfaces]`
    /// entry of each face (`spec_of(f)`), the report pairs, and the
    /// one-sided stencils of their cells; `fixed` per boundary face: true
    /// where the face value is a fixed temperature or a convection
    /// condition (the stencil point is the face centroid, else its
    /// projection on the normal, as in the device's least squares).
    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    pub fn new(
        mesh: &Mesh,
        geo: &Geometry,
        disc: &foehncraft_fv::disc::Discretisation,
        region: &[u32],
        names: &[String],
        spec_of: &dyn Fn(usize) -> Option<usize>,
        specs: Vec<Resistance>,
        fixed: &[bool],
    ) -> Result<Self> {
        let ni = mesh.n_interior_faces();
        let nc = mesh.n_cells;
        let mut faces = Vec::new();
        let mut slot = vec![u32::MAX; ni];
        let mut pairs: Vec<(u32, u32, Pair)> = Vec::new();
        let xc = &geo.cell_centroid;
        for f in 0..ni {
            let (o, n) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
            let (ro, rn) = (region[o], region[n]);
            if ro == rn {
                continue;
            }
            let (lo, hi) = (ro.min(rn), ro.max(rn));
            let pair = match pairs.iter().position(|p| p.0 == lo && p.1 == hi) {
                Some(i) => i,
                None => {
                    let name = if lo == 0 {
                        crate::cht::interface_name(&names[hi as usize])
                    } else {
                        format!(
                            "{}{}-{}",
                            crate::cht::INTERFACE,
                            names[lo as usize],
                            names[hi as usize]
                        )
                    };
                    pairs.push((lo, hi, Pair { name }));
                    pairs.len() - 1
                }
            };
            let (nrm, area) = geo.normal_area(f);
            let xf = geo.face_centroid[f];
            let rp = sub(xf, xc[o]);
            let rn_ = sub(xf, xc[n]);
            let dp = dot(rp, nrm);
            let dn = -dot(rn_, nrm);
            slot[f] = faces.len() as u32;
            faces.push(IFace {
                f: f as u32,
                o: o as u32,
                n: n as u32,
                area,
                normal: nrm,
                tp: axpy(-dp, nrm, rp),
                tn: axpy(dn, nrm, rn_),
                kvec: [disc.k.x[f], disc.k.y[f], disc.k.z[f]],
                dgeo: disc.dgeo[f],
                w: disc.weight[f],
                spec: spec_of(f),
                pair,
                owner_into: ro == lo,
                r_p: 0.0,
                r_c: 0.0,
                r_n: 0.0,
                a_p: [0.0; 3],
                a_n: [0.0; 3],
            });
        }
        // the one-sided stencils of the cells next to an interface
        let mut is_iface_cell = vec![false; nc];
        for fa in &faces {
            is_iface_cell[fa.o as usize] = true;
            is_iface_cell[fa.n as usize] = true;
        }
        let mut cells = Vec::new();
        for c in 0..nc {
            if !is_iface_cell[c] {
                continue;
            }
            let mut pts: Vec<(Point, [f64; 3], f64)> = Vec::new();
            for (f, is_nb) in mesh.cell_face_iter(c) {
                if f < ni {
                    let other = if is_nb {
                        mesh.owner[f] as usize
                    } else {
                        mesh.neighbour[f] as usize
                    };
                    if slot[f] != u32::MAX {
                        let d = sub(geo.face_centroid[f], xc[c]);
                        pts.push((Point::Face(slot[f], !is_nb), d, 1.0 / dot(d, d)));
                    } else {
                        let d = sub(xc[other], xc[c]);
                        pts.push((Point::Cell(other as u32), d, 1.0 / dot(d, d)));
                    }
                } else {
                    let d = sub(geo.face_centroid[f], xc[c]);
                    let d = if fixed[f - ni] {
                        d
                    } else {
                        let (nrm, _) = geo.normal_area(f);
                        let s = dot(d, nrm);
                        [s * nrm[0], s * nrm[1], s * nrm[2]]
                    };
                    pts.push((Point::Boundary((f - ni) as u32), d, 1.0 / dot(d, d)));
                }
            }
            let mut m = [[0.0; 3]; 3];
            for (_, d, w) in &pts {
                for i in 0..3 {
                    for j in 0..3 {
                        m[i][j] += w * d[i] * d[j];
                    }
                }
            }
            let Some(inv) = inverse3(m) else {
                anyhow::bail!(
                    "cell {c} next to a region interface: its one-sided least-squares stencil \
                     is degenerate"
                );
            };
            let pts = pts
                .into_iter()
                .map(|(p, d, w)| {
                    let mut v = [0.0; 3];
                    for (r, vr) in v.iter_mut().enumerate() {
                        *vr = w * (inv[r][0] * d[0] + inv[r][1] * d[1] + inv[r][2] * d[2]);
                    }
                    (p, v)
                })
                .collect();
            cells.push(OneCell { c: c as u32, pts });
        }
        let nf = faces.len();
        let mut cell_slot = vec![u32::MAX; nc];
        for (i, oc) in cells.iter().enumerate() {
            cell_slot[oc.c as usize] = i as u32;
        }
        Ok(Self {
            cell_slot,
            faces,
            slot,
            pairs: pairs.into_iter().map(|p| p.2).collect(),
            specs,
            grad: Vec::new(),
            cells,
            cross: Vec::new(),
            heat: vec![0.0; nf],
            t_o: vec![0.0; nf],
            t_n: vec![0.0; nf],
            boundary_cross: vec![0.0; mesh.n_boundary_faces()],
            phi: Vec::new(),
            g: [Vec::new(), Vec::new(), Vec::new()],
            bc: (Vec::new(), Vec::new()),
            n_interior: ni,
            robin: None,
        })
    }

    /// True when the energy equation needs the corrections (interfaces or
    /// orthotropic faces).
    pub fn is_active(&self) -> bool {
        !self.faces.is_empty() || !self.cross.is_empty()
    }

    /// The flux q |S| (W, from P to N) and the face temperatures of the two
    /// sides of interface face `i`, with the gradients `gp`, `gn`.
    #[inline]
    fn face(&self, i: usize, gp: [f64; 3], gn: [f64; 3]) -> (f64, f64, f64) {
        let fa = &self.faces[i];
        let tp = self.phi[fa.o as usize] + dot(gp, fa.tp);
        let tn = self.phi[fa.n as usize] + dot(gn, fa.tn);
        let cp = -dot(fa.a_p, gp);
        let cn = -dot(fa.a_n, gn);
        if let Some(r) = self.robin.as_ref().and_then(|r| r[i]) {
            // the solid phase of the cycle loop: the Robin pair on the
            // solid side, the face temperature of the solid side through
            // its half cell; the other side shows the same face value
            let (ts, rs) = if r.solid_owner {
                (tp + cp * fa.r_p, fa.r_p)
            } else {
                (tn - cn * fa.r_n, fa.r_n)
            };
            let q_in = r.g * (r.t_ref - ts);
            let tw = ts + rs * q_in / fa.area;
            let q_pn = if r.solid_owner { -q_in } else { q_in };
            return (q_pn, tw, tw);
        }
        let rt = fa.r_p + fa.r_c + fa.r_n;
        let q = (tp - tn + cp * fa.r_p + cn * fa.r_n) / rt;
        (q * fa.area, tp - fa.r_p * (q - cp), tn + fa.r_n * (q - cn))
    }

    /// The gradient of cell `c` (one-sided where it has a stencil).
    #[inline]
    fn grad_of(&self, cell_slot: &[u32], c: usize) -> [f64; 3] {
        match cell_slot.get(c) {
            Some(&s) if s != u32::MAX => self.grad[s as usize],
            _ => [self.g[0][c], self.g[1][c], self.g[2][c]],
        }
    }

    /// Replaces the device gradient `grad` of T in the interface cells by
    /// the one-sided gradient (§9.2), iterated `PASSES` times with the face
    /// temperatures from the last sweep's gradients.
    pub fn one_sided(
        &mut self,
        dev: &mut dyn FlowDevice,
        phi: Buf,
        grad: Buf3,
        bc_a: Buf,
        bc_b: Buf,
    ) -> Result<()> {
        const PASSES: usize = 3;
        self.phi = dev.download(phi)?;
        self.g = [
            dev.download(grad.x)?,
            dev.download(grad.y)?,
            dev.download(grad.z)?,
        ];
        self.bc = (dev.download(bc_a)?, dev.download(bc_b)?);
        if self.cells.is_empty() {
            return Ok(());
        }
        let cell_slot = std::mem::take(&mut self.cell_slot);
        if self.grad.len() != self.cells.len() {
            self.grad = self
                .cells
                .iter()
                .map(|oc| {
                    let c = oc.c as usize;
                    [self.g[0][c], self.g[1][c], self.g[2][c]]
                })
                .collect();
        }
        let mut tf = vec![(0.0, 0.0); self.faces.len()];
        for _ in 0..PASSES {
            for (i, t) in tf.iter_mut().enumerate() {
                let fa = &self.faces[i];
                let gp = self.grad_of(&cell_slot, fa.o as usize);
                let gn = self.grad_of(&cell_slot, fa.n as usize);
                let (_, a, b) = self.face(i, gp, gn);
                *t = (a, b);
            }
            let mut new = Vec::with_capacity(self.cells.len());
            for oc in &self.cells {
                let c = oc.c as usize;
                let pc = self.phi[c];
                let mut g = [0.0; 3];
                for &(p, v) in &oc.pts {
                    let val = match p {
                        Point::Cell(n) => self.phi[n as usize],
                        Point::Face(s, own) => {
                            let t = tf[s as usize];
                            if own { t.0 } else { t.1 }
                        }
                        Point::Boundary(fb) => {
                            let fb = fb as usize;
                            self.bc.0[fb] * pc + self.bc.1[fb]
                        }
                    };
                    g = axpy(val - pc, v, g);
                }
                new.push(g);
            }
            self.grad = new;
        }
        self.cell_slot = cell_slot;
        for (oc, g) in self.cells.iter().zip(&self.grad) {
            let c = oc.c as usize;
            for (d, v) in g.iter().enumerate() {
                self.g[d][c] = *v;
            }
        }
        dev.upload(grad.x, &self.g[0])?;
        dev.upload(grad.y, &self.g[1])?;
        dev.upload(grad.z, &self.g[2])?;
        Ok(())
    }

    /// Adds the corrections to the face constants `bf` of the face assembly
    /// (after it, with the gradients of [`Conduction::one_sided`]): the
    /// exact interface flux in place of the assembly's two-sided
    /// non-orthogonal part, and the cross terms of orthotropic regions.
    pub fn correct(&mut self, dev: &mut dyn FlowDevice, gamma: Buf, bf: Buf) -> Result<()> {
        let gam = dev.download(gamma)?;
        let mut b = dev.download(bf)?;
        let ni = self.n_interior;
        let cs = std::mem::take(&mut self.cell_slot);
        for i in 0..self.faces.len() {
            let fa = &self.faces[i];
            let (f, o, n, w) = (fa.f as usize, fa.o as usize, fa.n as usize, fa.w);
            let gp = self.grad_of(&cs, o);
            let gn = self.grad_of(&cs, n);
            let gf = [
                w * gp[0] + (1.0 - w) * gn[0],
                w * gp[1] + (1.0 - w) * gn[1],
                w * gp[2] + (1.0 - w) * gn[2],
            ];
            let g_f = gam[f] * fa.dgeo;
            let cp = -dot(fa.a_p, gp);
            let cn = -dot(fa.a_n, gn);
            // the exact flux minus the implicit part G (T_P − T_N), in
            // place of the assembly's −Γ k·∇T_f (G = Γ |S|²/(d·S) =
            // |S|/(R_P + R″_c + R_N) by the face conductivity)
            b[f] += g_f * (dot(gp, fa.tp) - dot(gn, fa.tn) + cp * fa.r_p + cn * fa.r_n)
                + gam[f] * dot(fa.kvec, gf);
            let (q, a, bb) = self.face(i, gp, gn);
            self.heat[i] = q;
            self.t_o[i] = a;
            self.t_n[i] = bb;
        }
        for v in self.boundary_cross.iter_mut() {
            *v = 0.0;
        }
        for cf in &self.cross {
            let f = cf.f as usize;
            let go = self.grad_of(&cs, cf.o as usize);
            if f < ni {
                // an interior face of one orthotropic region: −t·∇T_f with
                // the interpolated gradient
                let gn = self.grad_of(&cs, cf.n as usize);
                let gf = [
                    cf.w * go[0] + (1.0 - cf.w) * gn[0],
                    cf.w * go[1] + (1.0 - cf.w) * gn[1],
                    cf.w * go[2] + (1.0 - cf.w) * gn[2],
                ];
                b[f] -= dot(cf.t, gf);
            } else {
                // a boundary face: with the weight (1 − a) of its face
                // value, as the assembly's non-orthogonal part
                let fb = f - ni;
                let x = (1.0 - self.bc.0[fb]) * dot(cf.t, go);
                b[f] -= x;
                self.boundary_cross[fb] = x;
            }
        }
        self.cell_slot = cs;
        dev.upload(bf, &b)?;
        Ok(())
    }

    /// The heat flows and face temperatures of the interface faces from
    /// the temperature `phi` after the solve (with the gradients of the
    /// last sweep), so that the report is not one solve behind.
    pub fn settle(&mut self, dev: &mut dyn FlowDevice, phi: Buf) -> Result<()> {
        if self.faces.is_empty() {
            return Ok(());
        }
        self.phi = dev.download(phi)?;
        let cs = std::mem::take(&mut self.cell_slot);
        for i in 0..self.faces.len() {
            let fa = &self.faces[i];
            let gp = self.grad_of(&cs, fa.o as usize);
            let gn = self.grad_of(&cs, fa.n as usize);
            let (q, a, b) = self.face(i, gp, gn);
            self.heat[i] = q;
            self.t_o[i] = a;
            self.t_n[i] = b;
        }
        self.cell_slot = cs;
        Ok(())
    }

    /// The Robin sources of the solid phase (D204) per solid cell of a
    /// Robin face: (cell, the explicit heat flow Ḡ (T̄_ref − δT) (W) with
    /// δT the offset and cross term of the face from the last sweep's
    /// gradient, Ḡ (W/K) on the diagonal).
    pub fn robin_sources(&self) -> Vec<(u32, f64, f64)> {
        let Some(rb) = self.robin.as_ref() else {
            return Vec::new();
        };
        let cs = &self.cell_slot;
        let mut out = Vec::new();
        for (fa, r) in self.faces.iter().zip(rb) {
            let Some(r) = r else { continue };
            // T_s* − T_P: the offset and the cross term (−a·∇T r on the
            // owner's side, +a·∇T r on the neighbour's, as in `face`)
            let (c, off, a, rr) = if r.solid_owner {
                (fa.o, fa.tp, fa.a_p, -fa.r_p)
            } else {
                (fa.n, fa.tn, fa.a_n, fa.r_n)
            };
            let g = if self.phi.is_empty() {
                [0.0; 3]
            } else {
                self.grad_of(cs, c as usize)
            };
            let dt = dot(g, off) + dot(a, g) * rr;
            out.push((c, r.g * (r.t_ref - dt), r.g));
        }
        out
    }

    /// Per interface face of the last solve: (area, the conductance
    /// |S|/(R_P + R″_c + R_N) (W/K), the conductance of the fluid side
    /// |S|/R_fluid (W/K), the heat flow into the solid (W), the solid
    /// cell's temperature at the face T_s* (K), the face temperature on
    /// the solid side (K)), for the fluid–solid faces with `Some`; `fluid`
    /// tells per face whether the owner is the fluid (`None`: two solids).
    pub fn face_states(&self, fluid_owner: &[Option<bool>]) -> Vec<Option<[f64; 6]>> {
        self.faces
            .iter()
            .enumerate()
            .zip(fluid_owner)
            .map(|((i, fa), fo)| {
                let fo = (*fo)?;
                let rt = fa.r_p + fa.r_c + fa.r_n;
                let (q_in, tw, rs, rf) = if fo {
                    (self.heat[i], self.t_n[i], fa.r_n, fa.r_p)
                } else {
                    (-self.heat[i], self.t_o[i], fa.r_p, fa.r_n)
                };
                let ts = tw - rs * q_in / fa.area;
                Some([fa.area, fa.area / rt, fa.area / rf, q_in, ts, tw])
            })
            .collect()
    }

    /// The face temperatures and heat flows of the last solve, per
    /// report pair: (name, heat flow into the pair's `into` region (W),
    /// area, area-mean face temperature on that side).
    pub fn report(&self) -> Vec<(String, f64, f64, f64)> {
        let mut out: Vec<(String, f64, f64, f64)> = self
            .pairs
            .iter()
            .map(|p| (p.name.clone(), 0.0, 0.0, 0.0))
            .collect();
        for (i, fa) in self.faces.iter().enumerate() {
            let r = &mut out[fa.pair];
            // q from P to N: into the owner's side is −q
            let (q, t) = if fa.owner_into {
                (-self.heat[i], self.t_o[i])
            } else {
                (self.heat[i], self.t_n[i])
            };
            r.1 += q;
            r.2 += fa.area;
            r.3 += t * fa.area;
        }
        for r in &mut out {
            if r.2 > 0.0 {
                r.3 /= r.2;
            }
        }
        out
    }
}

/// A face inside an orthotropic region with its cross vector.
#[derive(Debug, Clone)]
pub(crate) struct CrossFace {
    /// The face, its owner and (interior faces) neighbour.
    pub f: u32,
    pub o: u32,
    pub n: u32,
    /// The owner weight of the gradient interpolation.
    pub w: f64,
    /// t = K S − (n̂·K·n̂) S (W/K per unit gradient).
    pub t: [f64; 3],
}
