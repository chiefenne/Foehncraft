//! The two temperatures of porous regions (D271, porous plan stage PM2):
//! the gas and the matrix of a porous cell each have their own temperature
//! (local thermal non-equilibrium; Schumann 1929; Kaviany 1995, ch. 3;
//! Nield & Bejan 2017, sec. 2.2.2),
//!
//!   gas:    ε ρ_g c_g ∂T_g/∂t + ∇·(Ċ T_g) − T_g ∇·Ċ = ∇·(ε k_g ∇T_g) + h a (T_s − T_g)
//!   matrix: (1 − ε) ρ_s ∂h_s/∂t = ∇·(K_s ∇T_s) + h a (T_g − T_s).
//!
//! The gas temperature is the temperature of the energy equation, one field
//! with the clear gas around the region; its region has a material of its
//! own with ε ρ_g and ε k_g ([`crate::scalar`]). The matrix temperature adds
//! one row per porous cell to the same matrix, coupled to the gas row of
//! its cell by h a V and to the matrix rows of its porous neighbours by the
//! conduction of the effective tensor K_s (D272), so gas and matrix are
//! solved in one linear system every outer iteration (monolithic, as the
//! solids of D158): the exchange is exact within the iteration however
//! large h a is. The device assembles the energy rows as for any case;
//! this module takes the matrix as CSR, adds the matrix rows and solves
//! the extended system with BiCGStab and Jacobi on the host (one process,
//! serial in this stage; ranks and the GPU in PM5). The storage of the
//! matrix is conservative in its enthalpy h_s(T) (D150): linearised at the
//! current T_s, exact at convergence.
//!
//! The matrix is adiabatic at the faces of the region (to the clear gas, to
//! a solid region, at the boundaries): it takes and gives heat through the
//! gas only. Its contact with a housing (D274) is a later step.

use crate::device::{Buf, FlowDevice};
use crate::energy::{Enthalpy, Material};
use crate::materials::Curve;
use crate::porous::{SolidK, Transfer};
use anyhow::Result;
use foehncraft_fv::disc::Discretisation;

/// One porous region with two temperatures in the report of a solve.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct PorousHeat {
    /// The region.
    pub region: String,
    /// Volume mean of the gas temperature over the region (K).
    pub gas_mean: f64,
    /// Volume mean of the matrix temperature (K).
    pub matrix_mean: f64,
    /// Lowest and highest matrix temperature (K).
    pub matrix_min: f64,
    /// See `matrix_min`.
    pub matrix_max: f64,
    /// The sensible heat held by the matrix, ∫ (1 − ε) ρ_s h_s(T_s) dV
    /// (J, h_s = 0 at 298.15 K).
    pub matrix_heat: f64,
    /// The heat from the gas into the matrix, Σ h a V (T_g − T_s) (W).
    pub exchange: f64,
}

/// The matrix of one porous region.
#[derive(Debug, Clone)]
struct Region {
    name: String,
    /// The material of the matrix (ρ_s, c_s(T), k_s(T)).
    solid: Material,
    /// h_s(T) of the matrix.
    h: Enthalpy,
    /// ε.
    eps: f64,
    /// a (1/m).
    area: f64,
    transfer: Transfer,
    /// k_g(T) of the gas (the Nusselt form of h).
    k_gas: Curve,
    /// e₁, e₂, e₃ of the region.
    axes: [[f64; 3]; 3],
    k: SolidK,
}

impl Region {
    /// h (W/(m² K)) at the gas temperature `tg`.
    fn h_at(&self, tg: f64) -> f64 {
        match self.transfer {
            Transfer::Fixed(h) => h,
            Transfer::Nusselt { nu, dh } => nu * self.k_gas.at(tg) / dh,
        }
    }

    /// The effective conductivity of the matrix along the unit normal `n`
    /// at `t` (W/(m K)).
    fn k_normal(&self, t: f64, n: [f64; 3]) -> f64 {
        let k = match self.k {
            SolidK::Values(v) => v,
            SolidK::Factors(f) => {
                let ks = self.solid.k.at(t);
                f.map(|x| x * ks)
            }
        };
        (0..3)
            .map(|i| {
                let c = (0..3).map(|d| self.axes[i][d] * n[d]).sum::<f64>();
                k[i] * c * c
            })
            .sum()
    }
}

/// A face between two matrix rows.
#[derive(Debug, Clone, Copy)]
struct Face {
    rows: [u32; 2],
    /// |S|²/(d·S).
    dgeo: f64,
    /// Unit normal.
    n: [f64; 3],
    /// Owner weight of the face interpolation.
    w: f64,
}

/// The name of the matrix temperature in the outputs and the restart file.
pub const FIELD: &str = "matrix_temperature";

/// The matrix temperatures of the porous regions of an energy mesh.
#[derive(Debug, Clone)]
pub struct Matrix {
    regions: Vec<Region>,
    /// The energy cell of each matrix row.
    cell: Vec<u32>,
    /// The region of each row.
    region: Vec<u32>,
    /// The volume of each row's cell.
    vol: Vec<f64>,
    faces: Vec<Face>,
    /// The CSR pattern of the energy matrix in its sliced-ELL storage:
    /// per cell `ptr[c]..ptr[c+1]` into (`pos`, `col`).
    ptr: Vec<usize>,
    pos: Vec<u32>,
    col: Vec<u32>,
    /// T_s per row (K).
    pub ts: Vec<f64>,
    /// The old levels n, n−1 and the one before (transient).
    old: Option<[Vec<f64>; 3]>,
    /// The last valid values (D63).
    valid: Option<Vec<f64>>,
    /// h a V (T_g − T_s) per row of the last solve (W): the heat into the
    /// matrix.
    pub exchange: Vec<f64>,
}

/// One porous region with its cells, for [`Matrix::new`].
pub struct RegionSpec {
    /// The region.
    pub name: String,
    /// The material of its matrix.
    pub solid: Material,
    /// Whether c_s follows T.
    pub dependent: bool,
    /// The resolved thermal values.
    pub thermal: crate::porous::Thermal,
    /// k_g(T) of its gas (unscaled).
    pub k_gas: Curve,
    /// Its frame.
    pub axes: [[f64; 3]; 3],
    /// Its cells of the energy mesh.
    pub cells: Vec<u32>,
}

impl Matrix {
    /// The matrix rows of `specs` on the energy discretisation `disc`, with
    /// the start temperatures `t0` of the energy cells.
    pub fn new(
        disc: &Discretisation,
        specs: Vec<RegionSpec>,
        t0: &[f64],
        transient: bool,
    ) -> Result<Self> {
        let nc = disc.n_cells;
        let mut row = vec![u32::MAX; nc];
        let (mut cell, mut region, mut vol) = (Vec::new(), Vec::new(), Vec::new());
        let mut regions = Vec::with_capacity(specs.len());
        for (ri, s) in specs.into_iter().enumerate() {
            for &c in &s.cells {
                anyhow::ensure!(
                    row[c as usize] == u32::MAX,
                    "porous region `{}`: cell {c} is in two porous regions",
                    s.name
                );
                row[c as usize] = cell.len() as u32;
                cell.push(c);
                region.push(ri as u32);
                vol.push(disc.cell_volume[c as usize]);
            }
            let h = if s.dependent {
                Enthalpy::of(&s.solid.cp)
            } else {
                Enthalpy::constant(s.solid.cp.at(s.solid.t_ref))
            };
            regions.push(Region {
                name: s.name,
                h,
                eps: s.thermal.porosity,
                area: s.thermal.area,
                transfer: s.thermal.transfer,
                k_gas: s.k_gas,
                axes: s.axes,
                k: s.thermal.solid_k,
                solid: s.solid,
            });
        }
        let mut faces = Vec::new();
        for f in 0..disc.n_interior {
            let (o, n) = (disc.owner[f] as usize, disc.neighbour[f] as usize);
            if row[o] == u32::MAX || row[n] == u32::MAX {
                continue;
            }
            let s = [
                disc.face_area.x[f],
                disc.face_area.y[f],
                disc.face_area.z[f],
            ];
            let a = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
            faces.push(Face {
                rows: [row[o], row[n]],
                dgeo: disc.dgeo[f],
                n: s.map(|x| x / a),
                w: disc.weight[f],
            });
        }
        // the storage positions of the CSR rows (as Discretisation::matrix_csr)
        let sell = &disc.sell;
        let mut ptr = Vec::with_capacity(nc + 1);
        ptr.push(0);
        let (mut pos, mut col) = (Vec::new(), Vec::new());
        for c in 0..nc {
            let lo = disc.cell_ptr[c] as usize;
            let hi = disc.cell_ptr[c + 1] as usize;
            let len = 1 + disc.cell_faces[lo..hi]
                .iter()
                .filter(|&&e| foehncraft_kernels::unpack_face(e).0 < disc.n_interior)
                .count();
            for j in 0..len {
                let k = sell.pos(c, j);
                pos.push(k as u32);
                col.push(sell.col[k]);
            }
            ptr.push(pos.len());
        }
        let ts: Vec<f64> = cell.iter().map(|&c| t0[c as usize]).collect();
        let n = ts.len();
        Ok(Self {
            regions,
            cell,
            region,
            vol,
            faces,
            ptr,
            pos,
            col,
            old: transient.then(|| [ts.clone(), ts.clone(), ts.clone()]),
            valid: None,
            exchange: vec![0.0; n],
            ts,
        })
    }

    /// The number of matrix rows.
    pub fn len(&self) -> usize {
        self.ts.len()
    }

    /// True without matrix rows.
    pub fn is_empty(&self) -> bool {
        self.ts.is_empty()
    }

    /// The region names.
    pub fn names(&self) -> Vec<String> {
        self.regions.iter().map(|r| r.name.clone()).collect()
    }

    /// The matrix temperature on the energy cells: T_s in the porous
    /// cells, `t` (the temperature of the energy equation) elsewhere.
    pub fn field(&self, t: &[f64]) -> Vec<f64> {
        let mut v = t.to_vec();
        for (r, &c) in self.cell.iter().enumerate() {
            v[c as usize] = self.ts[r];
        }
        v
    }

    /// The report per region, with the gas temperature `tg` of the energy
    /// cells.
    pub fn report(&self, tg: &[f64]) -> Vec<PorousHeat> {
        let mut out: Vec<PorousHeat> = self
            .regions
            .iter()
            .map(|r| PorousHeat {
                region: r.name.clone(),
                matrix_min: f64::INFINITY,
                matrix_max: f64::NEG_INFINITY,
                ..PorousHeat::default()
            })
            .collect();
        let mut vol = vec![0.0; self.regions.len()];
        for r in 0..self.len() {
            let i = self.region[r] as usize;
            let g = &self.regions[i];
            let v = self.vol[r];
            let o = &mut out[i];
            vol[i] += v;
            o.gas_mean += v * tg[self.cell[r] as usize];
            o.matrix_mean += v * self.ts[r];
            o.matrix_min = o.matrix_min.min(self.ts[r]);
            o.matrix_max = o.matrix_max.max(self.ts[r]);
            o.matrix_heat += (1.0 - g.eps) * g.solid.rho * v * g.h.at(self.ts[r]);
            o.exchange += self.exchange[r];
        }
        for (o, v) in out.iter_mut().zip(vol) {
            if v > 0.0 {
                o.gas_mean /= v;
                o.matrix_mean /= v;
            }
        }
        out
    }

    /// Per region: the heat held by the matrix relative to `t_ref`,
    /// ∫ (1 − ε) ρ_s (h_s(T_s) − h_s(t_ref)) dV (J), and the heat into it
    /// over the last solve, Σ h a V (T_g − T_s) (W).
    pub fn heat(&self, t_ref: f64) -> Vec<(String, f64, f64)> {
        let mut out: Vec<(String, f64, f64)> = self
            .regions
            .iter()
            .map(|r| (r.name.clone(), 0.0, 0.0))
            .collect();
        for r in 0..self.len() {
            let g = &self.regions[self.region[r] as usize];
            let m = (1.0 - g.eps) * g.solid.rho * self.vol[r];
            let o = &mut out[self.region[r] as usize];
            o.1 += m * (g.h.at(self.ts[r]) - g.h.at(t_ref));
            o.2 += self.exchange[r];
        }
        out
    }

    /// Per region: the heat stored by the matrix in the last time step
    /// (W), in the form of the time scheme: Σ (1 − ε) ρ_s V (a₀ h_s(T_s) +
    /// a₁ h_s(T_s,n) + a₂ h_s(T_s,n−1)) / Δt.
    pub fn storage(&self, a: [f64; 3], dt: f64) -> Vec<f64> {
        let mut out = vec![0.0; self.regions.len()];
        let Some(old) = &self.old else {
            return out;
        };
        for r in 0..self.len() {
            let g = &self.regions[self.region[r] as usize];
            let mut l = a[0] * g.h.at(self.ts[r]) + a[1] * g.h.at(old[0][r]);
            if a[2] != 0.0 {
                l += a[2] * g.h.at(old[1][r]);
            }
            out[self.region[r] as usize] += (1.0 - g.eps) * g.solid.rho * self.vol[r] * l / dt;
        }
        out
    }

    /// The label of a region's matrix in the heat reports.
    pub fn label(region: &str) -> String {
        format!("{region} matrix")
    }

    /// Start of a time step: T_s,n−1 ← T_s,n, T_s,n ← T_s (`has_old`), and
    /// the start value extrapolated with `extrapolate`, as the energy field.
    pub fn step_start(&mut self, has_old: bool, extrapolate: Option<f64>) {
        let Some(old) = &mut self.old else {
            return;
        };
        old[2] = old[1].clone();
        if has_old {
            old[1] = old[0].clone();
        }
        old[0] = self.ts.clone();
        if has_old && let Some(w) = extrapolate {
            for ((t, a), b) in self.ts.iter_mut().zip(&old[0]).zip(&old[1]) {
                *t = (1.0 + w) * a - w * b;
            }
        }
    }

    /// Takes the last step back exactly.
    pub fn step_back(&mut self) {
        if let Some(old) = &mut self.old {
            self.ts = old[0].clone();
            old[0] = old[1].clone();
            old[1] = old[2].clone();
        }
    }

    /// Keeps the values as the last valid state (D63).
    pub fn keep_valid(&mut self) {
        self.valid = Some(self.ts.clone());
    }

    /// Puts the last valid values back.
    pub fn rewind(&mut self) {
        if let Some(v) = &self.valid {
            self.ts = v.clone();
        }
    }

    /// After a rewind of a time step: the level n−1 back into level n.
    pub fn rewind_levels(&mut self) {
        if let Some(old) = &mut self.old {
            old[0] = old[1].clone();
        }
    }

    /// The level before the current one on the energy cells (restart: the
    /// level n−1 of the next step), as [`Matrix::field`].
    pub fn old_field(&self, t: &[f64]) -> Option<Vec<f64>> {
        let old = self.old.as_ref()?;
        let mut v = t.to_vec();
        for (r, &c) in self.cell.iter().enumerate() {
            v[c as usize] = old[0][r];
        }
        Some(v)
    }

    /// Sets T_s (and the level before) from fields on the energy cells
    /// (restart).
    pub fn load(&mut self, field: &[f64], old: Option<&[f64]>) -> Result<()> {
        let nc = self.ptr.len() - 1;
        anyhow::ensure!(
            field.len() == nc,
            "the matrix temperature of {} cells, {nc} expected",
            field.len()
        );
        for (r, &c) in self.cell.iter().enumerate() {
            self.ts[r] = field[c as usize];
        }
        if let (Some(lv), Some(o)) = (&mut self.old, old)
            && o.len() == nc
        {
            for (r, &c) in self.cell.iter().enumerate() {
                lv[0][r] = o[c as usize];
            }
        }
        Ok(())
    }

    /// Solves the energy rows `val`, `rhs` (the device assembly, relaxed
    /// with `alpha`; `diag` the relaxed diagonal) together with the matrix
    /// rows; `phi` holds the gas temperature (start value in, solution
    /// out). `time` = (1/Δt, a₀, a₁, a₂) of a time step. Returns the
    /// iterations, ‖b − A x₀‖ and the residual scale of the extended
    /// system (as [`foehncraft_kernels::scalar::residual_parts`]).
    #[allow(clippy::too_many_arguments)]
    pub fn solve(
        &mut self,
        dev: &mut dyn FlowDevice,
        val: Buf,
        rhs: Buf,
        diag: Buf,
        phi: Buf,
        alpha: f64,
        time: Option<(f64, [f64; 3])>,
        tol: f64,
        max_iter: usize,
    ) -> Result<(usize, f64, f64)> {
        let v = dev.download(val)?;
        let mut b = dev.download(rhs)?;
        let dg = dev.download(diag)?;
        let tg = dev.download(phi)?;
        let nc = self.ptr.len() - 1;
        let ns = self.len();
        let n = nc + ns;
        // h a V per row
        let hv: Vec<f64> = (0..ns)
            .map(|r| {
                let g = &self.regions[self.region[r] as usize];
                g.h_at(tg[self.cell[r] as usize]) * g.area * self.vol[r]
            })
            .collect();
        // the rows of the extended matrix: (column, value), the diagonal first
        let mut rows: Vec<Vec<(u32, f64)>> = Vec::with_capacity(n);
        let mut a_unrel = vec![0.0; n];
        for c in 0..nc {
            let mut r: Vec<(u32, f64)> = (self.ptr[c]..self.ptr[c + 1])
                .map(|k| (self.col[k], v[self.pos[k] as usize]))
                .collect();
            a_unrel[c] = dg[c] * alpha;
            rows.push(std::mem::take(&mut r));
        }
        for (s, &h) in hv.iter().enumerate() {
            let c = self.cell[s] as usize;
            if let Some(d) = rows[c].iter_mut().find(|e| e.0 as usize == c) {
                d.1 += h;
            }
            rows[c].push(((nc + s) as u32, -h));
            a_unrel[c] += h;
        }
        // the matrix rows: exchange, conduction, storage
        let mut srow: Vec<Vec<(u32, f64)>> = (0..ns)
            .map(|s| vec![((nc + s) as u32, hv[s]), (self.cell[s], -hv[s])])
            .collect();
        let mut sb = vec![0.0; ns];
        for f in &self.faces {
            let [p, q] = f.rows.map(|x| x as usize);
            let t = 0.5 * (self.ts[p] + self.ts[q]);
            let (gp, gq) = (
                &self.regions[self.region[p] as usize],
                &self.regions[self.region[q] as usize],
            );
            let kp = gp.k_normal(t, f.n);
            let kq = gq.k_normal(t, f.n);
            // series of the two half cells (the owner weight w is the share
            // of the neighbour's distance)
            let k = if self.region[p] == self.region[q] {
                kp
            } else if kp > 0.0 && kq > 0.0 {
                1.0 / ((1.0 - f.w) / kp + f.w / kq)
            } else {
                0.0
            };
            let g = k * f.dgeo;
            srow[p][0].1 += g;
            srow[q][0].1 += g;
            srow[p].push(((nc + q) as u32, -g));
            srow[q].push(((nc + p) as u32, -g));
        }
        if let (Some((inv_dt, a)), Some(old)) = (time, &self.old) {
            for s in 0..ns {
                let g = &self.regions[self.region[s] as usize];
                let m = (1.0 - g.eps) * g.solid.rho * self.vol[s] * inv_dt;
                let t = self.ts[s];
                let cp = if g.h.is_constant() {
                    g.solid.cp.at(g.solid.t_ref)
                } else {
                    g.solid.cp.at(t)
                };
                srow[s][0].1 += m * a[0] * cp;
                sb[s] -= m
                    * (a[0] * (g.h.at(t) - cp * t)
                        + a[1] * g.h.at(old[0][s])
                        + a[2] * g.h.at(old[1][s]));
            }
        }
        for s in 0..ns {
            let a = srow[s][0].1;
            a_unrel[nc + s] = a;
            srow[s][0].1 = a / alpha;
            sb[s] += (1.0 / alpha - 1.0) * a * self.ts[s];
        }
        rows.extend(srow);
        b.extend(sb);
        let a = foehncraft_fv::amg::Csr::from_rows(n, rows);
        let mut x: Vec<f64> = tg[..nc].iter().chain(&self.ts).copied().collect();
        // the residual scale ‖a x‖ + ‖b − (a/α − a) x‖
        let ad = diagonal(&a);
        let (mut sv, mut sw) = (0.0, 0.0);
        for i in 0..n {
            let p = a_unrel[i] * x[i];
            let q = b[i] - (ad[i] - a_unrel[i]) * x[i];
            sv += p * p;
            sw += q * q;
        }
        let scale = (sv.sqrt() + sw.sqrt()).max(f64::MIN_POSITIVE);
        let (its, r0) = bicgstab(&a, &ad, &b, &mut x, tol, max_iter);
        for (s, t) in self.ts.iter_mut().enumerate() {
            *t = x[nc + s];
        }
        for s in 0..ns {
            self.exchange[s] = hv[s] * (x[self.cell[s] as usize] - self.ts[s]);
        }
        dev.upload(phi, &x[..nc])?;
        Ok((its, r0, scale))
    }
}

fn diagonal(a: &foehncraft_fv::amg::Csr) -> Vec<f64> {
    (0..a.n_rows)
        .map(|i| {
            a.row(i)
                .find(|e| e.0 == i)
                .map_or(1.0, |e| if e.1 != 0.0 { e.1 } else { 1.0 })
        })
        .collect()
}

fn spmv(a: &foehncraft_fv::amg::Csr, x: &[f64], y: &mut [f64]) {
    for (i, yi) in y.iter_mut().enumerate() {
        *yi = a.row(i).map(|(j, v)| v * x[j]).sum();
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// BiCGStab with Jacobi on the host (deterministic: sequential sums), to
/// ‖r‖ ≤ `tol` ‖r₀‖. Returns the iterations and ‖r₀‖.
fn bicgstab(
    a: &foehncraft_fv::amg::Csr,
    d: &[f64],
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iter: usize,
) -> (usize, f64) {
    let n = b.len();
    let mut r = vec![0.0; n];
    spmv(a, x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }
    let r0 = dot(&r, &r).sqrt();
    let target = tol * r0;
    if r0 == 0.0 {
        return (0, 0.0);
    }
    let rh = r.clone();
    let mut p = r.clone();
    let (mut ph, mut vv, mut s, mut sh, mut t) = (
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
    );
    let mut rho = dot(&rh, &r);
    for it in 1..=max_iter {
        if rho == 0.0 {
            return (it - 1, r0);
        }
        for i in 0..n {
            ph[i] = p[i] / d[i];
        }
        spmv(a, &ph, &mut vv);
        let rv = dot(&rh, &vv);
        if rv == 0.0 {
            return (it - 1, r0);
        }
        let alpha = rho / rv;
        for i in 0..n {
            s[i] = r[i] - alpha * vv[i];
        }
        if dot(&s, &s).sqrt() <= target {
            for i in 0..n {
                x[i] += alpha * ph[i];
            }
            return (it, r0);
        }
        for i in 0..n {
            sh[i] = s[i] / d[i];
        }
        spmv(a, &sh, &mut t);
        let tt = dot(&t, &t);
        let omega = if tt > 0.0 { dot(&t, &s) / tt } else { 0.0 };
        for i in 0..n {
            x[i] += alpha * ph[i] + omega * sh[i];
            r[i] = s[i] - omega * t[i];
        }
        if dot(&r, &r).sqrt() <= target || omega == 0.0 {
            return (it, r0);
        }
        let rho_new = dot(&rh, &r);
        let beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * vv[i]);
        }
    }
    (max_iter, r0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bicgstab_solves_a_small_nonsymmetric_system() {
        let a = foehncraft_fv::amg::Csr::from_rows(
            3,
            vec![
                vec![(0, 4.0), (1, -1.0)],
                vec![(0, -2.0), (1, 5.0), (2, -1.0)],
                vec![(1, -1.5), (2, 3.0)],
            ],
        );
        let xe = [1.0, 2.0, 3.0];
        let mut b = vec![0.0; 3];
        spmv(&a, &xe, &mut b);
        let mut x = vec![0.0; 3];
        let d = diagonal(&a);
        bicgstab(&a, &d, &b, &mut x, 1e-14, 50);
        for i in 0..3 {
            assert!((x[i] - xe[i]).abs() < 1e-10, "{x:?}");
        }
    }
}
