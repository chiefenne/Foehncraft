//! The deforming mesh (D91, stage MS0 of the moving solids plan): the
//! nodes of the mesh move in time by a prescribed displacement, the cells
//! change their shape and volume, and the equations are those of a moving
//! control volume (arbitrary Lagrangian–Eulerian form). Case file: the
//! `[mesh_motion]` table.
//!
//! ```toml
//! [constants]
//! stroke = 0.03
//! omega = 6.283185307179586
//!
//! [mesh_motion]   # an evenly compressed cylinder of height 0.04 m (L22)
//! displacement = [0.0, "-y / 0.04 * stroke / 2 * (1 - cos(omega * t))", 0.0]
//! ```
//!
//! # The method
//!
//! - **Geometry per step.** Before every time step the nodes are put at
//!   x₀ + **d**(x₀, tⁿ⁺¹), with x₀ the node in the mesh file, and the face
//!   areas, centroids and cell volumes of the moved mesh are computed. The
//!   solver is built again on the moved mesh (as for the imprinted sliding
//!   mesh, D182 variant B) and takes over the fields and the old time
//!   levels: the cells keep their numbers.
//! - **Swept volumes.** Every face gets the volume δV_f that it sweeps
//!   between the two node positions
//!   ([`foehncraft_mesh::geometry::swept_volume`]), computed from the same
//!   fan of triangles as the cell volumes, so that the swept volumes of a
//!   cell add up to the change of its volume to round-off: the space
//!   conservation law (Demirdžić & Perić 1988) holds exactly, for any face
//!   shape.
//! - **Grid flux.** With the backward-difference coefficients (a₀, −a₁, −a₂)
//!   of the step, the volume rate of a face is Ḟ_f = (a₀ δV_fⁿ⁺¹ +
//!   a₂ δV_fⁿ)/Δt (δVⁿ: the sweep of the step before); then Σ_f ±Ḟ_f =
//!   (a₀ Vⁿ⁺¹ − a₁ Vⁿ − a₂ Vⁿ⁻¹)/Δt in every cell, the discrete dV/dt of the
//!   time terms. The grid mass flux ρ_f Ḟ_f is subtracted from the absolute
//!   face mass flux ṁ_f wherever a flux convects (momentum, scalars,
//!   energy, turbulence): the convecting flux is the relative one,
//!   ṁ_f − ρ_f Ḟ_f, as for the turning zones ([`crate::frame`]).
//! - **Time terms.** Each time level enters with its own volume:
//!   d(ρ V u)/dt = (a₀ ρⁿ⁺¹Vⁿ⁺¹uⁿ⁺¹ − a₁ ρⁿVⁿuⁿ − a₂ ρⁿ⁻¹Vⁿ⁻¹uⁿ⁻¹)/Δt.
//!   The rows keep the volume of the new level and take the old
//!   levels with the densities ρᵏ Vᵏ/Vⁿ⁺¹.
//! - **Continuity.** d(ρV)/dt + Σ_f (ṁ_f − ρ_f Ḟ_f) = 0. With a constant
//!   density the space conservation law makes this Σ_f ṁ_f = 0, the
//!   continuity of a fixed mesh. With the variable density of stage 7c the
//!   continuity rows take the mass source d(ρV)/dt − Σ_f ±ρ_f Ḟ_f, and
//!   p₀(t) of a closed volume follows from its mass with the volumes of the
//!   moved mesh: a gas compressed by a piston heats up with the pressure
//!   work V dp₀/dt (L22).
//! - **Walls.** A wall moves with its nodes: its velocity normal to itself
//!   is the grid velocity Ḟ_f/|S_f| (the given wall velocity keeps its
//!   tangential part), so the relative flux through a wall is zero to
//!   round-off. Every other boundary may move only in its own plane (nodes
//!   that slide along a symmetry plane, a wedge or an inlet): a face of
//!   such a boundary that sweeps a volume stops the run (MOV003).
//!
//! A uniform flow on a deforming mesh stays uniform to round-off
//! (free-stream preservation, check MS0.1), and a gas compressed slowly
//! and adiabatically keeps p₀ V^γ (check MS0.2, the L22 piston).
//!
//! # Limits (MS0)
//!
//! One process (no ranks yet), no restart, no turning zones, conjugate
//! regions, sprays, immersed bodies, porous regions, periodic boundaries
//! or error-controlled time steps; the displacement is prescribed (the
//! mesh motion solver for moving walls is stage MS1). The solver is built
//! again in every step, with its multigrid setup: a deforming case costs
//! about one setup per step more than a fixed one.

use crate::case::FlowCase;
use crate::device::Buf;
use crate::diag::{self, Diagnostic};
use crate::expr::{Value, Vars, Vector, vector_value};
use crate::simplec::Flow;
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::geometry::swept_volume;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// `[mesh_motion]`: the prescribed motion of the mesh nodes (stage MS0,
/// D91).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MeshMotionSection {
    /// Displacement of every mesh node from its place in the mesh file (m):
    /// three numbers or expressions of x, y, z (the node's place in the
    /// mesh file, m) and t (s) (D53). The walls move with their nodes;
    /// every other boundary may move only in its own plane.
    pub displacement: Vector,
}

/// The deforming mesh at the end of a run (`result.json`).
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct MotionSummary {
    /// The volume of the mesh at the end (m³).
    pub volume: f64,
    /// The smallest ratio of a cell volume to its volume in the mesh file
    /// during the run.
    pub volume_ratio_min: f64,
    /// The largest one.
    pub volume_ratio_max: f64,
    /// The largest closure |Σ_f ±δV_f − ΔV|/V of a cell in a step (the space
    /// conservation law; round-off).
    pub closure_max: f64,
}

/// The time-step coefficients and the swept volumes of the step that a
/// deforming mesh gives the flow solver: what the time terms, the
/// continuity and the convecting flux need.
#[derive(Debug, Clone)]
pub struct Grid {
    /// Vⁿ and Vⁿ⁻¹ per cell (Vⁿ⁻¹ empty before the second step).
    vol_old: [Vec<f64>; 2],
    /// δVⁿ⁺¹ and δVⁿ per face (δVⁿ empty before the second step).
    swept: [Vec<f64>; 2],
    /// Wall faces: face, unit normal (out of the domain) and area.
    walls: Vec<(usize, [f64; 3], f64)>,
    /// Volume rate Ḟ_f per face of the current step (m³/s; zero before
    /// the step starts).
    pub rate: Vec<f64>,
    /// The device buffers of the grid mass flux ρ_f Ḟ_f and of the
    /// relative flux.
    pub(crate) flux: Option<(Buf, Buf)>,
}

impl Grid {
    /// The volume rates Ḟ_f of a step with the coefficients `a` (a₀, −a₁,
    /// −a₂) and the step `dt`.
    pub fn set_step(&mut self, a: [f64; 3], dt: f64) {
        let two = a[2] != 0.0 && !self.swept[1].is_empty();
        self.rate = (0..self.swept[0].len())
            .map(|f| {
                let mut r = a[0] * self.swept[0][f];
                if two {
                    r += a[2] * self.swept[1][f];
                }
                r / dt
            })
            .collect();
    }

    /// The densities ρᵏ Vᵏ/Vⁿ⁺¹ of the old levels k = n, n−1 that the time
    /// terms of the rows (volume Vⁿ⁺¹ = `vol`) take, from the old cell
    /// densities `rho_old` (level n, n−1).
    pub fn old_densities(&self, vol: &[f64], rho_old: [&[f64]; 2]) -> [Vec<f64>; 2] {
        let lv = |k: usize| -> Vec<f64> {
            let v = if self.vol_old[k].is_empty() {
                &self.vol_old[0]
            } else {
                &self.vol_old[k]
            };
            (0..vol.len())
                .map(|c| rho_old[k][c] * v[c] / vol[c])
                .collect()
        };
        [lv(0), lv(1)]
    }

    /// The rate of change d(ρV)/dt per cell with the coefficients `a`, the
    /// step `dt`, the new density `rho` and the old ones.
    pub fn mass_rate(
        &self,
        vol: &[f64],
        rho: &[f64],
        rho_old: [&[f64]; 2],
        a: [f64; 3],
        dt: f64,
    ) -> Vec<f64> {
        let [o1, o2] = self.old_densities(vol, rho_old);
        (0..vol.len())
            .map(|c| {
                let mut d = a[0] * rho[c] - a[1] * o1[c];
                if a[2] != 0.0 {
                    d -= a[2] * o2[c];
                }
                d * vol[c] / dt
            })
            .collect()
    }

    /// The grid mass flux ρ_f Ḟ_f per face with the face densities
    /// `rho_face`.
    pub fn mass_flux(&self, rho_face: &[f64]) -> Vec<f64> {
        self.rate
            .iter()
            .zip(rho_face)
            .map(|(r, rho)| r * rho)
            .collect()
    }

    /// Σ_f ±ρ_f Ḟ_f per cell (out of the cell) of the grid mass flux `m`.
    pub fn cell_sum(disc: &foehncraft_fv::Discretisation, m: &[f64]) -> Vec<f64> {
        let mut s = vec![0.0; disc.n_cells];
        for (f, &mf) in m.iter().enumerate() {
            s[disc.owner[f] as usize] += mf;
            if f < disc.n_interior {
                s[disc.neighbour[f] as usize] -= mf;
            }
        }
        s
    }

    /// The boundary velocities `ub` with the walls moving with the mesh:
    /// the normal part of a wall's velocity is the grid velocity Ḟ_f/|S_f|.
    pub fn wall_velocities(&self, n_interior: usize, ub: &[Vec<f64>; 3]) -> [Vec<f64>; 3] {
        let mut out = ub.clone();
        if self.rate.is_empty() {
            return out;
        }
        for &(f, n, area) in &self.walls {
            let fb = f - n_interior;
            let u = [ub[0][fb], ub[1][fb], ub[2][fb]];
            let un = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
            let w = self.rate[f] / area;
            for i in 0..3 {
                out[i][fb] = u[i] + (w - un) * n[i];
            }
        }
        out
    }
}

/// The deforming mesh of a run: the nodes of the mesh file, the moved
/// mesh of the current time level and the volumes of the levels before.
pub struct Deforming {
    displacement: [Value; 3],
    base: Vec<[f64; 3]>,
    /// The cell volumes of the mesh file.
    base_vol: Vec<f64>,
    /// The mesh at the current time level.
    pub mesh: Mesh,
    /// Its geometry.
    pub geo: Geometry,
    /// Vⁿ⁻¹ (empty before the first step).
    vol_prev: Vec<f64>,
    /// δVⁿ of the last step (empty before the first step).
    swept_prev: Vec<f64>,
    /// Faces of the walls, and of the other boundaries with their names.
    walls: Vec<usize>,
    fixed: Vec<(usize, String)>,
    /// A domain without a pressure boundary and with a constant density: its
    /// volume must not change (MOV004).
    closed_constant: bool,
    /// The time of the current level.
    pub time: f64,
    /// The smallest ratio of a cell volume to its volume in the mesh file
    /// so far, and the largest.
    pub volume_ratio: (f64, f64),
}

/// The checks of `[mesh_motion]` that need no mesh (MOV001).
pub fn validate(case: &FlowCase) -> Result<()> {
    let Some(m) = &case.mesh_motion else {
        return Ok(());
    };
    let fail = |msg: String| -> Result<()> {
        Err(diag::err(
            Diagnostic::new(diag::MOTION_CONFIG, msg).key("mesh_motion"),
        ))
    };
    if case.time.is_none() {
        return fail(
            "[mesh_motion] moves the mesh in time: the case needs a [time] table (a transient \
             run)"
                .into(),
        );
    }
    let what = if !case.mrf.is_empty() || !case.sliding.is_empty() {
        Some("turning zones ([[mrf]], [[sliding]])")
    } else if !case.regions.is_empty() {
        Some("regions (conjugate heat transfer, porous media)")
    } else if !case.spray.is_empty() {
        Some("sprays")
    } else if !case.immersed.is_empty() {
        Some("immersed bodies")
    } else if case.boundary.values().any(|b| b.is_periodic()) {
        Some("periodic boundaries")
    } else if case.run.restart.is_some() {
        Some("a restart")
    } else if !case.run.ranks.is_empty() && case.run.ranks.iter().any(|&r| r > 1) {
        Some("several ranks")
    } else if case
        .time
        .is_some_and(|t| t.error_tolerance.is_some() && t.adaptive)
    {
        Some("an error-controlled time step (a rejected step would need the mesh of another time)")
    } else if case.verification.levels.len() > 1 {
        Some("a sweep over levels")
    } else {
        None
    };
    if let Some(w) = what {
        return fail(format!(
            "the deforming mesh (stage MS0) does not work with {w} yet"
        ));
    }
    let env = case.env()?;
    vector_value(&m.displacement, &env, Vars::ALL, "mesh_motion.displacement")?;
    Ok(())
}

impl Deforming {
    /// The deforming mesh of `case` on the mesh `mesh` (`None` without
    /// `[mesh_motion]`), with the nodes at the start time `t0`.
    pub fn new(case: &FlowCase, mesh: &Mesh, t0: f64) -> Result<Option<Self>> {
        let Some(m) = &case.mesh_motion else {
            return Ok(None);
        };
        validate(case)?;
        let env = case.env()?;
        let displacement =
            vector_value(&m.displacement, &env, Vars::ALL, "mesh_motion.displacement")?;
        let ni = mesh.n_interior_faces();
        let mut walls = Vec::new();
        let mut fixed = Vec::new();
        for p in &mesh.boundaries {
            let wall = matches!(
                case.boundary.get(&p.name),
                Some(crate::bc::BoundaryCondition::Wall { .. })
            );
            for f in p.range() {
                debug_assert!(f >= ni);
                if wall {
                    walls.push(f);
                } else {
                    fixed.push((f, p.name.clone()));
                }
            }
        }
        let open = mesh.boundaries.iter().any(|p| {
            case.boundary
                .get(&p.name)
                .is_some_and(crate::bc::BoundaryCondition::is_pressure)
        });
        let closed_constant = !open && !crate::density::wanted(case);
        let geo = Geometry::compute(mesh);
        let mut d = Self {
            displacement,
            base: mesh.points.clone(),
            base_vol: geo.cell_volume.clone(),
            mesh: mesh.clone(),
            geo,
            vol_prev: Vec::new(),
            swept_prev: Vec::new(),
            walls,
            fixed,
            closed_constant,
            time: t0,
            volume_ratio: (1.0, 1.0),
        };
        d.mesh.points = d.points_at(t0);
        d.geo = Geometry::compute(&d.mesh);
        d.check_volumes(t0)?;
        Ok(Some(d))
    }

    /// The node coordinates at time `t`.
    fn points_at(&self, t: f64) -> Vec<[f64; 3]> {
        self.base
            .iter()
            .map(|&x| {
                [
                    x[0] + self.displacement[0].at(x, t),
                    x[1] + self.displacement[1].at(x, t),
                    x[2] + self.displacement[2].at(x, t),
                ]
            })
            .collect()
    }

    /// Fails when a cell of the current mesh has no positive volume
    /// (MOV002); keeps the range of the volume ratios.
    fn check_volumes(&mut self, t: f64) -> Result<()> {
        let n = self.mesh.n_cells;
        let base = &self.base_vol;
        let mut worst = (f64::INFINITY, 0);
        for (c, (v, b)) in self.geo.cell_volume[..n].iter().zip(base).enumerate() {
            let r = v / b;
            if r < worst.0 {
                worst = (r, c);
            }
            self.volume_ratio.0 = self.volume_ratio.0.min(r);
            self.volume_ratio.1 = self.volume_ratio.1.max(r);
        }
        // (a NaN ratio never becomes the worst, so `worst` stays infinite)
        if worst.0 <= 0.0 || !worst.0.is_finite() {
            let x = self.geo.cell_centroid[worst.1];
            return Err(diag::err(
                Diagnostic::new(
                    diag::MOTION_INVERTED,
                    format!(
                        "the deforming mesh turns cell {} inside out at t = {t:.6e} s (volume \
                         {:.3e} of its volume in the mesh file, near ({:.4e}, {:.4e}, {:.4e}) m): \
                         the displacement of [mesh_motion] folds the mesh there",
                        worst.1, worst.0, x[0], x[1], x[2]
                    ),
                )
                .key("mesh_motion.displacement"),
            ));
        }
        Ok(())
    }

    /// Moves the mesh to time `t1` (the end of the next step): the new mesh
    /// and geometry, and the [`Grid`] of the step.
    pub fn advance(&mut self, t1: f64) -> Result<Grid> {
        let old = std::mem::take(&mut self.mesh.points);
        let new = self.points_at(t1);
        let swept: Vec<f64> = (0..self.mesh.n_faces())
            .map(|f| swept_volume(&self.mesh, &old, &new, f))
            .collect();
        let vol_now = self.geo.cell_volume[..self.mesh.n_cells].to_vec();
        self.mesh.points = new;
        self.geo = Geometry::compute(&self.mesh);
        self.check_volumes(t1)?;
        // a boundary other than a wall moves in its own plane only
        for (f, name) in &self.fixed {
            let c = self.mesh.owner[*f] as usize;
            let v = self.geo.cell_volume[c];
            if swept[*f].abs() > 1e-9 * v {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::MOTION_BOUNDARY,
                        format!(
                            "the deforming mesh moves boundary `{name}` out of its plane between \
                             t = {:.6e} s and {t1:.6e} s (a face sweeps {:.3e} of the volume of its \
                             cell): only walls move with the mesh; the nodes of other boundaries \
                             may slide along them",
                            self.time,
                            swept[*f].abs() / v
                        ),
                    )
                    .boundary(name),
                ));
            }
        }
        if self.closed_constant {
            let v0: f64 = vol_now.iter().copied().exact_sum();
            let v1 = self.volume();
            if (v1 - v0).abs() > 1e-9 * v0 {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::MOTION_CLOSED,
                        format!(
                            "the deforming mesh changes the volume of a closed domain by {:.3e}                              between t = {:.6e} s and {t1:.6e} s, and the density is constant:                              an incompressible fluid cannot be compressed. Give the domain a                              pressure boundary, or let the density follow the pressure ([energy]                              density = \"variable\" with a closed volume)",
                            (v1 - v0) / v0,
                            self.time
                        ),
                    )
                    .key("mesh_motion.displacement"),
                ));
            }
        }
        let walls = self
            .walls
            .iter()
            .map(|&f| {
                let (n, a) = self.geo.normal_area(f);
                (f, n, a)
            })
            .collect();
        let g = Grid {
            vol_old: [vol_now.clone(), std::mem::take(&mut self.vol_prev)],
            swept: [swept.clone(), std::mem::take(&mut self.swept_prev)],
            walls,
            rate: Vec::new(),
            flux: None,
        };
        self.vol_prev = vol_now;
        self.swept_prev = swept;
        self.time = t1;
        Ok(g)
    }

    /// The total volume of the current mesh (m³).
    pub fn volume(&self) -> f64 {
        self.geo.cell_volume[..self.mesh.n_cells]
            .iter()
            .copied()
            .exact_sum()
    }

    /// The largest |Σ_f ±δV_f − ΔV| / V of the cells for the last step: the
    /// closure of the space conservation law (round-off).
    pub fn closure(&self, grid: &Grid) -> f64 {
        let m = &self.mesh;
        let mut worst = 0.0_f64;
        for c in 0..m.n_cells {
            let mut s = 0.0;
            for (f, is_nb) in m.cell_face_iter(c) {
                s += if is_nb {
                    -grid.swept[0][f]
                } else {
                    grid.swept[0][f]
                };
            }
            let v = self.geo.cell_volume[c];
            let dv = v - grid.vol_old[0][c];
            worst = worst.max((s - dv).abs() / v);
        }
        worst
    }
}

impl Flow {
    /// Gives the flow solver the grid of a deforming mesh for the next
    /// step (stage MS0).
    pub fn attach_grid(&mut self, mut g: Grid) -> Result<()> {
        let nf = self.disc.n_faces;
        let flux = self.dev.alloc(nf)?;
        let conv = self.dev.alloc(nf)?;
        self.dev.upload(flux, &vec![0.0; nf])?;
        g.flux = Some((flux, conv));
        self.grid = Some(Box::new(g));
        Ok(())
    }

    /// The start of a time step on a deforming mesh: the volume rates of
    /// the faces with the coefficients `a` and the step `dt`, the walls
    /// moving with the mesh, the grid mass flux and the old-level
    /// densities of the time terms (a constant density; the variable one
    /// sets them in [`Flow::refresh_density`]).
    pub(crate) fn grid_step_start(&mut self, a: [f64; 3], dt: f64) -> Result<()> {
        let Some(g) = self.grid.as_deref_mut() else {
            return Ok(());
        };
        g.set_step(a, dt);
        let ub = self.bc.ub.clone();
        self.set_boundary_velocity(&ub)?;
        if self.density.is_none() {
            let rho = self.phys.rho;
            let g = self.grid.as_deref().expect("grid");
            let nf = self.disc.n_faces;
            let m: Vec<f64> = g.rate.iter().map(|r| rho * r).collect();
            debug_assert_eq!(m.len(), nf);
            let nc = self.disc.n_cells;
            let vol = &self.disc.cell_volume[..nc];
            let r = vec![rho; nc];
            let [o1, o2] = g.old_densities(vol, [&r, &r]);
            let (flux, _) = g.flux.expect("grid buffers");
            self.dev.upload(flux, &m)?;
            self.dev.set_variable_density(None, Some(&o1), Some(&o2))?;
            self.phys.grid_flux = Some(std::sync::Arc::new(m));
        }
        Ok(())
    }

    /// The relative flux of `m` on a deforming mesh (`m` itself without
    /// one).
    pub(crate) fn grid_relative(&mut self, m: Buf) -> Result<Buf> {
        let Some((flux, conv)) = self.grid.as_deref().and_then(|g| g.flux) else {
            return Ok(m);
        };
        self.dev.lincomb(conv, 1.0, m, -1.0, flux)?;
        Ok(conv)
    }
}
