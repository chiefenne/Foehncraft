//! Restart files (D42): the solution and the solver state of a flow run in
//! one HDF5 file, so that a run continued from it gives bitwise the result
//! of the run without the interruption (same mesh, same decomposition and
//! device).
//!
//! Layout (format `foehncraft-restart`, version 1):
//!
//! - attributes `format`, `version`, `case`, `iteration`, `amg_setups`,
//!   `transient`, and for a transient run `time`, `step`, `dt_prev`,
//!   `dt_prev2`, `levels`, `dt_next`, `scheme`;
//! - per cell (local order, with the global identities in `cell_id`):
//!   `u_x u_y u_z p`, the gradients `grad_u_<i><j>` and `grad_p_<j>`, the
//!   limiters `psi_<i>`, and for a transient run the old levels `u_old_<i>`
//!   and `g_old_<k>` (velocity gradients and limiters of the level
//!   before, for the extrapolation), with the error estimate of the time
//!   step (D65) also `u_old2_<i>` (the level before that);
//! - per face (with `face_id`): `mass_flux`, and transient `m_old`, `x_cur`,
//!   `x_old`; per boundary face `ub_<i>`, `pb`, and with convective outlets
//!   (D29) the boundary velocity of the level before, `ub_old_<i>`;
//! - the passive scalars (D51): attribute `scalars` (the names), per cell
//!   `scalar_<name>` and, transient, `scalar_old_<name>` (the level before),
//!   with `limiter = "mlp"` also `scalar_psi_<name>`, the limiter value
//!   that its relaxation carries between outer iterations (D237);
//! - the periodic pairs (D29): `periodic_force`, the driving force per
//!   pair;
//! - the turbulence model (stage 6): attribute `turbulence` (the family:
//!   `k-omega-sst`, or `k-epsilon` when `turb_omega` holds ε), per cell
//!   `turb_k`, `turb_omega`, `turb_nu_t`, `turb_p_modified` (the
//!   solver's p̂; `p` is the static pressure), `turb_wall_distance`,
//!   transient `turb_k_old`, `turb_omega_old`; per face `turb_mu_eff`; per
//!   boundary face `turb_k_b`, `turb_omega_b`, `turb_u_tau`, `turb_y_plus`,
//!   `turb_pb_modified`; k-ε-ζ-f adds per cell `turb_zeta`, `turb_f`,
//!   `turb_time_scale` (transient `turb_zeta_old`) and per boundary face
//!   `turb_zeta_b`, `turb_f_b`; `turb_amg_matrix_<slot>` the matrices of
//!   the last multigrid setups of the turbulence equations (slot 1 f,
//!   2 k, 3 ε or ω, 4 ζ; with their multigrid only); `turb_fr1`, f_r1
//!   of the curvature correction of SST, and per boundary face
//!   `turb_wall_blending`, the viscous weight of the wall treatment (both
//!   output only, recomputed in every iteration);
//! - `amg_matrix`: the pressure matrix of the last multigrid setup (the
//!   hierarchy is rebuilt from it, so the preconditioner is the same);
//! - the coupled algorithm with the multigrid momentum sweep (D179):
//!   `momentum_amg_use`, the hierarchy (0 to 2, −1 none) of each velocity
//!   component, and `momentum_amg_matrix_<h>`, the momentum matrix
//!   hierarchy h was built from;
//! - a conjugate case (stage 7): `cht_cell_id`, `cht_temperature` and,
//!   transient, `cht_temperature_old` on the full mesh (fluid and solids),
//!   with the energy limiter `mlp` also `cht_psi` (its limiter value, D237),
//!   and with the cycle-averaged loop (D204) `cht_cycle_state`, the bytes of
//!   its state as JSON (counters, accumulators, update history), with the
//!   frozen-flow march (D118) `cht_frozen_state` (its phase, settling
//!   windows and counters), and with a viscosity that follows T
//!   `cht_viscosity` (μ(T) per face of the fluid view as the energy equation
//!   last passed it to the flow);
//! - the sprays (L12): `spray_state`, the bytes of the JSON state of the
//!   parcels (their positions, the injector counters and the tallies);
//! - optional, for meshes that move (D91, L22): `point_id`,
//!   `points_x|y|z` and the cell volumes of the time levels `volume_<k>`
//!   (attribute `volume_levels`), the face swept volumes of the steps
//!   before `swept_<k>` (attribute `swept_levels`, stage MS0 part 2) and
//!   the attribute `volume_ratio` (the smallest and largest ratio of a cell
//!   volume to its volume in the mesh file so far); a fixed mesh leaves
//!   them out.
//!
//! A file whose identities are in another order (another renumbering) is
//! mapped by identity; the multigrid matrix is then dropped and the
//! hierarchy rebuilt, so the continuation is no longer bit-identical.

use crate::simplec::State;
use crate::transient::{Clock, OldLevels, TimeScheme};
use anyhow::Result;
#[cfg(feature = "hdf5")]
use anyhow::{Context, bail, ensure};
use foehncraft_fv::disc::Soa3;
use std::path::Path;

/// The content of a restart file.
#[derive(Debug, Clone)]
pub struct Restart {
    /// Case name.
    pub case: String,
    /// Global cell identities in local order.
    pub cell_id: Vec<u64>,
    /// Global face identities in local order.
    pub face_id: Vec<u64>,
    /// The fields.
    pub state: State,
    /// Outer iterations done.
    pub iteration: usize,
    /// Multigrid setups done.
    pub amg_setups: usize,
    /// Pressure matrix of the last multigrid setup.
    pub amg_matrix: Option<Vec<f64>>,
    /// The kept hierarchies of the momentum sweeps of the coupled
    /// algorithm (D179).
    pub sweep: crate::coupled::SweepState,
    /// Transient state: old levels, clock, next time step, scheme.
    pub transient: Option<(OldLevels, Clock, f64, TimeScheme)>,
    /// The geometry of the state, for meshes that move or deform (D91,
    /// L22): `None` for a fixed mesh, whose mesh file defines it.
    pub geometry: Option<RestartGeometry>,
    /// The state of the sprays (L12; JSON of `foehncraft_particles::Spray::state`),
    /// dataset `spray_state`; parcels follow their positions, not the cell order.
    pub spray: Option<Vec<u8>>,
    /// The temperature of a conjugate case on the full mesh (fluid and
    /// solids, stage 7) and the state of its cycle-averaged loop (D204).
    pub cht: Option<ChtRestart>,
}

/// The energy equation of a conjugate case in a restart file: datasets
/// `cht_cell_id` (the global identities of the full mesh), `cht_temperature`
/// (K) and, transient, `cht_temperature_old` (the level before), and with
/// the cycle-averaged loop `cht_cycle_state` (bytes of its JSON state).
#[derive(Debug, Clone, Default)]
pub struct ChtRestart {
    /// Global cell identities of the full mesh in local order.
    pub cell_id: Vec<u64>,
    /// The temperature per cell of the full mesh (K).
    pub t: Vec<f64>,
    /// Transient: the level before.
    pub t_old: Option<Vec<f64>>,
    /// The limiter value ψ per cell of the full mesh with the energy
    /// limiter `mlp` (D237; its relaxation carries ψ between iterations).
    pub psi: Option<Vec<f64>>,
    /// The state of the cycle-averaged loop (JSON), if the run has one.
    pub cycle: Option<Vec<u8>>,
    /// The state of the frozen-flow march (JSON, D118), if the run has one.
    pub frozen: Option<Vec<u8>>,
    /// The face viscosity μ(T) of the fluid view as the energy equation
    /// last passed it to the flow, when the viscosity follows T.
    pub viscosity: Option<Vec<f64>>,
}

/// Point coordinates and cell volumes of the time levels (D91): optional
/// datasets `point_id`, `points_x|y|z`, `volume_<k>` (k = 0 the current
/// level, 1, 2 … the older ones). Fixed meshes leave them out; a deforming
/// mesh (L22) writes them, since its mesh file no longer holds the state.
/// Units: points m, volumes m³.
#[derive(Debug, Clone, Default)]
pub struct RestartGeometry {
    /// Global point identities in local order.
    pub point_id: Vec<u64>,
    /// Point coordinates (m).
    pub points: Soa3,
    /// Cell volumes per time level, the current first (m³).
    pub volumes: Vec<Vec<f64>>,
    /// The volume δV_f each face swept in the steps before, the last
    /// first (m³; the second-order time schemes need the sweep of the
    /// step before, stage MS0 part 2).
    pub swept: Vec<Vec<f64>>,
    /// The smallest and the largest ratio of a cell volume to its volume
    /// in the mesh file so far.
    pub volume_ratio: Option<[f64; 2]>,
}

#[cfg(feature = "hdf5")]
mod io {
    use super::*;
    use hdf5_metno as h5;

    fn ds(g: &h5::File, name: &str, data: &[f64]) -> Result<()> {
        g.new_dataset::<f64>()
            .shape([data.len()])
            .create(name)?
            .write_raw(data)?;
        Ok(())
    }

    fn ds_u64(g: &h5::File, name: &str, data: &[u64]) -> Result<()> {
        g.new_dataset::<u64>()
            .shape([data.len()])
            .create(name)?
            .write_raw(data)?;
        Ok(())
    }

    fn attr<T: h5::H5Type>(g: &h5::File, name: &str, v: &T) -> Result<()> {
        g.new_attr::<T>().create(name)?.write_scalar(v)?;
        Ok(())
    }

    fn sattr(g: &h5::File, name: &str, v: &str) -> Result<()> {
        let s: h5::types::VarLenUnicode = v.parse().map_err(|e| anyhow::anyhow!("{e:?}"))?;
        g.new_attr::<h5::types::VarLenUnicode>()
            .create(name)?
            .write_scalar(&s)?;
        Ok(())
    }

    fn rd(g: &h5::File, name: &str) -> Result<Vec<f64>> {
        g.dataset(name)
            .with_context(|| format!("restart file: dataset `{name}`"))?
            .read_raw::<f64>()
            .map_err(Into::into)
    }

    fn rattr<T: h5::H5Type>(g: &h5::File, name: &str) -> Result<T> {
        g.attr(name)
            .with_context(|| format!("restart file: attribute `{name}`"))?
            .read_scalar::<T>()
            .map_err(Into::into)
    }

    fn rsattr(g: &h5::File, name: &str) -> Result<String> {
        Ok(rattr::<h5::types::VarLenUnicode>(g, name)?
            .as_str()
            .to_owned())
    }

    fn put3(g: &h5::File, name: &str, v: &Soa3) -> Result<()> {
        ds(g, &format!("{name}x"), &v.x)?;
        ds(g, &format!("{name}y"), &v.y)?;
        ds(g, &format!("{name}z"), &v.z)
    }

    fn get3(g: &h5::File, name: &str) -> Result<Soa3> {
        Ok(Soa3 {
            x: rd(g, &format!("{name}x"))?,
            y: rd(g, &format!("{name}y"))?,
            z: rd(g, &format!("{name}z"))?,
        })
    }

    /// Writes `r` to `path` (through a temporary file, so a crash never
    /// leaves a broken restart behind).
    pub fn write(path: &Path, r: &Restart) -> Result<()> {
        let tmp = path.with_extension("h5.tmp");
        {
            // no object timestamps: the same state gives the same bytes (D84)
            let g = h5::File::with_options()
                .with_fcpl(|p| p.obj_track_times(false))
                .create(&tmp)
                .with_context(|| format!("creating {}", tmp.display()))?;
            sattr(&g, "format", "foehncraft-restart")?;
            attr(&g, "version", &1i64)?;
            sattr(&g, "foehncraft_version", foehncraft_core::VERSION)?;
            sattr(&g, "units", UNITS)?;
            sattr(&g, "case", &r.case)?;
            attr(&g, "iteration", &(r.iteration as u64))?;
            attr(&g, "amg_setups", &(r.amg_setups as u64))?;
            attr(&g, "transient", &u8::from(r.transient.is_some()))?;
            ds_u64(&g, "cell_id", &r.cell_id)?;
            ds_u64(&g, "face_id", &r.face_id)?;
            let s = &r.state;
            put3(&g, "u_", &s.u)?;
            ds(&g, "p", &s.p)?;
            ds(&g, "mass_flux", &s.mass_flux)?;
            put3(&g, "ub_", &s.ub)?;
            ds(&g, "pb", &s.pb)?;
            for (i, c) in ["x", "y", "z"].iter().enumerate() {
                put3(&g, &format!("grad_u_{c}"), &s.grad_u[i])?;
            }
            put3(&g, "grad_p_", &s.grad_p)?;
            put3(&g, "psi_", &s.psi)?;
            // the passive scalars (D51): `scalar_<name>`, transient also
            // `scalar_old_<name>`; attribute `scalars` lists the names
            if !s.scalars.is_empty() {
                let names: Vec<&str> = s.scalars.iter().map(|x| x.name.as_str()).collect();
                sattr(&g, "scalars", &names.join(" "))?;
                for x in &s.scalars {
                    ds(&g, &format!("scalar_{}", x.name), &x.phi)?;
                    if let Some(o) = &x.old {
                        ds(&g, &format!("scalar_old_{}", x.name), o)?;
                    }
                    if let Some(p) = &x.psi {
                        ds(&g, &format!("scalar_psi_{}", x.name), p)?;
                    }
                }
            }
            if !s.periodic_force.is_empty() {
                ds(&g, "periodic_force", &s.periodic_force)?;
            }
            // the turbulence model (stage 6): `turb_*`, attribute
            // `turbulence`
            if let Some(t) = &s.turbulence {
                // the second variable keeps its dataset name `turb_omega`;
                // the attribute says whether it holds ω or ε
                sattr(
                    &g,
                    "turbulence",
                    if t.epsilon {
                        "k-epsilon"
                    } else {
                        "k-omega-sst"
                    },
                )?;
                ds(&g, "turb_k", &t.k)?;
                ds(&g, "turb_omega", &t.omega)?;
                ds(&g, "turb_nu_t", &t.nu_t)?;
                ds(&g, "turb_mu_eff", &t.mu_eff)?;
                ds(&g, "turb_k_b", &t.k_b)?;
                ds(&g, "turb_omega_b", &t.omega_b)?;
                ds(&g, "turb_u_tau", &t.u_tau)?;
                ds(&g, "turb_y_plus", &t.y_plus)?;
                if let Some(b) = &t.wall_blending {
                    ds(&g, "turb_wall_blending", b)?;
                }
                ds(&g, "turb_p_modified", &t.p_modified)?;
                ds(&g, "turb_pb_modified", &t.pb_modified)?;
                ds(&g, "turb_wall_distance", &t.wall_distance)?;
                if let Some(o) = &t.k_old {
                    ds(&g, "turb_k_old", o)?;
                }
                if let Some(o) = &t.omega_old {
                    ds(&g, "turb_omega_old", o)?;
                }
                if let Some(z) = &t.zeta {
                    ds(&g, "turb_zeta", &z.zeta)?;
                    ds(&g, "turb_f", &z.f)?;
                    ds(&g, "turb_zeta_b", &z.zeta_b)?;
                    ds(&g, "turb_f_b", &z.f_b)?;
                    ds(&g, "turb_time_scale", &z.t)?;
                    if let Some(o) = &z.zeta_old {
                        ds(&g, "turb_zeta_old", o)?;
                    }
                }
                for (slot, m) in &t.amg {
                    ds(&g, &format!("turb_amg_matrix_{slot}"), m)?;
                }
                if let Some(f) = &t.fr1 {
                    ds(&g, "turb_fr1", f)?;
                }
            }
            if let Some(m) = &r.amg_matrix {
                ds(&g, "amg_matrix", m)?;
            }
            if r.sweep.slot.iter().any(Option::is_some) {
                let used = r.sweep.slot.map(|s| s.map_or(-1.0, |s| s as f64));
                ds(&g, "momentum_amg_use", &used)?;
                for (h, m) in r.sweep.matrix.iter().enumerate() {
                    if let Some(m) = m {
                        ds(&g, &format!("momentum_amg_matrix_{h}"), m)?;
                    }
                }
            }
            if let Some((o, clock, dt_next, scheme)) = &r.transient {
                attr(&g, "time", &clock.time)?;
                attr(&g, "step", &(clock.step as u64))?;
                attr(&g, "dt_prev", &clock.dt_prev)?;
                attr(&g, "dt_prev2", &clock.dt_prev2)?;
                attr(&g, "levels", &(clock.levels as u64))?;
                attr(&g, "dt_next", dt_next)?;
                sattr(&g, "scheme", scheme_name(*scheme))?;
                ds(&g, "u_old_x", &o.u_old[0])?;
                ds(&g, "u_old_y", &o.u_old[1])?;
                ds(&g, "u_old_z", &o.u_old[2])?;
                ds(&g, "m_old", &o.m_old)?;
                ds(&g, "x_cur", &o.x_cur)?;
                ds(&g, "x_old", &o.x_old)?;
                for (i, v) in o.g_old.iter().enumerate() {
                    ds(&g, &format!("g_old_{i}"), v)?;
                }
                if let Some(u) = &o.u_old2 {
                    ds(&g, "u_old2_x", &u[0])?;
                    ds(&g, "u_old2_y", &u[1])?;
                    ds(&g, "u_old2_z", &u[2])?;
                }
                if let Some(u) = &o.ub_old {
                    ds(&g, "ub_old_x", &u[0])?;
                    ds(&g, "ub_old_y", &u[1])?;
                    ds(&g, "ub_old_z", &u[2])?;
                }
            }
            if let Some(c) = &r.cht {
                ds_u64(&g, "cht_cell_id", &c.cell_id)?;
                ds(&g, "cht_temperature", &c.t)?;
                if let Some(o) = &c.t_old {
                    ds(&g, "cht_temperature_old", o)?;
                }
                if let Some(p) = &c.psi {
                    ds(&g, "cht_psi", p)?;
                }
                if let Some(b) = &c.cycle {
                    g.new_dataset::<u8>()
                        .shape([b.len()])
                        .create("cht_cycle_state")?
                        .write_raw(b)?;
                }
                if let Some(b) = &c.frozen {
                    g.new_dataset::<u8>()
                        .shape([b.len()])
                        .create("cht_frozen_state")?
                        .write_raw(b)?;
                }
                if let Some(m) = &c.viscosity {
                    ds(&g, "cht_viscosity", m)?;
                }
            }
            if let Some(b) = &r.spray {
                g.new_dataset::<u8>()
                    .shape([b.len()])
                    .create("spray_state")?
                    .write_raw(b)?;
            }
            if let Some(geo) = &r.geometry {
                ds_u64(&g, "point_id", &geo.point_id)?;
                put3(&g, "points_", &geo.points)?;
                attr(&g, "volume_levels", &(geo.volumes.len() as u64))?;
                for (k, v) in geo.volumes.iter().enumerate() {
                    ds(&g, &format!("volume_{k}"), v)?;
                }
                attr(&g, "swept_levels", &(geo.swept.len() as u64))?;
                for (k, v) in geo.swept.iter().enumerate() {
                    ds(&g, &format!("swept_{k}"), v)?;
                }
                if let Some(r) = geo.volume_ratio {
                    g.new_attr::<f64>()
                        .shape([2])
                        .create("volume_ratio")?
                        .write_raw(&r)?;
                }
            }
        }
        std::fs::rename(&tmp, path).with_context(|| format!("writing {}", path.display()))?;
        Ok(())
    }

    fn scheme_name(s: TimeScheme) -> &'static str {
        match s {
            TimeScheme::Bdf1 => "bdf1",
            TimeScheme::Bdf2 => "bdf2",
            TimeScheme::CrankNicolson => "crank-nicolson",
        }
    }

    /// The passive scalars of a restart file (none in files without them).
    fn read_scalars(g: &h5::File, transient: bool) -> Result<Vec<crate::scalar::ScalarState>> {
        if g.attr("scalars").is_err() {
            return Ok(Vec::new());
        }
        let names = rsattr(g, "scalars")?;
        names
            .split_whitespace()
            .map(|n| {
                let old = format!("scalar_old_{n}");
                let psi = format!("scalar_psi_{n}");
                Ok(crate::scalar::ScalarState {
                    name: n.to_string(),
                    phi: rd(g, &format!("scalar_{n}"))?,
                    old: if transient && g.link_exists(&old) {
                        Some(rd(g, &old)?)
                    } else {
                        None
                    },
                    psi: if g.link_exists(&psi) {
                        Some(rd(g, &psi)?)
                    } else {
                        None
                    },
                    wall: None,
                })
            })
            .collect()
    }

    /// The turbulence state of a restart file (none without one).
    fn read_turbulence(g: &h5::File) -> Result<Option<crate::turbulence::TurbState>> {
        if g.attr("turbulence").is_err() {
            return Ok(None);
        }
        let opt = |n: &str| -> Result<Option<Vec<f64>>> {
            Ok(if g.link_exists(n) {
                Some(rd(g, n)?)
            } else {
                None
            })
        };
        let epsilon = rsattr(g, "turbulence")? == "k-epsilon";
        Ok(Some(crate::turbulence::TurbState {
            epsilon,
            k: rd(g, "turb_k")?,
            omega: rd(g, "turb_omega")?,
            nu_t: rd(g, "turb_nu_t")?,
            k_old: opt("turb_k_old")?,
            omega_old: opt("turb_omega_old")?,
            mu_eff: rd(g, "turb_mu_eff")?,
            k_b: rd(g, "turb_k_b")?,
            omega_b: rd(g, "turb_omega_b")?,
            u_tau: rd(g, "turb_u_tau")?,
            y_plus: rd(g, "turb_y_plus")?,
            wall_blending: if g.link_exists("turb_wall_blending") {
                Some(rd(g, "turb_wall_blending")?)
            } else {
                None
            },
            p_modified: rd(g, "turb_p_modified")?,
            pb_modified: rd(g, "turb_pb_modified")?,
            wall_distance: rd(g, "turb_wall_distance")?,
            zeta: if g.link_exists("turb_zeta") {
                Some(crate::turbulence::ZetaState {
                    zeta: rd(g, "turb_zeta")?,
                    f: rd(g, "turb_f")?,
                    zeta_old: opt("turb_zeta_old")?,
                    zeta_b: rd(g, "turb_zeta_b")?,
                    f_b: rd(g, "turb_f_b")?,
                    t: rd(g, "turb_time_scale")?,
                })
            } else {
                None
            },
            amg: (1..=4)
                .filter_map(|s| {
                    let n = format!("turb_amg_matrix_{s}");
                    g.link_exists(&n).then(|| rd(g, &n).map(|m| (s, m)))
                })
                .collect::<Result<_>>()?,
            fr1: if g.link_exists("turb_fr1") {
                Some(rd(g, "turb_fr1")?)
            } else {
                None
            },
        }))
    }

    /// Kind and progress of a restart file, from its attributes.
    pub fn progress(path: &Path) -> Result<(bool, u64)> {
        let g = h5::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let fmt = rsattr(&g, "format")?;
        ensure!(
            fmt == "foehncraft-restart",
            "{} is not a Foehncraft restart file (format `{fmt}`)",
            path.display()
        );
        let transient = rattr::<u8>(&g, "transient")? == 1;
        let n: u64 = if transient {
            rattr(&g, "step")?
        } else {
            rattr(&g, "iteration")?
        };
        Ok((transient, n))
    }

    /// Reads a restart file.
    pub fn read(path: &Path) -> Result<Restart> {
        let g = h5::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let fmt = rsattr(&g, "format")?;
        ensure!(
            fmt == "foehncraft-restart",
            "{} is not a Foehncraft restart file (format `{fmt}`)",
            path.display()
        );
        let version: i64 = rattr(&g, "version")?;
        ensure!(
            version == 1,
            "restart file version {version} is not supported"
        );
        let transient = rattr::<u8>(&g, "transient")? == 1;
        let state = State {
            u: get3(&g, "u_")?,
            p: rd(&g, "p")?,
            mass_flux: rd(&g, "mass_flux")?,
            ub: get3(&g, "ub_")?,
            pb: rd(&g, "pb")?,
            grad_u: [
                get3(&g, "grad_u_x")?,
                get3(&g, "grad_u_y")?,
                get3(&g, "grad_u_z")?,
            ],
            grad_p: get3(&g, "grad_p_")?,
            psi: get3(&g, "psi_")?,
            scalars: read_scalars(&g, transient)?,
            periodic_force: if g.link_exists("periodic_force") {
                rd(&g, "periodic_force")?
            } else {
                Vec::new()
            },
            turbulence: read_turbulence(&g)?,
        };
        let amg_matrix = if g.link_exists("amg_matrix") {
            Some(rd(&g, "amg_matrix")?)
        } else {
            None
        };
        let mut sweep = crate::coupled::SweepState::default();
        if g.link_exists("momentum_amg_use") {
            let used = rd(&g, "momentum_amg_use")?;
            anyhow::ensure!(used.len() == 3, "momentum_amg_use: three values expected");
            for (i, u) in used.iter().enumerate() {
                if (0.0..3.0).contains(u) {
                    let h = *u as usize;
                    let name = format!("momentum_amg_matrix_{h}");
                    anyhow::ensure!(g.link_exists(&name), "{name} is missing");
                    if sweep.matrix[h].is_none() {
                        sweep.matrix[h] = Some(rd(&g, &name)?);
                    }
                    sweep.slot[i] = Some(h);
                }
            }
        }
        let tr = if transient {
            let scheme = match rsattr(&g, "scheme")?.as_str() {
                "bdf1" => TimeScheme::Bdf1,
                "bdf2" => TimeScheme::Bdf2,
                "crank-nicolson" => TimeScheme::CrankNicolson,
                other => bail!("unknown time scheme `{other}` in the restart file"),
            };
            Some((
                OldLevels {
                    u_old: [rd(&g, "u_old_x")?, rd(&g, "u_old_y")?, rd(&g, "u_old_z")?],
                    m_old: rd(&g, "m_old")?,
                    x_cur: rd(&g, "x_cur")?,
                    x_old: rd(&g, "x_old")?,
                    g_old: (0..12)
                        .map(|i| rd(&g, &format!("g_old_{i}")))
                        .collect::<Result<_>>()?,
                    u_old2: if g.link_exists("u_old2_x") {
                        Some([
                            rd(&g, "u_old2_x")?,
                            rd(&g, "u_old2_y")?,
                            rd(&g, "u_old2_z")?,
                        ])
                    } else {
                        None
                    },
                    ub_old: if g.link_exists("ub_old_x") {
                        Some([
                            rd(&g, "ub_old_x")?,
                            rd(&g, "ub_old_y")?,
                            rd(&g, "ub_old_z")?,
                        ])
                    } else {
                        None
                    },
                },
                Clock {
                    time: rattr(&g, "time")?,
                    step: rattr::<u64>(&g, "step")? as usize,
                    dt_prev: rattr(&g, "dt_prev")?,
                    levels: rattr::<u64>(&g, "levels")? as usize,
                    // files of earlier versions: no uⁿ⁻² (D65)
                    dt_prev2: if g.attr("dt_prev2").is_ok() {
                        rattr(&g, "dt_prev2")?
                    } else {
                        0.0
                    },
                },
                rattr(&g, "dt_next")?,
                scheme,
            ))
        } else {
            None
        };
        let geometry = if g.link_exists("point_id") {
            let n: u64 = rattr(&g, "volume_levels")?;
            Some(RestartGeometry {
                point_id: g.dataset("point_id")?.read_raw::<u64>()?,
                points: get3(&g, "points_")?,
                volumes: (0..n)
                    .map(|k| rd(&g, &format!("volume_{k}")))
                    .collect::<Result<_>>()?,
                swept: if g.attr("swept_levels").is_ok() {
                    let n: u64 = rattr(&g, "swept_levels")?;
                    (0..n)
                        .map(|k| rd(&g, &format!("swept_{k}")))
                        .collect::<Result<_>>()?
                } else {
                    Vec::new()
                },
                volume_ratio: match g.attr("volume_ratio") {
                    Ok(a) => {
                        let v = a.read_raw::<f64>()?;
                        ensure!(v.len() == 2, "volume_ratio has {} values", v.len());
                        Some([v[0], v[1]])
                    }
                    Err(_) => None,
                },
            })
        } else {
            None
        };
        let cht = if g.link_exists("cht_temperature") {
            Some(ChtRestart {
                cell_id: g.dataset("cht_cell_id")?.read_raw::<u64>()?,
                t: rd(&g, "cht_temperature")?,
                t_old: if g.link_exists("cht_temperature_old") {
                    Some(rd(&g, "cht_temperature_old")?)
                } else {
                    None
                },
                psi: if g.link_exists("cht_psi") {
                    Some(rd(&g, "cht_psi")?)
                } else {
                    None
                },
                cycle: if g.link_exists("cht_cycle_state") {
                    Some(g.dataset("cht_cycle_state")?.read_raw::<u8>()?)
                } else {
                    None
                },
                frozen: if g.link_exists("cht_frozen_state") {
                    Some(g.dataset("cht_frozen_state")?.read_raw::<u8>()?)
                } else {
                    None
                },
                viscosity: if g.link_exists("cht_viscosity") {
                    Some(rd(&g, "cht_viscosity")?)
                } else {
                    None
                },
            })
        } else {
            None
        };
        let spray = if g.link_exists("spray_state") {
            Some(g.dataset("spray_state")?.read_raw::<u8>()?)
        } else {
            None
        };
        Ok(Restart {
            geometry,
            spray,
            cht,
            case: rsattr(&g, "case")?,
            cell_id: g.dataset("cell_id")?.read_raw::<u64>()?,
            face_id: g.dataset("face_id")?.read_raw::<u64>()?,
            state,
            iteration: rattr::<u64>(&g, "iteration")? as usize,
            amg_setups: rattr::<u64>(&g, "amg_setups")? as usize,
            amg_matrix,
            sweep,
            transient: tr,
        })
    }
}

/// Units of the fields in a restart file (attribute `units`, D75).
pub const UNITS: &str = "u m/s; p Pa; mass_flux kg/s; ub m/s; pb Pa; grad_u 1/s; grad_p Pa/m; \
     time s; dt_prev s; dt_prev2 s; dt_next s";

/// Whether a restart file is transient, and its step (transient) or outer
/// iteration (steady), read from the file (D75), never from its name.
pub fn progress(path: &Path) -> Result<(bool, u64)> {
    #[cfg(feature = "hdf5")]
    {
        io::progress(path)
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = path;
        anyhow::bail!("restart files need the `hdf5` feature")
    }
}

/// Restart files are available (the `hdf5` feature); without it a run
/// writes none.
pub const AVAILABLE: bool = cfg!(feature = "hdf5");

/// Writes a restart file.
pub fn write(path: &Path, r: &Restart) -> Result<()> {
    #[cfg(feature = "hdf5")]
    {
        io::write(path, r)
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = (path, r);
        anyhow::bail!("restart files need the `hdf5` feature")
    }
}

/// Reads a restart file.
pub fn read(path: &Path) -> Result<Restart> {
    #[cfg(feature = "hdf5")]
    {
        io::read(path)
    }
    #[cfg(not(feature = "hdf5"))]
    {
        let _ = path;
        anyhow::bail!("restart files need the `hdf5` feature")
    }
}

/// Position in the file of each local entry, by global identity; `None`
/// when the order is the same.
fn mapping(file_ids: &[u64], local_ids: &[u64]) -> Result<Option<Vec<usize>>> {
    anyhow::ensure!(
        file_ids.len() == local_ids.len(),
        "the restart file has {} entries, the mesh {} (another mesh?)",
        file_ids.len(),
        local_ids.len()
    );
    if file_ids == local_ids {
        return Ok(None);
    }
    let pos: std::collections::HashMap<u64, usize> =
        file_ids.iter().enumerate().map(|(i, &g)| (g, i)).collect();
    let m: Option<Vec<usize>> = local_ids.iter().map(|g| pos.get(g).copied()).collect();
    m.map(Some)
        .ok_or_else(|| anyhow::anyhow!("the restart file belongs to another mesh"))
}

fn remap(v: &[f64], m: &[usize]) -> Vec<f64> {
    m.iter().map(|&i| v[i]).collect()
}

/// What a restart's arrays become: every array per cell, per face and per
/// boundary face and every matrix (in the sliced-ELL storage of the
/// discretisation) goes through one of the four maps, always in the same
/// order ([`Restart::map_fields`]). A reordering by identity
/// ([`Restart::to_local`]), the gather of the ranks' parts into the global
/// order and the scatter of a global file to a rank's part
/// ([`crate::ranks`], MPI step W6) are such maps; on ranks every rank
/// calls them in the same order, so they may communicate.
pub trait FieldMap {
    /// An array per cell.
    fn cells(&mut self, v: &[f64]) -> Result<Vec<f64>>;
    /// An array per face.
    fn faces(&mut self, v: &[f64]) -> Result<Vec<f64>>;
    /// An array per boundary face.
    fn boundary(&mut self, v: &[f64]) -> Result<Vec<f64>>;
    /// A matrix; `None` drops it (its hierarchy is then rebuilt).
    fn matrix(&mut self, v: &[f64]) -> Result<Option<Vec<f64>>>;
}

struct Reorder {
    mc: Vec<usize>,
    mf: Vec<usize>,
    mb: Vec<usize>,
}

impl FieldMap for Reorder {
    fn cells(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        Ok(remap(v, &self.mc))
    }
    fn faces(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        Ok(remap(v, &self.mf))
    }
    fn boundary(&mut self, v: &[f64]) -> Result<Vec<f64>> {
        Ok(remap(v, &self.mb))
    }
    fn matrix(&mut self, _: &[f64]) -> Result<Option<Vec<f64>>> {
        // the matrices follow the old cell and face order: rebuilt
        Ok(None)
    }
}

/// The fields of `s` mapped by `m`; with `wall` the boundary coefficients
/// of the scalars too (the output of a run on ranks needs them, a restart
/// file does not keep them).
pub(crate) fn map_state(s: &State, m: &mut dyn FieldMap, wall: bool) -> Result<State> {
    let opt = |m: &mut dyn FieldMap,
               v: &Option<Vec<f64>>,
               f: fn(&mut dyn FieldMap, &[f64]) -> Result<Vec<f64>>|
     -> Result<Option<Vec<f64>>> { v.as_ref().map(|v| f(m, v)).transpose() };
    let c3 = |m: &mut dyn FieldMap, v: &Soa3| -> Result<Soa3> {
        Ok(Soa3 {
            x: m.cells(&v.x)?,
            y: m.cells(&v.y)?,
            z: m.cells(&v.z)?,
        })
    };
    let cells: fn(&mut dyn FieldMap, &[f64]) -> Result<Vec<f64>> = |m, v| m.cells(v);
    let bnd: fn(&mut dyn FieldMap, &[f64]) -> Result<Vec<f64>> = |m, v| m.boundary(v);
    Ok(State {
        u: c3(m, &s.u)?,
        p: m.cells(&s.p)?,
        mass_flux: m.faces(&s.mass_flux)?,
        ub: Soa3 {
            x: m.boundary(&s.ub.x)?,
            y: m.boundary(&s.ub.y)?,
            z: m.boundary(&s.ub.z)?,
        },
        pb: m.boundary(&s.pb)?,
        grad_u: [
            c3(m, &s.grad_u[0])?,
            c3(m, &s.grad_u[1])?,
            c3(m, &s.grad_u[2])?,
        ],
        grad_p: c3(m, &s.grad_p)?,
        psi: c3(m, &s.psi)?,
        scalars: s
            .scalars
            .iter()
            .map(|x| {
                Ok(crate::scalar::ScalarState {
                    name: x.name.clone(),
                    phi: m.cells(&x.phi)?,
                    old: opt(m, &x.old, cells)?,
                    psi: opt(m, &x.psi, cells)?,
                    wall: match (&x.wall, wall) {
                        (Some(w), true) => Some(crate::scalar::ScalarWall {
                            a: m.boundary(&w.a)?,
                            b: m.boundary(&w.b)?,
                            gamma: m.boundary(&w.gamma)?,
                        }),
                        _ => None,
                    },
                })
            })
            .collect::<Result<_>>()?,
        periodic_force: s.periodic_force.clone(),
        turbulence: match &s.turbulence {
            None => None,
            Some(t) => Some(crate::turbulence::TurbState {
                epsilon: t.epsilon,
                k: m.cells(&t.k)?,
                omega: m.cells(&t.omega)?,
                nu_t: m.cells(&t.nu_t)?,
                k_old: opt(m, &t.k_old, cells)?,
                omega_old: opt(m, &t.omega_old, cells)?,
                mu_eff: m.faces(&t.mu_eff)?,
                k_b: m.boundary(&t.k_b)?,
                omega_b: m.boundary(&t.omega_b)?,
                u_tau: m.boundary(&t.u_tau)?,
                y_plus: m.boundary(&t.y_plus)?,
                wall_blending: opt(m, &t.wall_blending, bnd)?,
                p_modified: m.cells(&t.p_modified)?,
                pb_modified: m.boundary(&t.pb_modified)?,
                wall_distance: m.cells(&t.wall_distance)?,
                zeta: match &t.zeta {
                    None => None,
                    Some(z) => Some(crate::turbulence::ZetaState {
                        zeta: m.cells(&z.zeta)?,
                        f: m.cells(&z.f)?,
                        zeta_old: opt(m, &z.zeta_old, cells)?,
                        zeta_b: m.boundary(&z.zeta_b)?,
                        f_b: m.boundary(&z.f_b)?,
                        t: m.cells(&z.t)?,
                    }),
                },
                amg: {
                    let mut kept = Vec::new();
                    for (slot, a) in &t.amg {
                        if let Some(a) = m.matrix(a)? {
                            kept.push((*slot, a));
                        }
                    }
                    kept
                },
                fr1: opt(m, &t.fr1, cells)?,
            }),
        },
    })
}

impl Restart {
    /// The content in the local order of a mesh with the identities
    /// `cell_id`, `face_id` (n_interior interior faces); false in the second
    /// value when the order differed (the multigrid matrix is then dropped).
    pub fn to_local(
        self,
        cell_id: &[u64],
        face_id: &[u64],
        n_interior: usize,
    ) -> Result<(Self, bool)> {
        let n_cells = cell_id.len();
        let mc = mapping(&self.cell_id, cell_id)?;
        let mf = mapping(&self.face_id, face_id)?;
        if mc.is_none() && mf.is_none() {
            return Ok((self, true));
        }
        let ident = |n: usize| (0..n).collect::<Vec<usize>>();
        let mc = mc.unwrap_or_else(|| ident(n_cells));
        let mf = mf.unwrap_or_else(|| ident(face_id.len()));
        anyhow::ensure!(
            mf[n_interior..].iter().all(|&i| i >= n_interior),
            "the restart file has another split of interior and boundary faces"
        );
        // boundary faces: the positions among the boundary faces
        let mb: Vec<usize> = mf[n_interior..].iter().map(|&i| i - n_interior).collect();
        let r = self.map_fields(&mut Reorder { mc, mf, mb }, cell_id, face_id)?;
        Ok((r, false))
    }

    /// Every array and matrix through `m` (see [`FieldMap`]); the result
    /// carries the identities `cell_id`, `face_id`. A matrix that `m`
    /// drops leaves its hierarchy to be rebuilt. The temperature of a
    /// conjugate case (on its own full mesh) and the points of a moving
    /// mesh are not mapped: the caller refuses them where they cannot be.
    pub fn map_fields(
        mut self,
        m: &mut dyn FieldMap,
        cell_id: &[u64],
        face_id: &[u64],
    ) -> Result<Self> {
        let state = map_state(&self.state, m, false)?;
        self.state = state;
        if let Some((o, ..)) = &mut self.transient {
            let c3a = |m: &mut dyn FieldMap, u: &[Vec<f64>; 3]| -> Result<[Vec<f64>; 3]> {
                Ok([m.cells(&u[0])?, m.cells(&u[1])?, m.cells(&u[2])?])
            };
            *o = OldLevels {
                u_old: c3a(m, &o.u_old)?,
                m_old: m.faces(&o.m_old)?,
                x_cur: m.faces(&o.x_cur)?,
                x_old: m.faces(&o.x_old)?,
                g_old: o.g_old.iter().map(|v| m.cells(v)).collect::<Result<_>>()?,
                u_old2: o.u_old2.as_ref().map(|u| c3a(m, u)).transpose()?,
                ub_old: match &o.ub_old {
                    None => None,
                    Some(u) => Some([m.boundary(&u[0])?, m.boundary(&u[1])?, m.boundary(&u[2])?]),
                },
            };
        }
        if let Some(geo) = &mut self.geometry {
            geo.volumes = geo
                .volumes
                .iter()
                .map(|v| m.cells(v))
                .collect::<Result<_>>()?;
            geo.swept = geo
                .swept
                .iter()
                .map(|v| m.faces(v))
                .collect::<Result<_>>()?;
        }
        self.amg_matrix = opt_matrix(m, &self.amg_matrix)?;
        // the momentum hierarchies (D179): `matrix` per hierarchy, `slot`
        // the hierarchy of each velocity component
        let mut sweep = crate::coupled::SweepState::default();
        for h in 0..3 {
            sweep.matrix[h] = opt_matrix(m, &self.sweep.matrix[h])?;
        }
        for i in 0..3 {
            sweep.slot[i] = self.sweep.slot[i].filter(|&h| sweep.matrix[h].is_some());
        }
        self.sweep = sweep;
        self.cell_id = cell_id.to_vec();
        self.face_id = face_id.to_vec();
        Ok(self)
    }
}

fn opt_matrix(m: &mut dyn FieldMap, a: &Option<Vec<f64>>) -> Result<Option<Vec<f64>>> {
    match a {
        Some(a) => m.matrix(a),
        None => Ok(None),
    }
}
