//! Torque, shaft power, angular momentum and the performance of rotating
//! machinery (D184), for the moving reference frame (`[[mrf]]`, D55) and the
//! sliding mesh (`[[sliding]]`, D182) alike.
//!
//! # Torque and shaft power
//!
//! The force of the fluid on a wall face is **F**_f = p_f **S**_f +
//! **F**_f^visc, as in the boundary forces, and its torque about the axis of
//! a turning zone T_f = **â**·((**x**_f − **x**₀) × **F**_f), with **â**
//! the unit vector along **Ω** (a positive torque turns the wall the way the
//! zone turns). The **rotor** of a zone is a set of walls: the `walls` of
//! its `[[mrf]]` or `[[sliding]]` table, by default every wall all of whose
//! faces lie in the zone. Its torque T is the sum over their faces, and its
//! **shaft power** is the power it gives the fluid,
//!
//! P = −**Ω**·**T** = −|Ω| T,
//!
//! positive for a pump, a fan, a compressor or a mixer and negative for a
//! turbine. For walls that turn with the zone, **u**_b = **Ω** × (**x**_f −
//! **x**₀), it is the work −Σ **F**_f·**u**_b of the walls on the fluid,
//! the work that the kinetic-energy ledger (D88) books for them; there it
//! is the shaft work of the rotor. `ROT001` warns when a wall of a rotor
//! does not turn with its zone (a stationary shroud in a rotating zone):
//! its torque is in that of the rotor, but it does no work.
//!
//! Every other wall gets its torque about the axis of the first turning
//! zone (the reaction torque of a stator, a volute, the baffles).
//!
//! # Angular momentum
//!
//! About the axis of each rotor, the torque of all boundaries on the fluid
//! balances the angular momentum that flows out, Σ ṁ_f^r **â**·((**x**_f −
//! **x**₀) × **u**_b), in a steady state; ṁ^r is the convecting flux (in a
//! turning zone relative to it, zero on a wall that turns with it). The
//! discrete equations conserve linear momentum, not angular momentum, so
//! the imbalance is a discretisation error that falls with the mesh. In a
//! closed vessel the torque on the stator balances that on the rotor; in a
//! machine the outflow is the Euler work divided by |Ω|. A transient run
//! has the storage of angular momentum in the imbalance; over a revolution
//! of a periodic state it averages out.
//!
//! # Performance of a machine
//!
//! `[[output.machines]]` names a rotor, its inlet and its outlet. With the
//! total pressure p_t = p + ½ρ|**u**|² (absolute velocity), mass-flow
//! weighted over a boundary list, and the static pressure p area weighted:
//!
//! - mass flow ṁ out of the outlet and volume flow Q = ṁ/ρ;
//! - total pressure rise Δp_t = p̄_t,out − p̄_t,in and total-to-static rise
//!   Δp_ts = p̄_out − p̄_t,in; head H = Δp_t / (ρ g);
//! - hydraulic power P_h = Σ_out p_t ṁ_f/ρ + Σ_in p_t ṁ_f/ρ (the flux of
//!   total pressure out of the domain; Q Δp_t when the flows balance);
//! - shaft power P of the rotor; hydraulic efficiency η = P_h / P for a
//!   machine that takes power (P > 0), η = P / P_h for a turbine, and the
//!   total-to-static efficiency with Q Δp_ts instead of P_h;
//! - Euler work P_E = ṁ |Ω| ((r c_u)_out − (r c_u)_in) with r c_u =
//!   **â**·((**x** − **x**₀) × **u**) mass-flow weighted (Euler's equation of
//!   turbomachinery: the work of the rotor when no other wall turns the
//!   flow), the Euler head H_E = P_E / (ṁ g), and the losses P − P_h.
//!
//! All values are those of the mesh: a two-dimensional case gives them per
//! its depth. Every sum is a pairwise sum in a fixed order (D31).

use crate::bc::{BoundaryCondition, Turning};
use crate::case::{FlowCase, wall_face_transposed, wall_face_values};
use crate::diag::{self, Diagnostic};
use crate::simplec::Solution;
use foehncraft_core::exact::exact_sum;
use foehncraft_mesh::geometry::{cross, dot, norm, sub};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// Standard gravity (m/s²), the default of [`MachineSection::gravity`].
pub const STANDARD_GRAVITY: f64 = 9.80665;

fn standard_gravity() -> f64 {
    STANDARD_GRAVITY
}

/// `[[output.machines]]`: the performance of a pump, fan or turbine (D184).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MachineSection {
    /// Name (report, `forces.csv`, log columns `machine.<name>.<quantity>`).
    pub name: String,
    /// The turning zone of the rotor (`zone` of an `[[mrf]]` or `[[sliding]]`
    /// table); default: the only turning zone.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rotor: Option<String>,
    /// The boundaries of the inlet (the suction side of a pump).
    pub inlet: Vec<String>,
    /// The boundaries of the outlet.
    pub outlet: Vec<String>,
    /// Gravity of the head H = Δp_t / (ρ g) (m/s²).
    #[serde(default = "standard_gravity")]
    pub gravity: f64,
}

/// The rotor of a turning zone on a mesh.
#[derive(Debug, Clone)]
pub struct Rotor {
    /// Name: the zone, or `rotor` when every cell turns.
    pub name: String,
    /// The zone.
    pub turning: Turning,
    /// Its walls (indices into `mesh.boundaries`).
    pub walls: Vec<usize>,
}

/// The name of the rotor of a zone.
pub fn rotor_name(zone: &Option<String>) -> String {
    zone.clone().unwrap_or_else(|| "rotor".into())
}

fn is_wall(case: &FlowCase, name: &str) -> bool {
    matches!(
        case.boundary.get(name),
        Some(BoundaryCondition::Wall { .. })
    )
}

/// The rotors of `case` on `mesh`, one per `[[mrf]]` and `[[sliding]]`
/// table, in that order (as [`FlowCase::turning`]).
pub fn rotors(case: &FlowCase, mesh: &Mesh) -> Vec<Rotor> {
    let turning = case.turning();
    let entries: Vec<(&Option<String>, &Vec<String>)> = case
        .mrf
        .iter()
        .map(|m| (&m.zone, &m.walls))
        .chain(case.sliding.iter().map(|s| (&s.zone, &s.walls)))
        .collect();
    if turning.len() != entries.len() {
        // an entry without a valid speed: validation reports it
        return Vec::new();
    }
    turning
        .into_iter()
        .zip(entries)
        .map(|(t, (zone, walls))| {
            let ws = mesh
                .boundaries
                .iter()
                .enumerate()
                .filter(|(_, p)| {
                    if walls.is_empty() {
                        is_wall(case, &p.name)
                            && p.len > 0
                            && p.range().all(|f| t.contains(mesh, mesh.owner[f] as usize))
                    } else {
                        walls.contains(&p.name)
                    }
                })
                .map(|(i, _)| i)
                .collect();
            Rotor {
                name: rotor_name(zone),
                turning: t,
                walls: ws,
            }
        })
        .collect()
}

/// Checks `walls` of the turning zones and `[[output.machines]]` against
/// the boundaries of the case.
pub fn validate(case: &FlowCase) -> anyhow::Result<()> {
    let lists = case
        .mrf
        .iter()
        .map(|m| ("mrf", &m.walls))
        .chain(case.sliding.iter().map(|s| ("sliding", &s.walls)));
    for (table, walls) in lists {
        for w in walls {
            crate::need!(
                is_wall(case, w),
                ROT_ENTRY,
                "[[{table}]] walls: `{w}` is not a wall of the case"
            );
        }
    }
    let names: Vec<String> = case
        .mrf
        .iter()
        .map(|m| rotor_name(&m.zone))
        .chain(case.sliding.iter().map(|s| rotor_name(&s.zone)))
        .collect();
    for (i, m) in case.output.machines.iter().enumerate() {
        crate::need!(
            !names.is_empty(),
            ROT_ENTRY,
            "[[output.machines]] `{}` needs a turning zone ([[mrf]] or [[sliding]])",
            m.name
        );
        match &m.rotor {
            Some(r) => crate::need!(
                names.contains(r),
                ROT_ENTRY,
                "[[output.machines]] `{}`: no turning zone `{r}` (zones: {})",
                m.name,
                names.join(", ")
            ),
            None => crate::need!(
                names.len() == 1,
                ROT_ENTRY,
                "[[output.machines]] `{}`: give `rotor`, the case has several turning zones ({})",
                m.name,
                names.join(", ")
            ),
        }
        crate::need!(
            !m.inlet.is_empty() && !m.outlet.is_empty(),
            ROT_ENTRY,
            "[[output.machines]] `{}` needs an inlet and an outlet",
            m.name
        );
        for b in m.inlet.iter().chain(&m.outlet) {
            crate::need!(
                case.boundary.contains_key(b),
                ROT_ENTRY,
                "[[output.machines]] `{}`: no boundary `{b}` in the case",
                m.name
            );
        }
        crate::need!(
            m.gravity > 0.0 && m.gravity.is_finite(),
            CASE_RANGE,
            "[[output.machines]] `{}`: gravity must be > 0",
            m.name
        );
        for o in &case.output.machines[i + 1..] {
            crate::need!(
                o.name != m.name,
                ROT_ENTRY,
                "[[output.machines]]: the name `{}` twice",
                m.name
            );
        }
    }
    Ok(())
}

/// The angular-momentum balance about the axis of a rotor (N m).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct AngularMomentum {
    /// Torque of the rotor walls on the fluid (−T).
    pub rotor: f64,
    /// Torque of the other walls on the fluid.
    pub walls: f64,
    /// Torque of the pressure on the other boundaries on the fluid.
    pub openings: f64,
    /// Angular momentum flowing out of the domain.
    pub outflow: f64,
    /// outflow − (rotor + walls + openings): a discretisation error in a
    /// steady state (and the storage in a transient one).
    pub imbalance: f64,
    /// The imbalance over the largest term (-).
    pub relative: f64,
}

/// Torque and power of one rotor.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct RotorValues {
    /// Name.
    pub name: String,
    /// |Ω| (rad/s).
    pub omega: f64,
    /// Its walls.
    pub walls: Vec<String>,
    /// Torque of the fluid on its walls about the axis along Ω (N m;
    /// negative for a pump).
    pub torque: f64,
    /// Shaft power into the fluid, −|Ω| T (W; negative for a turbine).
    pub power: f64,
    /// Work of its walls on the fluid, −Σ **F**·**u**_b (W): the shaft
    /// power when the walls turn with the zone.
    pub wall_work: f64,
    /// Of it the work of the transposed viscous stress (W): the part that
    /// the kinetic-energy ledger, which has the Laplacian form of the
    /// viscous term, does not book at the wall (it is part of the
    /// dissipation of the full stress).
    pub transposed_work: f64,
    /// Largest |**u**_b − **Ω** × (**x** − **x**₀)| on its walls over
    /// |Ω| times the largest distance from the axis (-): zero when the
    /// walls turn with the zone (ROT001).
    pub slip: f64,
    /// Angular momentum about its axis.
    pub angular_momentum: AngularMomentum,
}

/// The performance of one machine of `[[output.machines]]`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct MachineValues {
    /// Name.
    pub name: String,
    /// Its rotor.
    pub rotor: String,
    /// Mass flow out of the outlet (kg/s).
    pub mass_flow: f64,
    /// Volume flow (m³/s).
    pub volume_flow: f64,
    /// Mass-flow weighted total pressure of the inlet (Pa).
    pub inlet_total_pressure: f64,
    /// Mass-flow weighted total pressure of the outlet (Pa).
    pub outlet_total_pressure: f64,
    /// Area weighted static pressure of the outlet (Pa).
    pub outlet_static_pressure: f64,
    /// Δp_t (Pa).
    pub total_pressure_rise: f64,
    /// Δp_ts (Pa).
    pub total_to_static_pressure_rise: f64,
    /// H = Δp_t / (ρ g) (m).
    pub head: f64,
    /// Shaft power of the rotor (W).
    pub shaft_power: f64,
    /// Shaft torque, the torque that drives the rotor, −T (N m).
    pub shaft_torque: f64,
    /// Flux of total pressure out of the domain (W).
    pub hydraulic_power: f64,
    /// Hydraulic (total-to-total) efficiency (-).
    pub efficiency: f64,
    /// Total-to-static efficiency (-).
    pub total_to_static_efficiency: f64,
    /// Euler work ṁ |Ω| Δ(r c_u) (W).
    pub euler_work: f64,
    /// Euler head (m).
    pub euler_head: f64,
    /// Shaft power − hydraulic power (W).
    pub losses: f64,
}

/// The values of all rotors and machines of one state.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Values {
    /// One per turning zone.
    pub rotors: Vec<RotorValues>,
    /// One per `[[output.machines]]` entry.
    pub machines: Vec<MachineValues>,
    /// Torque of the fluid on every wall about the axis of its rotor, or of
    /// the first turning zone (N m), by boundary index; `None` for the other
    /// boundaries.
    #[serde(skip)]
    pub wall_torque: Vec<Option<f64>>,
}

impl Values {
    /// The rotor `name`.
    pub fn rotor(&self, name: &str) -> Option<&RotorValues> {
        self.rotors.iter().find(|r| r.name == name)
    }

    /// The machine `name`.
    pub fn machine(&self, name: &str) -> Option<&MachineValues> {
        self.machines.iter().find(|m| m.name == name)
    }
}

/// The frame flux ρ **S**·**V** of boundary face `f` of `mesh` (zero outside
/// the turning zones): the convecting flux there is ṁ_f minus it.
pub fn frame_flux(turning: &[Turning], mesh: &Mesh, geo: &Geometry, f: usize, rho: f64) -> f64 {
    let c = mesh.owner[f] as usize;
    match Turning::of(turning, mesh, c) {
        Some(t) => rho * dot(geo.face_area[f], t.at(geo.face_centroid[f])),
        None => 0.0,
    }
}

/// One boundary face: the force of the fluid on it, its velocity, absolute
/// and convecting mass flux, pressure.
struct Face {
    x: [f64; 3],
    force: [f64; 3],
    /// Of it the transposed part of the viscous stress (a moving wall).
    transposed: [f64; 3],
    ub: [f64; 3],
    flux: f64,
    conv: f64,
    p: f64,
    area: f64,
}

/// The torques, powers and machine values of the state `sol` (D184).
pub fn evaluate(sol: &Solution, mesh: &Mesh, geo: &Geometry, case: &FlowCase) -> Values {
    let rotors = rotors(case, mesh);
    if rotors.is_empty() {
        return Values::default();
    }
    let ni = mesh.n_interior_faces();
    let rho = sol.phys.rho;
    let turning: Vec<Turning> = rotors.iter().map(|r| r.turning.clone()).collect();
    // the faces of every boundary but the periodic ones (those are pairs of
    // the interior)
    let faces: Vec<Vec<Face>> = mesh
        .boundaries
        .iter()
        .map(|p| {
            let bc = &case.boundary[&p.name];
            if matches!(bc, BoundaryCondition::Periodic { .. }) {
                return Vec::new();
            }
            let wall = matches!(bc, BoundaryCondition::Wall { .. });
            p.range()
                .map(|f| {
                    let s = geo.face_area[f];
                    let (pf, fv) = if wall {
                        let (pf, fv, _) = wall_face_values(sol, geo, f, ni);
                        (pf, fv)
                    } else {
                        (sol.pb[f - ni], [0.0; 3])
                    };
                    let tr = if wall {
                        wall_face_transposed(sol, geo, f, ni)
                    } else {
                        [0.0; 3]
                    };
                    let flux = sol.mass_flux[f];
                    Face {
                        x: geo.face_centroid[f],
                        force: [
                            pf * s[0] + fv[0] + tr[0],
                            pf * s[1] + fv[1] + tr[1],
                            pf * s[2] + fv[2] + tr[2],
                        ],
                        transposed: tr,
                        ub: sol.ub.get(f - ni),
                        flux,
                        conv: flux - frame_flux(&turning, mesh, geo, f, rho),
                        p: pf,
                        area: norm(s),
                    }
                })
                .collect()
        })
        .collect();
    let walls: Vec<bool> = mesh
        .boundaries
        .iter()
        .map(|p| is_wall(case, &p.name))
        .collect();
    let torque = |t: &Turning, fs: &[Face]| -> f64 {
        exact_sum(
            &fs.iter()
                .map(|f| dot(t.axis, cross(sub(f.x, t.origin), f.force)))
                .collect::<Vec<_>>(),
        )
    };
    let mut wall_torque = vec![None; mesh.boundaries.len()];
    for (b, fs) in faces.iter().enumerate() {
        if !walls[b] {
            continue;
        }
        let r = rotors
            .iter()
            .find(|r| r.walls.contains(&b))
            .unwrap_or(&rotors[0]);
        wall_torque[b] = Some(torque(&r.turning, fs));
    }
    let values: Vec<RotorValues> = rotors
        .iter()
        .map(|r| {
            let t = &r.turning;
            let w = norm(t.omega);
            let mut on_rotor = Vec::new();
            let mut on_walls = Vec::new();
            let mut on_open = Vec::new();
            let mut out = Vec::new();
            let mut work = Vec::new();
            let mut transposed = Vec::new();
            let (mut slip, mut reach) = (0.0f64, 0.0f64);
            for (b, fs) in faces.iter().enumerate() {
                let own = r.walls.contains(&b);
                for f in fs {
                    let arm = sub(f.x, t.origin);
                    let m = dot(t.axis, cross(arm, f.force));
                    if own {
                        on_rotor.push(m);
                        work.push(-dot(f.force, f.ub));
                        transposed.push(-dot(f.transposed, f.ub));
                        let v = t.at(f.x);
                        slip = slip.max(norm(sub(f.ub, v)));
                        reach = reach.max(norm(v));
                    } else if walls[b] {
                        on_walls.push(-m);
                    } else {
                        // the pressure only: the viscous stress of an opening
                        // is small and that of a symmetry plane has no torque
                        on_open.push(-m);
                    }
                    if f.conv != 0.0 {
                        out.push(f.conv * dot(t.axis, cross(arm, f.ub)));
                    }
                }
            }
            let tq = exact_sum(&on_rotor);
            let mut am = AngularMomentum {
                rotor: -tq,
                walls: exact_sum(&on_walls),
                openings: exact_sum(&on_open),
                outflow: exact_sum(&out),
                ..Default::default()
            };
            am.imbalance = am.outflow - (am.rotor + am.walls + am.openings);
            let big = am
                .rotor
                .abs()
                .max(am.walls.abs())
                .max(am.openings.abs())
                .max(am.outflow.abs());
            am.relative = if big > 0.0 { am.imbalance / big } else { 0.0 };
            RotorValues {
                name: r.name.clone(),
                omega: w,
                walls: r
                    .walls
                    .iter()
                    .map(|&b| mesh.boundaries[b].name.clone())
                    .collect(),
                torque: tq,
                power: -w * tq,
                wall_work: exact_sum(&work),
                transposed_work: exact_sum(&transposed),
                slip: if reach > 0.0 { slip / reach } else { 0.0 },
                angular_momentum: am,
            }
        })
        .collect();
    let index = |names: &[String]| -> Vec<usize> {
        mesh.boundaries
            .iter()
            .enumerate()
            .filter(|(_, p)| names.contains(&p.name))
            .map(|(i, _)| i)
            .collect()
    };
    let machines = case
        .output
        .machines
        .iter()
        .filter_map(|m| {
            let ri = match &m.rotor {
                Some(n) => rotors.iter().position(|r| &r.name == n)?,
                None => 0,
            };
            let (r, rv) = (&rotors[ri], &values[ri]);
            let t = &r.turning;
            let w = rv.omega;
            let side = |bs: &[usize]| {
                let (mut m, mut pt, mut rcu, mut pa, mut a) =
                    (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
                for f in bs.iter().flat_map(|&b| &faces[b]) {
                    let q = 0.5 * rho * dot(f.ub, f.ub);
                    m.push(f.flux);
                    pt.push((f.p + q) * f.flux);
                    rcu.push(dot(t.axis, cross(sub(f.x, t.origin), f.ub)) * f.flux);
                    a.push(f.area);
                    pa.push(f.p * f.area);
                }
                (
                    exact_sum(&m),
                    exact_sum(&pt),
                    exact_sum(&rcu),
                    exact_sum(&pa),
                    exact_sum(&a),
                )
            };
            let (m_in, pt_in, l_in, _, _) = side(&index(&m.inlet));
            let (m_out, pt_out, l_out, pa_out, a_out) = side(&index(&m.outlet));
            let mean = |v: f64, by: f64| if by != 0.0 { v / by } else { f64::NAN };
            let p_in = mean(pt_in, m_in);
            let p_out = mean(pt_out, m_out);
            let ps_out = mean(pa_out, a_out);
            let dpt = p_out - p_in;
            let dpts = ps_out - p_in;
            let q = m_out / rho;
            let hydraulic = (pt_out + pt_in) / rho;
            let power = rv.power;
            let eff = |useful: f64| {
                if power > 0.0 {
                    useful / power
                } else if power < 0.0 && useful != 0.0 {
                    power / useful
                } else {
                    f64::NAN
                }
            };
            // Δ(r c_u) between the mass-flow weighted means
            let euler = m_out * w * (mean(l_out, m_out) - mean(l_in, m_in));
            Some(MachineValues {
                name: m.name.clone(),
                rotor: r.name.clone(),
                mass_flow: m_out,
                volume_flow: q,
                inlet_total_pressure: p_in,
                outlet_total_pressure: p_out,
                outlet_static_pressure: ps_out,
                total_pressure_rise: dpt,
                total_to_static_pressure_rise: dpts,
                head: dpt / (rho * m.gravity),
                shaft_power: power,
                shaft_torque: -rv.torque,
                hydraulic_power: hydraulic,
                efficiency: eff(hydraulic),
                total_to_static_efficiency: eff(q * dpts),
                euler_work: euler,
                euler_head: mean(euler, m_out * m.gravity),
                losses: power - hydraulic,
            })
        })
        .collect();
    Values {
        rotors: values,
        machines,
        wall_torque,
    }
}

/// `ROT001`: the walls of a rotor that do not turn with its zone.
pub fn check(v: &Values) -> Vec<Diagnostic> {
    v.rotors
        .iter()
        .filter(|r| r.slip > 1e-6)
        .map(|r| {
            Diagnostic::new(
                diag::ROT_SLIP,
                format!(
                    "a wall of rotor `{}` ({}) does not turn with its zone (largest slip {:.2e} of \
                     the zone speed): its torque is in that of the rotor, but the shaft power \
                     −Ω·T is not its work on the fluid; give the walls of the rotor in `walls`",
                    r.name,
                    r.walls.join(", "),
                    r.slip
                ),
            )
        })
        .collect()
}

/// The rotors and machines in the report (`flow.md`).
pub fn markdown(v: &Values) -> String {
    use std::fmt::Write;
    let mut s = String::from(
        "\n| Rotor | Walls | |Ω| (rad/s) | Torque (N m) | Shaft power (W) | Angular momentum imbalance |\n|---|---|---|---|---|---|\n",
    );
    for r in &v.rotors {
        let _ = writeln!(
            s,
            "| {} | {} | {:.6e} | {:.6e} | {:.6e} | {:.2e} |",
            r.name,
            r.walls.join(", "),
            r.omega,
            r.torque,
            r.power,
            r.angular_momentum.relative
        );
    }
    if !v.machines.is_empty() {
        s.push_str(
            "\n| Machine | Q (m³/s) | Δp_t (Pa) | Head (m) | Shaft power (W) | Hydraulic power (W) | η | η_ts | Euler work (W) |\n|---|---|---|---|---|---|---|---|---|\n",
        );
        for m in &v.machines {
            let _ = writeln!(
                s,
                "| {} | {:.5e} | {:.5e} | {:.5e} | {:.5e} | {:.5e} | {:.4} | {:.4} | {:.5e} |",
                m.name,
                m.volume_flow,
                m.total_pressure_rise,
                m.head,
                m.shaft_power,
                m.hydraulic_power,
                m.efficiency,
                m.total_to_static_efficiency,
                m.euler_work
            );
        }
    }
    s
}
