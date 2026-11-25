//! The sliding mesh with conformal stepping (D56, D182, variant A): a zone
//! of the mesh (the rotor) that turns about an axis through the time steps
//! of a transient run, past the cells of the rest of the mesh (the stator).
//! Case file: one `[[sliding]]` table.
//!
//! # The method
//!
//! The rotor mesh turns **rigidly**: it keeps its cells, volumes and faces,
//! so the geometric conservation law holds exactly. When it turns by a whole
//! number k of its sectors per time step (a rotor of N identical sectors,
//! Δt = k · 2π / (N |Ω|)), every rotor cell lands on the place of another
//! rotor cell, and the cells of the interface meet new neighbours without
//! any interpolation: the mesh stays **conformal**. The solver keeps the
//! mesh where it is: the fields are the fields at the places of the mesh
//! (Eulerian), also the old levels, from which the gradients, the
//! Rhie–Chow differences and the extrapolations are made. The only thing
//! that sees the motion is the **time term**: the old level of a cell that
//! moved is what that cell held, which is the Eulerian old level at the
//! place that it came from, k sectors back for the level n and 2k for the
//! level n−1. So the time term of the momentum and scalar rows of a rotor
//! cell reads the old levels with the rotor cells shifted by k and 2k
//! sectors ([`Plan::cell`]); everything else is the solver of a fixed mesh.
//!
//! # Equations
//!
//! The solved velocity is the absolute one **u** in the inertial basis
//! (the vectors of the shifted old levels do not turn with their cells: the
//! pressure gradient turns them, as it must), the equations are those of a
//! moving mesh,
//!
//! d/dt ∫ ρ **u** dV + ∮ ρ (**u** − **w**)·**S** **u** = −∮ p **S** + ∮ τ·**S**,
//!
//! with the grid velocity **w** = **Ω** × (**x** − **x**₀) of the rotor and
//! **no Coriolis source**. The convecting flux is the relative flux ṁ −
//! ρ **S**·**w**, the frame flux of the moving reference frame
//! ([`crate::frame`]), which is the grid flux of a rigid turn: Σ_f ρ
//! **S**·**w** = 0 for every closed rotor cell (a rigid turn keeps volumes),
//! and ρ **S**·**w** = 0 on an interface that is a surface of revolution
//! about the axis. Continuity, the Rhie–Chow flux and the pressure
//! correction work on the absolute flux. Scalars and turbulence variables
//! are convected with the relative flux; a scalar that is carried by the
//! rotor (relative flux 0) is carried exactly by the shift, without
//! numerical diffusion.
//!
//! # Variant B: the interface by imprinting
//!
//! With `rotor_interface` and `stator_interface` the mesh holds the two parts
//! with different angular pitch (or the same one, offset by any angle) as
//! separate meshes whose faces on the interface are two boundaries. The time
//! step is free. In every step the mesh of the new rotor angle is made
//! ([`foehncraft_mesh::imprint`]): the rotor nodes are turned, every
//! intersection of a rotor face with a stator face becomes an ordinary
//! internal face between a rotor cell and a stator cell (mass and momentum
//! flux conserved to round-off, no interpolation), and the cells, their
//! volumes (a rigid turn) and their numbering are those of the base mesh.
//! The flow solver is built again on that mesh (`rebuild_flow` in
//! `unsteady.rs`) and takes over the cell fields, the old time levels and
//! the face fluxes of the faces that stay; the intersections are new faces
//! with the interpolated flux. The scheme is fully implicit in the step (the
//! extrapolation of the explicit terms off, several outer iterations), since
//! an explicit term made of the vectors of the previous angle would be
//! stale on the turned mesh (the lesson of variant A). This costs a rebuild
//! of the discretisation and the multigrid hierarchy per step; the
//! variant is for the geometry and the accuracy, not yet for the speed.
//!
//! # Limits
//!
//! - The rotor mesh must be a ring of identical sectors: it maps onto itself
//!   when turned by 2π/N about the axis (cells, faces, volumes). The
//!   sectors are found from the geometry (`SLD004` says what failed).
//! - The step is fixed and Δt = k · 2π / (N |Ω|) (`SLD003`).
//! - The boundary faces of the rotor must be walls that turn with it, or
//!   symmetry planes (`SLD004`); a wall with a prescribed velocity expression
//!   (D53) of the rigid turn does that.
//! - One rotor per case in conformal stepping; several rotors, each at its
//!   own speed and with its own interface to the stator, by imprinting
//!   (D185). An interface between two rotors is not available yet.
//!
//! # Turbulence
//!
//! The turbulence variables are transported like the passive scalars: with
//! the relative flux, and in conformal stepping with the old levels of the
//! time terms of k, ε or ω and ζ moved with the rotor cells. The wall
//! distance does not change in conformal stepping (the rotor lands on
//! identical sectors). With imprinting the cells keep their numbers, so the
//! cell fields and old levels go along with the rebuilt solver ([`carry`]),
//! and the wall distance of every angle comes from the walls of the two
//! parts kept apart ([`Imprinting::wall_distance`]: only cells nearer to
//! the edge of their part than to its walls search the other part).
//! - On the CPU (the shifted old levels are made through the host); on
//!   ranks (MPI) the plan is made on the whole mesh and the moved values
//!   of the own cells are taken from the cells of all ranks ([`Shift`]),
//!   in conformal stepping only.

use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::frame::{MrfSection, assign_cells_in};
use crate::simplec::{Flow, State};
use crate::transient::{OldLevels, Transient};
use anyhow::Result;
use foehncraft_fv::Discretisation;
use foehncraft_fv::disc::Soa3;
use foehncraft_mesh::imprint::{FaceOrigin, Imprint, Interface, InterfaceSpec};
use foehncraft_mesh::query::TurningWallDistance;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::f64::consts::TAU;

/// `[[sliding]]`: a zone that turns about an axis through the time steps
/// (D56, D182).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct SlidingSection {
    /// The cell zone of the mesh that turns (a zone name of the mesh;
    /// without it all cells turn).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zone: Option<String>,
    /// A point on the axis of rotation (m).
    pub origin: [f64; 3],
    /// Direction of the axis (any length; the rotation is right-handed about
    /// it for a positive speed).
    pub axis: [f64; 3],
    /// Angular velocity (rad/s). Give this or `rpm`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub omega: Option<f64>,
    /// Angular velocity in revolutions per minute. Give this or `omega`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rpm: Option<f64>,
    /// Sectors of the rotor that it turns in one time step (default 1); the
    /// time step is `shift` × 2π / (sectors × |Ω|).
    pub shift: usize,
    /// Variant B of D182, a pitch ratio that is not one: the boundary of the
    /// mesh that holds the faces of the rotor on the interface (with
    /// `stator_interface`; `[mesh.annulus]` makes them with
    /// `n_theta_stator`, `[mesh.revolved]` always). The interface is any
    /// surface of revolution about the axis whose meridian is made of
    /// straight pieces: a cylinder, a disc, a cone, or a chain of them, also
    /// through the axis (D186). The interface is imprinted at the rotor angle of
    /// every time step, which is free then.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rotor_interface: Option<String>,
    /// The boundary of the faces of the stator on the interface (with
    /// `rotor_interface`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stator_interface: Option<String>,
    /// The walls of the rotor, whose torque gives the shaft power (D184);
    /// default: every wall all of whose faces lie in the zone.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub walls: Vec<String>,
}

impl Default for SlidingSection {
    fn default() -> Self {
        Self {
            zone: None,
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            omega: None,
            rpm: None,
            shift: 1,
            rotor_interface: None,
            stator_interface: None,
            walls: Vec::new(),
        }
    }
}

impl SlidingSection {
    /// The rotating zone as a frame of [`crate::frame`] (the same axis,
    /// speed and zone).
    pub fn frame(&self) -> MrfSection {
        MrfSection {
            zone: self.zone.clone(),
            origin: self.origin,
            axis: self.axis,
            omega: self.omega,
            rpm: self.rpm,
            walls: self.walls.clone(),
        }
    }

    /// The interface is imprinted in every step (variant B).
    pub fn imprinted(&self) -> bool {
        self.rotor_interface.is_some() || self.stator_interface.is_some()
    }

    /// The angular velocity vector (rad/s).
    pub fn omega_vector(&self) -> Result<[f64; 3]> {
        self.frame().omega_vector()
    }

    /// The time step that turns the rotor by `shift` of `sectors` sectors.
    pub fn time_step(&self, sectors: usize) -> Result<f64> {
        let o = self.omega_vector()?;
        Ok(self.shift as f64 * TAU / (sectors as f64 * norm(o)))
    }
}

/// Relative deviation of the time step from `shift` × 2π / (N |Ω|) that is
/// still accepted (a time step given with a few digits).
const STEP_TOLERANCE: f64 = 1e-6;

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn err(code: &'static str, key: &str, msg: String) -> anyhow::Error {
    let mut d = Diagnostic::new(code, msg);
    d.key = Some(key.into());
    diag::err(d)
}

/// The checks of the `[[sliding]]` tables that need no mesh (D182).
pub fn validate(case: &FlowCase) -> Result<()> {
    if case.sliding.is_empty() {
        return Ok(());
    }
    let Some(t) = case.time else {
        return Err(err(
            diag::SLIDING_STEADY,
            "sliding",
            "[[sliding]] moves the mesh in time (D182): the case needs a [time] table; a steady \
             case with a rotating zone uses [[mrf]] (D55)"
                .into(),
        ));
    };
    let several = case.sliding.len() > 1;
    for (i, s) in case.sliding.iter().enumerate() {
        validate_entry(case, i, s, several)?;
    }
    crate::need!(
        !t.adaptive && t.error_tolerance.is_none(),
        SLIDING_STEPPING,
        "[[sliding]]: the time step is fixed (shift × 2π / (sectors × Ω)); `adaptive` and \
         `error_tolerance` are not available"
    );
    let v = &case.verification;
    crate::need!(
        v.levels.len() <= 1 && v.time_steps.is_empty(),
        SLIDING_STEPPING,
        "[[sliding]]: the time step follows from the mesh, so a sweep over levels or time steps \
         is a series of runs, not a [verification] table"
    );
    Ok(())
}

/// The checks of the `[[sliding]]` table `i` (`several`: the case has more
/// than one, D185).
fn validate_entry(case: &FlowCase, i: usize, s: &SlidingSection, several: bool) -> Result<()> {
    let key = format!("sliding[{i}]");
    let need = |ok: bool, code: &'static str, msg: String| -> Result<()> {
        if ok {
            Ok(())
        } else {
            Err(err(code, &key, msg))
        }
    };
    if let Err(e) = s.omega_vector() {
        return Err(err(diag::SLIDING_ENTRY, &key, format!("[[sliding]]: {e}")));
    }
    need(
        s.omega_vector().is_ok_and(|o| norm(o) > 0.0),
        diag::SLIDING_ENTRY,
        "[[sliding]]: the angular velocity must not be 0 (a rotor at rest needs no sliding mesh)"
            .into(),
    )?;
    need(
        s.shift >= 1,
        diag::SLIDING_ENTRY,
        "[[sliding]]: `shift` (sectors per time step) must be at least 1".into(),
    )?;
    need(
        s.rotor_interface.is_some() == s.stator_interface.is_some(),
        diag::SLIDING_ENTRY,
        "[[sliding]]: `rotor_interface` and `stator_interface` go together (the two boundaries \
         of the interface of variant B)"
            .into(),
    )?;
    need(
        !several || s.imprinted(),
        diag::SLIDING_ENTRY,
        format!(
            "[[sliding]] entry {}: several rotors (D185) need the interface by imprinting \
             (`rotor_interface` and `stator_interface` in every table): their time steps of \
             whole sectors do not agree",
            i + 1
        ),
    )?;
    if several {
        for (j, o) in case.sliding.iter().enumerate().skip(i + 1) {
            need(
                s.zone != o.zone,
                diag::SLIDING_ENTRY,
                format!(
                    "[[sliding]] entries {} and {} turn the zone `{}` twice",
                    i + 1,
                    j + 1,
                    s.zone.as_deref().unwrap_or("")
                ),
            )?;
        }
    }
    if s.imprinted() {
        need(
            s.zone.is_some(),
            diag::SLIDING_ENTRY,
            "[[sliding]]: `zone` (the cell zone of the rotor) is needed with `rotor_interface`"
                .into(),
        )?;
        need(
            s.shift == 1,
            diag::SLIDING_ENTRY,
            "[[sliding]]: `shift` belongs to conformal stepping; with `rotor_interface` the time \
             step is free"
                .into(),
        )?;
        need(
            !case.solver.extrapolation && case.solver.outer_iterations >= 2,
            diag::SLIDING_STEPPING,
            "[[sliding]]: the sliding mesh by imprinting is fully implicit in the step: set \
             [solver] extrapolation = false and outer_iterations of at least 2 (3 to 6 is usual), \
             since an explicit term made of the vectors of the previous rotor angle would be \
             stale on the turned mesh"
                .into(),
        )?;
        need(
            case.run.restart.is_none(),
            diag::SLIDING_INTERFACE,
            "[[sliding]]: a restart of the sliding mesh by imprinting is not available yet".into(),
        )?;
    }
    Ok(())
}

/// A rotor zone of the mesh and the shift of its cells per time step.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    /// Identical sectors of the rotor (N).
    pub sectors: usize,
    /// Sectors per step (k).
    pub shift: usize,
    /// The angle turned in one step (rad), in the sense of the rotation.
    pub angle: f64,
    /// Cells of the rotor.
    pub rotor_cells: usize,
    /// Entry i of the shifted array is entry `cell[i]` of the old one: the
    /// cell that is now at the place of cell i was at the place of cell
    /// `cell[i]` one step ago (identity outside the rotor).
    pub cell: Vec<usize>,
}

/// How the data of the rotor cells move by one and by two time steps
/// (D182): entry i of an array moved by one step is entry `perm1[i]` of the
/// array, `perm2` for two steps. On ranks (MPI) the permutations are those
/// of the whole mesh: a rotor cell may land on a cell of another rank, so
/// the own cells of every rank are gathered into the order of the whole
/// mesh and each rank takes the moved values of its own cells (copies, so
/// bitwise those of one rank).
#[derive(Clone)]
pub struct Shift {
    perm1: Vec<usize>,
    perm2: Vec<usize>,
    ranks: Option<ShiftRanks>,
}

#[derive(Clone)]
struct ShiftRanks {
    comm: foehncraft_mesh::decompose::SharedComm,
    /// The own cells of every rank as positions in the whole mesh.
    all_cells: Vec<Vec<usize>>,
    /// This rank's own cells as positions in the whole mesh.
    own: Vec<usize>,
}

impl std::fmt::Debug for Shift {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Shift")
            .field("cells", &self.perm1.len())
            .field("ranks", &self.ranks.as_ref().map(|r| r.comm.size()))
            .finish()
    }
}

impl Shift {
    /// The shift of the cells of one rank: `perm` for one step, applied
    /// twice for two.
    pub fn new(perm: Vec<usize>) -> Self {
        let perm2 = perm.iter().map(|&i| perm[i]).collect();
        Self {
            perm1: perm,
            perm2,
            ranks: None,
        }
    }

    /// The shift `perm` of the whole mesh on the rank of `comm`, whose own
    /// cells are at the positions `own` of the whole mesh; `all_cells`: the
    /// own cells of every rank (in rank order).
    pub fn on_ranks(
        perm: Vec<usize>,
        comm: foehncraft_mesh::decompose::SharedComm,
        all_cells: Vec<Vec<usize>>,
        own: Vec<usize>,
    ) -> Self {
        Self {
            ranks: Some(ShiftRanks {
                comm,
                all_cells,
                own,
            }),
            ..Self::new(perm)
        }
    }

    /// The number of cells the shift is made for: those of the solver on
    /// one rank, the own cells on ranks.
    pub fn cells(&self) -> usize {
        self.ranks
            .as_ref()
            .map_or(self.perm1.len(), |r| r.own.len())
    }

    /// The own-cell values `v` moved by `steps` (1 or 2) steps. Collective
    /// on ranks.
    pub fn apply(&self, v: &[f64], steps: usize) -> Vec<f64> {
        let perm = if steps == 2 { &self.perm2 } else { &self.perm1 };
        match &self.ranks {
            None => perm.iter().map(|&i| v[i]).collect(),
            Some(r) => {
                let all = r.comm.all_gather_vec_f64(v);
                let mut whole = vec![0.0; perm.len()];
                for (cells, vals) in r.all_cells.iter().zip(&all) {
                    for (&p, &x) in cells.iter().zip(vals) {
                        whole[p] = x;
                    }
                }
                r.own.iter().map(|&p| whole[perm[p]]).collect()
            }
        }
    }

    /// `dst` = the cell buffer `src` moved by `steps` steps (`dst` may be
    /// `src`; on ranks the halo of `dst` comes from its owners).
    pub(crate) fn moved(
        &self,
        dev: &mut dyn crate::device::FlowDevice,
        src: crate::device::Buf,
        dst: crate::device::Buf,
        steps: usize,
    ) -> Result<()> {
        let v = dev.download(src)?;
        dev.upload(dst, &self.apply(&v, steps))
    }
}

/// Points that can be found again within a tolerance.
struct PointIndex<'a> {
    pts: &'a [[f64; 3]],
    bin: f64,
    tol: f64,
    map: HashMap<[i64; 3], Vec<u32>>,
}

impl<'a> PointIndex<'a> {
    fn new(pts: &'a [[f64; 3]], tol: f64) -> Self {
        let bin = 16.0 * tol;
        let mut map: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
        for (i, p) in pts.iter().enumerate() {
            map.entry(Self::key(*p, bin)).or_default().push(i as u32);
        }
        Self { pts, bin, tol, map }
    }

    fn key(p: [f64; 3], bin: f64) -> [i64; 3] {
        [
            (p[0] / bin).floor() as i64,
            (p[1] / bin).floor() as i64,
            (p[2] / bin).floor() as i64,
        ]
    }

    /// The nearest point within the tolerance of `q`.
    fn find(&self, q: [f64; 3]) -> Option<usize> {
        let k = Self::key(q, self.bin);
        let mut best: Option<(f64, usize)> = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(v) = self.map.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) {
                        for &i in v {
                            let d = norm(sub(self.pts[i as usize], q));
                            if d <= self.tol && best.is_none_or(|(b, _)| d < b) {
                                best = Some((d, i as usize));
                            }
                        }
                    }
                }
            }
        }
        best.map(|(_, i)| i)
    }
}

/// Rodrigues' rotation of `v` about the unit axis `a` by `angle`.
fn rotate(a: [f64; 3], angle: f64, v: [f64; 3]) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    let axv = cross(a, v);
    let adv = dot(a, v);
    std::array::from_fn(|i| v[i] * c + axv[i] * s + a[i] * adv * (1.0 - c))
}

/// The shift of the data of the cells `rotor` of `disc` by the turn of
/// `entry` (see the module documentation).
pub fn plan(entry: &SlidingSection, rotor: &[bool], disc: &Discretisation) -> Result<Plan> {
    let o = entry.origin;
    let omega = entry.omega_vector()?;
    let w = norm(omega);
    let a = omega.map(|x| x / w);
    let nc = disc.n_cells;
    let cells: Vec<usize> = (0..nc).filter(|&c| rotor[c]).collect();
    let fail = |msg: String| err(diag::SLIDING_ZONE, "sliding", msg);
    if cells.is_empty() {
        return Err(fail("[[sliding]]: the rotor has no cells".into()));
    }
    let centre = |c: usize| disc.cell_centroid.get(c);
    // the size of the rotor and the tolerance of the matches
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for &c in &cells {
        let x = centre(c);
        for i in 0..3 {
            lo[i] = lo[i].min(x[i]);
            hi[i] = hi[i].max(x[i]);
        }
    }
    let size = norm(sub(hi, lo)).max(f64::MIN_POSITIVE);
    let tol = 1e-6 * size;
    // the sectors: the cells at the radius and height of the outermost cell
    let cyl = |x: [f64; 3]| {
        let d = sub(x, o);
        let z = dot(d, a);
        let r = norm(sub(d, a.map(|ai| ai * z)));
        (z, r)
    };
    let c0 = *cells
        .iter()
        .max_by(|&&p, &&q| cyl(centre(p)).1.total_cmp(&cyl(centre(q)).1))
        .expect("a cell");
    let (z0, r0) = cyl(centre(c0));
    let sectors = cells
        .iter()
        .filter(|&&c| {
            let (z, r) = cyl(centre(c));
            (z - z0).abs() <= tol && (r - r0).abs() <= tol
        })
        .count();
    if sectors < 2 || r0 <= tol {
        return Err(fail(format!(
            "[[sliding]]: the rotor is not a ring of identical sectors about the axis {:?} through \
             {:?}: no second cell lies at the radius and height of the outermost cell",
            entry.axis, entry.origin
        )));
    }
    if entry.shift.is_multiple_of(sectors) {
        return Err(err(
            diag::SLIDING_ENTRY,
            "sliding[0].shift",
            format!(
                "[[sliding]]: a shift of {} sectors is a whole turn of the {sectors} sectors of the rotor",
                entry.shift
            ),
        ));
    }
    let angle = TAU * entry.shift as f64 / sectors as f64;
    let turn = |x: [f64; 3]| {
        let r = rotate(a, angle, sub(x, o));
        [o[0] + r[0], o[1] + r[1], o[2] + r[2]]
    };
    // cells: the image of each rotor cell
    let pts: Vec<[f64; 3]> = cells.iter().map(|&c| centre(c)).collect();
    let index = PointIndex::new(&pts, tol);
    let mut image = vec![usize::MAX; nc];
    let mut cell: Vec<usize> = (0..nc).collect();
    let mut missing = 0usize;
    let mut first_missing = None;
    let mut taken = vec![false; nc];
    for &c in &cells {
        let target = turn(centre(c));
        match index.find(target) {
            Some(i) => {
                let t = cells[i];
                let dv = (disc.cell_volume[t] - disc.cell_volume[c]).abs();
                if dv > 1e-8 * disc.cell_volume[c] || taken[t] {
                    missing += 1;
                    first_missing.get_or_insert(c);
                } else {
                    taken[t] = true;
                    image[c] = t;
                    cell[t] = c;
                }
            }
            None => {
                missing += 1;
                first_missing.get_or_insert(c);
            }
        }
    }
    if missing > 0 {
        let c = first_missing.expect("a cell");
        return Err(fail(format!(
            "[[sliding]]: {missing} of {} rotor cells have no image when the rotor is turned by \
             {:.4}° (1/{sectors} of a turn, the first one is the cell at {:?}): the rotor mesh is not \
             a ring of {sectors} identical sectors about the axis {:?} through {:?}",
            cells.len(),
            angle.to_degrees(),
            centre(c),
            entry.axis,
            entry.origin
        )));
    }
    Ok(Plan {
        sectors,
        shift: entry.shift,
        angle,
        rotor_cells: cells.len(),
        cell,
    })
}

/// What [`attach`] reports for the log and the diagnostics.
#[derive(Debug, Default)]
pub struct SlidingNotes {
    /// Log lines.
    pub lines: Vec<String>,
}

/// The sliding mesh of a run.
#[derive(Debug, Clone)]
pub struct Sliding {
    /// The shift of the data.
    pub plan: Plan,
    /// The time step of the run (s).
    pub dt: f64,
}

/// The sliding mesh of `case` on `disc`, or `None` without `[[sliding]]`;
/// `dt` is the time step of the run.
pub fn attach(
    case: &FlowCase,
    mesh: &Mesh,
    disc: &Discretisation,
    dt: f64,
) -> Result<Option<(Sliding, SlidingNotes)>> {
    validate(case)?;
    let Some(entry) = case.sliding.first() else {
        return Ok(None);
    };
    if entry.imprinted() {
        return Ok(None);
    }
    if let Some(name) = &entry.zone
        && !mesh.zones.contains(name)
    {
        let hint = crate::bc::did_you_mean(name, mesh.zones.iter().map(String::as_str))
            .map(|c| format!(" (did you mean `{c}`?)"))
            .unwrap_or_default();
        return Err(err(
            diag::SLIDING_ENTRY,
            "sliding[0].zone",
            format!(
                "[[sliding]]: the mesh has no cell zone `{name}`{hint}; its zones are {:?}",
                mesh.zones
            ),
        ));
    }
    let cell_zone = assign_cells_in(&[entry.frame()], mesh, "[[sliding]]")?;
    let rotor: Vec<bool> = cell_zone.iter().map(Option::is_some).collect();
    // the boundaries of the rotor: walls that turn with it, or symmetry
    for b in &mesh.boundaries {
        if !b
            .range()
            .any(|f| rotor.get(disc.owner[f] as usize).copied().unwrap_or(false))
        {
            continue;
        }
        let Some(bc) = case.boundary.get(&b.name) else {
            continue;
        };
        use crate::bc::BoundaryCondition as Bc;
        if !matches!(
            bc,
            Bc::Wall { .. } | Bc::SlipWall | Bc::Symmetry | Bc::Wedge
        ) {
            return Err(err(
                diag::SLIDING_ZONE,
                &format!("boundary.{}", b.name),
                format!(
                    "[[sliding]]: the boundary `{}` of the rotor is a `{}`; the boundary faces of \
                     the rotor must be walls that turn with it, slip walls or symmetry planes",
                    b.name,
                    bc.type_name()
                ),
            ));
        }
    }
    let plan = plan(entry, &rotor, disc)?;
    let want = entry.time_step(plan.sectors)?;
    if (dt - want).abs() > STEP_TOLERANCE * want {
        return Err(err(
            diag::SLIDING_STEPPING,
            "time.time_step",
            format!(
                "[[sliding]]: the rotor turns {} of its {} sectors in one step, so the time step \
                 must be {want:.10e} s, not {dt:.10e} s",
                plan.shift, plan.sectors
            ),
        ));
    }
    if let Some(t) = case.time
        && let Some(e) = t.end_time
    {
        let n = (e - t.start_time) / dt;
        if (n - n.round()).abs() > 1e-6 {
            return Err(err(
                diag::SLIDING_STEPPING,
                "time.end_time",
                format!(
                    "[[sliding]]: `end_time` must be a whole number of steps after `start_time` \
                     ({n:.6} steps of {dt:.6e} s); give `steps`, or {:.10e} s",
                    t.start_time + n.round() * dt
                ),
            ));
        }
    }
    let notes = SlidingNotes {
        lines: vec![format!(
            "sliding mesh (D182): {} ({} cells) turns about the axis {:?} through {:?} at {:.6e} \
             rad/s; {} sectors, {} per step = {:.4}° per step, conformal stepping",
            entry
                .zone
                .as_deref()
                .map_or("all cells".to_string(), |z| format!("zone `{z}`")),
            plan.rotor_cells,
            entry.axis,
            entry.origin,
            norm(entry.omega_vector()?),
            plan.sectors,
            plan.shift,
            plan.angle.to_degrees()
        )],
    };
    Ok(Some((Sliding { plan, dt }, notes)))
}

/// The sliding mesh by imprinting (variant B of D182) of a run: one or
/// several rotors (D185), each with its interface to the stator.
#[derive(Debug, Clone)]
pub struct Imprinting {
    imp: Imprint,
    /// Angular velocity of every rotor about the unit axis of its interface
    /// (rad/s).
    omega: Vec<f64>,
}

/// The mesh of one set of rotor angles.
#[derive(Debug, Clone)]
pub struct Working {
    /// The conformal mesh (the cells of the base mesh).
    pub mesh: Mesh,
    /// Its geometry (the cells as rigid bodies).
    pub geo: Geometry,
    /// The origin of every face.
    pub origin: Vec<FaceOrigin>,
    /// The angle of every rotor (rad).
    pub angles: Vec<f64>,
    /// Intersection polygons (faces between rotor and stator).
    pub pieces: usize,
}

/// The sliding mesh by imprinting of `case` on the base mesh, or `None` for
/// conformal stepping or without `[[sliding]]`.
pub fn attach_imprint(case: &FlowCase, mesh: &Mesh) -> Result<Option<(Imprinting, SlidingNotes)>> {
    validate(case)?;
    if !case.sliding.first().is_some_and(SlidingSection::imprinted) {
        return Ok(None);
    }
    let fail = |key: &str, msg: String| err(diag::SLIDING_INTERFACE, key, msg);
    let mut parts = Vec::with_capacity(case.sliding.len());
    let mut omega = Vec::with_capacity(case.sliding.len());
    let mut notes = SlidingNotes::default();
    for (i, entry) in case.sliding.iter().enumerate() {
        let (rp, sp) = (
            entry.rotor_interface.clone().unwrap_or_default(),
            entry.stator_interface.clone().unwrap_or_default(),
        );
        let zone = entry.zone.clone().unwrap_or_default();
        if !mesh.zones.contains(&zone) {
            let hint = crate::bc::did_you_mean(&zone, mesh.zones.iter().map(String::as_str))
                .map(|c| format!(" (did you mean `{c}`?)"))
                .unwrap_or_default();
            return Err(err(
                diag::SLIDING_ENTRY,
                &format!("sliding[{i}].zone"),
                format!(
                    "[[sliding]]: the mesh has no cell zone `{zone}`{hint}; its zones are {:?}",
                    mesh.zones
                ),
            ));
        }
        for name in [&rp, &sp] {
            if case.boundary.contains_key(name) {
                return Err(fail(
                    &format!("boundary.{name}"),
                    format!(
                        "[[sliding]]: `{name}` is an interface boundary of the sliding mesh; it \
                         gets no [boundary.{name}] table (its faces are replaced by the imprint)"
                    ),
                ));
            }
        }
        let omega_vec = entry.omega_vector()?;
        let iface = Interface::new(
            mesh,
            &InterfaceSpec {
                zone: zone.clone(),
                rotor_patch: rp.clone(),
                stator_patch: sp.clone(),
                origin: entry.origin,
                axis: entry.axis,
            },
        )
        .map_err(|e| {
            fail(
                &format!("sliding[{i}].rotor_interface"),
                format!("[[sliding]]: {e}"),
            )
        })?;
        let w = dot(omega_vec, iface.axis());
        // the other boundaries of the rotor: walls that turn with it, or
        // symmetry
        let rotor = iface.rotor_cells();
        let patches = iface.patches();
        for (bi, b) in mesh.boundaries.iter().enumerate() {
            if patches.contains(&bi) || !b.range().any(|f| rotor[mesh.owner[f] as usize]) {
                continue;
            }
            let Some(bc) = case.boundary.get(&b.name) else {
                continue;
            };
            use crate::bc::BoundaryCondition as Bc;
            if !matches!(
                bc,
                Bc::Wall { .. } | Bc::SlipWall | Bc::Symmetry | Bc::Wedge
            ) {
                return Err(err(
                    diag::SLIDING_ZONE,
                    &format!("boundary.{}", b.name),
                    format!(
                        "[[sliding]]: the boundary `{}` of the rotor `{zone}` is a `{}`; the \
                         boundary faces of a rotor must be walls that turn with it, slip walls \
                         or symmetry planes",
                        b.name,
                        bc.type_name()
                    ),
                ));
            }
        }
        let (nr, ns) = iface.faces();
        notes.lines.push(format!(
            "sliding mesh (D182, variant B): zone `{zone}` turns about the axis {:?} through {:?} \
             at {:.6e} rad/s; the interface ({}) has {nr} rotor and {ns} stator faces \
             and is imprinted in every step (the solver is built again on the mesh of the angle, \
             fully implicit)",
            entry.axis,
            entry.origin,
            w,
            iface.surface()
        ));
        parts.push(iface);
        omega.push(w);
    }
    let n = parts.len();
    let imp =
        Imprint::new(mesh, parts).map_err(|e| fail("sliding", format!("[[sliding]]: {e}")))?;
    if n > 1 {
        notes.lines.push(format!(
            "sliding mesh (D185): {n} rotors, each at its own speed, imprinted together on the \
             stator in every step"
        ));
    }
    Ok(Some((Imprinting { imp, omega }, notes)))
}

impl Imprinting {
    /// The angle of every rotor at time `t` (rad): the base mesh is the
    /// rotors at t = 0.
    pub fn angles(&self, t: f64) -> Vec<f64> {
        self.omega.iter().map(|w| w * t).collect()
    }

    /// The exact wall distance of the turbulence model (D48) for every set
    /// of rotor angles, kept per part of the base mesh (`mesh`, `geo`;
    /// [`TurningWallDistance`]); `None` in laminar flow, with the Poisson
    /// distance (computed again on the mesh of every step) or without walls.
    pub fn wall_distance(
        &self,
        case: &FlowCase,
        mesh: &Mesh,
        geo: &Geometry,
    ) -> Option<TurningWallDistance> {
        use crate::turbulence::WallDistanceMethod;
        if case.turbulence.is_laminar()
            || case.turbulence.wall_distance != WallDistanceMethod::Exact
        {
            return None;
        }
        let walls: Vec<usize> = mesh
            .boundaries
            .iter()
            .enumerate()
            .filter(|(_, b)| {
                matches!(
                    case.boundary.get(&b.name),
                    Some(crate::bc::BoundaryCondition::Wall { .. })
                )
            })
            .map(|(i, _)| i)
            .collect();
        if walls.is_empty() {
            return None;
        }
        let parts = self.imp.parts();
        let n = parts.len() as u32;
        let part: Vec<u32> = (0..mesh.n_cells)
            .map(|c| {
                parts
                    .iter()
                    .position(|p| p.rotor_cells()[c])
                    .map_or(n, |i| i as u32)
            })
            .collect();
        Some(TurningWallDistance::of_mesh_parts(
            mesh,
            geo,
            &part,
            &walls,
            parts.iter().map(|p| (p.origin(), p.axis())).collect(),
        ))
    }

    /// The mesh and the geometry of the rotors at time `t`.
    pub fn working(&self, mesh: &Mesh, base: &Geometry, t: f64) -> Result<Working> {
        let angles = self.angles(t);
        let a = self.imp.assemble(mesh, &angles).map_err(|e| {
            err(
                diag::SLIDING_INTERFACE,
                "sliding",
                format!("[[sliding]]: {e}"),
            )
        })?;
        let geo = self.imp.geometry(&a, base);
        Ok(Working {
            geo,
            origin: a.origin,
            pieces: a.pieces.len(),
            mesh: a.mesh,
            angles,
        })
    }

    /// The largest closure defect |Σ S|/Σ |S| of a cell of `w` (round-off
    /// for an extruded ring mesh).
    pub fn closure(&self, w: &Working) -> f64 {
        let mut worst = 0.0f64;
        for c in 0..w.mesh.n_cells {
            let s = w.geo.closure(&w.mesh, c);
            let tot: f64 = w
                .mesh
                .cell_face_iter(c)
                .map(|(f, _)| norm(w.geo.face_area[f]))
                .sum();
            worst = worst.max(norm(s) / tot);
        }
        worst
    }
}

/// Where the faces of one working mesh come from in the previous one.
#[derive(Debug, Clone)]
pub struct FaceMap {
    /// For every face of the new mesh the face of the previous mesh (a face
    /// of the base mesh), or `None` for an intersection polygon.
    pub face: Vec<Option<usize>>,
    /// The same for the boundary faces, counted among the boundary faces.
    pub boundary: Vec<Option<usize>>,
    /// The new faces without a predecessor.
    pub pieces: Vec<usize>,
}

/// The map from the faces of `prev` to those of `next`.
pub fn face_map(prev: &Working, next: &Working) -> FaceMap {
    let n_base = prev
        .origin
        .iter()
        .chain(&next.origin)
        .filter_map(|o| match o {
            FaceOrigin::Base(b) => Some(*b as usize + 1),
            FaceOrigin::Piece(_) => None,
        })
        .max()
        .unwrap_or(0);
    let mut pos = vec![None; n_base];
    for (f, o) in prev.origin.iter().enumerate() {
        if let FaceOrigin::Base(b) = o {
            pos[*b as usize] = Some(f);
        }
    }
    let (ni_prev, ni_next) = (prev.mesh.n_interior_faces(), next.mesh.n_interior_faces());
    let face: Vec<Option<usize>> = next
        .origin
        .iter()
        .map(|o| match o {
            FaceOrigin::Base(b) => pos[*b as usize],
            FaceOrigin::Piece(_) => None,
        })
        .collect();
    let boundary = face[ni_next..]
        .iter()
        .map(|p| p.map(|f| f - ni_prev))
        .collect();
    let pieces = face
        .iter()
        .enumerate()
        .filter(|(_, p)| p.is_none())
        .map(|(f, _)| f)
        .collect();
    FaceMap {
        face,
        boundary,
        pieces,
    }
}

fn remap(v: &[f64], m: &[Option<usize>]) -> Vec<f64> {
    m.iter().map(|o| o.map_or(0.0, |i| v[i])).collect()
}

fn remap3(v: &Soa3, m: &[Option<usize>]) -> Soa3 {
    Soa3 {
        x: remap(&v.x, m),
        y: remap(&v.y, m),
        z: remap(&v.z, m),
    }
}

/// The state and the old levels of the previous mesh on the next one: the
/// cell arrays are the same, the face arrays follow the faces of the base
/// mesh, and the intersection polygons start from the flux of the linearly
/// interpolated velocity (its Rhie–Chow differences 0). The turbulence
/// model's cell fields and old levels go along as they are (the cells keep
/// their numbers); its face values follow the faces, the face viscosity of
/// an intersection is μ + ρ ν_t of the linearly interpolated ν_t (`mu` is
/// the molecular viscosity), and the matrices of its multigrid setups are
/// dropped (the pattern of the matrix changed; the hierarchies are built
/// again).
pub fn carry(
    mut st: State,
    mut lv: OldLevels,
    map: &FaceMap,
    disc: &Discretisation,
    rho: f64,
    mu: f64,
) -> (State, OldLevels) {
    st.mass_flux = remap(&st.mass_flux, &map.face);
    st.ub = remap3(&st.ub, &map.boundary);
    st.pb = remap(&st.pb, &map.boundary);
    if let Some(t) = st.turbulence.as_mut() {
        let b = &map.boundary;
        t.mu_eff = remap(&t.mu_eff, &map.face);
        for &f in &map.pieces {
            let (o, n) = (disc.owner[f] as usize, disc.neighbour[f] as usize);
            let w = disc.weight[f];
            t.mu_eff[f] = mu + rho * (w * t.nu_t[o] + (1.0 - w) * t.nu_t[n]);
        }
        t.k_b = remap(&t.k_b, b);
        t.omega_b = remap(&t.omega_b, b);
        t.u_tau = remap(&t.u_tau, b);
        t.y_plus = remap(&t.y_plus, b);
        t.pb_modified = remap(&t.pb_modified, b);
        if let Some(z) = t.zeta.as_mut() {
            z.zeta_b = remap(&z.zeta_b, b);
            z.f_b = remap(&z.f_b, b);
        }
        t.amg.clear();
    }
    lv.m_old = remap(&lv.m_old, &map.face);
    lv.x_cur = remap(&lv.x_cur, &map.face);
    lv.x_old = remap(&lv.x_old, &map.face);
    lv.ub_old = lv.ub_old.as_ref().map(|u| {
        [
            remap(&u[0], &map.boundary),
            remap(&u[1], &map.boundary),
            remap(&u[2], &map.boundary),
        ]
    });
    for &f in &map.pieces {
        let (o, n) = (disc.owner[f] as usize, disc.neighbour[f] as usize);
        let w = disc.weight[f];
        let sv = disc.face_area.get(f);
        let flux = |u: [&[f64]; 3]| -> f64 {
            (0..3)
                .map(|i| rho * (w * u[i][o] + (1.0 - w) * u[i][n]) * sv[i])
                .sum()
        };
        st.mass_flux[f] = flux([&st.u.x, &st.u.y, &st.u.z]);
        lv.m_old[f] = flux([&lv.u_old[0], &lv.u_old[1], &lv.u_old[2]]);
        lv.x_cur[f] = 0.0;
        lv.x_old[f] = 0.0;
    }
    (st, lv)
}

impl Sliding {
    /// Tells the transient run how the old levels of the time terms move
    /// (call once, before the first step); on ranks (MPI) `on` is this
    /// rank's part and the plan is that of the whole mesh.
    pub(crate) fn install(
        &self,
        fl: &mut Flow,
        tr: &mut Transient,
        on: Option<&crate::ranks::OnRanks>,
    ) -> Result<()> {
        let shift = match on {
            Some(o) => o.shift(self.plan.cell.clone()),
            None => Shift::new(self.plan.cell.clone()),
        };
        tr.set_sliding(fl, shift)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_step_from_sectors_and_speed() {
        let s = SlidingSection {
            omega: Some(10.0),
            shift: 2,
            ..Default::default()
        };
        let dt = s.time_step(64).unwrap();
        assert!((dt - 2.0 * TAU / (64.0 * 10.0)).abs() < 1e-15);
        let rpm = SlidingSection {
            rpm: Some(60.0),
            ..Default::default()
        };
        assert!((rpm.time_step(36).unwrap() - 1.0 / 36.0).abs() < 1e-15);
        assert!(SlidingSection::default().time_step(8).is_err());
    }

    #[test]
    fn rotation_about_an_axis() {
        let r = rotate(
            [0.0, 0.0, 1.0],
            std::f64::consts::FRAC_PI_2,
            [1.0, 0.0, 0.5],
        );
        assert!(r[0].abs() < 1e-15 && (r[1] - 1.0).abs() < 1e-15 && (r[2] - 0.5).abs() < 1e-15);
    }

    #[test]
    fn points_are_found_within_the_tolerance() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0 - 1e-9]];
        let idx = PointIndex::new(&pts, 1e-6);
        assert_eq!(idx.find([1.0, 1e-8, 0.0]), Some(1));
        assert_eq!(idx.find([1.0, 1.0, 1.0]), Some(2));
        assert_eq!(idx.find([0.5, 0.5, 0.5]), None);
    }
}
