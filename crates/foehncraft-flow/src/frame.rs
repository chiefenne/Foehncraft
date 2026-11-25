//! The moving reference frame (MRF, D55): zones of the mesh that rotate
//! about an axis, for pumps, fans and mixers ("frozen rotor"), in steady
//! state and in a transient run (D185). Case file: one `[[mrf]]` table per
//! zone.
//!
//! # Equations
//!
//! The solved velocity is the **absolute** velocity **u** (in the inertial
//! basis); the convecting velocity of a rotating zone is the **relative**
//! one, **u**ᵣ = **u** − **Ω** × (**x** − **x**₀):
//!
//! ∇·(ρ **u**ᵣ ⊗ **u**) + ρ **Ω** × **u** = −∇p + ∇·τ(**u**)
//!
//! (a flow that is steady in the rotating frame has ∂**u**/∂t = **Ω** ×
//! **u** − (**Ω** × **r**)·∇**u** in the inertial frame; the first term is
//! the Coriolis-type source, the second moves the convecting velocity). The
//! stress is that of the absolute velocity: the symmetric part of ∇**u** is
//! the same in both frames, so a solid-body rotation is stress free. A
//! turbulence quantity or a passive scalar φ (a scalar of the rotating
//! frame) is convected with the relative flux:
//! ∇·(ρ **u**ᵣ φ) = …, and its production uses the strain rate, which is
//! the same in both frames (the vorticity of the Kato–Launder form is the
//! absolute one, TURB004).
//!
//! # Discretisation
//!
//! - The continuity equation, the Rhie–Chow flux and the pressure
//!   correction work on the **absolute** mass flux φ (ṁ_f), unchanged.
//! - The frame flux of a face is ρ **S**·**V**, with **V** = **Ω** × (**x**_f
//!   − **x**₀) in a rotating zone and 0 elsewhere; on a face between two
//!   zones the mean of the two values (it vanishes when the interface is a
//!   surface of revolution about the axis, **S** ⟂ **V**); on a boundary
//!   face the value of its cell. The relative flux ṁ_f − ρ **S**·**V** is
//!   the convecting flux of the momentum rows, of the turbulence variables
//!   and of the scalars. A wall that rotates with the frame (`velocity`
//!   expression of D53) has the absolute flux ρ **S**·**u**_wall = ρ **S**·**V**
//!   and so the relative flux 0.
//! - Σ_f ρ **S**·**V**(**x**_f) = 0 for every closed cell with planar faces
//!   (**V** is linear in **x** and divergence free), so the relative flux
//!   is divergence free like the absolute one wherever all faces of a cell
//!   have their exact frame flux: inside a zone, on boundaries and at an
//!   interface that is a surface of revolution. Next to any other
//!   interface it is not (the known limit of the method); MRF003 reports
//!   the worst cell.
//! - The source −ρ V_c **Ω** × **u**_c of the rotating cells is explicit
//!   (lagged one outer iteration).
//!
//! # Transient MRF (D185)
//!
//! In a transient run the equation gets the time term of the velocity
//! whose components are those of the rotating basis, at a point that
//! turns with the frame, u′(x′, t) = Rᵀ(t) **u**(R(t) x′, t): its time
//! derivative is Rᵀ (∂**u**/∂t + (**Ω** × **r**)·∇**u** − **Ω** × **u**), so
//! at the frozen instant (R = I)
//!
//! ∂(ρ **u**′)/∂t + ∇·(ρ **u**ᵣ ⊗ **u**) + ρ **Ω** × **u** = −∇p + ∇·τ(**u**),
//!
//! the steady equation with the ordinary time term of the solved velocity:
//! the transient frozen rotor. It is exact in a zone that is a body of
//! revolution about its axis with walls that turn with it or are bodies of
//! revolution too (a Taylor–Couette spin-up in the rotating frame is the
//! spin-up of the inertial frame); with blades it keeps the rotor at one
//! position relative to the stator, so the rotor–stator interaction is that
//! of the frozen position (the sliding mesh of D182, [`crate::sliding`],
//! moves the rotor instead: same frame flux as its grid flux, no Coriolis
//! source). A transient MRF zone and sliding rotors can share a case, in
//! different zones (`MRF001`).

use crate::case::FlowCase;
use crate::device::Buf;
use crate::device::Buf3;
use crate::diag::{self, Diagnostic};
use crate::simplec::Flow;
use anyhow::Result;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// `[[mrf]]`: a zone that rotates steadily about an axis (D55).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MrfSection {
    /// The cell zone of the mesh that rotates (a zone name of the mesh;
    /// without it all cells rotate, which needs a single `[[mrf]]` entry).
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
    /// The walls of the rotor, whose torque gives the shaft power (D184);
    /// default: every wall all of whose faces lie in the zone.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub walls: Vec<String>,
}

impl Default for MrfSection {
    fn default() -> Self {
        Self {
            zone: None,
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            omega: None,
            rpm: None,
            walls: Vec::new(),
        }
    }
}

impl MrfSection {
    /// The angular velocity vector Ω (rad/s).
    pub fn omega_vector(&self) -> Result<[f64; 3]> {
        let w = match (self.omega, self.rpm) {
            (Some(w), None) => w,
            (None, Some(r)) => r * std::f64::consts::TAU / 60.0,
            _ => anyhow::bail!("give exactly one of `omega` and `rpm`"),
        };
        anyhow::ensure!(w.is_finite(), "the angular velocity must be finite");
        let n = self.axis.iter().map(|a| a * a).sum::<f64>().sqrt();
        anyhow::ensure!(
            n > 0.0 && n.is_finite() && self.axis.iter().all(|a| a.is_finite()),
            "`axis` must be a non-zero vector"
        );
        Ok(self.axis.map(|a| w * a / n))
    }
}

/// The checks of the `[[mrf]]` tables that need no mesh (D55, D182).
pub fn validate(case: &FlowCase) -> Result<()> {
    if case.mrf.is_empty() {
        return Ok(());
    }
    // a zone is turned by one method: a frame or a sliding mesh
    for (i, m) in case.mrf.iter().enumerate() {
        for s in &case.sliding {
            crate::need!(
                m.zone.is_some() && s.zone.is_some() && m.zone != s.zone,
                MRF_SLIDING_ZONE,
                "[[mrf]] entry {} and a [[sliding]] table turn the same cells (zone {}): a zone \
                 is turned by a moving reference frame or by a sliding mesh, not both; give \
                 each its own `zone`",
                i + 1,
                m.zone
                    .as_deref()
                    .map_or("all cells".to_string(), |z| format!("`{z}`"))
            );
        }
    }
    for (i, m) in case.mrf.iter().enumerate() {
        if let Err(e) = m.omega_vector() {
            let mut d = Diagnostic::new(diag::MRF_ENTRY, format!("[[mrf]] entry {}: {e}", i + 1));
            d.key = Some(format!("mrf[{i}]"));
            return Err(diag::err(d));
        }
        crate::need!(
            m.zone.is_some() || case.mrf.len() == 1,
            MRF_ENTRY,
            "[[mrf]] entry {}: several entries need a `zone` each",
            i + 1
        );
        for (j, o) in case.mrf.iter().enumerate().skip(i + 1) {
            crate::need!(
                m.zone.is_none() || o.zone.is_none() || m.zone != o.zone,
                MRF_ENTRY,
                "[[mrf]] entries {} and {} name the zone `{}` twice",
                i + 1,
                j + 1,
                m.zone.as_deref().unwrap_or("")
            );
        }
    }
    Ok(())
}

/// A rotating zone on the device.
#[derive(Debug, Clone)]
pub(crate) struct FrameZone {
    /// Angular velocity vector (rad/s).
    pub omega: [f64; 3],
    /// ρ V_c in the cells of the zone, 0 elsewhere (a cell field).
    pub weight: Buf,
    /// The Coriolis source acts in the zone (the moving reference frame);
    /// not for the cells of a sliding mesh, whose mesh moves (D182).
    pub coriolis: bool,
}

/// The device data of the frames of a run.
#[derive(Debug, Clone)]
pub(crate) struct Frames {
    /// ρ **S**·**V** per face (kg/s).
    pub flux: Buf,
    /// The convecting flux ṁ − flux (a work field, rewritten on use).
    pub conv: Buf,
    /// One entry per `[[mrf]]` table.
    pub zones: Vec<FrameZone>,
    /// Work field for the Coriolis source.
    pub work: Buf3,
}

/// What [`attach`] reports for the log and the diagnostics.
#[derive(Debug, Default)]
pub struct FrameNotes {
    /// Log lines.
    pub lines: Vec<String>,
    /// Notes and warnings.
    pub diagnostics: Vec<Diagnostic>,
}

/// The frame flux and the cell imbalance, on the host (the testable part of
/// [`attach`]).
#[derive(Debug, Clone, PartialEq)]
pub struct FrameFluxes {
    /// ρ **S**·**V** per face (kg/s).
    pub flux: Vec<f64>,
    /// Zone entry (index into `[[mrf]]`) of each cell, if any.
    pub cell_zone: Vec<Option<usize>>,
    /// Largest |Σ_f ± flux_f| / Σ_f ρ |V_f| |S_f| over the cells (0 for
    /// exact closure).
    pub worst_imbalance: f64,
    /// Cells with an imbalance above 10⁻⁶.
    pub unbalanced_cells: usize,
}

/// Relative closure above which MRF003 warns.
const IMBALANCE_WARN: f64 = 1e-6;

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Assigns each cell to the `[[mrf]]` entry of its zone.
pub fn assign_cells(entries: &[MrfSection], mesh: &Mesh) -> Result<Vec<Option<usize>>> {
    assign_cells_in(entries, mesh, "[[mrf]]")
}

/// [`assign_cells`] for the entries of another table (`[[sliding]]`), named
/// `table` in the messages.
pub fn assign_cells_in(
    entries: &[MrfSection],
    mesh: &Mesh,
    table: &str,
) -> Result<Vec<Option<usize>>> {
    let nc = mesh.n_cells;
    let mut out: Vec<Option<usize>> = vec![None; nc];
    for (i, m) in entries.iter().enumerate() {
        let zone = match &m.zone {
            Some(name) => {
                let Some(z) = mesh.zones.iter().position(|n| n == name) else {
                    let hint = crate::bc::did_you_mean(name, mesh.zones.iter().map(String::as_str))
                        .map(|c| format!(" (did you mean `{c}`?)"))
                        .unwrap_or_default();
                    let mut d = Diagnostic::new(
                        diag::MRF_ENTRY,
                        format!(
                            "{table} entry {}: the mesh has no cell zone `{name}`{hint}; its zones are {:?}",
                            i + 1,
                            mesh.zones
                        ),
                    );
                    d.key = Some(format!("{}[{i}].zone", table.trim_matches(['[', ']'])));
                    return Err(diag::err(d));
                };
                Some(z as u32)
            }
            None => None,
        };
        for (c, slot) in out.iter_mut().enumerate() {
            if zone.is_none_or(|z| mesh.cell_zone[c] == z) {
                anyhow::ensure!(
                    slot.is_none(),
                    "{table} entries {} and {} both rotate cell {c}",
                    slot.unwrap_or(0) + 1,
                    i + 1
                );
                *slot = Some(i);
            }
        }
    }
    Ok(out)
}

/// The frame flux of every face and the closure of the cells (see the
/// module documentation).
pub fn frame_fluxes(
    entries: &[MrfSection],
    cell_zone: &[Option<usize>],
    disc: &foehncraft_fv::Discretisation,
    rho: f64,
) -> Result<FrameFluxes> {
    let omegas: Vec<[f64; 3]> = entries
        .iter()
        .map(MrfSection::omega_vector)
        .collect::<Result<_>>()?;
    let velocity = |c: usize, x: [f64; 3]| -> [f64; 3] {
        match cell_zone[c] {
            Some(z) => {
                let o = entries[z].origin;
                cross(omegas[z], [x[0] - o[0], x[1] - o[1], x[2] - o[2]])
            }
            None => [0.0; 3],
        }
    };
    let nf = disc.n_faces;
    let mut flux = vec![0.0; nf];
    // (on ranks the halo cells take the sums of the faces they share with
    // own cells; only the own cells are judged)
    let mut sum = vec![0.0; cell_zone.len().max(disc.n_cells)];
    let mut abs = vec![0.0; cell_zone.len().max(disc.n_cells)];
    for (f, out) in flux.iter_mut().enumerate() {
        let x = [
            disc.face_centroid.x[f],
            disc.face_centroid.y[f],
            disc.face_centroid.z[f],
        ];
        let s = [
            disc.face_area.x[f],
            disc.face_area.y[f],
            disc.face_area.z[f],
        ];
        let p = disc.owner[f] as usize;
        let mut v = velocity(p, x);
        if f < disc.n_interior {
            let n = disc.neighbour[f] as usize;
            let vn = velocity(n, x);
            for i in 0..3 {
                v[i] = 0.5 * (v[i] + vn[i]);
            }
        }
        let m = rho * (v[0] * s[0] + v[1] * s[1] + v[2] * s[2]);
        // the scale of the closure error: ρ |V| |S|, not |ρ V·S| (which
        // is rounding noise on a face that the frame velocity slides along)
        let scale = rho * norm(v) * norm(s);
        *out = m;
        sum[p] += m;
        abs[p] += scale;
        if f < disc.n_interior {
            let n = disc.neighbour[f] as usize;
            sum[n] -= m;
            abs[n] += scale;
        }
    }
    let mut worst = 0.0_f64;
    let mut count = 0;
    for c in 0..disc.n_cells {
        if abs[c] > 0.0 {
            let r = sum[c].abs() / abs[c];
            worst = worst.max(r);
            if r > IMBALANCE_WARN {
                count += 1;
            }
        }
    }
    Ok(FrameFluxes {
        flux,
        cell_zone: cell_zone.to_vec(),
        worst_imbalance: worst,
        unbalanced_cells: count,
    })
}

/// The rotating zones of `case` as frames: the `[[mrf]]` tables (with the
/// Coriolis source) and the `[[sliding]]` tables (without it, D182).
fn entries(case: &FlowCase) -> (Vec<MrfSection>, Vec<bool>) {
    let mut e = case.mrf.clone();
    let mut coriolis = vec![true; e.len()];
    for s in &case.sliding {
        e.push(s.frame());
        coriolis.push(false);
    }
    (e, coriolis)
}

/// Adds the frames of `case` to the flow solver `fl` (no-op without
/// `[[mrf]]` or `[[sliding]]`); call before loading a restart.
pub fn attach(
    fl: &mut Flow,
    case: &FlowCase,
    mesh: &Mesh,
    _geo: &Geometry,
) -> Result<Option<FrameNotes>> {
    let (entries, coriolis) = entries(case);
    crate::need!(
        !(case.turbulence.curvature_correction && !case.sliding.is_empty()),
        CASE_INCONSISTENT,
        "turbulence.curvature_correction does not work with [[sliding]] yet (the time part of DS/Dt \
         would need the strain rate that each rotor cell carries); use [[mrf]] for the rotor"
    );
    let spin = case.flow.frame_rotation;
    if spin != [0.0; 3] {
        return attach_rotating_domain(fl, case, spin).map(Some);
    }
    if entries.is_empty() {
        return Ok(None);
    }
    validate(case)?;
    let mut cell_zone = assign_cells_in(
        &entries,
        mesh,
        if case.sliding.is_empty() {
            "[[mrf]]"
        } else {
            "[[sliding]]"
        },
    )?;
    let comm = fl.dev.comm();
    if comm.is_some() {
        // on ranks (MPI): the zones of the halo cells from their owners, so
        // that a face between an own and a halo cell gets the frame flux of
        // the one-rank run
        cell_zone = halo_zones(fl, &cell_zone)?;
    }
    let rho = fl.phys.rho;
    let mut ff = frame_fluxes(&entries, &cell_zone, &fl.disc, rho)?;
    if let Some(c) = &comm {
        // the closure of the cells of all ranks
        ff.worst_imbalance = c.max(ff.worst_imbalance);
        let mut n = [ff.unbalanced_cells as i64];
        c.all_sum_i64(&mut n);
        ff.unbalanced_cells = n[0] as usize;
    }
    let nf = fl.disc.n_faces;
    let nc = fl.disc.n_cells;
    let flux = fl.dev.alloc(nf)?;
    fl.dev.upload(flux, &ff.flux)?;
    let conv = fl.dev.alloc(nf)?;
    let mut zones = Vec::new();
    let mut notes = FrameNotes::default();
    for (i, m) in entries.iter().enumerate() {
        let omega = m.omega_vector()?;
        let w: Vec<f64> = (0..nc)
            .map(|c| {
                if cell_zone[c] == Some(i) {
                    rho * fl.disc.cell_volume[c]
                } else {
                    0.0
                }
            })
            .collect();
        let mut cells = [w.iter().filter(|&&x| x > 0.0).count() as i64];
        if let Some(c) = &comm {
            c.all_sum_i64(&mut cells);
        }
        let cells = cells[0];
        let weight = fl.dev.alloc(nc)?;
        fl.dev.upload(weight, &w)?;
        zones.push(FrameZone {
            omega,
            weight,
            coriolis: coriolis[i],
        });
        if !coriolis[i] {
            // the sliding mesh reports itself (crate::sliding)
            continue;
        }
        let mag = omega.iter().map(|a| a * a).sum::<f64>().sqrt();
        let sign = if m.omega.is_some_and(|w| w < 0.0) || m.rpm.is_some_and(|r| r < 0.0) {
            -1.0
        } else {
            1.0
        };
        notes.lines.push(format!(
            "moving reference frame (D55): {} ({cells} cells) turns about the axis {:?} through {:?} at {:.6e} rad/s ({:.4} rpm)",
            m.zone
                .as_deref()
                .map_or("all cells".to_string(), |z| format!("zone `{z}`")),
            m.axis,
            m.origin,
            sign * mag,
            sign * mag * 60.0 / std::f64::consts::TAU
        ));
    }
    // the curvature correction of SST (D170): Ω^rot and the frame velocity
    // of the cells of the moving reference frames
    if case.turbulence.curvature_correction {
        let omegas: Vec<[f64; 3]> = entries
            .iter()
            .map(MrfSection::omega_vector)
            .collect::<Result<_>>()?;
        let c = &fl.disc.cell_centroid;
        let mut om = vec![[0.0; 3]; nc];
        let mut vel = vec![[0.0; 3]; nc];
        for i in 0..nc {
            if let Some(z) = cell_zone[i]
                && coriolis[z]
            {
                let o = entries[z].origin;
                om[i] = omegas[z];
                vel[i] = cross(omegas[z], [c.x[i] - o[0], c.y[i] - o[1], c.z[i] - o[2]]);
            }
        }
        if let Some(t) = fl.turbulence.as_mut() {
            t.set_curvature_frames(&mut *fl.dev, &om, &vel)?;
        }
    }
    let work = Buf3 {
        x: fl.dev.alloc(nc)?,
        y: fl.dev.alloc(nc)?,
        z: fl.dev.alloc(nc)?,
    };
    fl.frames = Some(Frames {
        flux,
        conv,
        zones,
        work,
    });
    if ff.unbalanced_cells > 0 {
        let what = if case.sliding.is_empty() {
            "the solution there is not that of a frozen rotor"
        } else {
            "the grid flux of the moving mesh does not close there (the geometric conservation law is violated)"
        };
        let m = format!(
            "the relative flux of {} cells does not close (worst: {:.2e} of ρ|V||S| summed over its faces): \
             the interface of a rotating zone is not a surface of revolution about its axis, or a boundary \
             face of the zone is not; {what}",
            ff.unbalanced_cells, ff.worst_imbalance
        );
        notes.lines.push(m.clone());
        let mut d = Diagnostic::new(diag::MRF_INTERFACE, m);
        d.key = Some("mrf".into());
        notes.diagnostics.push(d);
    }
    if case.turbulence.production_limiter == Some(crate::turbulence::ProductionLimiter::KatoLaunder)
    {
        let m = "the Kato–Launder production uses the vorticity of the absolute velocity in a rotating zone (D55); \
                 the strain rate is the same in both frames"
            .to_string();
        notes.lines.push(m.clone());
        let mut d = Diagnostic::new(diag::TURB_KATO_ROTATING, m);
        d.key = Some("turbulence.production_limiter".into());
        notes.diagnostics.push(d);
    }
    Ok(Some(notes))
}

/// The zone entries of the own cells `own` (local order) and of the halo
/// cells of a rank (MPI), these from their owners through the halo
/// exchange of the device.
fn halo_zones(fl: &mut Flow, own: &[Option<usize>]) -> Result<Vec<Option<usize>>> {
    let code: Vec<f64> = own.iter().map(|z| z.map_or(-1.0, |i| i as f64)).collect();
    let b = fl.dev.alloc(own.len())?;
    fl.dev.upload(b, &code)?;
    let all = fl.dev.download_with_halo(b)?;
    Ok(all
        .iter()
        .map(|&v| (v >= 0.0).then_some(v as usize))
        .collect())
}

/// The whole domain in a frame that turns at `spin` (`flow.frame_rotation`):
/// no frame flux (the solved velocity is the relative one), one zone over
/// all cells whose Coriolis source −ρ V (2 **Ω**) × **u** the momentum rows
/// take.
fn attach_rotating_domain(fl: &mut Flow, case: &FlowCase, spin: [f64; 3]) -> Result<FrameNotes> {
    crate::need!(
        case.mrf.is_empty() && case.sliding.is_empty(),
        CASE_INCONSISTENT,
        "flow.frame_rotation turns the whole domain; it does not combine with [[mrf]] or [[sliding]]"
    );
    crate::need!(
        spin.iter().all(|w| w.is_finite()),
        CASE_RANGE,
        "flow.frame_rotation must be finite"
    );
    let (nc, nf) = (fl.disc.n_cells, fl.disc.n_faces);
    let rho = fl.phys.rho;
    let flux = fl.dev.alloc(nf)?;
    fl.dev.fill(flux, 0.0)?;
    let conv = fl.dev.alloc(nf)?;
    let w: Vec<f64> = (0..nc).map(|c| rho * fl.disc.cell_volume[c]).collect();
    let weight = fl.dev.alloc(nc)?;
    fl.dev.upload(weight, &w)?;
    let work = Buf3 {
        x: fl.dev.alloc(nc)?,
        y: fl.dev.alloc(nc)?,
        z: fl.dev.alloc(nc)?,
    };
    fl.frames = Some(Frames {
        flux,
        conv,
        zones: vec![FrameZone {
            omega: spin.map(|x| 2.0 * x),
            weight,
            coriolis: true,
        }],
        work,
    });
    let mag = spin.iter().map(|a| a * a).sum::<f64>().sqrt();
    let mut notes = FrameNotes::default();
    notes.lines.push(format!(
        "rotating frame: the whole domain turns at {spin:?} rad/s (|Ω| = {mag:.6e} rad/s); the velocity \
         is the relative one, the Coriolis force −2ρ Ω × u acts, the pressure is p − ½ρ|Ω × r|²"
    ));
    Ok(notes)
}

impl Flow {
    /// The flux that convects momentum, turbulence variables and scalars:
    /// the mass flux, or without frames the relative flux ṁ − ρ **S**·**V**.
    pub(crate) fn convecting_flux(&mut self) -> Result<Buf> {
        self.relative_flux(self.f.mass_flux)
    }

    /// The relative flux of the face flux `m` (the work field of the frames
    /// holds it), or `m` itself without frames.
    pub(crate) fn relative_flux(&mut self, m: Buf) -> Result<Buf> {
        let Some(fr) = self.frames.as_ref() else {
            // (the grid flux of a deforming mesh, stage MS0)
            return self.grid_relative(m);
        };
        let (conv, flux) = (fr.conv, fr.flux);
        self.dev.lincomb(conv, 1.0, m, -1.0, flux)?;
        Ok(conv)
    }

    /// The frame flux ρ **S**·**V** per face, if there are frames.
    pub(crate) fn frame_flux(&self) -> Option<Buf> {
        match self.frames.as_ref() {
            Some(fr) => Some(fr.flux),
            // the grid mass flux of a deforming mesh (stage MS0)
            None => self.grid.as_deref().and_then(|g| g.flux).map(|(f, _)| f),
        }
    }

    /// Adds the Coriolis source −ρ V **Ω** × **u** of the rotating cells to
    /// the right-hand sides of the momentum rows (explicit, from the current
    /// velocity).
    pub(crate) fn add_coriolis(&mut self) -> Result<()> {
        let Some(fr) = self.frames.as_ref() else {
            return Ok(());
        };
        let work = fr.work;
        let zones: Vec<(Buf, [f64; 3])> = fr
            .zones
            .iter()
            .filter(|z| z.coriolis)
            .map(|z| (z.weight, z.omega))
            .collect();
        let (u, rhs) = (self.f.u, self.w.rhs);
        for (weight, om) in zones {
            self.dev.scale_vector(weight, u, work)?;
            let (t, r) = (work.all(), rhs.all());
            // (Ω × t)_i = Ω_j t_k − Ω_k t_j for (i, j, k) cyclic
            for (i, &ri) in r.iter().enumerate() {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                if om[j] != 0.0 {
                    self.dev.lincomb(ri, 1.0, ri, -om[j], t[k])?;
                }
                if om[k] != 0.0 {
                    self.dev.lincomb(ri, 1.0, ri, om[k], t[j])?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn omega_from_rpm_and_axis() {
        let m = MrfSection {
            rpm: Some(60.0),
            axis: [0.0, 0.0, 2.0],
            ..Default::default()
        };
        let o = m.omega_vector().unwrap();
        assert!((o[2] - std::f64::consts::TAU).abs() < 1e-12 && o[0] == 0.0 && o[1] == 0.0);
        let both = MrfSection {
            rpm: Some(1.0),
            omega: Some(1.0),
            ..Default::default()
        };
        assert!(both.omega_vector().is_err());
        let none = MrfSection::default();
        assert!(none.omega_vector().is_err());
        let zero_axis = MrfSection {
            omega: Some(1.0),
            axis: [0.0; 3],
            ..Default::default()
        };
        assert!(zero_axis.omega_vector().is_err());
    }
}
