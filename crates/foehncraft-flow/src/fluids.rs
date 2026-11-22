//! Several fluids and materials in one case (D194, stage MF-A).
//!
//! The flow domain is the union of the fluid regions. Fluid cells that are
//! connected through faces form a **fluid part**; parts that touch only
//! through solids (water in the tubes, air on the shell side of a heat
//! exchanger) are separate blocks of the flow matrix: their momentum and
//! pressure rows never couple, one AMG setup serves all of them.
//!
//! - Every fluid region may have a material of its own
//!   (`[regions.<name>] material`) and a temperature of its own at which
//!   that material is evaluated (`temperature`); its density, viscosity,
//!   c_p and k are per-cell (and per-face) fields of the solver.
//! - Fluids of different materials (or temperatures) must not share a
//!   face: that would be a two-phase interface, which needs a multiphase
//!   model (L3), so the case is refused with REG004. Regions of the same
//!   fluid may touch: they are one part.
//! - In a case of several parts every part needs its own pressure level:
//!   a pressure boundary (REG003); a closed part with its own reference
//!   cell and thermodynamic pressure comes with stage MF-B.
//! - The log lists the parts, their regions and fluids (REG005).

use crate::bc::BoundaryCondition;
use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::materials::RegionKind;
use anyhow::Result;
use foehncraft_mesh::Mesh;
use serde::Serialize;

/// The fluid parts of a flow mesh: the part of every cell.
#[derive(Debug, Clone)]
pub struct Parts {
    /// The part of each cell (0..n).
    pub part: Vec<u32>,
    /// Number of parts.
    pub n: usize,
}

/// The connected parts of `mesh` (cells joined by interior and periodic
/// faces), numbered in the order of their first cell.
pub fn parts(mesh: &Mesh) -> Parts {
    let nc = mesh.n_cells;
    let mut up: Vec<u32> = (0..nc as u32).collect();
    fn root(up: &mut [u32], mut c: u32) -> u32 {
        while up[c as usize] != c {
            let g = up[up[c as usize] as usize];
            up[c as usize] = g;
            c = g;
        }
        c
    }
    for f in 0..mesh.n_interior_faces() {
        let (a, b) = (
            root(&mut up, mesh.owner[f]),
            root(&mut up, mesh.neighbour[f]),
        );
        if a != b {
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            up[hi as usize] = lo;
        }
    }
    let mut id = vec![u32::MAX; nc];
    let mut part = vec![0u32; nc];
    let mut n = 0u32;
    for (c, pc) in part.iter_mut().enumerate() {
        let r = root(&mut up, c as u32) as usize;
        if id[r] == u32::MAX {
            id[r] = n;
            n += 1;
        }
        *pc = id[r];
    }
    Parts {
        part,
        n: n as usize,
    }
}

/// The fluid of every cell zone of `mesh` (0: the `[fluid]`; `i > 0`: the
/// i-th fluid region with a material or a temperature of its own) and the
/// label of each fluid.
pub fn zone_fluids(case: &FlowCase, mesh: &Mesh) -> (Vec<u32>, Vec<String>) {
    let mut labels = vec![fluid_label(case, None)];
    let mut of_zone = vec![0u32; mesh.zones.len()];
    for (name, r) in &case.regions {
        if r.kind != RegionKind::Fluid
            || !(crate::energy::own_fluid(case, name) || r.temperature.is_some())
        {
            continue;
        }
        labels.push(fluid_label(case, Some(name)));
        let i = (labels.len() - 1) as u32;
        for z in &r.zones {
            if let Some(k) = mesh.zones.iter().position(|y| y == z) {
                of_zone[k] = i;
            }
        }
    }
    (of_zone, labels)
}

fn fluid_label(case: &FlowCase, region: Option<&str>) -> String {
    let r = region.and_then(|n| case.regions.get(n));
    let mat = match (region, r.and_then(|r| r.material.clone())) {
        (Some(n), Some(m)) if crate::energy::own_fluid(case, n) => m,
        _ => case
            .fluid
            .material
            .clone()
            .unwrap_or_else(|| "fluid".to_string()),
    };
    let t = r
        .and_then(|r| r.temperature)
        .unwrap_or_else(|| case.fluid.conditions().temperature);
    format!("{mat} at {t} K")
}

/// True when the case has more than one fluid (D194): a fluid region with
/// a material or a temperature of its own.
pub fn several(case: &FlowCase) -> bool {
    case.regions.iter().any(|(n, r)| {
        r.kind == RegionKind::Fluid
            && (crate::energy::own_fluid(case, n) || r.temperature.is_some())
    })
}

/// The settings that do not yet work with several fluids (REG006).
pub fn validate(case: &FlowCase) -> Result<()> {
    if !several(case) {
        return Ok(());
    }
    let what = if crate::density::wanted(case) {
        Some("the variable density (per-fluid density curves come with stage MF-B)")
    } else if case.gravity.is_some() {
        Some("gravity and buoyancy (a buoyancy model per fluid comes with stage MF-C)")
    } else if !case.mrf.is_empty() || !case.sliding.is_empty() {
        Some("a rotating zone ([[mrf]] or [[sliding]])")
    } else if case.flow.frame_rotation != [0.0; 3] {
        Some("a rotating frame (flow.frame_rotation)")
    } else if !case.scalars.is_empty() {
        Some("passive scalars")
    } else if case.boundary.values().any(|b| b.is_convective_outlet()) {
        Some("a convective outlet")
    } else if case
        .boundary
        .values()
        .any(|b| matches!(b, BoundaryCondition::Periodic { .. }))
    {
        Some("a periodic boundary")
    } else {
        None
    };
    match what {
        Some(w) => Err(diag::err(
            Diagnostic::new(
                diag::REG_FLUIDS_UNSUPPORTED,
                format!(
                    "several fluids in one case (fluid regions with a material or temperature of \
                     their own) do not yet work with {w}"
                ),
            )
            .key("regions"),
        )),
        None => Ok(()),
    }
}

/// One fluid part, for the log and the report.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct PartInfo {
    /// Number of cells.
    pub cells: usize,
    /// The regions (`fluid`: cells outside the regions).
    pub regions: Vec<String>,
    /// The fluids (material at its temperature).
    pub fluids: Vec<String>,
    /// The boundaries that touch the part.
    pub boundaries: Vec<String>,
    /// A boundary of the part fixes the pressure level.
    pub pressure: bool,
}

/// The fluid parts of the flow mesh `mesh` of `case` (the fluid view of a
/// conjugate case): REG004 when two fluids share a face, REG003 when a
/// part of several has no pressure boundary.
pub fn check(case: &FlowCase, mesh: &Mesh) -> Result<Vec<PartInfo>> {
    let (of_zone, labels) = zone_fluids(case, mesh);
    let fluid_of = |c: usize| of_zone[mesh.cell_zone[c] as usize];
    let region_of_zone: Vec<String> = mesh
        .zones
        .iter()
        .map(|z| {
            case.regions
                .iter()
                .find(|(_, r)| r.kind == RegionKind::Fluid && r.zones.contains(z))
                .map_or_else(|| "fluid".to_string(), |(n, _)| n.clone())
        })
        .collect();
    // one fluid: nothing to check (the parts of a single fluid, such as a
    // rotor behind a sliding interface, are joined by their interfaces)
    if labels.len() <= 1 {
        return Ok(Vec::new());
    }
    for f in 0..mesh.n_interior_faces() {
        let (o, n) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
        let (a, b) = (fluid_of(o), fluid_of(n));
        if a != b {
            let ro = &region_of_zone[mesh.cell_zone[o] as usize];
            let rn = &region_of_zone[mesh.cell_zone[n] as usize];
            return Err(diag::err(
                Diagnostic::new(
                    diag::REG_FLUIDS_TOUCH,
                    format!(
                        "region `{ro}` ({}) and region `{rn}` ({}) share faces (one at {:?}): \
                         two different fluids may meet only through a solid region or a thin \
                         wall; a direct contact is a two-phase interface (multiphase, L3)",
                        labels[a as usize],
                        labels[b as usize],
                        face_centroid(mesh, f)
                    ),
                )
                .key(format!("regions.{}", if ro == "fluid" { rn } else { ro })),
            ));
        }
    }
    let p = parts(mesh);
    let mut out: Vec<PartInfo> = (0..p.n)
        .map(|_| PartInfo {
            cells: 0,
            regions: Vec::new(),
            fluids: Vec::new(),
            boundaries: Vec::new(),
            pressure: false,
        })
        .collect();
    for c in 0..mesh.n_cells {
        let i = &mut out[p.part[c] as usize];
        i.cells += 1;
        let r = &region_of_zone[mesh.cell_zone[c] as usize];
        if !i.regions.contains(r) {
            i.regions.push(r.clone());
        }
        let l = &labels[fluid_of(c) as usize];
        if !i.fluids.contains(l) {
            i.fluids.push(l.clone());
        }
    }
    for b in &mesh.boundaries {
        let pressure = case
            .boundary
            .get(&b.name)
            .is_some_and(BoundaryCondition::is_pressure);
        let mut seen = vec![false; p.n];
        for f in b.range() {
            seen[p.part[mesh.owner[f] as usize] as usize] = true;
        }
        for (k, s) in seen.iter().enumerate() {
            if *s {
                out[k].boundaries.push(b.name.clone());
                out[k].pressure |= pressure;
            }
        }
    }
    for i in &mut out {
        i.regions.sort();
        i.fluids.sort();
    }
    if p.n > 1
        && let Some((k, i)) = out.iter().enumerate().find(|(_, i)| !i.pressure)
    {
        return Err(diag::err(
            Diagnostic::new(
                diag::REG_PART_REFERENCE,
                format!(
                    "fluid part {} of {} (regions {}, {} cells, boundaries {}) has no pressure \
                     boundary: in a case of several fluid parts each part needs one (a pressure \
                     outlet or inlet) to fix its pressure level",
                    k + 1,
                    p.n,
                    i.regions.join(", "),
                    i.cells,
                    if i.boundaries.is_empty() {
                        "none".to_string()
                    } else {
                        i.boundaries.join(", ")
                    }
                ),
            )
            .key("boundary"),
        ));
    }
    Ok(out)
}

fn face_centroid(mesh: &Mesh, f: usize) -> [f64; 3] {
    let nodes = mesh.face(f);
    let mut x = [0.0; 3];
    for &v in nodes {
        for (d, xd) in x.iter_mut().enumerate() {
            *xd += mesh.points[v as usize][d];
        }
    }
    x.map(|v| v / nodes.len().max(1) as f64)
}

/// The REG005 note: one line per part, when the case has more than one
/// fluid or part.
pub fn notes(parts: &[PartInfo]) -> Vec<String> {
    let fluids: std::collections::BTreeSet<&String> =
        parts.iter().flat_map(|p| p.fluids.iter()).collect();
    if parts.len() < 2 && fluids.len() < 2 {
        return Vec::new();
    }
    let mut v = vec![format!(
        "[REG005] {} fluid part(s), {} fluid(s)",
        parts.len(),
        fluids.len()
    )];
    for (k, p) in parts.iter().enumerate() {
        v.push(format!(
            "[REG005]   part {}: {} cells, regions {}, fluid {}, boundaries {}",
            k + 1,
            p.cells,
            p.regions.join(", "),
            p.fluids.join(", "),
            p.boundaries.join(", ")
        ));
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_strip_is_one_part() {
        // one strip of four cells: one part
        let src = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [4.0, 1.0, 1.0]
cells = [4, 1, 1]
"#;
        let case: FlowCase = toml::from_str(src).unwrap();
        let (mesh, _, _) = crate::case::mesh_for(&case, None).unwrap();
        let p = parts(&mesh);
        assert_eq!(p.n, 1);
        assert!(p.part.iter().all(|&x| x == 0));
    }
}
