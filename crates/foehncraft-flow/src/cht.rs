//! Conjugate heat transfer (stage 7, PR 7d1; D112, D157): solid regions in
//! the energy equation.
//!
//! One mesh holds the fluid and the solid regions. The flow is solved on
//! the *fluid view*: the cells of the fluid, with every face between a
//! fluid and a solid cell turned into a no-slip wall boundary named
//! `interface-<solid region>`. The energy equation is solved on the *full
//! mesh*, all regions in one matrix: the mass flux of the fluid view is
//! copied onto the faces of the full mesh (zero in the solids), and the
//! faces between two materials take the series conductance
//!
//!   k_f = (δ_P + δ_N) / (δ_P/k_P + δ_N/k_N),
//!
//! with δ the normal distances of the two cell centroids to the face
//! (Patankar 1980, §4.2), so the heat flux is the same number on both sides
//! of the interface at every iteration (D117 (1)). The temperature-
//! dependent viscosity of the fluid goes back to the fluid view after each
//! energy solve.
//!
//! Case file: a `[regions.<name>]` table with `kind = "solid"`, the cell
//! zones and a material with density, specific heat and conductivity
//! (MAT006); the outer boundaries of a solid take the thermal conditions
//! of `[boundary.<name>.energy]` like any wall.

use crate::case::FlowCase;
use crate::materials::RegionKind;
use anyhow::{Result, anyhow};
use foehncraft_mesh::builder::MeshBuilder;
use foehncraft_mesh::{Geometry, Mesh};
use std::collections::BTreeMap;
use std::sync::Arc;

/// The prefix of the boundaries of the fluid view between the fluid and a
/// solid region.
pub const INTERFACE: &str = "interface-";

/// True for a case with the energy equation and at least one solid region.
pub fn is_cht(case: &FlowCase) -> bool {
    case.energy.is_some() && case.regions.values().any(|r| r.kind == RegionKind::Solid)
}

/// The full mesh behind a fluid view, and the maps between the two.
#[derive(Debug)]
pub struct FluidView {
    /// The case as written (all regions, all boundaries).
    pub full_case: FlowCase,
    /// The full mesh and its geometry.
    pub mesh: Mesh,
    /// See `mesh`.
    pub geo: Geometry,
    /// Per cell of the fluid view: the cell of the full mesh.
    pub cell: Vec<u32>,
    /// Per face of the fluid view: the face of the full mesh.
    pub face: Vec<u32>,
    /// Per face of the fluid view: true where its normal is the opposite
    /// of that of the full face (interface faces whose fluid cell is the
    /// neighbour in the full mesh).
    pub flip: Vec<bool>,
    /// The solid region of each cell of the full mesh (`None`: fluid).
    pub solid: Vec<Option<String>>,
    /// The interface boundaries of the fluid view, by solid region.
    pub interfaces: Vec<String>,
}

/// The name of the interface boundary of a solid region.
pub fn interface_name(region: &str) -> String {
    format!("{INTERFACE}{region}")
}

/// The solid region of every cell (by its cell zone).
fn solid_cells(case: &FlowCase, mesh: &Mesh) -> Vec<Option<String>> {
    let mut by_zone: BTreeMap<&str, &str> = BTreeMap::new();
    for (name, r) in &case.regions {
        if r.kind == RegionKind::Solid {
            for z in &r.zones {
                by_zone.insert(z.as_str(), name.as_str());
            }
        }
    }
    let zone_region: Vec<Option<String>> = mesh
        .zones
        .iter()
        .map(|z| by_zone.get(z.as_str()).map(|r| r.to_string()))
        .collect();
    (0..mesh.n_cells)
        .map(|c| zone_region[mesh.cell_zone[c] as usize].clone())
        .collect()
}

/// The cells of region `name` of the energy equation: `fluid` (every
/// fluid cell), a solid region (`mat_cell` indexes `regions`, the fluid
/// first) or a fluid region of `[regions]` (its zones); `None` for a name
/// the case does not have.
pub(crate) fn region_cells(
    case: &FlowCase,
    mesh: &Mesh,
    mat_cell: &[u16],
    regions: &[String],
    name: &str,
) -> Option<Vec<u32>> {
    let m = member(case, mesh, mat_cell, regions, name)?;
    Some(
        (0..mat_cell.len() as u32)
            .filter(|&c| m(c as usize))
            .collect(),
    )
}

/// Whether a cell is in region `name` (see [`region_cells`]), and whether
/// the region is a fluid.
#[allow(clippy::type_complexity)]
fn member<'a>(
    case: &'a FlowCase,
    mesh: &'a Mesh,
    mat_cell: &'a [u16],
    regions: &'a [String],
    name: &str,
) -> Option<Box<dyn Fn(usize) -> bool + 'a>> {
    if name == "fluid" {
        return Some(Box::new(move |c| mat_cell[c] == 0));
    }
    if let Some(i) = regions.iter().skip(1).position(|r| r == name) {
        let i = (i + 1) as u16;
        return Some(Box::new(move |c| mat_cell[c] == i));
    }
    let r = case.regions.get(name)?;
    if r.kind != RegionKind::Fluid {
        return None;
    }
    let zones: Vec<bool> = mesh.zones.iter().map(|z| r.zones.contains(z)).collect();
    Some(Box::new(move |c| {
        mat_cell[c] == 0 && zones[mesh.cell_zone[c] as usize]
    }))
}

/// The `[interfaces]` entries of a conjugate case (stage 7d2): their
/// resistances and, per interior face of the full mesh on one of them, the
/// entry; CHT004 for an entry that is not on a region interface.
#[allow(clippy::type_complexity)]
pub(crate) fn interface_specs(
    case: &FlowCase,
    env: &crate::expr::Env,
    mesh: &Mesh,
    mat_cell: &[u16],
    regions: &[String],
) -> Result<(
    Vec<crate::conduction::Resistance>,
    std::collections::HashMap<u32, usize>,
)> {
    use crate::conduction::Resistance;
    use crate::diag::{self, Diagnostic};
    let bad = |name: &str, msg: String| {
        diag::err(Diagnostic::new(diag::CHT_INTERFACE, msg).key(format!("interfaces.{name}")))
    };
    let ni = mesh.n_interior_faces();
    let mut res = Vec::new();
    let mut of_face = std::collections::HashMap::new();
    let is_fluid = |n: &str| {
        n == "fluid"
            || case
                .regions
                .get(n)
                .is_some_and(|r| r.kind == RegionKind::Fluid)
    };
    for (name, sec) in &case.interfaces {
        let [a, b] = &sec.regions;
        let known: Vec<&str> = std::iter::once("fluid")
            .chain(case.regions.keys().map(|k| k.as_str()))
            .collect();
        let ma = member(case, mesh, mat_cell, regions, a);
        let mb = member(case, mesh, mat_cell, regions, b);
        let (Some(ma), Some(mb)) = (ma, mb) else {
            return Err(bad(
                name,
                format!(
                    "[interfaces.{name}] regions = [\"{a}\", \"{b}\"]: a region the case does not \
                     have; the regions are {}",
                    known.join(", ")
                ),
            ));
        };
        if a == b || (is_fluid(a) && is_fluid(b)) {
            return Err(bad(
                name,
                format!(
                    "[interfaces.{name}] regions = [\"{a}\", \"{b}\"]: an interface is between two \
                     different regions, not two fluids (a wall between two fluids is a thin wall \
                     of D189)"
                ),
            ));
        }
        let r = match (&sec.contact_resistance, &sec.thin_wall) {
            (Some(r), None) => {
                if !(r.is_finite() && *r >= 0.0) {
                    return Err(bad(
                        name,
                        format!("[interfaces.{name}] contact_resistance must be ≥ 0 m²K/W, is {r}"),
                    ));
                }
                Resistance::Fixed(*r)
            }
            (None, Some(w)) => {
                if !(w.thickness.is_finite() && w.thickness > 0.0) {
                    return Err(bad(
                        name,
                        format!(
                            "[interfaces.{name}] thin_wall thickness must be > 0 m, is {}",
                            w.thickness
                        ),
                    ));
                }
                match (&w.conductivity, &w.material) {
                    (Some(k), None) if k.is_finite() && *k > 0.0 => {
                        Resistance::Fixed(w.thickness / k)
                    }
                    (None, Some(m)) => {
                        let sec_m = crate::materials::lookup(m, &case.materials)?;
                        let owner = format!("materials.{m}");
                        let Some(crate::materials::Conductivity::Isotropic(p)) =
                            &sec_m.conductivity
                        else {
                            return Err(bad(
                                name,
                                format!(
                                    "[interfaces.{name}] thin_wall material `{m}` needs one \
                                     conductivity (an isotropic k)"
                                ),
                            ));
                        };
                        let k = crate::materials::curve(
                            p,
                            "conductivity",
                            &owner,
                            &sec_m,
                            env,
                            case.fluid.conditions(),
                        )?;
                        Resistance::Wall(w.thickness, k)
                    }
                    _ => {
                        return Err(bad(
                            name,
                            format!(
                                "[interfaces.{name}] thin_wall needs one of `conductivity` (> 0 \
                                 W/(m K)) and `material`"
                            ),
                        ));
                    }
                }
            }
            _ => {
                return Err(bad(
                    name,
                    format!(
                        "[interfaces.{name}] needs one of `contact_resistance` (m²K/W) and \
                         `thin_wall`"
                    ),
                ));
            }
        };
        let idx = res.len();
        res.push(r);
        let mut count = 0usize;
        for f in 0..ni {
            let (o, n) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
            if (ma(o) && mb(n)) || (mb(o) && ma(n)) {
                if let Some(other) = of_face.insert(f as u32, idx) {
                    let other = case
                        .interfaces
                        .keys()
                        .nth(other)
                        .cloned()
                        .unwrap_or_default();
                    return Err(bad(
                        name,
                        format!(
                            "[interfaces.{name}] and [interfaces.{other}] share faces: one entry \
                             per pair of regions"
                        ),
                    ));
                }
                count += 1;
            }
        }
        if count == 0 {
            return Err(bad(
                name,
                format!(
                    "[interfaces.{name}]: regions `{a}` and `{b}` do not touch (no face between \
                     them)"
                ),
            ));
        }
    }
    Ok((res, of_face))
}

/// The fluid view of a case on its full mesh: the case of the flow (the
/// boundaries of the fluid cells and the interface walls, with the view
/// attached), its mesh and geometry. A case without solid regions comes
/// back as it is.
pub fn view(case: &FlowCase, mesh: Mesh, geo: Geometry) -> Result<(FlowCase, Mesh, Geometry)> {
    let full_case = match &case.fluid_view {
        Some(v) => v.full_case.clone(),
        None => case.clone(),
    };
    if !is_cht(&full_case) {
        crate::need!(
            full_case.interfaces.is_empty(),
            CHT_INTERFACE,
            "[interfaces] needs the energy equation and a solid region ([energy] and \
             [regions.<name>] kind = \"solid\")"
        );
        return Ok((full_case, mesh, geo));
    }
    for (name, r) in &full_case.regions {
        for z in &r.zones {
            crate::need!(
                mesh.zones.contains(z),
                REG_ZONE,
                "[regions.{name}] zone `{z}` is not a cell zone of the mesh; the mesh has {}",
                mesh.zones
                    .iter()
                    .map(|z| format!("`{z}`"))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }
    let solid = solid_cells(&full_case, &mesh);
    let nc = mesh.n_cells;
    let mut local = vec![u32::MAX; nc];
    let mut cell = Vec::new();
    for c in 0..nc {
        if solid[c].is_none() {
            local[c] = cell.len() as u32;
            cell.push(c as u32);
        }
    }
    crate::need!(
        !cell.is_empty(),
        ENE_UNSUPPORTED,
        "every cell is in a solid region: the flow needs fluid cells (a solid-only case is \
         not in this version)"
    );
    // nodes of the fluid cells, in the order of the full mesh
    let mut node = vec![u32::MAX; mesh.n_points()];
    let mut used = Vec::new();
    for &c in &cell {
        for (f, _) in mesh.cell_face_iter(c as usize) {
            for &v in mesh.face(f) {
                if node[v as usize] == u32::MAX {
                    node[v as usize] = 0;
                    used.push(v);
                }
            }
        }
    }
    used.sort_unstable();
    for (i, &v) in used.iter().enumerate() {
        node[v as usize] = i as u32;
    }
    let points: Vec<[f64; 3]> = used.iter().map(|&v| mesh.points[v as usize]).collect();
    let mut b = MeshBuilder::new(points, cell.len());
    b.set_zones(
        mesh.zones.clone(),
        cell.iter().map(|&c| mesh.cell_zone[c as usize]).collect(),
    );
    for p in &mesh.boundaries {
        if p.range().any(|f| solid[mesh.owner[f] as usize].is_none()) {
            b.declare_boundary(&p.name);
        }
    }
    let ni = mesh.n_interior_faces();
    crate::need!(
        mesh.periodic.is_none(),
        ENE_UNSUPPORTED,
        "periodic boundaries with solid regions are not in this version"
    );
    // (insertion order) full face and flip
    let mut inserted: Vec<(u32, bool)> = Vec::new();
    let mut interfaces = Vec::new();
    let remap = |nodes: &[u32], rev: bool| -> Vec<u32> {
        let mut v: Vec<u32> = nodes.iter().map(|&n| node[n as usize]).collect();
        if rev {
            v[1..].reverse();
        }
        v
    };
    for f in 0..mesh.n_faces() {
        let o = mesh.owner[f] as usize;
        if f < ni {
            let n = mesh.neighbour[f] as usize;
            match (&solid[o], &solid[n]) {
                (None, None) => {
                    b.add_interior(&remap(mesh.face(f), false), local[o], local[n]);
                    inserted.push((f as u32, false));
                }
                (None, Some(r)) | (Some(r), None) => {
                    let name = interface_name(r);
                    if mesh.boundary(&name).is_some() {
                        return Err(crate::diag::err(crate::diag::Diagnostic::new(
                            crate::diag::ENE_UNSUPPORTED,
                            format!(
                                "the mesh has a boundary `{name}`: the name of the interface \
                                 of solid region `{r}` is taken; rename the boundary"
                            ),
                        )));
                    }
                    if !interfaces.contains(&name) {
                        interfaces.push(name.clone());
                    }
                    let rev = solid[o].is_some();
                    let fluid_cell = if rev { n } else { o };
                    b.add_boundary(&remap(mesh.face(f), rev), local[fluid_cell], &name);
                    inserted.push((f as u32, rev));
                }
                (Some(_), Some(_)) => {}
            }
        } else if solid[o].is_none() {
            let p = mesh.boundary_of(f).expect("a boundary face");
            b.add_boundary(
                &remap(mesh.face(f), false),
                local[o],
                &mesh.boundaries[p].name,
            );
            inserted.push((f as u32, false));
        }
    }
    let (fluid, order) = b
        .finish_with_order()
        .map_err(|e| anyhow!("fluid view: {e}"))?;
    let face: Vec<u32> = order.iter().map(|&i| inserted[i as usize].0).collect();
    let flip: Vec<bool> = order.iter().map(|&i| inserted[i as usize].1).collect();
    let fgeo = Geometry::compute(&fluid);
    // the case of the flow: the boundaries of the fluid view, the
    // interfaces as walls at rest
    let mut fc = full_case.clone();
    fc.boundary
        .retain(|name, _| fluid.boundaries.iter().any(|p| &p.name == name));
    for name in &interfaces {
        let wall: crate::bc::BoundaryCondition =
            toml::from_str("type = \"wall\"").map_err(|e| anyhow!("interface wall: {e}"))?;
        fc.boundary.insert(name.clone(), wall);
    }
    fc.regions.retain(|_, r| r.kind != RegionKind::Solid);
    fc.fluid_view = Some(Arc::new(FluidView {
        full_case,
        mesh,
        geo,
        cell,
        face,
        flip,
        solid,
        interfaces,
    }));
    Ok((fc, fluid, fgeo))
}

impl FluidView {
    /// Face values of the fluid view on the faces of the full mesh (`sign`:
    /// flip the sign on reversed faces, for fluxes), `fill` elsewhere.
    pub fn faces_to_full(&self, v: &[f64], sign: bool, fill: f64) -> Vec<f64> {
        let mut out = vec![fill; self.mesh.n_faces()];
        for (i, (&f, &r)) in self.face.iter().zip(&self.flip).enumerate() {
            out[f as usize] = if sign && r { -v[i] } else { v[i] };
        }
        out
    }

    /// Cell values of the fluid view on the cells of the full mesh, `fill`
    /// in the solids.
    pub fn cells_to_full(&self, v: &[f64], fill: f64) -> Vec<f64> {
        let mut out = vec![fill; self.mesh.n_cells];
        for (i, &c) in self.cell.iter().enumerate() {
            out[c as usize] = v[i];
        }
        out
    }

    /// Face values of the full mesh on the faces of the fluid view.
    pub fn faces_from_full(&self, v: &[f64]) -> Vec<f64> {
        self.face.iter().map(|&f| v[f as usize]).collect()
    }

    /// Boundary-face values of the fluid view (in its boundary-face order)
    /// on the boundary faces of the full mesh, `fill` elsewhere.
    pub fn boundary_faces_to_full(&self, v: &[f64], fluid_ni: usize, fill: f64) -> Vec<f64> {
        let ni = self.mesh.n_interior_faces();
        let mut out = vec![fill; self.mesh.n_boundary_faces()];
        for (i, &f) in self.face[fluid_ni..].iter().enumerate() {
            let f = f as usize;
            if f >= ni {
                out[f - ni] = v[i];
            }
        }
        out
    }
}

/// The energy equation of a conjugate heat transfer case: a flow object
/// on the full mesh that only solves the energy equation, fed with the
/// mass flux, the velocity gradients and the viscosity of the fluid view.
pub struct Cht {
    /// The energy view (the full mesh).
    pub fe: crate::simplec::Flow,
    /// The maps between the views.
    pub view: Arc<FluidView>,
    /// The turbulent coupling of the energy view (stage 7d5), set up at
    /// the first energy solve of a turbulent flow.
    turb: Option<ChtTurb>,
}

/// The turbulence of the fluid view in the energy view (stage 7d5): the
/// turbulent conductivity c_p μ_t/Pr_t on the fluid faces, the thermal
/// wall function on the walls of the fluid that are boundaries of the full
/// mesh, and the y⁺ of the interface faces, where the wall function goes
/// in series with the solid half cell ([`crate::scalar::InterfaceWall`]).
struct ChtTurb {
    /// Laminar μ per face of the full mesh (1 in the solids).
    mu: crate::device::Buf,
    /// μ_eff per face of the full mesh for the turbulent conductivity: 0 on
    /// the interfaces and in the solids, so that no μ_t is added there.
    mu_eff: crate::device::Buf,
    /// 1 on the fluid's walls, per boundary face of the full mesh.
    wall: crate::device::Buf,
    /// y* per wall face of the full mesh.
    y_plus: crate::device::Buf,
    y_floor: f64,
    /// Per wall face of the full mesh: the wall face of the fluid view.
    outer: Vec<usize>,
    /// Per interface face: the wall face and the face of the fluid view.
    iface: Vec<(usize, usize)>,
    /// The interface faces with their geometry (y* and μ are set per
    /// iteration).
    walls: crate::scalar::InterfaceWall,
}

impl ChtTurb {
    /// The coupling of the walls of the fluid view `fl` (with its
    /// turbulence model) onto the full mesh of `view`.
    fn new(
        fl: &crate::simplec::Flow,
        fe: &mut crate::simplec::Flow,
        view: &FluidView,
    ) -> Result<Option<Self>> {
        let Some(t) = fl.turbulence.as_ref() else {
            return Ok(None);
        };
        let tc = t.scalar_coupling(fl.dev.as_ref());
        let (faces, y) = t.wall_faces();
        let (mesh, geo) = (&view.mesh, &view.geo);
        let ni = mesh.n_interior_faces();
        let mut outer: Vec<(u32, usize)> = Vec::new();
        let mut iface = Vec::new();
        let mut walls = crate::scalar::InterfaceWall {
            y_floor: tc.y_floor,
            ..Default::default()
        };
        for (j, &fv) in faces.iter().enumerate() {
            let fv = fv as usize;
            let full = view.face[fv];
            if (full as usize) >= ni {
                outer.push((full, j));
            } else {
                let f = full as usize;
                let (s, _) = geo.normal_area(f);
                let (o, n) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
                let dn: f64 = (0..3)
                    .map(|k| (geo.cell_centroid[n][k] - geo.cell_centroid[o][k]) * s[k])
                    .sum();
                iface.push((j, fv));
                walls.faces.push(full);
                walls.fluid_owner.push(!view.flip[fv]);
                walls.y.push(y[j]);
                walls.dn.push(dn);
            }
        }
        outer.sort_unstable();
        // the device arrays of the wall faces of the full mesh (face order)
        let mut stat = crate::device::WallStatic::default();
        let mut cells: Vec<(u32, u32)> = Vec::new();
        for (i, &(f, _)) in outer.iter().enumerate() {
            stat.faces.push(f);
            cells.push((mesh.owner[f as usize], i as u32));
        }
        cells.sort_unstable();
        for (k, &(c, i)) in cells.iter().enumerate() {
            if k == 0 || cells[k - 1].0 != c {
                stat.cells.push(c);
                stat.ptr.push(k as u32);
            }
            stat.slots.push(i);
        }
        stat.ptr.push(cells.len() as u32);
        fe.dev.wall_setup(&stat)?;
        let nb = mesh.n_boundary_faces();
        let mut mask = vec![0.0; nb];
        for &(f, _) in &outer {
            mask[f as usize - ni] = 1.0;
        }
        let nf = mesh.n_faces();
        let dev = fe.dev.as_mut();
        let wall = dev.alloc(nb)?;
        dev.upload(wall, &mask)?;
        let y_plus = dev.alloc(outer.len().max(1))?;
        Ok(Some(Self {
            mu: dev.alloc(nf)?,
            mu_eff: dev.alloc(nf)?,
            wall,
            y_plus,
            y_floor: tc.y_floor,
            outer: outer.iter().map(|&(_, j)| j).collect(),
            iface,
            walls,
        }))
    }
}

/// The temperature fields of the energy view of a conjugate case with the
/// device they live on (their time levels, kept and rewound states follow
/// those of the flow, stage 7d6).
pub(crate) fn energy_fields(
    fl: &mut crate::simplec::Flow,
) -> Option<(
    &mut crate::scalar::Scalars,
    &mut dyn crate::device::FlowDevice,
)> {
    let c = fl.cht.as_deref_mut()?;
    let fe = &mut c.fe;
    Some((fe.scalars.as_mut()?, fe.dev.as_mut()))
}

/// The energy equation of a conjugate case for a restart file: the
/// temperature on the full mesh, its level before (transient) and the
/// states of the cycle-averaged loop and of the frozen-flow march, if
/// given; `None` for a case without solids.
pub(crate) fn restart_part(
    fl: &mut crate::simplec::Flow,
    cycle: Option<Vec<u8>>,
    frozen: Option<Vec<u8>>,
) -> anyhow::Result<Option<crate::restart::ChtRestart>> {
    let Some(c) = fl.cht.as_deref() else {
        return Ok(None);
    };
    let cell_id = c.view.mesh.cell_id[..c.view.mesh.n_cells].to_vec();
    let viscosity = fl.phys.mu_face.as_deref().cloned().filter(|_| {
        c.fe.scalars
            .as_ref()
            .is_some_and(|s| s.energy_viscosity_follows())
    });
    let Some((sc, dev)) = energy_fields(fl) else {
        return Ok(None);
    };
    let Some(st) = sc
        .states(dev)?
        .into_iter()
        .find(|s| s.name == crate::energy::FIELD)
    else {
        return Ok(None);
    };
    Ok(Some(crate::restart::ChtRestart {
        cell_id,
        t: st.phi,
        t_old: st.old,
        cycle,
        frozen,
        viscosity,
    }))
}

/// Loads the temperature of a conjugate case from a restart file (mapped
/// by the identities of the full mesh); a conjugate case needs it, a case
/// without solids ignores it.
pub(crate) fn load_restart(
    fl: &mut crate::simplec::Flow,
    r: Option<&crate::restart::ChtRestart>,
    path: &std::path::Path,
) -> anyhow::Result<()> {
    let Some(c) = fl.cht.as_deref() else {
        return Ok(());
    };
    let Some(r) = r else {
        anyhow::bail!(
            "{} holds no temperatures of the solids (written by a run without solid regions or \
             by an older version); a conjugate case needs them",
            path.display()
        );
    };
    let ids = &c.view.mesh.cell_id[..c.view.mesh.n_cells];
    anyhow::ensure!(
        r.cell_id.len() == ids.len(),
        "{}: the restart file has {} cells of the full mesh, the case {}",
        path.display(),
        r.cell_id.len(),
        ids.len()
    );
    let map: Option<Vec<usize>> = if r.cell_id == ids {
        None
    } else {
        let pos: std::collections::HashMap<u64, usize> =
            r.cell_id.iter().enumerate().map(|(i, &g)| (g, i)).collect();
        Some(
            ids.iter()
                .map(|g| pos.get(g).copied())
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| {
                    anyhow!(
                        "{}: the restart file belongs to another mesh",
                        path.display()
                    )
                })?,
        )
    };
    let remap = |v: &[f64]| match &map {
        Some(m) => m.iter().map(|&i| v[i]).collect(),
        None => v.to_vec(),
    };
    let st = crate::scalar::ScalarState {
        name: crate::energy::FIELD.into(),
        phi: remap(&r.t),
        old: r.t_old.as_deref().map(remap),
        wall: None,
    };
    if let Some((sc, dev)) = energy_fields(fl) {
        sc.load(dev, &[st])?;
    }
    // μ(T) as the energy last passed it to the flow (else the first
    // momentum solve would take the viscosity of the case temperature)
    if let Some(mu) = &r.viscosity
        && mu.len() == fl.disc.n_faces
    {
        let target = match fl.turbulence.as_ref() {
            Some(t) => t.scalar_coupling(fl.dev.as_ref()).mu,
            None => fl.dev.face_viscosity(),
        };
        fl.dev.upload(target, mu)?;
        fl.phys.mu_face = Some(Arc::new(mu.clone()));
    }
    Ok(())
}

/// Attaches the energy equation on the full mesh to the flow of the fluid
/// view.
pub fn attach(fl: &mut crate::simplec::Flow, view: Arc<FluidView>) -> Result<()> {
    let full = &view.full_case;
    let (mesh, geo) = (&view.mesh, &view.geo);
    let phys = crate::case::physics(full)?;
    let values = full.compiled_boundaries()?;
    let at = crate::bc::At {
        analytic: None,
        values: &values,
        t: 0.0,
        turning: &full.turning(),
        hydrostatic: crate::buoyancy::hydrostatic(full),
    };
    let arrays = crate::bc::arrays(mesh, geo, &full.boundary, phys.rho, at)?;
    let disc = crate::simplec::discretisation(mesh, geo, &arrays)?;
    let phys = crate::case::physics_on(full, &phys, mesh, geo)?;
    // (the settings are those of the case, checked as transient or steady
    // with it; the energy view takes only its energy and scalar parts)
    let mut fe = crate::simplec::Flow::on(
        disc,
        arrays,
        phys,
        full.solver.clone(),
        full.time.is_some(),
        None,
        &crate::simplec::Target {
            device: full.run.backend.clone(),
            threads: full.run.threads,
            ptx: full.run.ptx.clone(),
        },
    )?;
    let nc = fe.disc.n_cells;
    fe.initialise(&foehncraft_fv::disc::Soa3::zeros(nc), &vec![0.0; nc])?;
    crate::scalar::attach(&mut fe, full, mesh, geo)?;
    fl.cht = Some(Box::new(Cht {
        fe,
        view,
        turb: None,
    }));
    Ok(())
}

/// The energy solve of an outer iteration: the fluxes of the fluid view
/// onto the full mesh, the solve there, μ(T) back to the fluid view.
pub fn solve(
    fl: &mut crate::simplec::Flow,
    step: Option<crate::scalar::StepTime>,
    res: &mut crate::simplec::Residuals,
) -> Result<()> {
    solve_with(fl, step, res, true)
}

/// [`solve`]; `back`: μ(T) goes back to the fluid view (false: a solve
/// that measures and leaves the flow alone, D213).
pub(crate) fn solve_with(
    fl: &mut crate::simplec::Flow,
    step: Option<crate::scalar::StepTime>,
    res: &mut crate::simplec::Residuals,
    back: bool,
) -> Result<()> {
    solve_with_flux(fl, step, res, back, None)
}

/// [`solve_with`]; `flux`: the face mass flux of the fluid view to convect
/// with in place of the current one (the cycle-mean flux of the mean-flow
/// model of the cycle loop, D204).
pub(crate) fn solve_with_flux(
    fl: &mut crate::simplec::Flow,
    step: Option<crate::scalar::StepTime>,
    res: &mut crate::simplec::Residuals,
    back: bool,
    flux: Option<&[f64]>,
) -> Result<()> {
    let mut c = fl.cht.take().expect("a conjugate heat transfer case");
    let r = (|| -> Result<()> {
        let v = c.view.clone();
        let fe = &mut c.fe;
        let m = match flux {
            Some(f) => f.to_vec(),
            None => {
                let conv = fl.convecting_flux()?;
                fl.dev.download(conv)?
            }
        };
        fe.dev
            .upload(fe.f.mass_flux, &v.faces_to_full(&m, true, 0.0))?;
        for i in 0..3 {
            for (src, dst) in fl.f.grad_u[i].all().into_iter().zip(fe.f.grad_u[i].all()) {
                let g = fl.dev.download(src)?;
                fe.dev.upload(dst, &v.cells_to_full(&g, 0.0))?;
            }
        }
        let mu = fl.dev.download(fl.dev.face_viscosity())?;
        let mu_full = v.faces_to_full(&mu, false, 0.0);
        let target = fe.dev.face_viscosity();
        fe.dev.upload(target, &mu_full)?;
        // a turbulent flow (stage 7d5): μ, μ_eff and y* of the fluid view
        // on the full mesh, the interface faces with their y* and μ
        if c.turb.is_none() && fl.turbulence.is_some() {
            c.turb = ChtTurb::new(fl, fe, &v)?;
        }
        let mut coupling = None;
        let mut lam_target = fl.dev.face_viscosity();
        if let (Some(ct), Some(t)) = (c.turb.as_mut(), fl.turbulence.as_ref()) {
            let tc = t.scalar_coupling(fl.dev.as_ref());
            lam_target = tc.mu;
            let lam = fl.dev.download(tc.mu)?;
            let ys = match tc.y_plus {
                Some(b) => fl.dev.download(b)?,
                None => Vec::new(),
            };
            fe.dev.upload(ct.mu, &v.faces_to_full(&lam, false, 1.0))?;
            let mut me = v.faces_to_full(&mu, false, 0.0);
            for &f in &ct.walls.faces {
                me[f as usize] = 0.0;
            }
            fe.dev.upload(ct.mu_eff, &me)?;
            if !ct.outer.is_empty() {
                let yo: Vec<f64> = ct.outer.iter().map(|&j| ys[j]).collect();
                fe.dev.upload(ct.y_plus, &yo)?;
            }
            ct.walls.y_star = ct.iface.iter().map(|&(j, _)| ys[j]).collect();
            ct.walls.mu = ct.iface.iter().map(|&(_, fv)| lam[fv]).collect();
            if let Some(sc) = fe.scalars.as_mut() {
                sc.set_energy_interfaces(ct.walls.clone());
            }
            coupling = Some(crate::scalar::TurbCoupling {
                mu: ct.mu,
                mu_eff: ct.mu_eff,
                wall: ct.wall,
                y_plus: (!ct.outer.is_empty()).then_some(ct.y_plus),
                y_floor: ct.y_floor,
            });
        }
        fe.solve_scalar_fields(step, res, coupling)?;
        // μ(T) back to the fluid view: the face viscosity of a laminar
        // flow, the laminar viscosity of a turbulence model (which adds μ_t)
        if let (true, Some(mu)) = (back, fe.scalars.as_ref().and_then(|s| s.energy_viscosity())) {
            let back = v.faces_from_full(mu);
            fl.dev.upload(lam_target, &back)?;
            fl.phys.mu_face = Some(Arc::new(back));
        }
        Ok(())
    })();
    fl.cht = Some(c);
    r
}

/// The scalar summaries of a conjugate heat transfer case: those of the
/// full mesh, with the interfaces as boundaries (the heat flow into the
/// fluid).
pub fn summaries(
    fl: &mut crate::simplec::Flow,
    dim: u8,
    t: f64,
) -> Result<Vec<crate::scalar::ScalarSummary>> {
    let c = fl.cht.as_mut().expect("a conjugate heat transfer case");
    let v = c.view.clone();
    let mut out = crate::scalar::summaries(&mut c.fe, &v.full_case, &v.mesh, &v.geo, dim, t)?;
    let Some(sc) = c.fe.scalars.as_ref() else {
        return Ok(out);
    };
    let faces = sc.energy_interfaces();
    if let Some(s) = out.iter_mut().find(|s| s.name == crate::energy::FIELD) {
        for (name, q, area, mean) in faces {
            s.boundaries.push(crate::scalar::ScalarBoundarySummary {
                boundary: name,
                mean,
                mixed_mean: None,
                flux: if area > 0.0 { q / area } else { 0.0 },
                heat_flow: Some(q),
                enthalpy_flow: None,
            });
        }
    }
    Ok(out)
}

/// The solution on the full mesh for the output: velocity, pressure and
/// mass flux of the fluid view (zero in the solids), the temperature of
/// the energy view.
pub fn full_state(
    fl: &mut crate::simplec::Flow,
    fluid: &crate::simplec::State,
) -> Result<crate::simplec::State> {
    let c = fl.cht.as_mut().expect("a conjugate heat transfer case");
    let v = &c.view;
    let mut s = c.fe.state()?;
    let fni = fluid.mass_flux.len() - fluid.pb.len();
    for (dst, src) in [
        (&mut s.u.x, &fluid.u.x),
        (&mut s.u.y, &fluid.u.y),
        (&mut s.u.z, &fluid.u.z),
        (&mut s.p, &fluid.p),
    ] {
        *dst = v.cells_to_full(src, 0.0);
    }
    for i in 0..3 {
        for (dst, src) in [
            (&mut s.grad_u[i].x, &fluid.grad_u[i].x),
            (&mut s.grad_u[i].y, &fluid.grad_u[i].y),
            (&mut s.grad_u[i].z, &fluid.grad_u[i].z),
        ] {
            *dst = v.cells_to_full(src, 0.0);
        }
    }
    for (dst, src) in [
        (&mut s.grad_p.x, &fluid.grad_p.x),
        (&mut s.grad_p.y, &fluid.grad_p.y),
        (&mut s.grad_p.z, &fluid.grad_p.z),
    ] {
        *dst = v.cells_to_full(src, 0.0);
    }
    s.mass_flux = v.faces_to_full(&fluid.mass_flux, true, 0.0);
    for (dst, src) in [
        (&mut s.ub.x, &fluid.ub.x),
        (&mut s.ub.y, &fluid.ub.y),
        (&mut s.ub.z, &fluid.ub.z),
        (&mut s.pb, &fluid.pb),
    ] {
        *dst = v.boundary_faces_to_full(src, fni, 0.0);
    }
    Ok(s)
}
