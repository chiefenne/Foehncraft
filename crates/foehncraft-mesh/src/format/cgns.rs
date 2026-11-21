//! Reader for CGNS files in the HDF5 encoding, unstructured 3D zones (D187).
//!
//! Written from the public description of the CGNS standard (SIDS and the
//! HDF5 mapping, <https://cgns.github.io>): every CGNS node is an HDF5
//! group whose attribute `label` names its type (`Zone_t`, `Elements_t`,
//! ...) and whose dataset ` data` (with the leading blank) holds its value.
//!
//! ```text
//! CGNSBase_t            data [cell dimension, physical dimension]
//!   Zone_t              data [vertices, cells, boundary vertices]
//!     ZoneType_t        "Unstructured"
//!     GridCoordinates_t CoordinateX, CoordinateY, CoordinateZ (DataArray_t)
//!     Elements_t        data [element type, boundary]; ElementRange (first,
//!                       last element number), ElementConnectivity (1-based
//!                       vertices); NGON_n and NFACE_n also ElementStartOffset
//!     ZoneBC_t
//!       BC_t            PointRange or PointList: element numbers of the faces
//! ```
//!
//! - One 3D unstructured zone. Cells: tetrahedra (`TETRA_4`), pyramids
//!   (`PYRA_5`), prisms (`PENTA_6`), hexahedra (`HEXA_8`), `MIXED` sections of
//!   them, and polyhedra as `NGON_n` faces with `NFACE_n` cells. Higher-order
//!   elements are refused. Element numbers run across all sections of the
//!   zone, as the standard has it.
//! - Each `BC_t` becomes the boundary of its name. Its faces are the 2D
//!   elements (`TRI_3`, `QUAD_4`, or `NGON_n` faces) listed in its
//!   `PointRange` or `PointList` (also the legacy `ElementRange` and
//!   `ElementList`); entries that are not faces (points and lines that a
//!   mesher also lists) are ignored. Boundary faces that no BC names go to
//!   `default`.
//! - A single zone is the cell zone `fluid`. With several zones (a rotor and a
//!   stator) all are read into one mesh, each `Zone_t` a cell zone of its
//!   name; the nodes of different zones stay apart (an interface is then two
//!   boundaries, as `[[sliding]]` wants it) unless `mesh.merge_nodes` joins
//!   coincident nodes. BCs of equal name in several zones are one boundary.
//!   2D bases and the legacy ADF encoding are refused with a message.

use super::{CellKind, ElementMesh, MeshFile, MeshFormat, ReadOptions};
use hdf5_metno as h5;
use std::path::Path;

/// The CGNS format in the registry of [`crate::format`].
#[derive(Debug, Clone, Copy)]
pub struct Cgns;

impl MeshFormat for Cgns {
    fn name(&self) -> &'static str {
        "cgns"
    }

    fn title(&self) -> &'static str {
        "CGNS (HDF5), unstructured 3D, linear elements and polyhedra"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &["cgns"]
    }

    fn sniff(&self, head: &[u8]) -> bool {
        head.starts_with(b"\x89HDF\r\n\x1a\n")
    }

    fn read_path(&self, path: &Path, options: &ReadOptions) -> Result<MeshFile, String> {
        read_file(path, options)
    }
}

// element types of the standard
const NODE: i32 = 2;
const BAR_2: i32 = 3;
const TRI_3: i32 = 5;
const QUAD_4: i32 = 7;
const TETRA_4: i32 = 10;
const PYRA_5: i32 = 12;
const PENTA_6: i32 = 14;
const HEXA_8: i32 = 17;
const MIXED: i32 = 20;
const NGON_N: i32 = 22;
const NFACE_N: i32 = 23;

fn type_name(t: i32) -> &'static str {
    match t {
        6 => "TRI_6",
        8 => "QUAD_8",
        9 => "QUAD_9",
        11 => "TETRA_10",
        13 => "PYRA_14",
        15 => "PENTA_15",
        16 => "PENTA_18",
        18 => "HEXA_20",
        19 => "HEXA_27",
        21 => "PYRA_13",
        _ => "a higher-order or unsupported element",
    }
}

/// Nodes of a fixed element type.
fn nodes_of(t: i32) -> Option<usize> {
    Some(match t {
        NODE => 1,
        BAR_2 => 2,
        TRI_3 => 3,
        QUAD_4 => 4,
        TETRA_4 => 4,
        PYRA_5 => 5,
        PENTA_6 => 6,
        HEXA_8 => 8,
        _ => return None,
    })
}

fn e<T: std::fmt::Display>(x: T) -> String {
    x.to_string()
}

/// A CGNS node: an HDF5 group.
struct Node(h5::Group);

impl Node {
    fn text_attr(&self, name: &str) -> Option<String> {
        let a = self.0.attr(name).ok()?;
        let s = a.read_scalar::<h5::types::FixedAscii<64>>().ok()?;
        Some(s.as_str().trim_end_matches('\0').to_string())
    }

    fn label(&self) -> String {
        self.text_attr("label").unwrap_or_default()
    }

    fn ints(&self) -> Result<Vec<i64>, String> {
        self.0
            .dataset(" data")
            .and_then(|d| d.read_raw::<i64>())
            .map_err(e)
    }

    fn floats(&self) -> Result<Vec<f64>, String> {
        self.0
            .dataset(" data")
            .and_then(|d| d.read_raw::<f64>())
            .map_err(e)
    }

    /// The value of a `C1` node as text.
    fn text(&self) -> String {
        self.0
            .dataset(" data")
            .and_then(|d| d.read_raw::<i8>())
            .map(|v| {
                v.iter()
                    .take_while(|&&c| c != 0)
                    .map(|&c| c as u8 as char)
                    .collect::<String>()
            })
            .unwrap_or_default()
    }

    fn children(&self) -> Vec<(String, Node)> {
        let mut out = Vec::new();
        for name in self.0.member_names().unwrap_or_default() {
            if name == " data" {
                continue;
            }
            if let Ok(g) = self.0.group(&name) {
                out.push((name, Node(g)));
            }
        }
        out
    }

    fn child(&self, name: &str) -> Option<Node> {
        self.0.group(name).ok().map(Node)
    }

    fn with_label(&self, label: &str) -> Vec<(String, Node)> {
        self.children()
            .into_iter()
            .filter(|(_, n)| n.label() == label)
            .collect()
    }
}

/// An `Elements_t` section, with every element's type and nodes found.
struct Section {
    first: i64,
    last: i64,
    conn: Vec<i64>,
    /// (type, start in `conn`, length) of each element.
    elems: Vec<(i32, usize, usize)>,
}

fn section(name: &str, n: &Node) -> Result<Section, String> {
    let d = n.ints()?;
    let ty = *d
        .first()
        .ok_or_else(|| format!("section '{name}' has no element type"))? as i32;
    let range = n
        .child("ElementRange")
        .ok_or_else(|| format!("section '{name}' has no ElementRange"))?
        .ints()?;
    if range.len() < 2 {
        return Err(format!("section '{name}': ElementRange needs two numbers"));
    }
    let (first, last) = (range[0], range[1]);
    let count = last.saturating_sub(first).saturating_add(1).max(0) as usize;
    let conn = match n.child("ElementConnectivity") {
        Some(c) => c.ints()?,
        None => Vec::new(),
    };
    // every element has one entry in the connectivity at least (D202)
    let count = super::declared(
        count,
        conn.len(),
        &format!("section '{name}': elements of the ElementRange"),
    )?;
    let mut elems = Vec::with_capacity(count);
    match ty {
        MIXED => {
            let mut at = 0usize;
            for _ in 0..count {
                let t = *conn
                    .get(at)
                    .ok_or_else(|| format!("section '{name}': connectivity ends early"))?
                    as i32;
                let Some(nn) = nodes_of(t) else {
                    return Err(format!(
                        "section '{name}': element type {t} ({}) is not read; only linear \
                         elements and polyhedra are (mesh with first-order elements)",
                        type_name(t)
                    ));
                };
                if at + 1 + nn > conn.len() {
                    return Err(format!("section '{name}': connectivity ends early"));
                }
                elems.push((t, at + 1, nn));
                at += 1 + nn;
            }
        }
        NGON_N | NFACE_N => {
            let off = n
                .child("ElementStartOffset")
                .ok_or_else(|| {
                    format!("section '{name}': NGON_n and NFACE_n need ElementStartOffset")
                })?
                .ints()?;
            if off.len() != count + 1 {
                return Err(format!(
                    "section '{name}': {} offsets for {count} elements",
                    off.len()
                ));
            }
            for k in 0..count {
                let (a, b) = (off[k] as usize, off[k + 1] as usize);
                if b < a || b > conn.len() {
                    return Err(format!("section '{name}': bad offsets"));
                }
                elems.push((ty, a, b - a));
            }
        }
        t => {
            let Some(nn) = nodes_of(t) else {
                return Err(format!(
                    "section '{name}': element type {t} ({}) is not read; only linear elements \
                     and polyhedra are (mesh with first-order elements)",
                    type_name(t)
                ));
            };
            if conn.len() != count * nn {
                return Err(format!(
                    "section '{name}': {} connectivity entries for {count} elements of {nn} nodes",
                    conn.len()
                ));
            }
            for k in 0..count {
                elems.push((t, k * nn, nn));
            }
        }
    }
    Ok(Section {
        first,
        last,
        conn,
        elems,
    })
}

struct Sections(Vec<Section>);

impl Sections {
    /// The element with number `id`: (type, entries of `conn`).
    fn element(&self, id: i64) -> Option<(i32, &[i64])> {
        let k = self.0.partition_point(|s| s.first <= id).checked_sub(1)?;
        let s = &self.0[k];
        if id > s.last {
            return None;
        }
        let (t, a, n) = *s.elems.get(id.checked_sub(s.first)? as usize)?;
        Some((t, &s.conn[a..a + n]))
    }
}

fn node_idx(v: i64, n_points: usize, off: u32) -> Result<u32, String> {
    if v < 1 || v as usize > n_points {
        return Err(format!("a vertex number {v} lies outside 1..={n_points}"));
    }
    Ok((v - 1) as u32 + off)
}

/// Reads a CGNS file (HDF5 encoding).
pub fn read_file(path: &Path, options: &ReadOptions) -> Result<MeshFile, String> {
    let mut sig = [0u8; 8];
    {
        use std::io::Read;
        let mut f = std::fs::File::open(path).map_err(e)?;
        let got = f.read(&mut sig).map_err(e)?;
        if got < 8 || &sig != b"\x89HDF\r\n\x1a\n" {
            return Err(
                "not an HDF5 file: only the HDF5 encoding of CGNS is read (convert an \
                        ADF file with `adf2hdf`)"
                    .into(),
            );
        }
    }
    let file = h5::File::open(path).map_err(e)?;
    let root = Node(file.as_group().map_err(e)?);
    let bases = root.with_label("CGNSBase_t");
    let Some((base_name, base)) = bases.first() else {
        return Err("the file has no CGNSBase_t".into());
    };
    if bases.len() > 1 {
        return Err(format!("{} bases: only one base is read", bases.len()));
    }
    let dims = base.ints()?;
    if dims.first() != Some(&3) {
        return Err(format!(
            "base '{base_name}' has cell dimension {}: only 3D meshes are read",
            dims.first().copied().unwrap_or(0)
        ));
    }
    let zones = base.with_label("Zone_t");
    if zones.is_empty() {
        return Err(format!("base '{base_name}' has no Zone_t"));
    }
    let single = zones.len() == 1;
    let mut em = ElementMesh {
        face_word: "boundary faces of the BCs",
        group_word: "BC",
        ..Default::default()
    };
    let mut warnings = Vec::new();
    for (zone_name, zone) in &zones {
        let cell_zone = if single { "fluid" } else { zone_name.as_str() };
        read_zone(&mut em, &mut warnings, base, zone_name, zone, cell_zone)?;
    }
    if !single {
        warnings.push(format!(
            "{} zones read as cell zones{}",
            zones.len(),
            if options.merge_nodes {
                ""
            } else {
                "; their nodes are kept apart (mesh.merge_nodes joins coincident nodes)"
            }
        ));
    }
    if options.merge_nodes {
        let removed = em.merge_points();
        warnings.push(format!("{removed} coincident nodes merged"));
    }
    let built = em.into_mesh()?;
    warnings.extend(built.warnings);
    Ok(MeshFile {
        mesh: built.mesh,
        reader: "cgns",
        format: "CGNS HDF5".into(),
        warnings,
    })
}

/// Adds one zone to `em`: its nodes (after those already there), cells in
/// cell zone `cell_zone` and the faces of its BCs.
fn read_zone(
    em: &mut ElementMesh,
    warnings: &mut Vec<String>,
    base: &Node,
    zone_name: &str,
    zone: &Node,
    cell_zone: &str,
) -> Result<(), String> {
    if let Some(t) = zone.child("ZoneType").map(|n| n.text())
        && t != "Unstructured"
    {
        return Err(format!(
            "zone '{zone_name}' is {t}: only Unstructured is read"
        ));
    }
    let zd = zone.ints()?;
    let n_points = *zd.first().ok_or("the zone has no sizes")? as usize;

    // coordinates
    let gc = zone
        .with_label("GridCoordinates_t")
        .into_iter()
        .next()
        .ok_or("the zone has no GridCoordinates")?
        .1;
    let coord = |name: &str| -> Result<Vec<f64>, String> {
        let v = gc
            .child(name)
            .ok_or_else(|| format!("GridCoordinates has no {name}"))?
            .floats()?;
        if v.len() != n_points {
            return Err(format!(
                "{name} has {} values for {n_points} vertices",
                v.len()
            ));
        }
        Ok(v)
    };
    let (x, y, z) = (
        coord("CoordinateX")?,
        coord("CoordinateY")?,
        coord("CoordinateZ")?,
    );
    let off = em.points.len() as u32;
    em.points.extend((0..n_points).map(|i| [x[i], y[i], z[i]]));

    // sections, by first element number
    let mut secs = Vec::new();
    for (name, n) in zone.with_label("Elements_t") {
        secs.push(section(&name, &n)?);
    }
    secs.sort_by_key(|s| s.first);
    let secs = Sections(secs);

    let z0 = em.zone(cell_zone);
    let mut n_cells = 0usize;
    for s in &secs.0 {
        for &(t, a, n) in &s.elems {
            let nodes = &s.conn[a..a + n];
            let kind = match t {
                TETRA_4 => CellKind::Tetrahedron,
                PYRA_5 => CellKind::Pyramid,
                PENTA_6 => CellKind::Prism,
                HEXA_8 => CellKind::Hexahedron,
                NFACE_N => {
                    let mut faces = Vec::with_capacity(n);
                    for &f in nodes {
                        // the sign says which side the face normal looks at; the
                        // orientation of the cell follows from its volume, so
                        // faces of one cell are made to agree with the sign
                        let (_, nodes) = secs.element(f.abs()).ok_or_else(|| {
                            format!("a polyhedron refers to the missing face {}", f.abs())
                        })?;
                        let mut v = nodes
                            .iter()
                            .map(|&x| node_idx(x, n_points, off))
                            .collect::<Result<Vec<_>, _>>()?;
                        if f < 0 {
                            v.reverse();
                        }
                        faces.push(v);
                    }
                    em.add_polyhedron(faces, z0);
                    n_cells += 1;
                    continue;
                }
                _ => continue,
            };
            let v = nodes
                .iter()
                .map(|&x| node_idx(x, n_points, off))
                .collect::<Result<Vec<_>, _>>()?;
            em.add_cell(kind, &v, z0);
            n_cells += 1;
        }
    }
    let n_elements: usize = secs.0.iter().map(|s| s.elems.len()).sum();
    if n_cells == 0 {
        return Err("zone has no 3D cells (TETRA_4, PYRA_5, PENTA_6, HEXA_8, NFACE_n)".into());
    }

    // boundary conditions
    let mut bcs: Vec<(i64, String, Vec<i64>)> = Vec::new();
    for zbc in zone.with_label("ZoneBC_t") {
        for (name, bc) in zbc.1.with_label("BC_t") {
            if let Some(g) = bc.child("GridLocation").map(|n| n.text())
                && g == "Vertex"
            {
                return Err(format!(
                    "BC '{name}' is given on vertices; the faces of a boundary must be given as \
                     elements (GridLocation FaceCenter)"
                ));
            }
            let ids = if let Some(r) = bc.child("PointRange").or_else(|| bc.child("ElementRange")) {
                let r = r.ints()?;
                if r.len() < 2 {
                    return Err(format!("BC '{name}': a range needs two numbers"));
                }
                let len = r[1].saturating_sub(r[0]).saturating_add(1).max(0) as usize;
                super::declared(len, n_elements, &format!("BC '{name}': the range"))?;
                (r[0]..=r[1]).collect()
            } else if let Some(l) = bc.child("PointList").or_else(|| bc.child("ElementList")) {
                l.ints()?
            } else {
                warnings.push(format!("BC '{name}' lists no elements and is ignored"));
                continue;
            };
            // Gmsh names the BC after its family and keeps the physical-group name
            // as a FamilyName_t child of that family; other writers use the BC name
            let name = bc
                .child("FamilyName")
                .map(|n| n.text())
                .and_then(|f| {
                    let (_, fam) = base.with_label("Family_t").into_iter().find(|x| x.0 == f)?;
                    fam.children().into_iter().find_map(|(n, c)| {
                        (c.label() == "FamilyName_t").then(|| {
                            let t = c.text();
                            if t.is_empty() { n } else { t }
                        })
                    })
                })
                .unwrap_or(name);
            bcs.push((ids.iter().copied().min().unwrap_or(0), name, ids));
        }
    }
    bcs.sort();
    for (_, name, ids) in &bcs {
        let b = em.boundary(name);
        for &id in ids {
            let Some((t, nodes)) = secs.element(id) else {
                return Err(format!("BC '{name}' lists the missing element {id}"));
            };
            if !matches!(t, TRI_3 | QUAD_4 | NGON_N) {
                continue; // a point, a line or a cell: not a boundary face
            }
            let v = nodes
                .iter()
                .map(|&x| node_idx(x, n_points, off))
                .collect::<Result<Vec<_>, _>>()?;
            em.add_face(&v, b);
        }
    }
    if bcs.is_empty() {
        warnings.push(format!(
            "zone '{zone_name}' has no ZoneBC: its boundary faces are in 'default'"
        ));
    }
    Ok(())
}
