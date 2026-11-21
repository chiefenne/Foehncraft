//! Reader for SU2 native mesh files (`.su2`, ASCII), 2D and 3D (D187).
//!
//! Written from the public description of the format in the SU2 user
//! documentation ("Mesh File"):
//!
//! ```text
//! % comment
//! NDIME= 3                 dimension, 2 or 3
//! NELEM= 2                 cells, one per line: VTK type, node indices
//! 10 0 1 2 3 0             (from 0), optionally the element index
//! ...
//! NPOIN= 5                 nodes: x y (z), optionally the node index;
//! 0.0 0.0 0.0 0            a second count (older parallel files) is ignored
//! ...
//! NMARK= 2                 boundary markers
//! MARKER_TAG= wall
//! MARKER_ELEMS= 4
//! 5 0 1 2                  boundary elements: VTK type, node indices
//! ...
//! ```
//!
//! - Cells: tetrahedron (VTK type 10), hexahedron (12), prism (13),
//!   pyramid (14); in 2D triangles (5) and quadrilaterals (9).
//! - Boundary elements: triangles and quadrilaterals in 3D, lines (3) in 2D.
//!   Each marker becomes the boundary of its name (`MARKER_TAG`); the case
//!   maps it to its boundary condition, renamed with `[mesh.boundaries]` if
//!   needed. Boundary faces of no marker go to `default`, with a warning.
//! - A 2D file is extruded one cell layer in z (D3), as the built-in 2D
//!   meshes: triangles become prisms, quadrilaterals hexahedra, the markers
//!   name the side faces and the two planes are the boundaries `zmin`
//!   (z = 0) and `zmax`. The layer thickness is `mesh.thickness` (in the
//!   units of the file), by default the median edge length of the 2D cells.
//! - One zone (`NZONE= 1` or none): the cells are the cell zone `fluid`.
//!   With `NZONE= n` and `IZONE=` blocks (3D only) every zone is a cell zone
//!   `zone1` .. `zone<n>` of one mesh (D195): the nodes of the zones stay apart,
//!   markers of equal name in several zones are one boundary, and
//!   `mesh.merge_nodes` joins coincident nodes. The interface connectivity
//!   entries (`MARKER_TAG= SEND_RECEIVE`) are not read.
//! - Higher-order elements are refused.

use super::{CellKind, ElementMesh, MeshFile, MeshFormat, ReadOptions, declared};

/// The SU2 format in the registry of [`crate::format`](mod@crate::format).
#[derive(Debug, Clone, Copy)]
pub struct Su2;

impl MeshFormat for Su2 {
    fn name(&self) -> &'static str {
        "su2"
    }

    fn title(&self) -> &'static str {
        "SU2 native mesh, ASCII, 2D (extruded one layer) and 3D"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &["su2"]
    }

    fn sniff(&self, head: &[u8]) -> bool {
        let text = String::from_utf8_lossy(head);
        text.lines()
            .map(|l| l.split('%').next().unwrap_or("").trim())
            .find(|l| !l.is_empty())
            .is_some_and(|l| {
                let k = l.split('=').next().unwrap_or("").trim();
                l.contains('=') && matches!(k, "NDIME" | "NZONE")
            })
    }

    fn parse(&self, bytes: &[u8], options: &ReadOptions) -> Result<MeshFile, String> {
        parse(bytes, options)
    }
}

/// VTK cell types of the SU2 format.
const LINE: u32 = 3;
const TRIANGLE: u32 = 5;
const QUADRILATERAL: u32 = 9;
const TETRAHEDRON: u32 = 10;
const HEXAHEDRON: u32 = 12;
const PRISM: u32 = 13;
const PYRAMID: u32 = 14;

fn nodes_of(t: u32) -> Result<usize, String> {
    Ok(match t {
        1 => 1,
        LINE => 2,
        TRIANGLE => 3,
        QUADRILATERAL => 4,
        TETRAHEDRON => 4,
        HEXAHEDRON => 8,
        PRISM => 6,
        PYRAMID => 5,
        21..=37 | 68..=72 => {
            return Err(format!(
                "element type {t} is a higher-order element; only linear elements are read"
            ));
        }
        _ => return Err(format!("unsupported element type {t}")),
    })
}

/// An element as read: VTK type and nodes.
struct Element {
    kind: u32,
    nodes: Vec<u32>,
}

struct Lines<'a> {
    it: std::iter::Enumerate<std::str::Lines<'a>>,
    line: usize,
    /// Size of the file in bytes (a count may not exceed it, D202).
    len: usize,
}

impl<'a> Lines<'a> {
    /// The next line with content, comments (`%`) removed.
    fn next(&mut self) -> Option<&'a str> {
        for (i, l) in self.it.by_ref() {
            let l = l.split('%').next().unwrap_or("").trim();
            if !l.is_empty() {
                self.line = i + 1;
                return Some(l);
            }
        }
        None
    }

    fn err(&self, msg: impl std::fmt::Display) -> String {
        format!("line {}: {msg}", self.line)
    }

    /// The next line, which must be `key= value`; returns the value.
    fn keyword(&mut self, key: &str) -> Result<&'a str, String> {
        let l = self
            .next()
            .ok_or_else(|| format!("unexpected end of file, expected {key}="))?;
        match split_keyword(l) {
            Some((k, v)) if k == key => Ok(v),
            _ => Err(self.err(format!("expected {key}=, found '{l}'"))),
        }
    }

    fn element(&mut self, what: &str) -> Result<Element, String> {
        let l = self
            .next()
            .ok_or_else(|| format!("unexpected end of file in the {what}"))?;
        let mut t = l.split_whitespace();
        let kind: u32 = parse_num(t.next().unwrap_or(""))
            .ok_or_else(|| self.err(format!("expected an element type, found '{l}'")))?;
        let n = nodes_of(kind).map_err(|e| self.err(e))?;
        let mut nodes = Vec::with_capacity(n);
        for _ in 0..n {
            let v = t.next().and_then(parse_num::<u32>).ok_or_else(|| {
                self.err(format!("element of type {kind} needs {n} nodes: '{l}'"))
            })?;
            nodes.push(v);
        }
        Ok(Element { kind, nodes })
    }
}

fn split_keyword(l: &str) -> Option<(&str, &str)> {
    let (k, v) = l.split_once('=')?;
    let k = k.trim();
    (!k.is_empty()
        && k.bytes()
            .all(|c| c.is_ascii_uppercase() || c == b'_' || c.is_ascii_digit()))
    .then_some((k, v.trim()))
}

fn parse_num<T: std::str::FromStr>(s: &str) -> Option<T> {
    s.parse().ok()
}

fn count(lines: &Lines, key: &str, v: &str) -> Result<usize, String> {
    let n = v
        .split_whitespace()
        .next()
        .and_then(parse_num)
        .ok_or_else(|| lines.err(format!("{key}= needs a count, found '{v}'")))?;
    declared(n, lines.len, &lines.err(format!("{key}=")))
}

/// Parses the bytes of a `.su2` file.
pub fn parse(bytes: &[u8], options: &ReadOptions) -> Result<MeshFile, String> {
    let text = std::str::from_utf8(bytes).map_err(|_| "an SU2 file is text (UTF-8 or ASCII)")?;
    let mut lines = Lines {
        it: text.lines().enumerate(),
        line: 0,
        len: text.len(),
    };
    let mut ndime: Option<usize> = None;
    let mut cells: Option<Vec<Element>> = None;
    let mut points: Option<Vec<[f64; 3]>> = None;
    let mut markers: Option<Vec<(String, Vec<Element>)>> = None;
    let mut warnings = Vec::new();
    let mut nzone = 1usize;
    let mut done: Vec<ZoneData> = Vec::new();

    while let Some(l) = lines.next() {
        let Some((key, value)) = split_keyword(l) else {
            return Err(lines.err(format!(
                "expected a keyword (NDIME=, NELEM=, ...), found '{l}'"
            )));
        };
        match key {
            "NDIME" => {
                let d = count(&lines, key, value)?;
                if d != 2 && d != 3 {
                    return Err(lines.err(format!("NDIME= {d}: the dimension is 2 or 3")));
                }
                ndime = Some(d);
            }
            "NZONE" => {
                nzone = count(&lines, key, value)?;
            }
            "IZONE" => {
                // the next zone starts: close the one before
                if cells.is_some() || points.is_some() {
                    done.push(finish_zone(
                        ndime,
                        cells.take(),
                        points.take(),
                        markers.take(),
                        &mut warnings,
                    )?);
                }
            }
            "NELEM" => {
                let n = count(&lines, key, value)?;
                let mut v = Vec::with_capacity(n);
                for _ in 0..n {
                    v.push(lines.element("cells (NELEM)")?);
                }
                cells = Some(v);
            }
            "NPOIN" => {
                let d = ndime.ok_or_else(|| lines.err("NPOIN= before NDIME="))?;
                let n = count(&lines, key, value)?;
                let mut v = Vec::with_capacity(n);
                for _ in 0..n {
                    let l = lines
                        .next()
                        .ok_or("unexpected end of file in the nodes (NPOIN)")?;
                    let mut t = l.split_whitespace();
                    let mut p = [0.0; 3];
                    for x in p.iter_mut().take(d) {
                        *x = t.next().and_then(parse_num).ok_or_else(|| {
                            lines.err(format!("a node needs {d} coordinates: '{l}'"))
                        })?;
                    }
                    v.push(p);
                }
                points = Some(v);
            }
            "NMARK" => {
                let n = count(&lines, key, value)?;
                let mut v = Vec::with_capacity(n);
                for _ in 0..n {
                    let name = lines.keyword("MARKER_TAG")?.to_string();
                    if name.is_empty() {
                        return Err(lines.err("MARKER_TAG= without a name"));
                    }
                    let k = lines.keyword("MARKER_ELEMS")?;
                    let k = count(&lines, "MARKER_ELEMS", k)?;
                    let mut e = Vec::with_capacity(k);
                    for _ in 0..k {
                        e.push(lines.element(&format!("marker '{name}'"))?);
                    }
                    v.push((name, e));
                }
                markers = Some(v);
            }
            _ => {
                if ndime.is_some() && cells.is_some() && points.is_some() && markers.is_some() {
                    warnings.push(format!(
                        "the file goes on after the mesh with {key}=; the rest is ignored"
                    ));
                    break;
                }
                return Err(lines.err(format!("unknown keyword {key}= inside the mesh")));
            }
        }
    }
    done.push(finish_zone(ndime, cells, points, markers, &mut warnings)?);
    let ndime = ndime.ok_or("no NDIME= (is this an SU2 mesh file?)")?;
    if nzone != done.len() {
        warnings.push(format!(
            "NZONE= {nzone} but the file has {} zones; all of them are read",
            done.len()
        ));
    }
    if done.len() > 1 {
        return multi_zone(done, ndime, options, warnings);
    }
    let ZoneData {
        cells,
        points,
        markers,
    } = done.pop().expect("one zone");
    let np = points.len() as u32;

    let mut em = ElementMesh {
        face_word: "marker elements",
        group_word: "marker",
        ..Default::default()
    };
    let z = em.zone("fluid");
    for (name, _) in &markers {
        em.boundary(name);
    }
    if ndime == 3 {
        em.points = points;
        for e in &cells {
            let kind = match e.kind {
                TETRAHEDRON => CellKind::Tetrahedron,
                HEXAHEDRON => CellKind::Hexahedron,
                PRISM => CellKind::Prism,
                PYRAMID => CellKind::Pyramid,
                t => return Err(format!("element type {t} is not a 3D cell (NDIME= 3)")),
            };
            em.add_cell(kind, &e.nodes, z);
        }
        for (b, (name, es)) in markers.iter().enumerate() {
            for e in es {
                if e.kind != TRIANGLE && e.kind != QUADRILATERAL {
                    return Err(format!(
                        "marker '{name}': element type {} is not a face of a 3D mesh \
                         (triangle 5, quadrilateral 9)",
                        e.kind
                    ));
                }
                em.add_face(&e.nodes, b as u32);
            }
        }
    } else {
        // one layer in z: node i at z = 0, node i + n at z = thickness
        let t = match options.thickness {
            Some(t) if t > 0.0 && t.is_finite() => t,
            Some(t) => return Err(format!("mesh.thickness must be positive, not {t}")),
            None => median_edge(&cells, &points)
                .ok_or("cannot find a layer thickness: the 2D cells have no edges of length > 0")?,
        };
        let n = np;
        let mut pts = points.clone();
        pts.extend(points.iter().map(|p| [p[0], p[1], t]));
        em.points = pts;
        let zmin = em.boundary("zmin");
        let zmax = em.boundary("zmax");
        let mut v = Vec::with_capacity(8);
        for e in &cells {
            let kind = match e.kind {
                TRIANGLE => CellKind::Prism,
                QUADRILATERAL => CellKind::Hexahedron,
                t => {
                    return Err(format!(
                        "element type {t} is not a 2D cell (NDIME= 2: triangle 5, quadrilateral 9)"
                    ));
                }
            };
            v.clear();
            v.extend_from_slice(&e.nodes);
            v.extend(e.nodes.iter().map(|&a| a + n));
            em.add_cell(kind, &v, z);
            em.add_face(&e.nodes, zmin);
            let top: Vec<u32> = e.nodes.iter().map(|&a| a + n).collect();
            em.add_face(&top, zmax);
        }
        for (b, (name, es)) in markers.iter().enumerate() {
            for e in es {
                if e.kind != LINE {
                    return Err(format!(
                        "marker '{name}': element type {} is not an edge of a 2D mesh (line 3)",
                        e.kind
                    ));
                }
                let (a, c) = (e.nodes[0], e.nodes[1]);
                em.add_face(&[a, c, c + n, a + n], b as u32);
            }
        }
        if markers.iter().any(|(m, _)| m == "zmin" || m == "zmax") {
            warnings.push(
                "a marker is called zmin or zmax, the names of the extruded planes; \
                 rename it with [mesh.boundaries]"
                    .into(),
            );
        }
    }
    let built = em.into_mesh()?;
    warnings.extend(built.warnings);
    Ok(MeshFile {
        mesh: built.mesh,
        reader: "su2",
        format: format!("SU2 {ndime}D ASCII"),
        warnings,
    })
}

/// Median edge length of the 2D cells.
fn median_edge(cells: &[Element], p: &[[f64; 3]]) -> Option<f64> {
    let mut l: Vec<f64> = Vec::new();
    for e in cells {
        let k = e.nodes.len();
        for i in 0..k {
            let (a, b) = (p[e.nodes[i] as usize], p[e.nodes[(i + 1) % k] as usize]);
            let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
            if d > 0.0 {
                l.push(d);
            }
        }
    }
    if l.is_empty() {
        return None;
    }
    l.sort_unstable_by(f64::total_cmp);
    Some(l[l.len() / 2])
}

/// Cells, nodes and markers of one zone of a file.
struct ZoneData {
    cells: Vec<Element>,
    points: Vec<[f64; 3]>,
    markers: Vec<(String, Vec<Element>)>,
}

/// Checks the sections of a zone and collects them.
fn finish_zone(
    ndime: Option<usize>,
    cells: Option<Vec<Element>>,
    points: Option<Vec<[f64; 3]>>,
    markers: Option<Vec<(String, Vec<Element>)>>,
    warnings: &mut Vec<String>,
) -> Result<ZoneData, String> {
    ndime.ok_or("no NDIME= (is this an SU2 mesh file?)")?;
    let cells = cells.ok_or("no NELEM= section")?;
    let points = points.ok_or("no NPOIN= section")?;
    let markers = markers.unwrap_or_else(|| {
        warnings.push("no NMARK= section: all boundary faces are in 'default'".into());
        Vec::new()
    });
    let np = points.len() as u32;
    let check = |e: &Element, what: &str| -> Result<(), String> {
        match e.nodes.iter().find(|&&n| n >= np) {
            Some(n) => Err(format!("{what} refers to node {n}, beyond the {np} nodes")),
            None => Ok(()),
        }
    };
    for e in &cells {
        check(e, "a cell")?;
    }
    for (name, es) in &markers {
        for e in es {
            check(e, &format!("an element of marker '{name}'"))?;
        }
    }
    Ok(ZoneData {
        cells,
        points,
        markers,
    })
}

/// Several zones (`NZONE= n`) into one mesh: the cells of zone k are the
/// cell zone `zone<k>`, the nodes of the zones are numbered one after the
/// other and stay apart unless `merge_nodes` joins coincident ones, markers
/// of the same name in several zones are one boundary.
fn multi_zone(
    zones: Vec<ZoneData>,
    ndime: usize,
    options: &ReadOptions,
    mut warnings: Vec<String>,
) -> Result<MeshFile, String> {
    if ndime != 3 {
        return Err(format!(
            "{} zones with NDIME= 2: only 3D multi-zone files are read",
            zones.len()
        ));
    }
    let mut em = ElementMesh {
        face_word: "marker elements",
        group_word: "marker",
        ..Default::default()
    };
    let n_zones = zones.len();
    for (k, zone) in zones.into_iter().enumerate() {
        let off = em.points.len() as u32;
        let z = em.zone(&format!("zone{}", k + 1));
        em.points.extend_from_slice(&zone.points);
        let mut v = Vec::with_capacity(8);
        for e in &zone.cells {
            let kind = match e.kind {
                TETRAHEDRON => CellKind::Tetrahedron,
                HEXAHEDRON => CellKind::Hexahedron,
                PRISM => CellKind::Prism,
                PYRAMID => CellKind::Pyramid,
                t => return Err(format!("element type {t} is not a 3D cell (NDIME= 3)")),
            };
            v.clear();
            v.extend(e.nodes.iter().map(|&a| a + off));
            em.add_cell(kind, &v, z);
        }
        for (name, es) in &zone.markers {
            let b = em.boundary(name);
            for e in es {
                if e.kind != TRIANGLE && e.kind != QUADRILATERAL {
                    return Err(format!(
                        "marker '{name}': element type {} is not a face of a 3D mesh \
                         (triangle 5, quadrilateral 9)",
                        e.kind
                    ));
                }
                v.clear();
                v.extend(e.nodes.iter().map(|&a| a + off));
                em.add_face(&v, b);
            }
        }
    }
    warnings.push(format!(
        "{n_zones} zones read as the cell zones zone1..zone{n_zones}{}",
        if options.merge_nodes {
            ""
        } else {
            "; their nodes are kept apart (mesh.merge_nodes joins coincident nodes)"
        }
    ));
    if options.merge_nodes {
        let removed = em.merge_points();
        warnings.push(format!("{removed} coincident nodes merged"));
    }
    let built = em.into_mesh()?;
    warnings.extend(built.warnings);
    Ok(MeshFile {
        mesh: built.mesh,
        reader: "su2",
        format: format!("SU2 {ndime}D ASCII, {n_zones} zones"),
        warnings,
    })
}
