//! Reader for the Fluent face-based `.msh` text format (D191).
//!
//! Fluent is a trademark of ANSYS, Inc.; the name is used here only to name
//! the file format. Foehncraft is not affiliated with or endorsed by ANSYS.
//!
//! Written from the public description of the format in the user guides of
//! the solvers that use it. A file is a list of sections in parentheses, the
//! numbers in hexadecimal:
//!
//! ```text
//! (0 "comment")
//! (2 3)                                  dimension
//! (10 (0 1 N 0 3))                       N nodes in all
//! (10 (zone first last type 3)(          nodes first..last
//!   x y z
//! ))
//! (12 (0 1 M 0 0))                       M cells in all
//! (12 (zone first last type elemtype))   cells first..last form a cell zone
//! (13 (0 1 F 0 0))                       F faces in all
//! (13 (zone first last bctype facetype)( faces first..last of a boundary zone
//!   n0 n1 .. c0 c1                       nodes, then the cells on both sides
//! ))
//! (45 (zone type name)())                the name of a zone
//! ```
//!
//! The mesh is face-based like ours: a face lists its nodes and the cells
//! `c0` and `c1` on its two sides (`c1` is 0 on the boundary), so polyhedral
//! cells are taken over as they are and the cell types (hexahedra,
//! tetrahedra, polyhedra, ...) do not matter. Face types 2 to 5 (fixed node
//! count) and 0 (the count leads each face, mixed and polygonal faces) are
//! read.
//!
//! - Zones: a face zone of the type `interior` (`bctype` 2) holds interior
//!   faces; every other face zone is a boundary of its name (`wall`,
//!   `velocity-inlet`, `symmetry`, ... all the same here: the boundary
//!   conditions of the case decide). A zone without a name from a `(45 ...)`
//!   section is called `zone<id>` (the numbers of `(45 ...)` are taken as hexadecimal
//!   first, else decimal). Cell zones are the cell sections
//!   (`fluid`, `solid`, named from `(45 ...)`).
//! - The side of a face is decided by geometry: the face is turned so that
//!   its normal points out of `c0`, so files with either orientation
//!   convention read alike.
//! - Only ASCII files are read (the binary sections 2010, 2012 and 2013 are
//!   refused with a hint), and 3D (2D files: extrude the mesh in the mesher).
//!   Periodic zones are ordinary boundaries.
//! - It shares the extension `.msh` with Gmsh: the first bytes decide.

use super::{MeshFile, MeshFormat, ReadOptions, declared};
use crate::builder::MeshBuilder;
use crate::geometry::{polygon_area_centroid, sub};
use std::collections::HashMap;

/// The format in the registry of [`crate::format`]; `mesh.format =
/// "fluent_msh"`.
#[derive(Debug, Clone, Copy)]
pub struct IndustryMsh;

impl MeshFormat for IndustryMsh {
    fn name(&self) -> &'static str {
        "fluent_msh"
    }

    fn title(&self) -> &'static str {
        "Fluent face-based .msh, ASCII, 3D, polyhedral"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &["msh"]
    }

    fn sniff(&self, head: &[u8]) -> bool {
        let start = head
            .iter()
            .position(|c| !c.is_ascii_whitespace())
            .unwrap_or(0);
        let h = &head[start..];
        h.starts_with(b"(0 ") || h.starts_with(b"(1 ") || h.starts_with(b"(2 ")
    }

    fn parse(&self, bytes: &[u8], _options: &ReadOptions) -> Result<MeshFile, String> {
        parse(bytes)
    }
}

/// The text between the parenthesis at `at` and its partner, and the
/// position after the partner. Quoted strings may hold parentheses.
fn balanced(t: &str, at: usize) -> Result<(&str, usize), String> {
    let b = t.as_bytes();
    let (mut depth, mut i, mut quoted) = (0usize, at, false);
    while i < b.len() {
        match b[i] {
            b'"' => quoted = !quoted,
            b'(' if !quoted => depth += 1,
            b')' if !quoted => {
                depth -= 1;
                if depth == 0 {
                    return Ok((&t[at + 1..i], i + 1));
                }
            }
            _ => {}
        }
        i += 1;
    }
    Err(format!("line {}: a section is not closed", line_of(t, at)))
}

fn line_of(t: &str, at: usize) -> usize {
    1 + t.as_bytes()[..at.min(t.len())]
        .iter()
        .filter(|&&c| c == b'\n')
        .count()
}

fn hex(tok: &str, what: &str) -> Result<usize, String> {
    usize::from_str_radix(tok, 16)
        .map_err(|_| format!("{what}: '{tok}' is not a hexadecimal number"))
}

/// A section: its index, its header list and its data list (if any).
struct Section<'a> {
    index: usize,
    header: Vec<&'a str>,
    data: Option<&'a str>,
    line: usize,
}

fn section(body: &str, line: usize) -> Result<Section<'_>, String> {
    let body = body.trim_start();
    let end = body
        .find(|c: char| c.is_ascii_whitespace() || c == '(')
        .unwrap_or(body.len());
    let index: usize = body[..end]
        .parse()
        .map_err(|_| format!("line {line}: a section needs a number first"))?;
    let rest = &body[end..];
    let header;
    let mut data = None;
    match rest.find('(') {
        Some(h) if index != 0 && index != 1 => {
            let (head, after) = balanced(rest, h)?;
            header = head.split_ascii_whitespace().collect();
            let tail = &rest[after..];
            if let Some(d) = tail.find('(') {
                data = Some(balanced(tail, d)?.0);
            }
        }
        _ => header = rest.split_ascii_whitespace().collect(),
    }
    Ok(Section {
        index,
        header,
        data,
        line,
    })
}

#[derive(Default)]
struct Face {
    zone: usize,
    nodes: Vec<u32>,
    c0: usize,
    c1: usize,
}

fn parse(bytes: &[u8]) -> Result<MeshFile, String> {
    let t = std::str::from_utf8(bytes)
        .map_err(|_| "not a text file: only the ASCII form of the format is read".to_string())?;
    let mut warnings = Vec::new();

    let (mut n_nodes, mut n_cells, mut n_faces) = (0usize, 0usize, 0usize);
    let mut dim = 3usize;
    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut cell_ranges: Vec<(usize, usize, usize)> = Vec::new(); // zone, first, last
    let mut face_zone_type: HashMap<usize, usize> = HashMap::new();
    let mut faces: Vec<Face> = Vec::new();
    let mut names: HashMap<String, String> = HashMap::new();
    let mut periodic = false;

    let mut pos = 0usize;
    while let Some(off) = t[pos..].find(|c: char| !c.is_ascii_whitespace()) {
        let at = pos + off;
        if t.as_bytes()[at] != b'(' {
            return Err(format!(
                "line {}: expected a section in parentheses",
                line_of(t, at)
            ));
        }
        let (body, next) = balanced(t, at)?;
        pos = next;
        let s = section(body, line_of(t, at))?;
        let line = s.line;
        match s.index {
            2 => {
                dim = s
                    .header
                    .first()
                    .and_then(|v| v.parse().ok())
                    .ok_or_else(|| format!("line {line}: the dimension is missing"))?;
                if dim != 3 {
                    return Err(format!(
                        "the file is {dim}D: only 3D meshes are read (extrude a 2D mesh in the \
                         mesher)"
                    ));
                }
            }
            10 => {
                let h = |i: usize| -> Result<usize, String> {
                    hex(
                        s.header.get(i).unwrap_or(&""),
                        &format!("line {line}: node section"),
                    )
                };
                let (zone, first, last) = (h(0)?, h(1)?, h(2)?);
                if let Some(nd) = s.header.get(4).and_then(|v| v.parse::<usize>().ok())
                    && nd != 3
                {
                    return Err(format!(
                        "line {line}: nodes with {nd} coordinates: only 3D meshes are read"
                    ));
                }
                let last = declared(last, t.len(), &format!("line {line}: nodes"))?;
                if zone == 0 {
                    n_nodes = last;
                    points.resize(n_nodes, [0.0; 3]);
                    continue;
                }
                let Some(d) = s.data else { continue };
                if last > points.len() {
                    if n_nodes == 0 {
                        points.resize(last, [0.0; 3]);
                    } else {
                        return Err(format!(
                            "line {line}: nodes up to {last} but only {n_nodes} are declared"
                        ));
                    }
                }
                let mut it = d.split_ascii_whitespace();
                for i in first.max(1)..=last {
                    let mut p = [0.0; 3];
                    for c in &mut p {
                        *c = it
                            .next()
                            .ok_or_else(|| format!("line {line}: too few node coordinates"))?
                            .parse()
                            .map_err(|_| {
                                format!("line {line}: a node coordinate is not a number")
                            })?;
                    }
                    points[i - 1] = p;
                }
            }
            12 => {
                let h = |i: usize| -> Result<usize, String> {
                    hex(
                        s.header.get(i).unwrap_or(&""),
                        &format!("line {line}: cell section"),
                    )
                };
                let (zone, first, last) = (h(0)?, h(1)?, h(2)?);
                let last = declared(last, t.len(), &format!("line {line}: cells"))?;
                if zone == 0 {
                    n_cells = last;
                    continue;
                }
                if let Some(ty) = s
                    .header
                    .get(4)
                    .and_then(|v| usize::from_str_radix(v, 16).ok())
                    && (ty == 1 || ty == 3)
                {
                    return Err(format!(
                        "line {line}: triangular or quadrilateral cells: only 3D meshes are read"
                    ));
                }
                cell_ranges.push((zone, first, last));
            }
            13 => {
                let h = |i: usize| -> Result<usize, String> {
                    hex(
                        s.header.get(i).unwrap_or(&""),
                        &format!("line {line}: face section"),
                    )
                };
                let (zone, _first, last) = (h(0)?, h(1)?, h(2)?);
                if zone == 0 {
                    n_faces = last;
                    continue;
                }
                let (bc, ftype) = (h(3)?, h(4)?);
                face_zone_type.insert(zone, bc);
                let Some(d) = s.data else { continue };
                let count = (last + 1)
                    .checked_sub(h(1)?)
                    .ok_or_else(|| format!("line {line}: a face section ends before it starts"))?;
                let count = declared(count, d.len(), &format!("line {line}: faces"))?;
                let mut it = d.split_ascii_whitespace();
                faces.reserve(count);
                let mut num = |what: &str| -> Result<usize, String> {
                    let tok = it
                        .next()
                        .ok_or_else(|| format!("line {line}: a face ends too early ({what})"))?;
                    hex(tok, &format!("line {line}"))
                };
                for _ in 0..count {
                    let nn = match ftype {
                        0 | 5 => num("node count")?,
                        2..=4 => ftype,
                        _ => {
                            return Err(format!("line {line}: unknown face type {ftype}"));
                        }
                    };
                    if !(2..=64).contains(&nn) {
                        return Err(format!("line {line}: a face with {nn} nodes"));
                    }
                    let mut nodes = Vec::with_capacity(nn);
                    for _ in 0..nn {
                        let n = num("node")?;
                        if n == 0 || n > u32::MAX as usize {
                            return Err(format!("line {line}: a face refers to node {n}"));
                        }
                        nodes.push((n - 1) as u32);
                    }
                    let (c0, c1) = (num("cell")?, num("cell")?);
                    faces.push(Face {
                        zone,
                        nodes,
                        c0,
                        c1,
                    });
                }
            }
            18 => periodic = true,
            2010 | 2012 | 2013 => {
                return Err(format!(
                    "line {line}: a binary section ({}): only ASCII files are read (write the \
                     mesh as ASCII)",
                    s.index
                ));
            }
            39 | 45 => {
                // (45 (id type name)()): the type may be left out
                if let (Some(id), Some(name)) = (s.header.first(), s.header.get(2)) {
                    names.insert(id.to_ascii_lowercase(), name.to_string());
                }
            }
            _ => {}
        }
    }
    let _ = dim;

    if faces.is_empty() {
        return Err("the file has no faces (section 13)".into());
    }
    if n_faces != 0 && n_faces != faces.len() {
        warnings.push(format!(
            "{n_faces} faces are declared but {} are given",
            faces.len()
        ));
    }
    let max_cell = faces.iter().map(|f| f.c0.max(f.c1)).max().unwrap_or(0);
    if n_cells == 0 {
        n_cells = declared(max_cell, t.len(), "cells")?;
    }
    if max_cell > n_cells {
        return Err(format!(
            "a face refers to cell {max_cell}, beyond the {n_cells} cells"
        ));
    }
    if n_cells == 0 {
        return Err("the mesh has no cells".into());
    }
    if let Some(f) = faces
        .iter()
        .find(|f| f.nodes.iter().any(|&n| n as usize >= points.len()))
    {
        return Err(format!(
            "a face of zone {} refers to a node beyond the {} nodes",
            f.zone,
            points.len()
        ));
    }

    // cell centres: the mean of the centres of the faces around a cell
    let fcentre: Vec<[f64; 3]> = faces
        .iter()
        .map(|f| {
            let mut c = [0.0; 3];
            for &n in &f.nodes {
                for (a, b) in c.iter_mut().zip(points[n as usize]) {
                    *a += b / f.nodes.len() as f64;
                }
            }
            c
        })
        .collect();
    let mut centre = vec![[0.0f64; 3]; n_cells];
    let mut count = vec![0usize; n_cells];
    for (f, fc) in faces.iter().zip(&fcentre) {
        for c in [f.c0, f.c1] {
            if c > 0 {
                for (a, b) in centre[c - 1].iter_mut().zip(fc) {
                    *a += b;
                }
                count[c - 1] += 1;
            }
        }
    }
    for (c, n) in centre.iter_mut().zip(&count) {
        if *n == 0 {
            continue;
        }
        for a in c.iter_mut() {
            *a /= *n as f64;
        }
    }
    if let Some(c) = count.iter().position(|&n| n == 0) {
        return Err(format!("cell {} has no faces", c + 1));
    }

    // the zone numbers of the sections 10 to 13 are hexadecimal; the files
    // write them in (45 ...) in either form (the same below 10)
    let name_of = |z: usize| -> Option<String> {
        names
            .get(&format!("{z:x}"))
            .or_else(|| names.get(&z.to_string()))
            .cloned()
    };
    let zone_name = |z: usize| -> String { name_of(z).unwrap_or_else(|| format!("zone{z}")) };

    // face zones in order of first use, so boundaries appear in file order
    let mut b = MeshBuilder::new(points.clone(), n_cells);
    let mut declared: Vec<usize> = Vec::new();
    for f in &faces {
        let boundary = f.c1 == 0 && f.c0 != 0 || face_zone_type.get(&f.zone) != Some(&2);
        if boundary && !declared.contains(&f.zone) {
            declared.push(f.zone);
            b.declare_boundary(&zone_name(f.zone));
        }
    }
    let mut baffles: Vec<usize> = Vec::new();
    for (f, fc) in faces.iter().zip(&fcentre) {
        let (mut c0, mut c1) = (f.c0, f.c1);
        if c0 == 0 && c1 == 0 {
            return Err(format!("a face of zone {} has no cell", f.zone));
        }
        if c0 == 0 {
            std::mem::swap(&mut c0, &mut c1);
        }
        if c0 == c1 {
            return Err(format!("a face has the cell {c0} on both sides"));
        }
        // turn the face so that its normal points out of c0
        let mut nodes = f.nodes.clone();
        let (n, _) = polygon_area_centroid(&points, &nodes);
        let out = sub(*fc, centre[c0 - 1]);
        if n[0] * out[0] + n[1] * out[1] + n[2] * out[2] < 0.0 {
            nodes.reverse();
        }
        let interior_zone = face_zone_type.get(&f.zone) == Some(&2);
        if c1 == 0 {
            if interior_zone {
                return Err(format!(
                    "a face of the interior zone '{}' has only one cell",
                    zone_name(f.zone)
                ));
            }
            b.add_boundary(&nodes, (c0 - 1) as u32, &zone_name(f.zone));
        } else {
            if !interior_zone && !baffles.contains(&f.zone) {
                baffles.push(f.zone);
            }
            b.add_interior(&nodes, (c0 - 1) as u32, (c1 - 1) as u32);
        }
    }
    for z in &baffles {
        warnings.push(format!(
            "zone '{}' has faces between two cells (a baffle or interface): they are interior \
             faces",
            zone_name(*z)
        ));
    }
    if periodic {
        warnings.push(
            "periodic zones are read as ordinary boundaries; pair them with `periodic` in the \
             case"
                .into(),
        );
    }

    // cell zones
    if !cell_ranges.is_empty() {
        let mut order: Vec<usize> = Vec::new();
        for r in &cell_ranges {
            if !order.contains(&r.0) {
                order.push(r.0);
            }
        }
        let mut zname: Vec<String> = order
            .iter()
            .map(|&z| {
                name_of(z).unwrap_or_else(|| {
                    if order.len() == 1 {
                        "fluid".into()
                    } else {
                        format!("zone{z}")
                    }
                })
            })
            .collect();
        let unset = u32::MAX;
        let mut cz = vec![unset; n_cells];
        for &(z, first, last) in &cell_ranges {
            let zi = order.iter().position(|&o| o == z).unwrap() as u32;
            if last > n_cells {
                return Err(format!("a cell zone reaches cell {last}, beyond {n_cells}"));
            }
            for c in first.max(1)..=last {
                cz[c - 1] = zi;
            }
        }
        let missing = cz.iter().filter(|&&z| z == unset).count();
        if missing > 0 {
            zname.push("default".into());
            let d = (zname.len() - 1) as u32;
            for z in &mut cz {
                if *z == unset {
                    *z = d;
                }
            }
            warnings.push(format!(
                "{missing} cells are in no cell zone and form the zone 'default'"
            ));
        }
        b.set_zones(zname, cz);
    }

    let mut mesh = b.finish()?;
    mesh.boundaries.retain(|p| p.len > 0);
    Ok(MeshFile {
        mesh,
        reader: "fluent_msh",
        format: "Fluent .msh ASCII 3D".into(),
        warnings,
    })
}
