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
//! - ASCII and binary files are read, 3D only (2D files: extrude the mesh in
//!   the mesher). Periodic zones are ordinary boundaries.
//! - Binary files: a section of index `2000 + i` (single precision) or
//!   `3000 + i` (double precision) is the binary form of the text section
//!   `i`. Its header is text as above; its data block holds the numbers as
//!   little-endian binary, coordinates as 32-bit or 64-bit reals, everything
//!   else (node and cell numbers, node counts of mixed faces, element types)
//!   as 32-bit integers; the words `End of Binary Section 20xx` close it:
//!
//!   ```text
//!   (3010 (1 1 N 1 3)(<N x 3 reals>)End of Binary Section 3010)
//!   (2012 (1 1 M 1 0)(<M element types>)End of Binary Section 2012)
//!   (2013 (3 1 F 2 0)(<n n0 n1 .. c0 c1, per face>)End of Binary Section 2013)
//!   ```
//!
//!   Sections 2010/3010, 2012/3012, 2013/3013 and 2018/3018 (periodic face
//!   pairs) are read by their counts; other binary sections are skipped up
//!   to their closing words. ASCII and binary sections may be mixed in one
//!   file.
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
        "Fluent face-based .msh, ASCII and binary, 3D, polyhedral"
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

/// The bytes between the parenthesis at `at` and its partner, and the
/// position after the partner. Quoted strings may hold parentheses.
fn balanced(t: &[u8], at: usize) -> Result<(&[u8], usize), String> {
    let (mut depth, mut i, mut quoted) = (0usize, at, false);
    while i < t.len() {
        match t[i] {
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

fn line_of(t: &[u8], at: usize) -> usize {
    1 + t[..at.min(t.len())].iter().filter(|&&c| c == b'\n').count()
}

fn hex(tok: &str, what: &str) -> Result<usize, String> {
    usize::from_str_radix(tok, 16)
        .map_err(|_| format!("{what}: '{tok}' is not a hexadecimal number"))
}

/// The leading decimal number of a section (its index) and the position
/// after it.
fn index_at(t: &[u8], at: usize) -> Option<(usize, usize)> {
    let mut i = at;
    while i < t.len() && t[i].is_ascii_whitespace() {
        i += 1;
    }
    let s = i;
    while i < t.len() && t[i].is_ascii_digit() {
        i += 1;
    }
    let n = std::str::from_utf8(&t[s..i]).ok()?.parse().ok()?;
    Some((n, i))
}

/// A text section: its index, its header list and its data list (if any).
struct Section<'a> {
    index: usize,
    header: Vec<&'a str>,
    data: Option<&'a str>,
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
            let (head, after) = balanced(rest.as_bytes(), h)?;
            header = std::str::from_utf8(head)
                .map_err(|_| format!("line {line}: the header is not text"))?
                .split_ascii_whitespace()
                .collect();
            let tail = &rest[after..];
            if let Some(d) = tail.find('(') {
                let (bytes, _) = balanced(tail.as_bytes(), d)?;
                data = Some(
                    std::str::from_utf8(bytes)
                        .map_err(|_| format!("line {line}: the data is not text"))?,
                );
            }
        }
        _ => header = rest.split_ascii_whitespace().collect(),
    }
    Ok(Section {
        index,
        header,
        data,
    })
}

#[derive(Default)]
struct Face {
    zone: usize,
    nodes: Vec<u32>,
    c0: usize,
    c1: usize,
}

/// Little-endian numbers of a binary data block, read one after the other.
struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
    what: String,
}

impl Cursor<'_> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N], String> {
        let end = self.pos + N;
        let b = self
            .bytes
            .get(self.pos..end)
            .ok_or_else(|| format!("{}: the binary data ends too early", self.what))?;
        self.pos = end;
        let mut a = [0u8; N];
        a.copy_from_slice(b);
        Ok(a)
    }

    /// A 32-bit integer, refused when negative.
    fn int(&mut self) -> Result<usize, String> {
        let v = i32::from_le_bytes(self.take()?);
        usize::try_from(v).map_err(|_| format!("{}: a negative number ({v})", self.what))
    }

    fn real(&mut self, double: bool) -> Result<f64, String> {
        Ok(if double {
            f64::from_le_bytes(self.take()?)
        } else {
            f32::from_le_bytes(self.take()?) as f64
        })
    }
}

/// What the sections give, gathered over the file.
#[derive(Default)]
struct Gathered {
    n_nodes: usize,
    n_cells: usize,
    n_faces: usize,
    points: Vec<[f64; 3]>,
    /// zone, first, last
    cell_ranges: Vec<(usize, usize, usize)>,
    face_zone_type: HashMap<usize, usize>,
    faces: Vec<Face>,
    names: HashMap<String, String>,
    periodic: bool,
}

/// The header numbers `zone first last type nd` of a node, cell or face
/// section.
fn header_hex(header: &[&str], i: usize, place: &str) -> Result<usize, String> {
    hex(header.get(i).copied().unwrap_or(""), place)
}

impl Gathered {
    /// Nodes `first..=last` of a node section; `None` when the section only
    /// declares the count (zone 0) and so has no data.
    fn node_range(
        &mut self,
        header: &[&str],
        len: usize,
        place: &str,
    ) -> Result<Option<(usize, usize)>, String> {
        let (zone, first, last) = (
            header_hex(header, 0, place)?,
            header_hex(header, 1, place)?,
            header_hex(header, 2, place)?,
        );
        if let Some(nd) = header.get(4).and_then(|v| v.parse::<usize>().ok())
            && nd != 3
        {
            return Err(format!(
                "{place}: nodes with {nd} coordinates: only 3D meshes are read"
            ));
        }
        let last = declared(last, len, &format!("{place}: nodes"))?;
        if zone == 0 {
            self.n_nodes = last;
            self.points.resize(self.n_nodes, [0.0; 3]);
            return Ok(None);
        }
        if last > self.points.len() {
            if self.n_nodes == 0 {
                self.points.resize(last, [0.0; 3]);
            } else {
                return Err(format!(
                    "{place}: nodes up to {last} but only {} are declared",
                    self.n_nodes
                ));
            }
        }
        Ok(Some((first.max(1), last)))
    }

    /// A cell section; returns its zone, its count and whether its cells
    /// are of mixed type (then the element types follow as data).
    fn cells(
        &mut self,
        header: &[&str],
        len: usize,
        place: &str,
    ) -> Result<(usize, usize, bool), String> {
        let (zone, first, last) = (
            header_hex(header, 0, place)?,
            header_hex(header, 1, place)?,
            header_hex(header, 2, place)?,
        );
        let last = declared(last, len, &format!("{place}: cells"))?;
        if zone == 0 {
            self.n_cells = last;
            return Ok((0, 0, false));
        }
        let ty = header
            .get(4)
            .and_then(|v| usize::from_str_radix(v, 16).ok());
        if matches!(ty, Some(1 | 3)) {
            return Err(format!(
                "{place}: triangular or quadrilateral cells: only 3D meshes are read"
            ));
        }
        self.cell_ranges.push((zone, first, last));
        Ok((zone, (last + 1).saturating_sub(first), ty == Some(0)))
    }

    /// The header of a face section: `(zone, count, face type)`, or `None`
    /// for the declaration of the count (zone 0).
    fn face_header(
        &mut self,
        header: &[&str],
        place: &str,
    ) -> Result<Option<(usize, usize, usize)>, String> {
        let (zone, first, last) = (
            header_hex(header, 0, place)?,
            header_hex(header, 1, place)?,
            header_hex(header, 2, place)?,
        );
        if zone == 0 {
            self.n_faces = last;
            return Ok(None);
        }
        let (bc, ftype) = (header_hex(header, 3, place)?, header_hex(header, 4, place)?);
        self.face_zone_type.insert(zone, bc);
        let count = (last + 1)
            .checked_sub(first)
            .ok_or_else(|| format!("{place}: a face section ends before it starts"))?;
        Ok(Some((zone, count, ftype)))
    }

    /// `count` faces of `zone`, their numbers from `num` (hexadecimal text
    /// or binary integers).
    fn faces(
        &mut self,
        zone: usize,
        count: usize,
        ftype: usize,
        place: &str,
        num: &mut dyn FnMut(&str) -> Result<usize, String>,
    ) -> Result<(), String> {
        self.faces.reserve(count);
        for _ in 0..count {
            let nn = match ftype {
                0 | 5 => num("node count")?,
                2..=4 => ftype,
                _ => return Err(format!("{place}: unknown face type {ftype}")),
            };
            if !(2..=64).contains(&nn) {
                return Err(format!("{place}: a face with {nn} nodes"));
            }
            let mut nodes = Vec::with_capacity(nn);
            for _ in 0..nn {
                let n = num("node")?;
                if n == 0 || n > u32::MAX as usize {
                    return Err(format!("{place}: a face refers to node {n}"));
                }
                nodes.push((n - 1) as u32);
            }
            let (c0, c1) = (num("cell")?, num("cell")?);
            self.faces.push(Face {
                zone,
                nodes,
                c0,
                c1,
            });
        }
        Ok(())
    }

    /// A text section.
    fn text(&mut self, s: &Section<'_>, len: usize, line: usize) -> Result<(), String> {
        let place = format!("line {line}");
        match s.index {
            2 => {
                let dim: usize = s
                    .header
                    .first()
                    .and_then(|v| v.parse().ok())
                    .ok_or_else(|| format!("{place}: the dimension is missing"))?;
                if dim != 3 {
                    return Err(format!(
                        "the file is {dim}D: only 3D meshes are read (extrude a 2D mesh in the \
                         mesher)"
                    ));
                }
            }
            10 => {
                let Some((first, last)) = self.node_range(&s.header, len, &place)? else {
                    return Ok(());
                };
                let Some(d) = s.data else { return Ok(()) };
                let mut it = d.split_ascii_whitespace();
                for i in first..=last {
                    let mut p = [0.0; 3];
                    for c in &mut p {
                        *c = it
                            .next()
                            .ok_or_else(|| format!("{place}: too few node coordinates"))?
                            .parse()
                            .map_err(|_| format!("{place}: a node coordinate is not a number"))?;
                    }
                    self.points[i - 1] = p;
                }
            }
            12 => {
                self.cells(&s.header, len, &place)?;
            }
            13 => {
                let Some((zone, count, ftype)) = self.face_header(&s.header, &place)? else {
                    return Ok(());
                };
                let Some(d) = s.data else { return Ok(()) };
                let count = declared(count, d.len(), &format!("{place}: faces"))?;
                let mut it = d.split_ascii_whitespace();
                let mut num = |what: &str| -> Result<usize, String> {
                    let tok = it
                        .next()
                        .ok_or_else(|| format!("{place}: a face ends too early ({what})"))?;
                    hex(tok, &place)
                };
                self.faces(zone, count, ftype, &place, &mut num)?;
            }
            18 => self.periodic = true,
            39 | 45 => {
                // (45 (id type name)()): the type may be left out
                if let (Some(id), Some(name)) = (s.header.first(), s.header.get(2)) {
                    self.names.insert(id.to_ascii_lowercase(), name.to_string());
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// A binary section `(2xxx (header)(data)End of Binary Section 2xxx)`
    /// (single precision) or `3xxx` (double precision) that starts at `at`
    /// with its index ending at `after`; returns the position after it.
    fn binary(&mut self, t: &[u8], at: usize, index: usize, after: usize) -> Result<usize, String> {
        let place = format!("section {index} at byte {at}");
        let double = index / 1000 == 3;
        let kind = index % 1000;
        // the header: text in parentheses
        let h = after
            + t[after..]
                .iter()
                .position(|c| !c.is_ascii_whitespace())
                .ok_or_else(|| format!("{place}: the section is not closed"))?;
        if t[h] != b'(' {
            return Err(format!("{place}: the header is missing"));
        }
        let (head, mut pos) = balanced(t, h)?;
        let header: Vec<&str> = std::str::from_utf8(head)
            .map_err(|_| format!("{place}: the header is not text"))?
            .split_ascii_whitespace()
            .collect();
        // the data block, if any: '(' bytes ')'
        while pos < t.len() && t[pos].is_ascii_whitespace() {
            pos += 1;
        }
        if pos < t.len() && t[pos] == b'(' {
            let mut c = Cursor {
                bytes: &t[pos + 1..],
                pos: 0,
                what: place.clone(),
            };
            match kind {
                10 => {
                    if let Some((first, last)) = self.node_range(&header, t.len(), &place)? {
                        let need =
                            (last + 1).saturating_sub(first) * 3 * if double { 8 } else { 4 };
                        declared(need, c.bytes.len(), &format!("{place}: nodes"))?;
                        for i in first..=last {
                            for x in &mut self.points[i - 1] {
                                *x = c.real(double)?;
                            }
                        }
                    }
                }
                12 => {
                    let (zone, count, mixed) = self.cells(&header, t.len(), &place)?;
                    if zone != 0 && mixed {
                        // the element type of every cell: not needed
                        let need = count
                            .checked_mul(4)
                            .ok_or_else(|| format!("{place}: too many cells"))?;
                        declared(need, c.bytes.len(), &format!("{place}: cells"))?;
                        c.pos = need;
                    }
                }
                13 => {
                    if let Some((zone, count, ftype)) = self.face_header(&header, &place)? {
                        // a face takes at least 4 integers of 4 bytes
                        declared(count, c.bytes.len() / 16, &format!("{place}: faces"))?;
                        self.faces(zone, count, ftype, &place, &mut |_| c.int())?;
                    }
                }
                18 => {
                    // periodic face pairs: two integers each
                    self.periodic = true;
                    let (first, last) = (
                        header_hex(&header, 0, &place)?,
                        header_hex(&header, 1, &place)?,
                    );
                    let need = (last + 1)
                        .saturating_sub(first)
                        .checked_mul(8)
                        .ok_or_else(|| format!("{place}: too many periodic faces"))?;
                    declared(need, c.bytes.len(), &format!("{place}: periodic faces"))?;
                    c.pos = need;
                }
                _ => {
                    // a section not needed here: its data end before the
                    // closing words
                    let mark = b"End of Binary Section";
                    let end = c
                        .bytes
                        .windows(mark.len())
                        .position(|w| w == mark)
                        .ok_or_else(|| {
                            format!("{place}: the end of the binary section is not found")
                        })?;
                    c.pos = c.bytes[..end]
                        .iter()
                        .rposition(|&b| b == b')')
                        .ok_or_else(|| format!("{place}: the data block is not closed"))?;
                }
            }
            pos += 1 + c.pos;
            if t.get(pos) != Some(&b')') {
                return Err(format!(
                    "{place}: the binary data do not end where the header says"
                ));
            }
            pos += 1;
        } else {
            match kind {
                10 => {
                    self.node_range(&header, t.len(), &place)?;
                }
                12 => {
                    self.cells(&header, t.len(), &place)?;
                }
                13 => {
                    self.face_header(&header, &place)?;
                }
                _ => {}
            }
        }
        // the closing words "End of Binary Section 2xxx" and ')'
        let close = t[pos..]
            .iter()
            .position(|&b| b == b')')
            .ok_or_else(|| format!("{place}: the section is not closed"))?;
        if t[pos..pos + close].contains(&b'(') {
            return Err(format!("{place}: the section is not closed"));
        }
        Ok(pos + close + 1)
    }
}

fn parse(t: &[u8]) -> Result<MeshFile, String> {
    let mut warnings = Vec::new();
    let mut g = Gathered::default();
    let mut binary = false;
    let mut double = false;

    let mut pos = 0usize;
    while let Some(off) = t[pos..].iter().position(|c| !c.is_ascii_whitespace()) {
        let at = pos + off;
        if t[at] != b'(' {
            return Err(format!(
                "line {}: expected a section in parentheses",
                line_of(t, at)
            ));
        }
        match index_at(t, at + 1) {
            Some((index, after)) if (2000..4000).contains(&index) => {
                binary = true;
                double |= index >= 3000;
                pos = g.binary(t, at, index, after)?;
            }
            Some((0 | 1, _)) => {
                // comments and headers: any bytes
                pos = balanced(t, at)?.1;
            }
            _ => {
                let (body, next) = balanced(t, at)?;
                pos = next;
                let line = line_of(t, at);
                let body = std::str::from_utf8(body)
                    .map_err(|_| format!("line {line}: the section is not text"))?;
                let s = section(body, line)?;
                g.text(&s, t.len(), line)?;
            }
        }
    }
    let Gathered {
        mut n_cells,
        n_faces,
        points,
        cell_ranges,
        face_zone_type,
        faces,
        names,
        periodic,
        ..
    } = g;

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
        format: match (binary, double) {
            (false, _) => "Fluent .msh ASCII 3D",
            (true, false) => "Fluent .msh binary single precision 3D",
            (true, true) => "Fluent .msh binary double precision 3D",
        }
        .into(),
        warnings,
    })
}
