//! Reader for OpenFOAM polyMesh directories (D187).
//!
//! Written from the public description of the format in the OpenFOAM user
//! guide ("Mesh description", "Boundaries"): a mesh is the directory
//! `constant/polyMesh` of a case with the files
//!
//! ```text
//! points      N points, (x y z) each
//! faces       M faces, k(n1 ... nk) each in ASCII; in binary the compact
//!             form (class faceCompactList): M+1 offsets, then the node labels
//! owner       M cell labels: the cell on the side of the face normal's tail
//! neighbour   one label per internal face (the cell the normal points to)
//! boundary    P patches, name { type ...; nFaces n; startFace s; }
//! cellZones   optional: name { type cellZone; cellLabels List<label> n(...); }
//! ```
//!
//! The faces come in the order internal (in `neighbour`), then patch after
//! patch (`startFace`, `nFaces`). The normal of an internal face points from
//! its owner to its neighbour, that of a boundary face out of the domain; the
//! mesh is polyhedral as it stands, so the cells are taken over as they are
//! (no cell shapes, no matching of faces: coincident faces of a baffle stay
//! two boundary faces).
//!
//! - ASCII and binary files (`format binary`, with the `arch` entry of the
//!   header: `label=32` or `label=64`, `scalar=32` or `scalar=64`, `LSB`
//!   or `MSB`), also mixed; comments in the files are skipped.
//! - The path can be the case directory, `constant`, `polyMesh` itself or a
//!   `<time>/polyMesh` of a moving mesh.
//! - Patches become boundaries of their names, whatever their `type`
//!   (`patch`, `wall`, `empty`, `symmetry`, `wedge`, `cyclic`): the boundary
//!   conditions of the case decide what a boundary is. A 2D OpenFOAM mesh,
//!   one cell thick with its `empty` patch, is read as the 3D mesh it is.
//!   Processor patches (a decomposed case) are refused: reconstruct first.
//! - `cellZones` give the cell zones (cells in none are in `default`).
//!   Compressed files (`points.gz`) are not read: `gunzip` them.

use super::{MeshFile, MeshFormat, ReadOptions, declared};
use crate::builder::MeshBuilder;
use std::path::{Path, PathBuf};

/// The OpenFOAM format in the registry of [`crate::format`].
#[derive(Debug, Clone, Copy)]
pub struct OpenFoam;

impl MeshFormat for OpenFoam {
    fn name(&self) -> &'static str {
        "openfoam"
    }

    fn title(&self) -> &'static str {
        "OpenFOAM polyMesh directory, ASCII and binary, polyhedral"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &[]
    }

    fn sniff(&self, _head: &[u8]) -> bool {
        false
    }

    fn accepts_dir(&self, dir: &Path) -> bool {
        poly_dir(dir).is_some()
    }

    fn read_path(&self, path: &Path, _options: &ReadOptions) -> Result<MeshFile, String> {
        let dir = poly_dir(path).ok_or_else(|| {
            format!(
                "{} holds no polyMesh (the files points, faces, owner, neighbour and boundary \
                 in the directory itself, in polyMesh/ or in constant/polyMesh/)",
                path.display()
            )
        })?;
        read_dir(&dir)
    }
}

/// The directory with the polyMesh files of `path`.
fn poly_dir(path: &Path) -> Option<PathBuf> {
    let has = |d: &Path| {
        ["points", "faces", "owner"]
            .iter()
            .all(|f| d.join(f).is_file() || d.join(format!("{f}.gz")).is_file())
    };
    [
        path.to_path_buf(),
        path.join("polyMesh"),
        path.join("constant").join("polyMesh"),
    ]
    .into_iter()
    .find(|d| has(d))
}

/// Header of a file: the `FoamFile` dictionary.
#[derive(Debug, Default)]
struct Header {
    binary: bool,
    class: String,
    label_bytes: usize,
    scalar_bytes: usize,
    big_endian: bool,
}

struct Scan<'a> {
    b: &'a [u8],
    pos: usize,
    head: Header,
    what: &'static str,
}

impl<'a> Scan<'a> {
    fn new(b: &'a [u8], what: &'static str) -> Result<Self, String> {
        let mut s = Scan {
            b,
            pos: 0,
            head: Header {
                label_bytes: 4,
                scalar_bytes: 8,
                ..Default::default()
            },
            what,
        };
        s.header()?;
        Ok(s)
    }

    fn err(&self, msg: impl std::fmt::Display) -> String {
        let line = 1 + self.b[..self.pos.min(self.b.len())]
            .iter()
            .filter(|&&c| c == b'\n')
            .count();
        format!("{} line {line}: {msg}", self.what)
    }

    /// Skips white space and comments.
    fn skip(&mut self) {
        loop {
            while self.pos < self.b.len() && self.b[self.pos].is_ascii_whitespace() {
                self.pos += 1;
            }
            if self.b[self.pos..].starts_with(b"//") {
                while self.pos < self.b.len() && self.b[self.pos] != b'\n' {
                    self.pos += 1;
                }
            } else if self.b[self.pos..].starts_with(b"/*") {
                match self.b[self.pos + 2..].windows(2).position(|w| w == b"*/") {
                    Some(k) => self.pos += k + 4,
                    None => self.pos = self.b.len(),
                }
            } else {
                return;
            }
        }
    }

    fn peek(&mut self) -> Option<u8> {
        self.skip();
        self.b.get(self.pos).copied()
    }

    fn eat(&mut self, c: u8) -> bool {
        if self.peek() == Some(c) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn expect(&mut self, c: u8) -> Result<(), String> {
        if self.eat(c) {
            Ok(())
        } else {
            let found = self
                .b
                .get(self.pos)
                .map_or("the end of the file".to_string(), |&x| {
                    format!("'{}'", x as char)
                });
            Err(self.err(format!("expected '{}', found {found}", c as char)))
        }
    }

    /// A word: up to white space or one of `(){};`; a quoted string is one
    /// word without its quotes.
    fn word(&mut self) -> Result<&'a str, String> {
        self.skip();
        let start = self.pos;
        if self.b.get(self.pos) == Some(&b'"') {
            self.pos += 1;
            let s = self.pos;
            while self.pos < self.b.len() && self.b[self.pos] != b'"' {
                self.pos += 1;
            }
            if self.pos >= self.b.len() {
                return Err(self.err("a quoted word without its closing quote"));
            }
            let w = std::str::from_utf8(&self.b[s..self.pos]).map_err(|_| self.err("bad text"))?;
            self.pos += 1;
            return Ok(w);
        }
        while self.pos < self.b.len()
            && !self.b[self.pos].is_ascii_whitespace()
            && !b"(){};".contains(&self.b[self.pos])
        {
            self.pos += 1;
        }
        if start == self.pos {
            return Err(self.err("unexpected end of file or punctuation"));
        }
        std::str::from_utf8(&self.b[start..self.pos]).map_err(|_| self.err("bad text"))
    }

    fn number<T: std::str::FromStr>(&mut self) -> Result<T, String> {
        let w = self.word()?;
        w.parse()
            .map_err(|_| self.err(format!("cannot read a number from '{w}'")))
    }

    fn header(&mut self) -> Result<(), String> {
        let w = self.word()?;
        if w != "FoamFile" {
            return Err(self.err(format!("expected FoamFile, found '{w}'")));
        }
        self.expect(b'{')?;
        while !self.eat(b'}') {
            let key = self.word()?;
            let mut value = String::new();
            while !self.eat(b';') {
                if !value.is_empty() {
                    value.push(' ');
                }
                value += self.word()?;
            }
            match key {
                "format" => match value.as_str() {
                    "ascii" => {}
                    "binary" => self.head.binary = true,
                    v => return Err(self.err(format!("format {v}"))),
                },
                "class" => self.head.class = value,
                "arch" => {
                    for part in value.split(';') {
                        match part.trim() {
                            "label=32" => self.head.label_bytes = 4,
                            "label=64" => self.head.label_bytes = 8,
                            "scalar=32" => self.head.scalar_bytes = 4,
                            "scalar=64" => self.head.scalar_bytes = 8,
                            "MSB" => self.head.big_endian = true,
                            _ => {}
                        }
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// A count, then the items inside `( )` or `{ }` after it. Returns the
    /// count and leaves the position behind the `(`; a uniform list
    /// `N{value}` is returned as `Err(Some(value))`.
    fn open_list(&mut self) -> Result<usize, String> {
        let n: usize = self.number()?;
        let n = declared(n, self.b.len(), &self.err("a list"))?;
        self.expect(b'(')?;
        Ok(n)
    }

    fn raw(&mut self, n: usize) -> Result<&'a [u8], String> {
        if self.pos + n > self.b.len() {
            return Err(self.err("unexpected end of binary data"));
        }
        let s = &self.b[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    fn label_from(&self, c: &[u8]) -> u64 {
        match (self.head.label_bytes, self.head.big_endian) {
            (4, false) => u32::from_le_bytes(c.try_into().unwrap()) as u64,
            (4, true) => u32::from_be_bytes(c.try_into().unwrap()) as u64,
            (_, false) => u64::from_le_bytes(c.try_into().unwrap()),
            (_, true) => u64::from_be_bytes(c.try_into().unwrap()),
        }
    }

    fn scalar_from(&self, c: &[u8]) -> f64 {
        match (self.head.scalar_bytes, self.head.big_endian) {
            (4, false) => f32::from_le_bytes(c.try_into().unwrap()) as f64,
            (4, true) => f32::from_be_bytes(c.try_into().unwrap()) as f64,
            (_, false) => f64::from_le_bytes(c.try_into().unwrap()),
            (_, true) => f64::from_be_bytes(c.try_into().unwrap()),
        }
    }

    /// A list of labels: `N(a b c)`, `N{a}` (ASCII), or `N(` raw `)`.
    fn labels(&mut self) -> Result<Vec<u32>, String> {
        let n: usize = self.number()?;
        let n = declared(n, self.b.len(), &self.err("a list of labels"))?;
        if self.eat(b'{') {
            let v: u32 = self.number()?;
            self.expect(b'}')?;
            return Ok(vec![v; n]);
        }
        self.expect(b'(')?;
        let mut out = Vec::with_capacity(n);
        if self.head.binary {
            let lb = self.head.label_bytes;
            for c in self.raw(n * lb)?.chunks_exact(lb) {
                out.push(to_u32(self.label_from(c), self)?);
            }
        } else {
            for _ in 0..n {
                out.push(self.number()?);
            }
        }
        self.expect(b')')?;
        Ok(out)
    }

    fn points(&mut self) -> Result<Vec<[f64; 3]>, String> {
        let n = self.open_list()?;
        let mut out = Vec::with_capacity(n);
        if self.head.binary {
            let sb = self.head.scalar_bytes;
            let data = self.raw(n * 3 * sb)?;
            for p in data.chunks_exact(3 * sb) {
                let x = |i: usize| self.scalar_from(&p[i * sb..(i + 1) * sb]);
                out.push([x(0), x(1), x(2)]);
            }
        } else {
            for _ in 0..n {
                self.expect(b'(')?;
                out.push([self.number()?, self.number()?, self.number()?]);
                self.expect(b')')?;
            }
        }
        self.expect(b')')?;
        Ok(out)
    }

    /// Faces as (offsets, nodes).
    fn faces(&mut self) -> Result<(Vec<u32>, Vec<u32>), String> {
        if self.head.class == "faceCompactList" {
            let offsets = self.labels()?;
            let nodes = self.labels()?;
            if offsets.last().map(|&o| o as usize) != Some(nodes.len()) {
                return Err(self.err("the last offset of the compact face list is not its length"));
            }
            return Ok((offsets, nodes));
        }
        if self.head.binary {
            return Err(self.err(
                "a binary faceList: only the compact form (faceCompactList) is read; \
                 write the mesh again, or convert it to ASCII (foamFormatConvert)",
            ));
        }
        let n = self.open_list()?;
        let mut offsets = Vec::with_capacity(n + 1);
        let mut nodes = Vec::new();
        offsets.push(0);
        for _ in 0..n {
            let k: usize = self.number()?;
            self.expect(b'(')?;
            for _ in 0..k {
                nodes.push(self.number()?);
            }
            self.expect(b')')?;
            offsets.push(nodes.len() as u32);
        }
        self.expect(b')')?;
        Ok((offsets, nodes))
    }
}

fn to_u32(v: u64, s: &Scan) -> Result<u32, String> {
    u32::try_from(v).map_err(|_| s.err(format!("label {v} does not fit 32 bits")))
}

struct Patch {
    name: String,
    kind: String,
    n_faces: usize,
    start: usize,
}

fn patches(s: &mut Scan) -> Result<Vec<Patch>, String> {
    let n = s.open_list()?;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let name = s.word()?.to_string();
        s.expect(b'{')?;
        let (mut kind, mut n_faces, mut start) = (String::new(), None, None);
        while !s.eat(b'}') {
            let key = s.word()?;
            if s.peek() == Some(b'{') {
                skip_balanced(s)?;
                continue;
            }
            let mut value: Vec<&str> = Vec::new();
            let mut depth = 0usize;
            loop {
                match s.peek() {
                    Some(b';') if depth == 0 => {
                        s.pos += 1;
                        break;
                    }
                    Some(b'(') => {
                        depth += 1;
                        s.pos += 1;
                    }
                    Some(b')') => {
                        depth = depth.saturating_sub(1);
                        s.pos += 1;
                    }
                    Some(_) => value.push(s.word()?),
                    None => return Err(s.err("unexpected end of file in the boundary")),
                }
            }
            let first = value.first().copied().unwrap_or("");
            match key {
                "type" => kind = first.to_string(),
                "nFaces" => n_faces = first.parse::<usize>().ok(),
                "startFace" => start = first.parse::<usize>().ok(),
                _ => {}
            }
        }
        let (Some(n_faces), Some(start)) = (n_faces, start) else {
            return Err(s.err(format!("patch '{name}' needs nFaces and startFace")));
        };
        out.push(Patch {
            name,
            kind,
            n_faces,
            start,
        });
    }
    Ok(out)
}

fn skip_balanced(s: &mut Scan) -> Result<(), String> {
    s.expect(b'{')?;
    let mut depth = 1;
    while depth > 0 {
        match s.peek() {
            Some(b'{') => depth += 1,
            Some(b'}') => depth -= 1,
            Some(_) => {}
            None => return Err(s.err("unbalanced braces")),
        }
        s.pos += 1;
    }
    Ok(())
}

/// The cell zones: (name, cells).
fn cell_zones(s: &mut Scan) -> Result<Vec<(String, Vec<u32>)>, String> {
    let n = s.open_list()?;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let name = s.word()?.to_string();
        s.expect(b'{')?;
        let mut cells = Vec::new();
        while !s.eat(b'}') {
            let key = s.word()?;
            if key == "cellLabels" {
                // `List<label> N(...)`, or `List<label> 0()`
                if s.peek() != Some(b'{') {
                    let _ = s.word()?; // List<label>
                }
                cells = s.labels()?;
                s.expect(b';')?;
            } else if s.peek() == Some(b'{') {
                skip_balanced(s)?;
            } else {
                while !s.eat(b';') {
                    let _ = s.word()?;
                }
            }
        }
        out.push((name, cells));
    }
    Ok(out)
}

fn read_file(dir: &Path, name: &str) -> Result<Vec<u8>, String> {
    let p = dir.join(name);
    if !p.is_file() && dir.join(format!("{name}.gz")).is_file() {
        return Err(format!(
            "{}.gz is compressed: unzip the polyMesh files (gunzip {}/*.gz)",
            p.display(),
            dir.display()
        ));
    }
    std::fs::read(&p).map_err(|e| format!("{}: {e}", p.display()))
}

/// Reads the polyMesh files of `dir`.
pub fn read_dir(dir: &Path) -> Result<MeshFile, String> {
    let points_b = read_file(dir, "points")?;
    let faces_b = read_file(dir, "faces")?;
    let owner_b = read_file(dir, "owner")?;
    let neighbour_b = read_file(dir, "neighbour")?;
    let boundary_b = read_file(dir, "boundary")?;
    let zones_b = std::fs::read(dir.join("cellZones")).ok();

    let mut s = Scan::new(&points_b, "points")?;
    let binary = s.head.binary;
    let arch = format!(
        "{}, label={}, scalar={}",
        if binary { "binary" } else { "ASCII" },
        8 * s.head.label_bytes,
        8 * s.head.scalar_bytes
    );
    let points = s.points()?;
    let mut s = Scan::new(&faces_b, "faces")?;
    let (offsets, nodes) = s.faces()?;
    let mut s = Scan::new(&owner_b, "owner")?;
    let owner = s.labels()?;
    let mut s = Scan::new(&neighbour_b, "neighbour")?;
    let neighbour = s.labels()?;
    let mut s = Scan::new(&boundary_b, "boundary")?;
    let patches = patches(&mut s)?;

    let Some(n_faces) = offsets.len().checked_sub(1) else {
        return Err("faces: an empty offset list".into());
    };
    if offsets[0] != 0 || offsets.windows(2).any(|w| w[0] > w[1]) {
        return Err("faces: the offsets of the compact face list do not rise from 0".into());
    }
    if offsets[n_faces] as usize > nodes.len() {
        return Err("faces: an offset beyond the node list".into());
    }
    let n_inner = neighbour.len();
    if owner.len() != n_faces {
        return Err(format!(
            "owner has {} entries for {n_faces} faces",
            owner.len()
        ));
    }
    if n_inner > n_faces {
        return Err(format!(
            "neighbour has {n_inner} entries for {n_faces} faces"
        ));
    }
    if let Some(&n) = nodes.iter().find(|&&n| n as usize >= points.len()) {
        return Err(format!(
            "a face refers to point {n}, beyond the {} points",
            points.len()
        ));
    }
    let n_cells = owner
        .iter()
        .chain(&neighbour)
        .copied()
        .max()
        .map_or(0, |m| m as usize + 1);
    if n_cells == 0 {
        return Err("the mesh has no cells".into());
    }
    // every cell owns a face or is the neighbour of one
    if n_cells > n_faces + n_inner {
        return Err(format!(
            "owner or neighbour refers to cell {}, but {n_faces} faces bound at most {} cells",
            n_cells - 1,
            n_faces + n_inner
        ));
    }

    // patches must tile the faces behind the internal ones
    let mut next = n_inner;
    for p in &patches {
        if p.kind == "processor" || p.kind == "processorCyclic" {
            return Err(format!(
                "patch '{}' is a processor patch: the case is decomposed; \
                 reconstruct it first (reconstructParMesh)",
                p.name
            ));
        }
        if p.start != next {
            return Err(format!(
                "patch '{}' starts at face {}, expected {next}: the patches must follow the \
                 internal faces in order",
                p.name, p.start
            ));
        }
        next = next
            .checked_add(p.n_faces)
            .filter(|&n| n <= n_faces)
            .ok_or_else(|| {
                format!(
                    "patch '{}' with {} faces reaches beyond the {n_faces} faces",
                    p.name, p.n_faces
                )
            })?;
    }
    if next != n_faces {
        return Err(format!(
            "the patches cover {} of the {} boundary faces",
            next - n_inner,
            n_faces - n_inner
        ));
    }

    let mut b = MeshBuilder::new(points, n_cells);
    for p in &patches {
        b.declare_boundary(&p.name);
    }
    let face = |f: usize| &nodes[offsets[f] as usize..offsets[f + 1] as usize];
    for f in 0..n_inner {
        if owner[f] == neighbour[f] {
            return Err(format!("face {f} has the cell {} on both sides", owner[f]));
        }
        b.add_interior(face(f), owner[f], neighbour[f]);
    }
    #[allow(clippy::needless_range_loop)] // `f` also selects the face
    for p in &patches {
        for f in p.start..p.start + p.n_faces {
            b.add_boundary(face(f), owner[f], &p.name);
        }
    }

    // cell zones
    let mut warnings = Vec::new();
    if let Some(zb) = zones_b {
        let mut s = Scan::new(&zb, "cellZones")?;
        let mut zones = cell_zones(&mut s)?;
        zones.retain(|z| !z.1.is_empty());
        if !zones.is_empty() {
            let mut names: Vec<String> = zones.iter().map(|z| z.0.clone()).collect();
            let unset = u32::MAX;
            let mut cz = vec![unset; n_cells];
            for (zi, (name, cells)) in zones.iter().enumerate() {
                for &c in cells {
                    let c = c as usize;
                    if c >= n_cells {
                        return Err(format!(
                            "cell zone '{name}' lists cell {c}, beyond {n_cells}"
                        ));
                    }
                    cz[c] = zi as u32;
                }
            }
            let missing = cz.iter().filter(|&&z| z == unset).count();
            if missing > 0 {
                names.push("default".into());
                let d = (names.len() - 1) as u32;
                for z in &mut cz {
                    if *z == unset {
                        *z = d;
                    }
                }
                warnings.push(format!(
                    "{missing} cells are in no cell zone and form the zone 'default'"
                ));
            }
            b.set_zones(names, cz);
        }
    }

    let mut mesh = b.finish()?;
    mesh.boundaries.retain(|p| p.len > 0);
    let empties: Vec<&str> = patches
        .iter()
        .filter(|p| p.kind == "empty" && p.n_faces > 0)
        .map(|p| p.name.as_str())
        .collect();
    if !empties.is_empty() {
        warnings.push(format!(
            "patch {} is of type empty (a 2D OpenFOAM mesh): it is an ordinary boundary here; \
             give it the boundary condition `symmetry`",
            empties.join(", ")
        ));
    }
    Ok(MeshFile {
        mesh,
        reader: "openfoam",
        format: format!("OpenFOAM polyMesh {arch}"),
        warnings,
    })
}
