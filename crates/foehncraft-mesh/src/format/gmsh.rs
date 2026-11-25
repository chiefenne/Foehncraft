//! Reader for Gmsh `.msh` files, format version 4.1, ASCII and binary (D33).
//!
//! - Linear 3D elements become cells: tetrahedra, hexahedra, prisms,
//!   pyramids. Higher-order elements are refused with a message (mesh with
//!   `-order 1`).
//! - 2D elements (triangles, quadrilaterals) on the boundary name the
//!   boundary faces: the boundary is the name of the physical group of their
//!   surface (or `surface_<tag>` without one). Boundary faces not covered by
//!   any 2D element go to the boundary `default`, with a warning.
//! - Volume physical groups become cell zones (or `volume_<tag>`).
//! - `$Periodic` links are read and kept for the periodic boundary condition.
//! - Nodes not used by any cell are dropped.
//!
//! The format is described in the Gmsh reference manual, section "MSH file
//! format" (<https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format>).

use super::{CellKind, ElementMesh, MeshFile, MeshFormat, ReadOptions, declared};
use crate::mesh::Mesh;
use std::collections::HashMap;
use std::path::Path;

/// A periodic link between two entities.
#[derive(Debug, Clone, Default)]
pub struct PeriodicLink {
    /// Dimension of the entities.
    pub dim: i32,
    /// Tag of the slave entity.
    pub tag: i32,
    /// Tag of the master entity.
    pub master: i32,
    /// Affine transformation (16 values, row-major 4×4), if given.
    pub affine: Vec<f64>,
    /// Pairs (slave node, master node) as indices into the mesh nodes;
    /// pairs with a node that no cell uses are left out.
    pub nodes: Vec<(u32, u32)>,
}

/// Result of reading a Gmsh file.
#[derive(Debug, Clone)]
pub struct GmshMesh {
    /// The mesh.
    pub mesh: Mesh,
    /// Periodic links.
    pub periodic: Vec<PeriodicLink>,
    /// Format version and encoding, for reports.
    pub format: String,
    /// Warnings.
    pub warnings: Vec<String>,
}

/// Reads a `.msh` file.
pub fn read(path: impl AsRef<Path>) -> Result<GmshMesh, String> {
    let path = path.as_ref();
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&bytes).map_err(|e| format!("{}: {e}", path.display()))
}

struct Cursor<'a> {
    b: &'a [u8],
    pos: usize,
    binary: bool,
    size_t: usize,
    swap: bool,
}

impl<'a> Cursor<'a> {
    fn skip_ws(&mut self) {
        while self.pos < self.b.len() && self.b[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn token(&mut self) -> Result<&'a str, String> {
        self.skip_ws();
        let start = self.pos;
        while self.pos < self.b.len() && !self.b[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        if start == self.pos {
            return Err("unexpected end of file".into());
        }
        std::str::from_utf8(&self.b[start..self.pos]).map_err(|_| "invalid text".to_string())
    }

    fn line(&mut self) -> Result<&'a str, String> {
        self.skip_ws();
        let start = self.pos;
        while self.pos < self.b.len() && self.b[self.pos] != b'\n' {
            self.pos += 1;
        }
        let s = std::str::from_utf8(&self.b[start..self.pos]).map_err(|_| "invalid text")?;
        Ok(s.trim_end_matches('\r'))
    }

    fn take<const N: usize>(&mut self) -> Result<[u8; N], String> {
        if self.pos + N > self.b.len() {
            return Err("unexpected end of binary data".into());
        }
        let mut a = [0u8; N];
        a.copy_from_slice(&self.b[self.pos..self.pos + N]);
        if self.swap {
            a.reverse();
        }
        self.pos += N;
        Ok(a)
    }

    fn ascii<T: std::str::FromStr>(&mut self) -> Result<T, String> {
        let t = self.token()?;
        t.parse()
            .map_err(|_| format!("cannot read number from '{t}'"))
    }

    fn int(&mut self) -> Result<i64, String> {
        if self.binary {
            Ok(i32::from_le_bytes(self.take::<4>()?) as i64)
        } else {
            self.ascii()
        }
    }

    fn size(&mut self) -> Result<u64, String> {
        if self.binary {
            Ok(match self.size_t {
                4 => u32::from_le_bytes(self.take::<4>()?) as u64,
                _ => u64::from_le_bytes(self.take::<8>()?),
            })
        } else {
            self.ascii()
        }
    }

    fn double(&mut self) -> Result<f64, String> {
        if self.binary {
            Ok(f64::from_le_bytes(self.take::<8>()?))
        } else {
            self.ascii()
        }
    }

    /// Moves past `$End<name>`.
    fn end(&mut self, name: &str) -> Result<(), String> {
        let tag = format!("$End{name}");
        let hay = &self.b[self.pos..];
        let at = hay
            .windows(tag.len())
            .position(|w| w == tag.as_bytes())
            .ok_or_else(|| format!("missing {tag}"))?;
        self.pos += at + tag.len();
        Ok(())
    }

    fn expect_end(&mut self, name: &str) -> Result<(), String> {
        let t = self.token()?;
        if t != format!("$End{name}") {
            return Err(format!("expected $End{name}, found '{t}'"));
        }
        Ok(())
    }
}

/// Element type: (number of nodes, dimension, order), or the reason it is
/// refused. Second-order types are read by [`parse_raw`] (the structural
/// solver's quadratic elements); [`parse`] refuses them.
fn element_info(t: i64) -> Result<(usize, u8, u8), String> {
    Ok(match t {
        15 => (1, 0, 1),
        1 => (2, 1, 1),
        2 => (3, 2, 1),
        3 => (4, 2, 1),
        4 => (4, 3, 1),
        5 => (8, 3, 1),
        6 => (6, 3, 1),
        7 => (5, 3, 1),
        8 => (3, 1, 2),
        9 => (6, 2, 2),
        10 => (9, 2, 2),
        11 => (10, 3, 2),
        12 => (27, 3, 2),
        13 => (18, 3, 2),
        14 => (14, 3, 2),
        16 => (8, 2, 2),
        17 => (20, 3, 2),
        18 => (15, 3, 2),
        19 => (13, 3, 2),
        20..=25 | 26..=31 | 36..=99 => {
            return Err(format!(
                "element type {t} is a higher-order element; mesh with order 1 or 2 \
                 (Gmsh option -order or Mesh.ElementOrder)"
            ));
        }
        _ => return Err(format!("unsupported element type {t}")),
    })
}

/// The message for a second-order element where only linear ones are read.
fn refuse_order(t: i64) -> String {
    format!(
        "element type {t} is a higher-order element; mesh with order 1 \
         (Gmsh option -order 1 or Mesh.ElementOrder = 1)"
    )
}

/// A block of elements of one type on one entity, as read.
#[derive(Debug, Clone)]
pub struct RawBlock {
    /// Dimension of the entity (2 or 3: lower dimensions are skipped).
    pub dim: u8,
    /// Entity tag.
    pub entity: i64,
    /// Gmsh element type (4 tetrahedron, 11 ten-node tetrahedron, …).
    pub etype: i64,
    /// Nodes per element.
    pub nodes_per: usize,
    /// Node tags, `nodes_per` per element in Gmsh order.
    pub nodes: Vec<u64>,
}

/// A `.msh` file as read, before any mesh is built: nodes, element blocks
/// of dimension 2 and 3 (any order up to 2), physical names and periodic
/// links.
#[derive(Debug, Clone, Default)]
pub struct RawGmsh {
    /// Format version and encoding, for reports.
    pub format: String,
    /// Node tags.
    pub node_tags: Vec<u64>,
    /// Node coordinates, in the order of `node_tags`.
    pub nodes: Vec<[f64; 3]>,
    /// Element blocks of dimension 2 and 3.
    pub blocks: Vec<RawBlock>,
    /// Physical names ((dimension, tag), name), in file order.
    pub physical_names: Vec<((i64, i64), String)>,
    /// Physical tags of each entity, by (dimension, entity tag).
    pub entity_physicals: HashMap<(u8, i64), Vec<i64>>,
    /// Periodic links (dimension, tag, master tag, affine, node tag pairs).
    pub periodic: Vec<RawPeriodic>,
    /// Warnings.
    pub warnings: Vec<String>,
}

/// A periodic link as read: (dimension, tag, master tag, affine transform,
/// node tag pairs).
pub type RawPeriodic = (i64, i64, i64, Vec<f64>, Vec<(u64, u64)>);

impl RawGmsh {
    /// The index of each node tag: a closure that maps a tag to its index
    /// in `nodes`, or an error for an unknown tag.
    pub fn indexer(&self) -> impl Fn(u64) -> Result<u32, String> + '_ {
        // a table for dense tags, a map for sparse ones: the same indices
        let max_tag = self.node_tags.iter().copied().max().unwrap_or(0) as usize;
        let dense = max_tag <= 4 * self.node_tags.len() + 1024;
        let mut tag_index = vec![u32::MAX; if dense { max_tag + 1 } else { 0 }];
        let mut tag_map: HashMap<u64, u32> = HashMap::new();
        for (i, &t) in self.node_tags.iter().enumerate() {
            if dense {
                tag_index[t as usize] = i as u32;
            } else {
                tag_map.insert(t, i as u32);
            }
        }
        move |t: u64| -> Result<u32, String> {
            if dense {
                tag_index.get(t as usize).copied()
            } else {
                tag_map.get(&t).copied()
            }
            .filter(|&i| i != u32::MAX)
            .ok_or_else(|| format!("element refers to unknown node {t}"))
        }
    }

    /// The name of a physical group, or `physical_<tag>`.
    pub fn physical_name(&self, dim: i64, tag: i64) -> String {
        self.physical_names
            .iter()
            .find(|(k, _)| *k == (dim, tag))
            .map(|(_, n)| n.clone())
            .unwrap_or_else(|| format!("physical_{tag}"))
    }

    /// The name of an entity: its first physical group, or
    /// `surface_<tag>` / `volume_<tag>`.
    pub fn entity_name(&self, dim: u8, entity: i64) -> String {
        match self.entity_physicals.get(&(dim, entity)).and_then(|v| v.first()) {
            Some(&p) => self.physical_name(dim as i64, p),
            None if dim == 2 => format!("surface_{entity}"),
            None => format!("volume_{entity}"),
        }
    }

    /// The names of every physical group of an entity (its own name if it
    /// has none).
    pub fn entity_names(&self, dim: u8, entity: i64) -> Vec<String> {
        match self.entity_physicals.get(&(dim, entity)) {
            Some(v) if !v.is_empty() => v
                .iter()
                .map(|&p| self.physical_name(dim as i64, p))
                .collect(),
            _ => vec![self.entity_name(dim, entity)],
        }
    }
}

/// Reads a `.msh` file without building a mesh: second-order elements are
/// kept (see [`RawGmsh`]).
pub fn read_raw(path: impl AsRef<Path>) -> Result<RawGmsh, String> {
    let path = path.as_ref();
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse_raw(&bytes).map_err(|e| format!("{}: {e}", path.display()))
}

/// Parses the bytes of a `.msh` file without building a mesh.
pub fn parse_raw(bytes: &[u8]) -> Result<RawGmsh, String> {
    let mut c = Cursor {
        b: bytes,
        pos: 0,
        binary: false,
        size_t: 8,
        swap: false,
    };
    let mut format = String::new();
    let mut phys_names: Vec<((i64, i64), String)> = Vec::new();
    let mut entity_phys: HashMap<(u8, i64), Vec<i64>> = HashMap::new();
    let mut node_tag: Vec<u64> = Vec::new();
    let mut node_xyz: Vec<[f64; 3]> = Vec::new();
    let mut blocks: Vec<RawBlock> = Vec::new();
    let mut periodic: Vec<RawPeriodic> = Vec::new();
    let mut warnings = Vec::new();

    loop {
        c.skip_ws();
        if c.pos >= bytes.len() {
            break;
        }
        let head = c.token()?;
        let Some(name) = head.strip_prefix('$') else {
            return Err(format!("expected a section, found '{head}'"));
        };
        let name = name.to_string();
        if c.binary {
            // binary data start right after the newline of the section line
            if c.b.get(c.pos) == Some(&b'\r') {
                c.pos += 1;
            }
            if c.b.get(c.pos) == Some(&b'\n') {
                c.pos += 1;
            }
        }
        match name.as_str() {
            "MeshFormat" => {
                let version: String = c.token()?.to_string();
                let ft: i64 = c.ascii()?;
                let ds: usize = c.ascii()?;
                if !version.starts_with("4.1") {
                    return Err(format!(
                        "format version {version}: only 4.1 is read (save with \
                         Mesh.MshFileVersion = 4.1)"
                    ));
                }
                c.binary = ft == 1;
                c.size_t = ds;
                if c.binary {
                    // one newline, then the integer 1 in the writer's byte order
                    while c.pos < c.b.len() && c.b[c.pos] != b'\n' {
                        c.pos += 1;
                    }
                    c.pos += 1;
                    let one = c.take::<4>()?;
                    if i32::from_le_bytes(one) != 1 {
                        if i32::from_be_bytes(one) == 1 {
                            c.swap = true;
                        } else {
                            return Err("cannot determine the byte order".into());
                        }
                    }
                }
                format = format!(
                    "MSH {version} {}",
                    if c.binary { "binary" } else { "ASCII" }
                );
                c.expect_end("MeshFormat")?;
            }
            "PhysicalNames" => {
                // always ASCII
                let n: usize = c.ascii()?;
                for _ in 0..n {
                    let dim: i64 = c.ascii()?;
                    let tag: i64 = c.ascii()?;
                    let rest = c.line()?;
                    let nm = rest.trim().trim_matches('"').to_string();
                    phys_names.retain(|(k, _)| *k != (dim, tag));
                    phys_names.push(((dim, tag), nm));
                }
                c.expect_end("PhysicalNames")?;
            }
            "Entities" => {
                let counts = [c.size()?, c.size()?, c.size()?, c.size()?];
                for (dim, &cnt) in counts.iter().enumerate() {
                    for _ in 0..cnt {
                        let tag = c.int()?;
                        let nbox = if dim == 0 { 3 } else { 6 };
                        for _ in 0..nbox {
                            c.double()?;
                        }
                        let np = c.size()?;
                        let mut tags = Vec::with_capacity(declared(
                            np as usize,
                            bytes.len(),
                            "physical tags of an entity",
                        )?);
                        for _ in 0..np {
                            tags.push(c.int()?);
                        }
                        entity_phys.insert((dim as u8, tag), tags);
                        if dim > 0 {
                            let nb = c.size()?;
                            for _ in 0..nb {
                                c.int()?;
                            }
                        }
                    }
                }
                c.expect_end("Entities")?;
            }
            "Nodes" => {
                let nblocks = c.size()?;
                let nnodes = declared(c.size()? as usize, bytes.len(), "nodes")?;
                c.size()?;
                c.size()?;
                node_tag.reserve(nnodes);
                node_xyz.reserve(nnodes);
                for _ in 0..nblocks {
                    let dim = c.int()?;
                    c.int()?;
                    let parametric = c.int()?;
                    let n = c.size()? as usize;
                    for _ in 0..n {
                        node_tag.push(c.size()?);
                    }
                    let extra = if parametric != 0 { dim as usize } else { 0 };
                    for _ in 0..n {
                        let p = [c.double()?, c.double()?, c.double()?];
                        for _ in 0..extra {
                            c.double()?;
                        }
                        node_xyz.push(p);
                    }
                }
                c.expect_end("Nodes")?;
            }
            "Elements" => {
                let nblocks = c.size()?;
                c.size()?;
                c.size()?;
                c.size()?;
                for _ in 0..nblocks {
                    let dim = c.int()? as u8;
                    let entity = c.int()?;
                    let etype = c.int()?;
                    let n = declared(c.size()? as usize, bytes.len(), "elements of a block")?;
                    let (nn, edim, _) = element_info(etype)?;
                    if edim != dim {
                        return Err(format!(
                            "an element block of dimension {dim} holds elements of type {etype} \
                             (dimension {edim})"
                        ));
                    }
                    let mut nodes = Vec::with_capacity(if dim >= 2 { n * nn } else { 0 });
                    for _ in 0..n {
                        c.size()?;
                        for _ in 0..nn {
                            let t = c.size()?;
                            if dim >= 2 {
                                nodes.push(t);
                            }
                        }
                    }
                    if dim >= 2 {
                        blocks.push(RawBlock {
                            dim,
                            entity,
                            etype,
                            nodes_per: nn,
                            nodes,
                        });
                    }
                }
                c.expect_end("Elements")?;
            }
            "Periodic" => {
                let n = c.size()?;
                for _ in 0..n {
                    let dim = c.int()?;
                    let tag = c.int()?;
                    let master = c.int()?;
                    let na = c.size()?;
                    let mut aff =
                        Vec::with_capacity(declared(na as usize, bytes.len(), "affine values")?);
                    for _ in 0..na {
                        aff.push(c.double()?);
                    }
                    let nn = c.size()?;
                    let mut pairs =
                        Vec::with_capacity(declared(nn as usize, bytes.len(), "periodic nodes")?);
                    for _ in 0..nn {
                        pairs.push((c.size()?, c.size()?));
                    }
                    periodic.push((dim, tag, master, aff, pairs));
                }
                c.expect_end("Periodic")?;
            }
            _ => {
                if name.starts_with("End") {
                    return Err(format!("unexpected ${name}"));
                }
                if name == "PartitionedEntities" {
                    warnings.push(
                        "the file is partitioned; partitions are ignored (Foehncraft \
                         partitions the mesh itself)"
                            .into(),
                    );
                }
                c.end(&name)?;
            }
        }
    }
    if format.is_empty() {
        return Err("no $MeshFormat section".into());
    }
    Ok(RawGmsh {
        format,
        node_tags: node_tag,
        nodes: node_xyz,
        blocks,
        physical_names: phys_names,
        entity_physicals: entity_phys,
        periodic,
        warnings,
    })
}

/// Parses the bytes of a `.msh` file.
pub fn parse(bytes: &[u8]) -> Result<GmshMesh, String> {
    let raw = parse_raw(bytes)?;
    if let Some(b) = raw
        .blocks
        .iter()
        .find(|b| element_info(b.etype).is_ok_and(|i| i.2 > 1))
    {
        return Err(refuse_order(b.etype));
    }
    let idx = raw.indexer();
    let blocks = &raw.blocks;
    if !blocks.iter().any(|b| b.dim == 3 && !b.nodes.is_empty()) {
        return Err(
            "the file has no 3D cells. For a 2D case, extrude the surface mesh one \
             layer in z (Extrude {0, 0, dz} { Surface{...}; Layers{1}; Recombine; })"
                .into(),
        );
    }

    let mut em = ElementMesh {
        points: raw.nodes.clone(),
        face_word: "surface elements",
        group_word: "surface physical group",
        ..Default::default()
    };
    // cells: volume physical groups are the zones
    let mut v = Vec::with_capacity(8);
    for b in blocks.iter().filter(|b| b.dim == 3) {
        let z = em.zone(&raw.entity_name(3, b.entity));
        let kind = match b.nodes_per {
            4 => CellKind::Tetrahedron,
            8 => CellKind::Hexahedron,
            6 => CellKind::Prism,
            5 => CellKind::Pyramid,
            _ => unreachable!(),
        };
        for e in b.nodes.chunks(b.nodes_per) {
            v.clear();
            for &t in e {
                v.push(idx(t)?);
            }
            em.add_cell(kind, &v, z);
        }
    }
    // boundary faces named by the 2D elements: boundaries in the order of
    // the physical names, then entity order
    for ((d, _), name) in &raw.physical_names {
        if *d == 2 {
            em.boundaries.push(name.clone());
        }
    }
    for b in blocks.iter().filter(|b| b.dim == 2) {
        let bi = em.boundary(&raw.entity_name(2, b.entity));
        for e in b.nodes.chunks(b.nodes_per) {
            v.clear();
            for &t in e {
                v.push(idx(t)?);
            }
            em.add_face(&v, bi);
        }
    }
    let built = em.into_mesh()?;
    let mut warnings = raw.warnings.clone();
    warnings.extend(built.warnings);
    let node = |t: u64| -> Result<u32, String> { Ok(built.node_map[idx(t)? as usize]) };

    let periodic = raw
        .periodic
        .iter()
        .cloned()
        .map(|(dim, tag, master, affine, pairs)| PeriodicLink {
            dim: dim as i32,
            tag: tag as i32,
            master: master as i32,
            affine,
            nodes: pairs
                .into_iter()
                .filter_map(|(a, b)| {
                    let (a, b) = (node(a).ok()?, node(b).ok()?);
                    (a != u32::MAX && b != u32::MAX).then_some((a, b))
                })
                .collect(),
        })
        .collect();

    Ok(GmshMesh {
        mesh: built.mesh,
        periodic,
        format: raw.format.clone(),
        warnings,
    })
}

/// The Gmsh format in the registry of [`crate::format`](mod@crate::format).
#[derive(Debug, Clone, Copy)]
pub struct Gmsh;

impl MeshFormat for Gmsh {
    fn name(&self) -> &'static str {
        "gmsh"
    }

    fn title(&self) -> &'static str {
        "Gmsh MSH 4.1, ASCII and binary"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &["msh"]
    }

    fn sniff(&self, head: &[u8]) -> bool {
        let start = head
            .iter()
            .position(|c| !c.is_ascii_whitespace())
            .unwrap_or(0);
        head[start..].starts_with(b"$MeshFormat")
    }

    fn parse(&self, bytes: &[u8], _options: &ReadOptions) -> Result<MeshFile, String> {
        let g = parse(bytes)?;
        Ok(MeshFile {
            mesh: g.mesh,
            reader: self.name(),
            format: g.format,
            warnings: g.warnings,
        })
    }
}
