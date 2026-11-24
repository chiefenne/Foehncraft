//! Mesh file formats (D33, D187): one reader interface for every format.
//!
//! A format is a [`MeshFormat`]: a name (`gmsh`, `su2`, `fluent_msh`, `vtk`,
//! `cgns`, `vtkhdf`, `openfoam`),
//! the file extensions it owns (an OpenFOAM mesh is a directory), a test of
//! the first bytes of a file, and a reader.
//! The readers are listed in [`FORMATS`]; [`read`] picks the reader from
//! the `mesh.format` key, else from the file extension, else from the
//! content, and returns a [`MeshFile`] with the same [`Mesh`] type for every
//! format.
//!
//! Most readers only parse: they hand their cells (tetrahedra, hexahedra,
//! prisms, pyramids, polyhedra) and their named boundary faces to an
//! [`ElementMesh`],
//! and [`ElementMesh::into_mesh`] builds the face-based mesh the same way
//! for every format (unused nodes dropped, faces matched, boundary faces
//! named, faces without a name in the boundary `default`). A new format is
//! a parser that fills an [`ElementMesh`], one [`MeshFormat`] value and one
//! line in [`FORMATS`].
//!
//! Boundary names come from the file (Gmsh physical groups, SU2 markers);
//! the case maps them to the names of its boundary conditions with
//! `[mesh.boundaries]` (`old = "new"`), whatever the format.

#[cfg(feature = "hdf5")]
pub mod cgns;
pub mod fluent_msh;
pub mod gmsh;
mod inflate;
pub mod openfoam;
pub mod su2;
pub mod vtk;
#[cfg(feature = "hdf5")]
pub mod vtkhdf;

use crate::builder::CellSoup;
use crate::generate::{hex_faces, tet_faces};
use crate::mesh::Mesh;
use std::collections::HashMap;
use std::path::Path;

/// A mesh read from a file.
#[derive(Debug, Clone)]
pub struct MeshFile {
    /// The mesh.
    pub mesh: Mesh,
    /// Name of the format (the `mesh.format` value, e.g. `gmsh`).
    pub reader: &'static str,
    /// Version, encoding and dimension as read, for reports
    /// (e.g. `MSH 4.1 ASCII`, `SU2 2D ASCII`).
    pub format: String,
    /// Warnings.
    pub warnings: Vec<String>,
}

/// Options of a reader that the file does not give.
#[derive(Debug, Clone, Copy, Default)]
pub struct ReadOptions {
    /// Cell layer thickness in z when a 2D mesh is extruded into one layer
    /// (in the units of the file; `None`: the median edge length of the 2D
    /// cells).
    pub thickness: Option<f64>,
    /// Join nodes that coincide (within 1e-8 of the bounding-box diagonal)
    /// after reading, so that the zones of a multi-zone file that touch share
    /// their interface nodes (SU2 and CGNS with several zones).
    pub merge_nodes: bool,
}

/// A mesh file format. A format that is a file implements [`Self::parse`]
/// (the default [`Self::read_path`] reads the file and calls it); one that
/// is a directory (OpenFOAM) or must be opened by a library (CGNS, HDF5)
/// implements [`Self::read_path`].
pub trait MeshFormat: Sync {
    /// Name, as given in `mesh.format`.
    fn name(&self) -> &'static str;
    /// One line for lists and docs.
    fn title(&self) -> &'static str;
    /// File extensions (lower case, without the dot).
    fn extensions(&self) -> &'static [&'static str];
    /// Whether the first bytes of a file look like this format.
    fn sniff(&self, head: &[u8]) -> bool;
    /// Reads a file from its bytes.
    fn parse(&self, _bytes: &[u8], _options: &ReadOptions) -> Result<MeshFile, String> {
        Err(format!("{} is not read from memory", self.name()))
    }
    /// Whether a directory holds a mesh of this format.
    fn accepts_dir(&self, _dir: &Path) -> bool {
        false
    }
    /// Reads the file or directory `path`.
    fn read_path(&self, path: &Path, options: &ReadOptions) -> Result<MeshFile, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        self.parse(&bytes, options)
    }
}

/// The formats Foehncraft reads.
pub static FORMATS: &[&dyn MeshFormat] = &[
    &gmsh::Gmsh,
    &su2::Su2,
    &fluent_msh::IndustryMsh,
    &vtk::Vtk,
    #[cfg(feature = "hdf5")]
    &cgns::Cgns,
    #[cfg(feature = "hdf5")]
    &vtkhdf::VtkHdf,
    &openfoam::OpenFoam,
];

/// A count, size or index that a file declares, checked before anything
/// is allocated for it (D202): every item a file declares takes at least
/// one byte of it, so a count above the size of the file (`len` bytes) is a
/// broken or hostile file, refused with a message instead of an allocation
/// that aborts the program.
pub(crate) fn declared(n: usize, len: usize, what: &str) -> Result<usize, String> {
    if n > len {
        Err(format!(
            "{what}: {n} declared, more than the file of {len} bytes can hold (a damaged file?)"
        ))
    } else {
        Ok(n)
    }
}

/// The format of a name.
pub fn by_name(name: &str) -> Result<&'static dyn MeshFormat, String> {
    // `industry_msh` is the former name of `fluent_msh`, still accepted
    let name = if name == "industry_msh" {
        "fluent_msh"
    } else {
        name
    };
    FORMATS
        .iter()
        .copied()
        .find(|f| f.name() == name)
        .ok_or_else(|| format!("unknown mesh format '{name}' (formats: {})", names()))
}

/// The format that owns the extension of `path`, if any.
pub fn by_extension(path: &Path) -> Option<&'static dyn MeshFormat> {
    let ext = path.extension()?.to_str()?.to_ascii_lowercase();
    FORMATS
        .iter()
        .copied()
        .find(|f| f.extensions().contains(&ext.as_str()))
}

/// The formats that own the extension of `path` (`.msh` is shared by Gmsh
/// and the Fluent format: the first bytes tell them apart).
fn owners(path: &Path) -> Vec<&'static dyn MeshFormat> {
    let Some(ext) = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
    else {
        return Vec::new();
    };
    FORMATS
        .iter()
        .copied()
        .filter(|f| f.extensions().contains(&ext.as_str()))
        .collect()
}

/// Whether `path` is a mesh of a known format: a file with the extension of
/// a format, or a directory that holds a mesh (OpenFOAM).
pub fn is_mesh_file(path: &Path) -> bool {
    if path.is_dir() {
        return FORMATS.iter().any(|f| f.accepts_dir(path));
    }
    by_extension(path).is_some()
}

fn names() -> String {
    FORMATS
        .iter()
        .map(|f| f.name())
        .collect::<Vec<_>>()
        .join(", ")
}

/// The format of a file: by `format` if given, else by extension, else by
/// content.
pub fn detect(
    path: &Path,
    bytes: &[u8],
    format: Option<&str>,
) -> Result<&'static dyn MeshFormat, String> {
    if let Some(name) = format {
        return by_name(name);
    }
    let own = owners(path);
    if own.len() > 1 {
        let head = &bytes[..bytes.len().min(4096)];
        if let Some(f) = own.iter().copied().find(|f| f.sniff(head)) {
            return Ok(f);
        }
    }
    if let Some(f) = own.first() {
        return Ok(*f);
    }
    if !cfg!(feature = "hdf5")
        && let Some(ext) = path.extension().and_then(|e| e.to_str())
        && (ext.eq_ignore_ascii_case("cgns") || ext.eq_ignore_ascii_case("vtkhdf"))
    {
        return Err(format!(
            "{} files are read through HDF5: build with the feature `hdf5`",
            if ext.eq_ignore_ascii_case("cgns") {
                "CGNS"
            } else {
                "VTKHDF"
            }
        ));
    }
    let head = &bytes[..bytes.len().min(4096)];
    if bytes.is_empty() {
        return Err("cannot tell the mesh format from the name; set `mesh.format`".into());
    }
    FORMATS
        .iter()
        .copied()
        .find(|f| f.sniff(head))
        .ok_or_else(|| {
            format!(
                "cannot tell the mesh format from the name or the content; set `mesh.format` \
                 (formats: {})",
                names()
            )
        })
}

/// Reads a mesh file with the default options; `format` forces the format.
pub fn read(path: impl AsRef<Path>, format: Option<&str>) -> Result<MeshFile, String> {
    read_with(path, format, &ReadOptions::default())
}

/// Reads a mesh file.
pub fn read_with(
    path: impl AsRef<Path>,
    format: Option<&str>,
    options: &ReadOptions,
) -> Result<MeshFile, String> {
    let path = path.as_ref();
    let err = |e: String| format!("{}: {e}", path.display());
    if path.is_dir() {
        let f = match format {
            Some(name) => by_name(name),
            None => FORMATS
                .iter()
                .copied()
                .find(|f| f.accepts_dir(path))
                .ok_or_else(|| {
                    "a directory that holds no mesh (an OpenFOAM mesh is a case directory with \
                     constant/polyMesh, or the polyMesh directory itself)"
                        .to_string()
                }),
        }
        .map_err(err)?;
        return f.read_path(path, options).map_err(err);
    }
    // sniff the first bytes only when the name does not decide
    let head = if format.is_some() || owners(path).len() == 1 {
        Vec::new()
    } else {
        use std::io::Read;
        let mut v = Vec::new();
        std::fs::File::open(path)
            .and_then(|f| f.take(4096).read_to_end(&mut v))
            .map_err(|e| err(e.to_string()))?;
        v
    };
    let f = detect(path, &head, format).map_err(err)?;
    f.read_path(path, options).map_err(err)
}

/// A linear 3D cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellKind {
    /// 4 nodes.
    Tetrahedron,
    /// 8 nodes: one quadrilateral (0-3), then the opposite one (4-7) with
    /// node 4 across from node 0.
    Hexahedron,
    /// 6 nodes: one triangle (0-2), then the opposite one (3-5) with node 3
    /// across from node 0.
    Prism,
    /// 5 nodes: the quadrilateral base (0-3), then the apex.
    Pyramid,
    /// Any closed set of polygons: the next polyhedron of
    /// [`ElementMesh::poly`].
    Polyhedron,
}

impl CellKind {
    /// Number of nodes.
    pub fn nodes(self) -> usize {
        match self {
            Self::Tetrahedron => 4,
            Self::Hexahedron => 8,
            Self::Prism => 6,
            Self::Pyramid => 5,
            Self::Polyhedron => 0,
        }
    }
}

/// Cells and named boundary faces, as most file formats store a mesh.
/// Orientation does not matter (the builder orients every cell from its
/// volume); only the node numbering of each kind ([`CellKind`]) does.
#[derive(Debug, Clone, Default)]
pub struct ElementMesh {
    /// Nodes.
    pub points: Vec<[f64; 3]>,
    /// Kind of each cell.
    pub cell_kind: Vec<CellKind>,
    /// Nodes of the cells, concatenated (none for a polyhedron).
    pub cell_nodes: Vec<u32>,
    /// The faces (node lists) of the polyhedra, one polyhedron per
    /// [`CellKind::Polyhedron`] cell in cell order.
    pub poly: Polyhedra,
    /// Zone of each cell (index into `zones`).
    pub cell_zone: Vec<u32>,
    /// Zone names (`fluid` if empty).
    pub zones: Vec<String>,
    /// Start of each boundary face in `face_nodes` (length faces + 1, or
    /// empty for none).
    pub face_ptr: Vec<u32>,
    /// Nodes of the boundary faces (triangles, quadrilaterals).
    pub face_nodes: Vec<u32>,
    /// Boundary of each face (index into `boundaries`). A face given twice
    /// takes the boundary of its last entry.
    pub face_boundary: Vec<u32>,
    /// Boundary names in the order they should appear; `default` is added.
    pub boundaries: Vec<String>,
    /// What the file calls its boundary faces in warnings
    /// (`surface elements`, `marker elements`).
    pub face_word: &'static str,
    /// What names a boundary in the file, in warnings
    /// (`surface physical group`, `marker`).
    pub group_word: &'static str,
}

/// The faces of the polyhedra of an [`ElementMesh`], kept as VTK face
/// streams in one flat array: for each polyhedron its number of faces,
/// then each face as its number of nodes and the nodes. A mesh of millions
/// of polyhedra is so a few large arrays of `u32` instead of one
/// allocation per face (a VTK reader moves the face stream of the file in
/// without a copy).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Polyhedra {
    stream: Vec<u32>,
    /// Start of each polyhedron in `stream`.
    start: Vec<usize>,
}

impl Polyhedra {
    /// Polyhedra whose face streams are in `stream`, starting at `start`;
    /// the caller has checked each with [`check_face_stream`].
    pub(crate) fn from_parts(stream: Vec<u32>, start: Vec<usize>) -> Self {
        Self { stream, start }
    }

    /// Adds a polyhedron from its faces.
    pub fn push<I>(&mut self, faces: I)
    where
        I: IntoIterator,
        I::Item: AsRef<[u32]>,
    {
        let at = self.stream.len();
        self.start.push(at);
        self.stream.push(0);
        let mut n = 0u32;
        for f in faces {
            let f = f.as_ref();
            self.stream.push(f.len() as u32);
            self.stream.extend_from_slice(f);
            n += 1;
        }
        self.stream[at] = n;
    }

    /// Number of polyhedra.
    pub fn len(&self) -> usize {
        self.start.len()
    }

    /// True without polyhedra.
    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    /// The faces of polyhedron `p` (none for an index out of range).
    pub fn faces(&self, p: usize) -> Faces<'_> {
        let at = self.start.get(p).copied().unwrap_or(self.stream.len());
        Faces {
            s: &self.stream,
            at: at + 1,
            left: self.stream.get(at).copied().unwrap_or(0),
        }
    }

    /// Calls `f` on every node of every face.
    fn nodes_mut(&mut self, mut f: impl FnMut(&mut u32)) {
        let Self { stream, start } = self;
        for &s in start.iter() {
            let nf = stream.get(s).copied().unwrap_or(0);
            let mut at = s + 1;
            for _ in 0..nf {
                let Some(&k) = stream.get(at) else { break };
                let end = (at + 1 + k as usize).min(stream.len());
                stream[at + 1..end].iter_mut().for_each(&mut f);
                at = end;
            }
        }
    }
}

/// The faces of one polyhedron of [`Polyhedra`].
#[derive(Debug, Clone)]
pub struct Faces<'a> {
    s: &'a [u32],
    at: usize,
    left: u32,
}

impl<'a> Iterator for Faces<'a> {
    type Item = &'a [u32];

    fn next(&mut self) -> Option<&'a [u32]> {
        if self.left == 0 {
            return None;
        }
        let k = *self.s.get(self.at)? as usize;
        let f = self.s.get(self.at + 1..self.at + 1 + k)?;
        self.at += 1 + k;
        self.left -= 1;
        Some(f)
    }
}

/// Checks the face stream of polyhedron `cell` (`n_faces, n0, nodes,
/// n1, nodes, ...`): at least one face, at least three nodes a face, and
/// no face beyond the end of `s`.
pub(crate) fn check_face_stream(s: &[u32], cell: usize) -> Result<(), String> {
    let nf = s.first().copied().unwrap_or(0);
    if nf == 0 {
        return Err(format!("polyhedron {cell} has no faces"));
    }
    let mut at = 1usize;
    for _ in 0..nf {
        let k = *s.get(at).ok_or("a polyhedron face stream ends too early")? as usize;
        if k < 3 {
            return Err(format!("polyhedron {cell} has a face of {k} points"));
        }
        at += 1 + k;
        if at > s.len() {
            return Err("a polyhedron face stream ends too early".into());
        }
    }
    Ok(())
}

/// A mesh built from an [`ElementMesh`].
#[derive(Debug, Clone)]
pub struct Built {
    /// The mesh.
    pub mesh: Mesh,
    /// Warnings.
    pub warnings: Vec<String>,
    /// New index of every input node (`u32::MAX`: used by no cell, dropped).
    pub node_map: Vec<u32>,
}

impl ElementMesh {
    /// Adds a cell.
    pub fn add_cell(&mut self, kind: CellKind, nodes: &[u32], zone: u32) {
        debug_assert_eq!(nodes.len(), kind.nodes());
        self.cell_kind.push(kind);
        self.cell_nodes.extend_from_slice(nodes);
        self.cell_zone.push(zone);
    }

    /// Adds a polyhedral cell from its faces (the orientation of the faces
    /// need not be consistent with the cell: it is made outward from the
    /// volume, but all faces of one cell must agree).
    pub fn add_polyhedron<I>(&mut self, faces: I, zone: u32)
    where
        I: IntoIterator,
        I::Item: AsRef<[u32]>,
    {
        self.cell_kind.push(CellKind::Polyhedron);
        self.poly.push(faces);
        self.cell_zone.push(zone);
    }

    /// Joins nodes that coincide, within 1e-8 of the bounding-box diagonal,
    /// into one and renumbers the cells and faces. Returns the number of
    /// nodes removed. Faces of two zones that touch then become interior
    /// faces of the one mesh.
    pub fn merge_points(&mut self) -> usize {
        let n = self.points.len();
        if n == 0 {
            return 0;
        }
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for p in &self.points {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let diag = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
        let tol = (1e-8 * diag).max(f64::MIN_POSITIVE);
        let cell = |p: &[f64; 3]| -> [i64; 3] {
            std::array::from_fn(|k| ((p[k] - lo[k]) / tol).floor() as i64)
        };
        let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
        let mut map = vec![0u32; n];
        let mut kept: Vec<[f64; 3]> = Vec::new();
        for (i, p) in self.points.iter().enumerate() {
            let c = cell(p);
            let mut found = None;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(list) = grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]) else {
                            continue;
                        };
                        for &j in list {
                            let q = kept[j as usize];
                            let d2: f64 = (0..3).map(|k| (p[k] - q[k]).powi(2)).sum();
                            if d2 <= tol * tol {
                                found = Some(j);
                                break 'search;
                            }
                        }
                    }
                }
            }
            map[i] = match found {
                Some(j) => j,
                None => {
                    kept.push(*p);
                    let j = (kept.len() - 1) as u32;
                    grid.entry(c).or_default().push(j);
                    j
                }
            };
        }
        let removed = n - kept.len();
        if removed == 0 {
            return 0;
        }
        self.points = kept;
        for v in &mut self.cell_nodes {
            *v = map[*v as usize];
        }
        self.poly.nodes_mut(|v| {
            if let Some(&m) = map.get(*v as usize) {
                *v = m;
            }
        });
        for v in &mut self.face_nodes {
            *v = map[*v as usize];
        }
        removed
    }

    /// Adds a boundary face.
    pub fn add_face(&mut self, nodes: &[u32], boundary: u32) {
        if self.face_ptr.is_empty() {
            self.face_ptr.push(0);
        }
        self.face_nodes.extend_from_slice(nodes);
        self.face_ptr.push(self.face_nodes.len() as u32);
        self.face_boundary.push(boundary);
    }

    /// Index of a boundary name, added if new.
    pub fn boundary(&mut self, name: &str) -> u32 {
        index_of(&mut self.boundaries, name)
    }

    /// Index of a zone name, added if new.
    pub fn zone(&mut self, name: &str) -> u32 {
        index_of(&mut self.zones, name)
    }

    /// Builds the face-based mesh: drops the nodes no cell uses (the others
    /// keep their order of first use by the cells), matches the faces,
    /// names the boundary faces and drops empty boundaries.
    pub fn into_mesh(self) -> Result<Built, String> {
        // the inputs are dropped as soon as they are used: a reader of a
        // large mesh peaks here
        let ElementMesh {
            points: input,
            cell_kind,
            cell_nodes,
            poly,
            cell_zone,
            zones,
            face_ptr,
            face_nodes,
            face_boundary,
            boundaries,
            face_word,
            group_word,
        } = self;
        let n_cells = cell_kind.len();
        if n_cells == 0 {
            return Err("the file has no 3D cells".into());
        }
        let np = input.len();
        let mut new_index = vec![u32::MAX; np];
        let mut points: Vec<[f64; 3]> = Vec::new();
        let mut use_node = |n: u32| -> Result<(), String> {
            let i = n as usize;
            if i >= np {
                return Err(format!("a cell refers to node {n}, beyond the {np} nodes"));
            }
            if new_index[i] == u32::MAX {
                new_index[i] = points.len() as u32;
                points.push(input[i]);
            }
            Ok(())
        };
        {
            let (mut at, mut ip) = (0usize, 0usize);
            for &kind in &cell_kind {
                if kind == CellKind::Polyhedron {
                    for &n in poly.faces(ip).flatten() {
                        use_node(n)?;
                    }
                    ip += 1;
                } else {
                    for &n in &cell_nodes[at..at + kind.nodes()] {
                        use_node(n)?;
                    }
                    at += kind.nodes();
                }
            }
        }
        drop(input);

        let mut soup = CellSoup::new(points);
        let (mut at, mut ip) = (0usize, 0usize);
        // the renumbered faces of one polyhedron
        let mut buf: Vec<u32> = Vec::new();
        let mut ends: Vec<usize> = Vec::new();
        for (c, &kind) in cell_kind.iter().enumerate() {
            if kind == CellKind::Polyhedron {
                buf.clear();
                ends.clear();
                for f in poly.faces(ip) {
                    buf.extend(f.iter().map(|&n| new_index[n as usize]));
                    ends.push(buf.len());
                }
                ip += 1;
                let mut s = 0;
                soup.add_cell_in_zone(
                    ends.iter().map(|&e| {
                        let f = &buf[s..e];
                        s = e;
                        f
                    }),
                    cell_zone[c],
                );
                continue;
            }
            let e = &cell_nodes[at..at + kind.nodes()];
            at += kind.nodes();
            let mut v = [0u32; 8];
            for (k, &n) in e.iter().enumerate() {
                v[k] = new_index[n as usize];
            }
            let z = cell_zone[c];
            match kind {
                CellKind::Tetrahedron => {
                    let f = tet_faces(&[v[0], v[1], v[2], v[3]]);
                    soup.add_cell_in_zone(f.iter().map(|x| &x[..]), z);
                }
                CellKind::Hexahedron => {
                    let f = hex_faces(&v);
                    soup.add_cell_in_zone(f.iter().map(|x| &x[..]), z);
                }
                CellKind::Prism => {
                    let f: [&[u32]; 5] = [
                        &[v[0], v[2], v[1]],
                        &[v[3], v[4], v[5]],
                        &[v[0], v[1], v[4], v[3]],
                        &[v[1], v[2], v[5], v[4]],
                        &[v[2], v[0], v[3], v[5]],
                    ];
                    soup.add_cell_in_zone(f, z);
                }
                CellKind::Pyramid => {
                    let f: [&[u32]; 5] = [
                        &[v[0], v[3], v[2], v[1]],
                        &[v[0], v[1], v[4]],
                        &[v[1], v[2], v[4]],
                        &[v[2], v[3], v[4]],
                        &[v[3], v[0], v[4]],
                    ];
                    soup.add_cell_in_zone(f, z);
                }
                CellKind::Polyhedron => unreachable!(),
            }
        }
        drop((poly, cell_nodes, cell_kind, cell_zone, buf, ends));
        let zones = if zones.is_empty() {
            vec!["fluid".to_string()]
        } else {
            zones
        };
        soup.set_zone_names(zones);

        // boundary faces by their sorted nodes
        let mut warnings = Vec::new();
        let mut named: HashMap<Vec<u32>, u32> = HashMap::new();
        let mut skipped = 0usize;
        for (f, &b) in face_boundary.iter().enumerate() {
            let nodes = &face_nodes[face_ptr[f] as usize..face_ptr[f + 1] as usize];
            let mut k = Vec::with_capacity(nodes.len());
            let mut ok = true;
            for &n in nodes {
                let m = new_index.get(n as usize).copied().ok_or_else(|| {
                    format!("a boundary face refers to node {n}, beyond the {np} nodes")
                })?;
                if m == u32::MAX {
                    ok = false;
                }
                k.push(m);
            }
            if !ok {
                skipped += 1;
                continue;
            }
            k.sort_unstable();
            named.insert(k, b);
        }
        if skipped > 0 {
            warnings.push(format!(
                "{skipped} {} are not faces of any cell and are ignored",
                face_word
            ));
        }
        let mut boundary_names = boundaries;
        let default = boundary_names.len() as u32;
        boundary_names.push("default".into());
        let unnamed = std::cell::Cell::new(0usize);
        let used = std::cell::Cell::new(0usize);
        let n_face_boundary = named.len();
        let names: Vec<&str> = boundary_names.iter().map(|s| s.as_str()).collect();
        let mut mesh = soup.finish(&names, |bf| {
            let mut k = bf.nodes.to_vec();
            k.sort_unstable();
            let b = match named.get(&k) {
                Some(&b) => {
                    used.set(used.get() + 1);
                    b
                }
                None => {
                    unnamed.set(unnamed.get() + 1);
                    default
                }
            };
            boundary_names[b as usize].clone()
        })?;
        if unnamed.get() > 0 {
            warnings.push(format!(
                "{} boundary faces belong to no {} and are in the boundary 'default'",
                unnamed.get(),
                group_word
            ));
        }
        if used.get() < n_face_boundary {
            warnings.push(format!(
                "{} {} lie inside the domain (internal surfaces) and are ignored",
                n_face_boundary - used.get(),
                face_word
            ));
        }
        mesh.boundaries.retain(|p| p.len > 0);
        Ok(Built {
            mesh,
            warnings,
            node_map: new_index,
        })
    }
}

fn index_of(list: &mut Vec<String>, name: &str) -> u32 {
    match list.iter().position(|x| x == name) {
        Some(i) => i as u32,
        None => {
            list.push(name.to_string());
            (list.len() - 1) as u32
        }
    }
}
