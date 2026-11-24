//! Reader for VTKHDF files as meshes: unstructured grids, whole or in
//! partitions, the blocked layout that KoflCraft hands over included.
//!
//! Written from the VTK documentation ("VTK File Formats", section
//! "VTKHDF File Format", version 2.4) and the blocked layout that
//! MurmelCraft and KoflCraft share (recipe 134). The file is opened
//! through HDF5, so any valid HDF5 layout reads (contiguous or chunked,
//! compressed or not).
//!
//! ```text
//! /VTKHDF            Version [1 or 2, minor], Type UnstructuredGrid or
//!                    PartitionedDataSetCollection
//!   internal         Type UnstructuredGrid: the cells, in one partition or many
//!     NumberOfPoints, NumberOfCells, NumberOfConnectivityIds   one per partition
//!     Points [Σ points, 3], Types, Offsets (cells + 1 per partition), Connectivity
//!     NumberOfFaces, NumberOfPolyhedronToFaceIds, NumberOfFaceConnectivityIds,
//!     FaceOffsets, FaceConnectivity, PolyhedronOffsets, PolyhedronToFaces
//!     CellData/      zone (any number type), vtkGhostType, OriginalCellIds
//!     PointData/     GlobalPointIds
//!     FieldData/     zone_names, boundary_names ("id=name")
//!   <patch>          Type PolyData: the faces of one boundary, in Polygons
//!     PointData/     GlobalPointIds
//! ```
//!
//! - Cells: tetrahedra, voxels, hexahedra, wedges, pyramids and polyhedra
//!   (the 2.4 arrays); triangles, quadrilaterals and polygons in the grid
//!   are boundary faces, as in a `.vtu`.
//! - Cells marked as ghosts (`CellData/vtkGhostType` with bit 1, the
//!   halo of each block) are dropped, so every cell is read once.
//! - Points of different partitions are joined by `PointData/GlobalPointIds`;
//!   a file of several partitions without them is joined by coordinates.
//! - Cells keep the order of `CellData/OriginalCellIds` when the file has
//!   them (the writer's own cell number), else their order in the file: a
//!   KoflCraft `foehncraft.vtkhdf` reads as the same mesh, bit for bit, as
//!   the `foehncraft.vtu` of the same export.
//! - Zones from the cell array `zone` (or the other region names of the
//!   VTK reader), named by `FieldData/zone_names`.
//! - Boundaries: each `PolyData` block in the order of its `Index`, named
//!   after the block; its faces find their nodes by `GlobalPointIds`, or by
//!   their coordinates when the block has none (a file without them must
//!   not have coincident distinct nodes on a boundary). Roles and the
//!   writer's root attributes are not part of the mesh and are left out.
//! - Memory: the partitions are read twice, one at a time. The first pass
//!   reads the counts of each cell and the cell order; the second puts the
//!   nodes of each cell straight where they belong in the arrays of the
//!   mesh. The file is never read whole.

use super::vtk::{REGION_NAMES, id_names};
use super::{
    CellKind, ElementMesh, MeshFile, MeshFormat, Polyhedra, ReadOptions, check_face_stream,
    declared,
};
use hdf5_metno as h5;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// The VTKHDF format in the registry of [`crate::format`].
#[derive(Debug, Clone, Copy)]
pub struct VtkHdf;

impl MeshFormat for VtkHdf {
    fn name(&self) -> &'static str {
        "vtkhdf"
    }

    fn title(&self) -> &'static str {
        "VTKHDF (HDF5), unstructured grids in one or many partitions, linear cells and polyhedra"
    }

    fn extensions(&self) -> &'static [&'static str] {
        &["vtkhdf"]
    }

    /// An HDF5 file is a CGNS file unless its name says otherwise.
    fn sniff(&self, _head: &[u8]) -> bool {
        false
    }

    fn read_path(&self, path: &Path, options: &ReadOptions) -> Result<MeshFile, String> {
        read_file(path, options)
    }
}

fn e<T: std::fmt::Display>(x: T) -> String {
    x.to_string()
}

/// The strings of a dataset or attribute (fixed or variable length, ASCII
/// or UTF-8), without trailing NULs.
fn texts(c: &h5::Container) -> Result<Vec<String>, String> {
    use h5::types::{FixedAscii, FixedUnicode, TypeDescriptor as T, VarLenAscii, VarLenUnicode};
    let d = c.dtype().and_then(|t| t.to_descriptor()).map_err(e)?;
    let v: Vec<String> = match d {
        T::FixedAscii(n) if n <= 1024 => c
            .read_raw::<FixedAscii<1024>>()
            .map_err(e)?
            .iter()
            .map(|s| s.as_str().to_string())
            .collect(),
        T::FixedUnicode(n) if n <= 1024 => c
            .read_raw::<FixedUnicode<1024>>()
            .map_err(e)?
            .iter()
            .map(|s| s.as_str().to_string())
            .collect(),
        T::VarLenAscii => c
            .read_raw::<VarLenAscii>()
            .map_err(e)?
            .iter()
            .map(|s| s.as_str().to_string())
            .collect(),
        T::VarLenUnicode => c
            .read_raw::<VarLenUnicode>()
            .map_err(e)?
            .iter()
            .map(|s| s.as_str().to_string())
            .collect(),
        other => return Err(format!("not a string ({other:?})")),
    };
    Ok(v.into_iter()
        .map(|s| s.trim_end_matches('\0').to_string())
        .collect())
}

fn text_attr(g: &h5::Group, name: &str) -> Option<String> {
    let a = g.attr(name).ok()?;
    texts(&a).ok()?.into_iter().next()
}

fn int_attr(g: &h5::Group, name: &str) -> Option<i64> {
    g.attr(name).ok()?.read_raw::<i64>().ok()?.first().copied()
}

/// A dataset of a group, if there is one.
fn dataset(g: &h5::Group, path: &str) -> Option<h5::Dataset> {
    let (dir, name) = path.rsplit_once('/').unwrap_or(("", path));
    if !dir.is_empty() && !g.link_exists(dir) {
        return None;
    }
    let parent = if dir.is_empty() {
        g.clone()
    } else {
        g.group(dir).ok()?
    };
    if !parent.link_exists(name) {
        return None;
    }
    parent.dataset(name).ok()
}

fn need(g: &h5::Group, path: &str) -> Result<h5::Dataset, String> {
    dataset(g, path).ok_or_else(|| format!("no dataset {path}"))
}

/// Rows `a..b` of a 1-D dataset, converted to `T` by HDF5.
fn rows<T: h5::H5Type>(d: &h5::Dataset, a: usize, b: usize, what: &str) -> Result<Vec<T>, String> {
    if b > d.size() {
        return Err(format!(
            "{what}: rows {a}..{b} asked, the dataset has {}",
            d.size()
        ));
    }
    if a == b {
        return Ok(Vec::new());
    }
    d.read_slice_1d::<T, _>(a..b)
        .map(|v| v.into_raw_vec_and_offset().0)
        .map_err(|x| format!("{what}: {x}"))
}

/// A count per partition (`NumberOf...`), each checked against the file.
fn counts(g: &h5::Group, name: &str, len: usize) -> Result<Vec<usize>, String> {
    let v: Vec<i64> = need(g, name)?.read_raw().map_err(e)?;
    v.into_iter()
        .map(|n| {
            usize::try_from(n)
                .map_err(|_| format!("{name}: negative count {n}"))
                .and_then(|n| declared(n, len, name))
        })
        .collect()
}

/// An index read from the file, checked to be in `0..n`.
fn index(v: i64, n: usize, what: &str) -> Result<usize, String> {
    usize::try_from(v)
        .ok()
        .filter(|&i| i < n)
        .ok_or_else(|| format!("{what} {v} is out of range (0..{n})"))
}

/// Starts of the rows of each partition: `per` rows each, plus `extra`
/// rows per partition (the `+ 1` of the offsets).
fn starts(n: &[usize], extra: usize) -> Vec<usize> {
    let mut s = Vec::with_capacity(n.len() + 1);
    let mut at = 0;
    s.push(0);
    for &k in n {
        at += k + extra;
        s.push(at);
    }
    s
}

/// The datasets of an unstructured grid and the rows of each partition.
struct Grid {
    g: h5::Group,
    n_cells: Vec<usize>,
    n_points: Vec<usize>,
    n_conn: Vec<usize>,
    n_faces: Vec<usize>,
    n_p2f: Vec<usize>,
    n_fconn: Vec<usize>,
    cell0: Vec<usize>,
    point0: Vec<usize>,
    conn0: Vec<usize>,
    face0: Vec<usize>,
    p2f0: Vec<usize>,
    fconn0: Vec<usize>,
    poly: bool,
    region: Option<h5::Dataset>,
    ghost: Option<h5::Dataset>,
    original: Option<h5::Dataset>,
    point_ids: Option<h5::Dataset>,
}

/// The arrays of one partition.
#[derive(Default)]
struct Part {
    types: Vec<u8>,
    offsets: Vec<i64>,
    ghost: Vec<u8>,
    region: Vec<i64>,
    original: Vec<i64>,
    point_ids: Vec<i64>,
    conn: Vec<i64>,
    poly_offsets: Vec<i64>,
    p2f: Vec<i64>,
    face_offsets: Vec<i64>,
    fconn: Vec<i64>,
    points: Vec<f64>,
}

/// What to read of a partition.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Pass {
    /// Counts, order and point ids.
    Sizes,
    /// Everything.
    Nodes,
}

impl Grid {
    fn open(g: h5::Group, len: usize) -> Result<Self, String> {
        let n_cells = counts(&g, "NumberOfCells", len)?;
        let np = n_cells.len();
        let n_points = counts(&g, "NumberOfPoints", len)?;
        let n_conn = counts(&g, "NumberOfConnectivityIds", len)?;
        if n_points.len() != np || n_conn.len() != np {
            return Err(format!(
                "{np} partitions in NumberOfCells, {} in NumberOfPoints, {} in \
                 NumberOfConnectivityIds",
                n_points.len(),
                n_conn.len()
            ));
        }
        let poly = g.link_exists("PolyhedronOffsets");
        let (n_faces, n_p2f, n_fconn) = if poly {
            let f = counts(&g, "NumberOfFaces", len)?;
            let t = counts(&g, "NumberOfPolyhedronToFaceIds", len)?;
            let q = counts(&g, "NumberOfFaceConnectivityIds", len)?;
            if f.len() != np || t.len() != np || q.len() != np {
                return Err("the polyhedron counts do not have one value per partition".into());
            }
            (f, t, q)
        } else {
            (vec![0; np], vec![0; np], vec![0; np])
        };
        let total = |v: &[usize], what: &str| declared(v.iter().sum(), len, what);
        total(&n_cells, "cells")?;
        total(&n_points, "points")?;
        total(&n_conn, "connectivity ids")?;
        total(&n_faces, "polyhedron faces")?;
        total(&n_p2f, "polyhedron face ids")?;
        total(&n_fconn, "face connectivity ids")?;
        let region = REGION_NAMES.iter().find_map(|want| {
            let cd = g.group("CellData").ok()?;
            let name = cd
                .member_names()
                .ok()?
                .into_iter()
                .find(|n| n.eq_ignore_ascii_case(want))?;
            cd.dataset(&name).ok()
        });
        Ok(Self {
            cell0: starts(&n_cells, 0),
            point0: starts(&n_points, 0),
            conn0: starts(&n_conn, 0),
            face0: starts(&n_faces, 0),
            p2f0: starts(&n_p2f, 0),
            fconn0: starts(&n_fconn, 0),
            ghost: dataset(&g, "CellData/vtkGhostType"),
            original: dataset(&g, "CellData/OriginalCellIds"),
            point_ids: dataset(&g, "PointData/GlobalPointIds"),
            region,
            poly,
            n_cells,
            n_points,
            n_conn,
            n_faces,
            n_p2f,
            n_fconn,
            g,
        })
    }

    fn partitions(&self) -> usize {
        self.n_cells.len()
    }

    fn read(&self, p: usize, pass: Pass) -> Result<Part, String> {
        let g = &self.g;
        let (c0, nc) = (self.cell0[p], self.n_cells[p]);
        // the offsets of partition p start after the nc + 1 rows of each
        // partition before it
        let o0 = c0 + p;
        let mut part = Part {
            types: rows(&need(g, "Types")?, c0, c0 + nc, "Types")?,
            offsets: rows(&need(g, "Offsets")?, o0, o0 + nc + 1, "Offsets")?,
            ..Default::default()
        };
        if let Some(d) = &self.ghost {
            part.ghost = rows(d, c0, c0 + nc, "vtkGhostType")?;
        }
        if let Some(d) = &self.region {
            // read as reals: HDF5 converts any number type to them, and
            // region numbers are small integers
            let v: Vec<f64> = rows(d, c0, c0 + nc, "the region array")?;
            part.region = v.into_iter().map(|x| x as i64).collect();
        }
        if let Some(d) = &self.original {
            part.original = rows(d, c0, c0 + nc, "OriginalCellIds")?;
        }
        let (p0, npt) = (self.point0[p], self.n_points[p]);
        if let Some(d) = &self.point_ids {
            part.point_ids = rows(d, p0, p0 + npt, "GlobalPointIds")?;
        }
        let poly = self.poly && part.types.contains(&42);
        if poly {
            let (f0, nf) = (self.face0[p] + p, self.n_faces[p]);
            let (t0, nt) = (self.p2f0[p], self.n_p2f[p]);
            part.poly_offsets = rows(
                &need(g, "PolyhedronOffsets")?,
                o0,
                o0 + nc + 1,
                "PolyhedronOffsets",
            )?;
            part.p2f = rows(
                &need(g, "PolyhedronToFaces")?,
                t0,
                t0 + nt,
                "PolyhedronToFaces",
            )?;
            part.face_offsets = rows(&need(g, "FaceOffsets")?, f0, f0 + nf + 1, "FaceOffsets")?;
        } else if part.types.contains(&42) {
            return Err("polyhedra need the arrays of VTKHDF 2.4 (PolyhedronOffsets, ...)".into());
        }
        if pass == Pass::Nodes {
            let (k0, nk) = (self.conn0[p], self.n_conn[p]);
            part.conn = rows(&need(g, "Connectivity")?, k0, k0 + nk, "Connectivity")?;
            if poly {
                let (q0, nq) = (self.fconn0[p], self.n_fconn[p]);
                part.fconn = rows(
                    &need(g, "FaceConnectivity")?,
                    q0,
                    q0 + nq,
                    "FaceConnectivity",
                )?;
            }
            let d = need(g, "Points")?;
            if d.ndim() != 2 || d.shape()[1] != 3 || d.shape()[0] < p0 + npt {
                return Err(format!("Points has the shape {:?}", d.shape()));
            }
            if npt > 0 {
                part.points = d
                    .read_slice_2d::<f64, _>((p0..p0 + npt, ..))
                    .map(|v| v.into_raw_vec_and_offset().0)
                    .map_err(|x| format!("Points: {x}"))?;
            }
        }
        Ok(part)
    }
}

impl Part {
    fn owned(&self, i: usize) -> bool {
        self.ghost.get(i).is_none_or(|g| g & 1 == 0)
    }

    /// The nodes of cell `i` in the connectivity.
    fn range(&self, i: usize, n: usize, what: &str) -> Result<(usize, usize), String> {
        let a = index(self.offsets[i], n + 1, what)?;
        let b = index(self.offsets[i + 1], n + 1, what)?;
        if a > b {
            return Err(format!("{what} decrease at cell {i}"));
        }
        Ok((a, b))
    }

    /// The faces of polyhedron `i`: their (start, end) in the face
    /// connectivity, each checked.
    fn faces(
        &self,
        i: usize,
        nt: usize,
        nf: usize,
        nq: usize,
    ) -> Result<Vec<(usize, usize)>, String> {
        let a = index(self.poly_offsets[i], nt + 1, "PolyhedronOffsets")?;
        let b = index(self.poly_offsets[i + 1], nt + 1, "PolyhedronOffsets")?;
        if a >= b {
            return Err(format!("polyhedron {i} of its partition has no faces"));
        }
        self.p2f[a..b]
            .iter()
            .map(|&f| {
                let f = index(f, nf, "a face of PolyhedronToFaces")?;
                let s = index(self.face_offsets[f], nq + 1, "FaceOffsets")?;
                let t = index(self.face_offsets[f + 1], nq + 1, "FaceOffsets")?;
                if t < s + 3 {
                    return Err(format!(
                        "polyhedron {i} of its partition has a face of {} points",
                        t.saturating_sub(s)
                    ));
                }
                Ok((s, t))
            })
            .collect()
    }
}

/// The kind of a VTK cell type, `None` for a boundary face and `Err` for a
/// type that is neither.
fn kind(t: u8) -> Result<Option<CellKind>, String> {
    Ok(Some(match t {
        10 => CellKind::Tetrahedron,
        11 | 12 => CellKind::Hexahedron,
        13 => CellKind::Prism,
        14 => CellKind::Pyramid,
        42 => CellKind::Polyhedron,
        5 | 7 | 9 => return Ok(None),
        other => {
            return Err(format!(
                "a cell of type {other}: only linear tetrahedra, hexahedra, wedges, pyramids \
                 and polyhedra (and triangles, quadrilaterals and polygons as boundary faces) \
                 are read"
            ));
        }
    }))
}

/// The blocks of a `PartitionedDataSetCollection`: the unstructured grid
/// and the boundary patches in the order of their `Index`.
fn blocks(root: &h5::Group) -> Result<(h5::Group, Vec<(String, h5::Group)>), String> {
    let mut grids = Vec::new();
    let mut patches = Vec::new();
    for name in root.member_names().map_err(e)? {
        if name == "Assembly" {
            continue;
        }
        let Ok(g) = root.group(&name) else { continue };
        let at = int_attr(&g, "Index").unwrap_or(i64::MAX);
        match text_attr(&g, "Type").as_deref() {
            Some("UnstructuredGrid") => grids.push((at, name, g)),
            Some("PolyData") => patches.push((at, name, g)),
            _ => {}
        }
    }
    grids.sort_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
    patches.sort_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
    let grid = match grids.len() {
        0 => return Err("no UnstructuredGrid block".into()),
        1 => grids.pop().unwrap().2,
        _ => {
            let names: Vec<&str> = grids.iter().map(|g| g.1.as_str()).collect();
            match grids.iter().position(|g| g.1 == "internal") {
                Some(i) => grids.swap_remove(i).2,
                None => {
                    return Err(format!(
                        "{} UnstructuredGrid blocks ({}): one is read, named `internal`",
                        names.len(),
                        names.join(", ")
                    ));
                }
            }
        }
    };
    Ok((grid, patches.into_iter().map(|(_, n, g)| (n, g)).collect()))
}

/// Dense node numbers from the file's point ids.
enum Ids {
    /// Partition start plus local id (no `GlobalPointIds`).
    Local,
    /// The ids themselves (all below the number of points of the file).
    Direct,
    /// Position in the sorted list of the ids.
    Sorted(Vec<i64>),
}

impl Ids {
    fn of(&self, id: i64) -> u32 {
        match self {
            Ids::Local | Ids::Direct => id as u32,
            Ids::Sorted(v) => v.binary_search(&id).unwrap_or(0) as u32,
        }
    }
}

/// Reads a VTKHDF file as a mesh.
pub fn read_file(path: &Path, options: &ReadOptions) -> Result<MeshFile, String> {
    let len = std::fs::metadata(path)
        .map_err(|x| format!("{}: {x}", path.display()))?
        .len() as usize;
    let file = h5::File::open(path).map_err(e)?;
    if !file.link_exists("VTKHDF") {
        return Err("no group /VTKHDF: not a VTKHDF file".into());
    }
    let root = file.group("VTKHDF").map_err(e)?;
    let version = root
        .attr("Version")
        .and_then(|a| a.read_raw::<i64>())
        .map_err(|_| "/VTKHDF has no attribute Version".to_string())?;
    let (major, minor) = (
        version.first().copied().unwrap_or(0),
        version.get(1).copied().unwrap_or(0),
    );
    if !(1..=2).contains(&major) {
        return Err(format!(
            "VTKHDF version {major}.{minor}: versions 1 and 2 are read"
        ));
    }
    let ty = text_attr(&root, "Type").unwrap_or_default();
    let (grid, patches) = match ty.as_str() {
        "UnstructuredGrid" => (root.clone(), Vec::new()),
        "PartitionedDataSetCollection" => blocks(&root)?,
        other => {
            return Err(format!(
                "a VTKHDF file of type `{other}`: meshes are UnstructuredGrid or \
                 PartitionedDataSetCollection files"
            ));
        }
    };
    if grid.link_exists("Steps") {
        return Err("a time-dependent VTKHDF file (group Steps) is not read as a mesh".into());
    }
    let grid = Grid::open(grid, len)?;
    let np = grid.partitions();
    let mut warnings = Vec::new();

    // names of zones and of boundaries given as polygons in the grid
    let field_names = |name: &str| -> HashMap<i64, String> {
        let list = dataset(&grid.g, &format!("FieldData/{name}")).and_then(|d| texts(&d).ok());
        id_names(list.as_ref())
    };
    let zone_names = field_names("zone_names");
    let boundary_names = field_names("boundary_names");

    // pass 1: the owned cells, their sizes, zones and order; the point ids
    let mut kinds: Vec<CellKind> = Vec::new();
    let mut sizes: Vec<u32> = Vec::new();
    let mut zone_id: Vec<i64> = Vec::new();
    let mut original: Vec<i64> = Vec::new();
    let mut zone_ids: BTreeMap<i64, ()> = BTreeMap::new();
    let mut bound_ids: BTreeMap<i64, ()> = BTreeMap::new();
    let (mut id_min, mut id_max) = (i64::MAX, i64::MIN);
    let all_points: usize = grid.n_points.iter().sum();
    for p in 0..np {
        let part = grid.read(p, Pass::Sizes)?;
        let nc = grid.n_cells[p];
        for &v in &part.point_ids {
            id_min = id_min.min(v);
            id_max = id_max.max(v);
        }
        for i in 0..nc {
            if !part.owned(i) {
                continue;
            }
            let region = part.region.get(i).copied();
            let Some(k) = kind(part.types[i])? else {
                bound_ids.insert(region.unwrap_or(i64::MIN), ());
                continue;
            };
            let size = if k == CellKind::Polyhedron {
                let f = part.faces(i, grid.n_p2f[p], grid.n_faces[p], grid.n_fconn[p])?;
                1 + f.iter().map(|&(s, t)| 1 + t - s).sum::<usize>()
            } else {
                let (a, b) = part.range(i, grid.n_conn[p], "Offsets")?;
                let want = if part.types[i] == 11 { 8 } else { k.nodes() };
                if b - a != want {
                    return Err(format!(
                        "cell {i} of partition {p} (type {}) has {} points, expected {want}",
                        part.types[i],
                        b - a
                    ));
                }
                want
            };
            kinds.push(k);
            sizes.push(u32::try_from(size).map_err(|_| "a polyhedron of 2^32 ids".to_string())?);
            zone_id.push(region.unwrap_or(0));
            zone_ids.insert(region.unwrap_or(0), ());
            if let Some(&o) = part.original.get(i) {
                original.push(o);
            }
        }
    }
    let n = kinds.len();
    if n == 0 {
        return Err("the grid has no 3D cells".into());
    }

    // dense point numbers
    let ids = match &grid.point_ids {
        None => {
            if np > 1 {
                warnings.push(
                    "the partitions have no PointData/GlobalPointIds: their points are joined \
                     by their coordinates"
                        .to_string(),
                );
            }
            Ids::Local
        }
        Some(_) if id_min >= 0 && (id_max as u64) < all_points as u64 => Ids::Direct,
        Some(d) => {
            let mut v: Vec<i64> = rows(d, 0, all_points, "GlobalPointIds")?;
            v.sort_unstable();
            v.dedup();
            Ids::Sorted(v)
        }
    };
    let n_points = match &ids {
        Ids::Local => all_points,
        Ids::Direct => (id_max + 1).max(0) as usize,
        Ids::Sorted(v) => v.len(),
    };

    // the order of the cells: their original numbers, else the file
    let at: Vec<usize> = if original.len() == n {
        let mut seen = vec![false; n];
        let permutation = original.iter().all(|&o| {
            usize::try_from(o)
                .ok()
                .filter(|&o| o < n)
                .is_some_and(|o| !std::mem::replace(&mut seen[o], true))
        });
        if permutation {
            original.iter().map(|&o| o as usize).collect()
        } else {
            let mut by: Vec<usize> = (0..n).collect();
            by.sort_by_key(|&i| (original[i], i));
            let mut at = vec![0; n];
            for (k, &i) in by.iter().enumerate() {
                at[i] = k;
            }
            at
        }
    } else {
        (0..n).collect()
    };
    drop(original);

    // where each cell's nodes go
    let single = zone_ids.len() == 1;
    let mut em = ElementMesh {
        points: vec![[0.0; 3]; n_points],
        face_word: "boundary faces",
        group_word: "patch",
        ..Default::default()
    };
    let zone_of: HashMap<i64, u32> = zone_ids
        .keys()
        .map(|&k| {
            let name = zone_names.get(&k).cloned().unwrap_or_else(|| {
                if single {
                    "fluid".to_string()
                } else {
                    format!("zone{k}")
                }
            });
            (k, em.zone(&name))
        })
        .collect();
    let mut inverse = vec![0usize; n];
    for (i, &k) in at.iter().enumerate() {
        inverse[k] = i;
    }
    em.cell_kind = inverse.iter().map(|&i| kinds[i]).collect();
    em.cell_zone = inverse.iter().map(|&i| zone_of[&zone_id[i]]).collect();
    drop((zone_id, kinds));
    // start of each cell (in file order) in `cell_nodes` or in the
    // polyhedron stream
    let mut start = vec![0usize; n];
    let (mut nodes_len, mut stream_len) = (0usize, 0usize);
    let mut poly_starts = Vec::new();
    for &i in &inverse {
        let s = sizes[i] as usize;
        if em.cell_kind[at[i]] == CellKind::Polyhedron {
            start[i] = stream_len;
            poly_starts.push(stream_len);
            stream_len += s;
        } else {
            start[i] = nodes_len;
            nodes_len += s;
        }
    }
    drop((sizes, inverse));
    em.cell_nodes = vec![0; nodes_len];
    let mut stream = vec![0u32; stream_len];
    let bound_of: HashMap<i64, u32> = bound_ids
        .keys()
        .map(|&k| {
            let name = boundary_names.get(&k).cloned().unwrap_or_else(|| {
                if k == i64::MIN {
                    "default".to_string()
                } else {
                    format!("boundary{k}")
                }
            });
            (k, em.boundary(&name))
        })
        .collect();

    // pass 2: points and nodes
    let mut j = 0usize; // owned volume cell in file order
    for p in 0..np {
        let part = grid.read(p, Pass::Nodes)?;
        let npt = grid.n_points[p];
        let node = |local: i64| -> Result<u32, String> {
            let l = index(local, npt, "a point id")?;
            Ok(match &ids {
                Ids::Local => (grid.point0[p] + l) as u32,
                _ => ids.of(part.point_ids[l]),
            })
        };
        for l in 0..npt {
            let v = node(l as i64)? as usize;
            em.points[v] = [
                part.points[3 * l],
                part.points[3 * l + 1],
                part.points[3 * l + 2],
            ];
        }
        for i in 0..grid.n_cells[p] {
            if !part.owned(i) {
                continue;
            }
            let t = part.types[i];
            let Some(k) = kind(t)? else {
                let (a, b) = part.range(i, grid.n_conn[p], "Offsets")?;
                if b - a < 3 {
                    return Err(format!("face {i} of partition {p} has {} points", b - a));
                }
                let f: Vec<u32> = part.conn[a..b]
                    .iter()
                    .map(|&v| node(v))
                    .collect::<Result<_, _>>()?;
                let region = part.region.get(i).copied().unwrap_or(i64::MIN);
                em.add_face(&f, bound_of[&region]);
                continue;
            };
            if k == CellKind::Polyhedron {
                let faces = part.faces(i, grid.n_p2f[p], grid.n_faces[p], grid.n_fconn[p])?;
                let mut w = start[j];
                stream[w] = faces.len() as u32;
                w += 1;
                for (s, t) in faces {
                    stream[w] = (t - s) as u32;
                    w += 1;
                    for &v in &part.fconn[s..t] {
                        stream[w] = node(v)?;
                        w += 1;
                    }
                }
            } else {
                let (a, b) = part.range(i, grid.n_conn[p], "Offsets")?;
                let s = start[j];
                for (m, &v) in part.conn[a..b].iter().enumerate() {
                    // a voxel's nodes in hexahedron order
                    let m = if t == 11 {
                        [0, 1, 3, 2, 4, 5, 7, 6][m]
                    } else {
                        m
                    };
                    em.cell_nodes[s + m] = node(v)?;
                }
            }
            j += 1;
        }
    }
    drop(start);
    for (k, &s) in poly_starts.iter().enumerate() {
        let end = poly_starts.get(k + 1).copied().unwrap_or(stream.len());
        check_face_stream(&stream[s..end], k)?;
    }
    em.poly = Polyhedra::from_parts(stream, poly_starts);

    // boundary patches
    let mut by_coordinates: Option<HashMap<[u64; 3], u32>> = None;
    for (name, g) in &patches {
        let b = em.boundary(name);
        read_patch(
            name,
            g,
            b,
            &ids,
            len,
            &mut em,
            &mut by_coordinates,
            &mut warnings,
        )?;
    }
    drop(by_coordinates);
    if options.merge_nodes || (matches!(ids, Ids::Local) && np > 1) {
        em.merge_points();
    }

    let blocked = int_attr(&root, "murmelcraft_blocks")
        .map(|v| format!(", blocked layout {v}"))
        .unwrap_or_default();
    let how = format!(
        "VTKHDF {major}.{minor}, {np} partition{}{blocked}",
        if np == 1 { "" } else { "s" }
    );
    let built = em.into_mesh()?;
    warnings.extend(built.warnings);
    Ok(MeshFile {
        mesh: built.mesh,
        reader: "vtkhdf",
        format: how,
        warnings,
    })
}

/// Reads the faces of one `PolyData` block into boundary `b`.
#[allow(clippy::too_many_arguments)]
fn read_patch(
    name: &str,
    g: &h5::Group,
    b: u32,
    ids: &Ids,
    len: usize,
    em: &mut ElementMesh,
    by_coordinates: &mut Option<HashMap<[u64; 3], u32>>,
    warnings: &mut Vec<String>,
) -> Result<(), String> {
    let err = |x: String| format!("block {name}: {x}");
    let n_points = counts(g, "NumberOfPoints", len).map_err(err)?;
    let poly = g
        .group("Polygons")
        .map_err(|_| err("no group Polygons".into()))?;
    let n_cells = counts(&poly, "NumberOfCells", len).map_err(err)?;
    let n_conn = counts(&poly, "NumberOfConnectivityIds", len).map_err(err)?;
    if n_cells.len() != n_points.len() || n_conn.len() != n_points.len() {
        return Err(err("the counts do not have one value per partition".into()));
    }
    for other in ["Strips", "Lines", "Vertices"] {
        if let Ok(o) = g.group(other)
            && let Ok(c) = counts(&o, "NumberOfCells", len)
            && c.iter().any(|&k| k > 0)
        {
            warnings.push(format!(
                "block {name}: its {} are not boundary faces and are left out",
                other.to_lowercase()
            ));
        }
    }
    let global = dataset(g, "PointData/GlobalPointIds");
    if global.is_none() && !matches!(ids, Ids::Local) {
        return Err(err(
            "no PointData/GlobalPointIds, which the cells have: its faces cannot be joined to them"
                .into(),
        ));
    }
    let (p0, c0, k0) = (
        starts(&n_points, 0),
        starts(&n_cells, 0),
        starts(&n_conn, 0),
    );
    let points = need(g, "Points").map_err(err)?;
    let offsets = need(&poly, "Offsets").map_err(err)?;
    let conn = need(&poly, "Connectivity").map_err(err)?;
    for p in 0..n_points.len() {
        let npt = n_points[p];
        // the mesh node of each point of the partition
        let node: Vec<u32> = match &global {
            Some(d) => rows::<i64>(d, p0[p], p0[p] + npt, "GlobalPointIds")
                .map_err(err)?
                .into_iter()
                .map(|v| match ids {
                    Ids::Sorted(s) => s
                        .binary_search(&v)
                        .map(|i| i as u32)
                        .map_err(|_| err(format!("point id {v} is not a point of the cells"))),
                    _ => index(v, em.points.len(), "a point id")
                        .map(|i| i as u32)
                        .map_err(err),
                })
                .collect::<Result<_, _>>()?,
            None => {
                let map = by_coordinates.get_or_insert_with(|| {
                    let mut m = HashMap::with_capacity(em.points.len());
                    for (i, q) in em.points.iter().enumerate() {
                        m.entry(q.map(f64::to_bits)).or_insert(i as u32);
                    }
                    m
                });
                let xyz: Vec<f64> = if npt == 0 {
                    Vec::new()
                } else {
                    points
                        .read_slice_2d::<f64, _>((p0[p]..p0[p] + npt, ..))
                        .map(|v| v.into_raw_vec_and_offset().0)
                        .map_err(|x| err(format!("Points: {x}")))?
                };
                xyz.as_chunks::<3>()
                    .0
                    .iter()
                    .map(|q| {
                        map.get(&[q[0], q[1], q[2]].map(f64::to_bits))
                            .copied()
                            .ok_or_else(|| {
                                err(format!("the point {q:?} is not a point of the cells"))
                            })
                    })
                    .collect::<Result<_, _>>()?
            }
        };
        let o: Vec<i64> =
            rows(&offsets, c0[p] + p, c0[p] + p + n_cells[p] + 1, "Offsets").map_err(err)?;
        let c: Vec<i64> = rows(&conn, k0[p], k0[p] + n_conn[p], "Connectivity").map_err(err)?;
        let mut f = Vec::new();
        for i in 0..n_cells[p] {
            let a = index(o[i], c.len() + 1, "an offset").map_err(err)?;
            let z = index(o[i + 1], c.len() + 1, "an offset").map_err(err)?;
            if z < a + 3 {
                return Err(err(format!("face {i} has {} points", z.saturating_sub(a))));
            }
            f.clear();
            for &v in &c[a..z] {
                f.push(node[index(v, npt, "a point id").map_err(err)?]);
            }
            em.add_face(&f, b);
        }
    }
    Ok(())
}
