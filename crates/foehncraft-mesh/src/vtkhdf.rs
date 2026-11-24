//! VTKHDF output (D41).
//!
//! One HDF5 file in the VTKHDF format of VTK, read directly by ParaView
//! (open the `.vtkhdf` file; reader "VTKHDF"), in the blocked layout that
//! MurmelCraft reads in parts (layout 1, shared with KoflCraft's hand-off).
//! The file is a `PartitionedDataSetCollection`:
//!
//! - block `internal`: the cells as an `UnstructuredGrid` in partitions of
//!   [`BLOCK_CELLS`] owned cells, the cells sorted along a Hilbert curve
//!   through their centroids (Skilling 2004), each partition with one layer
//!   of ghost cells (the cells of other blocks that share a point with its
//!   own, `vtkGhostType` 1, carrying the same field values); tetrahedra,
//!   hexahedra, prisms and pyramids as VTK's linear cells, everything else
//!   as `VTK_POLYHEDRON` with its faces; cell fields, `zone`,
//!   `GlobalCellIds` (the cell's place in the file), `OriginalCellIds` (its
//!   global identity) under `CellData`, `GlobalPointIds` under `PointData`,
//!   and per block the owned cells, bounds and field ranges under
//!   `FieldData`;
//! - one `PolyData` block per boundary, named after the boundary, with
//!   the boundary index, `OwnerCell` and the boundary face fields
//!   ([`FaceField`], such as the wall shear stress) as cell fields;
//! - an assembly `internal`, `boundary/<boundary>` for ParaView's block
//!   selection.
//!
//! The order does not depend on a renumbering of the cells: blocks follow
//! the centroids (ties by the global identities, input order, R2), the
//! points of a block and the boundary faces follow their identities.
//!
//! The writer streams: the datasets are created at their full size from a
//! first pass that counts each block, then each block's rows are built and
//! written in turn, so the file never exists in memory whole. A million
//! polyhedra with 6.2 million points take 0.26 GB beside the mesh and the
//! fields while they are written (1.26 GB of whole-file arrays before).
//! Points that no cell uses are not written.
//!
//! The format is described in the VTK documentation, "VTK File Formats",
//! section "VTKHDF File Format". This writer is serial; the parallel writer
//! (all ranks into one file through MPI-IO) follows with the solver stages.

use crate::mesh::{CellShape, Mesh};
use anyhow::{Context, Result};
use h5::types::FixedAscii;
use hdf5_metno as h5;
use std::path::Path;

/// VTKHDF format version written.
pub const VERSION: [i64; 2] = [2, 4];

const VTK_TETRA: u8 = 10;
const VTK_HEXAHEDRON: u8 = 12;
const VTK_WEDGE: u8 = 13;
const VTK_PYRAMID: u8 = 14;
const VTK_POLYHEDRON: u8 = 42;

/// A cell field to write.
pub struct CellField<'a> {
    /// Name shown in ParaView.
    pub name: &'a str,
    /// One value per local cell.
    pub values: &'a [f64],
}

/// A field on the boundary faces, written to the boundary blocks.
pub struct FaceField<'a> {
    /// Name shown in ParaView.
    pub name: &'a str,
    /// One value per local boundary face (boundary face `f` of the mesh at
    /// `f - n_interior_faces`).
    pub values: &'a [f64],
}

/// A file attribute of [`write_with`].
#[derive(Debug, Clone, Copy)]
pub enum Attr<'a> {
    /// Text.
    Str(&'a str),
    /// Real number.
    F64(f64),
    /// Count.
    U64(u64),
}

/// Writes the mesh and cell fields to `path`.
pub fn write(path: impl AsRef<Path>, mesh: &Mesh, fields: &[CellField]) -> Result<()> {
    write_with(path, mesh, fields, &[])
}

/// As [`write()`], with attributes of the file's root group that describe it
/// (D75: format, version, case, step, time, units). VTK reads only the
/// `/VTKHDF` group, so they do not disturb ParaView. The file holds no
/// object timestamps, so the same content gives the same bytes (D84).
pub fn write_with(
    path: impl AsRef<Path>,
    mesh: &Mesh,
    fields: &[CellField],
    attrs: &[(&str, Attr)],
) -> Result<()> {
    write_all(path, mesh, fields, &[], attrs)
}

/// Owned cells per block of the blocked layout (the ghosts come on top).
pub const BLOCK_CELLS: usize = 65_536;

/// As [`write_with`], with fields on the boundary faces, which go to the
/// boundary blocks.
pub fn write_all(
    path: impl AsRef<Path>,
    mesh: &Mesh,
    fields: &[CellField],
    face_fields: &[FaceField],
    attrs: &[(&str, Attr)],
) -> Result<()> {
    write_blocks(path, mesh, fields, face_fields, attrs, BLOCK_CELLS)
}

/// As [`write_all`] with `block` owned cells per block (small blocks for
/// tests; files for MurmelCraft keep [`BLOCK_CELLS`]).
#[doc(hidden)]
pub fn write_blocks(
    path: impl AsRef<Path>,
    mesh: &Mesh,
    fields: &[CellField],
    face_fields: &[FaceField],
    attrs: &[(&str, Attr)],
    block: usize,
) -> Result<()> {
    let path = path.as_ref();
    anyhow::ensure!(block > 0, "blocks of 0 cells");
    let nb = mesh.n_boundary_faces();
    for f in face_fields {
        anyhow::ensure!(
            f.values.len() == nb,
            "face field {} has {} values for {nb} boundary faces",
            f.name,
            f.values.len(),
        );
    }
    for f in fields {
        anyhow::ensure!(
            f.values.len() == mesh.n_cells,
            "field {} has {} values for {} cells",
            f.name,
            f.values.len(),
            mesh.n_cells
        );
    }
    let file = h5::File::with_options()
        .with_fcpl(|p| p.obj_track_times(false))
        .create(path)
        .with_context(|| format!("creating {}", path.display()))?;
    for (name, a) in attrs {
        match a {
            Attr::Str(v) => {
                let v: h5::types::VarLenUnicode =
                    v.parse().map_err(|e| anyhow::anyhow!("{e:?}"))?;
                file.new_attr::<h5::types::VarLenUnicode>()
                    .create(*name)?
                    .write_scalar(&v)?;
            }
            Attr::F64(v) => file.new_attr::<f64>().create(*name)?.write_scalar(v)?,
            Attr::U64(v) => file.new_attr::<u64>().create(*name)?.write_scalar(v)?,
        }
    }
    let root = group(&file, "VTKHDF")?;
    ascii_attr(&root, "Type", "PartitionedDataSetCollection")?;
    version_attr(&root)?;
    root.new_attr::<i64>()
        .create("murmelcraft_blocks")?
        .write_scalar(&1)?;
    ascii_attr(&root, "block_order", "hilbert")?;

    let layout = Layout::new(mesh, block);
    let internal = group(&root, "internal")?;
    write_internal(&internal, mesh, fields, &layout)?;
    index_attr(&internal, 0)?;
    internal
        .new_attr::<i64>()
        .create("block_cells")?
        .write_scalar(&(block as i64))?;

    let assembly = group(&root, "Assembly")?;
    assembly.link_soft("/VTKHDF/internal", "internal")?;
    let boundary = group(&assembly, "boundary")?;
    for (i, p) in mesh.boundaries.iter().enumerate() {
        let name = block_name(&p.name);
        let g = group(&root, &name)?;
        write_boundary(&g, mesh, p.range(), i as f64, face_fields, &layout.file_id)?;
        index_attr(&g, i as i64 + 1)?;
        boundary.link_soft(&format!("/VTKHDF/{name}"), &name)?;
    }
    Ok(())
}

/// A field on points (parcels): one value or one 3-vector per point.
pub struct PointField<'a> {
    /// Name shown in ParaView.
    pub name: &'a str,
    /// `n` values (scalar) or `3 n` values (vector, x y z of point 0 first).
    pub values: &'a [f64],
    /// 1 or 3.
    pub components: usize,
}

/// Writes points (the parcels of the sprays, L12) as a VTKHDF `PolyData`
/// file: one vertex cell per point, the fields as point data. The file
/// holds no object timestamps, so the same content gives the same bytes
/// (D84).
pub fn write_points(
    path: impl AsRef<Path>,
    points: &[[f64; 3]],
    fields: &[PointField],
    attrs: &[(&str, Attr)],
) -> Result<()> {
    let path = path.as_ref();
    let n = points.len();
    for f in fields {
        anyhow::ensure!(
            (f.components == 1 || f.components == 3) && f.values.len() == n * f.components,
            "point field {} has {} values for {n} points of {} components",
            f.name,
            f.values.len(),
            f.components
        );
    }
    let file = h5::File::with_options()
        .with_fcpl(|p| p.obj_track_times(false))
        .create(path)
        .with_context(|| format!("creating {}", path.display()))?;
    for (name, a) in attrs {
        match a {
            Attr::Str(v) => {
                let v: h5::types::VarLenUnicode =
                    v.parse().map_err(|e| anyhow::anyhow!("{e:?}"))?;
                file.new_attr::<h5::types::VarLenUnicode>()
                    .create(*name)?
                    .write_scalar(&v)?;
            }
            Attr::F64(v) => file.new_attr::<f64>().create(*name)?.write_scalar(v)?,
            Attr::U64(v) => file.new_attr::<u64>().create(*name)?.write_scalar(v)?,
        }
    }
    let g = group(&file, "VTKHDF")?;
    ascii_attr(&g, "Type", "PolyData")?;
    version_attr(&g)?;
    ds(&g, "NumberOfPoints", &[n as i64])?;
    points_ds(&g, points)?;
    let conn: Vec<i64> = (0..n as i64).collect();
    let offsets: Vec<i64> = (0..=n as i64).collect();
    for (name, cells, o, c) in [
        ("Vertices", n, offsets, conn),
        ("Lines", 0, vec![0i64], vec![]),
        ("Strips", 0, vec![0], vec![]),
        ("Polygons", 0, vec![0], vec![]),
    ] {
        let sg = group(&g, name)?;
        ds(&sg, "NumberOfCells", &[cells as i64])?;
        ds(&sg, "NumberOfConnectivityIds", &[c.len() as i64])?;
        ds(&sg, "Offsets", &o)?;
        ds(&sg, "Connectivity", &c)?;
    }
    let pd = group(&g, "PointData")?;
    for f in fields {
        if f.components == 1 {
            ds(&pd, f.name, f.values)?;
        } else {
            pd.new_dataset::<f64>()
                .shape([n, 3])
                .create(f.name)?
                .write_raw(f.values)?;
        }
    }
    group(&g, "CellData")?;
    group(&g, "FieldData")?;
    Ok(())
}

/// The size of the file that [`write()`] writes for `mesh` with `n_fields`
/// cell fields, from the arrays it stores (an estimate within a few per
/// cent: HDF5 adds headers of some kilobytes). The blocks and their ghosts
/// are counted as the writer lays them out.
pub fn estimate_bytes(mesh: &Mesh, n_fields: usize) -> f64 {
    let layout = Layout::new(mesh, BLOCK_CELLS);
    let n = layout.counts(&sizes(mesh));
    let nblocks = n.cells.len();
    let sum = |v: &[usize]| v.iter().sum::<usize>();
    let (cells, points) = (sum(&n.cells), sum(&n.points));
    // per cell: type, two offsets, zone, two ids, ghost flag, fields
    let per_cell = 1 + 8 + 8 + 4 + 8 + 8 + 1 + 8 * n_fields;
    // per point: coordinates, id, ghost flag
    let per_point = 24 + 8 + 1;
    // connectivity, the faces of the polyhedra (offset and reference) and
    // their connectivity
    let ids = sum(&n.conn) + 2 * sum(&n.faces) + sum(&n.fconn);
    // per block: six counts, three closing offsets, owned cells, bounds,
    // field ranges
    let index = nblocks * (9 * 8 + 8 + 48 + 16 * n_fields);
    // boundary blocks: points with their ids, the face connectivity, per
    // face an offset, the boundary index and the owner cell
    let mut boundary = 0usize;
    for p in &mesh.boundaries {
        let mut points = std::collections::HashSet::new();
        let mut nodes = 0usize;
        for f in p.range() {
            let fnodes = mesh.face(f);
            nodes += fnodes.len();
            points.extend(fnodes.iter().copied());
        }
        boundary += points.len() * (24 + 8) + nodes * 8 + p.len * (8 + 8 + 8);
    }
    (cells * per_cell + points * per_point + ids * 8 + index + boundary) as f64
        + 8192.0 * (3 + mesh.boundaries.len()) as f64
}

/// Block names must not clash with the reserved group names.
fn block_name(boundary: &str) -> String {
    match boundary {
        "internal" | "Assembly" => format!("boundary_{boundary}"),
        _ => boundary.replace('/', "_"),
    }
}

/// A group that tracks and indexes the creation order of its links: VTK
/// iterates over blocks and assembly nodes in creation order.
fn group(parent: &h5::Group, name: &str) -> Result<h5::Group> {
    use h5::plist::group_create::LinkCreationOrder;
    Ok(parent
        .create_group_builder()
        .with_gcpl(|p| {
            p.link_creation_order(LinkCreationOrder::Indexed)
                .obj_track_times(false)
        })
        .create(name)?)
}

/// A fixed-length ASCII attribute of exactly the text's length, as VTK
/// writes `Type`.
fn ascii_attr(g: &h5::Group, name: &str, t: &str) -> Result<()> {
    macro_rules! fixed {
        ($($n:literal)*) => {
            match t.len() {
                $($n => {
                    let v = FixedAscii::<$n>::from_ascii(t.as_bytes())?;
                    g.new_attr::<FixedAscii<$n>>().create(name)?.write_scalar(&v)?;
                })*
                _ => anyhow::bail!("unexpected attribute text {t}"),
            }
        };
    }
    fixed!(7 8 16 28);
    Ok(())
}

fn version_attr(g: &h5::Group) -> Result<()> {
    g.new_attr::<i64>()
        .shape([2])
        .create("Version")?
        .write_raw(&VERSION)?;
    Ok(())
}

fn index_attr(g: &h5::Group, i: i64) -> Result<()> {
    g.new_attr::<i64>().create("Index")?.write_scalar(&i)?;
    Ok(())
}

fn ds<T: h5::H5Type>(g: &h5::Group, name: &str, data: &[T]) -> Result<()> {
    g.new_dataset::<T>()
        .shape([data.len()])
        .create(name)?
        .write_raw(data)?;
    Ok(())
}

fn points_ds(g: &h5::Group, pts: &[[f64; 3]]) -> Result<()> {
    let flat: Vec<f64> = pts.iter().flatten().copied().collect();
    g.new_dataset::<f64>()
        .shape([pts.len(), 3])
        .create("Points")?
        .write_raw(&flat)?;
    Ok(())
}

/// Nodes of face `f` ordered with the normal out of cell `c`.
fn outward(mesh: &Mesh, f: usize, is_nb: bool) -> Vec<i64> {
    let n = mesh.face(f);
    if is_nb {
        let mut v: Vec<i64> = Vec::with_capacity(n.len());
        v.push(n[0] as i64);
        v.extend(n[1..].iter().rev().map(|&x| x as i64));
        v
    } else {
        n.iter().map(|&x| x as i64).collect()
    }
}

/// VTK node order of a standard cell, or `None` for a polyhedron.
fn standard_cell(mesh: &Mesh, c: usize) -> Option<(u8, Vec<i64>)> {
    let shape = mesh.cell_shape(c);
    if shape == CellShape::Polyhedron {
        return None;
    }
    let faces: Vec<Vec<i64>> = mesh
        .cell_face_iter(c)
        .map(|(f, nb)| outward(mesh, f, nb))
        .collect();
    // a cell whose faces do not close on the nodes of its shape (a face of an
    // imprinted interface, D182, has corners of its own) is a polyhedron
    let nodes: std::collections::HashSet<i64> = faces.iter().flatten().copied().collect();
    let expected = match shape {
        CellShape::Tetrahedron => 4,
        CellShape::Hexahedron => 8,
        CellShape::Prism => 6,
        _ => 5,
    };
    if nodes.len() != expected {
        return None;
    }
    // node adjacent to `n` along an edge, not in `base`
    let up = |n: i64, base: &[i64]| -> i64 {
        for f in &faces {
            let k = f.len();
            for i in 0..k {
                if f[i] == n {
                    for m in [f[(i + 1) % k], f[(i + k - 1) % k]] {
                        if !base.contains(&m) {
                            return m;
                        }
                    }
                }
            }
        }
        unreachable!("standard cell without a side edge")
    };
    let first_with = |len: usize| faces.iter().find(|f| f.len() == len).unwrap().clone();
    Some(match shape {
        CellShape::Tetrahedron => {
            let t = &faces[0];
            let apex = up(t[0], t);
            (VTK_TETRA, vec![t[0], t[2], t[1], apex])
        }
        CellShape::Hexahedron => {
            let q = &faces[0];
            let base = vec![q[0], q[3], q[2], q[1]];
            let mut v = base.clone();
            v.extend(base.iter().map(|&n| up(n, &base)));
            (VTK_HEXAHEDRON, v)
        }
        CellShape::Prism => {
            // base normal towards the opposite triangle, as for the
            // hexahedron: VTK's wedge lists (0, 2, 1) as its outward base
            // (vtkWedge, VTK 9.7), so (0, 1, 2) points towards (3, 4, 5).
            // Older VTK documentation said the opposite; the volumes VTK
            // computes and `tests::standard_cells_in_vtk_order` follow the
            // face table.
            let t = first_with(3);
            let base = vec![t[0], t[2], t[1]];
            let mut v = base.clone();
            v.extend(base.iter().map(|&n| up(n, &base)));
            (VTK_WEDGE, v)
        }
        CellShape::Pyramid => {
            let q = first_with(4);
            let base = vec![q[0], q[3], q[2], q[1]];
            let apex = up(base[0], &base);
            let mut v = base;
            v.push(apex);
            (VTK_PYRAMID, v)
        }
        CellShape::Polyhedron => unreachable!(),
    })
}

/// The 3-D Hilbert index of a point of the 2²¹ grid (63 bits), from John
/// Skilling, "Programming the Hilbert curve", AIP Conference Proceedings
/// 707, 381 (2004): his `AxestoTranspose`, then the bits of the transposed
/// index read from the top bit plane down, x, y, z within each plane.
fn hilbert(x: [u32; 3]) -> u64 {
    hilbert_bits(x, 21)
}

/// [`hilbert`] on the grid of `bits` bits per axis.
fn hilbert_bits(mut x: [u32; 3], bits: u32) -> u64 {
    let m = 1u32 << (bits - 1);
    // inverse undo
    let mut q = m;
    while q > 1 {
        let p = q - 1;
        for i in 0..3 {
            if x[i] & q != 0 {
                x[0] ^= p; // invert
            } else {
                let t = (x[0] ^ x[i]) & p; // exchange
                x[0] ^= t;
                x[i] ^= t;
            }
        }
        q >>= 1;
    }
    // Gray encode
    for i in 1..3 {
        x[i] ^= x[i - 1];
    }
    let mut t = 0;
    let mut q = m;
    while q > 1 {
        if x[2] & q != 0 {
            t ^= q - 1;
        }
        q >>= 1;
    }
    for v in &mut x {
        *v ^= t;
    }
    let mut h = 0u64;
    for b in (0..bits).rev() {
        for v in x {
            h = (h << 1) | u64::from((v >> b) & 1);
        }
    }
    h
}

/// What a cell writes: a linear cell's VTK type and nodes in VTK's order,
/// or a polyhedron's distinct nodes (in the order its faces first use
/// them) and its faces with their normals out of the cell.
enum Record {
    Linear(u8, Vec<i64>),
    Poly(Vec<i64>, Vec<Vec<i64>>),
}

impl Record {
    fn of(mesh: &Mesh, c: usize) -> Record {
        if let Some((t, v)) = standard_cell(mesh, c) {
            return Record::Linear(t, v);
        }
        let faces: Vec<Vec<i64>> = mesh
            .cell_face_iter(c)
            .map(|(f, nb)| outward(mesh, f, nb))
            .collect();
        let mut nodes: Vec<i64> = Vec::with_capacity(16);
        for &n in faces.iter().flatten() {
            if !nodes.contains(&n) {
                nodes.push(n);
            }
        }
        Record::Poly(nodes, faces)
    }

    /// Connectivity ids, faces and face connectivity ids of the cell.
    fn sizes(&self) -> [u32; 3] {
        match self {
            Record::Linear(_, v) => [v.len() as u32, 0, 0],
            Record::Poly(n, f) => [
                n.len() as u32,
                f.len() as u32,
                f.iter().map(Vec::len).sum::<usize>() as u32,
            ],
        }
    }
}

/// [`Record::sizes`] of every cell.
fn sizes(mesh: &Mesh) -> Vec<[u32; 3]> {
    use rayon::prelude::*;
    (0..mesh.n_cells)
        .into_par_iter()
        .map(|c| Record::of(mesh, c).sizes())
        .collect()
}

/// One partition of `internal`: its cells (the owned ones in file order,
/// then the ghosts in file order) and its points (by identity).
struct Part {
    cells: Vec<u32>,
    owned: usize,
    points: Vec<u32>,
    /// Minimum then maximum of the points of the owned cells.
    bounds: [f64; 6],
}

/// Rows of each partition.
#[derive(Default)]
struct Counts {
    points: Vec<usize>,
    cells: Vec<usize>,
    conn: Vec<usize>,
    faces: Vec<usize>,
    fconn: Vec<usize>,
}

/// The blocked layout of a mesh: the cells in Hilbert order of their
/// centroids, cut into runs of `size` owned cells, each run with one layer
/// of ghost cells (the cells of other blocks that share a point with it).
struct Layout {
    /// File cell id of each local cell.
    file_id: Vec<u32>,
    /// The block that owns each point: the lowest block whose owned cells
    /// use it.
    point_block: Vec<u32>,
    parts: Vec<Part>,
}

impl Layout {
    fn new(mesh: &Mesh, size: usize) -> Layout {
        use rayon::prelude::*;
        let n = mesh.n_cells;
        let np = mesh.n_points();
        // each axis of the bounding box scaled to 0 … 2²¹ − 1 on its own
        let (lo, hi) = mesh.bounding_box();
        let top = f64::from((1u32 << 21) - 1);
        let scale: [f64; 3] = std::array::from_fn(|k| {
            let d = hi[k] - lo[k];
            if d > 0.0 && d.is_finite() {
                top / d
            } else {
                0.0
            }
        });
        let key: Vec<u64> = (0..n)
            .into_par_iter()
            .map(|c| {
                // the centroid summed in the order of the point identities,
                // so that a renumbering does not move a cell
                let mut nodes = mesh.cell_nodes(c);
                nodes.sort_unstable_by_key(|&p| mesh.point_id[p as usize]);
                let mut s = [0.0; 3];
                for &p in &nodes {
                    for (k, x) in s.iter_mut().enumerate() {
                        *x += mesh.points[p as usize][k];
                    }
                }
                let m = nodes.len().max(1) as f64;
                hilbert(std::array::from_fn(|k| {
                    let v = ((s[k] / m - lo[k]) * scale[k]).floor();
                    if v >= 0.0 { v.min(top) as u32 } else { 0 }
                }))
            })
            .collect();
        // ties by the global identities
        let mut order: Vec<u32> = (0..n as u32).collect();
        order.par_sort_unstable_by_key(|&c| (key[c as usize], mesh.cell_id[c as usize]));
        drop(key);
        let mut file_id = vec![0u32; n];
        for (k, &c) in order.iter().enumerate() {
            file_id[c as usize] = k as u32;
        }
        // the cells around each point
        let mut ptr = vec![0usize; np + 1];
        for c in 0..n {
            for p in mesh.cell_nodes(c) {
                ptr[p as usize + 1] += 1;
            }
        }
        for p in 0..np {
            ptr[p + 1] += ptr[p];
        }
        let mut around = vec![0u32; ptr[np]];
        let mut next = ptr[..np].to_vec();
        for c in 0..n {
            for p in mesh.cell_nodes(c) {
                let p = p as usize;
                around[next[p]] = c as u32;
                next[p] += 1;
            }
        }
        drop(next);
        let point_block: Vec<u32> = (0..np)
            .map(|p| {
                around[ptr[p]..ptr[p + 1]]
                    .iter()
                    .map(|&c| file_id[c as usize] / size as u32)
                    .min()
                    .unwrap_or(u32::MAX)
            })
            .collect();
        // the partitions: owned cells, their point neighbours, the points
        let mut mark = vec![0u32; n];
        let mut seen = vec![0u32; np];
        let mut parts = Vec::with_capacity(n.div_ceil(size).max(1));
        for (b, own) in order.chunks(size).enumerate() {
            let tag = b as u32 + 1;
            for &c in own {
                mark[c as usize] = tag;
            }
            let mut points = Vec::new();
            let mut ghosts = Vec::new();
            let mut bounds = [f64::INFINITY; 6];
            bounds[3..].fill(f64::NEG_INFINITY);
            for &c in own {
                for p in mesh.cell_nodes(c as usize) {
                    let p = p as usize;
                    if seen[p] == tag {
                        continue;
                    }
                    seen[p] = tag;
                    points.push(p as u32);
                    for k in 0..3 {
                        bounds[k] = bounds[k].min(mesh.points[p][k]);
                        bounds[k + 3] = bounds[k + 3].max(mesh.points[p][k]);
                    }
                    for &d in &around[ptr[p]..ptr[p + 1]] {
                        if mark[d as usize] != tag {
                            mark[d as usize] = tag;
                            ghosts.push(d);
                        }
                    }
                }
            }
            ghosts.sort_unstable_by_key(|&c| file_id[c as usize]);
            for &g in &ghosts {
                for p in mesh.cell_nodes(g as usize) {
                    if seen[p as usize] != tag {
                        seen[p as usize] = tag;
                        points.push(p);
                    }
                }
            }
            points.sort_unstable_by_key(|&p| mesh.point_id[p as usize]);
            let mut cells = Vec::with_capacity(own.len() + ghosts.len());
            cells.extend_from_slice(own);
            cells.extend(ghosts);
            parts.push(Part {
                cells,
                owned: own.len(),
                points,
                bounds,
            });
        }
        if parts.is_empty() {
            parts.push(Part {
                cells: Vec::new(),
                owned: 0,
                points: Vec::new(),
                bounds: [f64::NAN; 6],
            });
        }
        Layout {
            file_id,
            point_block,
            parts,
        }
    }

    /// The rows of each partition, from the [`sizes`] of the cells.
    fn counts(&self, sizes: &[[u32; 3]]) -> Counts {
        let mut n = Counts::default();
        for part in &self.parts {
            let mut s = [0usize; 3];
            for &c in &part.cells {
                for (k, x) in s.iter_mut().enumerate() {
                    *x += sizes[c as usize][k] as usize;
                }
            }
            n.points.push(part.points.len());
            n.cells.push(part.cells.len());
            n.conn.push(s[0]);
            n.faces.push(s[1]);
            n.fconn.push(s[2]);
        }
        n
    }
}

/// A contiguous dataset of `n` rows, written block by block with [`put`].
fn dataset<T: h5::H5Type>(g: &h5::Group, name: &str, n: usize) -> Result<h5::Dataset> {
    Ok(g.new_dataset::<T>().shape([n]).create(name)?)
}

/// A contiguous dataset of `n` rows of `cols` values.
fn dataset2<T: h5::H5Type>(
    g: &h5::Group,
    name: &str,
    n: usize,
    cols: usize,
) -> Result<h5::Dataset> {
    Ok(g.new_dataset::<T>().shape([n, cols]).create(name)?)
}

/// Writes `v` to the rows from `at`.
fn put<T: h5::H5Type>(d: &h5::Dataset, at: usize, v: &[T]) -> Result<()> {
    if !v.is_empty() {
        d.write_slice(v, at..at + v.len())?;
    }
    Ok(())
}

/// Writes the rows of `cols` values in `v` from row `at`.
fn put2<T: h5::H5Type>(d: &h5::Dataset, at: usize, v: &[T], cols: usize) -> Result<()> {
    let rows = v.len() / cols;
    if rows > 0 {
        let a = ndarray::ArrayView2::from_shape((rows, cols), v)?;
        d.write_slice(a, (at..at + rows, ..))?;
    }
    Ok(())
}

/// Zone names as `id=name` entries in fixed-length strings of the longest
/// entry (null-padded), as KoflCraft writes them.
fn names_ds(g: &h5::Group, name: &str, entries: &[String]) -> Result<()> {
    use h5::types::{FixedUnicode, VarLenUnicode};
    let len = entries.iter().map(String::len).max().unwrap_or(1);
    macro_rules! fixed {
        ($($n:literal)*) => {
            $(if len <= $n {
                let v = entries
                    .iter()
                    .map(|e| e.parse::<FixedUnicode<$n>>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| anyhow::anyhow!("{e:?}"))?;
                return ds(g, name, &v);
            })*
        };
    }
    fixed!(8 16 32 64 128 256);
    let v = entries
        .iter()
        .map(|e| e.parse::<VarLenUnicode>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| anyhow::anyhow!("{e:?}"))?;
    ds(g, name, &v)
}

fn write_internal(g: &h5::Group, mesh: &Mesh, fields: &[CellField], layout: &Layout) -> Result<()> {
    use rayon::prelude::*;
    ascii_attr(g, "Type", "UnstructuredGrid")?;
    version_attr(g)?;
    let n = layout.counts(&sizes(mesh));
    let nb = layout.parts.len();
    let sum = |v: &[usize]| v.iter().sum::<usize>();
    let (tp, tc, tf) = (sum(&n.points), sum(&n.cells), sum(&n.faces));
    let i64s = |v: &[usize]| v.iter().map(|&x| x as i64).collect::<Vec<_>>();
    ds(g, "NumberOfPoints", &i64s(&n.points))?;
    ds(g, "NumberOfCells", &i64s(&n.cells))?;
    ds(g, "NumberOfConnectivityIds", &i64s(&n.conn))?;
    let points = dataset2::<f64>(g, "Points", tp, 3)?;
    let types = dataset::<u8>(g, "Types", tc)?;
    let offsets = dataset::<i64>(g, "Offsets", tc + nb)?;
    let conn = dataset::<i64>(g, "Connectivity", sum(&n.conn))?;
    ds(g, "NumberOfFaces", &i64s(&n.faces))?;
    ds(g, "NumberOfFaceConnectivityIds", &i64s(&n.fconn))?;
    ds(g, "NumberOfPolyhedronToFaceIds", &i64s(&n.faces))?;
    let face_offsets = dataset::<i64>(g, "FaceOffsets", tf + nb)?;
    let face_conn = dataset::<i64>(g, "FaceConnectivity", sum(&n.fconn))?;
    let poly_offsets = dataset::<i64>(g, "PolyhedronOffsets", tc + nb)?;
    let poly_faces = dataset::<i64>(g, "PolyhedronToFaces", tf)?;
    let cd = group(g, "CellData")?;
    let values = fields
        .iter()
        .map(|f| dataset::<f64>(&cd, f.name, tc))
        .collect::<Result<Vec<_>>>()?;
    let zone = if fields.iter().any(|f| f.name == "zone") {
        None
    } else {
        Some(dataset::<u32>(&cd, "zone", tc)?)
    };
    let cell_ids = dataset::<i64>(&cd, "GlobalCellIds", tc)?;
    let original = dataset::<i64>(&cd, "OriginalCellIds", tc)?;
    let cell_ghost = dataset::<u8>(&cd, "vtkGhostType", tc)?;
    let pd = group(g, "PointData")?;
    let point_ids = dataset::<i64>(&pd, "GlobalPointIds", tp)?;
    let point_ghost = dataset::<u8>(&pd, "vtkGhostType", tp)?;

    let mut ranges = vec![Vec::with_capacity(2 * nb); fields.len()];
    let mut local = vec![0u32; mesh.n_points()];
    let (mut at_p, mut at_c, mut at_k, mut at_f, mut at_fk) = (0, 0, 0, 0, 0);
    for (b, part) in layout.parts.iter().enumerate() {
        let b32 = b as u32;
        for (k, &p) in part.points.iter().enumerate() {
            local[p as usize] = k as u32;
        }
        let xyz: Vec<f64> = part
            .points
            .iter()
            .flat_map(|&p| mesh.points[p as usize])
            .collect();
        put2(&points, at_p, &xyz, 3)?;
        drop(xyz);
        let ids: Vec<i64> = part
            .points
            .iter()
            .map(|&p| mesh.point_id[p as usize] as i64)
            .collect();
        put(&point_ids, at_p, &ids)?;
        let ghost: Vec<u8> = part
            .points
            .iter()
            .map(|&p| u8::from(layout.point_block[p as usize] != b32))
            .collect();
        put(&point_ghost, at_p, &ghost)?;

        let records: Vec<Record> = part
            .cells
            .par_iter()
            .map(|&c| Record::of(mesh, c as usize))
            .collect();
        let nc = part.cells.len();
        let l = |v: &i64| local[*v as usize] as i64;
        let mut t = Vec::with_capacity(nc);
        let mut off = Vec::with_capacity(nc + 1);
        let mut cn = Vec::with_capacity(n.conn[b]);
        let mut poff = Vec::with_capacity(nc + 1);
        let mut pf = Vec::with_capacity(n.faces[b]);
        let mut foff = Vec::with_capacity(n.faces[b] + 1);
        let mut fc = Vec::with_capacity(n.fconn[b]);
        off.push(0i64);
        poff.push(0i64);
        foff.push(0i64);
        for r in &records {
            match r {
                Record::Linear(ty, v) => {
                    t.push(*ty);
                    cn.extend(v.iter().map(l));
                }
                Record::Poly(v, faces) => {
                    t.push(VTK_POLYHEDRON);
                    cn.extend(v.iter().map(l));
                    for f in faces {
                        pf.push(foff.len() as i64 - 1);
                        fc.extend(f.iter().map(l));
                        foff.push(fc.len() as i64);
                    }
                }
            }
            off.push(cn.len() as i64);
            poff.push(pf.len() as i64);
        }
        drop(records);
        debug_assert_eq!(cn.len(), n.conn[b]);
        debug_assert_eq!(fc.len(), n.fconn[b]);
        put(&types, at_c, &t)?;
        put(&offsets, at_c + b, &off)?;
        put(&conn, at_k, &cn)?;
        put(&poly_offsets, at_c + b, &poff)?;
        put(&poly_faces, at_f, &pf)?;
        put(&face_offsets, at_f + b, &foff)?;
        put(&face_conn, at_fk, &fc)?;

        for ((f, d), range) in fields.iter().zip(&values).zip(&mut ranges) {
            let v: Vec<f64> = part.cells.iter().map(|&c| f.values[c as usize]).collect();
            put(d, at_c, &v)?;
            let (lo, hi) = v[..part.owned]
                .iter()
                .filter(|x| x.is_finite())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &x| {
                    (lo.min(x), hi.max(x))
                });
            range.extend(if lo <= hi { [lo, hi] } else { [f64::NAN; 2] });
        }
        if let Some(d) = &zone {
            let v: Vec<u32> = part
                .cells
                .iter()
                .map(|&c| mesh.cell_zone[c as usize])
                .collect();
            put(d, at_c, &v)?;
        }
        let v: Vec<i64> = part
            .cells
            .iter()
            .map(|&c| layout.file_id[c as usize] as i64)
            .collect();
        put(&cell_ids, at_c, &v)?;
        let v: Vec<i64> = part
            .cells
            .iter()
            .map(|&c| mesh.cell_id[c as usize] as i64)
            .collect();
        put(&original, at_c, &v)?;
        let v: Vec<u8> = (0..nc).map(|k| u8::from(k >= part.owned)).collect();
        put(&cell_ghost, at_c, &v)?;

        at_p += n.points[b];
        at_c += nc;
        at_k += n.conn[b];
        at_f += n.faces[b];
        at_fk += n.fconn[b];
    }

    let fd = group(g, "FieldData")?;
    let owned: Vec<i64> = layout.parts.iter().map(|p| p.owned as i64).collect();
    ds(&fd, "BlockOwnedCells", &owned)?;
    let bounds: Vec<f64> = layout.parts.iter().flat_map(|p| p.bounds).collect();
    put2(&dataset2::<f64>(&fd, "BlockBounds", nb, 6)?, 0, &bounds, 6)?;
    for (f, r) in fields.iter().zip(&ranges) {
        let d = dataset2::<f64>(&fd, &format!("BlockRange_{}", f.name), nb, 2)?;
        put2(&d, 0, r, 2)?;
    }
    if zone.is_some() && !mesh.zones.is_empty() {
        let names: Vec<String> = mesh
            .zones
            .iter()
            .enumerate()
            .map(|(i, z)| format!("{i}={z}"))
            .collect();
        names_ds(&fd, "zone_names", &names)?;
    }
    Ok(())
}

fn write_boundary(
    g: &h5::Group,
    mesh: &Mesh,
    faces: std::ops::Range<usize>,
    boundary: f64,
    face_fields: &[FaceField],
    file_id: &[u32],
) -> Result<()> {
    ascii_attr(g, "Type", "PolyData")?;
    version_attr(g)?;
    let mut map = std::collections::HashMap::new();
    let mut pts: Vec<[f64; 3]> = Vec::new();
    let mut ids: Vec<i64> = Vec::new();
    let mut conn: Vec<i64> = Vec::new();
    let mut offsets: Vec<i64> = vec![0];
    let mut order: Vec<usize> = faces.clone().collect();
    order.sort_by_key(|&f| mesh.face_id[f]);
    for &f in &order {
        for &n in mesh.face(f) {
            let id = *map.entry(n).or_insert_with(|| {
                pts.push(mesh.points[n as usize]);
                ids.push(mesh.point_id[n as usize] as i64);
                pts.len() as i64 - 1
            });
            conn.push(id);
        }
        offsets.push(conn.len() as i64);
    }
    ds(g, "NumberOfPoints", &[pts.len() as i64])?;
    points_ds(g, &pts)?;
    for (name, n, o, c) in [
        ("Vertices", 0, vec![0i64], vec![]),
        ("Lines", 0, vec![0], vec![]),
        ("Strips", 0, vec![0], vec![]),
        ("Polygons", faces.len(), offsets, conn),
    ] {
        let sg = group(g, name)?;
        ds(&sg, "NumberOfCells", &[n as i64])?;
        ds(&sg, "NumberOfConnectivityIds", &[c.len() as i64])?;
        ds(&sg, "Offsets", &o)?;
        ds(&sg, "Connectivity", &c)?;
    }
    let cd = group(g, "CellData")?;
    ds(&cd, "boundary", &vec![boundary; faces.len()])?;
    let owner: Vec<i64> = order
        .iter()
        .map(|&f| file_id[mesh.owner[f] as usize] as i64)
        .collect();
    ds(&cd, "OwnerCell", &owner)?;
    let ni = mesh.n_interior_faces();
    for ff in face_fields {
        let v: Vec<f64> = order.iter().map(|&f| ff.values[f - ni]).collect();
        ds(&cd, ff.name, &v)?;
    }
    let pd = group(g, "PointData")?;
    ds(&pd, "GlobalPointIds", &ids)?;
    group(g, "FieldData")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::{CellKind, ElementMesh};
    use crate::generate::{BoxSpec, Kind};

    /// The faces of VTK's linear cells, each with its normal out of the
    /// cell, as `GetFace` lists them (vtkTetra, vtkHexahedron, vtkWedge,
    /// vtkPyramid; VTK 9.7).
    fn vtk_faces(t: u8) -> &'static [&'static [usize]] {
        match t {
            VTK_TETRA => &[&[0, 1, 3], &[1, 2, 3], &[2, 0, 3], &[0, 2, 1]],
            VTK_HEXAHEDRON => &[
                &[0, 4, 7, 3],
                &[1, 2, 6, 5],
                &[0, 1, 5, 4],
                &[3, 7, 6, 2],
                &[0, 3, 2, 1],
                &[4, 5, 6, 7],
            ],
            VTK_WEDGE => &[
                &[0, 2, 1],
                &[3, 4, 5],
                &[0, 1, 4, 3],
                &[1, 2, 5, 4],
                &[2, 0, 3, 5],
            ],
            VTK_PYRAMID => &[
                &[0, 3, 2, 1],
                &[0, 1, 4],
                &[1, 2, 4],
                &[2, 3, 4],
                &[3, 0, 4],
            ],
            _ => unreachable!(),
        }
    }

    /// `a` and `b` are the same cycle of nodes (same orientation).
    fn same_cycle(a: &[i64], b: &[i64]) -> bool {
        a.len() == b.len()
            && (0..b.len()).any(|r| (0..a.len()).all(|k| a[k] == b[(k + r) % b.len()]))
    }

    /// Volume enclosed by faces, by the divergence theorem over triangles
    /// fanned from each face's mean point; positive when they point out.
    fn volume(points: &[[f64; 3]], faces: &[Vec<i64>]) -> f64 {
        let mut v = 0.0;
        for f in faces {
            let p: Vec<[f64; 3]> = f.iter().map(|&n| points[n as usize]).collect();
            let mut m = [0.0; 3];
            for q in &p {
                for d in 0..3 {
                    m[d] += q[d] / p.len() as f64;
                }
            }
            for k in 0..p.len() {
                let (a, b) = (p[k], p[(k + 1) % p.len()]);
                v += a[0] * (b[1] * m[2] - b[2] * m[1]) - a[1] * (b[0] * m[2] - b[2] * m[0])
                    + a[2] * (b[0] * m[1] - b[1] * m[0]);
            }
        }
        v / 6.0
    }

    fn one_cell(kind: CellKind, points: &[[f64; 3]], nodes: &[u32]) -> Mesh {
        let mut em = ElementMesh {
            points: points.to_vec(),
            ..Default::default()
        };
        em.add_cell(kind, nodes, 0);
        em.into_mesh().unwrap().mesh
    }

    /// Invariant: every tetrahedron, hexahedron, prism and pyramid is
    /// written in VTK's node order, so VTK's own faces of the cell are the
    /// cell's faces, pointing out, and VTK finds a positive volume. Cells
    /// given inside out and perturbed generated meshes (faces seen from
    /// both sides) included.
    #[test]
    fn standard_cells_in_vtk_order() {
        let cube = [
            [0.0, 0.0, 0.0],
            [1.1, 0.1, 0.0],
            [1.0, 0.9, 0.1],
            [-0.1, 1.0, 0.0],
            [0.1, 0.0, 1.0],
            [1.0, 0.1, 1.2],
            [1.1, 1.0, 0.9],
            [0.0, 1.1, 1.0],
        ];
        let apex = [0.4, 0.6, 1.3];
        let mut meshes = Vec::new();
        for flip in [false, true] {
            let o = |v: &[u32]| -> Vec<u32> {
                let mut v = v.to_vec();
                if flip {
                    // the same cell given with its base the other way round
                    let h = v.len() / 2;
                    match v.len() {
                        4 => v.swap(1, 2),
                        5 => v.swap(1, 3),
                        _ => {
                            v[..h].reverse();
                            v[h..].reverse();
                        }
                    }
                }
                v
            };
            meshes.push(one_cell(
                CellKind::Hexahedron,
                &cube,
                &o(&[0, 1, 2, 3, 4, 5, 6, 7]),
            ));
            meshes.push(one_cell(CellKind::Prism, &cube, &o(&[0, 1, 3, 4, 5, 7])));
            meshes.push(one_cell(CellKind::Tetrahedron, &cube, &o(&[0, 1, 3, 4])));
            let mut p = cube[..4].to_vec();
            p.push(apex);
            meshes.push(one_cell(CellKind::Pyramid, &p, &o(&[0, 1, 2, 3, 4])));
        }
        for (kind, dim) in [
            (Kind::Hex, 3),
            (Kind::Tet, 2),
            (Kind::Tet, 3),
            (Kind::Mixed, 3),
        ] {
            let mut s = BoxSpec::new(kind, dim, [3, 3, 2]);
            s.perturb = 0.2;
            meshes.push(s.generate().unwrap());
        }
        let mut seen = std::collections::BTreeSet::new();
        for mesh in &meshes {
            for c in 0..mesh.n_cells {
                let Some((t, v)) = standard_cell(mesh, c) else {
                    continue;
                };
                seen.insert(t);
                let own: Vec<Vec<i64>> = mesh
                    .cell_face_iter(c)
                    .map(|(f, nb)| outward(mesh, f, nb))
                    .collect();
                let vtk: Vec<Vec<i64>> = vtk_faces(t)
                    .iter()
                    .map(|f| f.iter().map(|&k| v[k]).collect())
                    .collect();
                assert_eq!(vtk.len(), own.len(), "type {t}, cell {c}");
                for f in &vtk {
                    assert!(
                        own.iter().any(|g| same_cycle(f, g)),
                        "type {t}, cell {c}: VTK face {f:?} is not an outward face {own:?}"
                    );
                }
                let (a, b) = (volume(&mesh.points, &vtk), volume(&mesh.points, &own));
                assert!(
                    a > 0.0 && (a - b).abs() <= 1e-12 * b,
                    "type {t}, cell {c}: {a} {b}"
                );
            }
        }
        let all = [VTK_TETRA, VTK_HEXAHEDRON, VTK_WEDGE, VTK_PYRAMID];
        assert_eq!(seen.into_iter().collect::<Vec<_>>(), all);
    }

    /// Invariant: the Hilbert index visits every point of a 2^b grid once
    /// (a bijection onto 0 … 8^b − 1), each step to a neighbour along one
    /// axis, which is what keeps a block compact; the file's 21-bit index
    /// starts at the origin and stays within 63 bits.
    #[test]
    fn hilbert_index_walks_the_grid_in_unit_steps() {
        for bits in 1..=4u32 {
            let n = 1u32 << bits;
            let mut at = vec![None; (n * n * n) as usize];
            for x in 0..n {
                for y in 0..n {
                    for z in 0..n {
                        let h = hilbert_bits([x, y, z], bits) as usize;
                        assert!(at[h].is_none(), "index {h} twice at {bits} bits");
                        at[h] = Some([x, y, z]);
                    }
                }
            }
            for w in at.windows(2) {
                let (a, b) = (w[0].unwrap(), w[1].unwrap());
                let d: u32 = (0..3).map(|k| a[k].abs_diff(b[k])).sum();
                assert_eq!(d, 1, "{a:?} then {b:?} at {bits} bits");
            }
        }
        let top = (1u32 << 21) - 1;
        assert_eq!(hilbert([0, 0, 0]), 0);
        assert!(hilbert([top, top, top]) < 1 << 63);
    }
}
