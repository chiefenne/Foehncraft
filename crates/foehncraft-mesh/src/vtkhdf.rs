//! VTKHDF output (D41).
//!
//! One HDF5 file in the VTKHDF format of VTK, read directly by ParaView
//! (open the `.vtkhdf` file; reader "VTKHDF"). The file is a
//! `PartitionedDataSetCollection`:
//!
//! - block `internal`: the cells as an `UnstructuredGrid`; tetrahedra,
//!   hexahedra, prisms and pyramids as VTK's linear cells, everything else
//!   as `VTK_POLYHEDRON` with its faces; cell fields under `CellData`;
//! - one `PolyData` block per boundary boundary, named after the boundary, with
//!   the boundary index and the boundary face fields ([`FaceField`], such as
//!   the wall shear stress) as cell fields;
//! - an assembly `internal`, `boundary/<boundary>` for ParaView's block
//!   selection.
//!
//! Cells, nodes and boundary faces are written in the order of their global
//! identities (input order, R2), so the file does not depend on a
//! renumbering of the cells.
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

/// As [`write_with`], with fields on the boundary faces, which go to the
/// boundary blocks.
pub fn write_all(
    path: impl AsRef<Path>,
    mesh: &Mesh,
    fields: &[CellField],
    face_fields: &[FaceField],
    attrs: &[(&str, Attr)],
) -> Result<()> {
    let path = path.as_ref();
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
    type_attr(&root, "PartitionedDataSetCollection")?;
    version_attr(&root)?;

    let internal = group(&root, "internal")?;
    write_internal(&internal, mesh, fields)?;
    index_attr(&internal, 0)?;

    let assembly = group(&root, "Assembly")?;
    assembly.link_soft("/VTKHDF/internal", "internal")?;
    let boundary = group(&assembly, "boundary")?;
    for (i, p) in mesh.boundaries.iter().enumerate() {
        let name = block_name(&p.name);
        let g = group(&root, &name)?;
        write_boundary(&g, mesh, p.range(), i as f64, face_fields)?;
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
    type_attr(&g, "PolyData")?;
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
/// cent: HDF5 adds headers of some kilobytes).
pub fn estimate_bytes(mesh: &Mesh, n_fields: usize) -> f64 {
    let nc = mesh.n_cells;
    let mut ids = 0usize;
    for c in 0..nc {
        match mesh.cell_shape(c) {
            CellShape::Tetrahedron => ids += 4,
            CellShape::Hexahedron => ids += 8,
            CellShape::Prism => ids += 6,
            CellShape::Pyramid => ids += 5,
            CellShape::Polyhedron => {
                ids += mesh.cell_nodes(c).len();
                for (f, _) in mesh.cell_face_iter(c) {
                    // face connectivity, face offset, polyhedron-to-face id
                    ids += (mesh.face_ptr[f + 1] - mesh.face_ptr[f]) as usize + 2;
                }
            }
        }
    }
    // boundary blocks: points of each boundary and the face connectivity
    let mut boundary = 0usize;
    for p in &mesh.boundaries {
        let mut points = std::collections::HashSet::new();
        let mut nodes = 0usize;
        for f in p.range() {
            let fnodes = &mesh.face_nodes[mesh.face_ptr[f] as usize..mesh.face_ptr[f + 1] as usize];
            nodes += fnodes.len();
            points.extend(fnodes.iter().copied());
        }
        boundary += points.len() * 24 + nodes * 8 + p.len * (8 + 1 + 8 + 8);
    }
    (mesh.n_points() * 24 + nc * (8 * (n_fields + 1) + 8 + 1 + 8 + 8) + ids * 8 + boundary) as f64
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

fn type_attr(g: &h5::Group, t: &str) -> Result<()> {
    macro_rules! fixed {
        ($n:literal) => {{
            let v = FixedAscii::<$n>::from_ascii(t.as_bytes())?;
            g.new_attr::<FixedAscii<$n>>()
                .create("Type")?
                .write_scalar(&v)?;
        }};
    }
    match t.len() {
        8 => fixed!(8),
        16 => fixed!(16),
        28 => fixed!(28),
        _ => anyhow::bail!("unexpected type name {t}"),
    }
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

/// The nodes in the order of their identities, and the output position of
/// each local node.
fn points_in_id_order(mesh: &Mesh) -> (Vec<[f64; 3]>, Vec<u32>) {
    let mut order: Vec<usize> = (0..mesh.n_points()).collect();
    order.sort_by_key(|&v| mesh.point_id[v]);
    let mut pos = vec![0u32; mesh.n_points()];
    for (k, &v) in order.iter().enumerate() {
        pos[v] = k as u32;
    }
    (order.iter().map(|&v| mesh.points[v]).collect(), pos)
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
            // base normal towards the opposite triangle, as for the hexahedron
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

fn write_internal(g: &h5::Group, mesh: &Mesh, fields: &[CellField]) -> Result<()> {
    type_attr(g, "UnstructuredGrid")?;
    version_attr(g)?;
    let nc = mesh.n_cells;
    let (pts, pos) = points_in_id_order(mesh);
    let mut cells: Vec<usize> = (0..nc).collect();
    cells.sort_by_key(|&c| mesh.cell_id[c]);
    let mut types = Vec::with_capacity(nc);
    let mut offsets: Vec<i64> = Vec::with_capacity(nc + 1);
    let mut conn: Vec<i64> = Vec::with_capacity(8 * nc);
    let mut poly_offsets: Vec<i64> = Vec::with_capacity(nc + 1);
    let mut poly_faces: Vec<i64> = Vec::new();
    let mut face_offsets: Vec<i64> = vec![0];
    let mut face_conn: Vec<i64> = Vec::new();
    offsets.push(0);
    poly_offsets.push(0);
    let p = |n: i64| pos[n as usize] as i64;
    for &c in &cells {
        match standard_cell(mesh, c) {
            Some((t, v)) => {
                types.push(t);
                conn.extend(v.into_iter().map(p));
            }
            None => {
                types.push(VTK_POLYHEDRON);
                conn.extend(
                    mesh.cell_nodes(c)
                        .into_iter()
                        .map(|n| pos[n as usize] as i64),
                );
                for (f, nb) in mesh.cell_face_iter(c) {
                    poly_faces.push(face_offsets.len() as i64 - 1);
                    face_conn.extend(outward(mesh, f, nb).into_iter().map(p));
                    face_offsets.push(face_conn.len() as i64);
                }
            }
        }
        offsets.push(conn.len() as i64);
        poly_offsets.push(poly_faces.len() as i64);
    }
    ds(g, "NumberOfPoints", &[mesh.n_points() as i64])?;
    ds(g, "NumberOfCells", &[nc as i64])?;
    ds(g, "NumberOfConnectivityIds", &[conn.len() as i64])?;
    points_ds(g, &pts)?;
    ds(g, "Types", &types)?;
    ds(g, "Offsets", &offsets)?;
    ds(g, "Connectivity", &conn)?;
    ds(g, "NumberOfFaces", &[face_offsets.len() as i64 - 1])?;
    ds(g, "NumberOfFaceConnectivityIds", &[face_conn.len() as i64])?;
    ds(g, "NumberOfPolyhedronToFaceIds", &[poly_faces.len() as i64])?;
    ds(g, "FaceOffsets", &face_offsets)?;
    ds(g, "FaceConnectivity", &face_conn)?;
    ds(g, "PolyhedronOffsets", &poly_offsets)?;
    ds(g, "PolyhedronToFaces", &poly_faces)?;
    let cd = group(g, "CellData")?;
    for f in fields {
        let v: Vec<f64> = cells.iter().map(|&c| f.values[c]).collect();
        ds(&cd, f.name, &v)?;
    }
    let zone: Vec<f64> = cells.iter().map(|&c| mesh.cell_zone[c] as f64).collect();
    if !fields.iter().any(|f| f.name == "zone") {
        ds(&cd, "zone", &zone)?;
    }
    group(g, "PointData")?;
    group(g, "FieldData")?;
    Ok(())
}

fn write_boundary(
    g: &h5::Group,
    mesh: &Mesh,
    faces: std::ops::Range<usize>,
    boundary: f64,
    face_fields: &[FaceField],
) -> Result<()> {
    type_attr(g, "PolyData")?;
    version_attr(g)?;
    let mut map = std::collections::HashMap::new();
    let mut pts: Vec<[f64; 3]> = Vec::new();
    let mut conn: Vec<i64> = Vec::new();
    let mut offsets: Vec<i64> = vec![0];
    let mut order: Vec<usize> = faces.clone().collect();
    order.sort_by_key(|&f| mesh.face_id[f]);
    for &f in &order {
        for &n in mesh.face(f) {
            let id = *map.entry(n).or_insert_with(|| {
                pts.push(mesh.points[n as usize]);
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
    let ni = mesh.n_interior_faces();
    for ff in face_fields {
        let v: Vec<f64> = order.iter().map(|&f| ff.values[f - ni]).collect();
        ds(&cd, ff.name, &v)?;
    }
    group(g, "PointData")?;
    group(g, "FieldData")?;
    Ok(())
}
