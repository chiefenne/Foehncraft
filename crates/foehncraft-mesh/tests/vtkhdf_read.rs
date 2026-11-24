//! VTKHDF meshes (`format::vtkhdf`): blocked files as KoflCraft hands them
//! over (polyhedra in partitions with a ghost layer, `GlobalPointIds`,
//! `OriginalCellIds`, patches as `PolyData` blocks) read as the same mesh
//! as the `.vtu` of the same cells, and FoehnCraft's own `fields.vtkhdf`
//! reads back as its mesh.
//!
//! The blocked files are written here by a small writer of the layout,
//! independent of the reader, and by KoflCraft (the fixtures of
//! `tests/data/vtkhdf`, written by KoflCraft PR #38 with the `.vtu` of the
//! same export and its `mesh-evidence.json`).
#![cfg(feature = "hdf5")]

use foehncraft_mesh::format::vtk::write_vtu;
use foehncraft_mesh::format::{self, MeshFile};
use foehncraft_mesh::quality::QualityReport;
use foehncraft_mesh::vtkhdf::{self, CellField};
use foehncraft_mesh::{Geometry, Mesh};
use h5::types::{FixedAscii, VarLenUnicode};
use hdf5_metno as h5;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

fn data(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

fn dir() -> PathBuf {
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join("vtkhdf_read");
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn same(a: &Mesh, b: &Mesh, what: &str) {
    assert_eq!(a.n_cells, b.n_cells, "{what}: cells");
    assert_eq!(a.points, b.points, "{what}: points");
    assert_eq!(a.face_ptr, b.face_ptr, "{what}: face_ptr");
    assert_eq!(a.face_nodes, b.face_nodes, "{what}: face_nodes");
    assert_eq!(a.owner, b.owner, "{what}: owner");
    assert_eq!(a.neighbour, b.neighbour, "{what}: neighbour");
    assert_eq!(a.cell_zone, b.cell_zone, "{what}: cell_zone");
    assert_eq!(a.zones, b.zones, "{what}: zones");
    let names = |m: &Mesh| -> Vec<(String, usize, usize)> {
        m.boundaries
            .iter()
            .map(|b| (b.name.clone(), b.start, b.len))
            .collect()
    };
    assert_eq!(names(a), names(b), "{what}: boundaries");
}

/// The faces of each cell as `write_vtu` lists them (face order, the
/// neighbour's copy reversed).
fn cell_faces(m: &Mesh) -> Vec<Vec<Vec<u32>>> {
    let mut out = vec![Vec::new(); m.n_cells];
    for f in 0..m.n_faces() {
        let n = m.face(f).to_vec();
        if f < m.n_interior_faces() {
            let mut r = n.clone();
            r.reverse();
            out[m.neighbour[f] as usize].push(r);
        }
        out[m.owner[f] as usize].push(n);
    }
    out
}

fn ds<T: h5::H5Type>(g: &h5::Group, name: &str, v: &[T]) {
    g.new_dataset::<T>()
        .shape([v.len()])
        .create(name)
        .unwrap()
        .write_raw(v)
        .unwrap();
}

fn text_attr(g: &h5::Group, name: &str, v: &str) {
    let s = FixedAscii::<64>::from_ascii(v.as_bytes()).unwrap();
    g.new_attr::<FixedAscii<64>>()
        .create(name)
        .unwrap()
        .write_scalar(&s)
        .unwrap();
}

fn int_attr(g: &h5::Group, name: &str, v: i64) {
    g.new_attr::<i64>()
        .create(name)
        .unwrap()
        .write_scalar(&v)
        .unwrap();
}

fn points(g: &h5::Group, p: &[[f64; 3]]) {
    let flat: Vec<f64> = p.iter().flatten().copied().collect();
    g.new_dataset::<f64>()
        .shape([p.len(), 3])
        .create("Points")
        .unwrap()
        .write_raw(&flat)
        .unwrap();
}

/// How a blocked file is written.
struct Layout {
    /// Owned cells per block.
    block: usize,
    /// The mesh cells in file order.
    order: Vec<usize>,
    /// Write `GlobalPointIds` (cells and patches).
    ids: bool,
    /// Write `OriginalCellIds`.
    original: bool,
}

/// Writes `m` in the blocked VTKHDF layout: every cell a polyhedron with
/// its faces wound outward, partitions of `block` owned cells plus their
/// point neighbours as ghosts, one `PolyData` block per boundary.
fn write_blocked(path: &Path, m: &Mesh, l: &Layout) {
    let faces = cell_faces(m);
    let cell_points: Vec<Vec<u32>> = faces
        .iter()
        .map(|fs| {
            let mut v: Vec<u32> = fs.iter().flatten().copied().collect();
            v.sort_unstable();
            v.dedup();
            v
        })
        .collect();
    let mut point_cells = vec![Vec::new(); m.points.len()];
    for (c, ps) in cell_points.iter().enumerate() {
        for &p in ps {
            point_cells[p as usize].push(c);
        }
    }
    let blocks: Vec<&[usize]> = l.order.chunks(l.block).collect();
    let mut block_of = vec![0usize; m.n_cells];
    let mut file_id = vec![0i64; m.n_cells];
    for (b, cells) in blocks.iter().enumerate() {
        for &c in *cells {
            block_of[c] = b;
        }
    }
    for (k, &c) in l.order.iter().enumerate() {
        file_id[c] = k as i64;
    }
    let mut point_owner = vec![usize::MAX; m.points.len()];
    for (b, cells) in blocks.iter().enumerate() {
        for &c in *cells {
            for &p in &cell_points[c] {
                point_owner[p as usize] = point_owner[p as usize].min(b);
            }
        }
    }

    let file = h5::File::create(path).unwrap();
    let root = file.create_group("VTKHDF").unwrap();
    root.new_attr::<i64>()
        .shape([2])
        .create("Version")
        .unwrap()
        .write_raw(&[2i64, 4])
        .unwrap();
    text_attr(&root, "Type", "PartitionedDataSetCollection");
    int_attr(&root, "murmelcraft_blocks", 1);
    let v: VarLenUnicode = "0.9.0".parse().unwrap();
    root.new_attr::<VarLenUnicode>()
        .create("koflcraft_version")
        .unwrap()
        .write_scalar(&v)
        .unwrap();
    root.create_group("Assembly").unwrap();
    let g = root.create_group("internal").unwrap();
    text_attr(&g, "Type", "UnstructuredGrid");
    int_attr(&g, "Index", 0);
    int_attr(&g, "block_cells", l.block as i64);

    let (mut n_points, mut n_cells, mut n_conn) = (vec![], vec![], vec![]);
    let (mut n_faces, mut n_p2f, mut n_fconn) = (vec![], vec![], vec![]);
    let (mut pts, mut gpid, mut pghost) = (vec![], vec![], vec![]);
    let (mut types, mut offsets, mut conn) = (vec![], vec![], vec![]);
    let (mut foff, mut fconn, mut poff, mut p2f) = (vec![], vec![], vec![], vec![]);
    let (mut zone, mut gcid, mut ocid, mut cghost) = (vec![], vec![], vec![], vec![]);
    let mut owned_counts = vec![];
    for (b, owned) in blocks.iter().enumerate() {
        let mut ghosts = BTreeSet::new();
        for &c in *owned {
            for &p in &cell_points[c] {
                for &d in &point_cells[p as usize] {
                    if block_of[d] != b {
                        ghosts.insert((file_id[d], d));
                    }
                }
            }
        }
        let cells: Vec<(usize, u8)> = owned
            .iter()
            .map(|&c| (c, 0u8))
            .chain(ghosts.iter().map(|&(_, d)| (d, 1u8)))
            .collect();
        let mut local = std::collections::HashMap::new();
        let mut order = Vec::new();
        for &(c, _) in &cells {
            for f in &faces[c] {
                for &p in f {
                    local.entry(p).or_insert_with(|| {
                        order.push(p);
                        order.len() as i64 - 1
                    });
                }
            }
        }
        n_points.push(order.len() as i64);
        n_cells.push(cells.len() as i64);
        owned_counts.push(owned.len() as i64);
        for &p in &order {
            pts.push(m.points[p as usize]);
            gpid.push(p as i64);
            pghost.push(u8::from(point_owner[p as usize] != b));
        }
        let (k0, f0, q0, t0) = (conn.len(), foff.len(), fconn.len(), p2f.len());
        offsets.push(0i64);
        poff.push(0i64);
        foff.push(0i64);
        let mut nf = 0i64;
        for &(c, ghost) in &cells {
            types.push(42u8);
            conn.extend(cell_points[c].iter().map(|p| local[p]));
            offsets.push((conn.len() - k0) as i64);
            for f in &faces[c] {
                p2f.push(nf);
                nf += 1;
                fconn.extend(f.iter().map(|p| local[p]));
                foff.push((fconn.len() - q0) as i64);
            }
            poff.push((p2f.len() - t0) as i64);
            zone.push(m.cell_zone[c]);
            gcid.push(file_id[c]);
            ocid.push(c as i64);
            cghost.push(ghost);
        }
        n_conn.push((conn.len() - k0) as i64);
        n_faces.push(nf);
        n_p2f.push((p2f.len() - t0) as i64);
        n_fconn.push((fconn.len() - q0) as i64);
        let _ = f0;
    }
    ds(&g, "NumberOfPoints", &n_points);
    ds(&g, "NumberOfCells", &n_cells);
    ds(&g, "NumberOfConnectivityIds", &n_conn);
    ds(&g, "NumberOfFaces", &n_faces);
    ds(&g, "NumberOfPolyhedronToFaceIds", &n_p2f);
    ds(&g, "NumberOfFaceConnectivityIds", &n_fconn);
    points(&g, &pts);
    ds(&g, "Types", &types);
    ds(&g, "Offsets", &offsets);
    ds(&g, "Connectivity", &conn);
    ds(&g, "FaceOffsets", &foff);
    ds(&g, "FaceConnectivity", &fconn);
    ds(&g, "PolyhedronOffsets", &poff);
    ds(&g, "PolyhedronToFaces", &p2f);
    let cd = g.create_group("CellData").unwrap();
    ds(&cd, "zone", &zone);
    ds(&cd, "GlobalCellIds", &gcid);
    if l.original {
        ds(&cd, "OriginalCellIds", &ocid);
    }
    ds(&cd, "vtkGhostType", &cghost);
    let pd = g.create_group("PointData").unwrap();
    if l.ids {
        ds(&pd, "GlobalPointIds", &gpid);
    }
    ds(&pd, "vtkGhostType", &pghost);
    let fd = g.create_group("FieldData").unwrap();
    let names: Vec<FixedAscii<64>> = m
        .zones
        .iter()
        .enumerate()
        .map(|(i, z)| FixedAscii::from_ascii(format!("{i}={z}").as_bytes()).unwrap())
        .collect();
    ds(&fd, "zone_names", &names);
    ds(&fd, "BlockOwnedCells", &owned_counts);

    for (i, b) in m.boundaries.iter().enumerate() {
        let pg = root.create_group(&b.name).unwrap();
        text_attr(&pg, "Type", "PolyData");
        int_attr(&pg, "Index", i as i64 + 1);
        text_attr(&pg, "role", "wall");
        let mut local = std::collections::HashMap::new();
        let mut order = Vec::new();
        let (mut off, mut con, mut owner) = (vec![0i64], vec![], vec![]);
        for f in b.range() {
            for &p in m.face(f) {
                let k = *local.entry(p).or_insert_with(|| {
                    order.push(p);
                    order.len() as i64 - 1
                });
                con.push(k);
            }
            off.push(con.len() as i64);
            owner.push(file_id[m.owner[f] as usize]);
        }
        ds(&pg, "NumberOfPoints", &[order.len() as i64]);
        let p: Vec<[f64; 3]> = order.iter().map(|&p| m.points[p as usize]).collect();
        points(&pg, &p);
        let pd = pg.create_group("PointData").unwrap();
        if l.ids {
            let ids: Vec<i64> = order.iter().map(|&p| p as i64).collect();
            ds(&pd, "GlobalPointIds", &ids);
        }
        for (name, n, o, c) in [
            ("Vertices", 0, vec![0i64], vec![]),
            ("Lines", 0, vec![0], vec![]),
            ("Strips", 0, vec![0], vec![]),
            ("Polygons", b.len, off, con),
        ] {
            let sg = pg.create_group(name).unwrap();
            ds(&sg, "NumberOfCells", &[n as i64]);
            ds(&sg, "NumberOfConnectivityIds", &[c.len() as i64]);
            ds(&sg, "Offsets", &o);
            ds(&sg, "Connectivity", &c);
        }
        let cd = pg.create_group("CellData").unwrap();
        ds(&cd, "OwnerCell", &owner);
    }
}

/// A polyhedral mesh in two zones, and the mesh its `.vtu` reads as.
fn source(name: &str) -> (Mesh, Mesh) {
    let mut m = format::read(data("openfoam/dual_polyhedra"), None)
        .unwrap()
        .mesh;
    m.zones = vec!["fluid".into(), "rotor".into()];
    for c in 0..m.n_cells {
        m.cell_zone[c] = u32::from(c % 3 == 0);
    }
    let vtu = dir().join(format!("{name}.vtu"));
    write_vtu(&m, &vtu, true).unwrap();
    let r = format::read(&vtu, None).unwrap().mesh;
    (m, r)
}

/// A file order unlike the mesh's: odd cells backwards, then even ones.
fn scrambled(n: usize) -> Vec<usize> {
    (0..n)
        .rev()
        .filter(|c| c % 2 == 1)
        .chain((0..n).filter(|c| c % 2 == 0))
        .collect()
}

fn read(p: &Path) -> MeshFile {
    format::read(p, None).unwrap_or_else(|e| panic!("{e}"))
}

/// Invariant: a blocked file (polyhedra, ghosts, global ids, original cell
/// numbers) reads as the mesh of the `.vtu` of the same cells, bit for
/// bit, whatever the block size and the order of the cells in the file.
#[test]
fn blocked_file_reads_as_its_vtu() {
    let (m, vtu) = source("blocked_file_reads_as_its_vtu");
    assert!(m.n_cells > 20, "{} cells", m.n_cells);
    for (block, order) in [
        (m.n_cells, (0..m.n_cells).collect::<Vec<_>>()),
        (7, scrambled(m.n_cells)),
        (5, (0..m.n_cells).rev().collect()),
    ] {
        let p = dir().join(format!("blocked_{block}.vtkhdf"));
        write_blocked(
            &p,
            &m,
            &Layout {
                block,
                order,
                ids: true,
                original: true,
            },
        );
        let f = read(&p);
        assert_eq!(f.reader, "vtkhdf");
        let parts = m.n_cells.div_ceil(block);
        let s = if parts == 1 { "" } else { "s" };
        assert_eq!(
            f.format,
            format!("VTKHDF 2.4, {parts} partition{s}, blocked layout 1")
        );
        assert!(f.warnings.is_empty(), "{:?}", f.warnings);
        same(&vtu, &f.mesh, &format!("block {block}"));
    }
}

/// Invariant: without `GlobalPointIds` the partitions and the patches are
/// joined by their coordinates, and the mesh is the same.
#[test]
fn blocked_file_without_point_ids_joins_by_coordinates() {
    let (m, vtu) = source("blocked_file_without_point_ids_joins_by_coordinates");
    let p = dir().join("no_ids.vtkhdf");
    write_blocked(
        &p,
        &m,
        &Layout {
            block: 6,
            order: scrambled(m.n_cells),
            ids: false,
            original: true,
        },
    );
    let f = read(&p);
    assert_eq!(f.warnings.len(), 1, "{:?}", f.warnings);
    assert!(f.warnings[0].contains("joined by their coordinates"));
    same(&vtu, &f.mesh, "no ids");
}

/// Invariant: without `OriginalCellIds` the cells keep the file's order:
/// the same cells, zones and boundaries, numbered as in the file.
#[test]
fn blocked_file_without_original_ids_keeps_the_file_order() {
    let (m, vtu) = source("blocked_file_without_original_ids_keeps_the_file_order");
    let order = scrambled(m.n_cells);
    let p = dir().join("no_original.vtkhdf");
    write_blocked(
        &p,
        &m,
        &Layout {
            block: 9,
            order: order.clone(),
            ids: true,
            original: false,
        },
    );
    let f = read(&p).mesh;
    assert_eq!(f.n_cells, vtu.n_cells);
    assert_eq!(f.n_faces(), vtu.n_faces());
    let (ga, gb) = (Geometry::compute(&vtu), Geometry::compute(&f));
    for (k, &c) in order.iter().enumerate() {
        let (a, b) = (ga.cell_volume[c], gb.cell_volume[k]);
        assert!((a - b).abs() <= 1e-12 * a, "cell {c}: {a} {b}");
        assert_eq!(vtu.cell_zone[c], f.cell_zone[k]);
    }
    let count = |m: &Mesh| -> Vec<(String, usize)> {
        m.boundaries
            .iter()
            .map(|b| (b.name.clone(), b.len))
            .collect()
    };
    assert_eq!(count(&vtu), count(&f));
}

/// Invariant: FoehnCraft's own `fields.vtkhdf` (one partition, linear
/// cells and polyhedra, a real-valued zone array, boundary blocks without
/// point ids) reads back as its mesh: the same cells in the order of their
/// identities, the same volumes, zones and boundaries.
#[test]
fn own_output_reads_back_as_its_mesh() {
    use foehncraft_mesh::generate::{BoxSpec, Kind};
    let mut s = BoxSpec::new(Kind::Mixed, 3, [4, 3, 3]);
    s.perturb = 0.2;
    let mut m = s.generate().unwrap();
    m.zones = vec!["fluid".into(), "solid".into()];
    for c in 0..m.n_cells {
        m.cell_zone[c] = u32::from(c % 4 == 0);
    }
    let g = Geometry::compute(&m);
    let p = dir().join("own.vtkhdf");
    vtkhdf::write(
        &p,
        &m,
        &[CellField {
            name: "volume",
            values: &g.cell_volume,
        }],
    )
    .unwrap();
    let f = read(&p);
    assert_eq!(f.format, "VTKHDF 2.4, 1 partition");
    let r = f.mesh;
    r.validate().unwrap();
    assert_eq!(r.n_cells, m.n_cells);
    let gr = Geometry::compute(&r);
    let mut by_id: Vec<usize> = (0..m.n_cells).collect();
    by_id.sort_by_key(|&c| m.cell_id[c]);
    for (k, &c) in by_id.iter().enumerate() {
        let (a, b) = (g.cell_volume[c], gr.cell_volume[k]);
        assert!((a - b).abs() <= 1e-12 * a, "cell {c}: {a} {b}");
        // no zone names in the file: the zones are named by number
        assert_eq!(
            r.zones[r.cell_zone[k] as usize],
            format!("zone{}", m.cell_zone[c])
        );
    }
    let count = |m: &Mesh| -> Vec<(String, usize)> {
        m.boundaries
            .iter()
            .map(|b| (b.name.clone(), b.len))
            .collect()
    };
    assert_eq!(count(&m), count(&r));
}

/// Invariant: counts beyond the file and ids beyond the points are refused
/// with a message, not a panic or an allocation of the declared size.
#[test]
fn broken_counts_are_refused() {
    let (m, _) = source("broken_counts_are_refused");
    let p = dir().join("broken.vtkhdf");
    write_blocked(
        &p,
        &m,
        &Layout {
            block: 8,
            order: (0..m.n_cells).collect(),
            ids: true,
            original: true,
        },
    );
    let file = h5::File::open_rw(&p).unwrap();
    let g = file.group("VTKHDF/internal").unwrap();
    let d = g.dataset("NumberOfCells").unwrap();
    let mut v: Vec<i64> = d.read_raw().unwrap();
    let keep = v[0];
    v[0] = 1 << 40;
    d.write_raw(&v).unwrap();
    drop((d, g, file));
    let e = format::read(&p, None).unwrap_err();
    assert!(
        e.contains("NumberOfCells") && e.contains("more than the file"),
        "{e}"
    );

    let file = h5::File::open_rw(&p).unwrap();
    let g = file.group("VTKHDF/internal").unwrap();
    let d = g.dataset("NumberOfCells").unwrap();
    v[0] = keep;
    d.write_raw(&v).unwrap();
    let d = g.dataset("FaceConnectivity").unwrap();
    let mut c: Vec<i64> = d.read_raw().unwrap();
    c[3] = 1 << 20;
    d.write_raw(&c).unwrap();
    drop((d, g, file));
    let e = format::read(&p, None).unwrap_err();
    assert!(e.contains("point id") && e.contains("out of range"), "{e}");
}

/// Invariant: KoflCraft's blocked files read as the mesh of the `.vtu` of
/// the same export, bit for bit (the `koflcraft_*` root strings, the
/// patch roles and the block index are not part of the mesh).
#[test]
fn koflcraft_files_read_as_their_vtu() {
    for (h, v, parts) in [
        ("vtkhdf/box4.vtkhdf", "vtkhdf/box4.vtu", 1),
        ("vtkhdf/channel.vtkhdf", "vtkhdf/channel.vtu", 1),
    ] {
        let a = read(Path::new(&data(h)));
        let b = read(Path::new(&data(v)));
        let s = if parts == 1 { "" } else { "s" };
        assert_eq!(
            a.format,
            format!("VTKHDF 2.4, {parts} partition{s}, blocked layout 1")
        );
        assert!(a.warnings.is_empty(), "{h}: {:?}", a.warnings);
        same(&b.mesh, &a.mesh, h);
    }
}

/// Invariant: the channel file reproduces the summary of its
/// `mesh-evidence.json`: counts, volumes, patches with their faces and
/// areas, zones.
#[test]
fn koflcraft_channel_matches_its_evidence() {
    let text = std::fs::read_to_string(data("vtkhdf/channel-mesh-evidence.json")).unwrap();
    let ev: serde_json::Value = serde_json::from_str(&text).unwrap();
    let s = &ev["summary"];
    let m = read(Path::new(&data("vtkhdf/channel.vtkhdf"))).mesh;
    m.validate().unwrap();
    let n = |k: &str| s[k].as_u64().unwrap() as usize;
    assert_eq!(m.n_cells, n("cells"));
    assert_eq!(m.n_faces(), n("faces"));
    assert_eq!(m.n_interior_faces(), n("interior_faces"));
    assert_eq!(m.n_boundary_faces(), n("boundary_faces"));
    assert_eq!(m.points.len(), n("points"));
    let g = Geometry::compute(&m);
    let close = |a: f64, b: f64, what: &str| {
        assert!((a - b).abs() <= 1e-12 * b.abs(), "{what}: {a} against {b}");
    };
    let total: f64 = g.cell_volume[..m.n_cells].iter().sum();
    close(total, s["total_volume"].as_f64().unwrap(), "total volume");
    let (q, _) = QualityReport::compute(&m, &g);
    let patches = s["patches"].as_array().unwrap();
    assert_eq!(q.boundaries.len(), patches.len());
    for (b, p) in q.boundaries.iter().zip(patches) {
        assert_eq!(b.0, p["name"].as_str().unwrap());
        assert_eq!(b.1, p["faces"].as_u64().unwrap() as usize, "{}", b.0);
        close(b.2, p["area"].as_f64().unwrap(), &b.0);
    }
    let zones = s["zones"].as_array().unwrap();
    let names: Vec<&str> = zones.iter().map(|z| z["name"].as_str().unwrap()).collect();
    assert_eq!(m.zones, names);
}
