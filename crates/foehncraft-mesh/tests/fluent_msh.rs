//! Fluent `.msh` reader (D191): hexahedral and polyhedral meshes, zones,
//! both face orientations, `.msh` shared with Gmsh, and clear errors.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::case::{MeshSection, build_mesh};
use foehncraft_mesh::format::{self, MeshFile};
use foehncraft_mesh::quality::QualityReport;

fn path(name: &str) -> String {
    format!(
        "{}/tests/data/fluent_msh/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn check(f: &MeshFile, volume: f64) -> QualityReport {
    f.mesh.validate().unwrap_or_else(|e| panic!("{e}"));
    let geo = Geometry::compute(&f.mesh);
    let (q, _) = QualityReport::compute(&f.mesh, &geo);
    assert!(q.ok(), "{}", q.to_markdown("fluent_msh"));
    assert!(
        ((q.total_volume - volume) / volume).abs() < 1e-12,
        "{}",
        q.total_volume
    );
    assert!(((q.boundary_volume - volume) / volume).abs() < 1e-12);
    q
}

fn area(q: &QualityReport, name: &str) -> f64 {
    q.boundaries
        .iter()
        .find(|b| b.0 == name)
        .unwrap_or_else(|| panic!("no boundary {name}"))
        .2
}

fn names(f: &MeshFile) -> Vec<&str> {
    f.mesh.boundaries.iter().map(|b| b.name.as_str()).collect()
}

#[test]
fn hexahedral_channel() {
    let f = format::read(path("channel_hex.msh"), None).unwrap();
    assert_eq!(f.reader, "fluent_msh");
    let q = check(&f, 2.0);
    assert_eq!(f.mesh.n_cells, 16);
    assert_eq!(names(&f), ["inlet", "outlet", "walls", "sides"]);
    assert!((area(&q, "inlet") - 1.0).abs() < 1e-12);
    assert!((area(&q, "walls") - 4.0).abs() < 1e-12);
    assert_eq!(f.mesh.zones, ["right", "left"]);
    assert_eq!(f.mesh.cell_zone.iter().filter(|&&z| z == 0).count(), 8);
    assert!(f.warnings.is_empty(), "{:?}", f.warnings);
}

#[test]
fn face_orientation_does_not_matter() {
    let a = format::read(path("channel_hex.msh"), None).unwrap();
    let b = format::read(path("channel_hex_flipped.msh"), None).unwrap();
    check(&b, 2.0);
    assert_eq!(names(&a), names(&b));
    assert_eq!(a.mesh.n_cells, b.mesh.n_cells);
    assert_eq!(a.mesh.owner.len(), b.mesh.owner.len());
}

#[test]
fn polyhedral_cells_and_mixed_faces() {
    let f = format::read(path("dual_polyhedra.msh"), None).unwrap();
    let q = check(&f, 1.0);
    assert_eq!(f.mesh.n_cells, 81);
    assert_eq!(names(&f), ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]);
    for n in names(&f) {
        assert!((area(&q, n) - 1.0).abs() < 1e-12, "{n}");
    }
    // the same mesh as the OpenFOAM file it was made from
    let o = format::read(
        format!(
            "{}/tests/data/openfoam/dual_polyhedra",
            env!("CARGO_MANIFEST_DIR")
        ),
        None,
    )
    .unwrap();
    assert_eq!(o.mesh.n_cells, f.mesh.n_cells);
    assert_eq!(o.mesh.owner.len(), f.mesh.owner.len());
}

#[test]
fn shares_the_extension_with_gmsh() {
    // Gmsh files with the same extension still go to the Gmsh reader
    for n in [
        "box_tet_ascii.msh",
        "box_tet_binary.msh",
        "channel_mixed_ascii.msh",
    ] {
        let p = format!("{}/tests/data/{n}", env!("CARGO_MANIFEST_DIR"));
        assert_eq!(format::read(&p, None).unwrap().reader, "gmsh", "{n}");
    }
    assert_eq!(
        format::read(path("channel_hex.msh"), None).unwrap().reader,
        "fluent_msh"
    );
    assert_eq!(
        format::read(path("channel_hex.msh"), Some("fluent_msh"))
            .unwrap()
            .reader,
        "fluent_msh"
    );
    assert!(format::read(path("channel_hex.msh"), Some("gmsh")).is_err());
}

#[test]
fn case_path_with_renamed_zones_and_scale() {
    let sec = MeshSection {
        file: Some(path("channel_hex.msh").into()),
        scale: 0.1,
        boundaries: [("walls", "wall")]
            .into_iter()
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect(),
        ..Default::default()
    };
    let (mesh, source, _) = build_mesh(&sec).unwrap();
    assert!(source.contains("Fluent"), "{source}");
    let n: Vec<&str> = mesh.boundaries.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(n, ["inlet", "outlet", "wall", "sides"]);
    let v: f64 = Geometry::compute(&mesh).cell_volume.iter().sum();
    assert!((v - 2.0e-3).abs() < 1e-15, "{v}");
}

#[test]
fn bad_files_give_clear_errors() {
    let e = |s: &str| {
        format::fluent_msh::IndustryMsh
            .parse(s.as_bytes(), &Default::default())
            .unwrap_err()
    };
    use format::MeshFormat;
    assert!(e("(2 2)\n").contains("3D"));
    assert!(e("(2 3)\n(2010 (1 1 4 1 3)(\n").contains("not closed"));
    assert!(e("(2 3)\n(2010 (1 1 4 1 3)())\n").contains("binary"));
    assert!(e("(2 3)\nhello\n").contains("line 2"));
    assert!(e("(2 3)\n(0 \"c\")\n").contains("no faces"));
    // a face of an unknown cell
    let s = "(2 3)\n(10 (0 1 3 0 3))\n(10 (1 1 3 1 3)(\n0 0 0\n1 0 0\n0 1 0\n))\n\
             (12 (0 1 1 0 0))\n(13 (3 1 1 3 3)(\n1 2 3 5 0\n))\n";
    assert!(e(s).contains("beyond"), "{}", e(s));
}
