//! VTKHDF files for every generated mesh type; read back with VTK by
//! `tools/check/vtkhdf.py` (stage-1 acceptance: the files load in ParaView).
#![cfg(feature = "hdf5")]

use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::vtkhdf::{self, CellField};

#[test]
fn write_all_types() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("vtkhdf");
    std::fs::create_dir_all(&dir).unwrap();
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for dim in [2u8, 3] {
            let mut s = BoxSpec::new(kind, dim, [4, 3, 3]);
            s.perturb = if std::env::var("FOEHNCRAFT_NO_PERTURB").is_ok() {
                0.0
            } else {
                0.2
            };
            let m = s.generate().unwrap();
            let g = Geometry::compute(&m);
            let p = dir.join(format!("{}_{}d.vtkhdf", kind.name(), dim));
            vtkhdf::write(
                &p,
                &m,
                &[CellField {
                    name: "volume",
                    values: &g.cell_volume,
                }],
            )
            .unwrap();
            assert!(p.exists());
        }
    }
    // a renumbered mesh is written in input order
    let m = BoxSpec::new(Kind::Poly, 3, [4, 3, 3]).generate().unwrap();
    let g = Geometry::compute(&m);
    let mut r = m.clone();
    foehncraft_mesh::renumber::renumber(&mut r, foehncraft_mesh::renumber::Ordering::Shuffle, 5)
        .unwrap();
    let gr = Geometry::compute(&r);
    vtkhdf::write(
        dir.join("poly_input_order.vtkhdf"),
        &m,
        &[CellField {
            name: "volume",
            values: &g.cell_volume,
        }],
    )
    .unwrap();
    vtkhdf::write(
        dir.join("poly_shuffled.vtkhdf"),
        &r,
        &[CellField {
            name: "volume",
            values: &gr.cell_volume,
        }],
    )
    .unwrap();
    eprintln!("written to {}", dir.display());
}

/// Invariant: the size estimate of the run plan (`estimate_bytes`) is the
/// written file within 1 %, here for a mesh of two blocks with ghosts.
#[test]
fn size_estimate_matches_the_file() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("vtkhdf");
    std::fs::create_dir_all(&dir).unwrap();
    let m = BoxSpec::new(Kind::Hex, 3, [41, 41, 41]).generate().unwrap();
    assert!(m.n_cells > vtkhdf::BLOCK_CELLS);
    let a: Vec<f64> = (0..m.n_cells).map(|c| c as f64).collect();
    let p = dir.join("two_blocks.vtkhdf");
    let fields = [
        CellField {
            name: "a",
            values: &a,
        },
        CellField {
            name: "b",
            values: &a,
        },
    ];
    vtkhdf::write(&p, &m, &fields).unwrap();
    let size = std::fs::metadata(&p).unwrap().len() as f64;
    let est = vtkhdf::estimate_bytes(&m, fields.len());
    assert!(
        (est - size).abs() <= 0.01 * size,
        "estimate {est}, file {size}"
    );
}
