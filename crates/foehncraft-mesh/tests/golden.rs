//! Golden values for the test mesher and the partitioner (R4): a matrix of
//! mesh type × dimension × perturbation × refinement × 1, 2 and 4 ranks.
//! Counts must match exactly, real numbers to a relative 1e-9 (the libm
//! functions of the grading may differ in the last bit between platforms).
//! Every run also checks the invariants: closed cells, total volume, owned
//! volumes add up, every halo cell has exactly one owner.
//!
//! An unexplained difference is a bug report, not noise: regenerate the
//! file only for an intended change, with
//!
//! ```text
//! FOEHNCRAFT_BLESS=1 cargo test -p foehncraft-mesh --test mesh golden::
//! ```
//!
//! and say why in the commit message.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::decompose::decompose;
use foehncraft_mesh::generate::{BoxSpec, Grading, Kind};
use foehncraft_mesh::partition::{Method, PartitionStats, partition};
use foehncraft_mesh::quality::QualityReport;
use serde_json::{Map, Value, json};
use std::collections::BTreeMap;

const FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/mesher.json");

fn case_name(s: &BoxSpec) -> String {
    format!(
        "{}-{}d-p{}-r{}",
        s.kind.name(),
        s.dimension,
        s.perturb,
        s.refine
    )
}

fn matrix() -> Vec<BoxSpec> {
    let mut v = Vec::new();
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for dimension in [2u8, 3] {
            for perturb in [0.0, 0.2] {
                for refine in [0u32, 1] {
                    let mut s = BoxSpec::new(kind, dimension, [4, 3, 3]);
                    s.upper = [2.0, 1.0, 0.5];
                    s.grading.y = Grading {
                        ratio: 3.0,
                        symmetric: true,
                    };
                    s.perturb = perturb;
                    s.refine = refine;
                    v.push(s);
                }
            }
        }
    }
    v
}

/// Golden values of one mesh, after checking its invariants.
fn values(s: &BoxSpec) -> BTreeMap<String, Value> {
    let name = case_name(s);
    let m = s.generate().unwrap();
    m.validate().unwrap();
    let g = Geometry::compute(&m);
    let (q, _) = QualityReport::compute(&m, &g);
    assert!(q.ok(), "{name}");
    let vol = 2.0 * 1.0 * 0.5;
    assert!((q.total_volume - vol).abs() < 1e-12, "{name}");
    assert!(q.closure.max < 1e-13, "{name}");
    let mut v = BTreeMap::new();
    v.insert("cells".into(), json!(q.cells));
    v.insert("faces".into(), json!(q.faces));
    v.insert("interior_faces".into(), json!(q.interior_faces));
    v.insert("points".into(), json!(q.points));
    v.insert("total_volume".into(), json!(q.total_volume));
    v.insert(
        "non_orthogonality_max".into(),
        json!(q.non_orthogonality.max),
    );
    v.insert(
        "non_orthogonality_mean".into(),
        json!(q.non_orthogonality.mean),
    );
    v.insert("skewness_max".into(), json!(q.skewness.max));
    v.insert("skewness_mean".into(), json!(q.skewness.mean));
    v.insert("warp_max".into(), json!(q.warp.max));
    v.insert("aspect_ratio_max".into(), json!(q.aspect_ratio.max));
    v.insert("volume_ratio_max".into(), json!(q.volume_ratio.max));
    v.insert("faces_per_cell_mean".into(), json!(q.faces_per_cell.mean));
    for ranks in [1usize, 2, 4] {
        for method in [Method::Rcb, Method::Metis] {
            if method == Method::Metis && !cfg!(feature = "metis") {
                continue;
            }
            let part = partition(&m, &g, ranks, method).unwrap();
            let st = PartitionStats::compute(&m, &part, ranks);
            let locals = decompose(&m, &part, ranks, 1).unwrap();
            let mut owned_vol = 0.0;
            let mut owner = vec![usize::MAX; m.n_cells];
            for l in &locals {
                l.mesh.validate().unwrap();
                let lg = Geometry::compute(&l.mesh);
                owned_vol += lg.cell_volume[..l.mesh.n_cells].iter().sum::<f64>();
                for &c in &l.mesh.cell_id[..l.mesh.n_cells] {
                    assert_eq!(
                        owner[c as usize],
                        usize::MAX,
                        "{name}: cell {c} owned twice"
                    );
                    owner[c as usize] = l.rank;
                }
            }
            assert!((owned_vol - vol).abs() < 1e-12, "{name}: {ranks} ranks");
            let mut halo_total = 0;
            for l in &locals {
                for nb in &l.halo.neighbours {
                    for h in nb.recv.clone() {
                        // exactly one owner, and it is the rank the halo
                        // block names
                        assert_eq!(owner[l.mesh.cell_id[h] as usize], nb.rank, "{name}");
                    }
                }
                // R1: inner cells have no halo neighbour, outer cells do
                for c in 0..l.mesh.n_cells {
                    let touches = l
                        .mesh
                        .cell_face_iter(c)
                        .any(|(f, _)| l.mesh.across(f, c).is_some_and(|d| d >= l.mesh.n_cells));
                    assert_eq!(touches, c >= l.mesh.n_inner_cells, "{name}");
                }
                halo_total += l.mesh.n_halo;
            }
            if method == Method::Rcb {
                v.insert(format!("rcb{ranks}_cut_faces"), json!(st.cut_faces));
                v.insert(format!("rcb{ranks}_halo_cells"), json!(halo_total));
            }
        }
    }
    v
}

fn close(a: &Value, b: &Value) -> bool {
    match (a.as_u64(), b.as_u64()) {
        (Some(x), Some(y)) => x == y,
        _ => {
            let (x, y) = (a.as_f64().unwrap(), b.as_f64().unwrap());
            (x - y).abs() <= 1e-9 * x.abs().max(y.abs()) + 1e-14
        }
    }
}

fn command_output(cmd: &str, args: &[&str]) -> String {
    std::process::Command::new(cmd)
        .args(args)
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map_or("unknown".into(), |s| s.trim().to_string())
}

#[test]
fn golden_values_of_mesher_and_partitioner() {
    let mut cases = Map::new();
    for s in matrix() {
        let v: Map<String, Value> = values(&s).into_iter().collect();
        cases.insert(case_name(&s), Value::Object(v));
    }
    if std::env::var("FOEHNCRAFT_BLESS").is_ok() {
        let doc = json!({
            "metadata": {
                "what": "golden values of the stage-1 test mesher and partitioner (R4)",
                "commit": command_output("git", &["rev-parse", "--short", "HEAD"]),
                "os": format!("{} {}", std::env::consts::OS, std::env::consts::ARCH),
                "rust": command_output("rustc", &["--version"]),
            },
            "cases": cases,
        });
        std::fs::write(FILE, serde_json::to_string_pretty(&doc).unwrap() + "\n").unwrap();
        eprintln!("wrote {FILE}");
        return;
    }
    let text = std::fs::read_to_string(FILE).expect("golden file (bless it first)");
    let doc: Value = serde_json::from_str(&text).unwrap();
    let golden = doc["cases"].as_object().unwrap();
    let mut diffs = Vec::new();
    for (name, v) in &cases {
        let Some(g) = golden.get(name) else {
            diffs.push(format!("{name}: not in the golden file"));
            continue;
        };
        for (k, x) in v.as_object().unwrap() {
            match g.get(k) {
                Some(y) if close(x, y) => {}
                Some(y) => diffs.push(format!("{name}.{k}: {x} (golden {y})")),
                None => diffs.push(format!("{name}.{k}: not in the golden file")),
            }
        }
    }
    assert_eq!(golden.len(), cases.len(), "cases added or removed");
    assert!(
        diffs.is_empty(),
        "golden differences:\n{}",
        diffs.join("\n")
    );
}
