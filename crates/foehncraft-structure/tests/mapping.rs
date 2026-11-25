//! Verification of the temperature mapping from a flow result onto the
//! structure (MS6): `[[temperature_fields]]` reads the cell temperatures
//! of the solid zones of a VTKHDF result and fits a linear field to the
//! nearest cell centres at each structural node.
//!
//! - A linear temperature field is mapped exactly on hexahedral,
//!   tetrahedral and polyhedral cells, and a free body under it expands
//!   without stress: ε_th = α (T − T₀) I is compatible when T is linear
//!   (Timoshenko & Goodier, *Theory of Elasticity*, 1970, §154), and the
//!   quadratic elements hold the quadratic displacement exactly.
//! - Plates and beams take the value at the mid-surface or axis and the
//!   slope across the section: a field T₀ + ΔT + G z gives what
//!   `[[temperatures]]` with `change = ΔT` and `difference = G h` gives.
//! - A smooth field converges at second order with the cell size.
//! - A 2D flow case (one layer of cells) maps a field linear in x and y
//!   exactly; the field is constant across the layer.

#![allow(clippy::needless_range_loop)]

use foehncraft_mesh::generate::{BoxSpec, Kind, ZoneBox};
use foehncraft_mesh::geometry::Geometry;
use foehncraft_mesh::vtkhdf::{self, CellField};
use foehncraft_structure::case::{self, Case};
use foehncraft_structure::codes;
use std::path::PathBuf;

const E: f64 = 200e9;
const ALPHA: f64 = 1.2e-5;
const T0: f64 = 293.15;

/// A flow mesh of the unit box with the zone `metal` over `lo..hi`,
/// written with `field(centroid)` as `temperature`.
fn flow_file(
    name: &str,
    kind: Kind,
    dim: u8,
    cells: [usize; 3],
    upper: [f64; 3],
    zone: ([f64; 3], [f64; 3]),
    field: impl Fn([f64; 3]) -> f64,
) -> PathBuf {
    let mut spec = BoxSpec::new(kind, dim, cells);
    spec.upper = upper;
    spec.zones = vec![ZoneBox {
        lower: zone.0,
        upper: zone.1,
        name: "metal".into(),
    }];
    let mesh = spec.generate().unwrap();
    let geo = Geometry::compute(&mesh);
    // the fluid cells get a value far off: they must not be used
    let t: Vec<f64> = (0..mesh.n_cells)
        .map(|c| {
            if mesh.zones[mesh.cell_zone[c] as usize] == "metal" {
                field(geo.cell_centroid[c])
            } else {
                1e6
            }
        })
        .collect();
    let dir = std::env::temp_dir().join(format!("foehncraft-mapping-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("{name}.vtkhdf"));
    vtkhdf::write(
        &path,
        &mesh,
        &[CellField {
            name: "temperature",
            values: &t,
        }],
    )
    .unwrap();
    path
}

const MATERIAL: &str = r#"
[analysis]
static = true
[materials.steel]
young = 200e9
poisson = 0.3
density = 7800.0
expansion = 1.2e-5
"#;

fn linear(x: [f64; 3]) -> f64 {
    350.0 + 60.0 * x[0] - 25.0 * x[1] + 40.0 * x[2]
}

/// The largest |(r_p − r_q)·d|/|d|² over pairs of nodes, with r = u − the
/// exact thermal displacement: zero when r is a rigid motion.
fn not_rigid(nodes: &[[f64; 3]], u: &[f64], exact: impl Fn([f64; 3]) -> [f64; 3]) -> f64 {
    let r: Vec<[f64; 3]> = nodes
        .iter()
        .enumerate()
        .map(|(n, &x)| {
            let e = exact(x);
            [0, 1, 2].map(|d| u[6 * n + d] - e[d])
        })
        .collect();
    let mut worst = 0.0f64;
    let step = (nodes.len() / 97).max(1);
    for p in (0..nodes.len()).step_by(step) {
        for q in (p + 1..nodes.len()).step_by(step) {
            let d = [0, 1, 2].map(|i| nodes[p][i] - nodes[q][i]);
            let dd: f64 = d.iter().map(|v| v * v).sum();
            let s: f64 = (0..3).map(|i| (r[p][i] - r[q][i]) * d[i]).sum();
            worst = worst.max((s / dd).abs());
        }
    }
    worst
}

#[test]
fn linear_field_is_mapped_exactly_and_a_free_solid_expands_without_stress() {
    let zone = ([0.2, 0.2, 0.2], [0.8, 0.6, 0.7]);
    for (kind, tag) in [(Kind::Hex, "hex"), (Kind::Tet, "tet"), (Kind::Poly, "poly")] {
        let mut spec_cells = [10, 10, 10];
        if kind == Kind::Hex {
            spec_cells = [12, 11, 13];
        }
        let file = flow_file(
            &format!("linear-{tag}"),
            kind,
            3,
            spec_cells,
            [1.0, 1.0, 1.0],
            zone,
            linear,
        );
        for element in ["hex20", "tet10"] {
            let text = format!(
                r#"{MATERIAL}
[[solids]]
name = "b"
material = "steel"
block = {{ origin = [0.2, 0.2, 0.2], size = [0.6, 0.4, 0.5], divisions = [3, 2, 2], element = "{element}" }}
[[supports]]
at = [0.2, 0.2, 0.2]
fix = ["ux", "uy", "uz"]
[[supports]]
at = [0.8, 0.2, 0.2]
fix = ["uy", "uz"]
[[supports]]
at = [0.2, 0.6, 0.2]
fix = ["uz"]
[[temperature_fields]]
parts = ["b"]
file = '{}'   # a literal string: Windows paths have backslashes
zones = ["metal"]
reference = {T0}
"#,
                file.display()
            );
            let out = case::run(&Case::from_toml(&text).unwrap()).unwrap();
            let m = &out.built.model;
            for (n, &x) in m.nodes.iter().enumerate() {
                let t = out.built.node_temperature[n];
                assert!(
                    (t - linear(x)).abs() < 1e-9 * linear(x),
                    "{tag}/{element}: node {x:?}: {t} against {}",
                    linear(x)
                );
            }
            let u = &out.statics.as_ref().unwrap().u;
            // u = α [(T − T₀) x − ½ g |x|²] plus a rigid motion
            let g = [60.0, -25.0, 40.0];
            let exact = |x: [f64; 3]| {
                let x2: f64 = x.iter().map(|v| v * v).sum();
                [0, 1, 2].map(|d| ALPHA * ((linear(x) - T0) * x[d] - 0.5 * g[d] * x2))
            };
            let scale = ALPHA * 100.0;
            let err = not_rigid(&m.nodes, u, exact);
            assert!(err < 1e-9 * scale, "{tag}/{element}: {err:e}");
            let vm = out.solid_von_mises(u);
            let big = vm.iter().cloned().fold(0.0, f64::max);
            assert!(big < 1e-8 * E * scale, "{tag}/{element}: von Mises {big:e}");
            let rep = &out.to_json("t")["temperature_fields"][0];
            assert_eq!(rep["nodes"].as_u64().unwrap() as usize, m.nodes.len());
            assert!(rep["source_range"][1].as_f64().unwrap() < 1e3, "{rep}");
        }
    }
}

#[test]
fn plates_and_beams_take_the_value_and_the_slope_across_the_section() {
    // T = T₀ + 40 + 3000 (z − z_c): 40 K above the stress-free temperature
    // at the mid-plane z_c, 6 K hotter on top than below (h = 2 mm)
    let zc = 0.5;
    let field = move |x: [f64; 3]| T0 + 40.0 + 3000.0 * (x[2] - zc);
    let file = flow_file(
        "slope",
        Kind::Hex,
        3,
        [10, 10, 8],
        [0.04, 0.02, 1.0],
        ([0.0, 0.0, 0.4], [0.04, 0.02, 0.6]),
        field,
    );
    let h = 2e-3;
    for element in ["mitc4", "mitc3+"] {
        let plate = |temps: &str| {
            format!(
                r#"{MATERIAL}
[[plates]]
name = "p"
origin = [0.005, 0.004, {zc}]
edge_a = [0.03, 0.0, 0.0]
edge_b = [0.0, 0.012, 0.0]
divisions = [6, 3]
thickness = {h}
material = "steel"
element = "{element}"
[[supports]]
at = [0.005, 0.004, {zc}]
fix = ["all"]
[[supports]]
at = [0.035, 0.004, {zc}]
fix = ["uy"]
{temps}
"#
            )
        };
        let mapped = case::run(
            &Case::from_toml(&plate(&format!(
                "[[temperature_fields]]\nparts = [\"p\"]\nfile = '{}'\nzones = [\"metal\"]\nreference = {T0}\n",
                file.display()
            )))
            .unwrap(),
        )
        .unwrap();
        let given = case::run(
            &Case::from_toml(&plate(&format!(
                "[[temperatures]]\npart = \"p\"\nchange = 40.0\ndifference = {}\n",
                3000.0 * h
            )))
            .unwrap(),
        )
        .unwrap();
        let a = &mapped.statics.as_ref().unwrap().u;
        let b = &given.statics.as_ref().unwrap().u;
        let big = b.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(big > 1e-6);
        for i in 0..a.len() {
            assert!((a[i] - b[i]).abs() < 1e-9 * big, "{element}: {i}");
        }
    }
    let beam = |temps: &str| {
        format!(
            r#"{MATERIAL}
[sections.s]
shape = "rectangle"
width = 4e-3
thickness = {h}
[[beams]]
name = "r"
from = [0.005, 0.01, {zc}]
to = [0.035, 0.01, {zc}]
elements = 8
section = "s"
material = "steel"
[[supports]]
at = [0.005, 0.01, {zc}]
fix = ["all"]
{temps}
"#
        )
    };
    let mapped = case::run(
        &Case::from_toml(&beam(&format!(
            "[[temperature_fields]]\nparts = [\"r\"]\nfile = '{}'\nzones = [\"metal\"]\nreference = {T0}\n",
            file.display()
        )))
        .unwrap(),
    )
    .unwrap();
    let given = case::run(
        &Case::from_toml(&beam(&format!(
            "[[temperatures]]\npart = \"r\"\nchange = 40.0\ndifference = {}\n",
            3000.0 * h
        )))
        .unwrap(),
    )
    .unwrap();
    let a = &mapped.statics.as_ref().unwrap().u;
    let b = &given.statics.as_ref().unwrap().u;
    let big = b.iter().map(|v| v.abs()).fold(0.0, f64::max);
    for i in 0..a.len() {
        assert!((a[i] - b[i]).abs() < 1e-9 * big, "beam: {i}");
    }
}

#[test]
fn smooth_field_converges_at_second_order() {
    let pi = std::f64::consts::PI;
    let smooth =
        move |x: [f64; 3]| 400.0 + 50.0 * (pi * x[0]).sin() * (pi * x[1]).cos() * x[2].exp();
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut spec_cells = [n, n, n];
        spec_cells[2] = n / 2;
        let file = flow_file(
            &format!("smooth-{n}"),
            Kind::Hex,
            3,
            spec_cells,
            [1.0, 1.0, 0.5],
            ([0.0, 0.0, 0.0], [1.0, 1.0, 0.5]),
            smooth,
        );
        let text = format!(
            r#"{MATERIAL}
[[solids]]
name = "b"
material = "steel"
block = {{ origin = [0.1, 0.1, 0.05], size = [0.8, 0.8, 0.4], divisions = [4, 4, 2], element = "hex20" }}
[[supports]]
part = "b"
surface = "zmin"
fix = ["ux", "uy", "uz"]
[[temperature_fields]]
parts = ["b"]
file = '{}'   # a literal string: Windows paths have backslashes
reference = {T0}
"#,
            file.display()
        );
        let out = case::run(&Case::from_toml(&text).unwrap()).unwrap();
        let m = &out.built.model;
        let err = m
            .nodes
            .iter()
            .enumerate()
            .map(|(k, &x)| (out.built.node_temperature[k] - smooth(x)).abs())
            .fold(0.0, f64::max);
        errors.push(err);
    }
    for w in errors.windows(2) {
        let order = (w[0] / w[1]).log2();
        assert!(order > 1.8, "{errors:?}");
    }
    // 32 cells: within 0.4 K (0.5 % of the 78 K swing of the field); the
    // largest errors sit on the surface, where the fit extrapolates
    assert!(errors[2] < 0.4, "{errors:?}");
}

#[test]
fn a_two_dimensional_case_maps_a_field_linear_in_x_and_y() {
    // one layer of cells, 0.01 thick: the centres span x and y only
    let f2 = |x: [f64; 3]| 320.0 + 80.0 * x[0] + 30.0 * x[1];
    for kind in [Kind::Hex, Kind::Tet] {
        let file = flow_file(
            "flat",
            kind,
            2,
            [30, 12, 1],
            [0.3, 0.06, 0.01],
            ([0.05, 0.0, 0.0], [0.25, 0.02, 0.01]),
            f2,
        );
        let text = format!(
            r#"{MATERIAL}
[[solids]]
name = "b"
material = "steel"
block = {{ origin = [0.05, 0.0, 0.0], size = [0.2, 0.02, 0.01], divisions = [10, 2, 1], element = "hex20" }}
[[supports]]
part = "b"
surface = "ymin"
fix = ["ux", "uy", "uz"]
[[temperature_fields]]
parts = ["b"]
file = '{}'   # a literal string: Windows paths have backslashes
zones = ["metal"]
reference = {T0}
"#,
            file.display()
        );
        let out = case::run(&Case::from_toml(&text).unwrap()).unwrap();
        for (k, &x) in out.built.model.nodes.iter().enumerate() {
            let t = out.built.node_temperature[k];
            assert!((t - f2(x)).abs() < 1e-9 * f2(x), "{kind:?} {x:?}: {t}");
        }
    }
}

#[test]
fn temperature_field_errors() {
    let file = flow_file(
        "errors",
        Kind::Hex,
        3,
        [6, 6, 6],
        [1.0, 1.0, 1.0],
        ([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]),
        linear,
    );
    let case_with = |block: &str, keys: &str, expansion: bool| {
        let mat = if expansion {
            MATERIAL.to_string()
        } else {
            MATERIAL.replace("expansion = 1.2e-5\n", "")
        };
        format!(
            r#"{mat}
[[solids]]
name = "b"
material = "steel"
block = {{ origin = {block}, size = [0.3, 0.3, 0.3], divisions = [1, 1, 1], element = "hex20" }}
[[supports]]
part = "b"
surface = "zmin"
fix = ["ux", "uy", "uz"]
[[temperature_fields]]
file = '{}'   # a literal string: Windows paths have backslashes
reference = {T0}
{keys}
"#,
            file.display()
        )
    };
    let inside = "[0.1, 0.1, 0.1]";
    let cases = [
        (
            case_with(inside, "parts = [\"b\"]\nzones = [\"steel\"]", true),
            "no zone `steel`",
        ),
        (
            case_with(inside, "parts = [\"b\"]\nfield = \"pressure\"", true),
            "no cell field `pressure`",
        ),
        (
            case_with(inside, "parts = [\"c\"]", true),
            "no beam line, plate or solid `c`",
        ),
        (case_with(inside, "parts = []", true), "names no part"),
        (
            case_with(inside, "parts = [\"b\"]", false),
            "needs `expansion`",
        ),
        (
            case_with(
                "[3.0, 0.1, 0.1]",
                "parts = [\"b\"]\nzones = [\"metal\"]",
                true,
            ),
            "does not lie in the solid region",
        ),
        (
            case_with(inside, "parts = [\"b\"]\nneighbours = 2", true),
            "`neighbours` must be 4 to 1000",
        ),
        (
            case_with(inside, "parts = [\"b\"]\nmax_distance = -1.0", true),
            "`max_distance` must be positive",
        ),
    ];
    for (text, want) in cases {
        let e = case::run(&Case::from_toml(&text).unwrap()).unwrap_err();
        assert_eq!(e.code, codes::CASE);
        assert!(e.message.contains(want), "{want}: {}", e.message);
    }
    // a file that is not there
    let text = case_with(inside, "parts = [\"b\"]", true)
        .replace(&file.display().to_string(), "/nonexistent/flow.vtkhdf");
    let e = case::run(&Case::from_toml(&text).unwrap()).unwrap_err();
    assert_eq!(e.code, codes::CASE);
    // the offset moves the structure into the zone
    let text = case_with(
        "[3.1, 0.1, 0.1]",
        "parts = [\"b\"]\nzones = [\"metal\"]\noffset = [-3.0, 0.0, 0.0]",
        true,
    );
    let out = case::run(&Case::from_toml(&text).unwrap()).unwrap();
    let x = out.built.model.nodes[0];
    let t = out.built.node_temperature[0];
    let want = linear([x[0] - 3.0, x[1], x[2]]);
    assert!((t - want).abs() < 1e-9 * want, "{t} {want}");
}
