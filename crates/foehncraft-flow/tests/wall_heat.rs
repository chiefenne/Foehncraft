//! The wall heat values of the energy equation (stage 7, derived fields
//! `wall-heat-flux`, `heat-transfer-coefficient`, `nusselt-number`): exact
//! on the Couette flow heated by its friction, the laminar isothermal flat
//! plate against the similarity solution of Pohlhausen (1921), `foehn post
//! derive` against the file of the run, and the checks of the case keys.
#![cfg(feature = "hdf5")]

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::derived;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "foehncraft_wall_heat_{name}_{}",
        std::process::id()
    ))
}

fn load(file: &str, sets: &[&str]) -> FlowCase {
    let sets: Vec<String> = sets.iter().map(|s| s.to_string()).collect();
    FlowCase::from_file_with(&cases().join(file), &sets).unwrap()
}

const HEAT: &str =
    r#"output.derived=["wall-heat-flux","heat-transfer-coefficient","nusselt-number"]"#;

/// The values of `name` on the faces of boundary `b`, with the face
/// centroids.
fn on(
    mesh: &foehncraft_mesh::Mesh,
    geo: &foehncraft_mesh::Geometry,
    f: &derived::Fields,
    name: &str,
    b: &str,
) -> Vec<([f64; 3], f64)> {
    let v = f.get(name).unwrap();
    let ni = mesh.n_interior_faces();
    let p = mesh.boundaries.iter().find(|p| p.name == b).unwrap();
    p.range()
        .map(|i| (geo.face_centroid[i], v[i - ni]))
        .collect()
}

/// Couette flow heated by its friction between two walls at T₀ = 300 K:
/// the dissipation μU²/h² per volume leaves through both walls, so the heat
/// flux into the fluid is q_w = −μU²/(2h) = −0.05 W/m² on each wall (the
/// discrete heat balance closes, so this holds on any mesh, not only to
/// second order). The heat transfer coefficient of the wall cell
/// h = q_w/(T_w − T_P) is k/d_n = 0.8 W/(m² K) with d_n = 1/16 m, half the
/// cell height, and Nu = h L/k(T_w) = 16 with L = 1 m and k = 0.05 W/(m K).
/// `foehn post derive` writes the same file as the run.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn couette_wall_heat_values_are_exact() {
    let case = load(
        "couette_viscous_heating.toml",
        &[
            HEAT,
            "verification.levels=[0]",
            "output.reference.length=1.0",
            "run.restart_final=true",
            "output.vtkhdf=true",
        ],
    );
    let dir = out("couette");
    let rep = run(&case, "couette", &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let restart = dir.join("restart.h5");
    let (mesh, geo, f) = derived::evaluate(&case, &restart, &case.output.derived).unwrap();
    for b in ["ymin", "ymax"] {
        for (n, exact) in [
            ("wall_heat_flux", -0.05),
            ("wall_temperature", 300.0),
            ("heat_transfer_coefficient", 0.8),
            ("nusselt_number", 16.0),
        ] {
            for (x, v) in on(&mesh, &geo, &f, n, b) {
                assert!(
                    (v - exact).abs() < 1e-7 * exact.abs().max(1.0),
                    "{b} {n} at {x:?}: {v}, exact {exact}"
                );
            }
        }
    }
    // the transfer coefficient of the wall treatment (D124 W4): k/d_n of
    // the wall cells, d_n half the cell height 1/16
    let (_, _, w) = derived::wall_exchange_from_restart(&case, &restart).unwrap();
    let h = w.heat.expect("with the energy equation");
    assert_eq!(h.len(), w.faces.len());
    assert!(h.iter().all(|h| (h - 0.05 * 16.0).abs() < 1e-9), "{h:?}");
    // post derive: the file of the run
    let post = dir.join("derived.vtkhdf");
    derived::from_restart(&case, &restart, &case.output.derived, &post).unwrap();
    let a = std::fs::read(dir.join("fields.vtkhdf")).unwrap();
    let b = std::fs::read(&post).unwrap();
    assert!(a == b, "the files differ");
    let _ = std::fs::remove_dir_all(&dir);
}

/// The laminar isothermal flat plate at Re_L = 1e4, Pr 0.71 on the coarsest
/// mesh: Nu_x = q_w x/(k ΔT) against the free stream (ΔT = 1 K) within 1 %
/// of the similarity solution 0.294165 Re_x^(1/2) (tools/verification/
/// pohlhausen.py) from x = 0.05 to 0.9 m; the coefficient of the wall cell
/// h = q_w/(T_w − T_P) = k/d_n and Nu = h L/k = L/d_n with the normal
/// distance d_n of the wall cell (L = 1 m); and the wall heat flux summed
/// over the plate equal to the heat flow of the plate in the run summary.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_plate_meets_the_similarity_solution() {
    let case = load(
        "plate_heated_laminar.toml",
        &["verification.levels=[0]", "run.restart_final=true"],
    );
    let dir = out("plate");
    let rep = run(&case, "plate", &dir, |_| {}).unwrap();
    assert!(rep.failures.is_empty(), "{:?}", rep.failures);
    let (mesh, geo, f) =
        derived::evaluate(&case, &dir.join("restart.h5"), &case.output.derived).unwrap();
    let k = 0.14084507042253522;
    let q = on(&mesh, &geo, &f, "wall_heat_flux", "plate");
    let mut checked = 0;
    for (x, q) in &q {
        if !(0.05..=0.9).contains(&x[0]) {
            continue;
        }
        let exact = 0.294165 * (1e4 * x[0]).sqrt();
        let nu = q * x[0] / k;
        assert!(
            (nu / exact - 1.0).abs() < 1e-2,
            "Nu_x at x = {}: {nu} against {exact}",
            x[0]
        );
        checked += 1;
    }
    assert!(checked > 40, "{checked} faces");
    // h = k/d_n and Nu = L/d_n of the wall cells
    let ni = mesh.n_interior_faces();
    let p = mesh.boundaries.iter().find(|p| p.name == "plate").unwrap();
    let h = f.get("heat_transfer_coefficient").unwrap();
    let nu = f.get("nusselt_number").unwrap();
    let mut total = 0.0;
    let qa = f.get("wall_heat_flux").unwrap();
    for i in p.range() {
        let s = geo.face_area[i];
        let a = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        let c = mesh.owner[i] as usize;
        let d: Vec<f64> = (0..3)
            .map(|j| geo.face_centroid[i][j] - geo.cell_centroid[c][j])
            .collect();
        let dn = (d[0] * s[0] + d[1] * s[1] + d[2] * s[2]).abs() / a;
        let b = i - ni;
        assert!((h[b] * dn / k - 1.0).abs() < 1e-9, "h {} d_n {dn}", h[b]);
        assert!((nu[b] * dn - 1.0).abs() < 1e-9, "Nu {} d_n {dn}", nu[b]);
        total += qa[b] * a;
    }
    let s = &rep.solves[0].scalars[0];
    let flow = s
        .boundaries
        .iter()
        .find(|b| b.boundary == "plate")
        .and_then(|b| b.heat_flow)
        .unwrap();
    assert!(
        (total / flow - 1.0).abs() < 1e-9,
        "Σ q_w A {total}, heat flow {flow}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn checks_of_the_wall_heat_keys() {
    let code = |e: anyhow::Error| foehncraft_flow::diag::of_error(&e).code;
    // without the energy equation
    let base = std::fs::read_to_string(cases().join("couette_hex.toml")).unwrap();
    let e = FlowCase::from_toml(
        &format!("{base}\n[output]\nderived = [\"wall-heat-flux\"]"),
        &cases(),
    )
    .unwrap_err();
    assert_eq!(code(e), "CASE005");
    // the Couette flow heated by its friction
    // line ends as on Linux (a Windows checkout has CRLF)
    let heated = std::fs::read_to_string(cases().join("couette_viscous_heating.toml"))
        .unwrap()
        .replace("\r\n", "\n");
    let with = |reference: &str, derived: &str| {
        let t = heated.replace(
            "[output]\n",
            &format!("{reference}\n\n[output]\nderived = {derived}\n"),
        );
        FlowCase::from_toml(&t, &cases())
    };
    with("", r#"["wall-heat-flux", "heat-transfer-coefficient"]"#).unwrap();
    // the Nusselt number needs a reference length
    let e = with("", r#"["nusselt-number"]"#).unwrap_err();
    let d = foehncraft_flow::diag::of_error(&e);
    assert_eq!(d.code, "CASE005");
    assert!(d.message.contains("reference length"), "{}", d.message);
    with("[output.reference]\nlength = 1.0", r#"["nusselt-number"]"#).unwrap();
    let e = with("[output.reference]\nlength = 0.0", "[]").unwrap_err();
    assert_eq!(code(e), "CASE004");
}
