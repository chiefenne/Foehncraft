//! Steam furnace stage SF1: species transport of an ideal-gas mixture.
//! Exact checks: the decay of a cosine mode of binary diffusion, Σ Y = 1
//! in every cell, the closure of the species and element accounts of the
//! ledger, and the refusals of the case reader.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_species_{name}_{}", std::process::id()))
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

/// The code and message a case text is refused with.
fn refused(text: &str) -> (String, String) {
    let e = match FlowCase::from_toml(text, &cases()) {
        Err(e) => e,
        Ok(c) => match run(&c, "refused", &out("refused"), |_| {}) {
            Err(e) => e,
            Ok(_) => panic!("the case was accepted"),
        },
    };
    (diag::of_error(&e).code.to_string(), format!("{e:#}"))
}

/// Binary diffusion of CH₄ in N₂ at rest: Y = 0.5 + 0.4 cos(πx)
/// exp(−Dπ²t) to t = 50 s with BDF2 on 100 cells, L2 error below 1e-4 of
/// the amplitude; the closure fraction is exactly 1 − Y_CH4 and the
/// accounts of the closed slab close (no flow through any boundary).
#[test]
fn binary_diffusion_decays_as_the_exact_mode() {
    let rep = run_case(&parse(&text("species_binary_diffusion.toml")), "binary");
    let s = &rep.solves[0];
    let y = s.scalars.iter().find(|x| x.name == "Y_CH4").unwrap();
    let l2 = y.l2.unwrap();
    assert!(l2 < 1e-4, "L2 error {l2}");
    assert!(y.min >= 0.0 && y.max <= 1.0);
    let b = s.balance.as_ref().unwrap();
    for q in [
        "species CH4",
        "species N2",
        "element C",
        "element H",
        "element N",
    ] {
        let a = b
            .last
            .ledger
            .iter()
            .find(|a| a.quantity == q)
            .unwrap_or_else(|| panic!("no account {q}"));
        assert!(a.imbalance.abs() < 1e-12, "{q}: imbalance {}", a.imbalance);
    }
}

/// The coaxial jet of the gallery: the mixture closes its element
/// accounts to the solver tolerance, every mass fraction stays in [0, 1]
/// within the smooth-extremum allowance of the MLP limiter (ε = 0.01 of the
/// range, here a few 1e-6), and the jet density is that of methane,
/// p₀ M_CH4/(R T) within the same allowance.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn coaxial_jet_closes_the_element_accounts() {
    let rep = run_case(&parse(&text("species_coaxial_jet.toml")), "jet");
    let s = &rep.solves[0];
    assert!(s.converged, "not converged");
    for n in ["Y_CH4", "Y_O2"] {
        let y = s.scalars.iter().find(|x| x.name == n).unwrap();
        assert!(
            y.min >= -1e-5 && y.max <= 1.0 + 1e-5,
            "{n}: [{}, {}]",
            y.min,
            y.max
        );
    }
    let b = s.balance.as_ref().unwrap();
    for e in ["C", "H", "O", "N"] {
        let a = b
            .last
            .ledger
            .iter()
            .find(|a| a.quantity == format!("element {e}"))
            .unwrap();
        assert!(a.relative < 1e-6, "element {e}: {}", a.relative);
    }
    let d = s.density.as_ref().unwrap();
    let rho_ch4 = 101_325.0 * 0.016_042_46 / (8.314_462_618 * 300.0);
    assert!(
        // within the MLP allowance: Y_CH4 = 1 + 3e-6 at the nozzle lowers M by as much
        (d.rho_min / rho_ch4 - 1.0).abs() < 1e-5,
        "{} {rho_ch4}",
        d.rho_min
    );
}

const BASE: &str = r#"
case_format = "0.7"
[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 0.1, 0.01]
cells = [10, 2, 1]
[fluid]
density = 1.0
viscosity = 1.0e-3
pressure_values = "gauge"
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "slip-wall"
[boundary.ymax]
type = "slip-wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
"#;

#[test]
fn the_reader_refuses_broken_species_tables() {
    let sp = "[species]\nnames = [\"CH4\", \"N2\"]\nproperties = \"fluid\"\n";
    let ini = "initial = { mole_fractions = { N2 = 1.0 } }\n";
    let inlet = "[boundary.xmin.species]\nmole_fractions = { CH4 = 1.0 }\n";
    // an unknown species
    let (c, m) = refused(&format!(
        "{BASE}[species]\nnames = [\"CH4\", \"XY\"]\nproperties = \"fluid\"\n{ini}{inlet}"
    ));
    assert_eq!(c, "SPE001", "{m}");
    assert!(m.contains("XY"), "{m}");
    // an inlet without its composition
    let (c, m) = refused(&format!("{BASE}{sp}{ini}"));
    assert_eq!(c, "SPE001", "{m}");
    assert!(m.contains("xmin"), "{m}");
    // fractions above one
    let (c, m) = refused(&format!(
        "{BASE}{sp}{ini}[boundary.xmin.species]\nmole_fractions = {{ CH4 = 0.7, N2 = 0.7 }}\n"
    ));
    assert_eq!(c, "SPE001", "{m}");
    // the mixture properties need the energy equation
    let (c, m) = refused(&format!(
        "{BASE}[species]\nnames = [\"CH4\", \"N2\"]\n{ini}{inlet}"
    ));
    assert_eq!(c, "SPE001", "{m}");
    assert!(m.contains("energy"), "{m}");
    // a boundary table without [species]
    let (c, m) = refused(&format!("{BASE}{inlet}"));
    assert_eq!(c, "SPE001", "{m}");
    // a valid table is accepted
    assert!(FlowCase::from_toml(&format!("{BASE}{sp}{ini}{inlet}"), &cases()).is_ok());
}
