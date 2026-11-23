//! Start values per region and per fluid part (D238), docs/requirements.md
//! D238 and the codes REG007–REG010: region tables, `"inlet"` temperatures,
//! every fluid part at its own boundary pressure, and the report.

use foehncraft_flow::case::{FlowCase, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "foehncraft_start_values_{name}_{}",
        std::process::id()
    ))
}

fn two_fluids() -> String {
    std::fs::read_to_string(cases().join("cht_two_fluids_absolute.toml"))
        .unwrap()
        .replace("\r\n", "\n")
        .replace("max_iterations = 3000", "max_iterations = 5")
}

/// V7.45: the water and the oil start at their inlet temperatures, the
/// steel plate at its own; with gauge values and two outlets at 0 and
/// 5·10⁴ Pa each part starts at its own outlet pressure. The log (REG010)
/// and `flow.json` (`start_values`) say so.
#[test]
fn every_region_and_part_starts_on_its_own() {
    let text = two_fluids()
        .replace(
            "[regions.water]\nkind = \"fluid\"\nzones = [\"water\"]",
            "[regions.water]\nkind = \"fluid\"\nzones = [\"water\"]\n[regions.water.initial]\ntemperature = \"inlet\"",
        )
        .replace(
            "material = \"light-oil\"\n",
            "material = \"light-oil\"\n[regions.oil.initial]\ntemperature = \"inlet\"\n",
        )
        .replace(
            "material = \"plate-steel\"\n",
            "material = \"plate-steel\"\n[regions.plate.initial]\ntemperature = 320.0\n",
        )
        .replace("case_format = \"0.7\"", "case_format = \"0.7\"\n")
        .replace("pressure = 1.0e5 ", "pressure = 0.0 ")
        .replace("pressure = 1.5e5 ", "pressure = 5.0e4 ")
        .replace("[fluid]\n", "[fluid]\npressure_values = \"gauge\"\noperating_pressure = 1.0e5\n");
    let c = FlowCase::from_toml(&text, &cases()).unwrap();
    let d = out("two");
    run(&c, "two", &d, |_| {}).unwrap();
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(d.join("flow.json")).unwrap()).unwrap();
    let sv = j["start_values"].as_array().expect("start_values");
    let find = |scope: &str, q: &str| {
        sv.iter()
            .find(|s| s["scope"] == scope && s["quantity"].as_str().unwrap().starts_with(q))
            .unwrap_or_else(|| panic!("{scope} {q}: {sv:?}"))
            .clone()
    };
    assert_eq!(find("region water", "temperature")["value"], "340");
    assert_eq!(find("region oil", "temperature")["value"], "300");
    assert_eq!(find("region plate", "temperature")["value"], "320");
    let oil_p = sv
        .iter()
        .find(|s| s["quantity"].as_str().unwrap().starts_with("pressure"))
        .unwrap_or_else(|| panic!("{sv:?}"));
    assert_eq!(oil_p["value"], "50000");
    assert_eq!(oil_p["source"], "boundary oil_out");
    let log = std::fs::read_to_string(d.join("flow.log")).unwrap();
    assert!(
        log.contains("[REG010] start values, region water: temperature 340"),
        "{log}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

/// REG007: a solid region takes a start temperature only; REG009: an
/// `"inlet"` temperature needs an inlet with a temperature.
#[test]
fn region_tables_are_checked() {
    let text = two_fluids().replace(
        "material = \"plate-steel\"\n",
        "material = \"plate-steel\"\n[regions.plate.initial]\nvelocity = [1.0, 0.0, 0.0]\n",
    );
    let e = FlowCase::from_toml(&text, &cases()).unwrap_err();
    assert_eq!(diag::of_error(&e).code, "REG007", "{e:#}");
    let text = two_fluids()
        .replace(
            "material = \"light-oil\"\n",
            "material = \"light-oil\"\n[regions.oil.initial]\ntemperature = \"inlet\"\n",
        )
        .replace(
            "[boundary.oil_in]\ntype = \"velocity-inlet\"\nvelocity = [\"-0.036759001803*((y - 0.0025)/0.002)*(1 - (y - 0.0025)/0.002)\", \"0\", \"0\"]",
            "[boundary.oil_in]\ntype = \"wall\"",
        );
    let c = FlowCase::from_toml(&text, &cases());
    let e = match c {
        Err(e) => e,
        Ok(c) => run(&c, "inlet", &out("inlet"), |_| {}).unwrap_err(),
    };
    let _ = std::fs::remove_dir_all(out("inlet"));
    assert_eq!(diag::of_error(&e).code, "REG009", "{e:#}");
}
