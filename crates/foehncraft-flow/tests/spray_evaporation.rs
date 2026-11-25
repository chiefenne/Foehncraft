//! Evaporating sprays (SP4, SP5; checks of
//! `docs/requirements/step28-sprays.md`).
//!
//! - The books of a transient evaporating duct: the ledger pairs (vapour,
//!   momentum, heat) are zero to the bit, the parcel accounts (mass,
//!   momentum, enthalpy) and the mass account of the gas close to
//!   round-off, and the vapour the gas receives is the liquid the parcels
//!   lose (D244, D254 (1)).
//! - V28.1 with evaporation: the same bits on 1 and 4 threads.
//! - The setting check: the vapour the liquid names must be a species of
//!   the case.

use foehncraft_flow::case::{FlowCase, run};
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "foehncraft_evaporation_{name}_{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn examples() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

/// The evaporating duct of the examples, coarser and shorter.
fn duct(end: f64) -> String {
    std::fs::read_to_string(examples().join("spray_evaporation.toml"))
        .unwrap()
        .replace("cells = [100, 20, 1]", "cells = [40, 8, 1]")
        .replace("end_time = 0.4", &format!("end_time = {end}"))
        .replace("write_every = 50", "write_every = 0")
        .replace(
            "[solver.scalar]",
            "[output]\nvtkhdf = false\nwall_csv = false\n\n[solver.scalar]",
        )
}

struct Books {
    pairs: Vec<(u32, f64)>,
    /// (quantity, region, relative imbalance)
    accounts: Vec<(String, String, f64)>,
    gas_vapour: f64,
    evaporated: f64,
    csv: String,
}

fn go(name: &str, text: &str) -> Books {
    let d = out(name);
    let case = FlowCase::from_toml(text, &examples()).unwrap();
    let r = run(&case, name, &d, |_| {}).unwrap();
    let s = r.solves.last().unwrap();
    let bal = s.balance.as_ref().expect("balance");
    let run = bal.run.as_ref().expect("run ledger");
    let pairs = run.pairs.iter().map(|p| (p.pair, p.residual)).collect();
    let accounts = run
        .ledger
        .iter()
        .map(|a| (a.quantity.clone(), a.region.clone(), a.relative))
        .collect();
    let term = |q: &str, t: &str| {
        run.ledger
            .iter()
            .find(|a| a.quantity == q && a.region == "parcels")
            .and_then(|a| a.terms.iter().find(|x| x.name == t))
            .map_or(f64::NAN, |x| x.value)
    };
    let evaporated = -term("parcel mass", "to the gas");
    let gas_vapour = run
        .ledger
        .iter()
        .find(|a| a.quantity == "mass" && a.region == "domain")
        .and_then(|a| a.terms.iter().find(|x| x.name == "spray vapour"))
        .map_or(f64::NAN, |x| x.value);
    let csv = std::fs::read_to_string(d.join("spray.csv"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n");
    let _ = std::fs::remove_dir_all(&d);
    Books {
        pairs,
        accounts,
        gas_vapour,
        evaporated,
        csv,
    }
}

#[test]
fn the_evaporating_duct_closes_its_books() {
    let b = go("books", &duct(0.04));
    // vapour 2805, momentum 2801 to 2803, heat 2804
    for id in [2801, 2802, 2803, 2804, 2805] {
        let r = b.pairs.iter().find(|p| p.0 == id).map(|p| p.1);
        assert_eq!(r, Some(0.0), "pair {id}: {r:?}");
    }
    for q in [
        "parcel mass",
        "parcel momentum x",
        "parcel momentum y",
        "parcel enthalpy",
    ] {
        let a = b
            .accounts
            .iter()
            .find(|a| a.0 == q && a.1 == "parcels")
            .unwrap_or_else(|| panic!("no account {q}"));
        assert!(a.2.abs() < 1e-12, "{q}: relative imbalance {:e}", a.2);
    }
    let m = b
        .accounts
        .iter()
        .find(|a| a.0 == "mass" && a.1 == "domain")
        .unwrap();
    assert!(m.2.abs() < 1e-9, "gas mass: relative imbalance {:e}", m.2);
    assert!(b.evaporated > 0.0, "nothing evaporated");
    assert_eq!(b.gas_vapour, b.evaporated);
}

#[test]
fn evaporation_gives_the_same_bits_on_any_number_of_threads() {
    let text = duct(0.02);
    let runs: Vec<String> = [1, 4]
        .iter()
        .map(|&n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| go(&format!("threads{n}"), &text).csv)
        })
        .collect();
    assert!(runs[0].lines().count() > 5);
    assert_eq!(runs[0], runs[1]);
}

#[test]
fn evaporation_needs_a_vapour_species_of_the_case() {
    // the liquid names a vapour the case does not solve
    let text = duct(0.02).replace(
        "names = [\"H2O\", \"O2\", \"N2\"]",
        "names = [\"CO2\", \"O2\", \"N2\"]",
    );
    let text = text.replace("H2O = 0.005", "CO2 = 0.005");
    let e = FlowCase::from_toml(&text, &examples())
        .and_then(|c| run(&c, "novapour", &out("novapour"), |_| {}).map(|_| ()))
        .expect_err("a case without the vapour species was accepted");
    assert!(e.to_string().contains("H2O"), "{e}");
}
