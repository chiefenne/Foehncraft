//! Solve a circuit case and write `result.json` and `profiles.csv` (the
//! same files as the tool `circuits.run`, without the registry):
//!
//! ```sh
//! cargo run --release -p foehncraft-circuits --example circuit_run -- CASE OUT_DIR
//! ```

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: circuit_run CASE OUT_DIR");
        std::process::exit(2);
    }
    let path = std::path::Path::new(&args[1]);
    let case = foehncraft_circuits::case::Case::from_file(path).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let sol = foehncraft_circuits::output::run(&case).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let name = path
        .file_stem()
        .map_or("circuit".into(), |s| s.to_string_lossy().into_owned());
    let out = std::path::Path::new(&args[2]);
    std::fs::create_dir_all(out).expect("output folder");
    let json = sol.to_json(&name, &case);
    std::fs::write(
        out.join("result.json"),
        serde_json::to_string_pretty(&json).expect("json") + "\n",
    )
    .expect("write");
    std::fs::write(out.join("profiles.csv"), sol.to_csv()).expect("write");
    for w in &sol.warnings {
        eprintln!("{}: {}", w.code, w.message);
    }
    eprintln!("{name}: {} Newton iterations", sol.history.len() - 1);
}
