//! Datums of a flow run (SP1b, D262, V28.13): the gas through a plane
//! moved to the faces and through the outlet boundary closes with the
//! boundary flows, a plane cut through the cells stays within its quoted
//! cut error, the outlet datum is the boundary's own sum, and the liquid
//! through the outlet datum is the mass book's escaped liquid.
#![cfg(feature = "hdf5")]

use foehncraft_flow::case::{FlowCase, run};
use std::path::{Path, PathBuf};

fn out(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("foehncraft_datums_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

const CASE: &str = r#"
[mesh.generate]
type = "hex"
dimension = 2
upper = [2.0, 0.4, 0.1]
cells = [40, 10, 1]
perturb = 0.15
[fluid]
density = 1.2
viscosity = 1.8e-5
[gravity]
vector = [0.0, -9.81, 0.0]
[initial]
velocity = [4.0, 0.0, 0.0]
[boundary.xmin]
type = "velocity-inlet"
velocity = [4.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[time]
time_step = 0.01
end_time = 0.6
[spray]
wall_mode = "stick"
[spray.liquid]
density = 998.0
[[spray.injector]]
name = "nozzle"
position = [0.1, 0.39, 0.05]
direction = [0.5, -1.0, 0.0]
speed = 5.0
half_angle = 10.0
mass_flow = 0.01
parcels_per_second = 2000.0
size = { distribution = "log-normal", median = 80e-6, sigma = 1.8, min = 10e-6, max = 300e-6 }
[[datum]]
name = "cut"
type = "plane"
point = [0.63, 0.2, 0.05]
normal = [1.0, 0.0, 0.0]
[[datum]]
name = "snapped"
type = "plane"
point = [0.63, 0.2, 0.05]
normal = [1.0, 0.0, 0.0]
snap = "faces"
[[datum]]
name = "outlet"
type = "boundary"
boundary = "xmax"
[[datum]]
name = "inlet"
type = "boundary"
boundary = "xmin"
[[datum]]
name = "middle"
type = "cylinder"
start = [0.8, 0.2, 0.0]
end = [0.8, 0.2, 0.1]
radius = 0.15
[output]
vtkhdf = false
wall_csv = false
monitor_every = 5
[solver]
tolerance = 1e-8
max_iterations = 200
[verification]
require_convergence = false
"#;

#[test]
fn datums_close_with_the_boundaries() {
    let case = FlowCase::from_toml(CASE, Path::new(".")).unwrap();
    let d = out("close");
    let r = run(&case, "datums", &d, |_| {}).unwrap();
    let s = r.solves.last().unwrap();
    let ds = s.datums.as_ref().expect("datums in the result");
    let get = |n: &str| ds.iter().find(|x| x.name == n).unwrap();
    let q = |n: &str, k: &str| get(n).gas[k];
    let outlet = s.boundaries.iter().find(|b| b.name == "xmax").unwrap();
    let inlet = s.boundaries.iter().find(|b| b.name == "xmin").unwrap();
    // the boundary datum is the boundary's own sum, bit for bit
    assert_eq!(q("outlet", "mass_flow"), outlet.mass_flow);
    assert_eq!(q("inlet", "mass_flow"), inlet.mass_flow);
    assert!(get("outlet").exact && get("snapped").exact && !get("cut").exact);
    // the face layer carries the inlet flow to the solver's tolerance
    let m_in = -inlet.mass_flow;
    let snapped = q("snapped", "mass_flow");
    assert!(
        (snapped - m_in).abs() < 1e-6 * m_in,
        "snapped {snapped} against the inlet {m_in}"
    );
    // the cut plane within its quoted error
    let cut = q("cut", "mass_flow");
    let e = q("cut", "cut_error");
    assert!(e < 1e-3, "cut error {e}");
    assert!(((cut - snapped) / snapped).abs() <= e * (1.0 + 1e-9) + 1e-15);
    assert!(q("middle", "volume") > 0.0);
    // the liquid through the outlet datum is the escaped liquid of the book
    let sp = s.spray.as_ref().expect("spray summary");
    let escaped: f64 = sp.tally.escaped_mass.iter().sum();
    let liq = get("outlet").liquid.as_ref().unwrap();
    assert_eq!(liq.mass_pos, Some(escaped));
    // and the files are there
    for f in [
        "datums.csv",
        "spray_hist_cut.csv",
        "series.json",
        "spray.csv",
    ] {
        assert!(d.join(f).exists(), "{f} missing");
    }
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn bad_datums_are_refused_with_their_code() {
    for (extra, what) in [
        ("[[datum]]\nname = \"a\"\ntype = \"blob\"\n", "unknown type"),
        (
            "[[datum]]\nname = \"a\"\ntype = \"boundary\"\nboundary = \"nowhere\"\n",
            "no boundary",
        ),
        (
            "[[datum]]\nname = \"a\"\ntype = \"plane\"\npoint = [5.0, 0.2, 0.05]\nnormal = [1.0, 0.0, 0.0]\nsnap = \"faces\"\n",
            "meets no face",
        ),
    ] {
        let text = CASE.replace("[output]", &format!("{extra}[output]"));
        let d = out("bad");
        let e = FlowCase::from_toml(&text, Path::new("."))
            .and_then(|c| {
                c.validate()?;
                run(&c, "bad", &d, |_| {})
            })
            .err()
            .unwrap_or_else(|| panic!("{what}: accepted"));
        let m = format!("{e:#}");
        assert!(m.contains(what), "{what}: {m}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
