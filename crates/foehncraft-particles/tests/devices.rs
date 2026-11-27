//! Injector devices (D340 to D342, checks VD.1, VD.2, VD.5 of the step 28
//! injector device page).

use foehncraft_particles::SpraySection;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn write(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("foehn_devices_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

/// VD.1: every injector file of the examples reads as before: positions
/// and directions are those written in the file, nothing is added.
#[test]
fn existing_injector_files_are_unchanged() {
    let mut n = 0;
    for e in std::fs::read_dir(cases()).unwrap() {
        let p = e.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        if !name.ends_with("_injectors.toml") || name.contains("device") {
            continue;
        }
        let raw: toml::Table = toml::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
        let sp = SpraySection::from_file(&p).unwrap_or_else(|e| panic!("{name}: {e:#}"));
        let list = raw["injector"].as_array().unwrap();
        assert_eq!(sp.injectors.len(), list.len(), "{name}");
        assert!(
            sp.devices.is_empty() && sp.device_warnings.is_empty(),
            "{name}"
        );
        for (i, t) in sp.injectors.iter().zip(list) {
            let v = |k: &str| -> Vec<f64> {
                t[k].as_array()
                    .unwrap()
                    .iter()
                    .map(|x| {
                        x.as_float()
                            .unwrap_or_else(|| x.as_integer().unwrap() as f64)
                    })
                    .collect()
            };
            assert_eq!(i.position.to_vec(), v("position"), "{name}");
            assert_eq!(i.direction.to_vec(), v("direction"), "{name}");
            assert!(i.device.is_none() && i.hole.is_none());
        }
        n += 1;
    }
    assert!(n >= 10, "only {n} injector files found");
}

const DEVICE: &str = r#"
[liquid]
name = "water"
density = 998.0

[[device]]
name = "g"
mount = [0.0, 0.0, 0.0]
axis = [1.0, 0.0, 0.0]
reference = [0.0, 1.0, 0.0]

[[device.pattern]]
count = 8
polar = 33.0
radius = 1.0e-3
setback = 1.0e-4
first_azimuth = -11.25
"#;

const KEYS: &str = r#"speed = 10.0
mass_flow = 1e-3
diameter = 1e-5
parcels_per_second = 1000.0
"#;

/// VD.2: Spray G hole 1 from the device, with the 2 µm nudge into the
/// sector as `offset`, is the hand-written injector to round-off; the file
/// itself, given the device, gives no warning.
#[test]
fn spray_g_hole_one_from_the_device() {
    let a = (11.25f64).to_radians();
    let off = [0.0, 2e-6 * a.sin(), 2e-6 * a.cos()];
    let text = format!(
        "{DEVICE}\n[[injector]]\nname = \"h\"\ndevice = \"g\"\nhole = \"1\"\noffset = [{:e}, {:e}, {:e}]\n{KEYS}",
        off[0], off[1], off[2]
    );
    let sp = SpraySection::from_file(&write("g1.toml", &text)).unwrap();
    let i = &sp.injectors[0];
    let pos = [1.000000e-04, 9.81175461e-04, -1.93128751e-04];
    let dir = [0.83867057, 0.53417395, -0.10625380];
    for k in 0..3 {
        assert!((i.position[k] - pos[k]).abs() < 1e-11, "{:?}", i.position);
        assert!((i.direction[k] - dir[k]).abs() < 1e-8, "{:?}", i.direction);
    }
    // the hand-written values with the device named: no warning
    let text = format!(
        "{DEVICE}\n[[injector]]\nname = \"h\"\ndevice = \"g\"\nhole = \"1\"\n\
         position = [1.000000e-04, 9.81175461e-04, -1.93128751e-04]\n\
         direction = [0.83867057, 0.53417395, -0.10625380]\n{KEYS}"
    );
    let sp = SpraySection::from_file(&write("g1_abs.toml", &text)).unwrap();
    assert!(sp.device_warnings.is_empty(), "{:?}", sp.device_warnings);
    assert_eq!(sp.injectors[0].position, pos);
    // the same values named against hole 2: SPR010
    let sp = SpraySection::from_file(&write(
        "g1_wrong.toml",
        &text.replace("hole = \"1\"", "hole = \"2\""),
    ))
    .unwrap();
    assert_eq!(sp.device_warnings.len(), 1);
    assert_eq!(sp.injectors[0].position, pos);
}

/// VD.5: `hole = "all"` is the list of one injector per hole, in the
/// order of the holes, named `<name>/<hole>`.
#[test]
fn all_holes_expand_in_order() {
    let all =
        format!("{DEVICE}\n[[injector]]\nname = \"s\"\ndevice = \"g\"\nhole = \"all\"\n{KEYS}");
    let a = SpraySection::from_file(&write("all.toml", &all)).unwrap();
    assert_eq!(a.injectors.len(), 8);
    let mut one = DEVICE.to_string();
    for k in 1..=8 {
        one.push_str(&format!(
            "\n[[injector]]\nname = \"s/{k}\"\ndevice = \"g\"\nhole = \"{k}\"\n{KEYS}"
        ));
    }
    let b = SpraySection::from_file(&write("one.toml", &one)).unwrap();
    assert_eq!(a.injectors, b.injectors);
    assert_eq!(a.injectors[3].name, "s/4");
}

#[test]
fn refusals() {
    let bad = |name: &str, inj: &str| {
        let t = format!("{DEVICE}\n[[injector]]\nname = \"x\"\n{inj}\n{KEYS}");
        let e = SpraySection::from_file(&write(name, &t)).unwrap_err();
        format!("{e:#}")
    };
    assert!(bad("r1.toml", "device = \"g\"").contains("needs hole"));
    assert!(bad("r2.toml", "device = \"g\"\nhole = \"9\"").contains("no hole 9"));
    assert!(bad("r3.toml", "device = \"q\"\nhole = \"1\"").contains("not a [[spray.device]]"));
    assert!(
        bad(
            "r4.toml",
            "hole = \"1\"\nposition = [0.0, 0.0, 0.0]\ndirection = [1.0, 0.0, 0.0]"
        )
        .contains("need device")
    );
    assert!(bad("r5.toml", "direction = [1.0, 0.0, 0.0]").contains("position is required"));
    assert!(
        bad(
            "r6.toml",
            "device = \"g\"\nhole = \"all\"\nposition = [0.0, 0.0, 0.0]\ndirection = [1.0, 0.0, 0.0]"
        )
        .contains("leave them out")
    );
}

/// The Spray G device file gives the injector of the hand-written file to
/// 1e-11 m and 1e-8 in the direction (D342), with the standard tip.
#[test]
fn spray_g_device_file_matches_the_hand_written_one() {
    let a = SpraySection::from_file(&cases().join("spray_g_device_injectors.toml")).unwrap();
    let b = SpraySection::from_file(&cases().join("spray_g_g1_injectors.toml")).unwrap();
    assert_eq!(a.injectors.len(), 1);
    let (i, j) = (&a.injectors[0], &b.injectors[0]);
    for k in 0..3 {
        assert!((i.position[k] - j.position[k]).abs() < 1e-11);
        assert!((i.direction[k] - j.direction[k]).abs() < 1e-8);
    }
    assert!(a.device_warnings.is_empty());
    let mut ii = i.clone();
    ii.position = j.position;
    ii.direction = j.direction;
    ii.device = None;
    ii.hole = None;
    ii.offset = None;
    assert_eq!(&ii, j, "every other key the same");
}

