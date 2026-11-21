//! Stage 7, PR 7d1: solid regions in the energy equation, conjugate heat
//! transfer (D112, D157), docs/requirements/stage7-energy-buoyancy-cht.md,
//! check V7.4.

use foehncraft_flow::case::{FlowCase, FlowReport, run};
use foehncraft_flow::diag;
use std::path::{Path, PathBuf};

fn cases() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cases/flow")
}

fn out(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("foehncraft_cht_{name}_{}", std::process::id()))
}

fn run_case(case: &FlowCase, name: &str) -> FlowReport {
    let dir = out(name);
    let rep = run(case, name, &dir, |_| {}).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    rep
}

fn parse(text: &str) -> FlowCase {
    FlowCase::from_toml(text, &cases()).unwrap()
}

fn text(file: &str) -> String {
    std::fs::read_to_string(cases().join(file))
        .unwrap()
        .replace("\r\n", "\n")
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

/// V7.4: a composite wall of two solids and a fluid layer at rest,
/// conductivity ratios 1 : 100 : 10⁴. The series heat flux and the
/// interface temperatures to round-off; the same flux at both walls and
/// through both interfaces.
#[test]
fn composite_wall_gives_the_series_heat_flux() {
    let r = run_case(&parse(&text("cht_composite_wall.toml")), "wall");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged);
    let sc = &s.scalars[0];
    let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
    let q = 100.0 / (0.4 / 200.0 + 0.4 / 2.0 + 0.4 / 20000.0);
    for (name, sign) in [
        ("ymin", -1.0),
        ("ymax", 1.0),
        ("interface-a", -1.0),
        ("interface-b", 1.0),
    ] {
        let f = b(name).flux;
        assert!(
            (f - sign * q).abs() < 1e-9 * q,
            "{name}: {f} against {}",
            sign * q
        );
    }
    let ta = 300.0 + q * 0.4 / 200.0;
    let tb = 400.0 - q * 0.4 / 20000.0;
    assert!(
        (b("interface-a").mean - ta).abs() < 1e-9,
        "{}",
        b("interface-a").mean
    );
    assert!(
        (b("interface-b").mean - tb).abs() < 1e-9,
        "{}",
        b("interface-b").mean
    );
    let e = sc.energy.as_ref().unwrap();
    assert!(e.relative < 1e-9, "{e:?}");
}

/// The heated wall under a channel flow (the gallery case, on a coarser
/// mesh): all of the heater's power crosses the interface into the fluid
/// and leaves with the flow; the balance of the whole domain closes.
#[test]
fn heated_wall_passes_its_power_to_the_fluid() {
    let t = text("cht_heated_wall.toml").replace("cells = [80, 30, 1]", "cells = [40, 15, 1]");
    let r = run_case(&parse(&t), "heated_wall");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = &r.solves[0];
    assert!(s.converged);
    let sc = &s.scalars[0];
    let b = |n: &str| sc.boundaries.iter().find(|b| b.boundary == n).unwrap();
    // 1000 W/m² over 4 m × 0.04 m
    let p = 160.0;
    assert!(
        (b("ymin").heat_flow.unwrap() - p).abs() < 1e-6,
        "{:?}",
        b("ymin")
    );
    assert!(
        (b("interface-wall").heat_flow.unwrap() - p).abs() < 1e-6 * p,
        "{:?}",
        b("interface-wall")
    );
    let e = sc.energy.as_ref().unwrap();
    assert!(e.relative < 1e-8, "{e:?}");
    // the wall is warmer than the fluid at the inlet, cooler than the
    // heated face
    assert!(b("interface-wall").mean > 300.0 && b("interface-wall").mean < b("ymin").mean);
}

/// The diagnostics of solid regions.
#[test]
fn solid_regions_are_checked() {
    let base = text("cht_heated_wall.toml");
    // a solid without a material
    let bare = base.replace("material = \"wall-solid\"\n", "");
    let (code, msg) = refused(&bare);
    assert_eq!(code, "MAT006", "{msg}");
    // a solid material without a conductivity
    let no_k = base.replace("conductivity = 140.84507042253522\n", "");
    let (code, msg) = refused(&no_k);
    assert_eq!(code, "MAT006", "{msg}");
    // a zone the mesh does not have
    let zone = base.replace("zones = [\"wall\"]", "zones = [\"wal\"]");
    let (code, msg) = refused(&zone);
    assert_eq!(code, "REG002", "{msg}");
    // solids without the energy equation
    let (a, b) = (
        base.find("[energy]").unwrap(),
        base.find("[boundary.inlet]").unwrap(),
    );
    let flow_only = format!("{}{}", &base[..a], &base[b..])
        .replace("[boundary.inlet.energy]\ntemperature = 300.0\n", "")
        .replace(
            "[boundary.ymin.energy]\nheat_flux = 1000.0                   # W/m² into the wall\n",
            "",
        );
    let (code, msg) = refused(&flow_only);
    assert_eq!(code, "REG001", "{msg}");
    // transient conjugate heat transfer is a later step
    let transient = base + "\n[time]\ntime_step = 0.01\nend_time = 0.1\n";
    let (code, msg) = refused(&transient);
    assert_eq!(code, "ENE006", "{msg}");
}

/// The turbulent channel over a conducting wall as the fluid alone: the
/// mesh of the fluid cells, the wall with the convection condition
/// h = k_s/t = 2 W/(m² K) to the 350 K of the wall's outer face, in series
/// with the thermal wall function of the fluid-only code (the convection
/// condition takes the final face conductance, D135).
fn turbulent_reference(base: &str) -> String {
    let mut s = base.replace(
        "upper = [20.0, 1.1, 0.1]\ncells = [40, 33, 1]",
        "lower = [0.0, 0.1, 0.0]\nupper = [20.0, 1.1, 0.1]\ncells = [40, 30, 1]",
    );
    // the text from `from` up to the next `to` (kept) replaced by `with`
    let cut = |s: &str, from: &str, to: &str, with: &str| {
        let a = s.find(from).unwrap();
        let b = a + s[a..].find(to).unwrap();
        format!("{}{with}{}", &s[..a], &s[b..])
    };
    s = cut(&s, "split = [", "[mesh]\n", "\n");
    s = cut(&s, "[materials.wall-solid]", "[turbulence]", "");
    s = cut(
        &s,
        "[boundary.ymin]",
        "[boundary.zmin]",
        "[boundary.ymin]\ntype = \"wall\"\n[boundary.ymin.energy]\n\
         heat_transfer_coefficient = 2.0\nambient_temperature = 350.0\n",
    );
    s = cut(&s, "lines = [", "\n]\n", "");
    s.replace("\n]\n", "\n")
        .replace("[boundary.inlet]", "[boundary.xmin]")
        .replace("[boundary.inlet.energy]", "[boundary.xmin.energy]")
        .replace("[boundary.outlet]", "[boundary.xmax]")
}

/// V7.36 (stage 7d5, specification §9.4): a turbulent channel flow over a
/// conducting wall, k-ε with scalable wall functions and SST with the
/// automatic treatment. The fluid side of the interface is the thermal
/// wall function in series with the solid half cell. Against the fluid
/// alone with the wall function in series with h = k_s/t (an independent
/// path of the fluid-only code): the same heat flow, interface temperature
/// and outlet temperature (the solid's axial conduction, absent in the
/// reference, is about 10⁻⁴ of the heat flow); the temperature drop
/// across the wall is q t/k_s; every watt of the outer face crosses the
/// interface.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn turbulent_interface_is_the_wall_function_in_series() {
    let base = text("cht_turbulent_channel.toml");
    for (name, model) in [
        ("keps", ""),
        (
            "sst",
            "model = \"k-omega-sst\"\nwall_treatment = \"automatic\"\n",
        ),
    ] {
        let set = |t: &str| {
            if model.is_empty() {
                t.to_string()
            } else {
                t.replace(
                    "model = \"k-epsilon\"\nwall_treatment = \"scalable\"\n",
                    model,
                )
            }
        };
        let a = run_case(&parse(&set(&base)), &format!("turb_{name}"));
        let b = run_case(
            &parse(&set(&turbulent_reference(&base))),
            &format!("turb_{name}_ref"),
        );
        for r in [&a, &b] {
            assert!(r.failures.is_empty(), "{name}: {:?}", r.failures);
            assert!(r.solves[0].converged, "{name}");
        }
        let sa = &a.solves[0].scalars[0];
        let sb = &b.solves[0].scalars[0];
        let ba = |n: &str| sa.boundaries.iter().find(|b| b.boundary == n).unwrap();
        let bb = |n: &str| sb.boundaries.iter().find(|b| b.boundary == n).unwrap();
        let qa = ba("interface-wall").heat_flow.unwrap();
        let qo = ba("ymin").heat_flow.unwrap();
        let qb = bb("ymin").heat_flow.unwrap();
        eprintln!(
            "{name}: interface {qa:.6} W (outer face {qo:.6} W), fluid alone {qb:.6} W; \
             T_i {:.4} K against {:.4} K",
            ba("interface-wall").mean,
            bb("ymin").mean
        );
        assert!((qa - qo).abs() < 1e-6 * qo, "{name}: {qa} against {qo}");
        assert!((qa / qb - 1.0).abs() < 1e-3, "{name}: {qa} against {qb}");
        let ti = ba("interface-wall").mean;
        assert!((ti - bb("ymin").mean).abs() < 0.02, "{name}: T_i {ti}");
        assert!(
            (ba("outlet").mixed_mean.unwrap() - bb("xmax").mixed_mean.unwrap()).abs() < 1e-3,
            "{name}: outlet"
        );
        // the 1D wall: T_o − T_i = q t/k_s (q per area of the 20 m × 0.1 m face)
        let q = qo / 2.0;
        let drop = 350.0 - ti;
        assert!(
            (drop / (q * 0.1 / 0.2) - 1.0).abs() < 1e-3,
            "{name}: {drop} K against {}",
            q * 0.5
        );
        let e = sa.energy.as_ref().unwrap();
        assert!(e.relative < 1e-8, "{name}: {e:?}");
    }
}
