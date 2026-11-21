//! Stage 7, PRs 7d1, 7d5 and 7d6: solid regions in the energy equation,
//! conjugate heat transfer (D112, D157), steady, turbulent and transient,
//! docs/requirements/stage7-energy-buoyancy-cht.md, checks V7.4, V7.6,
//! V7.36 and V7.37.

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
    // passive scalars would diffuse into the solids
    let dye = base + "\n[scalars.dye]\ndiffusivity = 1e-3\ninitial = 0.0\n";
    let (code, msg) = refused(&dye);
    assert_eq!(code, "ENE006", "{msg}");
}

/// The contact temperature of two semi-infinite bodies (Carslaw & Jaeger
/// 1959, §2.15).
fn contact_temperature(k: [f64; 2], rho_cp: [f64; 2], t: [f64; 2]) -> f64 {
    let e = [(k[0] * rho_cp[0]).sqrt(), (k[1] * rho_cp[1]).sqrt()];
    (e[0] * t[0] + e[1] * t[1]) / (e[0] + e[1])
}

/// The temperature summary of a report and the interface boundary of it.
fn temperature(r: &FlowReport) -> &foehncraft_flow::scalar::ScalarSummary {
    r.solves[0]
        .scalars
        .iter()
        .find(|s| s.name == "temperature")
        .unwrap()
}

fn boundary<'a>(
    s: &'a foehncraft_flow::scalar::ScalarSummary,
    name: &str,
) -> &'a foehncraft_flow::scalar::ScalarBoundarySummary {
    s.boundaries.iter().find(|b| b.boundary == name).unwrap()
}

/// V7.6: two semi-infinite bodies brought into contact (the solid at
/// 400 K, k = 10, the fluid at rest at 300 K, k = 1, both ρ c_p = 10⁶).
/// The interface temperature is T_i from the first steps on; at 500 s the
/// L2 error against the erf profiles of `[energy] exact` is below 0.01 K
/// (in the case), the heat flow through the interface is
/// e₁ (T₁ − T_i)/√(π t) per area, and the heat the solid loses is what the
/// fluid stores, in every step.
#[test]
fn contact_of_two_bodies_keeps_the_contact_temperature() {
    let base = text("cht_contact_transient.toml");
    let ti = contact_temperature([10.0, 1.0], [1e6, 1e6], [400.0, 300.0]);
    // after the first 10 steps: within 0.1 % of ΔT (the profiles are only
    // a few cells deep, the L2 error is not checked yet)
    let early = base
        .replace("end_time = 500.0", "end_time = 10.0")
        .replace("scalar_max_error = 0.01 ", "scalar_max_error = 1.0 ");
    let r = run_case(&parse(&early), "contact_early");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let t = boundary(temperature(&r), "interface-bar").mean;
    assert!((t - ti).abs() < 1e-3 * 100.0, "T_i {t} against {ti}");
    // at 500 s
    let r = run_case(&parse(&base), "contact");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = temperature(&r);
    let b = boundary(s, "interface-bar");
    assert!((b.mean - ti).abs() < 1e-6, "T_i {} against {ti}", b.mean);
    let q = (10.0f64 * 1e6).sqrt() * (400.0 - ti) / (std::f64::consts::PI * 500.0).sqrt() * 1e-4;
    let qi = b.heat_flow.unwrap();
    assert!(
        (qi / q - 1.0).abs() < 2e-3,
        "interface {qi} W against {q} W"
    );
    assert!(s.l2.unwrap() < 6e-3, "L2 error {:?}", s.l2);
    let e = s.energy.as_ref().unwrap();
    let st: Vec<f64> = e.storage_regions.iter().map(|r| r.storage).collect();
    assert_eq!(st.len(), 2, "{e:?}");
    assert!((st[0] + st[1]).abs() < 1e-9 * st[0].abs(), "{e:?}");
    assert!((st[0] / q - 1.0).abs() < 2e-3, "{e:?}");
    let h = r.solves[0].transient.as_ref().unwrap();
    assert!(h.heat_relative_max.unwrap() < 1e-9, "{h:?}");
}

/// V7.6 variants: copper against air (the contact temperature stays at
/// that of the copper), and two equal bodies (T_i the mean, by symmetry
/// to round-off).
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn contact_of_copper_and_air_and_of_equal_bodies() {
    let base = text("cht_contact_transient.toml");
    // copper and air at 300 K (Incropera et al. 2007, tables A.1 and A.4)
    let cu = (401.0, 8933.0, 385.0);
    let air = (0.0263, 1.1614, 1007.0, 1.846e-5);
    let copper = base
        .replace("k1 = 10.0 ", "k1 = 401.0 ")
        .replace("k2 = 1.0\n", "k2 = 0.0263\n")
        .replace("c1 = 1.0e6 ", &format!("c1 = {} ", cu.1 * cu.2))
        .replace("c2 = 1.0e6\n", &format!("c2 = {}\n", air.1 * air.2))
        .replace(
            "density = 1000.0\nviscosity = 1.0e-3\nspecific_heat = 1000.0\nconductivity = 1.0\n",
            &format!(
                "density = {}\nviscosity = {}\nspecific_heat = {}\nconductivity = {}\n",
                air.1, air.3, air.2, air.0
            ),
        )
        .replace(
            "density = 1000.0\nspecific_heat = 1000.0\nconductivity = 10.0\n",
            &format!(
                "density = {}\nspecific_heat = {}\nconductivity = {}\n",
                cu.1, cu.2, cu.0
            ),
        )
        .replace("end_time = 500.0", "end_time = 100.0");
    assert!(
        copper.contains("conductivity = 401"),
        "the copper replacements"
    );
    assert!(
        copper.contains("conductivity = 0.0263"),
        "the air replacements"
    );
    let ti = contact_temperature([cu.0, air.0], [cu.1 * cu.2, air.1 * air.2], [400.0, 300.0]);
    let r = run_case(&parse(&copper), "contact_copper");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let t = boundary(temperature(&r), "interface-bar").mean;
    assert!(
        (t - ti).abs() < 1e-3 * 100.0,
        "copper: T_i {t} against {ti}"
    );
    assert!(ti > 399.98, "{ti}");
    // equal bodies
    let equal = base
        .replace("k1 = 10.0 ", "k1 = 1.0 ")
        .replace("conductivity = 10.0\n", "conductivity = 1.0\n");
    let r = run_case(&parse(&equal), "contact_equal");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let t = boundary(temperature(&r), "interface-bar").mean;
    assert!((t - 350.0).abs() < 1e-9, "equal bodies: T_i {t}");
}

/// V7.6 orders: the error against the exact solution over three meshes
/// (second order in space), the differences of successive solutions over
/// four time steps (BDF2, second order in time); the case checks both
/// against [1.8, 2.2].
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn contact_converges_at_second_order_in_space_and_time() {
    let base = text("cht_contact_transient.toml");
    let space = base
        .replace("cells = [400, 1, 1]", "cells = [100, 1, 1]")
        .replace("[verification]\n", "[verification]\nlevels = [0, 1, 2]\n");
    let r = run_case(&parse(&space), "contact_space");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let (_, o) = r
        .scalar_orders
        .iter()
        .find(|(n, _)| n == "temperature")
        .expect("space orders");
    assert!(
        o.len() == 2 && o.iter().all(|o| (1.8..=2.2).contains(o)),
        "space {o:?}"
    );
    let time = base.replace(
        "[verification]\n",
        "[verification]\ntime_steps = [20.0, 10.0, 5.0, 2.5]\n",
    );
    let r = run_case(&parse(&time), "contact_time");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let (_, o) = r
        .scalar_orders
        .iter()
        .find(|(n, _)| n == "temperature (time)")
        .expect("time orders");
    assert!(
        o.len() == 2 && o.iter().all(|o| (1.8..=2.2).contains(o)),
        "time {o:?}"
    );
}

/// V7.37: the heated wall with the heater switched on at t = 0 (the
/// gallery case on a coarser mesh). The heat balance closes in every step
/// (heat.csv), the heat goes first into the wall and at the end out with
/// the flow, and the end state is the steady solution.
#[test]
#[cfg_attr(debug_assertions, ignore = "release builds only (CI runs it)")]
fn heated_wall_switched_on_reaches_the_steady_state() {
    let coarse = |t: &str| t.replace("cells = [80, 30, 1]", "cells = [40, 15, 1]");
    let steady = run_case(
        &parse(&coarse(&text("cht_heated_wall.toml"))),
        "wall_steady",
    );
    assert!(steady.failures.is_empty(), "{:?}", steady.failures);
    let case = parse(&coarse(&text("cht_heated_wall_transient.toml")));
    let dir = out("wall_transient");
    let r = run(&case, "wall_transient", &dir, |_| {}).unwrap();
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let (a, b) = (temperature(&steady), temperature(&r));
    let dt = a.max - a.min;
    for n in ["ymin", "interface-wall", "outlet"] {
        let (x, y) = (boundary(a, n).mean, boundary(b, n).mean);
        assert!((x - y).abs() < 1e-3 * dt, "{n}: {y} against the steady {x}");
    }
    let h = r.solves[0].transient.as_ref().unwrap();
    assert!(h.heat_relative_max.unwrap() < 1e-9, "{h:?}");
    // the heat history: every step balanced; the wall stores most of the
    // heat at first, almost nothing at the end
    let text = std::fs::read_to_string(dir.join("heat.csv")).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let (cw, cc, ci, cr) = (
        col("wall_storage"),
        col("conduction"),
        col("interface-wall_heat"),
        col("relative"),
    );
    let rows: Vec<Vec<f64>> = lines
        .map(|l| l.split(',').map(|v| v.parse().unwrap()).collect())
        .collect();
    assert_eq!(rows.len(), 600);
    assert!(rows.iter().all(|r| r[cr] < 1e-9));
    let (first, last) = (&rows[0], rows.last().unwrap());
    assert!(first[cw] > 0.95 * first[cc], "{first:?}");
    assert!(last[cw].abs() < 1e-4 * last[cc], "{last:?}");
    // all of the heater's 160 W through the interface
    assert!((last[ci] / 160.0 - 1.0).abs() < 1e-4, "{last:?}");
    // the probes: the temperature falls from the heated face to the
    // channel; no velocity in the solid
    let text = std::fs::read_to_string(dir.join("probes.csv")).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |n: &str| head.iter().position(|h| *h == n).unwrap();
    let last: Vec<f64> = lines
        .next_back()
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect();
    let t = |n: &str| last[col(&format!("{n}_T"))];
    assert!(t("heated_face") > t("interface") && t("interface") > t("channel"));
    assert_eq!(last[col("heated_face_u")], 0.0);
    let _ = std::fs::remove_dir_all(&dir);
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

/// Stage 7d2, V7.4 with resistance: the composite wall with a contact
/// resistance at the lower interface and a thin wall (t/k) at the upper
/// one; the series heat flux with both and the temperature jumps q R″.
#[test]
fn composite_wall_with_contact_resistance_and_thin_wall() {
    let (rc, t, k) = (2.0e-3, 1.0e-3, 0.5);
    let case = text("cht_composite_wall.toml").replace(
        "[energy]",
        &format!(
            "[interfaces.a-air]\nregions = [\"a\", \"fluid\"]\ncontact_resistance = {rc}\n\n\
             [interfaces.b-air]\nregions = [\"fluid\", \"b\"]\n\
             thin_wall = {{ thickness = {t}, conductivity = {k} }}\n\n[energy]"
        ),
    );
    let r = run_case(&parse(&case), "wall_contact");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = temperature(&r);
    let q = 100.0 / (0.4 / 200.0 + 0.4 / 2.0 + 0.4 / 20000.0 + rc + t / k);
    for (name, sign) in [
        ("ymin", -1.0),
        ("ymax", 1.0),
        ("interface-a", -1.0),
        ("interface-b", 1.0),
    ] {
        let f = boundary(s, name).flux;
        assert!(
            (f - sign * q).abs() < 1e-9 * q,
            "{name}: {f} against {}",
            sign * q
        );
    }
    // the fluid-side face temperatures: above the solid's by q R″ at a,
    // below it by q t/k at b
    let ta = 300.0 + q * (0.4 / 200.0 + rc);
    let tb = 400.0 - q * (0.4 / 20000.0 + t / k);
    let (ma, mb) = (
        boundary(s, "interface-a").mean,
        boundary(s, "interface-b").mean,
    );
    assert!((ma - ta).abs() < 1e-9, "{ma} against {ta}");
    assert!((mb - tb).abs() < 1e-9, "{mb} against {tb}");
}

/// Stage 7d2: the patch test of `cht_orthotropic_patch.toml` (an
/// orthotropic laminate with axes turned by 30°, a fluid at rest, a block
/// behind a contact resistance, on a mesh with nodes moved by 25 %):
/// piecewise-linear T exact to round-off, the same heat flux 200 W/m²
/// through both outer walls and both interfaces, and a closed balance.
/// The same with the axes of the region in place of the material's.
#[test]
fn orthotropic_patch_is_exact_on_a_skewed_mesh() {
    let base = text("cht_orthotropic_patch.toml");
    let region_axes = base
        .replace(
            ", axes = [[0.8660254037844387, 0.5, 0.0], [-0.5, 0.8660254037844387, 0.0]] }",
            " }",
        )
        .replace(
            "material = \"fibre-composite\"\n",
            "material = \"fibre-composite\"\naxes = [[0.8660254037844387, 0.5, 0.0], [-0.5, 0.8660254037844387, 0.0]]\n",
        );
    for (name, text) in [("patch", &base), ("patch_region_axes", &region_axes)] {
        let r = run_case(&parse(text), name);
        assert!(r.failures.is_empty(), "{name}: {:?}", r.failures);
        let s = temperature(&r);
        assert!(s.linf.unwrap() < 1e-9, "{name}: max error {:?} K", s.linf);
        for (b, sign) in [
            ("xmin", -1.0),
            ("xmax", 1.0),
            ("interface-laminate", -1.0),
            ("interface-block", 1.0),
        ] {
            let f = boundary(s, b).flux;
            assert!(
                (f - sign * 200.0).abs() < 1e-8 * 200.0,
                "{name} {b}: {f} W/m²"
            );
        }
        // the fluid side of the contact: 2 K below the block
        let t = boundary(s, "interface-block").mean;
        assert!((t - (320.0 + 50.0 * 0.5)).abs() < 1e-9, "{name}: {t}");
        let e = s.energy.as_ref().unwrap();
        assert!(e.relative < 1e-10, "{name}: {e:?}");
    }
}

/// Stage 7d2: `[interfaces]` and `[[energy.sources]]` are checked.
#[test]
fn interfaces_and_region_sources_are_checked() {
    let base = text("cht_composite_wall.toml");
    let with = |extra: &str| base.replace("[energy]", &format!("{extra}\n[energy]"));
    for (extra, code) in [
        // an unknown region
        (
            "[interfaces.x]\nregions = [\"a\", \"c\"]\ncontact_resistance = 1e-3\n",
            "CHT004",
        ),
        // two solids that do not touch
        (
            "[interfaces.x]\nregions = [\"a\", \"b\"]\ncontact_resistance = 1e-3\n",
            "CHT004",
        ),
        // neither a resistance nor a thin wall
        ("[interfaces.x]\nregions = [\"a\", \"fluid\"]\n", "CHT004"),
        // a negative resistance
        (
            "[interfaces.x]\nregions = [\"a\", \"fluid\"]\ncontact_resistance = -1.0\n",
            "CHT004",
        ),
        // the same region twice
        (
            "[interfaces.x]\nregions = [\"a\", \"a\"]\ncontact_resistance = 1e-3\n",
            "CHT004",
        ),
    ] {
        let (c, msg) = refused(&with(extra));
        assert_eq!(c, code, "{extra}: {msg}");
    }
    let src = |s: &str| base.replace("initial = 350.0\n", &format!("initial = 350.0\n{s}"));
    for s in [
        "[[energy.sources]]\nregion = \"c\"\npower = 1.0\n",
        "[[energy.sources]]\nregion = \"a\"\n",
        "[[energy.sources]]\nregion = \"a\"\npower = 1.0\npower_density = 2.0\n",
    ] {
        let (c, msg) = refused(&src(s));
        assert_eq!(c, "ENE008", "{s}: {msg}");
    }
}

/// Stage 7d2: a region source of 5 W in the lower solid of the composite
/// wall: the balance closes with it, and the heat leaving through the two
/// walls is 5 W more than without it.
#[test]
fn region_source_enters_the_balance() {
    let base = text("cht_composite_wall.toml");
    let r0 = run_case(&parse(&base), "src0");
    let heated = base.replace(
        "initial = 350.0\n",
        "initial = 350.0\n\n[[energy.sources]]\nregion = \"a\"\npower = 5.0\n",
    );
    let r = run_case(&parse(&heated), "src");
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let s = temperature(&r);
    let e = s.energy.as_ref().unwrap();
    assert!((e.source - 5.0).abs() < 1e-12, "{e:?}");
    assert!(e.relative < 1e-9, "{e:?}");
    let out = |s: &foehncraft_flow::scalar::ScalarSummary| {
        boundary(s, "ymin").heat_flow.unwrap() + boundary(s, "ymax").heat_flow.unwrap()
    };
    let d = out(temperature(&r0)) - out(s);
    assert!((d - 5.0).abs() < 1e-8, "{d}");
}

/// Stage 7d2: the two bodies of V7.6 with a contact resistance between
/// them (R″_c = 10⁻³ m²K/W, the conductance 1000 W/(m² K) of the gap
/// against √(kρc/t) ≈ 140 at 500 s): the heat balance closes in every
/// step, what the solid loses the fluid stores, and the two sides of the
/// contact differ by q R″_c.
#[test]
fn transient_contact_resistance_closes_the_balance() {
    let base = text("cht_contact_transient.toml")
        .replace("end_time = 500.0", "end_time = 50.0")
        .replace("scalar_max_error = 0.01 ", "scalar_max_error = 100.0 ")
        .replace(
            "\n[energy]\n",
            "\n[interfaces.gap]\nregions = [\"bar\", \"fluid\"]\ncontact_resistance = 1e-3\n\n[energy]\n",
        );
    let r = run_case(&parse(&base), "contact_gap");
    let s = temperature(&r);
    let e = s.energy.as_ref().unwrap();
    let st: Vec<f64> = e.storage_regions.iter().map(|r| r.storage).collect();
    assert!((st[0] + st[1]).abs() < 1e-9 * st[0].abs(), "{e:?}");
    let h = r.solves[0].transient.as_ref().unwrap();
    assert!(h.heat_relative_max.unwrap() < 1e-9, "{h:?}");
    // the fluid side of the gap is below the contact temperature without it
    let ti = contact_temperature([10.0, 1.0], [1e6, 1e6], [400.0, 300.0]);
    let t = boundary(s, "interface-bar").mean;
    assert!(t < ti - 1.0, "{t} against {ti}");
}
