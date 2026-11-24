//! IAPWS-IF97 against the verification tables of the release (IAPWS
//! R7-97(2012), Tables 5, 15, 33, 35, 36, 42 and the B23 values), the
//! viscosity against IAPWS R12-08 Table 4 (μ̄₂ = 1), the conductivity
//! against IAPWS R15-11 Table 4 at the points where the critical
//! enhancement is negligible, and the surface tension against IAPWS
//! R1-76(2014). The tables give 9 significant digits; the band is half a
//! unit of the last digit (5e-9 relative), the rounding of the table.

use foehncraft_kernels::steam::{self, Phase};

fn close(what: &str, got: f64, want: f64, rel: f64) {
    let e = ((got - want) / want).abs();
    assert!(
        e <= rel,
        "{what}: got {got:.12e}, want {want:.12e}, rel. error {e:.2e} > {rel:.1e}"
    );
}

const TAB: f64 = 5e-9;

/// (T K, p MPa, v m³/kg, h kJ/kg, u kJ/kg, s kJ/(kg K), c_p kJ/(kg K), w m/s)
type Row = (f64, f64, f64, f64, f64, f64, f64, f64);

fn check_rows(name: &str, rows: &[Row], f: impl Fn(f64, f64) -> steam::State) {
    for &(t, p, v, h, u, s, cp, w) in rows {
        let st = f(p * 1e6, t);
        let at = format!("{name} T = {t} K, p = {p} MPa");
        close(&format!("{at} v"), st.v(), v, TAB);
        close(&format!("{at} h"), st.h, h * 1e3, TAB);
        close(&format!("{at} u"), st.u, u * 1e3, TAB);
        close(&format!("{at} s"), st.s, s * 1e3, TAB);
        close(&format!("{at} cp"), st.cp, cp * 1e3, TAB);
        close(&format!("{at} w"), st.w, w, TAB);
    }
}

#[test]
fn region1_table5() {
    check_rows(
        "region 1",
        &[
            (
                300.0,
                3.0,
                0.100215168e-2,
                0.115331273e3,
                0.112324818e3,
                0.392294792,
                0.417301218e1,
                0.150773921e4,
            ),
            (
                300.0,
                80.0,
                0.971180894e-3,
                0.184142828e3,
                0.106448356e3,
                0.368563852,
                0.401008987e1,
                0.163469054e4,
            ),
            (
                500.0,
                3.0,
                0.120241800e-2,
                0.975542239e3,
                0.971934985e3,
                0.258041912e1,
                0.465580682e1,
                0.124071337e4,
            ),
        ],
        steam::region1,
    );
}

#[test]
fn region2_table15() {
    check_rows(
        "region 2",
        &[
            (
                300.0,
                0.0035,
                0.394913866e2,
                0.254991145e4,
                0.241169160e4,
                0.852238967e1,
                0.191300162e1,
                0.427920172e3,
            ),
            (
                700.0,
                0.0035,
                0.923015898e2,
                0.333568375e4,
                0.301262819e4,
                0.101749996e2,
                0.208141274e1,
                0.644289068e3,
            ),
            (
                700.0,
                30.0,
                0.542946619e-2,
                0.263149474e4,
                0.246861076e4,
                0.517540298e1,
                0.103505092e2,
                0.480386523e3,
            ),
        ],
        steam::region2,
    );
}

#[test]
fn region5_table42() {
    check_rows(
        "region 5",
        &[
            (
                1500.0,
                0.5,
                0.138455090e1,
                0.521976855e4,
                0.452749310e4,
                0.965408875e1,
                0.261609445e1,
                0.917068690e3,
            ),
            (
                1500.0,
                30.0,
                0.230761299e-1,
                0.516723514e4,
                0.447495124e4,
                0.772970133e1,
                0.272724317e1,
                0.928548002e3,
            ),
            (
                2000.0,
                30.0,
                0.311385219e-1,
                0.657122604e4,
                0.563707038e4,
                0.853640523e1,
                0.288569882e1,
                0.106736948e4,
            ),
        ],
        steam::region5,
    );
}

#[test]
fn region3_table33() {
    // (T K, ρ kg/m³, p MPa, h, u, s, c_p, w)
    for (t, rho, p, h, u, s, cp, w) in [
        (
            650.0,
            500.0,
            0.255837018e2,
            0.186343019e4,
            0.181226279e4,
            0.405427273e1,
            0.138935717e2,
            0.502005554e3,
        ),
        (
            650.0,
            200.0,
            0.222930643e2,
            0.237512401e4,
            0.226365868e4,
            0.485438792e1,
            0.446579342e2,
            0.383444594e3,
        ),
        (
            750.0,
            500.0,
            0.783095639e2,
            0.225868845e4,
            0.210206932e4,
            0.446971906e1,
            0.634165359e1,
            0.760696041e3,
        ),
    ] {
        let st = steam::region3(rho, t);
        let at = format!("region 3 T = {t} K, rho = {rho}");
        close(&format!("{at} p"), st.p, p * 1e6, TAB);
        close(&format!("{at} h"), st.h, h * 1e3, TAB);
        close(&format!("{at} u"), st.u, u * 1e3, TAB);
        close(&format!("{at} s"), st.s, s * 1e3, TAB);
        close(&format!("{at} cp"), st.cp, cp * 1e3, TAB);
        close(&format!("{at} w"), st.w, w, TAB);
        // the density back from (p, T) through the region selection
        let back = steam::props_pt(st.p, t).expect("region 3 state");
        assert_eq!(back.region, 3, "{at}");
        close(&format!("{at} rho(p, T)"), back.rho, rho, 1e-12);
    }
}

#[test]
fn region4_tables35_36_and_b23() {
    for (t, p) in [
        (300.0, 0.353658941e-2),
        (500.0, 0.263889776e1),
        (600.0, 0.123443146e2),
    ] {
        close(&format!("p_s({t})"), steam::p_sat(t), p * 1e6, TAB);
    }
    for (p, t) in [
        (0.1, 0.372755919e3),
        (1.0, 0.453035632e3),
        (10.0, 0.584149488e3),
    ] {
        close(&format!("T_s({p})"), steam::t_sat(p * 1e6), t, TAB);
    }
    close(
        "p_B23(623.15)",
        steam::p_b23(623.15),
        0.165291643e2 * 1e6,
        TAB,
    );
    close(
        "T_B23(16.5291643)",
        steam::t_b23(0.165291643e2 * 1e6),
        0.623150000e3,
        TAB,
    );
    close("p_s(273.15)", steam::p_sat(273.15), steam::P_SAT_MIN, 1e-9);
    assert!(steam::p_sat(273.15) >= steam::P_SAT_MIN);
}

#[test]
fn regions_are_selected_and_continuous() {
    // region numbers on a grid, and the jump of the properties across the
    // region boundaries (IF97 consistency requirements: below 0.2 % in c_p
    // and 1e-5 relative... as stated in the release, here bands from it)
    assert_eq!(steam::region_pt(3e6, 300.0), 1);
    assert_eq!(steam::region_pt(0.0035e6, 300.0), 2);
    assert_eq!(steam::region_pt(25e6, 650.0), 3);
    assert_eq!(steam::region_pt(30e6, 1500.0), 5);
    assert_eq!(steam::region_pt(60e6, 1500.0), 0);
    assert_eq!(steam::region_pt(101e6, 500.0), 0);
    assert_eq!(steam::region_pt(1e6, 270.0), 0);
    // across T = 623.15 K (regions 1 and 3) and B23 (regions 2 and 3)
    for p in [20e6, 40e6, 80e6] {
        let a = steam::region1(p, 623.15);
        let b = steam::region3(steam::rho3(p, 623.15, true), 623.15);
        close(&format!("rho 1|3 at {p}"), b.rho, a.rho, 1e-4);
        close(&format!("h 1|3 at {p}"), b.h, a.h, 1e-4);
        let t = steam::t_b23(p);
        let c = steam::region2(p, t);
        let d = steam::region3(steam::rho3(p, t, false), t);
        close(&format!("rho 2|3 at {p}"), d.rho, c.rho, 5e-4);
        close(&format!("h 2|3 at {p}"), d.h, c.h, 1e-4);
    }
}

#[test]
fn saturation_states_and_phase_equilibrium() {
    // the Gibbs energies of the two saturated phases agree to the
    // consistency of the region 4 equation with the basic equations of
    // regions 1 to 3 (IF97: |Δg|/(R T) up to about 7e-5, 11 J/kg at 0.1 MPa)
    for p in [0.01e6, 0.1e6, 1e6, 5e6, 10e6, 16e6, 18e6, 20e6, 22e6] {
        let s = steam::saturation_p(p).expect("saturation");
        let gl = s.liquid.h - s.t * s.liquid.s;
        let gv = s.vapour.h - s.t * s.vapour.s;
        let dg = (gl - gv) / (steam::R * s.t);
        assert!(dg.abs() < 1e-4, "g' = {gl}, g'' = {gv} at {p}: {dg:.2e}");
        assert!(s.vapour.rho < s.liquid.rho);
        assert!(s.h_fg() > 0.0);
    }
    // latent heat at 0.101325 MPa: 2256.5 kJ/kg (IAPWS-95: 2256.47)
    let s = steam::saturation_p(101_325.0).unwrap();
    close("h_fg at 1 atm", s.h_fg(), 2256.47e3, 2e-4);
    close("T_s at 1 atm", s.t, 373.124, 1e-5);
}

#[test]
fn ph_ps_px_round_trips() {
    for (p, t) in [
        (0.1e6, 300.0),
        (0.1e6, 500.0),
        (3e6, 500.0),
        (10e6, 700.0),
        (25e6, 650.0),
        (25e6, 700.0),
        (30e6, 1500.0),
        (80e6, 300.0),
        (12e6, 1100.0),
    ] {
        let s = steam::state_pt(p, t).unwrap();
        let a = steam::state_ph(p, s.h).unwrap();
        close(&format!("T(p, h) at {p}, {t}"), a.t, t, 1e-10);
        let b = steam::state_ps(p, s.s).unwrap();
        close(&format!("T(p, s) at {p}, {t}"), b.t, t, 1e-10);
    }
    // two-phase: (p, x) → h → (p, h) gives the same quality
    for p in [0.1e6, 3e6, 18e6] {
        let m = steam::state_px(p, 0.3).unwrap();
        assert_eq!(m.phase, Phase::TwoPhase);
        let back = steam::state_ph(p, m.h).unwrap();
        close("x(p, h)", back.x, 0.3, 1e-12);
        close("rho_m", back.rho, m.rho, 1e-12);
        assert!(back.alpha > 0.3 && back.alpha < 1.0);
    }
    // subcooled and superheated qualities
    let sub = steam::state_pt(7e6, 500.0).unwrap();
    assert!(sub.x < 0.0 && sub.phase == Phase::Liquid);
    let sup = steam::state_pt(7e6, 700.0).unwrap();
    assert!(sup.x > 1.0 && sup.phase == Phase::Vapour);
    let sc = steam::state_pt(25e6, 700.0).unwrap();
    assert!(sc.x.is_nan() && sc.phase == Phase::Supercritical);
}

#[test]
fn viscosity_r12_08_table4() {
    for (t, rho, mu) in [
        (298.15, 998.0, 889.735100),
        (298.15, 1200.0, 1437.649467),
        (373.15, 1000.0, 307.883622),
        (433.15, 1.0, 14.538324),
        (433.15, 1000.0, 217.685358),
        (873.15, 1.0, 32.619287),
        (873.15, 100.0, 35.802262),
        (873.15, 600.0, 77.430195),
        (1173.15, 1.0, 44.217245),
        (1173.15, 100.0, 47.640433),
        (1173.15, 400.0, 64.154608),
    ] {
        // the table gives 6 decimals in μPa s: half a unit of the last one
        let got = steam::viscosity(rho, t) * 1e6;
        assert!(
            (got - mu).abs() <= 5e-7 + 1e-12 * mu,
            "mu({t}, {rho}): got {got:.9}, want {mu}"
        );
    }
}

#[test]
fn conductivity_r15_11_table4() {
    // the dilute gas (ρ → 0) and two liquid states where λ̄₂ vanishes
    for (t, p, lam) in [(298.15, 1e-3, 18.4341883), (873.15, 1e-3, 79.1034659)] {
        let s = steam::state_pt(p, t).unwrap();
        close(
            &format!("lambda({t}, 0)"),
            steam::conductivity(&s),
            lam * 1e-3,
            1e-7,
        );
    }
    // (298.15 K, 1200 kg/m³ of the table lies far above 100 MPa, outside
    // IF97, and is left out)
    for (t, rho, lam) in [(298.15, 998.0, 607.712868)] {
        // the pressure of (ρ, T) in region 1 by bisection
        let (mut lo, mut hi) = (1e3, 1e9);
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            if steam::region1(m, t).rho < rho {
                lo = m
            } else {
                hi = m
            }
        }
        let s = steam::props_pt(0.5 * (lo + hi), t).unwrap();
        close(
            &format!("lambda({t}, {rho})"),
            steam::conductivity(&s),
            lam * 1e-3,
            1e-6,
        );
    }
}

#[test]
fn surface_tension_r1_76() {
    // Table 1 of the release: 0.01 °C 75.65, 100 °C 58.91, 300 °C 14.37 mN/m
    close(
        "sigma(273.16)",
        steam::surface_tension(273.16),
        75.65e-3,
        1e-3,
    );
    close(
        "sigma(373.15)",
        steam::surface_tension(373.15),
        58.91e-3,
        1e-3,
    );
    close(
        "sigma(573.15)",
        steam::surface_tension(573.15),
        14.37e-3,
        1e-3,
    );
    assert_eq!(steam::surface_tension(steam::T_CRIT), 0.0);
}

#[test]
fn outside_the_range_is_none() {
    assert!(steam::state_pt(120e6, 500.0).is_none());
    assert!(steam::state_pt(1e6, 200.0).is_none());
    assert!(steam::state_ph(1e6, -1e7).is_none());
    assert!(steam::state_ph(1e6, 1e8).is_none());
    assert!(steam::state_px(30e6, 0.5).is_none());
    assert!(steam::p_sat(700.0).is_nan());
    let o = steam::state_outputs(steam::Input::PT, 1e6, 100.0);
    assert!(o.iter().all(|v| v.is_nan()));
}
