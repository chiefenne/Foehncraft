//! Evaporating droplets (SP4, SP5, D247 to D253) in a gas of one cell that
//! the test solves itself:
//!
//! - many parcels in one cell near saturation: the vapour slope (the chord
//!   to the saturation at the droplet surface, D250) keeps the gas from
//!   going above saturation or below zero vapour, the closing step keeps
//!   the water to round-off, and the results are the same bits on one and
//!   on four threads with the parcels stored in another order (D247, D249);
//! - a single droplet in hot dry air follows the d² law at the wet-bulb
//!   temperature (V28.7);
//! - a closed adiabatic box of hot dry air with more water than it can
//!   take ends at the adiabatic saturation state of the enthalpy balance,
//!   with water and energy kept to round-off (V28.5).

use foehncraft_kernels::evaporation::{surface_fraction, water};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, CouplingMode, Delta, Env, EvaporationModel, Gas, GasState,
    Injector, Parcels, Spray, SpraySection, Vapour,
};

const P0: f64 = 101_325.0;
const M_AIR: f64 = 0.028_965;
const CP_AIR: f64 = 1005.0;

/// The gas of the one cell: at rest, with the temperature and vapour mass
/// fraction the test holds.
struct Cell {
    t: f64,
    y: f64,
    rho: f64,
}

impl Gas for Cell {
    fn at(&self, _: usize, _: [f64; 3]) -> GasState {
        GasState {
            u: [0.0; 3],
            rho: self.rho,
            mu: 1.85e-5,
            t: self.t,
            cp: CP_AIR,
            k: 0.026,
            tke: 0.0,
            eps: 0.0,
            y_v: self.y,
            p0: P0,
            m_gas: M_AIR,
            diff: 2.5e-5,
        }
    }
}

/// Saturation mass fraction of water vapour in air at `t`.
fn ys(t: f64) -> f64 {
    surface_fraction(water::saturation_pressure(t), P0, water::MOLAR_MASS, M_AIR)
}

fn one_cell(edge: f64) -> (foehncraft_mesh::Mesh, Geometry) {
    let mut s = BoxSpec::new(Kind::Hex, 3, [1, 1, 1]);
    s.lower = [0.0; 3];
    s.upper = [edge; 3];
    let m = s.generate().unwrap();
    let g = Geometry::compute(&m);
    (m, g)
}

#[allow(clippy::too_many_arguments)]
fn spray(
    env: &Env,
    edge: f64,
    mass_flow: f64,
    pps: f64,
    d: f64,
    temp: f64,
    end: Option<f64>,
    two_way: bool,
) -> Spray {
    let mut cfg = SpraySection {
        injectors: vec![Injector {
            name: "n".into(),
            position: [0.5 * edge; 3],
            direction: [1.0, 0.0, 0.0],
            speed: 1e-4,
            half_angle: 0.0,
            mass_flow,
            diameter: d,
            size: None,
            parcels_per_second: pps,
            start: 0.0,
            end,
            temperature: Some(temp),
            ..Default::default()
        }],
        ..SpraySection::default()
    };
    if two_way {
        cfg.coupling.mode = CouplingMode::TwoWay;
    }
    cfg.liquid.vapour = Some("H2O".into());
    cfg.evaporation.model = EvaporationModel::RanzMarshall;
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Rebound,
            e_normal: 0.5,
            e_tangential: 0.5,
        };
        env.mesh.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, env, rules, [0.0; 3]).unwrap();
    sp.set_vapour(Some(Vapour::water()));
    sp
}

fn reversed(p: &Parcels) -> Parcels {
    let mut q = p.clone();
    q.id.reverse();
    q.x.reverse();
    q.v.reverse();
    q.d.reverse();
    q.count.reverse();
    q.born.reverse();
    q.temp.reverse();
    q.eddy.reverse();
    q.cell.reverse();
    q
}

/// The run of the owner's question: many parcels evaporating into one cell
/// whose gas starts at 95 % of saturation. Returns the vapour mass fraction
/// after each step and the parcels at the end (sorted by identity).
fn crowded_cell(threads: usize, shuffle: bool) -> (Vec<f64>, Vec<(u64, f64, f64)>, f64) {
    let edge = 0.01;
    let (m, geo) = one_cell(edge);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let t_gas = 300.0;
    let rho = P0 * M_AIR / (8.314_462_618 * t_gas);
    let rho_v = rho * geo.cell_volume[0];
    // 20 kg/s per m³ of liquid for 0.2 s: far more than the gas can take
    // (about 1.3e-9 kg of vapour between 95 % and saturation), in 4000
    // parcels, a little warmer than the gas
    let mut sp = spray(&env, edge, 2e-5, 20_000.0, 20e-6, 305.0, None, true);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let mut y = 0.95 * ys(t_gas);
    let mut ys_hist = Vec::new();
    let water0 = rho_v * y;
    let dt = 1e-3;
    for n in 0..200 {
        if shuffle {
            sp.parcels = reversed(&sp.parcels);
        }
        let gas = Cell { t: t_gas, y, rho };
        let pass = pool.install(|| sp.pass(&env, &gas, n as f64 * dt, dt).unwrap());
        let s = &pass.sources;
        // the gas of the cell: ρ V (Y₁ − Y) = E − β (Y₁ − Y*) with Y* = Y
        let (e, b) = match s.cells.first() {
            Some(_) => (s.vapour[0], s.vapour_slope[0]),
            None => (0.0, 0.0),
        };
        assert!(b >= 0.0 && b.is_finite(), "step {n}: slope {b}");
        let y1 = y + e / (rho_v + b);
        // the highest saturation any droplet surface or the gas can hold
        let t_max = sp
            .parcels
            .temp
            .iter()
            .copied()
            .fold(305.0f64.max(t_gas), f64::max);
        assert!(y1 >= 0.0, "step {n}: negative vapour {y1}");
        assert!(
            y1 <= ys(t_max) * (1.0 + 1e-12),
            "step {n}: Y = {y1} above the saturation {} of {t_max} K",
            ys(t_max)
        );
        let dy = y1 - y;
        let books = pool.install(|| {
            sp.commit(
                &env,
                pass,
                Some(&|_| Delta {
                    u: [0.0; 3],
                    t: 0.0,
                    y: dy,
                }),
            )
        });
        // what the gas received is what the parcels gave
        let to_gas = rho_v * dy;
        // (Y₁ − Y is a difference of two numbers near Y: its round-off is
        // that of Y)
        assert!(
            (to_gas + books.mass[2]).abs() <= 1e-12 * to_gas.abs() + 1e-14 * rho_v * y,
            "step {n}: gas {to_gas:e}, parcels {:e}",
            books.mass[2]
        );
        y = y1;
        ys_hist.push(y);
    }
    // the water of gas and liquid: what was there and what was injected
    let water = rho_v * y + sp.mass_in_flow_by_id();
    let injected = sp.tally.injected_mass;
    let err = (water - water0 - injected).abs() / (injected + water0);
    let mut parcels: Vec<(u64, f64, f64)> = (0..sp.parcels.len())
        .map(|i| (sp.parcels.id[i], sp.parcels.d[i], sp.parcels.temp[i]))
        .collect();
    parcels.sort_by_key(|p| p.0);
    for p in &parcels {
        assert!(p.1 >= 0.0 && p.1.is_finite() && p.2.is_finite(), "{p:?}");
    }
    (ys_hist, parcels, err)
}

#[test]
fn many_parcels_in_one_cell_near_saturation() {
    let (ya, pa, err) = crowded_cell(1, false);
    let (yb, pb, _) = crowded_cell(4, true);
    assert!(pa.len() > 3000, "{} parcels", pa.len());
    // bitwise the same on one and four threads, whatever the storage order
    assert_eq!(ya, yb);
    assert_eq!(pa, pb);
    // the water is kept to round-off
    assert!(err < 1e-12, "water book {err:e}");
    // the gas approaches saturation from below and stays there
    let sat = ys(300.0);
    let last = *ya.last().unwrap();
    assert!(last > 0.95 * sat, "{last} against {sat}");
    eprintln!(
        "crowded cell: {} parcels, Y from {:.6e} to {last:.6e} (saturation of the gas \
         {sat:.6e}, of the 305 K droplets {:.6e}), water book {err:.1e}",
        pa.len(),
        0.95 * sat,
        ys(305.0)
    );
}

/// The wet-bulb temperature of a droplet of diameter `d` in `gas`: the
/// temperature whose equilibrium temperature of the exact integrator is
/// itself.
fn wet_bulb(sp: &Spray, gas: &GasState, d: f64) -> f64 {
    let (mut lo, mut hi) = (250.0, gas.t);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let e = sp.evaporation(gas, d, mid, 0.0, 1e-9).unwrap();
        if e.t_eq > mid {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn a_single_droplet_follows_the_d2_law_at_the_wet_bulb() {
    let edge = 0.1;
    let (m, geo) = one_cell(edge);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let dt = 1e-3;
    // one droplet of 100 µm at 300 K in dry air at 500 K (one way)
    let mut sp = spray(&env, edge, 1e-9, 1.0 / dt, 100e-6, 300.0, Some(dt), false);
    let gas = Cell {
        t: 500.0,
        y: 0.0,
        rho: P0 * M_AIR / (8.314_462_618 * 500.0),
    };
    let mut hist = Vec::new();
    for n in 0..20_000 {
        sp.advance(&env, &gas, n as f64 * dt, dt).unwrap();
        if sp.parcels.is_empty() {
            break;
        }
        hist.push((
            (n + 1) as f64 * dt,
            sp.parcels.d[0] * sp.parcels.d[0],
            sp.parcels.temp[0],
        ));
    }
    assert_eq!(sp.tally.evaporated_parcels, 1, "the droplet evaporated");
    assert!(hist.len() > 100, "{} steps", hist.len());
    // after the heat-up d² falls linearly: the local slopes of the middle
    // half agree to 1 % (the wet bulb moves a little as the droplet
    // shrinks, through the film's Nusselt number only at rest: not at all)
    let n = hist.len();
    let mid = &hist[n / 4..3 * n / 4];
    let k: Vec<f64> = mid
        .windows(2)
        .map(|w| -(w[1].1 - w[0].1) / (w[1].0 - w[0].0))
        .collect();
    let kmax = k.iter().copied().fold(f64::MIN, f64::max);
    let kmin = k.iter().copied().fold(f64::MAX, f64::min);
    assert!((kmax - kmin) / kmax < 0.01, "K from {kmin:e} to {kmax:e}");
    // at the wet-bulb temperature
    let g = gas.at(0, [0.0; 3]);
    let (_, d2, tp) = hist[n / 2];
    let twb = wet_bulb(&sp, &g, d2.sqrt());
    assert!((tp - twb).abs() < 0.05, "droplet {tp} K, wet bulb {twb} K");
    // the life: the heat-up plus d0²/K
    let life = hist.last().unwrap().0;
    let d2_law = (100e-6f64).powi(2) / kmax;
    assert!(
        life > d2_law && life < 1.2 * d2_law,
        "life {life} s, d0²/K {d2_law} s"
    );
    eprintln!(
        "single droplet: K = {kmax:.4e} m²/s, wet bulb {twb:.3} K, droplet {tp:.3} K, life \
         {life:.3} s (d0²/K {d2_law:.3} s)"
    );
}

#[test]
fn a_closed_box_ends_at_the_adiabatic_saturation() {
    let edge = 0.01;
    let (m, geo) = one_cell(edge);
    let loc = Locator::new(&m, &geo);
    let env = Env {
        mesh: &m,
        geo: &geo,
        loc: &loc,
    };
    let vol = geo.cell_volume[0];
    let dt = 2e-3;
    let t0 = 400.0;
    let rho0 = P0 * M_AIR / (8.314_462_618 * t0);
    let m_air = rho0 * vol;
    // 1e-7 kg of water at 300 K in 100 parcels of 20 µm droplets, released
    // in the first step: more than the air can take
    let mut sp = spray(
        &env,
        edge,
        1e-7 / dt,
        100.0 / dt,
        20e-6,
        300.0,
        Some(dt),
        true,
    );
    // the vapour enthalpy of the parcels' reference, c_l T + L(T)
    let liq = sp.config().liquid.clone();
    let c_l = liq.specific_heat;
    let wv = Vapour::water();
    let hv = |t: f64| {
        c_l * t + liq.latent_heat + wv.h(t)
            - wv.h(liq.latent_heat_temperature)
            - c_l * (t - liq.latent_heat_temperature)
    };
    let cpv = |t: f64| wv.cp(t);
    // the gas: air and vapour, its enthalpy in the parcels' reference
    let mut m_v = 0.0f64;
    let mut e_g = m_air * CP_AIR * t0;
    let temp = |e: f64, mv: f64| {
        // E = m_a c_p T + m_v h_v(T), by Newton
        let mut t = 350.0;
        for _ in 0..50 {
            let f = m_air * CP_AIR * t + mv * hv(t) - e;
            let df = m_air * CP_AIR + mv * cpv(t);
            t -= f / df;
        }
        t
    };
    let (mut t, mut y) = (t0, 0.0);
    let mut first = None;
    for n in 0..3000 {
        let gas = Cell {
            t,
            y,
            rho: rho0 * t0 / t,
        };
        let pass = sp.pass(&env, &gas, n as f64 * dt, dt).unwrap();
        let s = pass.sources.clone();
        let (ev, bv, h, g, vh) = match s.cells.first() {
            Some(_) => (
                s.vapour[0],
                s.vapour_slope[0],
                s.heat[0],
                s.heat_slope[0],
                s.vapour_slope_h[0],
            ),
            None => (0.0, 0.0, 0.0, 0.0, 0.0),
        };
        // the gas at the end of the step with the value-plus-slope
        // sources, by fixed-point iteration
        let (mut y1, mut t1) = (y, t);
        for _ in 0..100 {
            let vap = ev - bv * (y1 - y);
            let heat = h - g * (t1 - t) - vh * (y1 - y);
            let mv1 = m_v + vap;
            y1 = mv1 / (m_air + mv1);
            t1 = temp(e_g + heat, mv1);
        }
        let (dy, dtemp) = (y1 - y, t1 - t);
        sp.commit(
            &env,
            pass,
            Some(&|_| Delta {
                u: [0.0; 3],
                t: dtemp,
                y: dy,
            }),
        );
        let vap = ev - bv * dy;
        let heat = h - g * dtemp - vh * dy;
        m_v += vap;
        e_g += heat;
        t = temp(e_g, m_v);
        y = m_v / (m_air + m_v);
        let liquid = sp.mass_in_flow_by_id();
        let water = m_v + liquid;
        let energy = e_g + sp.enthalpy();
        let (w0, e0) = *first.get_or_insert((water, energy));
        assert!(
            (water - w0).abs() < 1e-12 * w0,
            "step {n}: water {water:e} against {w0:e}"
        );
        assert!(
            (energy - e0).abs() < 1e-12 * e0,
            "step {n}: energy {energy:e} against {e0:e}"
        );
    }
    let (w0, e0) = first.unwrap();
    // the exact end state: air and vapour saturated at T, the rest liquid
    // at T, with the water and the energy of the start
    let end = |tt: f64| {
        let s = ys(tt);
        let mv = m_air * s / (1.0 - s);
        m_air * CP_AIR * tt + mv * hv(tt) + (w0 - mv) * c_l * tt - e0
    };
    let (mut lo, mut hi) = (280.0, 400.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if end(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let t_sat = 0.5 * (lo + hi);
    // the mass mean temperature of the liquid: each closing step gives a
    // droplet the vapour of its records' slopes back with the liquid
    // enthalpy and the share e_h of the latent heat (D250), so single
    // droplets scatter by about a kelvin around it while they keep
    // exchanging with the gas
    let tp = sp.enthalpy() / (sp.mass_in_flow_by_id() * c_l);
    let spread = sp
        .parcels
        .temp
        .iter()
        .map(|x| (x - t_sat).abs())
        .fold(0.0, f64::max);
    assert!(spread < 2.0, "a droplet {spread} K from the end state");
    assert!(!sp.parcels.is_empty(), "liquid is left");
    // the liquid holds a quarter of the heat capacity of the gas: a part
    // in a thousand of the vapour moves its mean by a tenth of a kelvin
    assert!(
        (t - t_sat).abs() < 0.01,
        "gas {t} K, adiabatic saturation {t_sat} K"
    );
    assert!((tp - t_sat).abs() < 0.1, "liquid {tp} K, {t_sat} K");
    assert!(
        (y - ys(t_sat)).abs() < 1e-3 * ys(t_sat),
        "Y {y}, saturation {}",
        ys(t_sat)
    );
    eprintln!(
        "closed box: gas {t:.4} K, liquid {tp:.4} K (droplets within {spread:.2} K), Y \
         {y:.6e}; adiabatic saturation {t_sat:.4} K, Y_s {:.6e}",
        ys(t_sat)
    );
}
