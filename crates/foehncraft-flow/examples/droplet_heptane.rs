//! Single n-heptane droplets in still nitrogen at 0.1 MPa through the
//! evaporation of the spray solver (D320, D321): the data of the gallery
//! image, against the measurements of Chauveau, Birouk, Halter and Gökalp
//! (2019) and Chauveau, Halter, Lalonde and Gökalp (2008).
//!
//! ```text
//! cargo run --release -p foehncraft-flow --example droplet_heptane -- output/droplet_heptane
//! python3 -m tools.plot.droplet_heptane output/droplet_heptane docs/figures/gallery
//! ```
//!
//! Writes `history.csv` (t/d0², d²/d0², droplet temperature per gas
//! temperature and film state) and `rates.csv` (K of the slope of d²/d0²
//! between 0.8 and 0.3 for both film states, the microgravity fit and the
//! fibre-free extrapolation of the measurements).

use std::fmt::Write as _;

use anyhow::{Context, Result};
use foehncraft_flow::species::database;
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, Env, EvaporationModel, Film, GasState, Spray, SpraySection, Vapour,
};

const P0: f64 = 101_325.0;
const D0: f64 = 370e-6;

fn spray(env: &Env, film: Film) -> Result<Spray> {
    let mut cfg = SpraySection::default();
    cfg.liquid.vapour = Some("NC7H16".into());
    cfg.liquid.library = Some("n-heptane".into());
    cfg.liquid.apply_library()?;
    cfg.evaporation.model = EvaporationModel::AbramzonSirignano;
    cfg.evaporation.film = film;
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Rebound,
            e_normal: 0.5,
            e_tangential: 0.5,
        };
        env.mesh.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, env, rules, [0.0; 3])?;
    let db = database();
    sp.set_vapour(Some(Vapour {
        species: db["NC7H16"].kernel(),
    }));
    sp.set_carrier(vec![(db["N2"].kernel(), 1.0)]);
    Ok(sp)
}

/// The history (t, d²/d0², T_p) of a droplet at 300 K in N2 at `tg`.
fn history(sp: &Spray, tg: f64) -> Vec<(f64, f64, f64)> {
    let db = database();
    let n2 = db["N2"].kernel();
    let mu = n2.viscosity(tg);
    let gs = GasState {
        u: [0.0; 3],
        rho: P0 * n2.molar_mass / (8.314_462_618 * tg),
        mu,
        t: tg,
        cp: n2.cp(tg),
        k: n2.conductivity(tg, mu),
        tke: 0.0,
        eps: 0.0,
        y_v: 0.0,
        p0: P0,
        m_gas: n2.molar_mass,
        diff: foehncraft_kernels::species::binary_diffusion(&db["NC7H16"].kernel(), &n2, tg, P0),
    };
    let (mut d, mut tp, mut t) = (D0, 300.0, 0.0);
    let mut h = vec![(0.0, 1.0, tp)];
    while d > 0.02 * D0 {
        let dt = 1e-4 * (d / D0).powi(2) + 1e-7;
        let Some(s) = sp.evaporation(&gs, d, tp, 0.0, dt) else {
            break;
        };
        d = s.d;
        tp = s.tp;
        t += dt;
        h.push((t, (d / D0).powi(2), tp));
        if d == 0.0 {
            break;
        }
    }
    h
}

fn slope(h: &[(f64, f64, f64)]) -> f64 {
    let a = h.iter().find(|x| x.1 < 0.8).unwrap();
    let b = h.iter().find(|x| x.1 < 0.3).unwrap();
    (a.1 - b.1) / (b.0 - a.0) * D0 * D0 * 1e6
}

fn main() -> Result<()> {
    let out = std::env::args().nth(1).unwrap_or_else(|| "output/droplet_heptane".into());
    std::fs::create_dir_all(&out).with_context(|| format!("create {out}"))?;
    let mut s = BoxSpec::new(Kind::Hex, 3, [1, 1, 1]);
    s.lower = [0.0; 3];
    s.upper = [1.0; 3];
    let mesh = s.generate().map_err(anyhow::Error::msg)?;
    let geo = Geometry::compute(&mesh);
    let loc = Locator::new(&mesh, &geo);
    let env = Env { mesh: &mesh, geo: &geo, loc: &loc };
    let films = [("one-third", spray(&env, Film::OneThird)?), ("far-field", spray(&env, Film::FarField)?)];
    let mut hist = String::from("T_gas,film,t_over_d0sq_s_per_mm2,d2_over_d02,T_drop\n");
    let mut rates = String::from("T_gas,K_one_third,K_far_field,K_microgravity,K_fibre_free\n");
    for tg in [473.0, 548.0, 623.0, 698.0, 748.0] {
        let mut k = [0.0; 2];
        for (j, (name, sp)) in films.iter().enumerate() {
            let h = history(sp, tg);
            k[j] = slope(&h);
            for (i, x) in h.iter().enumerate() {
                if i % 20 == 0 || i + 1 == h.len() {
                    writeln!(hist, "{tg},{name},{:.6},{:.6},{:.3}", x.0 / (D0 * D0 * 1e6), x.1, x.2)?;
                }
            }
        }
        writeln!(
            rates,
            "{tg},{:.5},{:.5},{:.5},{:.5}",
            k[0],
            k[1],
            3.6552e-4 * tg - 0.1078,
            4.485e-4 * tg - 0.1369
        )?;
    }
    std::fs::write(format!("{out}/history.csv"), hist)?;
    std::fs::write(format!("{out}/rates.csv"), rates)?;
    print!("{rates}");
    Ok(())
}
