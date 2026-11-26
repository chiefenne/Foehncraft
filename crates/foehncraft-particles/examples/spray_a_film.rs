//! The real-fluid gas film of an n-dodecane droplet in nitrogen at high
//! pressure (ECN Spray A, SA2, D313): the data of the gallery image.
//!
//! ```text
//! cargo run --release -p foehncraft-particles --example spray_a_film -- output/spray_a_film
//! python3 -m tools.plot.spray_a_film output/spray_a_film docs/figures/gallery
//! ```
//!
//! Writes into the output directory
//!
//! - `nitrogen.csv`: the dense-fluid factors η/η° and λ/λ° of Chung et al.
//!   (1988) for nitrogen at 6 and 15 MPa over the temperature, with the
//!   Peng-Robinson density, and `nist.csv` the NIST ratios to 0.1 MPa at
//!   the same temperatures (the points of `tests/dense_film.rs`);
//! - `film.csv`: the film of a droplet in 900 K nitrogen at 6 and 15 MPa
//!   over the droplet temperature (fuel fraction of the Peng-Robinson
//!   surface, 1/3 rule): density, viscosity, conductivity, c_p and
//!   diffusivity of the real-fluid film over those of the ideal-gas film;
//! - `droplet.csv`: the diameter and temperature of a 10 µm droplet of
//!   363 K n-dodecane in still 900 K nitrogen at 6 MPa with either film.

use std::fmt::Write as _;

use anyhow::{Context, Result};
use foehncraft_kernels::realfluid::{MAX, Mixture, Root, fuels::N_DODECANE, gases};
use foehncraft_kernels::species::{R, Species};
use foehncraft_mesh::Geometry;
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::locate::Locator;
use foehncraft_particles::{
    BoundaryMode, BoundaryRule, Env, Equilibrium, EvaporationModel, FilmFluid, Gas, GasState,
    Injector, Spray, SpraySection, Vapour, ideal_gas_species,
};

/// Nitrogen, NIST Chemistry WebBook fluid tables (Span et al. 2000;
/// Lemmon and Jacobsen 2004): T (K), p (MPa), η (Pa s), λ (W/(m K)).
const NIST: [(f64, f64, f64, f64); 11] = [
    (300.0, 0.1, 1.78899e-05, 0.0259682),
    (600.0, 0.1, 2.95770e-05, 0.0448405),
    (900.0, 0.1, 3.87800e-05, 0.0605199),
    (1200.0, 0.1, 4.67665e-05, 0.0746785),
    (300.0, 6.0, 1.89442e-05, 0.0287804),
    (600.0, 6.0, 2.99556e-05, 0.0460242),
    (900.0, 6.0, 3.89937e-05, 0.0612592),
    (1200.0, 6.0, 4.69097e-05, 0.0752132),
    (300.0, 15.0, 2.14870e-05, 0.0342397),
    (600.0, 15.0, 3.07601e-05, 0.0482899),
    (900.0, 15.0, 3.94277e-05, 0.0626566),
];

/// Still nitrogen at 900 K and `p0`: NIST at 6 MPa for the far field.
struct Still {
    p0: f64,
}

impl Gas for Still {
    fn at(&self, _: usize, _: [f64; 3]) -> GasState {
        GasState {
            u: [0.0; 3],
            rho: 21.97 * self.p0 / 6e6,
            mu: 3.899e-5,
            t: 900.0,
            cp: 1152.8,
            k: 0.0613,
            tke: 0.0,
            eps: 0.0,
            y_v: 0.0,
            p0: self.p0,
            m_gas: 0.028_013_4,
            diff: 2.5e-5,
        }
    }
}

fn nitrogen() -> Result<Species> {
    ideal_gas_species(
        include_str!("../../foehncraft-flow/data/species-gri30.toml"),
        "N2",
    )
}

/// The spray of one 10 µm droplet of 363 K n-dodecane with the
/// Peng-Robinson surface at `p0` and the film `film`.
fn spray(env: &Env, p0: f64, film: FilmFluid) -> Result<Spray> {
    let mut cfg = SpraySection {
        injectors: vec![Injector {
            name: "n".into(),
            position: [0.005; 3],
            direction: [1.0, 0.0, 0.0],
            speed: 1e-4,
            // one parcel of one droplet
            mass_flow: 703.7 * std::f64::consts::PI / 6.0 * 1e-15 / 1e-6,
            diameter: 10e-6,
            parcels_per_second: 1e6,
            end: Some(1e-6),
            temperature: Some(363.0),
            ..Injector::default()
        }],
        ..SpraySection::default()
    };
    cfg.liquid.library = Some("n-dodecane".into());
    cfg.liquid.temperature = 363.0;
    cfg.liquid.vapour = Some("NC12H26".into());
    cfg.liquid.apply_library()?;
    cfg.evaporation.model = EvaporationModel::AbramzonSirignano;
    cfg.evaporation.equilibrium = Equilibrium::PengRobinson;
    cfg.evaporation.peng_robinson.film = film;
    let rules = vec![
        BoundaryRule {
            mode: BoundaryMode::Rebound,
            e_normal: 0.5,
            e_tangential: 0.5,
        };
        env.mesh.boundaries.len()
    ];
    let mut sp = Spray::new(&cfg, env, rules, [0.0; 3])?;
    sp.set_vapour(Some(Vapour {
        species: N_DODECANE.ideal_gas,
    }));
    sp.set_carrier(vec![(nitrogen()?, 1.0)]);
    sp.set_surface_pressure(p0)?;
    Ok(sp)
}

fn main() -> Result<()> {
    let out = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "output/spray_a_film".into());
    std::fs::create_dir_all(&out).with_context(|| format!("create {out}"))?;
    let n2 = nitrogen()?;
    let m = gases::N2.molar_mass;

    // nitrogen: Chung's factors along the Peng-Robinson isobars
    let pure = Mixture::new(&[gases::N2]);
    let mut x = [0.0; MAX];
    x[0] = 1.0;
    let fluid = pure.chung(&x);
    let mut s = String::from("p_MPa,T,rho,eta_factor,lambda_factor\n");
    for p in [6e6, 15e6] {
        for k in 0..=36 {
            let t = 300.0 + 25.0 * f64::from(k);
            let cv = n2.cp(t) * m - R;
            let rm = 1.0 / pure.phase(&x, t, p, Root::Vapour).v;
            // to 0.1 MPa, as the NIST ratios
            let r0 = 1.0 / pure.phase(&x, t, 1e5, Root::Vapour).v;
            let fe = fluid.viscosity_factor(t, rm) / fluid.viscosity_factor(t, r0);
            let fl = fluid.conductivity_factor(t, rm, cv) / fluid.conductivity_factor(t, r0, cv);
            writeln!(s, "{},{t},{},{fe},{fl}", p / 1e6, rm * m)?;
        }
    }
    std::fs::write(format!("{out}/nitrogen.csv"), s)?;
    let mut s = String::from("p_MPa,T,eta_ratio,lambda_ratio\n");
    for (t, p, eta, lambda) in NIST {
        if p > 0.1 {
            let lo = NIST.iter().find(|r| r.0 == t && r.1 == 0.1).unwrap();
            writeln!(s, "{p},{t},{},{}", eta / lo.2, lambda / lo.3)?;
        }
    }
    std::fs::write(format!("{out}/nist.csv"), s)?;

    // the film of a droplet, real fluid over ideal gas
    let mut b = BoxSpec::new(Kind::Hex, 3, [1, 1, 1]);
    b.lower = [0.0; 3];
    b.upper = [0.01; 3];
    let mesh = b.generate().map_err(anyhow::Error::msg)?;
    let geo = Geometry::compute(&mesh);
    let loc = Locator::new(&mesh, &geo);
    let env = Env {
        mesh: &mesh,
        geo: &geo,
        loc: &loc,
    };
    let mf = N_DODECANE.component.molar_mass;
    let mut s = String::from("p_MPa,T_drop,Y_s,rho,eta,lambda,cp,D,T_end\n");
    for p in [6e6, 15e6] {
        let real = spray(&env, p, FilmFluid::RealFluid)?;
        let ideal = spray(&env, p, FilmFluid::IdealGas)?;
        let gas = Still { p0: p }.at(0, [0.0; 3]);
        let t_end = real.surface_table().context("table")?.t_end;
        let mut t = 370.0;
        while t < t_end {
            let st = real.surface_state(t).context("surface")?;
            let ys = st.y_fuel * mf / (st.y_fuel * mf + (1.0 - st.y_fuel) * m);
            let a = real.film_properties(&gas, t, ys);
            let c = ideal.film_properties(&gas, t, ys);
            write!(s, "{},{t},{ys}", p / 1e6)?;
            for k in 0..5 {
                write!(s, ",{}", a[k] / c[k])?;
            }
            writeln!(s, ",{t_end}")?;
            t += 5.0;
        }
        println!("{:4.0} MPa: surface table ends at {t_end:.1} K", p / 1e6);
    }
    std::fs::write(format!("{out}/film.csv"), s)?;

    // one droplet with either film
    let mut s = String::from("film,t_ms,d_um,T\n");
    let gas = Still { p0: 6e6 };
    for (name, film) in [
        ("real-fluid", FilmFluid::RealFluid),
        ("ideal-gas", FilmFluid::IdealGas),
    ] {
        let mut sp = spray(&env, 6e6, film)?;
        let dt = 2e-6;
        let mut t = 0.0;
        let mut life = f64::NAN;
        for _ in 0..20_000 {
            sp.advance(&env, &gas, t, dt)?;
            t += dt;
            let p = &sp.parcels;
            if p.d.is_empty() {
                if t > 1e-5 {
                    life = t;
                    break;
                }
                continue;
            }
            writeln!(s, "{name},{},{},{}", t * 1e3, p.d[0] * 1e6, p.temp[0])?;
        }
        println!("{name:10}: the droplet lives {:.3} ms", life * 1e3);
    }
    std::fs::write(format!("{out}/droplet.csv"), s)?;
    println!("wrote {out}/nitrogen.csv, nist.csv, film.csv, droplet.csv");
    Ok(())
}
