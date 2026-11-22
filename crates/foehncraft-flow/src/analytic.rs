//! Exact solutions and reference data of the stage-4 verification and
//! validation cases (D45):
//!
//! - uniform flow (the exact solution between slip walls, D29);
//! - plane Couette flow between two walls, the upper one moving in x;
//! - plane Poiseuille flow between two walls at rest, driven by a pressure
//!   drop in x;
//! - Hagen–Poiseuille flow in a pipe along the x axis (the axisymmetric
//!   check on wedge meshes, D29);
//! - the lid-driven cavity centre-line velocities of Ghia, Ghia & Shin
//!   (1982), tables I and II, Re = 100 and 1000;
//! - the decaying Taylor–Green vortex in two dimensions (Taylor & Green
//!   1937), an exact transient solution of the Navier–Stokes equations
//!   (stage 5);
//! - the reference ranges of the flow around a cylinder in a channel of
//!   Schäfer & Turek (1996), cases 2D-1 (Re 20, steady) and 2D-2 (Re 100,
//!   periodic).

use serde::{Deserialize, Serialize};

/// An exact solution, `[analytic]`. The walls are the planes y = `lower`
/// and y = `upper`; the flow is in x.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(tag = "type", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Analytic {
    /// Uniform flow: u = U, p constant (exact between slip walls).
    Uniform {
        /// The velocity U (m/s).
        #[serde(default = "unit_x")]
        velocity: [f64; 3],
        /// The constant pressure (Pa).
        #[serde(default)]
        pressure: f64,
        /// Length scale of the Reynolds number (m).
        #[serde(default = "one")]
        length: f64,
    },
    /// Plane Couette flow: u = U (y − y₀)/(y₁ − y₀), p constant.
    Couette {
        /// y of the wall at rest.
        #[serde(default)]
        lower: f64,
        /// y of the moving wall.
        #[serde(default = "one")]
        upper: f64,
        /// Speed of the upper wall U (m/s).
        #[serde(default = "one")]
        wall_velocity: f64,
        /// The constant pressure (Pa).
        #[serde(default)]
        pressure: f64,
    },
    /// Plane Couette flow with the viscosity μ = μ₀ (1 + a (y − y₀)) (an
    /// expression of y, D53): the shear stress is constant, so
    /// u = U ln(1 + a (y − y₀)) / ln(1 + a H), p constant. Verifies the
    /// viscosity per face of the property program (D94).
    CouetteViscosity {
        /// y of the wall at rest.
        #[serde(default)]
        lower: f64,
        /// y of the moving wall.
        #[serde(default = "one")]
        upper: f64,
        /// Speed of the upper wall U (m/s).
        #[serde(default = "one")]
        wall_velocity: f64,
        /// The relative viscosity gradient a (1/m).
        gradient: f64,
        /// The constant pressure (Pa).
        #[serde(default)]
        pressure: f64,
    },
    /// Plane Couette flow with transpiration: fluid enters through the
    /// lower wall and leaves through the upper one at the normal speed V
    /// (the wall mass flux of D124 W3), u = U (1 − e^{Re_w η})/(1 − e^{Re_w}),
    /// v = V, η = (y − y₀)/H, Re_w = V H/ν, p constant (an exact solution
    /// of the Navier–Stokes equations; Schlichting & Gersten, Boundary-Layer
    /// Theory, 2017, §5.2).
    CouetteTranspiration {
        /// y of the wall at rest (blowing).
        #[serde(default)]
        lower: f64,
        /// y of the moving wall (suction).
        #[serde(default = "one")]
        upper: f64,
        /// Speed of the upper wall U (m/s).
        #[serde(default = "one")]
        wall_velocity: f64,
        /// The transpiration speed V (m/s), > 0 from the lower to the
        /// upper wall.
        transpiration_velocity: f64,
        /// The constant pressure (Pa).
        #[serde(default)]
        pressure: f64,
        /// Kinematic viscosity ν (m²/s), set by the case from `[fluid]`.
        #[serde(default)]
        nu: f64,
    },
    /// Plane Poiseuille flow: u = 6 U_m η (1 − η), η = (y − y₀)/H,
    /// p = p_out + G (x_out − x) with G = 12 μ U_m / H².
    Poiseuille {
        /// y of the lower wall.
        #[serde(default)]
        lower: f64,
        /// y of the upper wall.
        #[serde(default = "one")]
        upper: f64,
        /// Mean velocity U_m (m/s).
        #[serde(default = "one")]
        mean_velocity: f64,
        /// Dynamic viscosity μ (Pa s), for the pressure gradient (set by the
        /// case from `[fluid]`).
        #[serde(default)]
        viscosity: f64,
        /// x of the outlet.
        #[serde(default = "one")]
        x_outlet: f64,
        /// Pressure at the outlet (Pa).
        #[serde(default)]
        pressure_outlet: f64,
    },
    /// Natural convection between two vertical walls of infinite height
    /// (stage 7b1): the hot wall at x = `lower`, the cold one at x =
    /// `upper`, gravity along −y, a linear temperature between them and the
    /// Boussinesq force about their mean temperature. The flow rises at the
    /// hot wall and falls at the cold one with no net flow:
    /// v = (g β ΔT D²/(48 ν)) (ξ³ − ξ), ξ = 2 (x − x₀)/D − 1, u = 0, and a
    /// constant reduced pressure (its level is not compared).
    BuoyantSlot {
        /// x of the hot wall.
        #[serde(default)]
        lower: f64,
        /// x of the cold wall.
        #[serde(default = "one")]
        upper: f64,
        /// g β ΔT (m/s²) of the wall temperature difference ΔT.
        buoyancy: f64,
        /// Kinematic viscosity ν (m²/s), set by the case from `[fluid]`.
        #[serde(default)]
        nu: f64,
    },
    /// The manufactured solution of the variable density (stage 7c, check
    /// V7c.4) on the unit square, an ideal gas with R/M = 1 at p₀ = 1:
    /// T = 1 + a sin πx sin πy, ρ = p₀/(R T/M) = 1/T, the mass flux from
    /// the stream function ψ = A sin²πx sin²πy (ρu = ∂ψ/∂y, ρv = −∂ψ/∂x,
    /// so ∇·(ρu) = 0 exactly and u = 0 on the walls), and p = P cos πx
    /// cos πy up to a constant. The momentum and heat sources that make it
    /// a solution are expressions of the case (written by
    /// `tools/verification/vd_mms.py`).
    VariableDensity {
        /// A (kg/(m s)), the amplitude of the stream function.
        flux: f64,
        /// a (K), the amplitude of the temperature.
        temperature: f64,
        /// P (Pa), the amplitude of the pressure.
        pressure: f64,
    },
    /// Hagen–Poiseuille flow in a pipe of radius R along the x axis:
    /// u = 2 U_m (1 − r²/R²), r² = y² + z², p = p_out + G (x_out − x) with
    /// G = 8 μ U_m / R².
    PipePoiseuille {
        /// Pipe radius R (m).
        #[serde(default = "one")]
        radius: f64,
        /// Mean velocity U_m (m/s).
        #[serde(default = "one")]
        mean_velocity: f64,
        /// Dynamic viscosity μ (Pa s), set by the case from `[fluid]`.
        #[serde(default)]
        viscosity: f64,
        /// x of the outlet.
        #[serde(default = "one")]
        x_outlet: f64,
        /// Pressure at the outlet (Pa).
        #[serde(default)]
        pressure_outlet: f64,
    },
    /// Laminar Taylor–Couette flow between coaxial cylinders about the z
    /// axis through the origin (a 2D case, or 3D with free ends): the
    /// azimuthal velocity u_θ = A r + B/r with
    /// A = (ω₂ r₂² − ω₁ r₁²)/(r₂² − r₁²), B = (ω₁ − ω₂) r₁² r₂²/(r₂² − r₁²),
    /// and the pressure from dp/dr = ρ u_θ²/r, p = ρ (A² r²/2 + 2 A B ln r
    /// − B²/(2 r²)) up to a constant (the level is free in a closed
    /// domain, so the error compares the volume means).
    TaylorCouette {
        /// Radius of the inner cylinder r₁ (m).
        #[serde(default = "inner_radius")]
        inner_radius: f64,
        /// Radius of the outer cylinder r₂ (m).
        #[serde(default = "outer_radius")]
        outer_radius: f64,
        /// Angular velocity of the inner cylinder ω₁ (rad/s).
        #[serde(default = "one")]
        inner_omega: f64,
        /// Angular velocity of the outer cylinder ω₂ (rad/s).
        #[serde(default)]
        outer_omega: f64,
        /// Density (kg/m³), set by the case from `[fluid]`.
        #[serde(skip, default)]
        rho: f64,
    },
    /// Decaying Taylor–Green vortex in the x–y plane, period L in x and y,
    /// k = 2π/L, F(t) = exp(−2 ν k² t):
    /// u = −U cos(kx) sin(ky) F, v = U sin(kx) cos(ky) F,
    /// p = p₀ − ¼ ρ U² (cos 2kx + cos 2ky) F².
    TaylorGreen {
        /// Velocity amplitude U (m/s).
        #[serde(default = "one")]
        velocity: f64,
        /// Period L in x and y (m).
        #[serde(default = "one")]
        length: f64,
        /// Mean pressure p₀ (Pa).
        #[serde(default)]
        pressure: f64,
        /// Kinematic viscosity ν (m²/s), set by the case from `[fluid]`.
        #[serde(skip, default)]
        nu: f64,
        /// Density (kg/m³), set by the case from `[fluid]`.
        #[serde(skip, default)]
        rho: f64,
        /// Time (s) at which the solution is evaluated, set by the solver.
        #[serde(skip, default)]
        time: f64,
    },
}

fn one() -> f64 {
    1.0
}

fn inner_radius() -> f64 {
    0.05
}

fn outer_radius() -> f64 {
    0.10
}

/// The constants A and B of u_θ = A r + B/r of Taylor–Couette flow.
fn taylor_couette_ab(r1: f64, r2: f64, w1: f64, w2: f64) -> (f64, f64) {
    let d = r2 * r2 - r1 * r1;
    (
        (w2 * r2 * r2 - w1 * r1 * r1) / d,
        (w1 - w2) * r1 * r1 * r2 * r2 / d,
    )
}

fn unit_x() -> [f64; 3] {
    [1.0, 0.0, 0.0]
}

impl Analytic {
    /// Exact velocity at `x`.
    pub fn velocity(&self, x: [f64; 3]) -> [f64; 3] {
        match *self {
            Self::Uniform { velocity, .. } => velocity,
            Self::Couette {
                lower,
                upper,
                wall_velocity,
                ..
            } => [wall_velocity * (x[1] - lower) / (upper - lower), 0.0, 0.0],
            Self::CouetteViscosity {
                lower,
                upper,
                wall_velocity,
                gradient,
                ..
            } => [
                wall_velocity * (gradient * (x[1] - lower)).ln_1p()
                    / (gradient * (upper - lower)).ln_1p(),
                0.0,
                0.0,
            ],
            Self::CouetteTranspiration {
                lower,
                upper,
                wall_velocity,
                transpiration_velocity: v,
                nu,
                ..
            } => {
                let h = upper - lower;
                let re = v * h / nu;
                let eta = (x[1] - lower) / h;
                [wall_velocity * (re * eta).exp_m1() / re.exp_m1(), v, 0.0]
            }
            Self::Poiseuille {
                lower,
                upper,
                mean_velocity,
                ..
            } => {
                let eta = (x[1] - lower) / (upper - lower);
                [6.0 * mean_velocity * eta * (1.0 - eta), 0.0, 0.0]
            }
            Self::PipePoiseuille {
                radius,
                mean_velocity,
                ..
            } => {
                let r2 = (x[1] * x[1] + x[2] * x[2]) / (radius * radius);
                [2.0 * mean_velocity * (1.0 - r2), 0.0, 0.0]
            }
            Self::BuoyantSlot {
                lower,
                upper,
                buoyancy,
                nu,
            } => {
                let d = upper - lower;
                let xi = 2.0 * (x[0] - lower) / d - 1.0;
                [
                    0.0,
                    buoyancy * d * d / (48.0 * nu) * (xi * xi * xi - xi),
                    0.0,
                ]
            }
            Self::TaylorCouette {
                inner_radius: r1,
                outer_radius: r2,
                inner_omega: w1,
                outer_omega: w2,
                ..
            } => {
                let (a, b) = taylor_couette_ab(r1, r2, w1, w2);
                let r2x = x[0] * x[0] + x[1] * x[1];
                // u_θ/r = A + B/r²
                let g = a + b / r2x;
                [-g * x[1], g * x[0], 0.0]
            }
            Self::VariableDensity {
                flux, temperature, ..
            } => {
                use std::f64::consts::PI;
                let (sx, sy) = ((PI * x[0]).sin(), (PI * x[1]).sin());
                let t = 1.0 + temperature * sx * sy;
                [
                    t * flux * PI * sx * sx * (2.0 * PI * x[1]).sin(),
                    -t * flux * PI * (2.0 * PI * x[0]).sin() * sy * sy,
                    0.0,
                ]
            }
            Self::TaylorGreen {
                velocity,
                length,
                nu,
                time,
                ..
            } => {
                let k = 2.0 * std::f64::consts::PI / length;
                let f = (-2.0 * nu * k * k * time).exp();
                let (kx, ky) = (k * x[0], k * x[1]);
                [
                    -velocity * kx.cos() * ky.sin() * f,
                    velocity * kx.sin() * ky.cos() * f,
                    0.0,
                ]
            }
        }
    }

    /// Exact pressure at `x`.
    pub fn pressure(&self, x: [f64; 3]) -> f64 {
        match *self {
            Self::Uniform { pressure, .. }
            | Self::Couette { pressure, .. }
            | Self::CouetteViscosity { pressure, .. }
            | Self::CouetteTranspiration { pressure, .. } => pressure,
            Self::BuoyantSlot { .. } => 0.0,
            Self::VariableDensity { pressure, .. } => {
                use std::f64::consts::PI;
                pressure * (PI * x[0]).cos() * (PI * x[1]).cos()
            }
            Self::Poiseuille {
                lower,
                upper,
                mean_velocity,
                viscosity,
                x_outlet,
                pressure_outlet,
            } => {
                let h = upper - lower;
                pressure_outlet + 12.0 * viscosity * mean_velocity / (h * h) * (x_outlet - x[0])
            }
            Self::PipePoiseuille {
                radius,
                mean_velocity,
                viscosity,
                x_outlet,
                pressure_outlet,
            } => {
                pressure_outlet
                    + 8.0 * viscosity * mean_velocity / (radius * radius) * (x_outlet - x[0])
            }
            Self::TaylorCouette {
                inner_radius: r1,
                outer_radius: r2,
                inner_omega: w1,
                outer_omega: w2,
                rho,
            } => {
                let (a, b) = taylor_couette_ab(r1, r2, w1, w2);
                let r = (x[0] * x[0] + x[1] * x[1]).sqrt();
                rho * (0.5 * a * a * r * r + 2.0 * a * b * r.ln() - 0.5 * b * b / (r * r))
            }
            Self::TaylorGreen {
                velocity,
                length,
                pressure,
                nu,
                rho,
                time,
            } => {
                let k = 2.0 * std::f64::consts::PI / length;
                let f = (-2.0 * nu * k * k * time).exp();
                pressure
                    - 0.25
                        * rho
                        * velocity
                        * velocity
                        * ((2.0 * k * x[0]).cos() + (2.0 * k * x[1]).cos())
                        * f
                        * f
            }
        }
    }

    /// True for a solution that changes in time.
    pub fn is_transient(&self) -> bool {
        matches!(self, Self::TaylorGreen { .. })
    }

    /// True when the pressure of the exact solution is known only up to a
    /// constant that the run does not fix (closed domain): the error then
    /// compares the volume means.
    pub fn free_pressure_level(&self) -> bool {
        matches!(
            self,
            Self::TaylorCouette { .. } | Self::BuoyantSlot { .. } | Self::VariableDensity { .. }
        )
    }

    /// The solution at time `t` (steady solutions ignore it).
    pub fn at_time(mut self, t: f64) -> Self {
        if let Self::TaylorGreen { time, .. } = &mut self {
            *time = t;
        }
        self
    }

    /// Sets the fluid properties the solution depends on.
    pub fn with_fluid(mut self, density: f64, viscosity: f64) -> Self {
        match &mut self {
            Self::Poiseuille { viscosity: v, .. } | Self::PipePoiseuille { viscosity: v, .. } => {
                *v = viscosity
            }
            Self::TaylorGreen { nu, rho, .. } => {
                *nu = viscosity / density;
                *rho = density;
            }
            Self::TaylorCouette { rho, .. } => *rho = density,
            Self::CouetteTranspiration { nu, .. } | Self::BuoyantSlot { nu, .. } => {
                *nu = viscosity / density
            }
            Self::Uniform { .. }
            | Self::Couette { .. }
            | Self::CouetteViscosity { .. }
            | Self::VariableDensity { .. } => {}
        }
        self
    }

    /// Length scale of the Reynolds number (channel height or period).
    pub fn length_scale(&self) -> f64 {
        match *self {
            Self::Couette { lower, upper, .. }
            | Self::CouetteViscosity { lower, upper, .. }
            | Self::CouetteTranspiration { lower, upper, .. }
            | Self::Poiseuille { lower, upper, .. }
            | Self::BuoyantSlot { lower, upper, .. } => upper - lower,
            Self::PipePoiseuille { radius, .. } => 2.0 * radius,
            Self::VariableDensity { .. } => 1.0,
            Self::Uniform { length, .. } => length,
            Self::TaylorGreen { length, .. } => length,
            Self::TaylorCouette {
                inner_radius,
                outer_radius,
                ..
            } => outer_radius - inner_radius,
        }
    }

    /// Velocity scale for the error norms (the wall or the largest speed).
    pub fn velocity_scale(&self) -> f64 {
        match *self {
            Self::Uniform { velocity: v, .. } => (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt(),
            Self::Couette { wall_velocity, .. }
            | Self::CouetteViscosity { wall_velocity, .. }
            | Self::CouetteTranspiration { wall_velocity, .. } => wall_velocity.abs(),
            Self::Poiseuille { mean_velocity, .. } => 1.5 * mean_velocity.abs(),
            Self::PipePoiseuille { mean_velocity, .. } => 2.0 * mean_velocity.abs(),
            Self::BuoyantSlot {
                lower,
                upper,
                buoyancy,
                nu,
            } => {
                // the largest speed, at ξ = ∓1/√3
                let d = upper - lower;
                (buoyancy * d * d / (48.0 * nu)).abs() * 2.0 / (3.0 * 3f64.sqrt())
            }
            Self::TaylorGreen { velocity, .. } => velocity.abs(),
            // the largest |ρu| = π A over the smallest ρ = 1/(1 + a)
            Self::VariableDensity {
                flux, temperature, ..
            } => std::f64::consts::PI * flux.abs() * (1.0 + temperature.abs()),
            Self::TaylorCouette {
                inner_radius,
                outer_radius,
                inner_omega,
                outer_omega,
                ..
            } => (inner_omega * inner_radius)
                .abs()
                .max((outer_omega * outer_radius).abs()),
        }
    }

    /// Sets the viscosity of the pressure gradient.
    pub fn with_viscosity(mut self, mu: f64) -> Self {
        if let Self::Poiseuille { viscosity, .. } | Self::PipePoiseuille { viscosity, .. } =
            &mut self
        {
            *viscosity = mu;
        }
        self
    }
}

/// Ghia et al. (1982), table I: y on the vertical centre line x = 0.5.
pub const GHIA_Y: [f64; 17] = [
    1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719,
    0.1016, 0.0703, 0.0625, 0.0547, 0.0000,
];
/// Table I: u / U_lid at `GHIA_Y`, Re = 100.
pub const GHIA_U_RE100: [f64; 17] = [
    1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090,
    -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000,
];
/// Table I: u / U_lid at `GHIA_Y`, Re = 1000.
pub const GHIA_U_RE1000: [f64; 17] = [
    1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648,
    -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000,
];
/// Ghia et al. (1982), table II: x on the horizontal centre line y = 0.5.
pub const GHIA_X: [f64; 17] = [
    1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563,
    0.0938, 0.0781, 0.0703, 0.0625, 0.0000,
];
/// Table II: v / U_lid at `GHIA_X`, Re = 100.
pub const GHIA_V_RE100: [f64; 17] = [
    0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454,
    0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000,
];
/// Table II: v / U_lid at `GHIA_X`, Re = 1000.
pub const GHIA_V_RE1000: [f64; 17] = [
    0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526,
    0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000,
];

/// Reference data of a case, `[verification] reference`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Reference {
    /// Ghia et al. 1982, Re = 100.
    GhiaRe100,
    /// Ghia et al. 1982, Re = 1000.
    GhiaRe1000,
}

impl Reference {
    /// (y, u) on the vertical and (x, v) on the horizontal centre line.
    pub fn centre_lines(self) -> ([f64; 17], [f64; 17], [f64; 17], [f64; 17]) {
        match self {
            Self::GhiaRe100 => (GHIA_Y, GHIA_U_RE100, GHIA_X, GHIA_V_RE100),
            Self::GhiaRe1000 => (GHIA_Y, GHIA_U_RE1000, GHIA_X, GHIA_V_RE1000),
        }
    }

    /// Reynolds number of the data.
    pub fn reynolds(self) -> f64 {
        match self {
            Self::GhiaRe100 => 100.0,
            Self::GhiaRe1000 => 1000.0,
        }
    }

    /// Label for reports.
    pub fn label(self) -> &'static str {
        match self {
            Self::GhiaRe100 => "Ghia et al. (1982), Re = 100",
            Self::GhiaRe1000 => "Ghia et al. (1982), Re = 1000",
        }
    }
}

/// A reference value of the cylinder benchmark: the range of the
/// participants that Schäfer & Turek (1996) give as the reference
/// solution.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Range {
    /// Quantity.
    pub name: &'static str,
    /// Lower bound.
    pub lower: f64,
    /// Upper bound.
    pub upper: f64,
}

impl Range {
    /// Relative distance of `v` from the range (0 inside), over the
    /// midpoint.
    pub fn deviation(&self, v: f64) -> f64 {
        let mid = 0.5 * (self.lower + self.upper);
        let d = if v < self.lower {
            self.lower - v
        } else if v > self.upper {
            v - self.upper
        } else {
            0.0
        };
        d / mid.abs().max(f64::MIN_POSITIVE)
    }
}

/// Benchmark data of a case, `[verification] benchmark`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum Benchmark {
    /// Schäfer & Turek (1996) 2D-1: channel 2.2 × 0.41 m, cylinder D = 0.1
    /// m at (0.2, 0.2), parabolic inflow with U_max = 0.3 m/s, ν = 10⁻³
    /// m²/s, Re = Ū D / ν = 20, steady.
    #[serde(rename = "schafer-turek-2d1")]
    SchaferTurek2d1,
    /// Schäfer & Turek (1996) 2D-2: as 2D-1 with U_max = 1.5 m/s, Re = 100,
    /// periodic vortex shedding.
    #[serde(rename = "schafer-turek-2d2")]
    SchaferTurek2d2,
}

impl Benchmark {
    /// Mean inflow velocity Ū = 2/3 U_max (m/s).
    pub fn mean_velocity(self) -> f64 {
        match self {
            Self::SchaferTurek2d1 => 0.2,
            Self::SchaferTurek2d2 => 1.0,
        }
    }

    /// Reynolds number Ū D / ν.
    pub fn reynolds(self) -> f64 {
        match self {
            Self::SchaferTurek2d1 => 20.0,
            Self::SchaferTurek2d2 => 100.0,
        }
    }

    /// Cylinder diameter (m).
    pub const DIAMETER: f64 = 0.1;

    /// The points of the pressure difference Δp = p(front) − p(back).
    pub const FRONT: [f64; 2] = [0.15, 0.2];
    /// See [`Self::FRONT`].
    pub const BACK: [f64; 2] = [0.25, 0.2];

    /// Reference ranges (Schäfer & Turek 1996, tables 3 and 4): 2D-1
    /// C_D, C_L, Δp; 2D-2 C_D,max, C_L,max, St, Δp at t₀ + 1/(2f) (t₀ the
    /// time of C_L,max).
    pub fn ranges(self) -> Vec<Range> {
        let r = |name, lower, upper| Range { name, lower, upper };
        match self {
            Self::SchaferTurek2d1 => vec![
                r("C_D", 5.57, 5.59),
                r("C_L", 0.0104, 0.0110),
                r("Δp (Pa)", 0.1172, 0.1176),
            ],
            Self::SchaferTurek2d2 => vec![
                r("C_D,max", 3.22, 3.24),
                r("C_L,max", 0.99, 1.01),
                r("St", 0.295, 0.305),
                r("Δp (Pa)", 2.46, 2.50),
            ],
        }
    }

    /// Label for reports.
    pub fn label(self) -> &'static str {
        match self {
            Self::SchaferTurek2d1 => "Schäfer & Turek (1996) 2D-1, Re = 20",
            Self::SchaferTurek2d2 => "Schäfer & Turek (1996) 2D-2, Re = 100",
        }
    }
}
