//! One-dimensional water and steam tube circuits of a boiler (steam furnace
//! stage SF7b of the
//! [development path](https://github.com/chiefenne/Foehncraft/blob/main/docs/assessments/steam-furnace-development-path.md)).
//!
//! A circuit is a graph of tubes between nodes: headers (junctions), fixed
//! pressures, mass-flow inlets and steam drums. Each tube is marched in
//! segments with the IAPWS-IF97 states of `foehncraft_kernels::steam`
//! (stage SF7a): mass, momentum (friction with the Friedel two-phase
//! multiplier or the homogeneous model, gravity with the drift-flux or the
//! homogeneous void, acceleration) and energy in the specific enthalpy. The
//! wall temperature follows from the heat flux by Gnielinski (single
//! phase) and Chen or Gungor–Winterton (subcooled and saturated flow
//! boiling); the void of subcooled boiling starts at the Saha–Zuber point
//! with the Levy profile of the true quality; the critical heat flux comes
//! from the Groeneveld look-up table (2006) and its margin is reported per
//! segment. The flows, pressures and enthalpies of the whole graph, the
//! natural circulation of a drum loop included, are solved together by a
//! damped Newton method.
pub mod case;
pub mod chf;
pub mod correlations;
pub mod coupling;
pub mod fluid;
pub mod linalg;
pub mod network;
pub mod output;

/// The diagnostic codes of the circuits (registered with the other codes
/// in `foehncraft_flow::diag`, D78).
pub mod codes {
    /// A circuit case that is not valid.
    pub const CASE: &str = "CIR001";
    /// The Newton solve of the circuit did not converge.
    pub const CONVERGENCE: &str = "CIR002";
    /// A state left the range of the water and steam properties.
    pub const RANGE: &str = "CIR003";
    /// The heat flux reaches the critical heat flux (margin below the
    /// limit of the case).
    pub const CHF: &str = "CIR004";
}

/// An error with its diagnostic code ([`codes`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Error {
    /// Diagnostic code.
    pub code: &'static str,
    /// What went wrong, for a person.
    pub message: String,
}

impl Error {
    /// A new error.
    pub fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}

impl std::error::Error for Error {}

/// Standard gravity (m/s²).
pub const G: f64 = 9.80665;
