//! A small linear structural finite element solver (L53, stage MS5 of the
//! [moving solids plan](https://github.com/chiefenne/Foehncraft/blob/main/docs/assessments/moving-solids-reed-valves-plan.md)).
pub mod analysis;
pub mod beam;
pub mod case;
pub mod clamp;
pub mod dynamics;
pub mod eigen;
pub mod ldlt;
pub mod model;
pub mod shell;
pub mod solid;
pub mod solid_mesh;
pub mod sparse;

/// The diagnostic codes of the structural solver (registered with the
/// other codes in `foehncraft_flow::diag`, D78).
pub mod codes {
    /// A structural case or model that is not valid.
    pub const CASE: &str = "STR001";
    /// A singular stiffness: the structure can move freely.
    pub const SINGULAR: &str = "STR002";
    /// The modal analysis may have missed modes or did not converge.
    pub const MODES: &str = "STR003";
    /// The clamp fit cannot reach the measured frequency.
    pub const FIT: &str = "STR004";
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
    pub fn new(code: &'static str, message: String) -> Self {
        Self { code, message }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}

impl std::error::Error for Error {}
