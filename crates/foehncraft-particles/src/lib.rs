//! Lagrangian parcels of the sprays (L12, version 1, D240 to D246).
//!
//! - [`batch`]: the parcel step for a batch of parcels, the CPU backends of the
//!   kernel that also runs on the GPU (D245);
//! - [`config`]: the `[spray]` table and its `include` file;
//! - [`equilibrium`]: the high-pressure phase equilibrium at the droplet
//!   surface, tabulated over the droplet temperature (Spray A, SA1, D311);
//! - [`datum`]: the datums of the sprays, the crossing records and the
//!   droplet statistics from moment sums (SP1b, D262 to D265);
//! - [`sizes`]: the droplet size distributions of the injectors (D242);
//! - [`spray`]: the parcels in structure-of-arrays form, the injectors,
//!   the tracking of each parcel through the polyhedral mesh (face walk of
//!   Macpherson, Nordin and Weller 2009, `foehncraft_mesh::locate`) with
//!   the exact drag relaxation of `foehncraft_kernels::particle`, the
//!   walls, the tallies, the state for a restart and the VTKHDF output.
//!
//! - [`coupling`]: the records, the sorted per-cell sums and the books of
//!   the two-way coupling (SP2, D247 to D252).
//!
//! One-way coupling by default; `[spray.coupling] mode = "two-way"` makes
//! the drag a source of the gas momentum (SP2). One process: parcel
//! migration between ranks is SP3.

pub mod batch;
pub mod config;
pub mod coupling;
pub mod datum;
pub mod equilibrium;
pub mod sizes;
pub mod spray;

pub use config::{
    BoundaryMode, Breakup, BreakupModel, Coupling, CouplingMode, Dispersion, DragLaw, Ecn,
    Equilibrium, Evaporation, EvaporationModel, Film, Injector, Liquid, Passes, PassesName,
    Penetration, PengRobinson, Saturation, SprayBoundary, SpraySection, Steady,
};
pub use coupling::{CellSources, Delta, Record, StepBooks, limiter_factor};
pub use datum::{Bound, Datum, DatumTally, HistSpec, Moments, Shape, Stats, Surface, Volume};
pub use equilibrium::{SurfacePoint, SurfaceTable, TableEnd, fit_mathias_copeman};
pub use sizes::SizeDistribution;
pub use spray::{
    BoundaryRule, BreakupState, Eddy, Env, Gas, GasState, InjectorStats, MomentumFlows, Parcels,
    Pass, Spray, StreamBooks, StreamPass, Streams, Tally, Vapour, injector_stats, moments_of,
};
