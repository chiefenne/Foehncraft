//! Lagrangian parcels of the sprays (L12, version 1, D240 to D246).
//!
//! - [`batch`]: the parcel step for a batch of parcels, the CPU backends of the
//!   kernel that also runs on the GPU (D245);
//! - [`config`]: the `[spray]` table and its `include` file;
//! - [`datum`]: the datums of the sprays, the crossing records and the
//!   droplet statistics from moment sums (SP1b, D262 to D265);
//! - [`sizes`]: the droplet size distributions of the injectors (D242);
//! - [`spray`]: the parcels in structure-of-arrays form, the injectors,
//!   the tracking of each parcel through the polyhedral mesh (face walk of
//!   Macpherson, Nordin and Weller 2009, `foehncraft_mesh::locate`) with
//!   the exact drag relaxation of `foehncraft_kernels::particle`, the
//!   walls, the tallies, the state for a restart and the VTKHDF output.
//!
//! One-way coupling in this version: the parcels feel the gas, the gas does
//! not feel them (two-way coupling is SP2). One process: parcel migration
//! between ranks is SP3.

pub mod batch;
pub mod config;
pub mod datum;
pub mod sizes;
pub mod spray;

pub use config::{BoundaryMode, Injector, Liquid, Penetration, SprayBoundary, SpraySection};
pub use datum::{Bound, Datum, DatumTally, HistSpec, Moments, Shape, Stats, Surface, Volume};
pub use sizes::SizeDistribution;
pub use spray::{
    BoundaryRule, Env, Gas, GasState, InjectorStats, Parcels, Spray, Tally, injector_stats,
    moments_of,
};
