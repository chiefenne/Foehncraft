//! Machine-readable diagnostics (D78): one registry of stable codes for
//! the validation of the case file (D39), the boundary-condition check
//! (D30), the run (stops, limits, re-reads, checks) and, with the health
//! checks, divergence (D63) and the conservation ledger (D88).
//!
//! A [`Diagnostic`] has a code from [`CODES`], a severity, the key path or
//! the boundary it concerns, the message (the same text the log shows) and,
//! where one exists, a fix as `key=value` changes that `--set` (D62) can
//! apply. A code is never reused for another meaning; a retired code stays
//! in the registry, marked as such. `docs/interfaces.md` lists the codes
//! (a test keeps the two in step).
//!
//! Errors of the case reader and the validation carry their diagnostic
//! inside the `anyhow::Error` ([`err`]), so a front end finds it with
//! [`of_error`] and the human message stays as it was.

use serde::Serialize;
use std::fmt;

/// Severity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, schemars::JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    /// Information.
    Note,
    /// The run goes on; the result may not be what was meant.
    Warning,
    /// The case cannot be run, or the run failed.
    Error,
}

impl Severity {
    /// Lower-case name.
    pub fn name(self) -> &'static str {
        match self {
            Self::Note => "note",
            Self::Warning => "warning",
            Self::Error => "error",
        }
    }
}

/// One entry of the registry.
#[derive(Debug, Clone, Copy, Serialize, schemars::JsonSchema)]
pub struct CodeInfo {
    /// The code.
    pub code: &'static str,
    /// Usual severity.
    pub severity: Severity,
    /// What it means.
    pub title: &'static str,
}

macro_rules! codes {
    ($($name:ident = $code:literal, $sev:ident, $title:literal;)*) => {
        $(
            #[doc = $title]
            pub const $name: &str = $code;
        )*
        /// Every code, in order.
        pub const CODES: &[CodeInfo] = &[
            $(CodeInfo { code: $code, severity: Severity::$sev, title: $title },)*
        ];
    };
}

codes! {
    // the case file (D39)
    CASE_SYNTAX = "CASE001", Error, "The case file is not valid TOML.";
    CASE_UNKNOWN_KEY = "CASE002", Error, "A key that the case file does not have (with the nearest key as a hint).";
    CASE_BAD_VALUE = "CASE003", Error, "A value of the wrong type or an unknown choice.";
    SCHED_ENTRY = "SCHED001", Error, "An entry of `[[schedule]]` that is not valid: its trigger, its order, a key that cannot change during a run, or a value that does not validate (D128).";
    CASE_RANGE = "CASE004", Error, "A value outside its range.";
    CASE_INCONSISTENT = "CASE005", Error, "Settings that do not fit together.";
    CASE_FILE = "CASE006", Error, "A file that cannot be read: the case, the mesh or the restart file.";
    CASE_OUTSIDE_ROOT = "CASE007", Error, "A path outside the read roots or the write roots (D82).";
    CASE_SET = "CASE008", Error, "A `--set key=value` that cannot be applied.";
    CASE_MISSING_KEY = "CASE009", Error, "A key that the case needs is missing.";
    CASE_REFUSED = "CASE010", Error, "A case file that its reader refuses for a reason without a more specific code yet (the message says which); D202.";
    // the boundary conditions (D30)
    BC_NO_CONDITION = "BC001", Error, "A mesh boundary without a boundary condition.";
    BC_NO_BOUNDARY = "BC002", Error, "A boundary condition for a boundary the mesh does not have.";
    BC_NOT_AVAILABLE = "BC003", Error, "A condition that a later stage brings.";
    BC_INLET_INPUT = "BC004", Error, "A velocity inlet with none or more than one of `velocity`, `normal_velocity` and `profile`.";
    BC_PROFILE_ANALYTIC = "BC005", Error, "`profile = \"analytic\"` without an `[analytic]` table.";
    BC_PROFILE_PARABOLIC = "BC006", Error, "`profile = \"parabolic\"` without `max_velocity` or `profile_axes`.";
    BC_PROFILE_UNUSED = "BC007", Warning, "`max_velocity` or `profile_axes` without a parabolic profile.";
    BC_NEGATIVE_MASS_FLOW = "BC008", Warning, "A mass-flow inlet with a negative mass flow.";
    BC_TURBULENCE_UNUSED = "BC009", Warning, "Turbulence values in a laminar case.";
    BC_WALL_NORMAL_VELOCITY = "BC010", Warning, "A wall velocity with a component normal to the wall.";
    BC_EXACT_OUTFLOW = "BC011", Note, "The exact solution leaves the domain through a velocity inlet.";
    BC_INLET_OUTFLOW = "BC012", Warning, "A velocity inlet whose velocity points out of the domain.";
    BC_SYMMETRY_NOT_PLANAR = "BC013", Warning, "A symmetry boundary that is not planar.";
    BC_UNBALANCED = "BC014", Error, "A net mass flow through the prescribed boundaries and no pressure boundary to balance it.";
    BC_REFERENCE_CELL = "BC015", Note, "No pressure boundary: the pressure level is fixed in the reference cell.";
    BC_OVER_SPECIFIED = "BC016", Error, "A flow path whose flow rate is over-specified.";
    BC_FLUID_AT_REST = "BC017", Warning, "Nothing drives a flow: the result is the fluid at rest.";
    BC_PRESSURE_DRIVEN = "BC018", Note, "No prescribed inflow: the flow rate follows from the pressure boundaries.";
    BC_TABLE = "BC019", Error, "A table inlet profile that is missing, unreadable or not a valid table.";
    BC_TABLE_RANGE = "BC020", Warning, "Faces of a table inlet outside the range of the table (they take the end values).";
    BC_CONVECTIVE_STEADY = "BC021", Error, "A convective outlet in a steady case.";
    BC_CONVECTIVE_VELOCITY = "BC022", Error, "A convective outlet without `convection_velocity` and no prescribed inflow to take it from.";
    BC_WEDGE = "BC023", Error, "Wedge boundaries that do not form one wedge about an axis.";
    BC_WEDGE_ANGLE = "BC024", Warning, "A wedge wider than 10° (the axisymmetric error grows with the angle).";
    BC_POTENTIAL_CONVECTIVE = "BC025", Error, "A potential-flow start with a convective outlet and no pressure boundary.";
    BC_WALL_MASS_CLOSED = "BC026", Error, "Walls that add or remove mass in a closed volume (no pressure boundary, no uniform thermodynamic pressure).";
    BC_ROUGHNESS = "BC027", Warning, "A wall roughness k_s larger than the first-cell distance on more than 10 % of the wall's area.";
    BC_TURBULENCE_UNUSED_VALUE = "BC028", Warning, "A turbulence value that the run does not use (on a boundary other than an inlet or a pressure outlet, whose backflow takes the values; `zeta` with a model other than k-ε-ζ-f, ζ-f₀ or the elliptic-blending k-ε); `omega` with an ε model and `epsilon` with SST are converted (ε = C_μ k ω); BC009 stays for laminar cases.";
    TURB_INLET_CONFLICT = "BC029", Error, "Two incompatible turbulence specifications on one inlet (for example `turbulence_length_scale` and `viscosity_ratio`, `k` with intensity, or `omega` with `epsilon`).";
    BC_PERIODIC = "BC030", Error, "A periodic boundary without a periodic partner that names it back, a pair whose faces do not match, or a driving given twice.";
    BC_PERIODIC_UNDRIVEN = "BC031", Warning, "A periodic pair with no driving (`mass_flow`, `pressure_gradient`) and no body force: the flow decays to rest.";
    // passive scalars (D51, D116)
    SCALAR_BOUNDARY = "SCAL001", Error, "A `[boundary.<name>.scalar.<scalar>]` table for an unknown scalar or boundary, without exactly one of `value` and `flux`, or with a condition its boundary does not take (a flux at an inlet, anything at an outlet, `periodic` away from a periodic pair, the heat-flux mode with fixed values).";
    // turbulence models (stage 6, D137)
    TURB_Y_PLUS = "TURB001", Warning, "First-cell y⁺ outside the claimed range of the wall treatment on more than 10 % of a wall's area: above 2 with `low-re`, above 300 (above the log layer) with `automatic`, below 30 with `scalable`, outside 30 to 300 with `standard` and `compound` of the k-ε family (per boundary with the area-weighted min, mean and max).";
    TURB_FLOOR = "TURB002", Warning, "A positivity floor of k, ω or ε acted on more than 0.1 % of the cells in an iteration of the final tenth of the solves (equation and cell count); the run record counts all iterations (`turbulence_floors`).";
    TURB_INLET_RANGE = "TURB003", Warning, "Inlet turbulence outside the usual range: I > 0.5, or a viscosity ratio below 10⁻² or above 10⁵.";
    TURB_KATO_ROTATING = "TURB004", Note, "Kato–Launder production in a rotating zone: Ω is the absolute-frame vorticity.";
    TURB_NO_WALL = "TURB005", Note, "No wall boundary: the wall distance is infinite (F₁ = 0, F₂ = 0; the ζ-f and EB near-wall terms vanish).";
    TURB_UNBOUNDED_SCHEME = "TURB006", Error, "A convection scheme or time scheme that is not bounded for a turbulence variable (central in a steady run, Crank–Nicolson).";
    TURB_CONSTANTS = "TURB007", Warning, "Model constants changed in `[turbulence.constants]`.";
    TURB_CAP = "TURB008", Warning, "The ν_t cap `max_viscosity_ratio` acted (cell count).";
    TURB_INLET_DEFAULT = "TURB009", Note, "Inlet turbulence defaulted to I = 0.05, r = 10 (boundary named).";
    TURB_CURVATURE_KATO = "TURB010", Note, "The curvature correction of SST together with Kato–Launder production: both act on stagnation and swirl regions (the paper uses the standard production).";
    // expressions (D53)
    EXPR_SYNTAX = "EXPR001", Error, "An expression that cannot be parsed.";
    EXPR_NAME = "EXPR002", Error, "An unknown or reserved name in an expression, a constant or a table (with the nearest name as a hint).";
    EXPR_VARIABLE = "EXPR003", Error, "An expression that uses a variable its key does not take (`t` in a steady case or in a property, x, y, z in a value for a whole boundary).";
    EXPR_TABLE = "EXPR004", Error, "A `[tables]` entry that is missing, unreadable or not a valid table.";
    EXPR_VALUE = "EXPR005", Error, "An expression or constant whose value is not finite where it is evaluated.";
    // materials and regions (D94-D102, D112)
    MAT_UNKNOWN = "MAT001", Error, "A material that neither the case nor the database defines (with the nearest name as a hint).";
    MAT_FORM = "MAT002", Error, "A property form that this version does not evaluate, or a property that the flow cannot take yet (a density that varies).";
    MAT_RANGE = "MAT003", Warning, "A property evaluated outside the validity range of its data (clamped to the range).";
    MAT_BASE = "MAT004", Error, "Materials whose `base` refer to each other in a cycle.";
    MAT_VALUE = "MAT005", Error, "A property value that is not positive and finite.";
    MAT_SOLID = "MAT006", Error, "A solid material without density, specific heat or thermal conductivity.";
    MAT_THERMAL = "MAT007", Error, "A fluid without specific heat or thermal conductivity in a run with the energy equation.";
    MAT_AXES = "MAT008", Error, "An orthotropic thermal conductivity without both axes, or with parallel axes.";
    MAT_FILE = "MAT009", Error, "A user material database file (`[materials] include`, `FOEHN_MATERIALS`) or a property data file (NASA-7 record, table file) that cannot be read or has the wrong format.";
    MAT_STEAM = "MAT010", Error, "A water or steam state (`materials.steam`, `foehn steam`) outside IAPWS-IF97, or asked with a set of inputs other than (p, T), (p, h), (p, s), (p, x), (T, x), p or T alone.";
    // sprays (L12, D240 to D246)
    SPRAY_CONFIG = "SPR001", Error, "A `[spray]` value that is out of range, an injector outside the mesh, a spray file (`[spray] include`) that cannot be read or defines `injector` or `liquid` twice, or a `[spray.boundary.<name>]` for a boundary the mesh does not have.";
    SPRAY_STEADY = "SPR002", Error, "The steady spray mode (D260): `[spray.steady]` in a case with `[time]`, or an injector with `start` or `end` in a steady case (its mass flow must not change in time).";
    SPRAY_UNSUPPORTED = "SPR003", Error, "A `[spray]` in a case with MRF zones, sliding interfaces or a periodic boundary, or two-way coupling with a setting it does not work with yet (solid regions, the frozen-flow march, an algorithm other than PISO), or the steady spray mode with solid regions or, two-way, with the coupled algorithm (not in this version).";
    SPRAY_LOST = "SPR004", Warning, "Parcels lost: the walk through the mesh failed or a parcel needed more than `[spray] max_substeps` sub-steps in a time step (counted in `spray.csv`).";
    SPRAY_BOOK = "SPR005", Error, "The mass book of the sprays does not close: released liquid differs from the liquid in the flow, stuck, left and lost by more than round-off.";
    SPRAY_COUPLING = "SPR006", Warning, "Two-way spray coupling (SP2, D252): the parcel passes of a time step ended at `[spray.coupling] max_passes` with the coupling lag above `tolerance` (fix: more passes or a smaller time step; `spray.csv` counts such steps).";
    SPRAY_UNFINISHED = "SPR007", Warning, "Steady spray mode (D260): streams still in the flow after `[spray.steady] max_time` (caught in a recirculation, or a domain the droplets need longer to cross); their mass flow is counted as unfinished in `spray.csv` and `result.json` (fix: a larger `max_time`).";
    SPRAY_SOURCES = "SPR008", Warning, "Steady spray mode (D260 (4)): the source residual of the final stream pass is above `[spray.steady] tolerance`, so the run is not reported as converged (fix: more outer iterations, or a smaller `[spray.steady] relaxation` if the residual swings).";
    SPRAY_LIMITER = "SPR009", Warning, "The saturation limiter of D253 held vapour back: the evaporation of a spray cell was above what brings the gas the cell holds and the gas that flows in to adiabatic saturation (transient: in a step; steady: in the final stream pass). Counted in `spray.csv` (`limiter_cells`, `limiter_kg_s`) and `result.json` (fix: more passes or a smaller time step; steady, more outer iterations or a smaller `[spray.steady] relaxation`).";
    SPRAY_DEVICE = "SPR010", Warning, "An injector names a device hole (D342) and gives its own `position` and `direction` as well, and they lie more than 5 % of the exit radius plus 10 µm from the hole exit (plus `offset`) or more than 0.5° from the hole axis. The injector's own values are used; check the device's `mount`, `axis`, `reference` and holes, or leave `position` and `direction` out.";
    // two fluids by the volume-of-fluid method (L3, stage MP1)
    VOF_CONFIG = "VOF001", Error, "A `[vof]` value that does not fit: not exactly two `[[vof.phase]]` tables, a density or viscosity that is not > 0, a negative surface tension, `thinc_beta` outside (0, 20], `max_courant` outside (0, 1], more than 10 smoothing passes, an α outside [0, 1] (`initial`, a region, an inflow), a region without the keys of its shape, or an inflow for a boundary the case does not have.";
    VOF_UNSUPPORTED = "VOF002", Error, "A `[vof]` case with a setting the first version (MP1) does not have: a steady case, a turbulence model, the energy equation, passive scalars or species, a `[buoyancy]` model, regions or several fluids, sprays, immersed bodies, the coupled algorithm, rotating zones or sliding interfaces, a moving mesh, periodic boundaries, ranks, or a restart.";
    VOF_CLIPPED = "VOF003", Warning, "The volume of phase 1 clipped to keep 0 ≤ α ≤ 1 is above 10⁻⁶ of the phase volume over the run (fix: a smaller `[vof] max_courant` or converged pressure corrections; `vof.csv` counts it per step).";
    // immersed bodies, the surface check (L13 stage 0, D225)
    IMM_SURFACE_EMPTY = "IMM001", Error, "A triangle surface (STL) with no triangles or with corners whose coordinates are not finite.";
    IMM_SURFACE_OPEN = "IMM002", Warning, "An open triangle surface: edges with only one triangle (holes, open borders). Inside and outside still come from the generalised winding number, which closes small holes smoothly; nothing is repaired.";
    IMM_SURFACE_NON_MANIFOLD = "IMM003", Warning, "A triangle surface with edges shared by more than two triangles (non-manifold: touching bodies, duplicate or internal triangles).";
    IMM_SURFACE_ORIENTATION = "IMM004", Warning, "A triangle surface whose orientation is not consistent: edges whose two triangles run them in the same direction (some triangles flipped).";
    IMM_SURFACE_INTERSECTING = "IMM005", Warning, "A self-intersecting triangle surface: pairs of triangles that cross or touch away from their shared corners and edges (exact orientation predicates, D225). The winding number still decides inside and outside.";
    IMM_SURFACE_DEGENERATE = "IMM006", Warning, "Degenerate triangles: corners that are repeated or exactly collinear (zero area).";
    IMM_SURFACE_DUPLICATE = "IMM007", Warning, "Duplicate triangles: two triangles with the same three corners.";
    IMM_SURFACE_INVERTED = "IMM008", Warning, "A closed surface whose normals point inward (negative enclosed volume): as given, the body is its complement (fix: flip the triangles in the CAD export).";
    IMM_SURFACE_UNITS = "IMM009", Note, "A surface hint: a bounding box larger than 20 m (a file in millimetres needs `scale = 1e-3`), or distinct corners closer than 1e-9 of the box diagonal (welding candidates, not welded).";
    // immersed bodies, Brinkman penalisation (L13 stage 1, D226)
    IMM_CONFIG = "IMM101", Error, "An `[immersed]` value that does not fit: a permeability that is not > 0, a body without the keys of its shape (`centre` and `radius`, `lower` and `upper`), a radius or length that is not > 0, a zero cylinder axis, two bodies with the same name, a `field` body without `fraction`, an `stl` body without `file` or with a `scale` that is not > 0, `file` or `scale` on a body of another shape, no body at all, or an `[immersed.topology]` value out of range (an empty design region, a fluid fraction not in (0, 1], a negative filter radius).";
    IMM_UNSUPPORTED = "IMM102", Error, "An `[immersed]` body in a case it does not fit (L13 stage 1, D226): a turbulence model (penalised walls are laminar and low-Re only; a wall function on a penalised wall is refused), the energy equation (no penalised conjugate heat transfer, Q9), passive scalars, gravity, regions or several fluids, rotating zones, or sprays; a topology search in a transient run.";
    IMM_RESOLUTION = "IMM103", Warning, "An immersed body that the mesh does not resolve: no cell lies more than half inside it, or its solid volume on the mesh differs from its exact volume inside the mesh by more than 5 %.";
    IMM_BODIES = "IMM104", Note, "The immersed bodies of the case: their solid volume on the mesh against the exact volume, the cells they touch and the sink coefficient μ/K; the design cells of a topology search.";
    // immersed bodies, cut cells (L13 stage 2, D227 to D229)
    IMM_CUT = "IMM201", Error, "The cut of `[immersed] method = \"cut\"` failed: the bodies fill the whole mesh, a body has the name of a mesh boundary, the mesh is a rank's part or periodic before the cut, or an `[immersed.cut]` value is out of range (`min_volume_fraction` not in [0, 0.9]).";
    IMM_CUT_UNSUPPORTED = "IMM202", Error, "Cut cells (L13 stage 2) with a setting they do not have yet: a cut body kept as a solid region without the energy equation, a topology search (it needs the penalised density field), rotating zones or a moving mesh (rigid motion is stage 5).";
    IMM_CUT_REPORT = "IMM203", Note, "The cut of the immersed bodies: cells cut, removed, split into pieces and merged, the smallest kept volume fraction, the wall area of every body against its exact area where simple (an STL body: the area of its triangles, when the surface is clean and lies inside the mesh), the closure of the cut cells, the cells cut with fewer wall and face points so that the discretisation takes their faces (the fallback), and the merged cells that share several faces with a neighbour (faces not in one plane are kept apart).";
    IMM_CUT_QUALITY = "IMM204", Warning, "Cut cells the merge could not fix: a small cut cell without a neighbour to merge into, or faces whose normal does not point from the cell centroid to its neighbour (the run refuses them). Fix: refine the mesh near the body, or move the body by a small fraction of a cell.";
    IMM_THIN_SOLID = "IMM205", Error, "A cut body kept as a solid region (L13 stage 2d, D229) is thinner than about two cells: its deepest solid cell lies less than one cell width inside the surface, or small solid cells found no solid neighbour to merge into. Thin solids are not resolved on cut cells (owner decision Q4). Fix: refine the mesh in the body, or model the wall as a thin wall on a region interface (`[interfaces.<name>] thin_wall`, D189).";
    // datums and their reports (SP1b, D262 to D267)
    DATUM_CONFIG = "DAT001", Error, "A `[[datum]]` table that cannot be used: a name used twice or with other characters than letters, digits, `_` and `-`, an unknown type, a missing or invalid key of its type (plane point and normal, bound, boundary, zones, STL file, box corners, cylinder, zone), or bad histogram bins.";
    DATUM_EMPTY = "DAT002", Error, "A `[[datum]]` that meets no face or cell of the mesh (a plane outside the domain, an interface between zones that do not touch, an empty volume).";
    ENE_UNSUPPORTED = "ENE006", Error, "An energy setting this version does not have: a passive scalar named `temperature` beside `[energy]`, solid regions with a sliding mesh, a restart with `[energy.acceleration]`, passive scalars or a periodic mesh, or with the variable density (stage 7c) passive scalars, a restart, a periodic mesh, a convective outlet, a wall mass flux or a rotating zone.";
    ENE_SOURCE = "ENE008", Error, "An `[[energy.sources]]` entry for a region the case does not have, or without exactly one of `power` and `power_density`.";
    CHT_INTERFACE = "CHT004", Error, "An `[interfaces.<name>]` table that is not on a region interface (regions unknown, the same, two fluids, or not touching), or without exactly one of `contact_resistance` (≥ 0) and `thin_wall` (thickness > 0 with `conductivity` or `material`).";
    ACC_SETTING = "ACC001", Error, "An `[energy.acceleration]` setting that cannot be used: a steady case, no solid region, a region that is not a solid, `inertia_scale` outside (0, 1] and not `\"auto\"`, `restore` neither `\"settled\"` nor a time ≥ 0, or a `period`, `window_steps`, `settle_tolerance` or `target_windows` out of range.";
    ACC_UNRESTORED = "ACC002", Warning, "The run ended before the true heat capacity of the accelerated solids was restored: the final state is that of the accelerated solids and the history is not physical.";
    ACC_SKIN = "ACC003", Note, "The accelerated solids' skin depth of the period, √(2k/(f ρ c_p ω)), is larger than the thickness V/A of a region (V its volume, A its interface area): while accelerated the swing of the period reaches through the wall; it goes away with the restore.";
    ACC_PHASE = "ACC004", Note, "A step of the thermal acceleration: the conductances and factors measured, the solids settled, the restore started or finished.";
    CYC_SETTING = "CYC001", Error, "An `[energy.cycle]` setting that cannot be used: a steady case, no solid region, together with `[energy.acceleration]`, an adaptive time step, a `period` that is not a whole number (≥ 2) of time steps, or a count, tolerance or relaxation out of range.";
    CYC_NOT_PERIODIC = "CYC002", Warning, "A fluid phase of the cycle-averaged loop ended at `max_cycles` while the mean interface heat flow still changed between the last two cycles by more than `cycle_tolerance`: the averages are not those of a periodic flow.";
    CYC_SKIN = "CYC003", Warning, "`skin = \"effusivity\"` with a skin depth of the period √(2k/(ρ c_p ω)) that is not small against the thickness V/A of a solid region (larger than a third of it): the swing reaches through the wall there and the semi-infinite solution does not hold; `resolved_cycles` gives the exact swing.";
    CYC_MAX_UPDATES = "CYC004", Warning, "The cycle-averaged loop stopped at `max_updates` before the wall temperatures and the interface heat flows settled within `tolerance` and `heat_tolerance`.";
    CYC_PHASE = "CYC005", Note, "A step of the cycle-averaged loop: a fluid phase averaged, a solid phase solved, the wall temperatures updated, the loop converged.";
    FRZ_SETTING = "FRZ001", Error, "An `[energy.frozen_flow]` setting that cannot be used: a steady case, no solid region, together with `[energy.acceleration]` or `[energy.cycle]`, the variable density, an adaptive time step, a `time_step` below the `[time]` step, or a `start`, window, tolerance or count out of range.";
    FRZ_NOT_SETTLED = "FRZ002", Warning, "The flow of a frozen-flow run did not settle by the window test (an unsteady flow, or a tolerance too tight): the run stayed fully coupled with the `[time]` step from the last freeze or refresh on.";
    FRZ_PHASE = "FRZ003", Note, "A step of the frozen-flow march: the flow settled and was held, or it was refreshed by coupled steps because the temperatures it depends on had moved.";
    FRZ_TIME = "FRZ004", Warning, "A frozen-flow run whose flow boundary values or body force depend on time: the held flow does not follow them between refreshes (`refresh_every`).";
    ENE_BALANCE = "ENE007", Warning, "A steady run that converged while its heat balance does not close (imbalance above 1 % of the largest heat flow): the temperature residual is small but the energy solve is not converged; tighten `[solver.energy] tolerance` and raise `max_iter`.";
    THM_LEVEL = "THM001", Error, "A steady run whose boundaries leave the temperature level undetermined: no inlet and no `temperature` or convection boundary, and no net heat input.";
    THM_DRIFT = "THM002", Error, "A steady run with a net heat input (heat flux, heat flow, source or viscous heating) and no way for the heat to leave: no inlet and no `temperature` or convection boundary.";
    THM_SIDE = "THM003", Error, "A thermal condition on a boundary that cannot carry it: a flux or convection at an inlet, any condition at an outlet, a periodic boundary, a symmetry plane or a wedge.";
    THM_RANGE = "THM004", Error, "A thermal value out of range: a temperature or ambient temperature ≤ 0 K, a heat-transfer coefficient ≤ 0, a heat flow that is not finite.";
    THM_ADIABATIC = "THM005", Warning, "Walls without a `[boundary.<name>.energy]` table in a run with the energy equation: they are adiabatic (listed once).";
    BUO_NO_GRAVITY = "BUO001", Error, "`[gravity]` with a vector of zero length, or a buoyancy model without `[gravity]`.";
    BUO_RANGE = "BUO002", Warning, "The Boussinesq model outside its range: β ΔT above 0.1 over the temperatures of the case (the density changes by more than 10 %).";
    BUO_NO_EXPANSION = "BUO003", Error, "The Boussinesq model without what it needs: no energy equation, no expansion coefficient (the fluid's density does not depend on T and `[buoyancy] expansion_coefficient` is not given), or a reference value out of range.";
    BUO_NONE = "BUO004", Warning, "`[gravity]` with the energy equation but no buoyancy model (the temperature drives no flow); as information: a turbulence model with buoyancy, whose production of k by buoyancy comes in stage 7e1.";
    BUO_CLOSED = "BUO005", Error, "A variable density in a closed domain (no inlet, no outlet) without `[energy] thermodynamic_pressure = \"closed-volume\"`: the mass in the domain would not be conserved.";
    VD_SETTINGS = "ENE009", Error, "Variable-density settings that do not fit together: a density that does not depend on T, `density = \"variable\"` with the Boussinesq model, a closed volume with an inlet or outlet, a closed volume of a liquid, or `closed-volume` without the variable density.";
    VD_MACH = "ENE004", Warning, "A Mach number above 0.3 in a run with the variable density: the low-Mach form loses its validity (the density does not follow the dynamic pressure).";
    BUO_TALL_OPENING = "BUO006", Warning, "A pressure boundary with `pressure_kind = \"static\"` over a height span whose hydrostatic difference ρ_ref g Δz is above the dynamic pressure: a uniform static value drives a flow through the opening.";
    THM_NO_ENERGY = "THM009", Error, "A `[boundary.<name>.energy]` table without `[energy]`: the energy equation is solved only with an `[energy]` table.";
    ENE_BOUNDARY = "THM010", Error, "A `[boundary.<name>.energy]` table for an unknown boundary, without exactly one condition (`temperature`, `heat_flux`, `heat_flow`, or `heat_transfer_coefficient` with `ambient_temperature`), or with `heat_transfer_coefficient` and `ambient_temperature` apart.";
    ENE_INLET = "THM011", Error, "An inlet without `[boundary.<name>.energy] temperature` in a run with the energy equation.";
    REG_NOT_FLUID = "REG001", Error, "A region that is not a fluid where it cannot be: a solid region without `[energy]` (solids are solved by the energy equation, stage 7d1), or a solid material as the `[fluid]`.";
    REG_ZONE = "REG002", Error, "A region that names a cell zone the mesh does not have, or a zone in two regions.";
    REG_PART_REFERENCE = "REG003", Error, "A fluid part (fluid cells connected without a solid between them) without a pressure boundary in a case of several fluid parts: give it a pressure outlet or inlet; a closed part of its own comes with stage MF-B.";
    REG_FLUIDS_TOUCH = "REG004", Error, "Two fluids of different materials (or temperatures) share a face: that is a two-phase interface (multiphase, L3); put a solid region or a thin wall between them, or give both regions the same fluid.";
    REG_PARTS = "REG005", Note, "The fluid parts of the case: their number, regions and fluids.";
    REG_FLUIDS_UNSUPPORTED = "REG006", Error, "A setting that does not yet work with several fluids in one case (D194): the variable density, gravity, a rotating zone or frame, passive scalars, a convective outlet or a periodic mass-flow drive.";
    REG_INITIAL_KEY = "REG007", Error, "A key of `[regions.<name>.initial]` that does not fit: a solid region takes `temperature` only, turbulence start values come as `k` with one of `omega` or `epsilon`, and a start temperature needs `[energy]`.";
    REG_INITIAL_VALUE = "REG008", Error, "A start temperature, k, omega or epsilon of `[regions.<name>.initial]` that is not > 0.";
    POROUS = "POR001", Error, "A porous region without its matrix (`[regions.<name>.porous]` with `viscous`/`inertial` or `honeycomb`), contradicting or invalid values, direction-dependent values without `axes`, a matrix in a region that is not porous, or with `[energy]` a matrix without `solid`, porosity, h or a.";
    POROUS_UNSUPPORTED = "POR002", Error, "A setting that does not yet work with porous regions (porous plan, D279): with the energy equation a cycle-averaged or accelerated energy mode or a density that varies, passive scalars, a turbulence model, a rotating zone, the GPU (stage PM5), or direction-dependent resistances whose axes are turned against the mesh axes.";
    SPECIES = "SPE001", Error, "Invalid `[species]` settings or compositions: an unknown species (the GRI-Mech 3.0 set), a species named twice, a closure species not in `names`, fractions that do not sum to 1 or leave the closure species negative, both or neither of `mole_fractions` and `mass_fractions`, an inlet without `[boundary.<name>.species]`, `properties = \"mixture\"` or `diffusion = \"lewis\"` without `[energy]`, or a variable density without an ideal-gas fluid with its molar mass.";
    SPECIES_UNSUPPORTED = "SPE002", Error, "A setting that does not yet work with species (stage SF1): regions (solids, several fluids, porous regions); conjugate walls come with SF5.";
    SPECIES_CLOSURE = "SPE003", Warning, "The transported mass fractions summed to more than 1 somewhere by more than 1e-4, so the closure species was held at zero there: the solve is not converged, the tolerance of `[solver.scalar]` is too loose, or the limiter lets extrema pass (`limiter_epsilon`).";
    CHEMISTRY = "CHM001", Error, "An invalid chemistry mechanism or reactor input: an unknown built-in mechanism, a mechanism file that does not parse, a reaction that does not conserve the elements, a species without thermodynamic data, an unsupported rate form (pressure-dependent Arrhenius, Chebyshev, SRI fall-off, a specific collider) or unit, or a reactor input with an invalid composition, temperature, pressure or time.";
    CHEMISTRY_INTEGRATOR = "CHM002", Warning, "The stiff chemistry integrator of a reactor did not reach the end of its interval within `max_steps`, or its step became too small: raise `max_steps`, loosen `rtol`, or check the mechanism and the state.";
    CHEMISTRY_UNSUPPORTED = "CHM003", Error, "A setting that does not work with chemistry in the flow: a steady laminar case, a turbulence model without a turbulence–chemistry interaction (`interaction = \"edc\"` or `\"edm\"`, SF3), an interaction without a turbulence model, sprays, or ranks.";
    CHEMISTRY_FINE_STRUCTURES = "CHM004", Warning, "The fine-structure reactor of the eddy dissipation concept (SF3) failed in some cells (too many steps, the step became too small, or a state out of range); those cells kept the rates of their last renewal (steady) or got no reaction source in that step (transient): raise `[chemistry.integrator] max_steps`, loosen `rtol`, or check the state of the cells.";
    REG_INITIAL_INLET = "REG009", Error, "`temperature = \"inlet\"` in `[regions.<name>.initial]` for a region whose fluid part has no inlet with `[boundary.<name>.energy] temperature`.";
    REG_START = "REG010", Note, "The start values of the regions and fluid parts and where they come from: given, a boundary value, a Laplace solve between the pressure boundaries of a part, or the inlets.";
    // pressure values absolute or gauge (D239)
    PRS_NOT_POSITIVE = "PRS001", Error, "An absolute pressure value (`[fluid] pressure_values = \"absolute\"`, the default from case format 0.7, D239) that is zero or negative: a boundary `pressure` or `total_pressure`, `[initial] pressure`, `[pressure] reference_value` or an operating pressure; set `pressure_values = \"gauge\"` for values relative to the operating pressure.";
    PRS_LOOKS_GAUGE = "PRS002", Warning, "An absolute pressure below 1 kPa in a gas (density below 50 kg/m³): it looks like a gauge value; set `[fluid] pressure_values = \"gauge\"` if it is one.";
    PRS_NO_OPERATING = "PRS003", Error, "A fluid part whose operating pressure cannot be taken from its pressure boundaries (absolute values, D239): their values are expressions only; give `[fluid] operating_pressure` or `[regions.<name>] operating_pressure`.";
    PRS_INITIAL_PARTS = "PRS004", Error, "A case-wide absolute `[initial] pressure` in a case whose fluid parts have different operating pressures: leave it out (each part starts at its own operating pressure) or give the start pressure per region.";
    PRS_FRAME = "PRS005", Note, "The pressure convention of the case file (absolute or gauge) and the operating pressure of every fluid part with where it comes from (D239).";
    // moving reference frame and sliding mesh (D55, D182)
    MRF_SLIDING_ZONE = "MRF001", Error, "An `[[mrf]]` table and a `[[sliding]]` table that turn the same cells (one of them without a `zone`, or both with the same one): a zone is turned by a moving reference frame (steady or transient, D185) or by a sliding mesh (D182), not by both.";
    MRF_ENTRY = "MRF002", Error, "An `[[mrf]]` table with no or two angular velocities (`omega`, `rpm`), a zero axis, a cell zone the mesh does not have (with the nearest name as a hint), no `zone` where there are several tables, or a zone that two tables rotate.";
    MRF_INTERFACE = "MRF003", Warning, "The interface or a boundary of a rotating zone is not a surface of revolution about its axis: the relative flux of the cells next to it does not close, and the solution there is not that of a frozen rotor (cell count and worst closure error).";
    SLIDING_STEADY = "SLD001", Error, "A `[[sliding]]` table in a steady case: the sliding mesh moves in time (D182); a steady case with a rotating zone uses `[[mrf]]` (D55).";
    SLIDING_ENTRY = "SLD002", Error, "A `[[sliding]]` table with no or two angular velocities (`omega`, `rpm`), a zero speed or axis, a `shift` of 0 or of a whole turn, a cell zone the mesh does not have (with the nearest name as a hint), or a second table (one rotor per case so far); with the interface imprinted: `rotor_interface` without `stator_interface`, no `zone`, or a `shift` (the time step is free).";
    SLIDING_STEPPING = "SLD003", Error, "The time stepping does not fit the sliding mesh: the time step is not `shift` × 2π/(sectors × Ω) (the right value is given), the step is adaptive or error controlled, `end_time` is not a whole number of steps, or `[verification]` asks for a sweep over levels or time steps; with the interface imprinted the scheme must be fully implicit (`extrapolation = false`, at least 2 `outer_iterations`).";
    SLIDING_ZONE = "SLD004", Error, "The rotor cannot be moved by whole sectors: its cells have no image when it is turned by one sector (it is not a ring of identical sectors about the axis), or a boundary of the rotor is not a wall, slip wall, symmetry plane or wedge.";
    SLIDING_INTERFACE = "SLD005", Error, "The interface of a sliding mesh by imprinting cannot be imprinted: a boundary or cell zone of `[[sliding]]` is missing or wrong, the two interface boundaries are not on one circular cylinder about the axis or do not cover the same surface, rotor and stator cells share a face or a node, a `[boundary]` table names an interface boundary, or a restart is asked for (not yet available).";
    SLIDING_CLOSURE = "SLD006", Warning, "The cells next to the imprinted interface of a sliding mesh do not close (|Σ S|/Σ |S| above 1e-8): the faces of the interface are not those of a mesh extruded along the axis, and the mesh flux of the rigid turn leaves a small mass source.";
    // the deforming mesh (D91, stage MS0)
    MOTION_CONFIG = "MOV001", Error, "A `[mesh_motion]` table that cannot be used: a steady case (the mesh moves in time), a displacement expression that does not compile, or a case with what the deforming mesh of stage MS0 does not work with yet (turning zones, regions, sprays, immersed bodies, periodic boundaries, the check of `[run] ranks`, an error-controlled time step, a sweep over levels), or a restart file whose nodes the displacement does not put there or whose node, cell or face order is another.";
    MOTION_INVERTED = "MOV002", Error, "The prescribed displacement of `[mesh_motion]` turns a cell inside out (its volume is no longer positive); the message gives the time, the cell and its place.";
    MOTION_BOUNDARY = "MOV003", Error, "The deforming mesh moves a boundary that is not a wall out of its own plane (a face sweeps a volume): only walls move with the mesh; the nodes of symmetry planes, wedges, inlets and outlets may only slide along them.";
    MOTION_CLOSED = "MOV004", Error, "The deforming mesh changes the volume of a closed domain (no pressure boundary) of constant density: an incompressible fluid cannot be compressed; give it a pressure boundary or the variable density of a closed volume.";
    // torque, shaft power and machine performance (D184)
    ROT_SLIP = "ROT001", Warning, "A wall of a rotor does not turn with its zone (a stationary shroud in a rotating zone): its torque is in that of the rotor, but the shaft power −Ω·T is not its work on the fluid; give the walls of the rotor in `walls` of the `[[mrf]]` or `[[sliding]]` table.";
    ROT_ENTRY = "ROT002", Error, "A `walls` entry of `[[mrf]]` or `[[sliding]]` that is not a wall of the case, or an `[[output.machines]]` table without a turning zone, with a `rotor` that is no turning zone (or none where there are several), without an inlet or an outlet, with a boundary the case does not have, or with a name used twice.";
    // the mesh quality (D50)
    MESH_SKEWED = "MESH001", Warning, "Faces whose skewness exceeds `solver.skewness_limit` (poor cells, typically tetrahedra of a mesh that was not optimised).";
    MESH_REMEDIATION = "MESH002", Note, "Cell-quality remediation (D50): the run uses other schemes than set because of poor cells.";
    MESH_OPERATOR_LINEAR = "MESH004", Warning, "The operator check (D236): the least-squares gradient of a linear field is not exact on some cells (a degenerate or badly conditioned stencil); the cells are listed worst first.";
    MESH_OPERATOR_HESSIAN = "MESH005", Warning, "The operator check (D236): with `turbulence.curvature_correction`, the velocity Hessian (DS/Dt) is not consistent for a quadratic field on some cells (error above `solver.operator_tolerance`), so the correction does not converge with the mesh there; the cells are listed worst first.";
    MESH_OPERATOR_GRADIENT = "MESH006", Warning, "The operator check (D236): with `turbulence.qcr`, the velocity gradient is first order on some cells (error above `solver.operator_tolerance` for a quadratic field), so the divergence of the explicit QCR stress does not converge with the mesh there; the cells are listed worst first.";
    MESH_READ = "MESH003", Error, "A mesh file that cannot be read: missing, of an unknown format, truncated or damaged (the reader says where); D202.";
    // the run (D40, D62, D64, D83)
    RUN_STOP_FILE = "RUN001", Note, "Stopped by the `.stop` control file.";
    RUN_SIGNAL = "RUN002", Note, "Stopped by a signal (SIGTERM, SIGINT or SIGUSR2).";
    RUN_LIMIT_WALL_TIME = "RUN003", Warning, "Stopped at the wall-time limit (`[run.limits] max_wall_time`).";
    RUN_LIMIT_STEPS = "RUN004", Warning, "Stopped at the step or iteration limit (`[run.limits] max_steps`, `max_iterations`).";
    RUN_LIMIT_OUTPUT = "RUN005", Warning, "Stopped at the output-volume limit (`[run.limits] max_output`).";
    RUN_NOT_CONVERGED = "RUN006", Error, "A steady solve that did not converge to the tolerance.";
    RUN_CHECK_FAILED = "RUN007", Error, "A verification check of the case failed (order, error, reference or benchmark band, mass balance).";
    RUN_NOT_APPLIED = "RUN008", Warning, "A change of the re-read case file that cannot be applied during a run.";
    RUN_REREAD_INVALID = "RUN009", Warning, "A re-read case file that does not validate; the run goes on with the settings it had.";
    RUN_RESTART_MAPPED = "RUN010", Warning, "A restart file in another cell or face order: the continuation is not bit-identical.";
    RUN_NO_RESTART = "RUN011", Error, "`--restart latest` found no restart file in the output directory.";
    RUN_NO_BENCHMARK = "RUN012", Note, "No benchmark of this hardware: the estimate has no wall time.";
    RUN_IO = "RUN013", Error, "An output file could not be written.";
    RUN_SCHEDULED = "RUN015", Note, "A setting changed by an entry of `[[schedule]]` (D128): key, old and new value.";
    RUN_DIVERGED = "RUN014", Error, "The solution diverged (a value that is not finite); see `divergence.md` (D63).";
    // iteration and study budgets (D198 (c), (d))
    RUN_DIAGNOSIS = "RUN016", Note, "The diagnosis of a steady run at `[run.limits] diagnose_at` (default 2000 iterations): the trend of each residual, the verdict and advice (D198).";
    RUN_DIAGNOSIS_STOP = "RUN017", Warning, "Stopped by the diagnosis: the steady run is not on its way to the tolerance within its iteration budget (`[run.limits] diagnosis = \"stop\"`).";
    RUN_STUDY_BUDGET = "RUN018", Warning, "Stopped at a study budget (`[run.limits] study_cell_hours`, `study_gpu_hours`, `study_output`).";
    RUN_STUDY_SPENT = "RUN019", Error, "A study budget is already spent by the earlier runs of the study (`study_budget.jsonl`): the run does not start.";
    RUN_STUDY_ESTIMATE = "RUN020", Warning, "The estimate of the run (D81) does not fit into what is left of a study budget.";
    RUN_LARGE_BUDGET = "RUN021", Warning, "A steady iteration budget above 10 000 iterations: a steady case needs hundreds, at most a few thousand; look at the run after 2000 to 3000 iterations instead.";
    // the live control port (L19)
    RUN_CONTROL_STOP = "RUN022", Note, "Stopped through the live control port (`run.stop`): restart file, final output, exit code 3 (L19).";
    RUN_CONTROL_SET = "RUN023", Note, "A setting changed through the live control port (`run.set`): key, old and new value, with the equivalent `[[schedule]]` entry (L19).";
    RUN_CONTROL_REJECTED = "RUN024", Warning, "A setting asked for through the live control port that cannot be applied (not valid for the case at that point); the run goes on with the settings it had.";
    RUN_REREAD = "RUN025", Note, "A setting changed by re-reading the case file (`.read_toml_N`, `[run] reread_every`): key, old and new value (D40, D214).";
    // the causes of a divergence (D63)
    DIV_TIME_STEP = "DIV001", Note, "Divergence cause: the Courant number near the failing cells is large; reduce the time step.";
    DIV_RELAXATION = "DIV002", Note, "Divergence cause: under-relaxation factors too large for the case.";
    DIV_CELL_QUALITY = "DIV003", Note, "Divergence cause: poor cells (skewness, non-orthogonality, volume ratio) at the failing cells.";
    DIV_BOUNDARY = "DIV004", Note, "Divergence cause: a boundary condition next to the failing cells.";
    DIV_INITIAL_STATE = "DIV005", Note, "Divergence cause: the run failed within its first steps or iterations, from the initial state.";
    DIV_LINEAR_SOLVER = "DIV006", Note, "Divergence cause: a linear solve that stopped at its iteration limit.";
    // the conservation ledger (D88, D89, D90)
    BAL_MASS_STEP = "BAL001", Warning, "A relative mass imbalance of a step or iteration above `[balance] tolerance_step` (D89).";
    BAL_MASS_RUN = "BAL002", Warning, "A relative mass imbalance over the run above `[balance] tolerance_run` (D89).";
    BAL_NEGATIVE_DISSIPATION = "BAL003", Warning, "A numerical dissipation of kinetic energy below `[balance] tolerance_energy`: the discretisation adds energy (D90).";
    // statistics of transient runs (D71)
    STAT_WINDOW = "STAT001", Warning, "The averaging window of a statistic of `[output.statistics]` is too short: the batch means are correlated (lag-1 autocorrelation above 0.3), so the confidence interval of the mean is too narrow; run longer or use fewer batches (D71).";
    // limitations of a run (D199): never a failure, the exit code stays
    LIM_GRADE = "LIM001", Note, "A feature in use is graded partial or unverified on the verification status page (docs/verification/status.md); the message names its first gap.";
    LIM_RANGE = "LIM002", Note, "A feature in use ran outside its validated range (first-cell y⁺, Reynolds, Prandtl or Courant number, skewness), as `range` of docs/verification/features.toml gives it.";
    LIM_NOT_SETTLED = "LIM003", Note, "A steady solve that stopped before its tolerance (iteration limit): the result may still change with more iterations.";
    LIM_ONE_MESH = "LIM004", Note, "One mesh only: the discretisation uncertainty is not estimated (a level sweep gives the GCI, D67).";
    LIM_DEVICE = "LIM005", Note, "The GPU backend ran on a device that is not in the list of verified devices (Tesla T4 so far).";
    LIM_ACCELERATED = "LIM007", Note, "A transient conjugate run with thermal acceleration of its solids (D213): the history before `physical_from` is not physical (the solids' time ran faster, their swing was too large); the end state and later cycle means are.";
    LIM_OPERATORS = "LIM008", Note, "An operator a feature in use needs is not exact on some cells (MESH005, MESH006, D236): the feature does not converge with the mesh there.";
    LIM_CYCLE = "LIM009", Note, "A cycle-averaged conjugate run (D204): the solids carry the steady temperatures of the cycle-mean heat transfer (with `solid = \"transient\"` a history under it, marched in large steps), and each fluid phase ran with the wall temperatures held: the swing of the wall temperature over a cycle is in the result only as the estimate of `skin` or from `resolved_cycles`.";
    LIM_FROZEN = "LIM010", Note, "A frozen-flow run (D118): the flow was held during the march, so its history is that of the settled flow at the temperatures of the last freeze or refresh; exact for a flow that does not depend on the temperature.";
    LIM_CUT_WALL = "LIM011", Note, "Wall functions on the walls of cut bodies (L13 stage 2c, D229): the first-cell distance is the exact distance of each cut or merged cell to the body, so the first-cell y⁺ varies with the cut fraction from cell to cell along the wall; the message gives the smallest, mean and largest y⁺ of every body (the map is the `y_plus` of the wall output).";
    LIM_REMEDIATION = "LIM006", Note, "Automatic mesh remediation (D50, MESH002) changed the schemes on poor faces: first order there.";
    // the tool registry (D196)
    TOOL_UNKNOWN = "TOOL001", Error, "A tool id that the registry does not have (with the nearest id as a hint); `foehn commands` lists them.";
    TOOL_INPUT = "TOOL002", Error, "A tool input that does not match the input schema of the tool (unknown field, wrong type, missing field).";
    TOOL_JOURNAL = "TOOL003", Warning, "The journal of tool calls (D197) could not be written: the calls run, but the journal misses them.";
    TOOL_SCRIPT = "TOOL004", Error, "A script to replay (`foehn run --script`) that cannot be read: not a file of JSON lines, or without calls.";
    TOOL_REPLAY_DIFFERS = "TOOL005", Warning, "A replayed call (`foehn run --script`) that ended with another status or exit code than in its record; the replay stops there unless `--keep-going`.";
    TOOL_SERIES = "TOOL006", Error, "A series of `result.page` that the run directory does not have, or a cursor that does not belong to it; the call without `series` lists the series.";
    // the MCP server and the approvals of agent calls (L18, D198)
    MCP_ROOTS = "MCP001", Error, "`foehn mcp` started without `--read-root` and `--write-root`, or with an action table (`--policy`) inside a write root, where a tool could rewrite it.";
    MCP_LIMITS = "MCP002", Error, "A run started through MCP without `[run.limits] max_wall_time` and `max_output`; the fix sets them with `set`.";
    MCP_REFUSED = "MCP003", Error, "A call that the action table of the project (`foehn-mcp.toml`) refuses: start a run, exceed the estimate caps, write outside the write roots, overwrite results, publish.";
    MCP_DECLINED = "MCP004", Error, "A call that the person declined when the server asked for approval.";
    TOOL_CONTROL = "TOOL007", Error, "A request to the live control port of a run (`run.status`, `run.stop`, `run.set`) that cannot be made or is refused: no run listens in the directory, the token is wrong, or the setting cannot be changed during a run (L19).";
    MCP_CANNOT_ASK = "MCP005", Error, "A call that needs the person's approval (`ask`), from a host that cannot ask (no MCP elicitation); the person allows it in `foehn-mcp.toml` or uses a host that asks.";
    MCP_POLICY = "MCP006", Error, "The action table `foehn-mcp.toml` cannot be read or written: not TOML, an unknown action kind or setting.";
    // the structural solver (L53, crate foehncraft-structure)
    STR_CASE = "STR001", Error, "A structural case or model that is not valid: an unknown key or a wrong value, a point that is not a node, a beam without length, a missing material or section, non-positive properties.";
    STR_SINGULAR = "STR002", Error, "The stiffness of a structural case is singular: the structure can move freely (a missing support or clamp, or a mechanism); the message names the first free unknown. For the modes of a free structure give a negative `shift`.";
    STR_MODES = "STR003", Warning, "A modal analysis may have missed modes: the Sturm count (the inertia of K − σM) finds more eigenvalues in the range than Lanczos returned, or Lanczos stopped above its tolerance.";
    STR_FIT = "STR004", Error, "The clamp fit (D290) cannot reach the measured frequency: it lies at or above the frequency of the rigidly clamped model (the model is too soft), or the mode does not depend on the spring.";
    CIR_CASE = "CIR001", Error, "A circuit case (SF7b) that is not valid: an unknown key or a wrong value, a tube to a node that does not exist, a circuit without a pressure or drum node, a drum above the critical pressure, a CHF table that is not complete.";
    CIR_CONVERGENCE = "CIR002", Error, "The Newton solve of a circuit did not converge or stalled: give start flows (`start_flow`) near the expected ones, check that every part of the circuit reaches a pressure or drum node, and that the flow can carry the heat.";
    CIR_RANGE = "CIR003", Warning, "A state of a circuit left the range of the IAPWS-IF97 properties (the pressure fell below 1 kPa, or the enthalpy of a tube reached 4.1 MJ/kg): too much flow for the pressure level, or too little flow for the heat.";
    CIR_CHF = "CIR004", Warning, "The critical heat flux ratio q_CHF/q of a tube segment is below the margin of the case (`chf_margin`, default 1.3); below 1 the wall is beyond the CHF and its temperature is that of vapour convection (Dougall-Rohsenow).";
    CIR_COUPLING = "CIR005", Error, "A `[circuits]` table (SF7c) that is not valid: no energy equation, a transient case (SF18), a boundary that is not a wall or carries its own energy table, a tube that is not in the circuit file, a tube axis that does not span the faces of its boundary.";
    CIR_NOT_COUPLED = "CIR006", Warning, "The water and steam circuit behind the walls (SF7c) is not converged with the flow: the coupling residual is above `[circuits] tolerance` at the end, or a circuit solve failed during the run and the water side of the solve before was kept.";
    // anything else
    GENERAL = "GEN001", Error, "An error without a more specific code.";
}

/// A call refused before it ran (an unknown tool, a bad input, a refusal or
/// a declined approval of the MCP server): a replay skips it (D197).
pub fn refused_before_run(code: &str) -> bool {
    matches!(
        code,
        "TOOL001" | "TOOL002" | "MCP002" | "MCP003" | "MCP004" | "MCP005"
    )
}

/// The registry entry of a code.
pub fn info(code: &str) -> Option<&'static CodeInfo> {
    CODES.iter().find(|c| c.code == code)
}

/// A fix: case keys and the values that resolve the finding, as
/// `--set key=value` arguments.
#[derive(Debug, Clone, PartialEq, Serialize, schemars::JsonSchema)]
pub struct Fix {
    /// `key=value` (TOML values), in order.
    pub set: Vec<String>,
}

/// One finding.
#[derive(Debug, Clone, PartialEq, Serialize, schemars::JsonSchema)]
pub struct Diagnostic {
    /// Code of [`CODES`].
    pub code: &'static str,
    /// Severity.
    pub severity: Severity,
    /// Key path in the case file (`fluid.density`), if one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub key: Option<String>,
    /// Boundary (boundary) it concerns, if one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub boundary: Option<String>,
    /// Line and column in the case file (1-based), if known.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub location: Option<[usize; 2]>,
    /// What is wrong and what to change.
    pub message: String,
    /// A fix that `--set` can apply.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fix: Option<Fix>,
}

impl Diagnostic {
    /// A diagnostic with the usual severity of `code`.
    pub fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            severity: info(code).map_or(Severity::Error, |i| i.severity),
            key: None,
            boundary: None,
            location: None,
            message: message.into(),
            fix: None,
        }
    }

    /// With the key path.
    pub fn key(mut self, k: impl Into<String>) -> Self {
        self.key = Some(k.into());
        self
    }

    /// With the boundary.
    pub fn boundary(mut self, b: impl Into<String>) -> Self {
        self.boundary = Some(b.into());
        self
    }

    /// With another severity.
    pub fn severity(mut self, s: Severity) -> Self {
        self.severity = s;
        self
    }

    /// With a fix.
    pub fn fix<S: Into<String>>(mut self, set: impl IntoIterator<Item = S>) -> Self {
        self.fix = Some(Fix {
            set: set.into_iter().map(Into::into).collect(),
        });
        self
    }
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Diagnostic {}

/// An error that carries `d`.
pub fn err(d: Diagnostic) -> anyhow::Error {
    anyhow::Error::new(d)
}

/// Fails with `code` and `message` (and the key path) unless `cond`.
pub fn ensure(
    cond: bool,
    code: &'static str,
    key: &str,
    message: impl FnOnce() -> String,
) -> anyhow::Result<()> {
    if cond {
        Ok(())
    } else {
        let mut d = Diagnostic::new(code, message());
        if !key.is_empty() {
            d = d.key(key);
        }
        Err(err(d))
    }
}

/// The key path that a validation message names: `[time] max_cfl must …`
/// gives `time.max_cfl`, `[time] needs …` gives `time`.
pub fn key_of(message: &str) -> Option<String> {
    let rest = message.strip_prefix('[')?;
    let (section, rest) = rest.split_once(']')?;
    let word: String = rest
        .trim_start()
        .chars()
        .take_while(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || *c == '_')
        .collect();
    let after = rest.trim_start()[word.len()..].chars().next();
    let is_key = !word.is_empty()
        && ![
            "needs",
            "need",
            "is",
            "has",
            "takes",
            "benchmark",
            "reference",
        ]
        .contains(&word.as_str())
        && matches!(after, None | Some(' ' | ',' | ':'));
    Some(if is_key {
        format!("{section}.{word}")
    } else {
        section.to_string()
    })
}

/// Fails with `code`, the message and the key path it names (see
/// [`key_of`]) unless `cond`.
pub fn check(
    cond: bool,
    code: &'static str,
    message: impl FnOnce() -> String,
) -> anyhow::Result<()> {
    check_fix(cond, code, &[], message)
}

/// As [`check`], with a fix (`key=value` for `--set`).
pub fn check_fix(
    cond: bool,
    code: &'static str,
    fix: &[&str],
    message: impl FnOnce() -> String,
) -> anyhow::Result<()> {
    if cond {
        return Ok(());
    }
    let m = message();
    let mut d = Diagnostic::new(code, m.clone());
    d.key = key_of(&m);
    if !fix.is_empty() {
        d = d.fix(fix.iter().copied());
    }
    Err(err(d))
}

/// `need!(cond, CODE, "format", args…)`: [`check`] with a formatted
/// message and `?`.
#[macro_export]
macro_rules! need {
    ($cond:expr, $code:ident, $($fmt:tt)+) => {
        $crate::diag::check($cond, $crate::diag::$code, || format!($($fmt)+))?
    };
}

/// The diagnostic of an error: the first one in its chain, with the full
/// message of the chain (context included). An error without one gets its
/// code from what failed (D202): a mesh file that cannot be read `MESH003`,
/// a TOML error of a case file `CASE001`, `CASE002`, `CASE003` or `CASE009`
/// (the readers of the mesh, scalar, Poisson and benchmark cases), a file
/// that cannot be read `CASE006`; `GEN001` otherwise.
pub fn of_error(e: &anyhow::Error) -> Diagnostic {
    let msg = format!("{e:#}");
    let mut d = e
        .chain()
        .find_map(|c| c.downcast_ref::<Diagnostic>())
        .cloned()
        .or_else(|| {
            e.chain()
                .any(|c| c.is::<foehncraft_mesh::case::ReadError>())
                .then(|| Diagnostic::new(MESH_READ, ""))
        })
        .or_else(|| {
            e.chain()
                .find_map(|c| c.downcast_ref::<toml::de::Error>())
                .map(|t| crate::case::toml_diagnostic("", t))
        })
        .or_else(|| {
            e.chain()
                .any(|c| c.is::<std::io::Error>())
                .then(|| Diagnostic::new(CASE_FILE, ""))
        })
        .unwrap_or_else(|| Diagnostic::new(GENERAL, ""));
    d.message = msg;
    d
}

/// An error of a case reader (flow, mesh, scalar, Poisson and benchmark
/// cases) with its diagnostic (D202): the one of [`of_error`], or
/// `CASE010` when that finds no specific code. The message is that of the
/// whole chain, as before.
pub fn case_error(e: anyhow::Error) -> anyhow::Error {
    let mut d = of_error(&e);
    if d.code == GENERAL {
        d.code = CASE_REFUSED;
        d.severity = Severity::Error;
        if d.key.is_none() {
            d.key = key_of(&d.message);
        }
    }
    err(d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_unique_and_well_formed() {
        let mut seen = std::collections::BTreeSet::new();
        for c in CODES {
            assert!(seen.insert(c.code), "{} twice", c.code);
            let (p, n) = c
                .code
                .split_at(c.code.find(|ch: char| ch.is_ascii_digit()).unwrap());
            assert!(
                [
                    "CASE", "SCHED", "BC", "SCAL", "RUN", "GEN", "DIV", "BAL", "EXPR", "MAT",
                    "REG", "MESH", "TURB", "MRF", "SLD", "ROT", "STAT", "ENE", "THM", "CHT",
                    "TOOL", "LIM", "MCP", "ACC", "CYC", "BUO", "FRZ", "PRS", "SPR", "DAT", "POR",
                    "IMM", "SPE", "CHM", "STR", "MOV", "CIR", "VOF"
                ]
                .contains(&p),
                "{}",
                c.code
            );
            assert_eq!(n.len(), 3, "{}", c.code);
        }
    }

    #[test]
    fn interfaces_doc_lists_every_code() {
        let doc = include_str!("../../../docs/interfaces.md");
        for c in CODES {
            let row = format!("| `{}` | {} | {} |", c.code, c.severity.name(), c.title);
            assert!(
                doc.contains(&row),
                "docs/interfaces.md lacks the row\n{row}"
            );
        }
    }

    #[test]
    fn keys_from_messages() {
        assert_eq!(
            key_of("[time] max_cfl must be > 0").as_deref(),
            Some("time.max_cfl")
        );
        assert_eq!(
            key_of("[time] needs `end_time` or `steps`").as_deref(),
            Some("time")
        );
        assert_eq!(
            key_of("[output.coefficients] velocity and area must be > 0").as_deref(),
            Some("output.coefficients.velocity")
        );
        assert_eq!(key_of("density must be > 0"), None);
    }

    #[test]
    fn diagnostics_travel_inside_errors() {
        let e = err(Diagnostic::new(CASE_RANGE, "density must be > 0").key("fluid.density"));
        let e = e.context("in case file a.toml");
        let d = of_error(&e);
        assert_eq!(d.code, "CASE004");
        assert_eq!(d.key.as_deref(), Some("fluid.density"));
        assert_eq!(d.message, "in case file a.toml: density must be > 0");
        assert_eq!(of_error(&anyhow::anyhow!("x")).code, "GEN001");
    }
}
