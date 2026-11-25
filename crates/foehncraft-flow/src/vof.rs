//! Two immiscible fluids by the volume-of-fluid method (L3, stage MP1).
//!
//! One set of momentum equations carries the mixture; the volume fraction
//! α of phase 1 (α = 1 in phase 1, 0 in phase 2) moves with the flow,
//!
//! ∂α/∂t + ∇·(α u) = 0,
//!
//! and sets the density ρ = α ρ₁ + (1 − α) ρ₂ and the viscosity
//! μ = α μ₁ + (1 − α) μ₂ per cell and face.
//!
//! **A time step.** At the start of each step α is advanced explicitly
//! with the volume fluxes F_f = ṁ_f/ρ_f of the step before, in sub-steps
//! that keep the outflow Courant number of every cell below
//! `max_courant`; the face values come from THINC (default), CICSAM or
//! upwind ([`foehncraft_kernels::vof`]). The mean face value of the
//! sub-steps gives the face density, so that the mass flux ṁ_f = ρ_f F_f
//! that convects the momentum is the one that moved α. The new densities
//! enter the time term of the momentum rows as level n + 1 (level n is the
//! density before the advection). The continuity of the two incompressible
//! phases is that of the volume, Σ_f F_f = 0: the pressure-correction rows
//! take the face coefficients and the fluxes divided by the face density
//! ([`FlowDevice::quotient`]); the corrections of the mass fluxes are the
//! same as for one fluid.
//!
//! **Surface tension** by the continuum surface force (Brackbill, Kothe &
//! Zemach 1992) in balanced-force form (Francois et al. 2006): the face
//! jump σ κ_f (α_N − α_P) is an offset of the reduced pressure jump, beside
//! the gravity offset ½ (b_P + b_N)·d of the balanced-force interpolation
//! (D154), so a static interface of constant curvature is in discrete
//! equilibrium with p′ = σ κ α. The curvature κ = −∇·(∇α̃/|∇α̃|) is taken
//! from α smoothed `smoothing` times; walls and symmetry planes have a
//! contact angle of 90°.
//!
//! **Ledger.** The update is conservative: the volume of phase 1 changes
//! by the boundary fluxes alone, except for the part clipped to keep
//! 0 ≤ α ≤ 1, which is booked and reported (`clipped_volume`).
//!
//! Case file:
//!
//! ```toml
//! [vof]
//! surface_tension = 0.072      # σ (N/m), default 0
//! scheme = "thinc"             # "thinc" (default) | "cicsam" | "upwind"
//! thinc_beta = 3.5             # steepness of THINC
//! max_courant = 0.25           # outflow Courant number of the α sub-steps
//! smoothing = 2                # passes of α̃ for the curvature
//! initial = 0.0                # α of the cells outside every region
//!
//! [[vof.phase]]                # phase 1 (α = 1), then phase 2
//! name = "water"
//! density = 998.2
//! viscosity = 1.0e-3
//!
//! [[vof.phase]]
//! name = "air"
//! density = 1.205
//! viscosity = 1.82e-5
//!
//! [[vof.region]]               # start regions, in order
//! shape = "box"                # "box" | "sphere" | "cylinder" | "field"
//! lower = [0.0, 0.0, 0.0]
//! upper = [0.146, 0.292, 1.0]
//! alpha = 1.0
//!
//! [vof.inflow]                 # α of the inflow per boundary (default 0)
//! inlet = 1.0
//! ```

use crate::device::{Buf, Buf3, FlowDevice, VofBufs, VofOp};
use crate::diag::{self, Diagnostic};
use anyhow::Result;
use foehncraft_core::exact::ExactSumIter as _;
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The face scheme of α.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, schemars::JsonSchema,
)]
#[serde(rename_all = "kebab-case")]
pub enum Scheme {
    /// THINC with the time-integrated face value (Xiao et al. 2005, 2011).
    #[default]
    Thinc,
    /// CICSAM (Ubbink & Issa 1999).
    Cicsam,
    /// First-order upwind (smears the interface; for checks).
    Upwind,
}

impl Scheme {
    fn code(self) -> u32 {
        use foehncraft_kernels::vof as k;
        match self {
            Self::Thinc => k::SCHEME_THINC,
            Self::Cicsam => k::SCHEME_CICSAM,
            Self::Upwind => k::SCHEME_UPWIND,
        }
    }

    /// The name of the case file.
    pub fn name(self) -> &'static str {
        match self {
            Self::Thinc => "thinc",
            Self::Cicsam => "cicsam",
            Self::Upwind => "upwind",
        }
    }
}

/// One `[[vof.phase]]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PhaseSection {
    /// Name of the phase (output, report).
    pub name: String,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Dynamic viscosity μ (Pa s).
    pub viscosity: f64,
}

/// One `[[vof.region]]`: α = `alpha` inside the shape (cell fractions by
/// the quadrature of the immersed bodies, D226), applied in order over the
/// cells' current α.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RegionSection {
    /// Shape.
    pub shape: crate::immersed::Shape,
    /// Cylinder: a point on the axis; sphere: the centre (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub centre: Option<[f64; 3]>,
    /// Radius of a cylinder or sphere (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    /// Axis of a cylinder (default z).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axis: Option<[f64; 3]>,
    /// Length of a cylinder (m); default endless.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub length: Option<f64>,
    /// Lower corner of a box (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lower: Option<[f64; 3]>,
    /// Upper corner of a box (m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper: Option<[f64; 3]>,
    /// The fraction (0 to 1) of a `field` region, an expression of x, y, z.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fraction: Option<crate::expr::Scalar>,
    /// α inside (default 1).
    #[serde(default = "one")]
    pub alpha: f64,
}

fn one() -> f64 {
    1.0
}

/// `[vof]`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VofSection {
    /// The two phases: phase 1 (α = 1) first.
    #[serde(default)]
    pub phase: Vec<PhaseSection>,
    /// Surface tension σ (N/m), default 0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub surface_tension: Option<f64>,
    /// Face scheme of α.
    #[serde(default)]
    pub scheme: Scheme,
    /// Steepness β of THINC (default 3.5).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thinc_beta: Option<f64>,
    /// The largest outflow Courant number of a cell in one α sub-step
    /// (default 0.25, at most 1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_courant: Option<f64>,
    /// Smoothing passes of α for the curvature (default 2).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub smoothing: Option<u32>,
    /// α outside every region at the start (default 0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial: Option<f64>,
    /// Start regions.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub region: Vec<RegionSection>,
    /// α of the inflow per boundary (default 0).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inflow: BTreeMap<String, f64>,
}

impl VofSection {
    /// True without `[vof]`.
    pub fn is_empty(&self) -> bool {
        self.phase.is_empty()
    }

    /// σ (N/m).
    pub fn sigma(&self) -> f64 {
        self.surface_tension.unwrap_or(0.0)
    }

    /// β of THINC.
    pub fn beta(&self) -> f64 {
        self.thinc_beta.unwrap_or(3.5)
    }

    /// The Courant number of a sub-step.
    pub fn courant(&self) -> f64 {
        self.max_courant.unwrap_or(0.25)
    }

    /// Smoothing passes.
    pub fn passes(&self) -> u32 {
        self.smoothing.unwrap_or(2)
    }
}

fn bad(key: &str, msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::VOF_CONFIG, msg).key(key.to_string()))
}

fn unsupported(msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(diag::VOF_UNSUPPORTED, msg))
}

/// Ranges of `[vof]` and the settings it does not work with yet (VOF001,
/// VOF002).
pub fn validate(case: &crate::case::FlowCase) -> Result<()> {
    let v = &case.vof;
    if v.is_empty() {
        return Ok(());
    }
    if v.phase.len() != 2 {
        return Err(bad(
            "vof.phase",
            format!(
                "[vof] needs exactly two [[vof.phase]] tables, not {}",
                v.phase.len()
            ),
        ));
    }
    for (i, p) in v.phase.iter().enumerate() {
        if !(p.density > 0.0
            && p.density.is_finite()
            && p.viscosity > 0.0
            && p.viscosity.is_finite())
        {
            return Err(bad(
                "vof.phase",
                format!(
                    "[[vof.phase]] {} (`{}`): density and viscosity must be > 0",
                    i + 1,
                    p.name
                ),
            ));
        }
    }
    let s = v.sigma();
    if !(s >= 0.0 && s.is_finite()) {
        return Err(bad(
            "vof.surface_tension",
            format!("[vof] surface_tension {s} must be >= 0"),
        ));
    }
    let b = v.beta();
    if !(b > 0.0 && b <= 20.0) {
        return Err(bad(
            "vof.thinc_beta",
            format!("[vof] thinc_beta {b} must be in (0, 20]"),
        ));
    }
    let c = v.courant();
    if !(c > 0.0 && c <= 1.0) {
        return Err(bad(
            "vof.max_courant",
            format!("[vof] max_courant {c} must be in (0, 1]"),
        ));
    }
    if v.passes() > 10 {
        return Err(bad(
            "vof.smoothing",
            "[vof] smoothing must be at most 10".into(),
        ));
    }
    let a0 = v.initial.unwrap_or(0.0);
    if !(0.0..=1.0).contains(&a0) {
        return Err(bad(
            "vof.initial",
            format!("[vof] initial {a0} must be in [0, 1]"),
        ));
    }
    for (i, r) in v.region.iter().enumerate() {
        if !(0.0..=1.0).contains(&r.alpha) {
            return Err(bad(
                "vof.region",
                format!(
                    "[[vof.region]] {}: alpha {} must be in [0, 1]",
                    i + 1,
                    r.alpha
                ),
            ));
        }
    }
    for (name, a) in &v.inflow {
        if !(0.0..=1.0).contains(a) {
            return Err(bad(
                "vof.inflow",
                format!("[vof.inflow] {name} = {a} must be in [0, 1]"),
            ));
        }
        if !case.boundary.contains_key(name) {
            return Err(bad(
                "vof.inflow",
                format!("[vof.inflow] names the boundary `{name}`, which the case does not have"),
            ));
        }
    }
    // what the first version does not have
    let mut no = Vec::new();
    if case.time.is_none() {
        no.push("a steady case (VOF needs [time])");
    }
    if case.turbulence.model != crate::turbulence::TurbulenceModel::Laminar {
        no.push("a turbulence model");
    }
    if case.energy.is_some() {
        no.push("the energy equation (stage MP2)");
    }
    if !case.scalars.is_empty() {
        no.push("passive scalars or species");
    }
    if case.buoyancy.model != crate::buoyancy::Model::None {
        no.push("a [buoyancy] model (gravity acts through [gravity] alone)");
    }
    if !case.regions.is_empty() {
        no.push("regions or several fluids");
    }
    if !case.spray.is_empty() {
        no.push("sprays");
    }
    if !case.immersed.is_empty() {
        no.push("immersed bodies");
    }
    if case.solver.algorithm == Some(crate::simplec::Algorithm::Coupled) {
        no.push("the coupled algorithm (PISO and SIMPLEC only)");
    }
    if !case.mrf.is_empty() || !case.sliding.is_empty() {
        no.push("rotating zones or sliding interfaces");
    }
    if case.mesh_motion.is_some() {
        no.push("a moving or deforming mesh");
    }
    if case.run.restart.is_some() {
        no.push("a restart");
    }
    if case.boundary.values().any(|b| b.is_periodic()) {
        no.push("periodic boundaries");
    }
    if !no.is_empty() {
        return Err(unsupported(format!(
            "[vof] does not work with {} in this version",
            no.join(", ")
        )));
    }
    Ok(())
}

/// The volume-of-fluid state of a flow.
#[derive(Debug)]
pub struct Vof {
    /// The settings.
    pub section: VofSection,
    rho: [f64; 2],
    mu: [f64; 2],
    /// α per cell.
    pub alpha: Buf,
    alpha_start: Buf,
    alpha_new: Buf,
    alpha_f: Buf,
    alpha_f_mean: Buf,
    flux: Buf,
    inv_rho_f: Buf,
    /// ρ per face on the device (the divisor of the volume fluxes).
    pub rho_f: Buf,
    grad: Buf3,
    gnorm: Buf,
    clipped: Buf,
    work: Buf,
    work2: Buf,
    kappa: Buf,
    tension: Buf,
    alpha_b: Buf,
    /// The raw face offsets o_f of gravity and surface tension (device).
    pub offset: Buf,
    /// ρ per cell (device).
    pub rho_c: Buf,
    /// The balanced body force B_P per unit volume (device).
    pub body: Buf3,
    /// Scratch of the volume-based pressure rows: c_f/ρ_f.
    pub cf_vol: Buf,
    /// Scratch: ṁ_f/ρ_f.
    pub m_vol: Buf,
    /// ρ per cell of the last mix (host).
    pub rho_cell: Vec<f64>,
    /// ρ per face of the last mix (host).
    pub rho_face: Vec<f64>,
    /// α per cell of the last mix (host).
    pub alpha_host: Vec<f64>,
    /// α per face of the last mix (host).
    alpha_face: Vec<f64>,
    /// ρ at the time levels n and n − 1.
    old: Option<[Vec<f64>; 2]>,
    saved: Option<Saved>,
    /// σ κ_f (α_N − α_P) per face of the last refresh.
    pub tension_offset: Vec<f64>,
    /// The totals of the run.
    pub totals: Totals,
    volume: Vec<f64>,
}

/// What [`Vof::step_back`] restores.
#[derive(Debug)]
struct Saved {
    old: Option<[Vec<f64>; 2]>,
    rho_cell: Vec<f64>,
    rho_face: Vec<f64>,
    alpha: Vec<f64>,
    alpha_face: Vec<f64>,
    totals: Totals,
}

/// The totals of a run (for `result.json` and `vof.csv`).
#[derive(Debug, Clone, Default, Serialize, Deserialize, schemars::JsonSchema, PartialEq)]
pub struct Totals {
    /// Volume of phase 1 at the start (m³).
    pub volume_start: f64,
    /// Volume of phase 1 now (m³).
    pub volume: f64,
    /// Volume of phase 1 that left through the boundaries (m³, net).
    pub volume_out: f64,
    /// Volume clipped to keep 0 ≤ α ≤ 1 (m³, signed sum).
    pub clipped_volume: f64,
    /// Sub-steps of the last step.
    pub substeps: usize,
    /// Most sub-steps in a step.
    pub max_substeps: usize,
    /// Smallest and largest α after the last step.
    pub alpha_range: [f64; 2],
}

/// The start field of α from `[vof] initial` and the regions.
pub fn start_alpha(case: &crate::case::FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<Vec<f64>> {
    let v = &case.vof;
    let nc = mesh.n_cells;
    let mut a = vec![v.initial.unwrap_or(0.0); nc];
    if v.region.is_empty() {
        return Ok(a);
    }
    let sec = crate::immersed::ImmersedSection {
        // (the permeability is not used: only the cell fractions)
        permeability: Some(1.0),
        samples: Some(6),
        body: v
            .region
            .iter()
            .enumerate()
            .map(|(i, r)| crate::immersed::BodySection {
                name: format!("vof-region-{}", i + 1),
                shape: r.shape,
                centre: r.centre,
                radius: r.radius,
                axis: r.axis,
                length: r.length,
                lower: r.lower,
                upper: r.upper,
                fraction: r.fraction.clone(),
            })
            .collect(),
        topology: None,
    };
    let bodies = crate::immersed::bodies_with(&sec, &case.env()?)
        .map_err(|e| bad("vof.region", format!("[[vof.region]]: {e:#}")))?;
    let fr = crate::immersed::fractions(&bodies, 6, mesh, geo);
    for (r, f) in v.region.iter().zip(&fr) {
        for c in 0..nc {
            a[c] = a[c] * (1.0 - f[c]) + r.alpha * f[c];
        }
    }
    Ok(a)
}

impl Vof {
    /// The VOF state of `case` with the start field `alpha` (per cell) on a
    /// flow device; `open` marks the boundary faces whose curvature takes
    /// the owner's normal, `alpha_b` is the inflow α per boundary face.
    pub fn new(
        section: VofSection,
        dev: &mut dyn FlowDevice,
        disc: &foehncraft_fv::disc::Discretisation,
        alpha: &[f64],
        open: &[u32],
        alpha_b: &[f64],
    ) -> Result<Self> {
        let (nc, nf) = (disc.n_cells, disc.n_faces);
        let cell = |dev: &mut dyn FlowDevice| dev.alloc(nc);
        let face = |dev: &mut dyn FlowDevice| dev.alloc(nf);
        let a = cell(dev)?;
        dev.upload(a, alpha)?;
        let grad = Buf3 {
            x: cell(dev)?,
            y: cell(dev)?,
            z: cell(dev)?,
        };
        let ab = dev.alloc(nf - disc.n_interior)?;
        dev.upload(ab, alpha_b)?;
        dev.vof_open(open)?;
        let p = &section.phase;
        let rho = [p[0].density, p[1].density];
        let mu = [p[0].viscosity, p[1].viscosity];
        let mut v = Self {
            rho,
            mu,
            alpha: a,
            alpha_start: cell(dev)?,
            alpha_new: cell(dev)?,
            alpha_f: face(dev)?,
            alpha_f_mean: face(dev)?,
            flux: face(dev)?,
            inv_rho_f: face(dev)?,
            rho_f: face(dev)?,
            grad,
            gnorm: cell(dev)?,
            clipped: cell(dev)?,
            work: cell(dev)?,
            work2: cell(dev)?,
            kappa: cell(dev)?,
            tension: face(dev)?,
            offset: face(dev)?,
            rho_c: cell(dev)?,
            body: Buf3 {
                x: cell(dev)?,
                y: cell(dev)?,
                z: cell(dev)?,
            },
            alpha_b: ab,
            cf_vol: face(dev)?,
            m_vol: face(dev)?,
            rho_cell: Vec::new(),
            rho_face: Vec::new(),
            alpha_host: alpha.to_vec(),
            alpha_face: Vec::new(),
            old: None,
            saved: None,
            tension_offset: vec![0.0; nf],
            totals: Totals::default(),
            volume: disc.cell_volume[..nc].to_vec(),
            section,
        };
        // the start face values: linear between the cells
        let af = linear_faces(disc, alpha);
        dev.upload(v.alpha_f_mean, &af)?;
        v.mix(alpha, &af);
        v.totals.volume_start = v.phase_volume(alpha);
        v.totals.volume = v.totals.volume_start;
        v.totals.alpha_range = range(alpha);
        Ok(v)
    }

    /// Sets the volume fluxes F_f (m³/s) of the next advance (a prescribed
    /// velocity: the advection checks).
    pub fn set_flux(&mut self, dev: &mut dyn FlowDevice, flux: &[f64]) -> Result<()> {
        dev.upload(self.flux, flux)
    }

    /// A face buffer for the density ratio of a step start (free outside
    /// the pressure rows).
    pub fn ratio_buf(&self) -> Buf {
        self.cf_vol
    }

    /// The buffer of the volume fluxes.
    pub fn flux_buf(&self) -> Buf {
        self.flux
    }

    /// α per face of the last mix (host).
    pub fn face_alpha(&self) -> Vec<f64> {
        self.alpha_face.clone()
    }

    /// The densities and viscosities of phase 1 and 2.
    pub fn properties(&self) -> ([f64; 2], [f64; 2]) {
        (self.rho, self.mu)
    }

    fn phase_volume(&self, a: &[f64]) -> f64 {
        a.iter().zip(&self.volume).map(|(a, v)| a * v).exact_sum()
    }

    /// ρ and μ from α per cell and per face (host).
    fn mix(&mut self, a: &[f64], af: &[f64]) {
        let [r1, r2] = self.rho;
        self.rho_cell = a.iter().map(|&x| x * r1 + (1.0 - x) * r2).collect();
        self.rho_face = af.iter().map(|&x| x * r1 + (1.0 - x) * r2).collect();
        self.alpha_host = a.to_vec();
        self.alpha_face = af.to_vec();
    }

    /// μ per face of the last mix.
    pub fn mu_face(&self, af: &[f64]) -> Vec<f64> {
        let [m1, m2] = self.mu;
        af.iter().map(|&x| x * m1 + (1.0 - x) * m2).collect()
    }

    /// Advances α over a step `dt` with the volume fluxes of the current
    /// mass flux `mass_flux` (divided by the face densities of the last
    /// mix): sub-steps, clipping, the mean face values. Returns the face
    /// values α_f (mean over the sub-steps).
    pub fn advance(
        &mut self,
        dev: &mut dyn FlowDevice,
        mass_flux: Buf,
        dt: f64,
        boundary_start: usize,
    ) -> Result<Vec<f64>> {
        dev.quotient(self.flux, 1.0, mass_flux, self.rho_f)?;
        self.advance_flux(dev, dt, boundary_start)
    }

    /// [`Vof::advance`] with the volume fluxes already in their buffer
    /// ([`Vof::set_flux`]).
    pub fn advance_flux(
        &mut self,
        dev: &mut dyn FlowDevice,
        dt: f64,
        boundary_start: usize,
    ) -> Result<Vec<f64>> {
        dev.copy(self.alpha, self.alpha_start)?;
        dev.vof(
            VofOp::OutflowRate,
            &VofBufs {
                flux: Some(self.flux),
                out: Some(self.work),
                ..Default::default()
            },
        )?;
        let rate = dev.min_max(self.work)?.1;
        let co = rate * dt;
        let n = if co.is_finite() && co > 0.0 {
            ((co / self.section.courant()).ceil() as usize).max(1)
        } else {
            1
        };
        anyhow::ensure!(
            n <= 100_000,
            "[vof]: {n} sub-steps of α in one time step (outflow Courant number {co:.3e}); \
             the time step is far too long for the velocity"
        );
        let h = dt / n as f64;
        dev.fill(self.alpha_f_mean, 0.0)?;
        let mut clip = 0.0;
        for _ in 0..n {
            dev.vof(
                VofOp::Gradient,
                &VofBufs {
                    alpha: Some(self.alpha),
                    out: Some(self.gnorm),
                    out3: Some(self.grad),
                    ..Default::default()
                },
            )?;
            dev.vof(
                VofOp::Faces,
                &VofBufs {
                    flux: Some(self.flux),
                    alpha: Some(self.alpha),
                    grad: Some(self.grad),
                    alpha_b: Some(self.alpha_b),
                    scheme: self.section.scheme.code(),
                    dt: h,
                    beta: self.section.beta(),
                    out: Some(self.alpha_f),
                    ..Default::default()
                },
            )?;
            dev.vof(
                VofOp::Update,
                &VofBufs {
                    flux: Some(self.flux),
                    alpha_f: Some(self.alpha_f),
                    alpha: Some(self.alpha),
                    dt: h,
                    out: Some(self.alpha_new),
                    out2: Some(self.clipped),
                    ..Default::default()
                },
            )?;
            dev.copy(self.alpha_new, self.alpha)?;
            dev.lincomb(
                self.alpha_f_mean,
                1.0,
                self.alpha_f_mean,
                1.0 / n as f64,
                self.alpha_f,
            )?;
            let c = dev.download(self.clipped)?;
            clip += c.iter().copied().exact_sum();
            // the volume of phase 1 through the boundary faces
            let af = dev.download(self.alpha_f)?;
            let fl = dev.download(self.flux)?;
            self.totals.volume_out += (boundary_start..fl.len())
                .map(|f| fl[f] * af[f] * h)
                .exact_sum();
        }
        self.totals.clipped_volume += clip;
        self.totals.substeps = n;
        self.totals.max_substeps = self.totals.max_substeps.max(n);
        let a = dev.download(self.alpha)?;
        let af = dev.download(self.alpha_f_mean)?;
        self.totals.volume = self.phase_volume(&a);
        self.totals.alpha_range = range(&a);
        self.mix(&a, &af);
        Ok(af)
    }

    /// Starts a step: the density of now becomes level n (level n − 1 the
    /// one before); keeps what [`Vof::step_back`] restores.
    pub fn start_step(&mut self) {
        let prev = self.old.take();
        self.saved = Some(Saved {
            old: prev.clone(),
            rho_cell: self.rho_cell.clone(),
            rho_face: self.rho_face.clone(),
            alpha: self.alpha_host.clone(),
            alpha_face: self.alpha_face.clone(),
            totals: self.totals.clone(),
        });
        let n1 = match prev {
            Some([n, _]) => n,
            None => self.rho_cell.clone(),
        };
        self.old = Some([self.rho_cell.clone(), n1]);
    }

    /// The old density levels of the time term (n, n − 1).
    pub fn old_levels(&self) -> Option<&[Vec<f64>; 2]> {
        self.old.as_ref()
    }

    /// Takes the last step back (a rejected step): α, the densities and
    /// the levels as before [`Vof::start_step`].
    pub fn step_back(&mut self, dev: &mut dyn FlowDevice) -> Result<()> {
        if let Some(sv) = self.saved.take() {
            self.old = sv.old;
            self.rho_cell = sv.rho_cell;
            self.rho_face = sv.rho_face;
            self.alpha_host = sv.alpha;
            self.alpha_face = sv.alpha_face;
            self.totals = sv.totals;
            dev.copy(self.alpha_start, self.alpha)?;
            // the volume fluxes of the start of the step (ṁ/ρ_f of the
            // restored densities follow in `upload_faces`)
        }
        Ok(())
    }

    /// Uploads 1/ρ_f and ρ_f of the last mix and rescales the mass flux to
    /// the new face densities: ṁ_f = ρ_f F_f with the volume flux of the
    /// advection.
    pub fn upload_faces(&mut self, dev: &mut dyn FlowDevice, mass_flux: Option<Buf>) -> Result<()> {
        let inv: Vec<f64> = self.rho_face.iter().map(|r| 1.0 / r).collect();
        dev.upload(self.inv_rho_f, &inv)?;
        dev.upload(self.rho_f, &self.rho_face)?;
        if let Some(m) = mass_flux {
            dev.quotient(m, 1.0, self.flux, self.inv_rho_f)?;
        }
        Ok(())
    }

    /// The surface-tension offsets σ κ_f (α_N − α_P) per face from the
    /// current α (zero for σ = 0), into [`Vof::tension_offset`].
    pub fn refresh_tension(&mut self, dev: &mut dyn FlowDevice, delta: f64) -> Result<()> {
        let sigma = self.section.sigma();
        if sigma == 0.0 {
            return Ok(());
        }
        // α̃: `smoothing` passes, alternating between two buffers
        let mut src = self.alpha;
        let bufs = [self.work, self.work2];
        for i in 0..self.section.passes() as usize {
            let dst = bufs[i % 2];
            dev.vof(
                VofOp::Smooth,
                &VofBufs {
                    alpha: Some(src),
                    out: Some(dst),
                    ..Default::default()
                },
            )?;
            src = dst;
        }
        dev.vof(
            VofOp::Gradient,
            &VofBufs {
                alpha: Some(src),
                out: Some(self.gnorm),
                out3: Some(self.grad),
                ..Default::default()
            },
        )?;
        dev.vof(
            VofOp::Curvature,
            &VofBufs {
                grad: Some(self.grad),
                s: delta,
                out: Some(self.kappa),
                ..Default::default()
            },
        )?;
        dev.vof(
            VofOp::Tension,
            &VofBufs {
                alpha: Some(self.alpha),
                kappa: Some(self.kappa),
                gnorm: Some(self.gnorm),
                s: sigma,
                out: Some(self.tension),
                ..Default::default()
            },
        )?;
        self.tension_offset = dev.download(self.tension)?;
        Ok(())
    }

    /// The curvature per cell of the last [`Vof::refresh_tension`].
    pub fn curvature(&self, dev: &mut dyn FlowDevice) -> Result<Vec<f64>> {
        dev.download(self.kappa)
    }
}

fn range(a: &[f64]) -> [f64; 2] {
    let lo = a.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = a.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    [lo, hi]
}

/// The linear face values of a cell field (the owner value on boundary
/// faces).
pub fn linear_faces(disc: &foehncraft_fv::disc::Discretisation, a: &[f64]) -> Vec<f64> {
    (0..disc.n_faces)
        .map(|f| {
            let p = disc.owner[f] as usize;
            if f < disc.n_interior {
                let n = disc.neighbour[f] as usize;
                let w = disc.weight[f];
                w * a[p] + (1.0 - w) * a[n]
            } else {
                a[p]
            }
        })
        .collect()
}

/// The open-boundary flags (1: inlets and outlets, 0: walls and symmetry)
/// and the inflow α per boundary face of a case.
pub fn boundary_arrays(
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    n_interior: usize,
    n_faces: usize,
) -> (Vec<u32>, Vec<f64>) {
    let ni = n_interior;
    let nb = n_faces - ni;
    let mut open = vec![0u32; nb];
    let mut ab = vec![0.0; nb];
    for b in &mesh.boundaries {
        let a = case.vof.inflow.get(&b.name).copied().unwrap_or(0.0);
        let o = case.boundary.get(&b.name).is_some_and(|c| c.is_open());
        for f in b.range() {
            if f >= ni && f - ni < nb {
                ab[f - ni] = a;
                open[f - ni] = u32::from(o);
            }
        }
    }
    (open, ab)
}

/// Attaches `[vof]` of `case` to the flow (after the start fields; nothing
/// without `[vof]`).
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &crate::case::FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
) -> Result<()> {
    if case.vof.is_empty() {
        return Ok(());
    }
    let alpha = start_alpha(case, mesh, geo)?;
    let (open, ab) = boundary_arrays(case, mesh, fl.disc.n_interior, fl.disc.n_faces);
    let gravity = crate::buoyancy::hydrostatic(case).map(|h| (h.g, h.rho_ref));
    fl.attach_vof(case.vof.clone(), &alpha, &open, &ab, gravity)
}

/// The result of [`advect`].
#[derive(Debug, Clone)]
pub struct Advected {
    /// The cell volumes.
    pub volume: Vec<f64>,
    /// α at the start.
    pub start: Vec<f64>,
    /// α at the end.
    pub end: Vec<f64>,
    /// The totals (volume, clipped volume, sub-steps).
    pub totals: Totals,
}

/// The advection of the start α of `case` alone, in the prescribed
/// velocity `u(x, t)` (m/s) for `steps` steps of `dt`, on the flow device
/// of the case (the kernels of the run; verification: Zalesak's disk, the
/// shear reversal). The face fluxes are u(x_f, t + Δt/2)·S_f.
pub fn advect(
    case: &crate::case::FlowCase,
    u: &dyn Fn([f64; 3], f64) -> [f64; 3],
    dt: f64,
    steps: usize,
    threads: Option<usize>,
) -> Result<Advected> {
    let (mesh, geo, _) = crate::case::mesh_for(case, None)?;
    let (mesh, geo) = (&mesh, &geo);
    let phys = crate::simplec::Physics::constant(case.rho(), case.mu(), [0.0; 3], true);
    let arrays = crate::bc::arrays(
        mesh,
        geo,
        &case.boundary,
        phys.rho,
        crate::bc::At::numbers(None),
    )?;
    let disc = crate::simplec::discretisation(mesh, geo, &arrays)?;
    let mut fl = crate::simplec::Flow::new(
        disc,
        arrays,
        phys,
        case.solver.clone(),
        [0.0; 3],
        0.0,
        Some(0),
        &crate::simplec::Target {
            device: "cpu".into(),
            threads,
            ptx: Default::default(),
            ranks: None,
        },
    )?;
    let alpha = start_alpha(case, mesh, geo)?;
    let (open, ab) = boundary_arrays(case, mesh, fl.disc.n_interior, fl.disc.n_faces);
    let mut v = Vof::new(
        case.vof.clone(),
        fl.dev.as_mut(),
        &fl.disc,
        &alpha,
        &open,
        &ab,
    )?;
    let nf = fl.disc.n_faces;
    let ni = fl.disc.n_interior;
    let mut t = 0.0;
    for _ in 0..steps {
        let tm = t + 0.5 * dt;
        let flux: Vec<f64> = (0..nf)
            .map(|f| {
                let x = geo.face_centroid[f];
                let s = geo.face_area[f];
                let w = u(x, tm);
                w[0] * s[0] + w[1] * s[1] + w[2] * s[2]
            })
            .collect();
        v.set_flux(fl.dev.as_mut(), &flux)?;
        v.advance_flux(fl.dev.as_mut(), dt, ni)?;
        t += dt;
    }
    Ok(Advected {
        volume: geo.cell_volume[..mesh.n_cells].to_vec(),
        start: alpha,
        end: v.alpha_host.clone(),
        totals: v.totals.clone(),
    })
}

/// Units of `vof.csv`.
pub const CSV_UNITS: &str = "step -, time s, volume m3 (phase 1), volume_out m3 (net through the \
    boundaries), clipped m3, substeps -, alpha_min alpha_max -, x_min..z_max m (extent of the \
    cells with alpha > 0.5)";

/// The header of `vof.csv`.
pub const CSV_HEAD: &str = "step,time,volume,volume_out,clipped,substeps,alpha_min,alpha_max,\
    x_min,x_max,y_min,y_max,z_min,z_max";

impl Vof {
    /// One row of `vof.csv` after step `n` at `time`, with the cell
    /// centroids `x`.
    pub fn csv_row(&self, n: usize, time: f64, x: &[[f64; 3]]) -> String {
        let t = &self.totals;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (c, &a) in self.alpha_host.iter().enumerate() {
            if a > 0.5 {
                for i in 0..3 {
                    lo[i] = lo[i].min(x[c][i]);
                    hi[i] = hi[i].max(x[c][i]);
                }
            }
        }
        format!(
            "{n},{time:e},{:e},{:e},{:e},{},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e}",
            t.volume,
            t.volume_out,
            t.clipped_volume,
            t.substeps,
            t.alpha_range[0],
            t.alpha_range[1],
            lo[0],
            hi[0],
            lo[1],
            hi[1],
            lo[2],
            hi[2]
        )
    }
}
