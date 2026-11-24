//! Porous regions (D270, porous plan stage PM1): a cell zone filled with a
//! fixed porous matrix (a honeycomb, a packed bed, a filter) that the gas
//! flows through. The velocity of the porous cells is the superficial
//! (Darcy) velocity, so the mass flux is continuous across the interface
//! to the clear fluid; the matrix adds the Darcy–Forchheimer sink
//!
//!   S = −(μ K⁻¹ + ½ ρ |u| F) · u
//!
//! per unit volume (Nield & Bejan 2017, ch. 1 and 2; Whitaker 1999), with
//! the tensors K⁻¹ (1/m²) and F (1/m) given by principal values along the
//! axes of the region (the material frame of L48: `[regions.<name>] axes`,
//! the first axis along the channels of a honeycomb). The convective and
//! time terms keep the form of the clear fluid (the superficial
//! formulation without the 1/ε and 1/ε² factors of the volume-averaged
//! equation): the sink dominates in the matrix and the two forms give the
//! same steady pressure drop.
//!
//! In the momentum rows the sink is implicit per component; the factors D
//! of the face fluxes take the smallest coefficient (along the channels).
//! Direction-dependent resistances need axes along the mesh axes in this
//! version (POR002): with turned axes the mesh-frame diagonal is dominated
//! by the blocked directions, the scalar D then under-corrects the flow
//! along the channels by the blocking ratio and the iteration stalls
//! (tried on a honeycomb turned by 30°); a tensor factor D of the face
//! fluxes would lift this. The code keeps the off-diagonal part, deferred
//! with its size on the diagonal (exact at convergence), for that step.
//!
//! `honeycomb = { channel, pitch }` computes the values of square channels
//! from laminar duct flow (Shah & London 1978): f_D Re = 56.9 gives
//! K⁻¹ along the channels = (f_D Re / 2)/(ε d_h²) with d_h the channel
//! width and ε = (channel/pitch)² the open area; across the channels the
//! walls are closed: `blocking` (default 10⁴) times the value along them.
//! With `length` (the length of one brick) the incremental pressure drop of
//! the developing flow at each brick entrance, K(∞) = 1.43, adds
//! F = K(∞)/(ε² length) along the channels.

use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::materials::RegionKind;
use anyhow::Result;
use serde::{Deserialize, Serialize};

/// f_D Re of fully developed laminar flow in a square duct (Shah & London
/// 1978).
pub const SQUARE_DUCT_FRE: f64 = 56.9;
/// The incremental pressure drop K(∞) of the developing laminar flow in a
/// square duct (Shah & London 1978).
pub const SQUARE_DUCT_K_INF: f64 = 1.43;
/// The default ratio of the resistance across the channels of a honeycomb
/// to the resistance along them.
pub const DEFAULT_BLOCKING: f64 = 1.0e4;

/// `[regions.<name>.porous]`: the matrix of a porous region.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PorousSection {
    /// Open volume fraction ε (0 < ε ≤ 1); with `honeycomb` it follows
    /// from the channel and the pitch.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub porosity: Option<f64>,
    /// Viscous resistance K⁻¹ (1/m²): three principal values along the
    /// axes of the region, or one value for every direction.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub viscous: Option<Principal>,
    /// Inertial resistance F (1/m), as `viscous`; the sink is ½ ρ |u| F u.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inertial: Option<Principal>,
    /// The values of a honeycomb of square channels (instead of `viscous`,
    /// `inertial` and `porosity`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub honeycomb: Option<Honeycomb>,
    /// Where the values come from.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
}

/// Principal values of a resistance: one number for all directions or
/// three along the axes of the region.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(untagged)]
pub enum Principal {
    /// The same value in every direction.
    One(f64),
    /// Along the first, second and third axis.
    Three([f64; 3]),
}

impl Principal {
    fn values(self) -> [f64; 3] {
        match self {
            Self::One(v) => [v; 3],
            Self::Three(v) => v,
        }
    }
}

/// `honeycomb = { channel, pitch, length, blocking }`: square channels
/// along the first axis of the region.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Honeycomb {
    /// Width of a square channel (m), also its hydraulic diameter.
    pub channel: f64,
    /// Distance between the channel axes (m): channel plus one wall.
    pub pitch: f64,
    /// Length of one brick along the channels (m): the entrance loss of
    /// each brick (optional; left out: fully developed flow only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub length: Option<f64>,
    /// Resistance across the channels over the resistance along them
    /// (default 10⁴: closed walls).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub blocking: Option<f64>,
}

/// The resolved values of a porous region in its own frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Resolved {
    /// ε.
    pub porosity: f64,
    /// K⁻¹ along the three axes (1/m²).
    pub viscous: [f64; 3],
    /// F along the three axes (1/m).
    pub inertial: [f64; 3],
}

impl PorousSection {
    /// The principal values; `Err` with the reason for missing or
    /// contradicting keys.
    pub fn resolve(&self) -> std::result::Result<Resolved, String> {
        if let Some(h) = &self.honeycomb {
            if self.viscous.is_some() || self.inertial.is_some() || self.porosity.is_some() {
                return Err(
                    "give either `honeycomb` or `viscous`, `inertial` and `porosity`, \
                            not both"
                        .into(),
                );
            }
            if !(h.channel > 0.0 && h.pitch > h.channel) {
                return Err(format!(
                    "honeycomb channel = {} m and pitch = {} m: the pitch is the channel plus \
                     one wall, so 0 < channel < pitch",
                    h.channel, h.pitch
                ));
            }
            let blocking = h.blocking.unwrap_or(DEFAULT_BLOCKING);
            if blocking.is_nan() || blocking < 1.0 {
                return Err(format!("honeycomb blocking = {blocking}: at least 1"));
            }
            if let Some(l) = h.length
                && (l.is_nan() || l <= 0.0)
            {
                return Err(format!("honeycomb length = {l} m: positive"));
            }
            let eps = (h.channel / h.pitch).powi(2);
            let along = 0.5 * SQUARE_DUCT_FRE / (eps * h.channel * h.channel);
            let f = h
                .length
                .map_or(0.0, |l| SQUARE_DUCT_K_INF / (eps * eps * l));
            return Ok(Resolved {
                porosity: eps,
                viscous: [along, blocking * along, blocking * along],
                inertial: [f, blocking * f, blocking * f],
            });
        }
        let viscous = self.viscous.map_or([0.0; 3], Principal::values);
        let inertial = self.inertial.map_or([0.0; 3], Principal::values);
        if self.viscous.is_none() && self.inertial.is_none() {
            return Err("give `viscous` (1/m²) and/or `inertial` (1/m), or `honeycomb`".into());
        }
        if viscous
            .iter()
            .chain(&inertial)
            .any(|v| !(v.is_finite() && *v >= 0.0))
        {
            return Err("the resistances are finite and not negative".into());
        }
        let porosity = self.porosity.unwrap_or(1.0);
        if !(porosity > 0.0 && porosity <= 1.0) {
            return Err(format!("porosity = {porosity}: 0 < ε ≤ 1"));
        }
        Ok(Resolved {
            porosity,
            viscous,
            inertial,
        })
    }
}

/// Whether the case has a porous region.
pub fn any(case: &FlowCase) -> bool {
    case.regions.values().any(|r| r.kind == RegionKind::Porous)
}

/// The checks of the porous regions (POR001).
pub fn validate(case: &FlowCase) -> Result<()> {
    for (name, r) in &case.regions {
        let key = format!("regions.{name}.porous");
        if r.kind != RegionKind::Porous {
            if r.porous.is_some() {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::POROUS,
                        format!(
                            "[regions.{name}.porous] in a region of kind \"{}\": only a \
                             porous region has a matrix (kind = \"porous\")",
                            match r.kind {
                                RegionKind::Fluid => "fluid",
                                RegionKind::Solid => "solid",
                                RegionKind::Porous => "porous",
                            }
                        ),
                    )
                    .key(&key),
                ));
            }
            continue;
        }
        let Some(p) = &r.porous else {
            return Err(diag::err(
                Diagnostic::new(
                    diag::POROUS,
                    format!(
                        "[regions.{name}] kind = \"porous\" needs its matrix: add \
                         [regions.{name}.porous] with `viscous` and `inertial`, or `honeycomb`"
                    ),
                )
                .key(&key),
            ));
        };
        let v = p.resolve().map_err(|m| {
            diag::err(
                Diagnostic::new(diag::POROUS, format!("[regions.{name}.porous] {m}")).key(&key),
            )
        })?;
        let isotropic = |a: [f64; 3]| a[1] == a[0] && a[2] == a[0];
        if r.axes.is_none() && !(isotropic(v.viscous) && isotropic(v.inertial)) {
            return Err(diag::err(
                Diagnostic::new(
                    diag::POROUS,
                    format!(
                        "[regions.{name}] has direction-dependent resistances but no axes: add \
                         `axes = [[..], [..]]` (the first along the channels)"
                    ),
                )
                .key(format!("regions.{name}.axes")),
            ));
        }
        if let Some(a) = r.axes {
            let e = frame(a).map_err(|m| {
                diag::err(
                    Diagnostic::new(diag::POROUS, format!("[regions.{name}] axes: {m}"))
                        .key(format!("regions.{name}.axes")),
                )
            })?;
            // axes turned against the mesh need the tensor factor of the
            // face fluxes (porous plan, D270): not in this version
            let along_mesh = e
                .iter()
                .all(|v| v.iter().filter(|x| x.abs() > 1e-9).count() == 1);
            if !(isotropic(v.viscous) && isotropic(v.inertial)) && !along_mesh {
                return Err(diag::err(
                    Diagnostic::new(
                        diag::POROUS_UNSUPPORTED,
                        format!(
                            "[regions.{name}] axes {a:?} are turned against the mesh axes: \
                             direction-dependent resistances need axes along x, y or z in \
                             this version (mesh the porous block along its channels)"
                        ),
                    )
                    .key(format!("regions.{name}.axes")),
                ));
            }
        }
    }
    if !any(case) {
        return Ok(());
    }
    // not in this version (porous plan: two temperatures PM2, ranks and GPU
    // PM5)
    let unsupported = if case.energy.is_some() {
        Some("the energy equation (stage PM2 of the porous plan)")
    } else if !case.scalars.is_empty() {
        Some("passive scalars")
    } else if !case.turbulence.is_laminar() {
        Some("a turbulence model")
    } else if !case.mrf.is_empty() || !case.sliding.is_empty() {
        Some("a rotating zone")
    } else if !matches!(case.run.backend.as_str(), "" | "cpu") {
        Some("the GPU (stage PM5 of the porous plan)")
    } else {
        None
    };
    if let Some(what) = unsupported {
        return Err(diag::err(
            Diagnostic::new(
                diag::POROUS_UNSUPPORTED,
                format!("porous regions do not yet work with {what}"),
            )
            .key("regions"),
        ));
    }
    Ok(())
}

/// The orthonormal frame e₁, e₂, e₃ of two axes (e₂ made normal to e₁), as
/// the orthotropic conductivity (L48).
pub fn frame(a: [[f64; 3]; 2]) -> std::result::Result<[[f64; 3]; 3], String> {
    let n = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let l1 = n(a[0]);
    if l1.is_nan() || l1 <= 0.0 {
        return Err("the first axis is zero".into());
    }
    let e1 = a[0].map(|x| x / l1);
    let p = a[1][0] * e1[0] + a[1][1] * e1[1] + a[1][2] * e1[2];
    let t = [
        a[1][0] - p * e1[0],
        a[1][1] - p * e1[1],
        a[1][2] - p * e1[2],
    ];
    let l2 = n(t);
    if l2.is_nan() || l2 <= 1e-12 * n(a[1]).max(f64::MIN_POSITIVE) {
        return Err("the second axis is zero or parallel to the first".into());
    }
    let e2 = t.map(|x| x / l2);
    let e3 = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    Ok([e1, e2, e3])
}

/// Σᵢ vᵢ eᵢ eᵢᵀ in the mesh frame.
pub fn tensor(v: [f64; 3], e: Option<[[f64; 3]; 3]>) -> [[f64; 3]; 3] {
    let e = e.unwrap_or([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let mut t = [[0.0; 3]; 3];
    for (k, ek) in e.iter().enumerate() {
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] += v[k] * ek[i] * ek[j];
            }
        }
    }
    t
}

/// One porous region on the cells of the flow.
#[derive(Debug, Clone)]
pub struct Block {
    /// The region name.
    pub name: String,
    /// Its cells (local indices of the flow).
    pub cells: Vec<usize>,
    /// The resolved values.
    pub values: Resolved,
    /// K⁻¹ in the mesh frame (1/m²).
    pub viscous: [[f64; 3]; 3],
    /// F in the mesh frame (1/m).
    pub inertial: [[f64; 3]; 3],
    /// Whether F is zero (no |u| needed).
    pub darcy_only: bool,
}

/// The porous regions of a flow.
#[derive(Debug, Clone, Default)]
pub struct Porous {
    /// The regions.
    pub blocks: Vec<Block>,
}

/// The porous regions of `case` on `mesh` (`None` without any).
pub fn setup(case: &FlowCase, mesh: &foehncraft_mesh::Mesh) -> Result<Option<Porous>> {
    if !any(case) {
        return Ok(None);
    }
    let mut blocks = Vec::new();
    for (name, r) in &case.regions {
        if r.kind != RegionKind::Porous {
            continue;
        }
        let p = r.porous.as_ref().expect("checked by validate");
        let values = p.resolve().map_err(anyhow::Error::msg)?;
        let e = r.axes.map(frame).transpose().map_err(anyhow::Error::msg)?;
        let cells: Vec<usize> = crate::start::cells_of(case, mesh, name)
            .into_iter()
            .filter(|&c| c < mesh.n_cells)
            .collect();
        blocks.push(Block {
            name: name.clone(),
            cells,
            viscous: tensor(values.viscous, e),
            inertial: tensor(values.inertial, e),
            darcy_only: values.inertial.iter().all(|&f| f == 0.0),
            values,
        });
    }
    Ok(Some(Porous { blocks }))
}

/// Attaches the porous regions of `case` to the flow.
pub fn attach(
    fl: &mut crate::simplec::Flow,
    case: &FlowCase,
    mesh: &foehncraft_mesh::Mesh,
) -> Result<()> {
    fl.porous = setup(case, mesh)?;
    Ok(())
}

impl crate::simplec::Flow {
    /// Whether the momentum rows have a diagonal of their own per component
    /// (symmetry planes, porous regions).
    pub(crate) fn per_component_diagonal(&self) -> bool {
        self.bc
            .mom_kind
            .contains(&foehncraft_kernels::flow::MOM_SYMMETRY)
            || self.porous.is_some()
    }

    /// Adds the Darcy–Forchheimer sink of the porous regions to the
    /// momentum rows just assembled (relaxation `alpha`): the diagonal of
    /// the sink tensor implicitly per component, the rest with the current
    /// velocity, and the smallest diagonal coefficient to D, D^s and D^C.
    pub(crate) fn add_porous(&mut self, alpha: f64, theta: f64) -> Result<()> {
        let Some(p) = self.porous.as_ref() else {
            return Ok(());
        };
        let (f, w) = (self.f, self.w);
        let u: Vec<Vec<f64>> =
            f.u.all()
                .into_iter()
                .map(|b| self.dev.download(b))
                .collect::<Result<_>>()?;
        let mut diag: Vec<Vec<f64>> = w
            .diag
            .all()
            .into_iter()
            .map(|b| self.dev.download(b))
            .collect::<Result<_>>()?;
        let mut rhs: Vec<Vec<f64>> = w
            .rhs
            .all()
            .into_iter()
            .map(|b| self.dev.download(b))
            .collect::<Result<_>>()?;
        let mut d = self.dev.download(w.d)?;
        let mut ds = self.dev.download(w.ds)?;
        let mut dc = self.dev.download(w.dc)?;
        let vol = &self.disc.cell_volume;
        let mu = self.phys.mu;
        for b in &p.blocks {
            for &c in &b.cells {
                let v = vol[c];
                let uc = [u[0][c], u[1][c], u[2][c]];
                let speed = if b.darcy_only {
                    0.0
                } else {
                    (uc[0] * uc[0] + uc[1] * uc[1] + uc[2] * uc[2]).sqrt()
                };
                let q = 0.5 * self.phys.rho_at(c) * speed;
                let a = |i: usize, j: usize| (mu * b.viscous[i][j] + q * b.inertial[i][j]) * v;
                let mut amin = f64::INFINITY;
                for i in 0..3 {
                    // the off-diagonal coupling of axes turned against the
                    // mesh is deferred: its size |a_ij| goes onto the
                    // diagonal and, with the current velocity, onto the
                    // right-hand side, so the row stays diagonally dominant
                    // (exact at convergence; nothing changes for axes
                    // along the mesh)
                    let mut lag = 0.0;
                    let mut r = 0.0;
                    for (j, &uj) in uc.iter().enumerate() {
                        if j != i {
                            let aij = a(i, j);
                            lag += aij.abs();
                            r -= aij * uj;
                        }
                    }
                    let aii = a(i, i) + lag;
                    amin = amin.min(aii);
                    diag[i][c] += aii / alpha;
                    r += (aii / alpha - a(i, i)) * uc[i];
                    rhs[i][c] += r;
                }
                let add = amin / (alpha * v);
                d[c] = 1.0 / (1.0 / d[c] + add);
                dc[c] = 1.0 / (1.0 / dc[c] + add);
                ds[c] = 1.0 / (1.0 / ds[c] + add / theta);
            }
        }
        for (b, x) in w.diag.all().into_iter().zip(&diag) {
            self.dev.upload(b, x)?;
        }
        for (b, x) in w.rhs.all().into_iter().zip(&rhs) {
            self.dev.upload(b, x)?;
        }
        self.dev.upload(w.d, &d)?;
        self.dev.upload(w.ds, &ds)?;
        self.dev.upload(w.dc, &dc)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn honeycomb_values_of_the_brick() {
        let h = PorousSection {
            honeycomb: Some(Honeycomb {
                channel: 6.9e-3,
                pitch: 150e-3 / 11.0,
                length: None,
                blocking: None,
            }),
            ..Default::default()
        };
        let r = h.resolve().unwrap();
        assert!((r.porosity - 121.0 * 6.9f64.powi(2) / 150f64.powi(2)).abs() < 1e-12);
        let along = 28.45 / (r.porosity * 6.9e-3 * 6.9e-3);
        assert!((r.viscous[0] - along).abs() < 1e-9 * along);
        assert_eq!(r.viscous[1], 1e4 * r.viscous[0]);
        assert_eq!(r.inertial, [0.0; 3]);
    }

    #[test]
    fn a_turned_frame_gives_a_symmetric_tensor() {
        let e = frame([[1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]).unwrap();
        let t = tensor([1.0, 100.0, 100.0], Some(e));
        for (i, row) in t.iter().enumerate() {
            for (j, x) in row.iter().enumerate() {
                assert!((x - t[j][i]).abs() < 1e-12);
            }
        }
        // along the first axis the value 1
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let along: f64 = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| [s, s, 0.0][i] * t[i][j] * [s, s, 0.0][j])
                    .sum::<f64>()
            })
            .sum();
        assert!((along - 1.0).abs() < 1e-12);
        assert!(frame([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]).is_err());
    }
}
