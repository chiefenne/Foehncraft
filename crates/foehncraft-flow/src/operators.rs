//! The operator exactness check of a run (D236): before the solve, the
//! velocity gradient and the gradient of the gradient (the Hessian of the
//! curvature correction) are applied to known linear and quadratic fields
//! on the mesh of the run ([`foehncraft_fv::exactness`]); the cells where
//! they are not exact are counted and listed, worst first, in the log,
//! `flow.md`, `flow.json`, the `case.check` report and the warnings:
//!
//! | Code | When | Meaning |
//! |---|---|---|
//! | `MESH004` | always | the gradient of a linear field is not exact: a degenerate or badly conditioned stencil |
//! | `MESH005` | `turbulence.curvature_correction` | the velocity Hessian (DS/Dt) is inconsistent: the correction does not converge with the mesh there |
//! | `MESH006` | `turbulence.qcr` | the velocity gradient is first order: the divergence of the explicit QCR stress does not converge there |
//!
//! With `MESH005` or `MESH006` the limitation `LIM008` names the feature.
//! The check changes no result.

use crate::case::{FlowCase, Log};
use crate::diag;
use crate::turbulence::QcrMode;
use anyhow::Result;
use foehncraft_fv::exactness::{OperatorErrors, QuadraticGradient};
use foehncraft_fv::problem::BoundaryKind;
use foehncraft_mesh::{Geometry, Mesh};
use serde::Serialize;
use std::collections::BTreeMap;

/// Error of the gradient of a linear field above which a stencil counts as
/// broken (round-off of a sound stencil is far below).
pub const LINEAR_TOLERANCE: f64 = 1e-6;

/// Cells listed per operator, worst first.
pub const LISTED: usize = 50;

/// A cell where an operator is not exact.
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct FlaggedCell {
    /// Cell index of the run (the order of the output fields).
    pub cell: usize,
    /// The dimensionless error.
    pub error: f64,
    /// Cell centroid (m).
    pub centroid: [f64; 3],
}

/// The errors of one operator over the mesh.
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct OperatorPart {
    /// The tolerance a cell is flagged above.
    pub tolerance: f64,
    /// Largest error.
    pub max: f64,
    /// Mean error over the cells.
    pub mean: f64,
    /// Cells above the tolerance.
    pub flagged: usize,
    /// Up to [`LISTED`] of them, worst first.
    pub worst: Vec<FlaggedCell>,
}

/// The operator exactness check of a mesh (D236).
#[derive(Debug, Clone, Serialize, schemars::JsonSchema)]
pub struct OperatorReport {
    /// Cells checked.
    pub cells: usize,
    /// Coordinates of the test fields (2: x and y of a two-dimensional case).
    pub dimensions: usize,
    /// The velocity gradient on linear fields, |∇ₕφ − ∇φ| L.
    pub linear: OperatorPart,
    /// The velocity gradient on quadratic fields, |∇ₕq − ∇q|/(h |H|):
    /// zero where the gradient is second order.
    pub gradient: OperatorPart,
    /// The gradient of the gradient on quadratic fields, |Hₕ − H|/|H|.
    pub hessian: OperatorPart,
    /// With the QCR (step 2 of D236): the error of the velocity gradient
    /// the stress takes, the quadratic least-squares fit on the cells
    /// flagged in `gradient` (the face-neighbour gradient where the fit is
    /// rank deficient): the cells still above the tolerance.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quadratic: Option<OperatorPart>,
    /// With the curvature correction (step 3 of D236): the error of the
    /// gradient of S_ij DS_ij/Dt takes, a quadratic fit of S_ij on the
    /// cells flagged in `hessian` (the gradient of the gradient where the
    /// fit or the velocity fit in its stencil is rank deficient): the cells
    /// still above the tolerance.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quadratic_hessian: Option<OperatorPart>,
    /// What the flagged cells mean for the run (the messages of the
    /// diagnostics, filled in a run).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub findings: Vec<String>,
}

fn part(v: &[f64], tolerance: f64, geo: &Geometry) -> OperatorPart {
    let mut idx: Vec<usize> = (0..v.len()).filter(|&c| v[c] > tolerance).collect();
    let flagged = idx.len();
    idx.sort_by(|&a, &b| v[b].total_cmp(&v[a]).then(a.cmp(&b)));
    idx.truncate(LISTED);
    OperatorPart {
        tolerance,
        max: v.iter().copied().fold(0.0, f64::max),
        mean: foehncraft_core::repro::pairwise_sum(v) / v.len().max(1) as f64,
        flagged,
        worst: idx
            .into_iter()
            .map(|c| FlaggedCell {
                cell: c,
                error: v[c],
                centroid: geo.cell_centroid[c],
            })
            .collect(),
    }
}

/// The check on `mesh` with the velocity stencil points `u_points` (as the
/// velocity gradient of the run takes them), `dims` coordinates and the
/// tolerance of `solver.operator_tolerance`.
pub fn check(
    mesh: &Mesh,
    geo: &Geometry,
    u_points: &BTreeMap<String, BoundaryKind>,
    dims: usize,
    tolerance: f64,
) -> Result<OperatorReport> {
    Ok(check_with_errors(mesh, geo, u_points, dims, tolerance)?.0)
}

/// [`check`] with the errors of every cell.
pub fn check_with_errors(
    mesh: &Mesh,
    geo: &Geometry,
    u_points: &BTreeMap<String, BoundaryKind>,
    dims: usize,
    tolerance: f64,
) -> Result<(OperatorReport, OperatorErrors)> {
    let kind = kinds(mesh, u_points);
    let e = foehncraft_fv::exactness::operator_errors(mesh, geo, &kind, dims)?;
    let r = OperatorReport {
        cells: mesh.n_cells,
        dimensions: dims,
        linear: part(&e.linear, LINEAR_TOLERANCE, geo),
        gradient: part(&e.gradient, tolerance, geo),
        hessian: part(&e.hessian, tolerance, geo),
        quadratic: None,
        quadratic_hessian: None,
        findings: Vec::new(),
    };
    Ok((r, e))
}

/// The boundary kind of every boundary face from the velocity stencil
/// points.
fn kinds(mesh: &Mesh, u_points: &BTreeMap<String, BoundaryKind>) -> Vec<BoundaryKind> {
    let ni = mesh.n_interior_faces();
    let mut kind = vec![BoundaryKind::FixedValue; mesh.n_faces() - ni];
    for p in &mesh.boundaries {
        let k = u_points
            .get(&p.name)
            .copied()
            .unwrap_or(BoundaryKind::FixedValue);
        for f in p.range() {
            kind[f - ni] = k;
        }
    }
    kind
}

/// The quadratic least-squares fits of D236 a run takes: the velocity
/// gradient on the cells whose face-neighbour gradient is first order
/// (`gradient` error above `solver.operator_tolerance`; with the QCR or
/// the curvature correction), and with the curvature correction the
/// gradient of S_ij on the cells whose gradient of the gradient is not
/// exact (`hessian` error above it), fitted over interior cells only (no
/// boundary values of S_ij are needed).
pub fn quadratic_fits(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    u_points: &BTreeMap<String, BoundaryKind>,
) -> Result<(QuadraticGradient, Option<QuadraticGradient>)> {
    let dims = crate::case::dimension(case) as usize;
    let kind = kinds(mesh, u_points);
    let e = foehncraft_fv::exactness::operator_errors(mesh, geo, &kind, dims)?;
    let tol = case.solver.operator_tolerance;
    let (qcr, curvature) = uses(case);
    Ok((
        velocity_fit(mesh, geo, &kind, dims, &e, tol, qcr || curvature),
        curvature.then(|| strain_fit(mesh, geo, dims, &e, tol)),
    ))
}

/// Which operators the run differentiates once more: (the QCR, the
/// curvature correction).
fn uses(case: &FlowCase) -> (bool, bool) {
    let t = &case.turbulence;
    (t.qcr != QcrMode::None, t.curvature_correction)
}

fn velocity_fit(
    mesh: &Mesh,
    geo: &Geometry,
    kind: &[BoundaryKind],
    dims: usize,
    e: &OperatorErrors,
    tolerance: f64,
    on: bool,
) -> QuadraticGradient {
    let cells: Vec<usize> = (0..e.gradient.len())
        .filter(|&c| on && e.gradient[c] > tolerance)
        .collect();
    foehncraft_fv::exactness::quadratic_gradient(mesh, geo, kind, dims, &cells)
}

/// Boundary kinds that leave every boundary face out of a fit.
fn interior_only(mesh: &Mesh) -> Vec<BoundaryKind> {
    vec![BoundaryKind::ZeroGradient; mesh.n_faces() - mesh.n_interior_faces()]
}

fn strain_fit(
    mesh: &Mesh,
    geo: &Geometry,
    dims: usize,
    e: &OperatorErrors,
    tolerance: f64,
) -> QuadraticGradient {
    let cells: Vec<usize> = (0..e.hessian.len())
        .filter(|&c| e.hessian[c] > tolerance)
        .collect();
    foehncraft_fv::exactness::quadratic_gradient(mesh, geo, &interior_only(mesh), dims, &cells)
}

/// Fills `quadratic` (with the QCR) and `quadratic_hessian` (with the
/// curvature correction) of `r`: the errors of the operators the run
/// takes with the quadratic fits on the flagged cells. Nothing without
/// either feature.
fn with_quadratic(
    r: &mut OperatorReport,
    e: &OperatorErrors,
    mesh: &Mesh,
    geo: &Geometry,
    u_points: &BTreeMap<String, BoundaryKind>,
    tolerance: f64,
    (qcr, curvature): (bool, bool),
) {
    let kind = kinds(mesh, u_points);
    let dims = r.dimensions;
    let q = velocity_fit(mesh, geo, &kind, dims, e, tolerance, qcr || curvature);
    let qe = foehncraft_fv::exactness::quadratic_gradient_errors(mesh, geo, &kind, &q);
    let mut g = e.gradient.clone();
    for (&c, &x) in q.cells.iter().zip(&qe) {
        g[c as usize] = x;
    }
    if qcr {
        r.quadratic = Some(part(&g, tolerance, geo));
    }
    if curvature {
        // the strain fit is exact for the linear strain of a quadratic
        // field; that strain is exact where the velocity gradient is (the
        // largest gradient error over the stencil counts as well)
        let qs = strain_fit(mesh, geo, dims, e, tolerance);
        let zero = interior_only(mesh);
        let se = foehncraft_fv::exactness::quadratic_linear_errors(mesh, geo, &zero, &qs);
        let mut v = e.hessian.clone();
        for (t, (&c, &x)) in qs.cells.iter().zip(&se).enumerate() {
            let (a, b) = (qs.ptr[t] as usize, qs.ptr[t + 1] as usize);
            let worst = qs.idx[a..b]
                .iter()
                .map(|&p| g[p as usize])
                .fold(g[c as usize], f64::max);
            v[c as usize] = x.max(worst);
        }
        r.quadratic_hessian = Some(part(&v, tolerance, geo));
    }
}

/// The file of a run with the errors of every cell (written when the check
/// reports something).
pub const CSV: &str = "operators.csv";

/// Writes `cell,x,y,z,linear,gradient,hessian`, one row per cell.
pub fn write_csv(path: &std::path::Path, geo: &Geometry, e: &OperatorErrors) -> Result<()> {
    use std::fmt::Write as _;
    let mut s = String::from("cell,x,y,z,linear,gradient,hessian\n");
    for c in 0..e.linear.len() {
        let x = geo.cell_centroid[c];
        let _ = writeln!(
            s,
            "{c},{:.9e},{:.9e},{:.9e},{:.6e},{:.6e},{:.6e}",
            x[0], x[1], x[2], e.linear[c], e.gradient[c], e.hessian[c]
        );
    }
    std::fs::write(path, s)?;
    Ok(())
}

/// The check of the first mesh of `case` (the operators of the velocity of
/// the case's boundary conditions).
pub fn check_case(case: &FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<OperatorReport> {
    let values = case.compiled_boundaries()?;
    let at = crate::bc::At {
        analytic: None,
        values: &values,
        t: 0.0,
        turning: &case.turning(),
        hydrostatic: crate::buoyancy::hydrostatic(case),
    };
    let arrays = crate::bc::arrays(mesh, geo, &case.boundary, case.rho(), at)?;
    let tol = case.solver.operator_tolerance;
    let (mut r, e) = check_with_errors(
        mesh,
        geo,
        &arrays.u_points,
        crate::case::dimension(case) as usize,
        tol,
    )?;
    with_quadratic(&mut r, &e, mesh, geo, &arrays.u_points, tol, uses(case));
    Ok(r)
}

fn cells_text(p: &OperatorPart) -> String {
    let list: Vec<String> = p.worst.iter().take(5).map(|c| c.cell.to_string()).collect();
    format!(
        "{} cells (largest {:.3}, worst cells {}{})",
        p.flagged,
        p.max,
        list.join(", "),
        if p.flagged > 5 { ", …" } else { "" }
    )
}

/// The diagnostics of a check for `case`: `MESH004` with a broken stencil,
/// `MESH005` and `MESH006` where a feature in use needs the operator.
pub fn diagnostics(r: &OperatorReport, case: &FlowCase) -> Vec<diag::Diagnostic> {
    let mut out = Vec::new();
    if r.linear.flagged > 0 {
        out.push(diag::Diagnostic::new(
            diag::MESH_OPERATOR_LINEAR,
            format!(
                "operator check (D236): the least-squares gradient of a linear field is not exact \
                 on {}: degenerate or badly conditioned cells; improve or remesh them",
                cells_text(&r.linear)
            ),
        ));
    }
    let t = &case.turbulence;
    let h = r.quadratic_hessian.as_ref().unwrap_or(&r.hessian);
    if t.curvature_correction && h.flagged > 0 {
        out.push(
            diag::Diagnostic::new(
                diag::MESH_OPERATOR_HESSIAN,
                format!(
                    "operator check (D236): the velocity Hessian of the curvature correction \
                     (DS/Dt) has an error above {} for a quadratic field{} on {}: the correction \
                     does not converge with the mesh on those cells",
                    h.tolerance,
                    if r.quadratic_hessian.is_some() {
                        "; the quadratic fit is rank deficient there"
                    } else {
                        ""
                    },
                    cells_text(h)
                ),
            )
            .key("turbulence.curvature_correction"),
        );
    }
    let g = r.quadratic.as_ref().unwrap_or(&r.gradient);
    if t.qcr != QcrMode::None && g.flagged > 0 {
        out.push(
            diag::Diagnostic::new(
                diag::MESH_OPERATOR_GRADIENT,
                format!(
                    "operator check (D236): the velocity gradient of the QCR stress is first order \
                     (error above {} for a quadratic field{}) on {}: the divergence of the explicit \
                     QCR stress does not converge with the mesh on those cells",
                    g.tolerance,
                    if r.quadratic.is_some() {
                        "; the quadratic fit is rank deficient there"
                    } else {
                        ""
                    },
                    cells_text(g)
                ),
            )
            .key("turbulence.qcr"),
        );
    }
    out
}

/// The check of a run on its mesh, with the log lines and warnings;
/// `None` (and nothing written) where every operator the run uses is exact
/// and no cell takes the quadratic fit, so that the reports of such runs
/// stay as they were.
#[allow(clippy::too_many_arguments)]
pub(crate) fn report<F: FnMut(&str)>(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    u_points: &BTreeMap<String, BoundaryKind>,
    out: &std::path::Path,
    lg: &mut Log<F>,
    warnings: &mut Vec<diag::Diagnostic>,
) -> Option<OperatorReport> {
    let (mut r, e) = check_with_errors(
        mesh,
        geo,
        u_points,
        crate::case::dimension(case) as usize,
        case.solver.operator_tolerance,
    )
    .ok()?;
    with_quadratic(
        &mut r,
        &e,
        mesh,
        geo,
        u_points,
        case.solver.operator_tolerance,
        uses(case),
    );
    if let Some(q) = &r.quadratic {
        let treated = r.gradient.flagged - q.flagged;
        if treated > 0 {
            lg.line(&format!(
                "QCR: quadratic least-squares velocity gradient on {treated} cells, skewness-corrected \
                 face stress (D236)"
            ));
        }
    }
    if let Some(q) = &r.quadratic_hessian {
        let treated = r.hessian.flagged - q.flagged;
        if treated > 0 {
            lg.line(&format!(
                "curvature correction: S_ij from the quadratic velocity fit on {} cells, DS_ij/Dt \
                 from a quadratic least-squares fit of S_ij on {treated} cells (D236)",
                r.gradient.flagged
            ));
        }
    }
    let ds = diagnostics(&r, case);
    // cells the quadratic fit treats are reported too (flow.json,
    // operators.csv), without a diagnostic
    let treated = (r.quadratic.is_some() && r.gradient.flagged > 0)
        || (r.quadratic_hessian.is_some() && r.hessian.flagged > 0);
    if ds.is_empty() && !treated {
        return None;
    }
    // every cell with its errors, for a plot or an agent
    if write_csv(&out.join(CSV), geo, &e).is_ok() {
        lg.line(&format!(
            "operator check: the errors of every cell in {CSV}"
        ));
    }
    for d in ds {
        lg.line(&d.message);
        r.findings.push(d.message.clone());
        if !warnings.iter().any(|w| w.code == d.code) {
            warnings.push(d);
        }
    }
    Some(r)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case(mesh: &str, extra: &str) -> FlowCase {
        let text = format!(
            r#"
[mesh.generate]
type = "{mesh}"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [8, 8, 1]
[fluid]
viscosity = 1e-3
[turbulence]
model = "k-omega-sst"
{extra}
[boundary.xmin]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
[boundary.xmax]
type = "static-pressure-outlet"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
"#
        );
        FlowCase::from_toml(&text, std::path::Path::new(".")).unwrap()
    }

    fn run(mesh: &str, extra: &str) -> (OperatorReport, Vec<diag::Diagnostic>) {
        let c = case(mesh, extra);
        let (m, g, _) = crate::case::mesh_for(&c, None).unwrap();
        let r = check_case(&c, &m, &g).unwrap();
        let d = diagnostics(&r, &c);
        (r, d)
    }

    #[test]
    fn regular_hexahedra_flag_only_the_boundary_cells() {
        let (r, d) = run("hex", "curvature_correction = true\nqcr = \"qcr2000\"");
        assert_eq!(r.linear.flagged, 0);
        assert!(r.gradient.flagged > 0 && r.gradient.flagged < r.cells);
        // inside: exact
        let inner = r.cells - 2 * 8 - 2 * 8 + 4;
        assert!(r.gradient.flagged <= r.cells - inner, "{:?}", r.gradient);
        assert!(r.hessian.flagged > r.gradient.flagged, "{:?}", r.hessian);
        // the quadratic fit on the flagged cells is exact for the QCR
        // stress and the curvature correction alike: nothing reported
        let q = r.quadratic.as_ref().unwrap();
        assert_eq!(q.flagged, 0, "{q:?}");
        let q = r.quadratic_hessian.as_ref().unwrap();
        assert_eq!(q.flagged, 0, "{q:?}");
        assert!(d.is_empty(), "{d:?}");
    }

    #[test]
    fn quadratic_fit_is_exact_on_tetrahedra() {
        let (r, d) = run("tet", "curvature_correction = true\nqcr = \"qcr2000\"");
        assert_eq!(
            r.quadratic.as_ref().unwrap().flagged,
            0,
            "{:?}",
            r.quadratic
        );
        assert_eq!(
            r.quadratic_hessian.as_ref().unwrap().flagged,
            0,
            "{:?}",
            r.quadratic_hessian
        );
        assert!(d.is_empty(), "{d:?}");
        // without the features no fit
        let (r, _) = run("tet", "");
        assert!(r.quadratic.is_none() && r.quadratic_hessian.is_none());
    }

    #[test]
    fn tetrahedra_flag_every_cell() {
        let (r, _) = run("tet", "");
        assert_eq!(r.linear.flagged, 0);
        assert_eq!(r.gradient.flagged, r.cells);
        assert!(r.hessian.flagged > r.cells / 2, "{:?}", r.hessian);
        assert!(r.worst_sorted());
    }

    #[test]
    fn nothing_is_reported_without_the_features() {
        let (_, d) = run("tet", "");
        assert!(d.is_empty(), "{d:?}");
    }

    impl OperatorReport {
        fn worst_sorted(&self) -> bool {
            [&self.gradient, &self.hessian]
                .iter()
                .all(|p| p.worst.windows(2).all(|w| w[0].error >= w[1].error))
        }
    }
}
