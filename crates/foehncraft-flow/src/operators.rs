//! The operator exactness check of a run (D235): before the solve, the
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
use foehncraft_fv::exactness::OperatorErrors;
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

/// The operator exactness check of a mesh (D235).
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
    let e = foehncraft_fv::exactness::operator_errors(mesh, geo, &kind, dims)?;
    let r = OperatorReport {
        cells: mesh.n_cells,
        dimensions: dims,
        linear: part(&e.linear, LINEAR_TOLERANCE, geo),
        gradient: part(&e.gradient, tolerance, geo),
        hessian: part(&e.hessian, tolerance, geo),
        findings: Vec::new(),
    };
    Ok((r, e))
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
    };
    let arrays = crate::bc::arrays(mesh, geo, &case.boundary, case.rho(), at)?;
    check(
        mesh,
        geo,
        &arrays.u_points,
        crate::case::dimension(case) as usize,
        case.solver.operator_tolerance,
    )
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
                "operator check (D235): the least-squares gradient of a linear field is not exact \
                 on {}: degenerate or badly conditioned cells; improve or remesh them",
                cells_text(&r.linear)
            ),
        ));
    }
    let t = &case.turbulence;
    if t.curvature_correction && r.hessian.flagged > 0 {
        out.push(
            diag::Diagnostic::new(
                diag::MESH_OPERATOR_HESSIAN,
                format!(
                    "operator check (D235): the velocity Hessian of the curvature correction \
                     (DS/Dt) has an error above {} on {}: the correction does not converge \
                     with the mesh on those cells",
                    r.hessian.tolerance,
                    cells_text(&r.hessian)
                ),
            )
            .key("turbulence.curvature_correction"),
        );
    }
    if t.qcr != QcrMode::None && r.gradient.flagged > 0 {
        out.push(
            diag::Diagnostic::new(
                diag::MESH_OPERATOR_GRADIENT,
                format!(
                    "operator check (D235): the velocity gradient is first order (error above {} \
                     for a quadratic field) on {}: the divergence of the explicit QCR stress \
                     does not converge with the mesh on those cells",
                    r.gradient.tolerance,
                    cells_text(&r.gradient)
                ),
            )
            .key("turbulence.qcr"),
        );
    }
    out
}

/// The check of a run on its mesh, with the log lines and warnings;
/// `None` (and nothing written) where every operator the run uses is
/// exact, so that the reports of such runs stay as they were.
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
    let ds = diagnostics(&r, case);
    if ds.is_empty() {
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
        assert_eq!(
            d.iter().map(|d| d.code).collect::<Vec<_>>(),
            [diag::MESH_OPERATOR_HESSIAN, diag::MESH_OPERATOR_GRADIENT]
        );
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
