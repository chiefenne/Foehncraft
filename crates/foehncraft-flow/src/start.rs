//! Start values per region and per fluid part (D238).
//!
//! In a case of several domains every region gets its own start values,
//! so that no domain starts far from its operating state:
//!
//! - `[regions.<name>.initial]` sets `velocity`, `pressure`, `temperature`,
//!   `k`, `omega` and `epsilon` of a fluid region (numbers or expressions of
//!   x, y, z, D53); a solid region takes `temperature` only (REG007). The
//!   region's table comes before the case-wide `[initial]` and
//!   `[energy] initial`. `temperature = "inlet"` starts a fluid region at
//!   the mean temperature of the inlets of its fluid part, weighted by their
//!   inflow.
//! - Without a start pressure, every fluid part starts from its own
//!   pressure boundaries: at their value when they all have one, else from
//!   the solution of a Laplace equation with those values on them and zero
//!   normal gradient elsewhere, so the start field runs smoothly from one
//!   boundary to the other (a total-pressure boundary counts with its total
//!   pressure: the fluid starts at rest). With absolute pressure values
//!   (D239) each part's own outlet is its operating pressure, so its gauge
//!   start is 0 there. Not with `[initial] pressure` (it stays the start of
//!   every cell), a reference cell (a closed case starts at its reference),
//!   or gravity (the hydrostatic start per part comes with MF-B).
//! - A restart replaces all of them; the log (REG010) and `flow.json`
//!   (`start_values`) say what each part and region started from.
//!
//! Volume selections inside a region (`[[initial.patch]]`) come later
//! (D238, with D73).

use crate::bc::{self, BoundaryCondition};
use crate::case::FlowCase;
use crate::diag::{self, Diagnostic};
use crate::expr::{Env, Scalar, Vars, Vector};
use crate::materials::RegionKind;
use anyhow::Result;
use foehncraft_mesh::geometry::{dot, norm, sub};
use foehncraft_mesh::{Geometry, Mesh};
use serde::{Deserialize, Serialize};

/// `[regions.<name>.initial]`: the start values of one region.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RegionInitial {
    /// Start velocity (m/s): three numbers or expressions of x, y, z.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub velocity: Option<Vector>,
    /// Start pressure (Pa), in the convention of `[fluid] pressure_values`
    /// (absolute by default): a number or an expression of x, y, z.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure: Option<Scalar>,
    /// Start temperature (K): a number, an expression of x, y, z, or for a
    /// fluid region `"inlet"`, the inflow-weighted mean temperature of the
    /// inlets of its fluid part. The only key of a solid region.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature: Option<Scalar>,
    /// Start k (m²/s²), with `omega` or `epsilon`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub k: Option<Scalar>,
    /// Start ω (1/s), with `k`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub omega: Option<Scalar>,
    /// Start ε (m²/s³), with `k`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub epsilon: Option<Scalar>,
}

impl RegionInitial {
    fn flow_keys(&self) -> Vec<&'static str> {
        let mut v = Vec::new();
        if self.velocity.is_some() {
            v.push("velocity");
        }
        if self.pressure.is_some() {
            v.push("pressure");
        }
        for (k, s) in [
            ("k", &self.k),
            ("omega", &self.omega),
            ("epsilon", &self.epsilon),
        ] {
            if s.is_some() {
                v.push(k);
            }
        }
        v
    }

    /// True for `temperature = "inlet"`.
    pub fn temperature_from_inlet(&self) -> bool {
        matches!(&self.temperature, Some(Scalar::Expr(s)) if s.trim() == "inlet")
    }
}

/// One line of the start report: what started where, from what.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct StartValue {
    /// `region <name>` or `part <n>`.
    pub scope: String,
    /// `velocity`, `pressure`, `temperature`, `turbulence`.
    pub quantity: String,
    /// The value (a number, a range or an expression).
    pub value: String,
    /// Where it comes from: `given`, `boundary <name>`, `Laplace solve
    /// between N boundaries`, `inlet`, `case-wide`.
    pub source: String,
}

fn err(code: &'static str, key: String, msg: String) -> anyhow::Error {
    diag::err(Diagnostic::new(code, msg).key(key))
}

/// The checks of the region start tables (REG007, REG008); run with the
/// case.
pub fn validate(case: &FlowCase) -> Result<()> {
    for (name, r) in &case.regions {
        let Some(i) = &r.initial else { continue };
        if r.kind != RegionKind::Fluid {
            let keys = i.flow_keys();
            if !keys.is_empty() || i.temperature_from_inlet() {
                return Err(err(
                    diag::REG_INITIAL_KEY,
                    format!("regions.{name}.initial"),
                    format!(
                        "[regions.{name}.initial] is a solid region's start table: it takes \
                         `temperature` (a number or an expression) only, not {}",
                        if keys.is_empty() {
                            "temperature = \"inlet\"".to_string()
                        } else {
                            keys.join(", ")
                        }
                    ),
                ));
            }
        }
        if i.k.is_some() != (i.omega.is_some() || i.epsilon.is_some())
            || (i.omega.is_some() && i.epsilon.is_some())
        {
            return Err(err(
                diag::REG_INITIAL_KEY,
                format!("regions.{name}.initial.k"),
                format!(
                    "[regions.{name}.initial] gives turbulence start values as `k` with one of \
                     `omega` or `epsilon`"
                ),
            ));
        }
        if i.temperature.is_some() && case.energy.is_none() {
            return Err(err(
                diag::REG_INITIAL_KEY,
                format!("regions.{name}.initial.temperature"),
                format!(
                    "[regions.{name}.initial] temperature needs the energy equation ([energy])"
                ),
            ));
        }
        let pos = |s: &Option<Scalar>, k: &str| -> Result<()> {
            if let Some(v) = s.as_ref().and_then(Scalar::number)
                && !(v > 0.0 && v.is_finite())
            {
                return Err(err(
                    diag::REG_INITIAL_VALUE,
                    format!("regions.{name}.initial.{k}"),
                    format!("[regions.{name}.initial] {k} must be > 0, is {v}"),
                ));
            }
            Ok(())
        };
        pos(&i.temperature, "temperature")?;
        pos(&i.k, "k")?;
        pos(&i.omega, "omega")?;
        pos(&i.epsilon, "epsilon")?;
    }
    Ok(())
}

/// The fluid region of every cell zone of `mesh` (`None`: no region).
fn zone_region<'a>(case: &'a FlowCase, mesh: &Mesh) -> Vec<Option<&'a str>> {
    mesh.zones
        .iter()
        .map(|z| {
            case.regions
                .iter()
                .find(|(_, r)| r.zones.contains(z))
                .map(|(n, _)| n.as_str())
        })
        .collect()
}

/// The cells of region `name` on `mesh`.
pub fn cells_of(case: &FlowCase, mesh: &Mesh, name: &str) -> Vec<usize> {
    let zr = zone_region(case, mesh);
    (0..mesh.n_cells)
        .filter(|&c| zr[mesh.cell_zone[c] as usize] == Some(name))
        .collect()
}

fn describe(s: &Scalar) -> String {
    match s {
        Scalar::Number(v) => format!("{v}"),
        Scalar::Expr(e) => format!("\"{e}\""),
    }
}

/// The start velocity and pressure of the flow cells: the region tables
/// and the start pressure of every fluid part from its own boundaries.
/// `u` and `p` hold the case-wide start values and are changed in place;
/// the result is `None` when nothing changed (the start of the case stays
/// bit for bit as before).
#[allow(clippy::too_many_arguments)]
pub fn flow_fields(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    env: &Env,
    values: &bc::Compiled,
    u: &mut foehncraft_fv::disc::Soa3,
    p: &mut [f64],
    reference: bool,
) -> Result<Option<Vec<StartValue>>> {
    let nc = mesh.n_cells;
    let mut rep = Vec::new();
    let mut fixed_p = vec![false; nc];
    // the region tables
    for (name, r) in &case.regions {
        let Some(i) = &r.initial else { continue };
        if r.kind != RegionKind::Fluid {
            continue;
        }
        let cells = cells_of(case, mesh, name);
        if cells.is_empty() {
            continue;
        }
        if let Some(v) = &i.velocity {
            let key = format!("regions.{name}.initial.velocity");
            let vv = crate::expr::vector_value(v, env, Vars::SPACE, &key)?;
            for &c in &cells {
                let x = geo.cell_centroid[c];
                u.set(c, [vv[0].at(x, 0.0), vv[1].at(x, 0.0), vv[2].at(x, 0.0)]);
            }
            rep.push(StartValue {
                scope: format!("region {name}"),
                quantity: "velocity".into(),
                value: format!(
                    "({}, {}, {})",
                    describe(&v[0]),
                    describe(&v[1]),
                    describe(&v[2])
                ),
                source: "given".into(),
            });
        }
        if let Some(s) = &i.pressure {
            let key = format!("regions.{name}.initial.pressure");
            let pv = crate::expr::value(s, env, Vars::SPACE, &key)?;
            for &c in &cells {
                p[c] = pv.at(geo.cell_centroid[c], 0.0);
                fixed_p[c] = true;
            }
            rep.push(StartValue {
                scope: format!("region {name}"),
                quantity: if case.pressure_frame.as_deref().is_some_and(|f| f.absolute()) {
                    "pressure (p - p_op of its part)".into()
                } else {
                    "pressure".into()
                },
                value: describe(s),
                source: "given".into(),
            });
        }
    }
    // the start pressure of every part from its own boundaries: not with a
    // case-wide start pressure or a reference cell, and not yet with
    // gravity (its hydrostatic start comes with MF-B)
    let given = case
        .pressure_frame
        .as_deref()
        .is_some_and(|f| f.initial_given);
    if !given && !reference && case.gravity.is_none() {
        let comp = crate::pressure_values::components(case, mesh);
        // the report speaks the file's convention: absolute values with
        // the operating pressure of the part added back
        let frame = case.pressure_frame.as_deref().filter(|f| f.absolute());
        let offsets = frame.and_then(|f| f.cell_offsets(case, mesh));
        let quantity = if frame.is_some() {
            "pressure (absolute)"
        } else {
            "pressure"
        };
        let first_cell = |k: usize| (0..nc).find(|&c| comp.of[c] == k as u32);
        let shown = |k: usize, v: f64| -> f64 {
            match (&offsets, first_cell(k)) {
                (Some(o), Some(c)) => v + o[c],
                _ => v,
            }
        };
        let ni = mesh.n_interior_faces();
        // the Dirichlet faces of every part: (face, value, boundary)
        let mut faces: Vec<Vec<(usize, f64)>> = vec![Vec::new(); comp.n];
        let mut names: Vec<Vec<String>> = vec![Vec::new(); comp.n];
        for b in &mesh.boundaries {
            let Some(bcond) = case.boundary.get(&b.name) else {
                continue;
            };
            if !bcond.is_pressure() || matches!(bcond, BoundaryCondition::Periodic { .. }) {
                continue;
            }
            let vals = values.get(&b.name);
            for f in b.range() {
                let k = comp.of[mesh.owner[f] as usize];
                if k == u32::MAX {
                    continue;
                }
                let v = bcond.pressure_at(vals, geo.face_centroid[f], 0.0);
                faces[k as usize].push((f, v));
                if !names[k as usize].contains(&b.name) {
                    names[k as usize].push(b.name.clone());
                }
            }
        }
        let mut laplace = vec![false; comp.n];
        let mut uniform: Vec<Option<f64>> = vec![None; comp.n];
        for k in 0..comp.n {
            let Some(&(_, v0)) = faces[k].first() else {
                continue;
            };
            if faces[k].iter().all(|&(_, v)| v.to_bits() == v0.to_bits()) {
                uniform[k] = Some(v0);
            } else {
                laplace[k] = true;
            }
        }
        // a part at one value: every cell there (no change at 0 Pa, the
        // uniform start of before)
        for c in 0..nc {
            let k = comp.of[c];
            if k == u32::MAX || fixed_p[c] {
                continue;
            }
            if let Some(v) = uniform[k as usize]
                && v.to_bits() != p[c].to_bits()
            {
                p[c] = v;
            }
        }
        for k in 0..comp.n {
            if let Some(v) = uniform[k]
                && v != 0.0
            {
                rep.push(StartValue {
                    scope: format!("part {}", k + 1),
                    quantity: quantity.into(),
                    value: format!("{}", shown(k, v)),
                    source: format!("boundary {}", names[k].join(", ")),
                });
            }
        }
        if laplace.iter().any(|&l| l) {
            let free: Vec<bool> = (0..nc)
                .map(|c| {
                    let k = comp.of[c];
                    k != u32::MAX && laplace[k as usize] && !fixed_p[c]
                })
                .collect();
            let mut dir = vec![None; mesh.n_faces() - ni];
            for list in &faces {
                for &(f, v) in list {
                    dir[f - ni] = Some(v);
                }
            }
            let (its, res) = laplace_solve(mesh, geo, &free, &dir, p);
            for k in (0..comp.n).filter(|&k| laplace[k]) {
                let (lo, hi) = (0..nc)
                    .filter(|&c| comp.of[c] == k as u32)
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), c| {
                        (a.min(p[c]), b.max(p[c]))
                    });
                rep.push(StartValue {
                    scope: format!("part {}", k + 1),
                    quantity: quantity.into(),
                    value: format!("{} to {}", shown(k, lo), shown(k, hi)),
                    source: format!(
                        "Laplace solve between {} boundaries ({}; {its} CG iterations, \
                         residual {res:.1e})",
                        names[k].len(),
                        names[k].join(", ")
                    ),
                });
            }
        }
    }
    Ok((!rep.is_empty()).then_some(rep))
}

/// ∇²p = 0 on the `free` cells with the values `dir` on boundary faces
/// (zero normal gradient elsewhere and the other cells held at `p`):
/// two-point fluxes |S|/|d|, Jacobi-preconditioned CG to 10⁻¹² of the
/// right-hand side. Returns the iterations and the relative residual.
fn laplace_solve(
    mesh: &Mesh,
    geo: &Geometry,
    free: &[bool],
    dir: &[Option<f64>],
    p: &mut [f64],
) -> (usize, f64) {
    let nc = mesh.n_cells;
    let ni = mesh.n_interior_faces();
    let w = |f: usize, a: [f64; 3], b: [f64; 3]| {
        let s = geo.face_area[f];
        let d = sub(b, a);
        let dd = dot(d, d).max(1e-300);
        // the orthogonal part |S|²/(S·d), bounded for skewed faces
        let sd = dot(s, d).abs().max(1e-3 * norm(s) * dd.sqrt());
        dot(s, s) / sd
    };
    let mut diag = vec![0.0; nc];
    let mut rhs = vec![0.0; nc];
    let mut off: Vec<(usize, usize, f64)> = Vec::new();
    for f in 0..ni {
        let (o, n) = (mesh.owner[f] as usize, mesh.neighbour[f] as usize);
        let a = w(f, geo.cell_centroid[o], geo.cell_centroid[n]);
        match (free[o], free[n]) {
            (true, true) => {
                diag[o] += a;
                diag[n] += a;
                off.push((o, n, a));
            }
            (true, false) => {
                diag[o] += a;
                rhs[o] += a * p[n];
            }
            (false, true) => {
                diag[n] += a;
                rhs[n] += a * p[o];
            }
            _ => {}
        }
    }
    for (fb, d) in dir.iter().enumerate() {
        let f = ni + fb;
        let o = mesh.owner[f] as usize;
        if let Some(v) = d
            && free[o]
        {
            let a = w(f, geo.cell_centroid[o], geo.face_centroid[f]);
            diag[o] += a;
            rhs[o] += a * v;
        }
    }
    let idx: Vec<usize> = (0..nc).filter(|&c| free[c] && diag[c] > 0.0).collect();
    let apply = |x: &[f64], y: &mut [f64]| {
        for &c in &idx {
            y[c] = diag[c] * x[c];
        }
        for &(o, n, a) in &off {
            y[o] -= a * x[n];
            y[n] -= a * x[o];
        }
    };
    let mut x: Vec<f64> = p.to_vec();
    let mut ax = vec![0.0; nc];
    apply(&x, &mut ax);
    let mut r = vec![0.0; nc];
    for &c in &idx {
        r[c] = rhs[c] - ax[c];
    }
    let norm_b = idx
        .iter()
        .map(|&c| rhs[c] * rhs[c])
        .sum::<f64>()
        .sqrt()
        .max(1e-300);
    let mut z = vec![0.0; nc];
    for &c in &idx {
        z[c] = r[c] / diag[c];
    }
    let mut d = z.clone();
    let mut rz: f64 = idx.iter().map(|&c| r[c] * z[c]).sum();
    let mut q = vec![0.0; nc];
    let mut its = 0;
    let mut res = idx.iter().map(|&c| r[c] * r[c]).sum::<f64>().sqrt() / norm_b;
    while res > 1e-12 && its < 20_000 {
        apply(&d, &mut q);
        let dq: f64 = idx.iter().map(|&c| d[c] * q[c]).sum();
        if dq <= 0.0 {
            break;
        }
        let alpha = rz / dq;
        for &c in &idx {
            x[c] += alpha * d[c];
            r[c] -= alpha * q[c];
            z[c] = r[c] / diag[c];
        }
        let rz1: f64 = idx.iter().map(|&c| r[c] * z[c]).sum();
        let beta = rz1 / rz;
        rz = rz1;
        for &c in &idx {
            d[c] = z[c] + beta * d[c];
        }
        its += 1;
        res = idx.iter().map(|&c| r[c] * r[c]).sum::<f64>().sqrt() / norm_b;
    }
    for &c in &idx {
        p[c] = x[c];
    }
    (its, res)
}

/// The start temperature per cell of the energy mesh `mesh` (the full mesh
/// of a conjugate case): the region tables over `[energy] initial` in
/// `phi`. Returns the report lines.
pub fn temperature(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    env: &Env,
    phi: &mut [f64],
) -> Result<Vec<StartValue>> {
    let mut rep = Vec::new();
    let nc = phi.len().min(mesh.n_cells);
    for (name, r) in &case.regions {
        let Some(t) = r.initial.as_ref().and_then(|i| i.temperature.as_ref()) else {
            continue;
        };
        let cells: Vec<usize> = cells_of(case, mesh, name)
            .into_iter()
            .filter(|&c| c < nc)
            .collect();
        if cells.is_empty() {
            continue;
        }
        let from_inlet = r
            .initial
            .as_ref()
            .is_some_and(RegionInitial::temperature_from_inlet);
        if from_inlet {
            let (t_in, inlets) = inlet_temperature(case, mesh, geo, env, &cells)?;
            for &c in &cells {
                phi[c] = t_in;
            }
            rep.push(StartValue {
                scope: format!("region {name}"),
                quantity: "temperature".into(),
                value: format!("{t_in}"),
                source: format!("inlet ({})", inlets.join(", ")),
            });
        } else {
            let key = format!("regions.{name}.initial.temperature");
            let v = crate::expr::value(t, env, Vars::SPACE, &key)?;
            for &c in &cells {
                let x = v.at(geo.cell_centroid[c], 0.0);
                if !(x > 0.0 && x.is_finite()) {
                    return Err(err(
                        diag::REG_INITIAL_VALUE,
                        key,
                        format!(
                            "[regions.{name}.initial] temperature is {x} K at the cell centroid \
                             {:?}; a temperature must be > 0 K",
                            geo.cell_centroid[c]
                        ),
                    ));
                }
                phi[c] = x;
            }
            rep.push(StartValue {
                scope: format!("region {name}"),
                quantity: "temperature".into(),
                value: describe(t),
                source: "given".into(),
            });
        }
    }
    Ok(rep)
}

/// The inflow-weighted mean temperature of the inlets of the fluid part(s)
/// of `cells` (REG009 when they have none with a temperature).
fn inlet_temperature(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    env: &Env,
    cells: &[usize],
) -> Result<(f64, Vec<String>)> {
    let comp = crate::pressure_values::components(case, mesh);
    let mut mine = vec![false; comp.n];
    for &c in cells {
        if comp.of[c] != u32::MAX {
            mine[comp.of[c] as usize] = true;
        }
    }
    let values = case.compiled_boundaries()?;
    let (mut sw, mut swt) = (0.0, 0.0);
    let (mut t_lo, mut t_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut names = Vec::new();
    for b in &mesh.boundaries {
        let Some(bcond) = case.boundary.get(&b.name) else {
            continue;
        };
        if !crate::turbulence::inlet::is_inlet(bcond) {
            continue;
        }
        let Some(t) = case
            .boundary_energy
            .get(&b.name)
            .and_then(|e| e.temperature.as_ref())
        else {
            continue;
        };
        let tv = crate::expr::value(
            t,
            env,
            Vars::SPACE,
            &format!("boundary.{}.energy.temperature", b.name),
        )?;
        let vals = values.get(&b.name);
        let mut used = false;
        for f in b.range() {
            let k = comp.of[mesh.owner[f] as usize];
            if k == u32::MAX || !mine[k as usize] {
                continue;
            }
            let x = geo.face_centroid[f];
            let s = geo.face_area[f];
            let a = norm(s);
            // the inflow: |u·n| A of a velocity inlet, the area otherwise
            let wgt = match (bcond, vals.and_then(|v| v.velocity.as_ref())) {
                (BoundaryCondition::VelocityInlet { .. }, Some(v)) => {
                    let uu = [v[0].at(x, 0.0), v[1].at(x, 0.0), v[2].at(x, 0.0)];
                    let un = dot(uu, s).abs();
                    if un > 0.0 { un } else { a }
                }
                _ => a,
            };
            let t = tv.at(x, 0.0);
            sw += wgt;
            swt += wgt * t;
            t_lo = t_lo.min(t);
            t_hi = t_hi.max(t);
            used = true;
        }
        if used {
            names.push(b.name.clone());
        }
    }
    if sw <= 0.0 {
        return Err(err(
            diag::REG_INITIAL_INLET,
            "regions".into(),
            "temperature = \"inlet\" in [regions.<name>.initial]: the fluid part of the region \
             has no inlet with `[boundary.<name>.energy] temperature`"
                .into(),
        ));
    }
    // one inlet temperature: exactly that value (no round-off of the mean)
    let t = if t_lo.to_bits() == t_hi.to_bits() {
        t_lo
    } else {
        (swt / sw).clamp(t_lo, t_hi)
    };
    Ok((t, names))
}

/// The region turbulence start values over `k` and `w` (ω or ε, as the
/// model solves; `eps` true for the k-ε family, `c_mu` to convert).
#[allow(clippy::too_many_arguments)]
pub fn turbulence(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    env: &Env,
    eps: bool,
    c_mu: f64,
    k: &mut [f64],
    w: &mut [f64],
) -> Result<()> {
    let nc = k.len();
    for (name, r) in &case.regions {
        let Some(i) = &r.initial else { continue };
        let (Some(ks), second) = (&i.k, (&i.omega, &i.epsilon)) else {
            continue;
        };
        let (ws, given_eps, key) = match second {
            (Some(o), _) => (o, false, "omega"),
            (None, Some(e)) => (e, true, "epsilon"),
            _ => continue,
        };
        let kv = crate::expr::value(ks, env, Vars::SPACE, &format!("regions.{name}.initial.k"))?;
        let wv = crate::expr::value(
            ws,
            env,
            Vars::SPACE,
            &format!("regions.{name}.initial.{key}"),
        )?;
        for c in cells_of(case, mesh, name).into_iter().filter(|&c| c < nc) {
            let x = geo.cell_centroid[c];
            let kk = kv.at(x, 0.0);
            let ww = wv.at(x, 0.0);
            k[c] = kk;
            w[c] = if given_eps == eps {
                ww
            } else if eps {
                c_mu * kk * ww
            } else {
                ww / (c_mu * kk).max(1e-300)
            };
        }
    }
    Ok(())
}

/// The turbulence start lines of the report.
pub fn turbulence_report(case: &FlowCase) -> Vec<StartValue> {
    case.regions
        .iter()
        .filter_map(|(name, r)| {
            let i = r.initial.as_ref()?;
            let k = i.k.as_ref()?;
            let (w, what) = match (&i.omega, &i.epsilon) {
                (Some(o), _) => (o, "omega"),
                (None, Some(e)) => (e, "epsilon"),
                _ => return None,
            };
            Some(StartValue {
                scope: format!("region {name}"),
                quantity: "turbulence".into(),
                value: format!("k = {}, {what} = {}", describe(k), describe(w)),
                source: "given".into(),
            })
        })
        .collect()
}

/// The REG010 lines of the log.
pub fn notes(rep: &[StartValue]) -> Vec<String> {
    rep.iter()
        .map(|s| {
            format!(
                "[REG010] start values, {}: {} {} ({})",
                s.scope, s.quantity, s.value, s.source
            )
        })
        .collect()
}

/// A start velocity and pressure field with its report lines.
pub type FlowStart = (foehncraft_fv::disc::Soa3, Vec<f64>, Vec<StartValue>);

/// The start velocity and pressure of a run that does not restart: the
/// case-wide values `u0`, `p0` with the region tables and the start
/// pressure of every part (`flow_fields`). `None` when they stay uniform
/// (the flow keeps the start it was built with).
pub fn initial_flow(
    case: &FlowCase,
    mesh: &Mesh,
    geo: &Geometry,
    u0: [f64; 3],
    p0: f64,
    reference: bool,
) -> Result<Option<FlowStart>> {
    let nc = mesh.n_cells;
    let mut u = foehncraft_fv::disc::Soa3 {
        x: vec![u0[0]; nc],
        y: vec![u0[1]; nc],
        z: vec![u0[2]; nc],
    };
    let mut p = vec![p0; nc];
    let env = case.env()?;
    let values = case.compiled_boundaries()?;
    Ok(
        flow_fields(case, mesh, geo, &env, &values, &mut u, &mut p, reference)?
            .map(|rep| (u, p, rep)),
    )
}

/// The temperature lines of the start report: the given value, or the
/// inlet temperature of the region's fluid part (on the flow mesh `mesh`).
pub fn temperature_report(case: &FlowCase, mesh: &Mesh, geo: &Geometry) -> Result<Vec<StartValue>> {
    let mut rep = Vec::new();
    let env = case.env()?;
    // the solid regions are in the full case of a conjugate run
    let full: &FlowCase = match &case.fluid_view {
        Some(v) => &v.full_case,
        None => case,
    };
    for (name, r) in &full.regions {
        let Some(i) = &r.initial else { continue };
        let Some(t) = &i.temperature else { continue };
        if i.temperature_from_inlet() {
            let cells = cells_of(case, mesh, name);
            if cells.is_empty() {
                continue;
            }
            let (t_in, inlets) = inlet_temperature(case, mesh, geo, &env, &cells)?;
            rep.push(StartValue {
                scope: format!("region {name}"),
                quantity: "temperature".into(),
                value: format!("{t_in}"),
                source: format!("inlet ({})", inlets.join(", ")),
            });
        } else {
            rep.push(StartValue {
                scope: format!("region {name}"),
                quantity: "temperature".into(),
                value: describe(t),
                source: "given".into(),
            });
        }
    }
    Ok(rep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn channel(inlet: &str, extra: &str) -> FlowCase {
        let text = format!(
            r#"case_format = "0.6"
[mesh.generate]
type = "hex"
dimension = 2
upper = [4.0, 1.0, 0.1]
cells = [16, 4, 1]

[fluid]
density = 1.0
viscosity = 0.05

[boundary.xmin]
{inlet}
[boundary.xmax]
type = "static-pressure-outlet"
pressure = 0.0
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
{extra}
"#
        );
        FlowCase::from_toml(&text, Path::new(".")).unwrap()
    }

    type Start = Option<(Vec<f64>, Vec<StartValue>)>;

    fn start(c: &FlowCase) -> (Mesh, Start) {
        let (mesh, geo, _) = crate::case::mesh_for(c, None).unwrap();
        let r = initial_flow(
            c,
            &mesh,
            &geo,
            c.initial.velocity,
            c.initial.pressure,
            false,
        )
        .unwrap()
        .map(|(_, p, rep)| (p, rep));
        (mesh, r)
    }

    /// V7.46: a channel between a total-pressure inlet at 10 Pa and an
    /// outlet at 0 Pa starts on the straight line between them (the
    /// Laplace solve on a uniform mesh is exact for a linear field).
    #[test]
    fn a_part_between_two_pressures_starts_on_the_laplace_solution() {
        let c = channel("type = \"total-pressure-inlet\"\ntotal_pressure = 10.0", "");
        let (mesh, r) = start(&c);
        let (p, rep) = r.expect("a Laplace start");
        let geo = foehncraft_mesh::Geometry::compute(&mesh);
        for (cell, &pc) in p.iter().enumerate() {
            let x = geo.cell_centroid[cell][0];
            assert!((pc - 10.0 * (1.0 - x / 4.0)).abs() < 1e-9, "{x}: {pc}");
        }
        assert_eq!(rep.len(), 1);
        assert!(
            rep[0]
                .source
                .starts_with("Laplace solve between 2 boundaries"),
            "{rep:?}"
        );
    }

    /// One pressure boundary at 0 Pa: the start stays the uniform start of
    /// before (no change, no report).
    #[test]
    fn one_outlet_at_zero_keeps_the_start() {
        let c = channel("type = \"velocity-inlet\"\nvelocity = [1.0, 0.0, 0.0]", "");
        assert!(start(&c).1.is_none());
        // a case-wide start pressure stays as given
        let c = channel(
            "type = \"total-pressure-inlet\"\ntotal_pressure = 10.0",
            "[initial]\npressure = 3.0",
        );
        assert!(start(&c).1.is_none());
    }
}
