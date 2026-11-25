//! Reverse sweep of a face kernel from local dual Jacobians (L17 stage A0,
//! decision D221 option A): the transposed Jacobian product of
//! [`crate::flow::momentum_face`] without a tape, without hand-derived
//! code and without atomics.
//!
//! The discrete adjoint (D220) needs v̄ ↦ (∂R/∂q)ᵀ v̄. A face kernel reads
//! a few values of its two cells and of the face, so its Jacobian is a
//! small dense block. This module
//!
//! 1. copies the inputs of one face into a two-cell local problem and runs
//!    the *same* kernel body on [`Dual<LOCAL>`], one seeded direction per
//!    local input ([`momentum_face_jacobian`]): the block comes from the
//!    kernel source by instantiation (D59), exact up to round-off, with the
//!    switches (upwind direction, `max`/`min`) frozen at the state;
//! 2. multiplies the transposed block with the cotangent of the face
//!    outputs ([`pullback`]), one face at a time, written by that face only;
//! 3. sums the cell parts of the face cotangents over the faces of each
//!    cell in their stored order ([`gather_cell`]): the reverse of a
//!    gather is a scatter, and the scatter is written as a gather over the
//!    inverse stencil (cell → its faces, `Topology::cell_faces`), so the
//!    result is bitwise the same for every thread count and on the GPU
//!    (D31, D93).
//!
//! The block can be stored per face (about 1 kB, fast) or recomputed per
//! product (no memory): both give bitwise the same numbers, since step 1
//! is deterministic.
//!
//! **Local inputs** ([`LOCAL`] = 26 directions, index constants below):
//! the face mass flux ṁ_f, the face viscosity μ_f, the owner velocity
//! u_P, the neighbour velocity u_N (on a boundary face: the boundary
//! velocity u_b), the owner gradient ∇u_P and the neighbour gradient ∇u_N
//! (row i is ∇u_i). The velocity of the central correction (`uc`) is the
//! velocity of a steady state, so it takes the seed of `u`. Frozen
//! (constants of the linearisation): the limiter ψ, the geometry (shape
//! derivatives are stage A2) and the scheme switches. Periodic faces are
//! not supported yet (their neighbour does not list them); the functions
//! return `None` for them.
//!
//! **Outputs** ([`OUT`] = 5): c_P, c_N, b_x, b_y, b_z of
//! [`crate::flow::momentum_face`].
//!
//! All functions here are safe: every index is checked while the local
//! copy is made, so a wrong index gives `None`, never an unchecked load.

use crate::flow::{self, FaceGeometry, MomentumInputs, VelocityGradient};
use crate::fv::{Topology, Vec3};
use crate::unpack_face;
use foehncraft_core::{Dual, Real};

/// Number of local inputs (seeded directions) of one face.
pub const LOCAL: usize = 26;
/// Number of outputs of one face.
pub const OUT: usize = 5;
/// Local input: face mass flux ṁ_f.
pub const IN_MASS_FLUX: usize = 0;
/// Local input: face viscosity μ_f.
pub const IN_MU: usize = 1;
/// Local inputs 2..5: owner velocity u_P.
pub const IN_U_OWNER: usize = 2;
/// Local inputs 5..8: neighbour velocity u_N, or the boundary velocity u_b
/// on a boundary face.
pub const IN_U_OTHER: usize = 5;
/// Local inputs 8..17: owner gradient, entry 8 + 3 i + j = ∂u_i/∂x_j.
pub const IN_GRAD_OWNER: usize = 8;
/// Local inputs 17..26: neighbour gradient (interior faces).
pub const IN_GRAD_NEIGHBOUR: usize = 17;
/// Number of cell cotangents gathered per cell: ū (3) and the gradient
/// ∇̄u (9, the layout of [`IN_GRAD_OWNER`]).
pub const CELL: usize = 12;

/// The Jacobian block of one face: `jac[o][j]` = ∂(output o)/∂(local
/// input j).
pub type LocalJacobian = [[f64; LOCAL]; OUT];

type D = Dual<LOCAL>;

fn mk<'a>(a: [&'a [D]; 3]) -> Vec3<'a, D> {
    Vec3 {
        x: a[0],
        y: a[1],
        z: a[2],
    }
}

fn cells3(a: &[[D; 2]; 3], nc: usize) -> [&[D]; 3] {
    [&a[0][..nc], &a[1][..nc], &a[2][..nc]]
}

fn at(s: &[f64], i: usize) -> Option<f64> {
    s.get(i).copied()
}

fn at3(v: &Vec3<f64>, i: usize) -> Option<[f64; 3]> {
    Some([at(v.x, i)?, at(v.y, i)?, at(v.z, i)?])
}

fn grad_at(g: &VelocityGradient<f64>, c: usize) -> Option<[[f64; 3]; 3]> {
    Some([at3(&g.x, c)?, at3(&g.y, c)?, at3(&g.z, c)?])
}

/// The Jacobian block and the values (c_P, c_N, b) of
/// [`flow::momentum_face`] on face `f`, from the same kernel body on dual
/// numbers. `None` when an index of `f` is outside its array or `f` is a
/// periodic face.
pub fn momentum_face_jacobian(
    f: usize,
    inp: &MomentumInputs<f64>,
) -> Option<(LocalJacobian, [f64; OUT])> {
    let t = &inp.topo;
    let g = &inp.geo;
    let interior = f < t.n_interior;
    if interior && f >= t.periodic {
        return None;
    }
    let p = *t.owner.get(f)? as usize;
    let n = if interior {
        Some(*t.neighbour.get(f)? as usize)
    } else {
        None
    };
    let fb = f.wrapping_sub(t.n_interior);

    // seeded local inputs
    let var = |v: f64, k: usize| D::variable(v, k);
    let c = D::constant;
    let m = [var(at(inp.mass_flux, f)?, IN_MASS_FLUX)];
    let mu = [var(at(inp.mu, f)?, IN_MU)];
    let cells: &[usize] = match n {
        Some(ref n) => &[p, *n],
        None => &[p],
    };
    let mut u = [[D::ZERO; 2]; 3];
    let mut uc = [[D::ZERO; 2]; 3];
    let mut gr = [[[D::ZERO; 2]; 3]; 3];
    let mut psi = [[D::ZERO; 2]; 3];
    for (slot, &cell) in cells.iter().enumerate() {
        let (ui, gi) = if slot == 0 {
            (IN_U_OWNER, IN_GRAD_OWNER)
        } else {
            (IN_U_OTHER, IN_GRAD_NEIGHBOUR)
        };
        let uv = at3(&inp.u, cell)?;
        let ucv = at3(&inp.uc, cell)?;
        let gv = grad_at(&inp.grad, cell)?;
        let pv = at3(&inp.psi, cell)?;
        for i in 0..3 {
            u[i][slot] = var(uv[i], ui + i);
            uc[i][slot] = var(ucv[i], ui + i);
            psi[i][slot] = c(pv[i]);
            for j in 0..3 {
                gr[i][j][slot] = var(gv[i][j], gi + 3 * i + j);
            }
        }
    }
    let mut ub = [[D::ZERO; 1]; 3];
    let mut kind = [0u32; 1];
    if !interior {
        let v = at3(&inp.ub, fb)?;
        for i in 0..3 {
            ub[i][0] = var(v[i], IN_U_OTHER + i);
        }
        kind[0] = *inp.mom_kind.get(fb)?;
    }

    // geometry of the face and its cells, as constants
    let s = at3(&g.area, f)?;
    let xf = at3(&g.face_centroid, f)?;
    let kf = at3(&g.k, f)?;
    let area = [[c(s[0])], [c(s[1])], [c(s[2])]];
    let fc = [[c(xf[0])], [c(xf[1])], [c(xf[2])]];
    let kk = [[c(kf[0])], [c(kf[1])], [c(kf[2])]];
    let dgeo = [c(at(g.dgeo, f)?)];
    let weight = [c(if interior { at(g.weight, f)? } else { 0.0 })];
    let mut cc = [[D::ZERO; 2]; 3];
    for (slot, &cell) in cells.iter().enumerate() {
        let x = at3(&g.cell_centroid, cell)?;
        for i in 0..3 {
            cc[i][slot] = c(x[i]);
        }
    }

    // the local problem: face 0, owner 0, neighbour 1
    let nc = cells.len();
    let owner = [0u32];
    let neighbour = [1u32];
    let topo = Topology {
        owner: &owner,
        neighbour: if interior { &neighbour } else { &[] },
        cell_ptr: &[],
        cell_faces: &[],
        n_interior: usize::from(interior),
        periodic: usize::from(interior),
    };

    let empty: &[D] = &[];
    let local = MomentumInputs {
        topo,
        geo: FaceGeometry {
            area: mk([&area[0], &area[1], &area[2]]),
            face_centroid: mk([&fc[0], &fc[1], &fc[2]]),
            cell_centroid: mk(cells3(&cc, nc)),
            dgeo: &dgeo,
            k: mk([&kk[0], &kk[1], &kk[2]]),
            weight: &weight,
            shift: Vec3 {
                x: empty,
                y: empty,
                z: empty,
            },
        },
        mass_flux: &m,
        mu: &mu,
        u: mk(cells3(&u, nc)),
        uc: mk(cells3(&uc, nc)),
        grad: VelocityGradient {
            x: mk(cells3(&gr[0], nc)),
            y: mk(cells3(&gr[1], nc)),
            z: mk(cells3(&gr[2], nc)),
        },
        psi: mk(cells3(&psi, nc)),
        ub: if interior {
            Vec3 {
                x: empty,
                y: empty,
                z: empty,
            }
        } else {
            mk([&ub[0], &ub[1], &ub[2]])
        },
        mom_kind: if interior { &[] } else { &kind },
        blend: c(inp.blend),
        central: c(inp.central),
        beta: c(inp.beta),
        nonorth: c(inp.nonorth),
        convection: c(inp.convection),
        transposed: c(inp.transposed),
        trace_free: c(inp.trace_free),
    };
    // SAFETY: the local problem has one face (index 0) whose owner (0) and
    // neighbour (1, interior faces only) index the local cell arrays of
    // length `nc`; the face arrays have one entry, the boundary arrays one
    // entry on a boundary face (index 0 − n_interior = 0); the face is not
    // periodic (`periodic` = `n_interior`), so `shift` is not read.
    let (cp, cn, b) = unsafe { flow::momentum_face(0, &local) };
    let out = [cp, cn, b[0], b[1], b[2]];
    let mut jac = [[0.0; LOCAL]; OUT];
    let mut val = [0.0; OUT];
    for (o, y) in out.iter().enumerate() {
        jac[o] = y.d;
        val[o] = y.v;
    }
    Some((jac, val))
}

/// x̄ = Jᵀ ȳ for one face: the cotangent of the local inputs from the
/// cotangent `ybar` of the outputs, summed over the outputs in order.
pub fn pullback(jac: &LocalJacobian, ybar: &[f64; OUT]) -> [f64; LOCAL] {
    let mut x = [0.0; LOCAL];
    for (row, &y) in jac.iter().zip(ybar) {
        for (a, &j) in x.iter_mut().zip(row) {
            *a += j * y;
        }
    }
    x
}

/// ẏ = J ẋ for one face (the tangent product, for the dot-product test).
pub fn push_forward(jac: &LocalJacobian, xdot: &[f64; LOCAL]) -> [f64; OUT] {
    let mut y = [0.0; OUT];
    for (a, row) in y.iter_mut().zip(jac) {
        for (&j, &x) in row.iter().zip(xdot) {
            *a += j * x;
        }
    }
    y
}

/// The cell cotangents of cell `c`: ū_c (entries 0..3) and ∇̄u_c (3..12)
/// summed over the faces of the cell in stored order from the face
/// cotangents `face_bar` (one [`pullback`] per face). The owner side of a
/// face gives its [`IN_U_OWNER`] and [`IN_GRAD_OWNER`] parts, the
/// neighbour side its [`IN_U_OTHER`] and [`IN_GRAD_NEIGHBOUR`] parts; the
/// boundary velocity part of a boundary face is not a cell value and stays
/// with the face. `None` when an index is outside its array.
pub fn gather_cell(c: usize, topo: &Topology, face_bar: &[[f64; LOCAL]]) -> Option<[f64; CELL]> {
    let lo = *topo.cell_ptr.get(c)? as usize;
    let hi = *topo.cell_ptr.get(c + 1)? as usize;
    let mut out = [0.0; CELL];
    for &e in topo.cell_faces.get(lo..hi)? {
        let (f, is_nb) = unpack_face(e);
        let x = face_bar.get(f)?;
        let (u, g) = if is_nb {
            (IN_U_OTHER, IN_GRAD_NEIGHBOUR)
        } else {
            (IN_U_OWNER, IN_GRAD_OWNER)
        };
        for i in 0..3 {
            out[i] += x[u + i];
        }
        for i in 0..9 {
            out[3 + i] += x[g + i];
        }
    }
    Some(out)
}
