//! Mesh geometry: face area vectors and centroids, cell volumes and
//! centroids.
//!
//! Faces may be non-planar polygons. A face is split into triangles fanning
//! out from the average of its nodes; its area vector is the sum of the
//! triangle area vectors, which depends only on the boundary loop of the
//! face, so the area vectors of a closed cell sum to zero up to round-off
//! for any face shape. The face centroid is the average of the triangle
//! centroids weighted by the triangle areas projected on the face normal.
//!
//! A cell is split into tetrahedra spanned by an interior point (the average
//! of its face centroids), the fan point of a face (its node average) and one
//! edge of the face, so cell volumes use the same face surfaces as the area
//! vectors and neighbouring cells tile the domain exactly. Volume and centroid are the sums over these tetrahedra (Gauss
//! divergence theorem). Geometry is always computed in float64, whatever the
//! floating-point type of the solver (D7).
//!
//! The functions are generic over [`Coord`], so the same source also runs
//! on dual numbers (`foehncraft_core::Dual`) and gives the exact derivatives
//! of areas, centroids and volumes with respect to the node coordinates
//! (D59, item 5: shape sensitivities, mesh morphing L4).
//!
//! The same decomposition is used by point location ([`crate::locate`]).

use crate::mesh::Mesh;
use foehncraft_core::{Dual, Real};
use rayon::prelude::*;

/// A coordinate type: [`Real`] with a square root (`f64`, or dual numbers
/// for derivatives with respect to the nodes).
pub trait Coord: Real {
    /// Square root.
    fn sqrt(self) -> Self;
}

impl Coord for f64 {
    #[inline(always)]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
}

impl<const N: usize> Coord for Dual<N> {
    /// √x with the derivatives d/(2√x).
    #[inline(always)]
    fn sqrt(self) -> Self {
        let r = self.v.sqrt();
        let s = 0.5 / r;
        let mut d = [0.0; N];
        for (o, i) in d.iter_mut().zip(self.d) {
            *o = s * i;
        }
        Dual::new(r, d)
    }
}

/// Vector sum.
#[inline(always)]
pub fn add<T: Coord>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Vector difference.
#[inline(always)]
pub fn sub<T: Coord>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scaled vector.
#[inline(always)]
pub fn scale<T: Coord>(a: [T; 3], s: T) -> [T; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
#[inline(always)]
pub fn dot<T: Coord>(a: [T; 3], b: [T; 3]) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline(always)]
pub fn cross<T: Coord>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm.
#[inline(always)]
pub fn norm<T: Coord>(a: [T; 3]) -> T {
    dot(a, a).sqrt()
}

/// Average of the nodes of a polygon.
#[inline]
pub fn node_average<T: Coord>(points: &[[T; 3]], nodes: &[u32]) -> [T; 3] {
    let mut c = [T::ZERO; 3];
    for &n in nodes {
        c = add(c, points[n as usize]);
    }
    scale(c, T::from_f64(1.0 / nodes.len() as f64))
}

/// Area vector and centroid of a (possibly non-planar) polygon.
pub fn polygon_area_centroid<T: Coord>(points: &[[T; 3]], nodes: &[u32]) -> ([T; 3], [T; 3]) {
    let k = nodes.len();
    if k == 3 {
        let (a, b, c) = (
            points[nodes[0] as usize],
            points[nodes[1] as usize],
            points[nodes[2] as usize],
        );
        let s = scale(cross(sub(b, a), sub(c, a)), T::from_f64(0.5));
        return (s, scale(add(add(a, b), c), T::from_f64(1.0 / 3.0)));
    }
    let m = node_average(points, nodes);
    let mut s = [T::ZERO; 3];
    let mut tri_s = [[T::ZERO; 3]; 32];
    let mut big: Vec<[T; 3]> = Vec::new();
    for i in 0..k {
        let p = points[nodes[i] as usize];
        let q = points[nodes[(i + 1) % k] as usize];
        let t = scale(cross(sub(p, m), sub(q, m)), T::from_f64(0.5));
        s = add(s, t);
        if k <= 32 {
            tri_s[i] = t;
        } else {
            big.push(t);
        }
    }
    let ns = norm(s);
    if ns == T::ZERO {
        return (s, m);
    }
    let n = scale(s, T::ONE / ns);
    let mut c = [T::ZERO; 3];
    let mut w = T::ZERO;
    let mut turned = false;
    for i in 0..k {
        let p = points[nodes[i] as usize];
        let q = points[nodes[(i + 1) % k] as usize];
        let t = if k <= 32 { tri_s[i] } else { big[i] };
        let a = dot(t, n);
        turned |= a.to_f64() <= 0.0;
        c = add(c, scale(add(add(p, q), m), a / T::from_f64(3.0)));
        w += a;
    }
    let centroid = scale(c, T::ONE / w);
    if turned && !in_node_box(points, nodes, centroid) {
        // a fan triangle turns against the mean normal and the centroid of
        // the signed weights has left the box of the nodes, where the true
        // centroid of a polygon always lies: a ribbon-thin face twisted
        // along its length, whose signed weights nearly cancel. Weigh by
        // the triangle areas instead, which keeps the centroid in the
        // convex hull of the nodes. Every other face, planar or warped,
        // convex or not, keeps the signed weights.
        let mut c = [T::ZERO; 3];
        let mut w = T::ZERO;
        for i in 0..k {
            let p = points[nodes[i] as usize];
            let q = points[nodes[(i + 1) % k] as usize];
            let t = if k <= 32 { tri_s[i] } else { big[i] };
            let a = norm(t);
            c = add(c, scale(add(add(p, q), m), a / T::from_f64(3.0)));
            w += a;
        }
        if w.to_f64() > 0.0 {
            return (s, scale(c, T::ONE / w));
        }
        return (s, m);
    }
    (s, centroid)
}

/// Relative margin of the node box in [`in_node_box`] (of its largest side).
const BOX_MARGIN: f64 = 1e-9;

/// `true` if `x` lies in the bounding box of the polygon's nodes, widened
/// by [`BOX_MARGIN`].
fn in_node_box<T: Coord>(points: &[[T; 3]], nodes: &[u32], x: [T; 3]) -> bool {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &v in nodes {
        let p = points[v as usize];
        for d in 0..3 {
            lo[d] = lo[d].min(p[d].to_f64());
            hi[d] = hi[d].max(p[d].to_f64());
        }
    }
    let side = (0..3).map(|d| hi[d] - lo[d]).fold(0.0, f64::max);
    let pad = BOX_MARGIN * side;
    (0..3).all(|d| {
        let v = x[d].to_f64();
        v >= lo[d] - pad && v <= hi[d] + pad
    })
}

/// Below this many faces [`Geometry::faces`] runs on one thread.
const SERIAL_FACES: usize = 50_000;

/// Geometry of a mesh, for the local cells and faces.
#[derive(Debug, Clone, Default)]
pub struct Geometry {
    /// Area vector of each face (normal out of the owner, length = area).
    pub face_area: Vec<[f64; 3]>,
    /// Centroid of each face.
    pub face_centroid: Vec<[f64; 3]>,
    /// Volume of each cell (local cells first; halo cells, if any, after
    /// them, filled by the halo exchange or by the decomposition).
    pub cell_volume: Vec<f64>,
    /// Centroid of each cell (local and halo).
    pub cell_centroid: Vec<[f64; 3]>,
}

impl Geometry {
    /// Computes the geometry of all faces and local cells. Halo entries are
    /// left at zero.
    pub fn compute(mesh: &Mesh) -> Self {
        let (face_area, face_centroid) = Self::faces(mesh);
        let mut cell_volume = vec![0.0; mesh.n_cells + mesh.n_halo];
        let mut cell_centroid = vec![[0.0; 3]; mesh.n_cells + mesh.n_halo];
        cell_volume[..mesh.n_cells]
            .par_iter_mut()
            .zip(cell_centroid[..mesh.n_cells].par_iter_mut())
            .enumerate()
            .for_each(|(c, (v, x))| {
                let (vol, cen) = cell_volume_centroid(mesh, &face_centroid, c);
                *v = vol;
                *x = cen;
            });
        Self {
            face_area,
            face_centroid,
            cell_volume,
            cell_centroid,
        }
    }

    /// The area vectors and centroids of all faces (the face part of
    /// [`Geometry::compute`]).
    /// Small meshes run on the calling thread: handing the work to the
    /// thread pool would cost more than it saves (same result).
    pub fn faces(mesh: &Mesh) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let one = |f: usize| face_area_centroid(mesh, &mesh.points, f);
        if mesh.n_faces() < SERIAL_FACES {
            (0..mesh.n_faces()).map(one).unzip()
        } else {
            (0..mesh.n_faces()).into_par_iter().map(one).unzip()
        }
    }

    /// Sum of the outward area vectors of local cell `c` (zero for a closed
    /// cell).
    pub fn closure(&self, mesh: &Mesh, c: usize) -> [f64; 3] {
        let mut s = [0.0; 3];
        for (f, is_nb) in mesh.cell_face_iter(c) {
            let a = self.face_area[f];
            s = if is_nb { sub(s, a) } else { add(s, a) };
        }
        s
    }

    /// Unit normal and area of face `f`.
    #[inline]
    pub fn normal_area(&self, f: usize) -> ([f64; 3], f64) {
        let a = norm(self.face_area[f]);
        (scale(self.face_area[f], 1.0 / a), a)
    }
}

/// Area vector (out of the local owner) and centroid of face `f`, evaluated
/// in the node order of the global mesh: a flipped face of a local mesh
/// gets the global bits with the sign of the area vector reversed, which is
/// exact (D93 (3)).
#[inline]
pub fn face_area_centroid<T: Coord>(mesh: &Mesh, points: &[[T; 3]], f: usize) -> ([T; 3], [T; 3]) {
    let (nodes, flipped) = mesh.face_global_order(f);
    let (s, c) = polygon_area_centroid(points, &nodes);
    if flipped {
        ([-s[0], -s[1], -s[2]], c)
    } else {
        (s, c)
    }
}

/// Interior point of a local cell used for its tetrahedral decomposition:
/// the average of its face centroids.
pub fn cell_apex<T: Coord>(mesh: &Mesh, face_centroid: &[[T; 3]], c: usize) -> [T; 3] {
    let mut x = [T::ZERO; 3];
    let mut k = T::ZERO;
    for (f, _) in mesh.cell_face_iter(c) {
        x = add(x, face_centroid[f]);
        k += T::ONE;
    }
    scale(x, T::ONE / k)
}

/// Volume and centroid of local cell `c`.
pub fn cell_volume_centroid(mesh: &Mesh, face_centroid: &[[f64; 3]], c: usize) -> (f64, [f64; 3]) {
    cell_volume_centroid_of(mesh, &mesh.points, face_centroid, c)
}

/// Volume and centroid of local cell `c` of `mesh` with the node
/// coordinates `points` (of any [`Coord`] type; `face_centroid` from the
/// same points).
pub fn cell_volume_centroid_of<T: Coord>(
    mesh: &Mesh,
    points: &[[T; 3]],
    face_centroid: &[[T; 3]],
    c: usize,
) -> (T, [T; 3]) {
    let apex = cell_apex(mesh, face_centroid, c);
    let mut vol = T::ZERO;
    let mut cen = [T::ZERO; 3];
    for (f, is_nb) in mesh.cell_face_iter(c) {
        // the face in global node order and with the global sense, so a
        // local mesh sums the same tetrahedra in the same order as the
        // global mesh (cell faces are listed in global face order)
        let (nodes, flipped) = mesh.face_global_order(f);
        let is_nb = is_nb != flipped;
        let nodes = &nodes[..];
        let k = nodes.len();
        // fan about the node average, the same surface as the area vector
        let fc = if k == 3 {
            face_centroid[f]
        } else {
            node_average(points, nodes)
        };
        for i in 0..k {
            let p = points[nodes[i] as usize];
            let q = points[nodes[(i + 1) % k] as usize];
            // tetrahedron (apex, fc, p, q), outward orientation of (fc, p, q)
            let mut v = dot(sub(fc, apex), cross(sub(p, apex), sub(q, apex))) / T::from_f64(6.0);
            if is_nb {
                v = -v;
            }
            vol += v;
            let tc = scale(add(add(apex, fc), add(p, q)), T::from_f64(0.25));
            cen = add(cen, scale(tc, v));
        }
    }
    (vol, scale(cen, T::ONE / vol))
}

/// The volume that face `f` sweeps when the nodes move on straight lines
/// from `old` to `new` (m³; positive when the face moves out of its local
/// owner), for the space conservation law of a deforming mesh (D91, stage
/// MS0; Demirdžić & Perić 1988).
///
/// The face is the fan of triangles about its node average that the area
/// vector and the cell volumes use. A triangle whose corners move linearly
/// in τ ∈ [0, 1] sweeps ∫ **d̄**·**S**(τ) dτ, with **d̄** the mean of the three
/// corner displacements (the velocity is linear over the triangle) and
/// **S**(τ) its area vector, which is quadratic in τ, so Simpson's rule is
/// exact. Since the fan point is a node average it moves on a straight line
/// too, and the swept volumes of the faces of a cell add up to the change
/// of its volume to round-off: Σ_f ±δV_f = V(new) − V(old), for any motion
/// of the nodes and any face shape.
pub fn swept_volume(mesh: &Mesh, old: &[[f64; 3]], new: &[[f64; 3]], f: usize) -> f64 {
    let (nodes, flipped) = mesh.face_global_order(f);
    let nodes = &nodes[..];
    let k = nodes.len();
    let half = |n: u32| scale(add(old[n as usize], new[n as usize]), 0.5);
    let (m0, m1) = (node_average(old, nodes), node_average(new, nodes));
    let mh = scale(add(m0, m1), 0.5);
    let dm = sub(m1, m0);
    let mut v = 0.0;
    for i in 0..k {
        let (a, b) = (nodes[i], nodes[(i + 1) % k]);
        let s0 = cross(sub(old[a as usize], m0), sub(old[b as usize], m0));
        let s1 = cross(sub(new[a as usize], m1), sub(new[b as usize], m1));
        let sh = cross(sub(half(a), mh), sub(half(b), mh));
        // ∫ S dτ = ½ (S₀ + 4 S½ + S₁)/6 with S = ½ (p − m) × (q − m)
        let s = scale(add(add(s0, scale(sh, 4.0)), s1), 1.0 / 12.0);
        let d = scale(
            add(
                add(dm, sub(new[a as usize], old[a as usize])),
                sub(new[b as usize], old[b as usize]),
            ),
            1.0 / 3.0,
        );
        v += dot(d, s);
    }
    if flipped { -v } else { v }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn swept_volumes_close_every_cell() {
        // a perturbed hex mesh and a prism mesh, moved by a smooth motion
        // with shear, stretch and turn: Σ ±δV_f = ΔV per cell to round-off
        for kind in [
            crate::generate::Kind::Hex,
            crate::generate::Kind::Tet,
            crate::generate::Kind::Poly,
            crate::generate::Kind::Mixed,
        ] {
            let spec = crate::generate::BoxSpec {
                kind,
                cells: [4, 3, 3],
                perturb: 0.2,
                ..Default::default()
            };
            let mesh = spec.generate().expect("mesh");
            let old = mesh.points.clone();
            let new: Vec<[f64; 3]> = old
                .iter()
                .map(|p| {
                    [
                        p[0] + 0.1 * (3.0 * p[1]).sin() * p[2],
                        p[1] * 0.8 + 0.05 * p[0] * p[0],
                        p[2] + 0.1 * p[0] * p[1],
                    ]
                })
                .collect();
            let g0 = Geometry::compute(&mesh);
            let mut moved = mesh.clone();
            moved.points = new.clone();
            let g1 = Geometry::compute(&moved);
            let swept: Vec<f64> = (0..mesh.n_faces())
                .map(|f| swept_volume(&mesh, &old, &new, f))
                .collect();
            for c in 0..mesh.n_cells {
                let mut s = 0.0;
                for (f, is_nb) in mesh.cell_face_iter(c) {
                    s += if is_nb { -swept[f] } else { swept[f] };
                }
                let dv = g1.cell_volume[c] - g0.cell_volume[c];
                assert!(
                    (s - dv).abs() <= 1e-14 * g0.cell_volume[c].max(1e-30) * 100.0,
                    "{kind:?} cell {c}: Σ δV = {s:e}, ΔV = {dv:e}"
                );
            }
            // a pure translation sweeps d·S on every face
            let d = [0.01, -0.02, 0.03];
            let shifted: Vec<[f64; 3]> = old.iter().map(|p| add(*p, d)).collect();
            for f in 0..mesh.n_faces() {
                let v = swept_volume(&mesh, &old, &shifted, f);
                let e = dot(d, g0.face_area[f]);
                assert!(
                    (v - e).abs() <= 1e-15 + 1e-12 * e.abs(),
                    "face {f}: {v:e} vs {e:e}"
                );
            }
        }
    }

    #[test]
    fn unit_square_and_warped_quad() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.2],
        ];
        let (s, c) = polygon_area_centroid(&pts, &[0, 1, 2, 3]);
        assert!((s[2] - 1.0).abs() < 1e-15 && s[0] == 0.0 && s[1] == 0.0);
        assert!((c[0] - 0.5).abs() < 1e-15 && (c[1] - 0.5).abs() < 1e-15);
        // a warped quad has the same projected area as its loop
        let (s, _) = polygon_area_centroid(&pts, &[0, 1, 4, 3]);
        assert!((s[2] - 1.0).abs() < 1e-15);
        // triangle centroid
        let (s, c) = polygon_area_centroid(&pts, &[0, 1, 3]);
        assert!((s[2] - 0.5).abs() < 1e-15);
        assert!((c[0] - 1.0 / 3.0).abs() < 1e-15);
    }

    /// D59, item 5: areas, centroids and volumes on dual numbers give the
    /// exact derivatives with respect to the node coordinates (against
    /// central differences), the values are the f64 values bit for bit, and
    /// moving an interior node leaves the total volume unchanged.
    #[test]
    fn shape_derivatives_of_the_geometry() {
        type D = Dual<3>;
        let sec: crate::case::MeshSection = toml::from_str(
            "[generate]\ntype = \"mixed\"\ndimension = 3\ncells = [3, 3, 2]\nperturb = 0.2\nseed = 5\n",
        )
        .unwrap();
        let (mesh, _, _) = crate::case::build_mesh(&sec).unwrap();
        let nf = mesh.n_faces();
        let geometry = |pts: &[[f64; 3]]| -> Vec<f64> {
            let fa: Vec<_> = (0..nf)
                .map(|f| polygon_area_centroid(pts, mesh.face(f)))
                .collect();
            let fc: Vec<_> = fa.iter().map(|a| a.1).collect();
            let mut out = Vec::new();
            for (s, c) in &fa {
                out.extend(s.iter().chain(c));
            }
            for c in 0..mesh.n_cells {
                let (v, x) = cell_volume_centroid_of(&mesh, pts, &fc, c);
                out.push(v);
                out.extend(x);
            }
            out
        };
        let base = geometry(&mesh.points);
        // the generic code on f64 is the code of Geometry::compute
        let g = Geometry::compute(&mesh);
        assert_eq!(base[..6], [g.face_area[0], g.face_centroid[0]].concat()[..]);
        let boundary: std::collections::BTreeSet<u32> = (mesh.n_interior_faces()..nf)
            .flat_map(|f| mesh.face(f).iter().copied())
            .collect();
        let h = 1e-6;
        let mut checked = 0;
        let mut worst = 0.0f64;
        for j in 0..mesh.points.len() {
            let pts: Vec<[D; 3]> = mesh
                .points
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    std::array::from_fn(|k| {
                        if i == j {
                            D::variable(p[k], k)
                        } else {
                            D::constant(p[k])
                        }
                    })
                })
                .collect();
            let fa: Vec<_> = (0..nf)
                .map(|f| polygon_area_centroid(&pts, mesh.face(f)))
                .collect();
            let fc: Vec<_> = fa.iter().map(|a| a.1).collect();
            let mut dual = Vec::new();
            for (s, c) in &fa {
                dual.extend(s.iter().chain(c).copied());
            }
            let mut dvol = [0.0; 3];
            for c in 0..mesh.n_cells {
                let (v, x) = cell_volume_centroid_of(&mesh, &pts, &fc, c);
                dual.push(v);
                dual.extend(x);
                for (s, d) in dvol.iter_mut().zip(v.d) {
                    *s += d;
                }
            }
            for (a, b) in dual.iter().zip(&base) {
                assert!(a.v == *b);
            }
            for k in 0..3 {
                let mut p = mesh.points.clone();
                p[j][k] += h;
                let gp = geometry(&p);
                p[j][k] -= 2.0 * h;
                let gm = geometry(&p);
                for (i, a) in dual.iter().enumerate() {
                    let fd = (gp[i] - gm[i]) / (2.0 * h);
                    // round-off of the differences: about 1e-10 (coordinates of size 1)
                    let scale = a.d[k].abs().max(fd.abs()).max(1e-3);
                    worst = worst.max((a.d[k] - fd).abs() / scale);
                    assert!(
                        (a.d[k] - fd).abs() < 1e-6 * scale,
                        "node {j}, direction {k}, output {i}: dual {:e}, central difference {fd:e}",
                        a.d[k]
                    );
                    checked += 1;
                }
                if !boundary.contains(&(j as u32)) {
                    assert!(
                        dvol[k].abs() < 1e-14,
                        "interior node {j}: dV/dx = {:e}",
                        dvol[k]
                    );
                }
            }
        }
        println!(
            "{} cells, {} nodes: {checked} derivatives, largest relative difference {worst:.1e}",
            mesh.n_cells,
            mesh.points.len()
        );
        assert!(checked > 10_000);
    }

    /// The twisted ribbon face found by the point-search benchmark (#235):
    /// thickness 1.3e-5, length 0.017, twisted by 1.7e-3 across. The signed
    /// weights put its centroid at z = 0.19; the area weights keep it on
    /// the face.
    #[test]
    fn twisted_ribbon_face_keeps_its_centroid_on_the_face() {
        let pts = vec![
            [0.43840603849876486, 0.0, 0.0],
            [0.4367040568405628, 1.2137783485223362e-5, 0.0],
            [
                0.4379280535783903,
                1.2897727829962324e-5,
                0.016268504775330917,
            ],
            [0.43617218170643896, 0.0, 0.016969859729795846],
        ];
        let (_, c) = polygon_area_centroid(&pts, &[0, 1, 2, 3]);
        for d in 0..3 {
            let lo = pts.iter().map(|p| p[d]).fold(f64::INFINITY, f64::min);
            let hi = pts.iter().map(|p| p[d]).fold(f64::NEG_INFINITY, f64::max);
            assert!(
                c[d] >= lo && c[d] <= hi,
                "axis {d}: {} not in [{lo}, {hi}]",
                c[d]
            );
        }
    }

    /// A planar non-convex polygon whose fan triangles from the node
    /// average turn both ways keeps the exact centroid of the signed
    /// weights.
    #[test]
    fn planar_non_convex_polygon_keeps_the_exact_centroid() {
        // a thin L-shaped hexagon in the plane z = 1: unit squares at
        // (0..3, 0..1) and (0..1, 1..3), area 5
        let pts = vec![
            [0.0, 0.0, 1.0],
            [3.0, 0.0, 1.0],
            [3.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 3.0, 1.0],
            [0.0, 3.0, 1.0],
        ];
        let (s, c) = polygon_area_centroid(&pts, &[0, 1, 2, 3, 4, 5]);
        assert!((s[2] - 5.0).abs() < 1e-14);
        // centroid: (3 * (1.5, 0.5) + 2 * (0.5, 2)) / 5 = (1.1, 1.1)
        assert!(
            (c[0] - 1.1).abs() < 1e-14 && (c[1] - 1.1).abs() < 1e-14,
            "{c:?}"
        );
        assert_eq!(c[2], 1.0);
    }
}
