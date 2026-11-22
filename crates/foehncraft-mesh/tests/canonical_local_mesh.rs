//! MPI plan step W1 (D93 (1), (3)): the local mesh of every rank is
//! canonical, so its geometry and its per-cell sums over faces are bitwise
//! those of the one-rank run on the global mesh, for 1 to 4 ranks, METIS and
//! recursive coordinate bisection, hexahedra, tetrahedra, polyhedra and
//! mixed meshes, 2D and 3D, one and two halo layers, and also when the
//! global mesh was renumbered before the decomposition.
//!
//! The decomposition is held in one process ([`HaloPlan::exchange_all`]);
//! the example `mpi_matrix` checks the same across MPI processes in CI.

use foehncraft_mesh::Geometry;
use foehncraft_mesh::Mesh;
use foehncraft_mesh::decompose::{HaloPlan, LocalMesh, canonical_check, decompose};
use foehncraft_mesh::generate::{BoxSpec, Kind};
use foehncraft_mesh::geometry::{dot, sub};
use foehncraft_mesh::partition::{Method, partition};
use foehncraft_mesh::renumber::{Ordering, renumber};
use std::collections::HashMap;

fn methods() -> Vec<Method> {
    if cfg!(feature = "metis") {
        vec![Method::Rcb, Method::Metis]
    } else {
        vec![Method::Rcb]
    }
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

fn same3(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|d| same(a[d], b[d]))
}

/// One row of a pure diffusion matrix (|S|²/(d·S) per face): the diagonal
/// summed over the faces in the order of the cell's face list, and the
/// off-diagonals keyed by the global identity of the neighbour cell.
fn rows(mesh: &Mesh, geo: &Geometry) -> HashMap<u64, (f64, HashMap<u64, f64>)> {
    let ni = mesh.n_interior_faces();
    let mut out = HashMap::new();
    for c in 0..mesh.n_cells {
        let mut diag = 0.0;
        let mut off = HashMap::new();
        for (f, _) in mesh.cell_face_iter(c) {
            // the local stored orientation: for a flipped face both the
            // area vector and the centroid difference change sign, which
            // leaves the coefficient bitwise unchanged
            let p = mesh.owner[f] as usize;
            let s = geo.face_area[f];
            let d = if f < ni {
                let q = mesh.neighbour[f] as usize;
                sub(geo.cell_centroid[q], geo.cell_centroid[p])
            } else {
                sub(geo.face_centroid[f], geo.cell_centroid[p])
            };
            let a = dot(s, s) / dot(d, s);
            if f < ni {
                let q = mesh.neighbour[f] as usize;
                off.insert(mesh.cell_id[if p == c { q } else { p }], -a);
            }
            diag += a;
        }
        out.insert(mesh.cell_id[c], (diag, off));
    }
    out
}

/// Geometry of every rank with the halo filled from the owners.
fn local_geometry(locals: &[LocalMesh]) -> Vec<Geometry> {
    let mut geos: Vec<Geometry> = locals.iter().map(|l| Geometry::compute(&l.mesh)).collect();
    let plans: Vec<&HaloPlan> = locals.iter().map(|l| &l.halo).collect();
    let mut vol: Vec<Vec<f64>> = geos.iter().map(|g| g.cell_volume.clone()).collect();
    HaloPlan::exchange_all(&plans, &mut vol);
    for d in 0..3 {
        let mut x: Vec<Vec<f64>> = geos
            .iter()
            .map(|g| g.cell_centroid.iter().map(|c| c[d]).collect())
            .collect();
        HaloPlan::exchange_all(&plans, &mut x);
        for (g, x) in geos.iter_mut().zip(&x) {
            for (c, v) in g.cell_centroid.iter_mut().zip(x) {
                c[d] = *v;
            }
        }
    }
    for (g, v) in geos.iter_mut().zip(vol) {
        g.cell_volume = v;
    }
    geos
}

/// Compares every rank with the global mesh; returns the number of
/// entities that differ in any bit.
fn differences(global: &Mesh, locals: &[LocalMesh], name: &str) -> usize {
    let ggeo = Geometry::compute(global);
    let grows = rows(global, &ggeo);
    let gface: HashMap<u64, usize> = global
        .face_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let gcell: HashMap<u64, usize> = global
        .cell_id
        .iter()
        .enumerate()
        .map(|(i, &g)| (g, i))
        .collect();
    let mut bad = 0;
    for (l, lgeo) in locals.iter().zip(local_geometry(locals)) {
        let m = &l.mesh;
        assert_eq!(m.face_flip.len(), m.n_faces(), "{name}");
        for f in 0..m.n_faces() {
            let g = gface[&m.face_id[f]];
            let flipped =
                global.cell_id[global.owner[g] as usize] != m.cell_id[m.owner[f] as usize];
            assert_eq!(flipped, m.is_flipped(f), "{name}: flip flag of face {f}");
            let s = ggeo.face_area[g];
            let s = if flipped { [-s[0], -s[1], -s[2]] } else { s };
            if !same3(s, lgeo.face_area[f]) || !same3(ggeo.face_centroid[g], lgeo.face_centroid[f])
            {
                bad += 1;
            }
        }
        // cell face lists in ascending global face position
        for c in 0..m.n_cells {
            let pos: Vec<usize> = m
                .cell_face_iter(c)
                .map(|(f, _)| gface[&m.face_id[f]])
                .collect();
            assert!(
                pos.windows(2).all(|w| w[0] < w[1]),
                "{name}: face order of cell {c}"
            );
        }
        for c in 0..m.n_cells + m.n_halo {
            let g = gcell[&m.cell_id[c]];
            if !same(ggeo.cell_volume[g], lgeo.cell_volume[c])
                || !same3(ggeo.cell_centroid[g], lgeo.cell_centroid[c])
            {
                bad += 1;
            }
        }
        for (id, (diag, off)) in rows(m, &lgeo) {
            let (gd, goff) = &grows[&id];
            if !same(diag, *gd)
                || off.len() != goff.len()
                || off
                    .iter()
                    .any(|(k, v)| goff.get(k).is_none_or(|w| !same(*v, *w)))
            {
                bad += 1;
            }
        }
    }
    bad
}

#[test]
fn local_geometry_is_bitwise_the_one_rank_geometry() {
    for kind in [Kind::Hex, Kind::Tet, Kind::Poly, Kind::Mixed] {
        for dim in [2u8, 3] {
            let mut s = BoxSpec::new(kind, dim, [9, 7, 5]);
            s.perturb = 0.25;
            let m = s.generate().unwrap();
            let g = Geometry::compute(&m);
            for method in methods() {
                for parts in 1usize..=4 {
                    let part = partition(&m, &g, parts, method).unwrap();
                    for layers in [1usize, 2] {
                        let name = format!(
                            "{} {dim}D {} {parts} ranks {layers} layers",
                            kind.name(),
                            method.name()
                        );
                        let locals = decompose(&m, &part, parts, layers).unwrap();
                        assert_eq!(differences(&m, &locals, &name), 0, "{name}");
                        // the library check used by `foehn mesh` agrees
                        for c in canonical_check(&m, &g, &locals) {
                            assert_eq!((c.differing_faces, c.differing_cells), (0, 0), "{name}");
                        }
                    }
                }
            }
        }
    }
}

/// The reference is the mesh handed to the decomposition: after a
/// renumbering (here a random order, the hardest case for face
/// orientation) the ranks reproduce the renumbered one-rank run.
#[test]
fn canonical_after_renumbering() {
    for kind in [Kind::Poly, Kind::Mixed] {
        let mut s = BoxSpec::new(kind, 3, [8, 6, 5]);
        s.perturb = 0.25;
        let mut m = s.generate().unwrap();
        renumber(&mut m, Ordering::Shuffle, 11).unwrap();
        let g = Geometry::compute(&m);
        for parts in [1usize, 3] {
            let part = partition(&m, &g, parts, Method::Rcb).unwrap();
            let locals = decompose(&m, &part, parts, 1).unwrap();
            let name = format!("{} shuffled {parts} ranks", kind.name());
            assert_eq!(differences(&m, &locals, &name), 0, "{name}");
        }
    }
}

/// A mesh that is not a rank's part has no flipped faces, so the geometry
/// of a run without a decomposition does not change.
#[test]
fn undecomposed_meshes_have_no_flips() {
    let m = BoxSpec::new(Kind::Mixed, 3, [4, 4, 3]).generate().unwrap();
    assert!(m.face_flip.is_empty());
    assert!((0..m.n_faces()).all(|f| !m.is_flipped(f)));
}
