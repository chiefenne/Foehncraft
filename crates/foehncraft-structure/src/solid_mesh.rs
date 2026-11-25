//! Meshes of quadratic solids: read from a Gmsh `.msh` 4.1 file (second
//! order, with physical groups) or generated as a block.
//!
//! - Gmsh: ten-node tetrahedra (type 11) and twenty-node hexahedra (type
//!   17, `Mesh.SecondOrderIncomplete = 1`) are the elements; six-node
//!   triangles (type 9) and eight-node quadrilaterals (type 16) on surfaces
//!   with a physical group are the named surfaces for supports and
//!   pressures. Volume physical groups name the parts of the mesh.
//! - Block: a box divided into hexahedra (`hex20`) or six tetrahedra per
//!   hexahedron (`tet10`, the Kuhn split, conforming); its sides are the
//!   surfaces `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`.
//!
//! Every surface face is stored with the node order of the element face it
//! lies on: counter-clockwise seen from outside, so its normal points out
//! of the solid.

use crate::solid::{self, Kind};
use foehncraft_mesh::gmsh::RawGmsh;
use std::collections::{BTreeMap, HashMap};

/// A mesh of quadratic solid elements.
#[derive(Debug, Clone, Default)]
pub struct SolidMesh {
    /// Node coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Elements: kind and nodes (VTK order).
    pub elements: Vec<(Kind, Vec<usize>)>,
    /// The volume group of each element: an index into `volumes`.
    pub volume_of: Vec<usize>,
    /// The names of each volume group (its physical groups).
    pub volumes: Vec<Vec<String>>,
    /// Named surfaces: faces (6 or 8 nodes) pointing out of the solid.
    pub surfaces: BTreeMap<String, Vec<Vec<usize>>>,
    /// Where it came from (`MSH 4.1 ASCII`, `block`).
    pub source: String,
}

impl SolidMesh {
    /// Keeps only the elements whose volume group has one of `names`
    /// (every element if `names` is empty), then drops the nodes and
    /// surface faces no element uses. Unknown names are an error.
    pub fn select(mut self, names: &[String]) -> Result<Self, String> {
        if names.is_empty() {
            return Ok(self);
        }
        for n in names {
            if !self.volumes.iter().any(|v| v.contains(n)) {
                let mut all: Vec<&String> = self.volumes.iter().flatten().collect();
                all.sort();
                all.dedup();
                return Err(format!(
                    "no volume `{n}` in the mesh; it has {}",
                    all.iter()
                        .map(|s| format!("`{s}`"))
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
            }
        }
        let keep: Vec<bool> = self
            .volume_of
            .iter()
            .map(|&v| self.volumes[v].iter().any(|n| names.contains(n)))
            .collect();
        let mut elements = Vec::new();
        let mut volume_of = Vec::new();
        for (i, e) in self.elements.into_iter().enumerate() {
            if keep[i] {
                elements.push(e);
                volume_of.push(self.volume_of[i]);
            }
        }
        self.elements = elements;
        self.volume_of = volume_of;
        // faces of the kept elements
        let mut faces = std::collections::HashSet::new();
        for (k, e) in &self.elements {
            for f in k.faces() {
                faces.insert(face_key(&f.iter().map(|&a| e[a]).collect::<Vec<_>>()));
            }
        }
        for v in self.surfaces.values_mut() {
            v.retain(|f| faces.contains(&face_key(f)));
        }
        self.surfaces.retain(|_, v| !v.is_empty());
        Ok(self.compact())
    }

    /// Drops the nodes no element uses and numbers the rest in order.
    fn compact(mut self) -> Self {
        let mut map = vec![usize::MAX; self.nodes.len()];
        let mut nodes = Vec::new();
        for (_, e) in &mut self.elements {
            for n in e.iter_mut() {
                if map[*n] == usize::MAX {
                    map[*n] = nodes.len();
                    nodes.push(self.nodes[*n]);
                }
                *n = map[*n];
            }
        }
        for v in self.surfaces.values_mut() {
            for f in v.iter_mut() {
                for n in f.iter_mut() {
                    *n = map[*n];
                }
            }
        }
        self.nodes = nodes;
        self
    }

    /// Points the surface faces out of the solid: each takes the node
    /// order of the element face it lies on. A face on no element is an
    /// error.
    fn orient(&mut self) -> Result<(), String> {
        let mut owner: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
        for (k, e) in &self.elements {
            for f in k.faces() {
                let nodes: Vec<usize> = f.iter().map(|&a| e[a]).collect();
                owner.entry(face_key(&nodes)).or_insert(nodes);
            }
        }
        for (name, faces) in &mut self.surfaces {
            for f in faces.iter_mut() {
                match owner.get(&face_key(f)) {
                    Some(o) if o.len() == f.len() => *f = o.clone(),
                    _ => {
                        return Err(format!(
                            "surface `{name}` has a face that is not a face of a solid element"
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

/// The sorted corners of a face: its identity.
fn face_key(f: &[usize]) -> Vec<usize> {
    let c = if f.len() == 8 { 4 } else { 3 };
    let mut k = f[..c].to_vec();
    k.sort_unstable();
    k
}

/// The solid mesh of a Gmsh file read by
/// [`foehncraft_mesh::gmsh::parse_raw`].
pub fn from_gmsh(raw: &RawGmsh) -> Result<SolidMesh, String> {
    let idx = raw.indexer();
    let mut m = SolidMesh {
        nodes: raw.nodes.clone(),
        source: raw.format.clone(),
        ..Default::default()
    };
    let mut groups: HashMap<i64, usize> = HashMap::new();
    for b in &raw.blocks {
        match (b.dim, b.etype) {
            (3, 11) | (3, 17) => {
                let (kind, perm): (Kind, &[usize]) = if b.etype == 11 {
                    (Kind::Tet10, &solid::GMSH_TET10)
                } else {
                    (Kind::Hex20, &solid::GMSH_HEX20)
                };
                let g = *groups.entry(b.entity).or_insert_with(|| {
                    m.volumes.push(raw.entity_names(3, b.entity));
                    m.volumes.len() - 1
                });
                for e in b.nodes.chunks(b.nodes_per) {
                    let v = perm
                        .iter()
                        .map(|&p| idx(e[p]).map(|i| i as usize))
                        .collect::<Result<Vec<_>, _>>()?;
                    m.elements.push((kind, v));
                    m.volume_of.push(g);
                }
            }
            (3, 4) | (3, 5) | (3, 6) | (3, 7) => {
                return Err(
                    "the solid mesh has linear elements; mesh with second order \
                     (Gmsh option -order 2, Mesh.ElementOrder = 2)"
                        .into(),
                );
            }
            (3, 12) => {
                return Err(
                    "the solid mesh has 27-node hexahedra; use the 20-node ones \
                     (Mesh.SecondOrderIncomplete = 1)"
                        .into(),
                );
            }
            (3, t) => {
                return Err(format!(
                    "element type {t} is not a solid element of the structural solver \
                     (ten-node tetrahedra and twenty-node hexahedra)"
                ));
            }
            (2, 9) | (2, 16) => {
                let names = raw.entity_names(2, b.entity);
                for e in b.nodes.chunks(b.nodes_per) {
                    let v = e
                        .iter()
                        .map(|&t| idx(t).map(|i| i as usize))
                        .collect::<Result<Vec<_>, _>>()?;
                    for n in &names {
                        m.surfaces.entry(n.clone()).or_default().push(v.clone());
                    }
                }
            }
            (2, 10) => {
                return Err(
                    "the solid mesh has 9-node quadrilaterals; use the 8-node ones \
                     (Mesh.SecondOrderIncomplete = 1)"
                        .into(),
                );
            }
            // linear surface elements (a mesh of first order) and the rest
            // are refused by the volume check above or are not needed
            _ => {}
        }
    }
    if m.elements.is_empty() {
        return Err(
            "the mesh has no solid elements (ten-node tetrahedra or twenty-node hexahedra)".into(),
        );
    }
    let mut m = m.compact();
    m.orient()?;
    Ok(m)
}

/// A block `origin + [0, size]` with `divisions` hexahedra along each
/// axis, as twenty-node hexahedra or as six ten-node tetrahedra each.
pub fn block(origin: [f64; 3], size: [f64; 3], divisions: [usize; 3], kind: Kind) -> SolidMesh {
    let [nx, ny, nz] = divisions;
    let mut m = SolidMesh {
        source: "block".into(),
        volumes: vec![vec!["block".into()]],
        ..Default::default()
    };
    // nodes on the lattice of half cells, created when first used
    let mut at: HashMap<[usize; 3], usize> = HashMap::new();
    let mut lattice: Vec<[usize; 3]> = Vec::new();
    let mut node = |m: &mut SolidMesh, p: [usize; 3]| -> usize {
        *at.entry(p).or_insert_with(|| {
            m.nodes.push([
                origin[0] + size[0] * p[0] as f64 / (2 * nx) as f64,
                origin[1] + size[1] * p[1] as f64 / (2 * ny) as f64,
                origin[2] + size[2] * p[2] as f64 / (2 * nz) as f64,
            ]);
            lattice.push(p);
            m.nodes.len() - 1
        })
    };
    let lat = |base: [usize; 3], off: [f64; 3]| -> [usize; 3] {
        [
            base[0] + off[0] as usize,
            base[1] + off[1] as usize,
            base[2] + off[2] as usize,
        ]
    };
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let base = [2 * i, 2 * j, 2 * k];
                match kind {
                    Kind::Hex20 => {
                        let v: Vec<usize> = Kind::Hex20
                            .natural_nodes()
                            .iter()
                            .map(|c| node(&mut m, lat(base, [c[0] + 1.0, c[1] + 1.0, c[2] + 1.0])))
                            .collect();
                        m.elements.push((Kind::Hex20, v));
                        m.volume_of.push(0);
                    }
                    Kind::Tet10 => {
                        // the Kuhn split: paths from corner 000 to 111 along
                        // the axes in each of the six orders
                        for perm in [
                            [0, 1, 2],
                            [0, 2, 1],
                            [1, 0, 2],
                            [1, 2, 0],
                            [2, 0, 1],
                            [2, 1, 0],
                        ] {
                            let mut c = [[0usize; 3]; 4];
                            for s in 1..4 {
                                c[s] = c[s - 1];
                                c[s][perm[s - 1]] = 2;
                            }
                            // positive volume: swap two corners if needed
                            let d = |a: usize, b: usize, x: usize| c[a][x] as f64 - c[b][x] as f64;
                            let e1 = [d(1, 0, 0), d(1, 0, 1), d(1, 0, 2)];
                            let e2 = [d(2, 0, 0), d(2, 0, 1), d(2, 0, 2)];
                            let e3 = [d(3, 0, 0), d(3, 0, 1), d(3, 0, 2)];
                            if crate::beam::dot3(crate::beam::cross(e1, e2), e3) < 0.0 {
                                c.swap(1, 2);
                            }
                            let mut v: Vec<usize> = c
                                .iter()
                                .map(|p| node(&mut m, [base[0] + p[0], base[1] + p[1], base[2] + p[2]]))
                                .collect();
                            for [a, b] in solid::TET_EDGES {
                                let p = [
                                    base[0] + (c[a][0] + c[b][0]) / 2,
                                    base[1] + (c[a][1] + c[b][1]) / 2,
                                    base[2] + (c[a][2] + c[b][2]) / 2,
                                ];
                                v.push(node(&mut m, p));
                            }
                            m.elements.push((Kind::Tet10, v));
                            m.volume_of.push(0);
                        }
                    }
                }
            }
        }
    }
    // the sides
    let top = [2 * nx, 2 * ny, 2 * nz];
    const SIDES: [[&str; 2]; 3] = [["xmin", "xmax"], ["ymin", "ymax"], ["zmin", "zmax"]];
    for (kind, e) in &m.elements {
        for f in kind.faces() {
            let nodes: Vec<usize> = f.iter().map(|&a| e[a]).collect();
            let corners = &nodes[..if nodes.len() == 8 { 4 } else { 3 }];
            for d in 0..3 {
                for (s, val) in [0, top[d]].into_iter().enumerate() {
                    if corners.iter().all(|&n| lattice[n][d] == val) {
                        m.surfaces
                            .entry(SIDES[d][s].to_string())
                            .or_default()
                            .push(nodes.clone());
                    }
                }
            }
        }
    }
    m
}
