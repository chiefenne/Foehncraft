//! The structural model: nodes with six unknowns each, beam elements,
//! supports, springs to ground and loads; assembly of the stiffness K,
//! the consistent mass M and the load vector f on the free unknowns.

use crate::beam::{self, BeamProps};
use crate::shell::{self, ShellProps};
use crate::sparse::{Assembly, SymCsc};

/// Unknowns per node: u, v, w, θx, θy, θz (global axes).
pub const DOFS: usize = 6;

/// Names of the unknowns of a node.
pub const DOF_NAMES: [&str; DOFS] = ["ux", "uy", "uz", "rx", "ry", "rz"];

/// A two-node beam element.
#[derive(Debug, Clone, Copy)]
pub struct Beam {
    /// Its nodes.
    pub nodes: [usize; 2],
    /// Index into [`Model::props`].
    pub props: usize,
    /// Rotation global → local (rows: local x, y, z).
    pub frame: [[f64; 3]; 3],
    /// Length.
    pub len: f64,
}

/// A uniform load per length on a beam element, global axes (N/m).
#[derive(Debug, Clone, Copy)]
pub struct LineLoad {
    /// Element.
    pub element: usize,
    /// Force per length.
    pub q: [f64; 3],
}

/// A shell element: four nodes (MITC4) or three (MITC3+, `tri`), counter-
/// clockwise about its normal.
#[derive(Debug, Clone, Copy)]
pub struct Shell {
    /// Its nodes (a triangle uses the first three).
    pub nodes: [usize; 4],
    /// A three-node element.
    pub tri: bool,
    /// Index into [`Model::shell_props`].
    pub props: usize,
}

/// A load per area on a shell element (N/m²): a vector in global axes
/// plus a pressure along the element's normal.
#[derive(Debug, Clone, Copy)]
pub struct AreaLoad {
    /// Element.
    pub element: usize,
    /// Force per area, global axes.
    pub q: [f64; 3],
    /// Pressure along the normal (positive along it).
    pub pressure: f64,
}

/// The model.
#[derive(Debug, Clone, Default)]
pub struct Model {
    /// Node coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Section and material sets.
    pub props: Vec<BeamProps>,
    /// Beam elements.
    pub beams: Vec<Beam>,
    /// Fixed unknowns (`DOFS` per node).
    pub fixed: Vec<bool>,
    /// Springs to ground on single unknowns (unknown, stiffness).
    pub springs: Vec<(usize, f64)>,
    /// Point masses (node, mass): translational only.
    pub masses: Vec<(usize, f64)>,
    /// Nodal forces and moments (`DOFS` per node).
    pub forces: Vec<f64>,
    /// Distributed loads.
    pub line_loads: Vec<LineLoad>,
    /// Shell thickness and material sets.
    pub shell_props: Vec<ShellProps>,
    /// Shell elements.
    pub shells: Vec<Shell>,
    /// Loads per area on shells.
    pub area_loads: Vec<AreaLoad>,
    /// Prescribed values of fixed unknowns (unknown, value); a fixed
    /// unknown not listed is zero.
    pub prescribed: Vec<(usize, f64)>,
}

/// The assembled system on the free unknowns.
#[derive(Debug, Clone)]
pub struct System {
    /// Stiffness on the free unknowns.
    pub k: SymCsc,
    /// Consistent mass on the free unknowns.
    pub m: SymCsc,
    /// Loads on the free unknowns.
    pub f: Vec<f64>,
    /// Free index of each unknown of the model (`None`: fixed).
    pub map: Vec<Option<usize>>,
    /// Model unknown of each free index.
    pub free: Vec<usize>,
    /// Full stiffness (every unknown), for the reactions.
    pub k_all: SymCsc,
    /// Full loads (every unknown).
    pub f_all: Vec<f64>,
    /// Values of every unknown at the fixed ones (zero elsewhere).
    pub prescribed: Vec<f64>,
}

impl System {
    /// A free vector expanded to every unknown (the prescribed values,
    /// zero by default, at fixed ones).
    pub fn expand(&self, x: &[f64]) -> Vec<f64> {
        let mut out = self.prescribed.clone();
        for (k, &d) in self.free.iter().enumerate() {
            out[d] = x[k];
        }
        out
    }

    /// A free vector expanded with zeros at the fixed unknowns (a mode).
    pub fn expand_zero(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.map.len()];
        for (k, &d) in self.free.iter().enumerate() {
            out[d] = x[k];
        }
        out
    }

    /// A full vector restricted to the free unknowns.
    pub fn restrict(&self, x: &[f64]) -> Vec<f64> {
        self.free.iter().map(|&d| x[d]).collect()
    }
}

impl Shell {
    /// The nodes in use.
    pub fn used(&self) -> &[usize] {
        if self.tri {
            &self.nodes[..3]
        } else {
            &self.nodes
        }
    }
}

impl Model {
    /// The unit normal of shell element `e` at its centre.
    pub fn shell_normal(&self, e: usize) -> [f64; 3] {
        let s = &self.shells[e];
        if s.tri {
            let x = [
                self.nodes[s.nodes[0]],
                self.nodes[s.nodes[1]],
                self.nodes[s.nodes[2]],
            ];
            let n = beam::cross(beam::sub(x[1], x[0]), beam::sub(x[2], x[0]));
            beam::scale(n, 1.0 / beam::norm(n))
        } else {
            shell::centre_normal(&s.nodes.map(|n| self.nodes[n]))
        }
    }

    /// Adds a structured grid of `na × nb` cells on the map `at(a, b)`,
    /// each cut into two three-node shells along the diagonal; as
    /// [`Model::add_shell_grid`].
    pub fn add_tri_grid(
        &mut self,
        at: impl Fn(f64, f64) -> [f64; 3],
        na: usize,
        nb: usize,
        props: usize,
    ) -> Vec<Vec<usize>> {
        let first = self.shells.len();
        let grid = self.add_shell_grid(at, na, nb, props);
        let quads: Vec<Shell> = self.shells.drain(first..).collect();
        for (c, q) in quads.iter().enumerate() {
            let n = q.nodes;
            // alternate the diagonal so the mesh has no preferred direction
            let (t1, t2) = if (c / nb + c % nb).is_multiple_of(2) {
                ([n[0], n[1], n[2]], [n[0], n[2], n[3]])
            } else {
                ([n[0], n[1], n[3]], [n[1], n[2], n[3]])
            };
            for t in [t1, t2] {
                self.shells.push(Shell {
                    nodes: [t[0], t[1], t[2], usize::MAX],
                    tri: true,
                    props,
                });
            }
        }
        grid
    }

    /// Adds a node and returns its index.
    pub fn add_node(&mut self, x: [f64; 3]) -> usize {
        self.nodes.push(x);
        self.fixed.extend([false; DOFS]);
        self.forces.extend([0.0; DOFS]);
        self.nodes.len() - 1
    }

    /// Adds a beam from node `a` to node `b`; `y_hint` sets the local y
    /// axis (`None`: [`beam::default_y`]). `None` if the beam has no
    /// length or the hint is parallel to it.
    pub fn add_beam(
        &mut self,
        a: usize,
        b: usize,
        props: usize,
        y_hint: Option<[f64; 3]>,
    ) -> Option<usize> {
        let (xa, xb) = (self.nodes[a], self.nodes[b]);
        let hint = y_hint.unwrap_or_else(|| beam::default_y(beam::sub(xb, xa)));
        let (frame, len) = beam::frame(xa, xb, hint)?;
        self.beams.push(Beam {
            nodes: [a, b],
            props,
            frame,
            len,
        });
        Some(self.beams.len() - 1)
    }

    /// Adds a structured grid of `na × nb` shell elements on the map
    /// `at(a, b)` of the unit square (a, b ∈ [0, 1]); returns the node
    /// indices by grid point, `grid[i][j]` at a = i/na, b = j/nb. The
    /// elements' normals are along ∂x/∂a × ∂x/∂b.
    pub fn add_shell_grid(
        &mut self,
        at: impl Fn(f64, f64) -> [f64; 3],
        na: usize,
        nb: usize,
        props: usize,
    ) -> Vec<Vec<usize>> {
        let grid: Vec<Vec<usize>> = (0..=na)
            .map(|i| {
                (0..=nb)
                    .map(|j| self.add_node(at(i as f64 / na as f64, j as f64 / nb as f64)))
                    .collect()
            })
            .collect();
        for i in 0..na {
            for j in 0..nb {
                self.shells.push(Shell {
                    nodes: [
                        grid[i][j],
                        grid[i + 1][j],
                        grid[i + 1][j + 1],
                        grid[i][j + 1],
                    ],
                    tri: false,
                    props,
                });
            }
        }
        grid
    }

    /// Number of unknowns.
    pub fn dofs(&self) -> usize {
        self.nodes.len() * DOFS
    }

    /// The unknowns that no element, spring or mass touches (they would
    /// make K singular): fixed by [`Model::assemble`].
    fn untouched(&self) -> Vec<bool> {
        let mut used = vec![false; self.dofs()];
        for b in &self.beams {
            for &n in &b.nodes {
                used[DOFS * n..DOFS * n + DOFS]
                    .iter_mut()
                    .for_each(|u| *u = true);
            }
        }
        for s in &self.shells {
            for &n in s.used() {
                used[DOFS * n..DOFS * n + DOFS]
                    .iter_mut()
                    .for_each(|u| *u = true);
            }
        }
        for &(d, _) in &self.springs {
            used[d] = true;
        }
        used.iter().map(|u| !u).collect()
    }

    /// Assembles K, M and f on the free unknowns.
    pub fn assemble(&self) -> System {
        let n = self.dofs();
        let untouched = self.untouched();
        let mut map = vec![None; n];
        let mut free = Vec::new();
        for d in 0..n {
            if !self.fixed[d] && !untouched[d] {
                map[d] = Some(free.len());
                free.push(d);
            }
        }
        let nf = free.len();
        let mut k = Assembly::new(nf);
        let mut m = Assembly::new(nf);
        let mut ka = Assembly::new(n);
        let mut f_all = self.forces.clone();
        for (e, b) in self.beams.iter().enumerate() {
            let p = &self.props[b.props];
            let (kl, ml) = beam::local_matrices(p, b.len);
            let kg = beam::to_global(&b.frame, &kl);
            let mg = beam::to_global(&b.frame, &ml);
            let all: Vec<usize> = b
                .nodes
                .iter()
                .flat_map(|&nd| (0..DOFS).map(move |d| DOFS * nd + d))
                .collect();
            let fr: Vec<Option<usize>> = all.iter().map(|&d| map[d]).collect();
            let al: Vec<Option<usize>> = all.iter().map(|&d| Some(d)).collect();
            k.add_block(&fr, &kg);
            m.add_block(&fr, &mg);
            ka.add_block(&al, &kg);
            for l in self.line_loads.iter().filter(|l| l.element == e) {
                let ql = beam::rotate(&b.frame, l.q);
                let fl = beam::local_line_load(p, b.len, ql);
                let fg = beam::vec_to_global(&b.frame, &fl);
                for (i, &d) in all.iter().enumerate() {
                    f_all[d] += fg[i];
                }
            }
        }
        let dirs = if self.shells.is_empty() {
            Vec::new()
        } else {
            self.directors()
        };
        for e in 0..self.shells.len() {
            let sh = self.shells[e];
            let p = &self.shell_props[sh.props];
            let normal = self.shell_normal(e);
            let used = sh.used();
            let (kg, mg) = if sh.tri {
                let geo = shell::TriGeometry {
                    x: [
                        self.nodes[used[0]],
                        self.nodes[used[1]],
                        self.nodes[used[2]],
                    ],
                    v: [dirs[used[0]], dirs[used[1]], dirs[used[2]]],
                    h: 0.5 * p.thickness,
                };
                shell::mitc3plus(&geo, p)
            } else {
                shell::mitc4(&self.shell_geometry(e, &dirs), p)
            };
            let all: Vec<usize> = used
                .iter()
                .flat_map(|&nd| (0..DOFS).map(move |d| DOFS * nd + d))
                .collect();
            let fr: Vec<Option<usize>> = all.iter().map(|&d| map[d]).collect();
            let al: Vec<Option<usize>> = all.iter().map(|&d| Some(d)).collect();
            k.add_block(&fr, &kg);
            m.add_block(&fr, &mg);
            ka.add_block(&al, &kg);
            for l in self.area_loads.iter().filter(|l| l.element == e) {
                let q = [
                    l.q[0] + l.pressure * normal[0],
                    l.q[1] + l.pressure * normal[1],
                    l.q[2] + l.pressure * normal[2],
                ];
                let fg: Vec<f64> = if sh.tri {
                    let x = [
                        self.nodes[used[0]],
                        self.nodes[used[1]],
                        self.nodes[used[2]],
                    ];
                    shell::tri_area_load(&x, q).0.to_vec()
                } else {
                    shell::area_load(&self.shell_geometry(e, &dirs), q)
                        .0
                        .to_vec()
                };
                for (i, &d) in all.iter().enumerate() {
                    f_all[d] += fg[i];
                }
            }
        }
        for &(d, s) in &self.springs {
            if let Some(i) = map[d] {
                k.add(i, i, s);
            }
            ka.add(d, d, s);
        }
        for &(nd, mass) in &self.masses {
            for d in 0..3 {
                if let Some(i) = map[DOFS * nd + d] {
                    m.add(i, i, mass);
                }
            }
        }
        let mut up = vec![0.0; n];
        for &(d, v) in &self.prescribed {
            up[d] = v;
        }
        let kup = if self.prescribed.is_empty() {
            vec![0.0; n]
        } else {
            ka.finish().apply(&up)
        };
        let f = free.iter().map(|&d| f_all[d] - kup[d]).collect();
        System {
            k: k.finish(),
            m: m.finish(),
            f,
            map,
            free,
            k_all: ka.finish(),
            f_all,
            prescribed: up,
        }
    }

    /// The unit directors of the shell nodes: the mean of the normals of
    /// the elements around each node (zero where there is no shell).
    pub fn directors(&self) -> Vec<[f64; 3]> {
        let mut v = vec![[0.0; 3]; self.nodes.len()];
        for (e, s) in self.shells.iter().enumerate() {
            let n = self.shell_normal(e);
            for &k in s.used() {
                for d in 0..3 {
                    v[k][d] += n[d];
                }
            }
        }
        for a in &mut v {
            let l = beam::norm(*a);
            if l > 0.0 {
                *a = beam::scale(*a, 1.0 / l);
            }
        }
        v
    }

    /// The geometry of shell element `e` with the nodal directors `v`.
    pub fn shell_geometry(&self, e: usize, v: &[[f64; 3]]) -> shell::Geometry {
        let s = &self.shells[e];
        shell::Geometry {
            x: s.nodes.map(|n| self.nodes[n]),
            v: s.nodes.map(|n| v[n]),
            h: 0.5 * self.shell_props[s.props].thickness,
        }
    }

    /// Total translational mass (elements and point masses).
    pub fn total_mass(&self) -> f64 {
        let beams: f64 = self
            .beams
            .iter()
            .map(|b| {
                let p = &self.props[b.props];
                p.density * p.area * b.len
            })
            .sum();
        let shells: f64 = (0..self.shells.len())
            .map(|e| {
                let s = &self.shells[e];
                let p = &self.shell_props[s.props];
                let area = if s.tri {
                    let x = [
                        self.nodes[s.nodes[0]],
                        self.nodes[s.nodes[1]],
                        self.nodes[s.nodes[2]],
                    ];
                    shell::tri_area_load(&x, [0.0; 3]).1
                } else {
                    let geo = shell::Geometry {
                        x: s.nodes.map(|n| self.nodes[n]),
                        v: [[0.0; 3]; 4],
                        h: 0.0,
                    };
                    shell::area_load(&geo, [0.0; 3]).1
                };
                p.density * p.thickness * area
            })
            .sum();
        beams + shells + self.masses.iter().map(|m| m.1).sum::<f64>()
    }

    /// The node closest to `x` and its distance.
    pub fn nearest(&self, x: [f64; 3]) -> Option<(usize, f64)> {
        self.nodes
            .iter()
            .enumerate()
            .map(|(i, p)| (i, beam::norm(beam::sub(*p, x))))
            .min_by(|a, b| a.1.total_cmp(&b.1))
    }

    /// Size of the model (diagonal of the bounding box).
    pub fn size(&self) -> f64 {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &self.nodes {
            for d in 0..3 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
        if self.nodes.is_empty() {
            0.0
        } else {
            beam::norm(beam::sub(hi, lo))
        }
    }
}
