//! The structural model: nodes with six unknowns each, beam elements,
//! supports, springs to ground and loads; assembly of the stiffness K,
//! the consistent mass M and the load vector f on the free unknowns.

use crate::beam::{self, BeamProps};
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
}

impl System {
    /// A free vector expanded to every unknown (zeros at fixed ones).
    pub fn expand(&self, x: &[f64]) -> Vec<f64> {
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

impl Model {
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
        let f = free.iter().map(|&d| f_all[d]).collect();
        System {
            k: k.finish(),
            m: m.finish(),
            f,
            map,
            free,
            k_all: ka.finish(),
            f_all,
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
        beams + self.masses.iter().map(|m| m.1).sum::<f64>()
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
