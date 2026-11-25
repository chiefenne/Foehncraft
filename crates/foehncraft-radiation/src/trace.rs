//! Paths of chosen bundles: where a bundle leaves each cell, meets a wall
//! and ends, with the quanta it still carries there (R1a, the furnace
//! showcase).
//!
//! [`Domain::trace`] emits the chosen bundles exactly as [`Domain::solve`]
//! does (the same emitters, weights and Sobol points for the same
//! [`Settings`]) and moves them with the same kernel body, whose tracing
//! hook ([`Tally::trace`]) does nothing in a call. A traced path is
//! therefore the path that bundle takes in the call, and the quanta it
//! leaves in the cells and on the walls are the bundle's share of the
//! call's tallies, to the quantum.

use crate::solve::{Domain, Settings, Wall};
use foehncraft_kernels::radiation::{Fate, Tally, Trace, transport};
use rayon::prelude::*;
use std::cell::{Cell, RefCell};

/// What emits a traced bundle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Source {
    /// A cell of the gas (index in the mesh).
    Cell(usize),
    /// A boundary face (index from the first boundary face).
    Wall(usize),
}

/// A point of a path: where the bundle stands and the quanta it carries
/// on from there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathPoint {
    /// What happens there.
    pub at: Trace,
    /// Position (m).
    pub x: [f64; 3],
    /// Quanta the bundle carries on (zero where it ends by absorption).
    pub w: i64,
}

/// The path of one bundle and the quanta it leaves.
#[derive(Debug, Clone, PartialEq)]
pub struct Path {
    /// What emitted it.
    pub source: Source,
    /// Index of the bundle among the emitter's bundles.
    pub index: u32,
    /// Quanta at emission.
    pub w0: i64,
    /// Point of emission (m).
    pub start: [f64; 3],
    /// The points after the start, in order; the last one is where it
    /// ends ([`Trace::Absorbed`], [`Trace::Escaped`] or [`Trace::Lost`]).
    pub points: Vec<PathPoint>,
    /// Quanta absorbed in each cell it crossed, in order (cell, quanta).
    pub cells: Vec<(u32, i64)>,
    /// Each wall it met (boundary face from the first boundary face,
    /// incident quanta, absorbed quanta).
    pub walls: Vec<(u32, i64, i64)>,
    /// Quanta that left through an open boundary.
    pub escaped: i64,
    /// Quanta Russian roulette removed (negative: added).
    pub roulette: i64,
    /// Quanta lost to degenerate geometry.
    pub lost: i64,
}

/// The traced bundles of a call.
#[derive(Debug, Clone, PartialEq)]
pub struct Traced {
    /// Power of one quantum (W), as in the call.
    pub quantum: f64,
    /// The paths, emitter by emitter (cells first, then boundary faces),
    /// bundle by bundle.
    pub paths: Vec<Path>,
}

/// Records one bundle.
#[derive(Default)]
struct Recorder {
    points: RefCell<Vec<PathPoint>>,
    cells: RefCell<Vec<(u32, i64)>>,
    walls: RefCell<Vec<(u32, i64, i64)>>,
    escaped: Cell<i64>,
    roulette: Cell<i64>,
    lost: Cell<i64>,
}

impl Tally for Recorder {
    fn cell(&self, _batch: u32, cell: usize, q: i64) {
        self.cells.borrow_mut().push((cell as u32, q));
    }
    fn wall(&self, _batch: u32, face: usize, incident: i64, absorbed: i64) {
        self.walls
            .borrow_mut()
            .push((face as u32, incident, absorbed));
    }
    fn escaped(&self, q: i64) {
        self.escaped.set(self.escaped.get() + q);
    }
    fn roulette(&self, q: i64) {
        self.roulette.set(self.roulette.get() + q);
    }
    fn lost(&self, q: i64) {
        self.lost.set(self.lost.get() + q);
    }
    fn trace(&self, at: Trace, x: [f64; 3], w: i64) {
        self.points.borrow_mut().push(PathPoint { at, x, w });
    }
}

impl Domain<'_> {
    /// Traces the first `pick(source)` bundles of every emitter of the call
    /// `solve(temperature, kappa, walls, s)` would make (zero: none; more
    /// than the emitter has: all). Runs on the CPU threads of one process;
    /// a rank's part of a mesh is refused, since a path would end at the
    /// processor faces.
    pub fn trace(
        &self,
        temperature: &[f64],
        kappa: &[f64],
        walls: &[Wall],
        s: &Settings,
        pick: impl Fn(Source) -> u32,
    ) -> anyhow::Result<Traced> {
        anyhow::ensure!(
            !self.on_ranks(),
            "bundle paths are traced on the whole mesh, not on a rank's part"
        );
        let mesh = self.mesh();
        anyhow::ensure!(
            temperature.len() == mesh.n_cells && kappa.len() == mesh.n_cells,
            "one temperature and one absorption coefficient per cell"
        );
        anyhow::ensure!(
            walls.len() == mesh.boundaries.len(),
            "one wall per boundary"
        );
        let Some(setup) = self.setup(temperature, kappa, walls, s) else {
            return Ok(Traced {
                quantum: 0.0,
                paths: Vec::new(),
            });
        };
        let medium = self.medium(kappa, &setup.kinds, &setup.emissivity);
        let wm = self.locator().walk_mesh();
        // emission here (cheap), the transport on the threads
        let mut bundles = Vec::new();
        for e in &setup.emitters {
            let source = if e.face == u32::MAX {
                Source::Cell(e.cell as usize)
            } else {
                Source::Wall(e.face as usize)
            };
            for k in 0..pick(source).min(e.bundles) {
                bundles.push((source, k, self.emit(e, k, s, &wm)));
            }
        }
        let paths = bundles
            .into_par_iter()
            .map(|(source, index, mut b)| {
                let (start, w0) = (b.walk.origin, b.w);
                let r = Recorder::default();
                let fate = transport(&wm, &medium, &mut b, s.call, &r);
                debug_assert_eq!(fate, Fate::Done, "one process has no processor faces");
                Path {
                    source,
                    index,
                    w0,
                    start,
                    points: r.points.into_inner(),
                    cells: r.cells.into_inner(),
                    walls: r.walls.into_inner(),
                    escaped: r.escaped.get(),
                    roulette: r.roulette.get(),
                    lost: r.lost.get(),
                }
            })
            .collect();
        Ok(Traced {
            quantum: setup.quantum,
            paths,
        })
    }
}
