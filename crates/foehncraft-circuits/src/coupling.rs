//! A circuit whose tube heat comes from outside, segment by segment: the
//! water side of the walls of a 3D furnace (steam furnace stage SF7c).
//!
//! The 3D side owns the geometry: it maps its wall faces to the segments
//! of the tubes, sums the heat each segment receives and calls
//! [`Coupled::solve`] with those heats. The circuit answers per segment
//! with the bulk temperature of the water and steam and the heat transfer
//! coefficient of the inner wall, which the 3D side turns into a Robin
//! condition h (T_water − T_wall) on its faces. The Newton solve starts
//! from the unknowns of the last solve, so the repeated solves of a
//! coupled run take one to three iterations each.

use crate::case::Case;
use crate::network::{Model, Segment};
use crate::output::Solution;
use crate::{Error, codes};
use std::f64::consts::PI;

/// The heat flux (W/m²) at which the heat transfer coefficient of a
/// segment is taken when its own heat flux is smaller: in boiling the
/// coefficient grows with the heat flux, and a wall that is not heated yet
/// (the first solve, a cold corner) still needs a finite one.
pub const Q_FLOOR: f64 = 1e3;

/// The water side of one segment for the 3D wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterSide {
    /// Bulk temperature (K).
    pub t_bulk: f64,
    /// Heat transfer coefficient of the inner wall q/(T_w − T_b)
    /// (W/(m² K)), at the segment's heat flux or [`Q_FLOOR`].
    pub htc: f64,
    /// Inner wall temperature at the segment's heat flux (K).
    pub t_wall: f64,
    /// Heat flux on the inner wall (W/m²).
    pub q: f64,
    /// CHF ratio q_CHF/q (NaN where not computed).
    pub chf_ratio: f64,
}

/// A circuit driven by the heat of the segments of some of its tubes.
#[derive(Clone, Debug)]
pub struct Coupled {
    /// The case.
    pub case: Case,
    /// The model; the heat flux of the coupled tubes is replaced at every
    /// solve.
    pub model: Model,
    /// The tubes whose heat comes from outside.
    pub coupled: Vec<bool>,
    /// The unknowns of the last solve (the warm start).
    x: Option<Vec<f64>>,
    /// The last solution.
    pub solution: Option<Solution>,
    /// Newton iterations of the last solve.
    pub iterations: usize,
}

impl Coupled {
    /// The circuit of a case; `tubes` names the tubes whose heat comes from
    /// outside (CIR001 for a name that is not a tube).
    pub fn new(case: Case, tubes: &[&str]) -> Result<Self, Error> {
        let model = Model::new(&case)?;
        let mut coupled = vec![false; model.tubes.len()];
        for name in tubes {
            let i = model
                .tubes
                .iter()
                .position(|t| t.name == *name)
                .ok_or_else(|| {
                    Error::new(codes::CASE, format!("no tube `{name}` in the circuit"))
                })?;
            coupled[i] = true;
        }
        Ok(Self {
            case,
            model,
            coupled,
            x: None,
            solution: None,
            iterations: 0,
        })
    }

    /// Index of a tube.
    pub fn tube(&self, name: &str) -> Option<usize> {
        self.model.tubes.iter().position(|t| t.name == name)
    }

    /// Segments of tube `i`.
    pub fn segments(&self, i: usize) -> usize {
        self.model.tubes[i].q.len()
    }

    /// The inner wall area (m²) of one segment of tube `i`, all its
    /// parallel tubes together.
    pub fn inner_area(&self, i: usize) -> f64 {
        let t = &self.model.tubes[i];
        t.count * PI * t.d * t.dz()
    }

    /// The heat (W) of each segment of tube `i`, all parallel tubes
    /// together, of the current model.
    pub fn heat(&self, i: usize) -> Vec<f64> {
        let a = self.inner_area(i);
        self.model.tubes[i].q.iter().map(|q| q * a).collect()
    }

    /// Solve with the heat (W, all parallel tubes of a segment together)
    /// of every segment of the coupled tubes, `heat[tube][segment]`; the
    /// tubes that are not coupled keep the heat of the case and may have
    /// an empty row. Returns the water side of every segment of every
    /// tube, `[tube][segment]` in the order from `from` to `to`.
    pub fn solve(&mut self, heat: &[Vec<f64>]) -> Result<Vec<Vec<WaterSide>>, Error> {
        let mut model = self.model.clone();
        for (i, t) in model.tubes.iter_mut().enumerate() {
            if !self.coupled[i] {
                continue;
            }
            let row = heat
                .get(i)
                .filter(|r| r.len() == t.q.len())
                .ok_or_else(|| {
                    Error::new(
                        codes::CASE,
                        format!("tube `{}`: one heat per segment is needed", t.name),
                    )
                })?;
            let a = t.count * PI * t.d * t.dz();
            for (q, w) in t.q.iter_mut().zip(row) {
                *q = w / a;
            }
        }
        // the warm start, and the start of the case when it fails (a loop
        // that stood still without heat has a singular Jacobian)
        let (x, marches, history) = match self.x.clone() {
            Some(x0) => match model.solve_from(x0) {
                Ok(s) => s,
                Err(_) => model.solve()?,
            },
            None => model.solve()?,
        };
        self.iterations = history.len() - 1;
        let sol = Solution::of(model.clone(), x.clone(), marches, history);
        let water = Self::water_sides(&model, &sol)?;
        self.model = model;
        self.x = Some(x);
        self.solution = Some(sol);
        Ok(water)
    }

    /// The water side of every segment of a solution.
    fn water_sides(model: &Model, sol: &Solution) -> Result<Vec<Vec<WaterSide>>, Error> {
        let mut out = Vec::with_capacity(model.tubes.len());
        for (ti, (t, mr)) in model.tubes.iter().zip(&sol.marches).enumerate() {
            let n = t.q.len();
            let dz = t.dz();
            let m = sol.x[t.im];
            let mut row = vec![
                WaterSide {
                    t_bulk: f64::NAN,
                    htc: f64::NAN,
                    t_wall: f64::NAN,
                    q: 0.0,
                    chf_ratio: f64::NAN,
                };
                n
            ];
            for g in &mr.segments {
                let i = ((g.s / dz).floor() as usize).min(n - 1);
                let at: Segment = if g.q >= Q_FLOOR && g.htc.is_finite() {
                    g.clone()
                } else {
                    model.wall_at(ti, i, g.p, g.h, m, Q_FLOOR)?
                };
                row[i] = WaterSide {
                    t_bulk: g.t,
                    htc: at.htc,
                    t_wall: g.t_wall,
                    q: g.q,
                    chf_ratio: g.chf_ratio,
                };
            }
            if let Some(w) = row.iter().find(|w| !(w.htc.is_finite() && w.htc > 0.0)) {
                return Err(Error::new(
                    codes::RANGE,
                    format!(
                        "tube `{}`: no positive heat transfer coefficient of the water side \
                         (T_b = {} K)",
                        t.name, w.t_bulk
                    ),
                ));
            }
            out.push(row);
        }
        Ok(out)
    }
}
