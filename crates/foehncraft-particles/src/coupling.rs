//! Two-way coupling of the parcels with the gas (SP2, D247 to D252,
//! `docs/requirements/step28-sprays.md` section 5).
//!
//! A *parcel pass* ([`crate::Spray::pass`]) is a pure function of the
//! parcels at the old time level, the frozen gas state and the keyed random
//! numbers (D247, D248). Each sub-step of a parcel writes one
//! [`Record`]: the cell, the parcel, the sub-step, the momentum the gas
//! receives (the parcel's own change minus its weight, D249) and the slope
//! of that exchange with respect to the cell's gas velocity (D250).
//!
//! The records are sorted by (cell, parcel, sub-step) and summed per cell in
//! that order ([`CellSources::from_records`]): no atomics, no dependence on
//! the parcel storage order or the thread count.
//!
//! The gas takes the sources in value-plus-slope form,
//! S = E − β (u − u*), with E = impulse / Δt and β = slope / Δt summed per
//! cell and u* the frozen cell velocity of the pass. The **closing step**
//! ([`crate::Spray::commit`]) gives each parcel the correction the gas
//! received: + slope_r (u_new − u*) on its momentum for every record r,
//! so both books hold the same numbers.
//!
//! **Slope of a record (own derivation, implementation note of D250).**
//! For a sub-step of length δt with response time τ the exact relaxation
//! gives the end velocity v₁ = a + (v₀ − a) e, e = exp(−δt/τ), a = u* +
//! g_eff τ, so ∂v₁/∂u* = 1 − e for that sub-step alone. A parcel takes
//! several sub-steps in one flow step; its velocity at the end of the step
//! depends on the gas of sub-step r through (1 − e_r) Π_{l>r} e_l. The
//! slope of record r is the parcel mass times that product, so that the
//! sum of the slopes of one parcel is m (1 − Π e) ≤ m: the closing step
//! moves a parcel by at most the largest change of the gas it saw, and it
//! cannot overshoot the gas. Summing the per-sub-step slopes m (1 − e_r)
//! instead would count the same relaxation several times (up to m times
//! the number of sub-steps). For a parcel that stays in one cell the slope
//! is the exact linearisation of the parcel's momentum change; for one
//! that crosses cells the later decay belongs, strictly, to the off-diagonal
//! coupling of the neighbour cells, which the diagonal form leaves to the
//! passes (D252).

/// One sub-step of one parcel: what the gas of cell `cell` receives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Record {
    /// Cell of the sub-step.
    pub cell: u32,
    /// Parcel identity.
    pub id: u64,
    /// Sub-step number of the parcel in this step (from 1).
    pub sub: u32,
    /// Momentum the gas receives in the sub-step (kg m/s): the drag impulse
    /// on the parcel with the opposite sign.
    pub impulse: [f64; 3],
    /// Slope of the impulse with respect to the gas velocity of the cell
    /// (kg), the parcel mass times (1 − e_r) Π_{l>r} e_l.
    pub slope: f64,
    /// Length of the sub-step (s).
    pub dt: f64,
    /// exp(−δt/τ) of the sub-step (kept for the slopes).
    pub(crate) decay: f64,
}

impl Record {
    /// The sort key of D249.
    pub fn key(&self) -> (u32, u64, u32) {
        (self.cell, self.id, self.sub)
    }
}

/// Sets the slopes of the records of one parcel (in sub-step order), of
/// mass `m`: slope_r = m (1 − e_r) Π_{l>r} e_l (module documentation).
pub(crate) fn set_slopes(records: &mut [Record], m: f64) {
    let mut later = 1.0;
    for r in records.iter_mut().rev() {
        r.slope = m * (1.0 - r.decay) * later;
        later *= r.decay;
    }
}

/// The per-cell sums of a pass: sorted by cell, each summed over its
/// records in the order (parcel, sub-step) (D249).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CellSources {
    /// The cells with records, ascending.
    pub cells: Vec<u32>,
    /// Momentum the gas of each cell receives over the step (kg m/s).
    pub impulse: Vec<[f64; 3]>,
    /// Summed slope of each cell (kg).
    pub slope: Vec<f64>,
    /// Records per cell.
    pub records: Vec<u32>,
}

impl CellSources {
    /// The sorted sums of `records` (any order on entry).
    pub fn from_records(records: &[Record]) -> Self {
        let mut ix: Vec<usize> = (0..records.len()).collect();
        ix.sort_unstable_by_key(|&i| records[i].key());
        let mut s = Self::default();
        for i in ix {
            let r = &records[i];
            if s.cells.last() != Some(&r.cell) {
                s.cells.push(r.cell);
                s.impulse.push(r.impulse);
                s.slope.push(r.slope);
                s.records.push(1);
            } else {
                let k = s.cells.len() - 1;
                for (a, b) in s.impulse[k].iter_mut().zip(r.impulse) {
                    *a += b;
                }
                s.slope[k] += r.slope;
                s.records[k] += 1;
            }
        }
        s
    }

    /// Position of `cell` in [`CellSources::cells`].
    pub fn find(&self, cell: u32) -> Option<usize> {
        self.cells.binary_search(&cell).ok()
    }

    /// True without records.
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }
}

/// The momentum books of one committed step (kg m/s over the step, before
/// division by Δt): what the ledger rows of D244 are made of.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct StepBooks {
    /// Momentum of the parcels in the flow at the start of the step.
    pub start: [f64; 3],
    /// At the end (after the closing step).
    pub end: [f64; 3],
    /// Received from the gas (the opposite of what the gas receives, per
    /// cell sums of the closed sources summed in cell order).
    pub from_gas: [f64; 3],
    /// Impulse of the weight (gravity less buoyancy) over the sub-steps.
    pub gravity: [f64; 3],
    /// Momentum of the parcels released in the step, at their release.
    pub injected: [f64; 3],
    /// Momentum of the parcels that left the flow (stuck, escaped, lost),
    /// with the sign of a loss (negative of what they carried).
    pub removed: [f64; 3],
    /// Impulse of the walls on rebounding and mirrored parcels.
    pub walls: [f64; 3],
}

impl StepBooks {
    /// end − start − Σ terms, per component: zero to round-off.
    pub fn imbalance(&self) -> [f64; 3] {
        let mut r = [0.0; 3];
        for (k, rk) in r.iter_mut().enumerate() {
            *rk = self.end[k]
                - self.start[k]
                - (self.from_gas[k]
                    + self.gravity[k]
                    + self.injected[k]
                    + self.removed[k]
                    + self.walls[k]);
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(cell: u32, id: u64, sub: u32, decay: f64) -> Record {
        Record {
            cell,
            id,
            sub,
            impulse: [1.0, 2.0, 3.0],
            slope: 0.0,
            dt: 1.0,
            decay,
        }
    }

    #[test]
    fn slopes_of_one_parcel_add_up_to_its_relaxation() {
        let mut r = vec![rec(0, 1, 1, 0.5), rec(0, 1, 2, 0.5), rec(1, 1, 3, 0.25)];
        set_slopes(&mut r, 2.0);
        let total: f64 = r.iter().map(|x| x.slope).sum();
        assert!((total - 2.0 * (1.0 - 0.5 * 0.5 * 0.25)).abs() < 1e-15);
        // the first sub-step is decayed by the two after it
        assert!((r[0].slope - 2.0 * 0.5 * 0.5 * 0.25).abs() < 1e-15);
        assert!((r[2].slope - 2.0 * 0.75).abs() < 1e-15);
    }

    #[test]
    fn cell_sums_do_not_depend_on_the_order_of_the_records() {
        let mut r: Vec<Record> = (0..50)
            .map(|i| {
                let mut x = rec((i * 7 % 5) as u32, i as u64, 1 + (i % 3) as u32, 0.3);
                x.impulse = [0.1 * i as f64, 1.0 / (1.0 + i as f64), -(i as f64).sqrt()];
                x.slope = 1.0 / (3.0 + i as f64);
                x
            })
            .collect();
        let a = CellSources::from_records(&r);
        r.reverse();
        r.swap(3, 17);
        let b = CellSources::from_records(&r);
        assert_eq!(a, b);
        assert_eq!(a.cells, vec![0, 1, 2, 3, 4]);
        assert_eq!(a.records.iter().sum::<u32>(), 50);
    }
}
