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
    /// Heat the gas receives in the sub-step (J): the parcel's own change
    /// of sensible enthalpy with the opposite sign (SP2b).
    pub heat: f64,
    /// Slope of the heat with respect to the gas temperature of the cell
    /// (J/K), the parcel's heat capacity times (1 − e_r) Π_{l>r} e_l of the
    /// thermal relaxation.
    pub heat_slope: f64,
    /// exp(−δt/τ_h) of the sub-step (1 without heat exchange).
    pub(crate) heat_decay: f64,
    /// Vapour mass the gas receives in the sub-step (kg; SP5): the parcel's
    /// own loss of liquid (negative: condensation). `heat` then holds the
    /// whole enthalpy the parcel gives up (sensible and that of the vapour,
    /// with the liquid enthalpy c_l T as the reference).
    pub vapour: f64,
    /// Slope of the vapour with respect to the vapour mass fraction of the
    /// gas of the cell (kg per unit mass fraction): the chord to the
    /// saturation state at the droplet surface, vapour / (Y_s − Y*) (D250),
    /// scaled so that the closing step cannot take more than the parcel
    /// holds (`cap_vapour_slopes`).
    pub vapour_slope: f64,
    /// The frozen vapour mass fraction Y* the sub-step saw.
    pub vapour_y: f64,
    /// Enthalpy of the liquid at the end of the sub-step (J/kg, c_l T, the
    /// reference of the parcels): what a closing correction of the vapour
    /// and vapour held back by the limiter of D253 carry. The mass the
    /// parcel keeps or regains stays liquid at its temperature, so the gas
    /// keeps the latent heat of it: the heat the droplet draws follows its
    /// evaporation (Q ≈ L β ΔY at the wet bulb), and an evaporating cell
    /// tends to the adiabatic saturation, not below it (SP5).
    pub liquid_h: f64,
}

/// The vapour slopes of one parcel (SP5, D250): the closing step gives the
/// parcel Σ β_r (Y_new − Y*_r) back, which takes at most Σ β_r Y*_r from it
/// (the gas cannot go below Y = 0). A parcel that evaporated completely
/// (`m_end` = 0) has no slopes: it cannot take vapour back or give more.
/// Otherwise the slopes are scaled down so that Σ β_r Y*_r ≤ m_end.
pub(crate) fn cap_vapour_slopes(records: &mut [Record], m_end: f64) {
    let held: f64 = records.iter().map(|r| r.vapour_slope * r.vapour_y).sum();
    let f = if m_end <= 0.0 {
        0.0
    } else if held > m_end {
        m_end / held
    } else {
        return;
    };
    for r in records.iter_mut() {
        r.vapour_slope *= f;
    }
}

impl Record {
    /// The sort key of D249.
    pub fn key(&self) -> (u32, u64, u32) {
        (self.cell, self.id, self.sub)
    }
}

/// Sets the slopes of the records of one parcel (in sub-step order), of
/// mass `m` and heat capacity `mc` (J/K): slope_r = m (1 − e_r) Π_{l>r}
/// e_l (module documentation), and the same for the heat with the thermal
/// decays.
pub(crate) fn set_slopes(records: &mut [Record], m: f64, mc: f64) {
    let mut later = 1.0;
    let mut later_h = 1.0;
    for r in records.iter_mut().rev() {
        r.slope = m * (1.0 - r.decay) * later;
        later *= r.decay;
        r.heat_slope = mc * (1.0 - r.heat_decay) * later_h;
        later_h *= r.heat_decay;
    }
}

/// The slopes of the records of one steady stream (D260), in sub-step
/// order: the local derivative of the exchange of each visit of a cell
/// with respect to that cell's gas, the inflow to the visit held fixed.
/// slope_r = m (1 − e_r) Π e_l over the later sub-steps l of the same
/// visit, so the slopes of a visit add up to m (1 − Π e) ≤ m. A stream
/// has no closing step (no state carried from pass to pass), so the decay
/// of later cells does not belong to the slope of this one; the transient
/// form of [`set_slopes`] would make the source of the first cells of a
/// long path nearly explicit.
pub(crate) fn set_slopes_local(records: &mut [Record], m: f64, mc: f64) {
    let mut later = 1.0;
    let mut later_h = 1.0;
    let mut next: Option<u32> = None;
    for r in records.iter_mut().rev() {
        if next != Some(r.cell) {
            later = 1.0;
            later_h = 1.0;
        }
        next = Some(r.cell);
        r.slope = m * (1.0 - r.decay) * later;
        later *= r.decay;
        r.heat_slope = mc * (1.0 - r.heat_decay) * later_h;
        later_h *= r.heat_decay;
    }
}

/// The factor of the saturation limiter of D253 for a cell whose
/// evaporating records sum to `evaporated` and whose net vapour is `net`
/// (kg, or kg/s): the factor f in (0, 1] of the evaporating records that
/// brings the net vapour down to `capacity`; 1 when it is not above.
pub fn limiter_factor(evaporated: f64, net: f64, capacity: f64) -> f64 {
    // NaN in net or capacity also leaves the records as they are
    if net.partial_cmp(&capacity) != Some(std::cmp::Ordering::Greater) || evaporated <= 0.0 {
        return 1.0;
    }
    let condensed = net - evaporated;
    ((capacity.max(0.0) - condensed) / evaporated).clamp(0.0, 1.0)
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
    /// Heat the gas of each cell receives over the step (J).
    pub heat: Vec<f64>,
    /// Summed heat slope of each cell (J/K).
    pub heat_slope: Vec<f64>,
    /// Vapour the gas of each cell receives over the step (kg; SP5).
    pub vapour: Vec<f64>,
    /// Summed vapour slope of each cell (kg per unit mass fraction).
    pub vapour_slope: Vec<f64>,
    /// Σ vapour slope × liquid enthalpy of each cell (J per unit mass
    /// fraction): the enthalpy that follows the vapour slope.
    pub vapour_slope_h: Vec<f64>,
    /// Σ of the evaporating (positive) vapour records of each cell (kg).
    pub evaporated: Vec<f64>,
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
                s.heat.push(r.heat);
                s.heat_slope.push(r.heat_slope);
                s.vapour.push(r.vapour);
                s.vapour_slope.push(r.vapour_slope);
                s.vapour_slope_h.push(r.vapour_slope * r.liquid_h);
                s.evaporated.push(r.vapour.max(0.0));
                s.records.push(1);
            } else {
                let k = s.cells.len() - 1;
                for (a, b) in s.impulse[k].iter_mut().zip(r.impulse) {
                    *a += b;
                }
                s.slope[k] += r.slope;
                s.heat[k] += r.heat;
                s.heat_slope[k] += r.heat_slope;
                s.vapour[k] += r.vapour;
                s.vapour_slope[k] += r.vapour_slope;
                s.vapour_slope_h[k] += r.vapour_slope * r.liquid_h;
                s.evaporated[k] += r.vapour.max(0.0);
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
    /// The enthalpy of the parcels (J, Σ m c_l T, the liquid at 0 K as
    /// the reference): [start, end, from the gas, injected, removed]
    /// (SP2b; with evaporation "from the gas" holds the enthalpy of the
    /// vapour as well, SP5).
    pub enthalpy: [f64; 5],
    /// The liquid mass of the parcels (kg): [start, end, from the gas
    /// (minus the vapour the gas received), injected, removed] (SP5).
    pub mass: [f64; 5],
}

/// The change of the gas of a cell against the frozen state of a pass, for
/// the closing step (D250).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Delta {
    /// u_new − u* (m/s).
    pub u: [f64; 3],
    /// T_new − T* (K).
    pub t: f64,
    /// Y_new − Y* of the vapour (SP5).
    pub y: f64,
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

    /// end − start − (from the gas + injected + removed) of the enthalpy:
    /// zero to round-off.
    pub fn heat_imbalance(&self) -> f64 {
        let [start, end, from_gas, injected, removed] = self.enthalpy;
        end - start - (from_gas + injected + removed)
    }

    /// The same for the liquid mass (SP5): zero to round-off.
    pub fn mass_imbalance(&self) -> f64 {
        let [start, end, from_gas, injected, removed] = self.mass;
        end - start - (from_gas + injected + removed)
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
            heat: 0.5,
            heat_slope: 0.0,
            heat_decay: decay,
            vapour: 0.0,
            vapour_slope: 0.0,
            vapour_y: 0.0,
            liquid_h: 0.0,
        }
    }

    #[test]
    fn slopes_of_one_parcel_add_up_to_its_relaxation() {
        let mut r = vec![rec(0, 1, 1, 0.5), rec(0, 1, 2, 0.5), rec(1, 1, 3, 0.25)];
        set_slopes(&mut r, 2.0, 3.0);
        let total: f64 = r.iter().map(|x| x.slope).sum();
        assert!((total - 2.0 * (1.0 - 0.5 * 0.5 * 0.25)).abs() < 1e-15);
        // the first sub-step is decayed by the two after it
        assert!((r[0].slope - 2.0 * 0.5 * 0.5 * 0.25).abs() < 1e-15);
        assert!((r[2].slope - 2.0 * 0.75).abs() < 1e-15);
        let heat: f64 = r.iter().map(|x| x.heat_slope).sum();
        assert!((heat - 3.0 * (1.0 - 0.5 * 0.5 * 0.25)).abs() < 1e-15);
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
