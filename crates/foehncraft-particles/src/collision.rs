//! Collisions of the drops (D331): the no-time-counter (NTC) sampling of
//! Schmidt and Rutland (2000, *A new droplet collision algorithm*, J.
//! Comput. Phys. 164, 62-80) with the collision outcomes of O'Rourke (1981,
//! PhD thesis, Princeton; Amsden, O'Rourke and Butler 1989, LA-11560-MS
//! §IV.D).
//!
//! A parcel stands for `count` identical drops. In a cell of volume V with
//! N parcels the expected number of drop collisions between parcels i and j
//! in a step Δt is q_i q_j σ_ij v_ij Δt / V, with σ_ij = π (r_i + r_j)² and
//! v_ij the relative speed. NTC draws
//!
//!   M = N (N − 1)/2 · q_max (σ v)_max Δt / V
//!
//! candidate pairs (rounded at random) and accepts a pair with the
//! probability P = q_m σ_ij v_ij / (q_max (σ v)_max), q_m = max(q_i, q_j);
//! an accepted pair is q_min drop collisions (each drop of the parcel with
//! fewer drops meets one drop of the other), so the expected number of drop
//! collisions is that of the pair, at a cost linear in the collisions
//! instead of quadratic in the parcels.
//!
//! Outcome (O'Rourke): with r_s ≤ r_l the radii, the collision Weber number
//! We = ρ_l v² r_s/σ_l and γ = r_l/r_s, the critical impact parameter is
//! b_cr² = (r_s + r_l)² min(1, 2.4 f(γ)/We), f(γ) = γ³ − 2.4 γ² + 2.7 γ. The
//! impact parameter b = (r_s + r_l) √ξ (ξ uniform): below b_cr the drops
//! coalesce, above it they separate with their sizes and the velocities
//!
//!   v₁' = (m₁ v₁ + m₂ v₂ + m₂ (v₁ − v₂) s)/(m₁ + m₂),  s = (b − b_cr)/(r₁ + r₂ − b_cr)
//!
//! (and the same with 1 and 2 swapped). Mass, momentum and liquid enthalpy
//! (constant c_p) are kept exactly; coalescence lowers the kinetic energy.
//!
//! The random numbers are keyed by the smallest parcel identity of the
//! cell, the time at the end of the step and the candidate (stream
//! `SPRAY_COLLISION`, D121); the parcels of a cell are taken in the order of
//! their identities, so a run gives the same drops for any thread count or
//! parcel order.

use crate::spray::Parcels;
use foehncraft_core::rng::{self, stream};
use rayon::prelude::*;

/// The most candidate pairs of one cell in one step per parcel of the cell;
/// a cell that would need more takes this many and is counted
/// ([`CollisionTally::capped_cells`]).
pub const MAX_CANDIDATES_PER_PARCEL: f64 = 20000.0;

/// What the collisions of one step did.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct CollisionTally {
    /// Accepted parcel pairs (each is q_min drop collisions).
    pub events: u64,
    /// Of them, coalescences.
    pub coalescences: u64,
    /// Drop collisions (the sum of q_min over the events).
    pub drop_collisions: f64,
    /// Parcels whose drops were all taken up by coalescence (removed).
    pub merged_parcels: u64,
    /// Cells whose candidate count was capped.
    pub capped_cells: u64,
}

/// The liquid and step data the collisions need.
#[derive(Debug, Clone, Copy)]
pub struct CollisionInput {
    /// Liquid density (kg/m³).
    pub rho_l: f64,
    /// Surface tension (N/m).
    pub sigma: f64,
    /// Time step (s).
    pub dt: f64,
    /// Time at the end of the step (s), a word of the random key.
    pub t1: f64,
    /// Coalescence on (otherwise every collision is a grazing separation).
    pub coalescence: bool,
}

#[derive(Debug, Clone, Copy)]
struct Drop {
    i: usize,
    r: f64,
    q: f64,
    v: [f64; 3],
    temp: f64,
    merged: bool,
}

/// O'Rourke's critical impact parameter squared over (r₁ + r₂)².
pub fn critical_impact(we: f64, gamma: f64) -> f64 {
    let f = gamma * gamma * gamma - 2.4 * gamma * gamma + 2.7 * gamma;
    if we <= 0.0 {
        1.0
    } else {
        (2.4 * f / we).min(1.0)
    }
}

/// Collides the drops of the parcels `p` in place (the cells of `p.cell`
/// with the volumes `volume`); parcels flagged in `skip` (the injected
/// blobs) take no part. Parcels left with no drops are removed.
pub fn collide(
    p: &mut Parcels,
    volume: &[f64],
    skip: &dyn Fn(usize) -> bool,
    inp: CollisionInput,
) -> CollisionTally {
    // the members of each cell in the order of their identities
    let mut members: Vec<(u32, u64, usize)> = (0..p.len())
        .filter(|&i| !skip(i) && p.count[i] > 0.0 && p.d[i] > 0.0)
        .map(|i| (p.cell[i], p.id[i], i))
        .collect();
    members.sort_unstable_by_key(|&(c, id, _)| (c, id));
    let mut groups: Vec<&[(u32, u64, usize)]> = Vec::new();
    let mut a = 0;
    while a < members.len() {
        let mut b = a + 1;
        while b < members.len() && members[b].0 == members[a].0 {
            b += 1;
        }
        if b - a >= 2 {
            groups.push(&members[a..b]);
        }
        a = b;
    }
    let pr: &Parcels = p;
    let done: Vec<(Vec<Drop>, CollisionTally)> = groups
        .par_iter()
        .map(|g| {
            let mut drops: Vec<Drop> = g
                .iter()
                .map(|&(_, _, i)| Drop {
                    i,
                    r: 0.5 * pr.d[i],
                    q: pr.count[i],
                    v: pr.v[i],
                    temp: pr.temp.get(i).copied().unwrap_or(0.0),
                    merged: false,
                })
                .collect();
            let t = cell(&mut drops, volume[g[0].0 as usize], g[0].1, inp);
            (drops, t)
        })
        .collect();
    let mut tally = CollisionTally::default();
    let mut gone = vec![false; p.len()];
    for (drops, t) in done {
        tally.events += t.events;
        tally.coalescences += t.coalescences;
        tally.drop_collisions += t.drop_collisions;
        tally.merged_parcels += t.merged_parcels;
        tally.capped_cells += t.capped_cells;
        for d in drops {
            let i = d.i;
            if d.q <= 0.0 {
                gone[i] = true;
                continue;
            }
            if d.merged && i < p.brk.len() {
                // a new drop: undistorted, no wave grown on it yet
                p.brk[i].y = 0.0;
                p.brk[i].ydot = 0.0;
                p.brk[i].rt = 0.0;
                p.brk[i].m0 = d.q * inp.rho_l * std::f64::consts::PI / 6.0 * 8.0 * d.r.powi(3);
            }
            p.d[i] = 2.0 * d.r;
            p.count[i] = d.q;
            p.v[i] = d.v;
            if i < p.temp.len() {
                p.temp[i] = d.temp;
            }
        }
    }
    if gone.iter().any(|&g| g) {
        p.retain(|i| !gone[i]);
    }
    tally
}

/// The collisions of the drops of one cell (in the order of the parcel
/// identities; `key` the smallest identity).
fn cell(drops: &mut [Drop], vol: f64, key: u64, inp: CollisionInput) -> CollisionTally {
    let mut t = CollisionTally::default();
    let n = drops.len();
    if n < 2 || vol <= 0.0 || inp.dt <= 0.0 {
        return t;
    }
    let q_max = drops.iter().map(|d| d.q).fold(0.0, f64::max);
    let r_max = drops.iter().map(|d| d.r).fold(0.0, f64::max);
    // an upper bound of the relative speed: twice the largest distance of
    // a velocity from the mean (the velocities stay inside their hull)
    let mut mean = [0.0; 3];
    for d in drops.iter() {
        for (m, v) in mean.iter_mut().zip(d.v) {
            *m += v / n as f64;
        }
    }
    let v_max = 2.0 * drops.iter().map(|d| dist(d.v, mean)).fold(0.0, f64::max);
    let sv_max = std::f64::consts::PI * 4.0 * r_max * r_max * v_max;
    if !(sv_max > 0.0 && q_max > 0.0) {
        return t;
    }
    let pairs = 0.5 * n as f64 * (n as f64 - 1.0);
    let mut m_exp = pairs * q_max * sv_max * inp.dt / vol;
    let cap = MAX_CANDIDATES_PER_PARCEL * n as f64;
    if m_exp > cap {
        m_exp = cap;
        t.capped_cells = 1;
    }
    let tw = inp.t1.to_bits();
    let u = |e: u64| rng::uniform_keyed([stream::SPRAY_COLLISION, key, tw, e]);
    let m = (m_exp + u(0)).floor() as u64;
    for c in 0..m {
        let e = 1 + 4 * c;
        let i = ((u(e) * n as f64) as usize).min(n - 1);
        let mut j = ((u(e + 1) * (n - 1) as f64) as usize).min(n - 2);
        if j >= i {
            j += 1;
        }
        let (a, b) = (drops[i], drops[j]);
        if a.q <= 0.0 || b.q <= 0.0 {
            continue;
        }
        let w = [a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]];
        let vr = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        let rs = a.r + b.r;
        let sv = std::f64::consts::PI * rs * rs * vr;
        let prob = (a.q.max(b.q) * sv / (q_max * sv_max)).min(1.0);
        if u(e + 2) >= prob {
            continue;
        }
        // the collector: the parcel with fewer drops, each of its drops
        // meets one drop of the other
        let (ia, ib) = if a.q <= b.q { (i, j) } else { (j, i) };
        let (ca, cb) = (drops[ia], drops[ib]);
        let qa = ca.q;
        let (r_s, r_l) = if ca.r <= cb.r {
            (ca.r, cb.r)
        } else {
            (cb.r, ca.r)
        };
        let we = inp.rho_l * vr * vr * r_s / inp.sigma;
        let bcr2 = critical_impact(we, r_l / r_s);
        let xi = u(e + 3);
        let ma = ca.r * ca.r * ca.r;
        let mb = cb.r * cb.r * cb.r;
        t.events += 1;
        t.drop_collisions += qa;
        if inp.coalescence && xi < bcr2 {
            t.coalescences += 1;
            let mt = ma + mb;
            let v: [f64; 3] = std::array::from_fn(|k| (ma * ca.v[k] + mb * cb.v[k]) / mt);
            drops[ia].r = mt.cbrt();
            drops[ia].v = v;
            drops[ia].temp = (ma * ca.temp + mb * cb.temp) / mt;
            drops[ia].merged = true;
            let left = cb.q - qa;
            if left <= 1e-12 * cb.q {
                drops[ib].q = 0.0;
                t.merged_parcels += 1;
            } else {
                drops[ib].q = left;
            }
        } else {
            // grazing separation: the sizes stay, the velocities move
            // toward the common one by the share 1 − s
            let b = xi.sqrt();
            let bc = bcr2.sqrt();
            let s = if inp.coalescence && bc < 1.0 {
                ((b - bc) / (1.0 - bc)).clamp(0.0, 1.0)
            } else {
                b
            };
            let mt = ma + mb;
            let p: [f64; 3] = std::array::from_fn(|k| ma * ca.v[k] + mb * cb.v[k]);
            drops[ia].v = std::array::from_fn(|k| (p[k] + mb * (ca.v[k] - cb.v[k]) * s) / mt);
            // qa of the drops of b changed
            let f = qa / cb.q;
            drops[ib].v = std::array::from_fn(|k| {
                let vb = (p[k] + ma * (cb.v[k] - ca.v[k]) * s) / mt;
                (1.0 - f) * cb.v[k] + f * vb
            });
        }
    }
    t
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spray::{BreakupState, Eddy};

    fn parcels(list: &[(u64, f64, f64, [f64; 3])]) -> Parcels {
        let mut p = Parcels::default();
        for &(id, d, q, v) in list {
            p.id.push(id);
            p.x.push([0.0; 3]);
            p.v.push(v);
            p.d.push(d);
            p.count.push(q);
            p.born.push(0.0);
            p.temp.push(300.0 + id as f64);
            p.eddy.push(Eddy::default());
            p.brk.push(BreakupState::default());
            p.cell.push(0);
        }
        p
    }

    fn totals(p: &Parcels) -> (f64, [f64; 3], f64) {
        let mut m = 0.0;
        let mut mv = [0.0; 3];
        let mut h = 0.0;
        for i in 0..p.len() {
            let mi = p.count[i] * p.d[i].powi(3);
            m += mi;
            for (m, v) in mv.iter_mut().zip(p.v[i]) {
                *m += mi * v;
            }
            h += mi * p.temp[i];
        }
        (m, mv, h)
    }

    fn input(t1: f64, coalescence: bool) -> CollisionInput {
        CollisionInput {
            rho_l: 700.0,
            sigma: 0.02,
            dt: 1e-5,
            t1,
            coalescence,
        }
    }

    /// Two opposed groups of parcels: the drop collisions over many steps
    /// match q_i q_j σ v Δt / V summed over the pairs (grazing only; the
    /// parcels are set up anew each step).
    #[test]
    fn ntc_rate_matches_the_pair_rate() {
        // (drops of the parcels moving +x, of those moving −x, volume,
        // steps): about 0.25 accepted pairs per step, so the velocities a
        // collision changes hardly meet again in the same step
        for (qa, qb, vol, steps) in [(1.0, 1.0, 1e-10, 40_000), (10.0, 1000.0, 1e-6, 400_000)] {
            let (r, v, n) = (10e-6, 10.0, 20usize);
            let mut sum = 0.0;
            for s in 0..steps {
                let list: Vec<_> = (0..n)
                    .map(|i| {
                        let (sign, q) = if i % 2 == 0 { (1.0, qa) } else { (-1.0, qb) };
                        (i as u64, 2.0 * r, q, [sign * v, 0.0, 0.0])
                    })
                    .collect();
                let mut p = parcels(&list);
                let t = collide(&mut p, &[vol], &|_| false, input(s as f64 * 1e-5, false));
                sum += t.drop_collisions;
            }
            let half = (n / 2) as f64;
            let expected =
                half * half * qa * qb * std::f64::consts::PI * (2.0 * r).powi(2) * (2.0 * v) * 1e-5
                    / vol
                    * steps as f64;
            let ratio = sum / expected;
            assert!(
                (ratio - 1.0).abs() < 0.03,
                "q {qa}/{qb}: {sum} vs {expected}: {ratio}"
            );
        }
    }

    /// Mass, momentum and liquid enthalpy stay through coalescence and
    /// grazing; slow drops coalesce (We small: b_cr = r₁ + r₂).
    #[test]
    fn collisions_keep_mass_momentum_and_enthalpy() {
        for coalescence in [true, false] {
            let list: Vec<_> = (0..30u64)
                .map(|i| {
                    let a = i as f64;
                    (
                        i,
                        (8.0 + (a * 0.37).sin() * 3.0) * 1e-6,
                        50.0 + 40.0 * (a * 1.3).cos().abs(),
                        [(a * 0.7).sin(), (a * 1.1).cos(), 0.3 * (a * 0.2).sin()],
                    )
                })
                .collect();
            let mut p = parcels(&list);
            let before = totals(&p);
            let mut tally = CollisionTally::default();
            for s in 0..50 {
                let t = collide(&mut p, &[1e-12], &|_| false, input(s as f64, coalescence));
                tally.events += t.events;
                tally.coalescences += t.coalescences;
            }
            let after = totals(&p);
            assert!(tally.events > 10, "{tally:?}");
            assert_eq!(tally.coalescences > 0, coalescence, "{tally:?}");
            let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
            assert!(rel(after.0, before.0) < 1e-12, "mass {before:?} {after:?}");
            assert!(rel(after.2, before.2) < 1e-12, "enthalpy");
            let scale = before.1.iter().map(|x| x.abs()).fold(0.0, f64::max) + before.0;
            for k in 0..3 {
                assert!(
                    (after.1[k] - before.1[k]).abs() < 1e-12 * scale,
                    "momentum {k}"
                );
            }
        }
    }

    /// The same drops whatever the order of the parcels.
    #[test]
    fn collisions_do_not_depend_on_the_parcel_order() {
        let list: Vec<_> = (0..40u64)
            .map(|i| {
                let a = i as f64;
                (
                    i,
                    (6.0 + (a * 0.9).sin() * 2.0) * 1e-6,
                    100.0,
                    [(a * 0.3).sin() * 5.0, 0.0, 0.0],
                )
            })
            .collect();
        let mut rev = list.clone();
        rev.reverse();
        let mut a = parcels(&list);
        let mut b = parcels(&rev);
        for s in 0..20 {
            collide(&mut a, &[1e-12], &|_| false, input(s as f64, true));
            collide(&mut b, &[1e-12], &|_| false, input(s as f64, true));
        }
        let key = |p: &Parcels| {
            let mut v: Vec<_> = (0..p.len())
                .map(|i| {
                    (
                        p.id[i],
                        p.d[i].to_bits(),
                        p.count[i].to_bits(),
                        p.v[i][0].to_bits(),
                    )
                })
                .collect();
            v.sort();
            v
        };
        assert_eq!(key(&a), key(&b));
    }

    /// O'Rourke's critical impact parameter: everything coalesces at small
    /// Weber numbers, and the limit 2.4 f(γ)/We for equal drops (f(1) = 1.3).
    #[test]
    fn critical_impact_parameter() {
        assert_eq!(critical_impact(0.1, 1.0), 1.0);
        assert!((critical_impact(31.2, 1.0) - 0.1).abs() < 1e-12);
    }
}
