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
//! Outcomes of Post and Abraham (2002, *Modeling the outcome of drop-drop
//! collisions in Diesel sprays*, Int. J. Multiphase Flow 28, 997-1019;
//! D349), with `post`: with Δ = r_s/r_l ≤ 1, B = b/(r_s + r_l) and the
//! Weber number of the small diameter We_d = 2 We,
//!
//! 1. bouncing (Estrade et al. 1999, Int. J. Heat Fluid Flow 20) when
//!    We_d < Δ (1 + Δ²)(4φ′ − 12)/(χ (1 − B²)), φ′ = 3.351, τ = (1 − B)(1 + Δ),
//!    χ = 1 − (2 − τ)²(1 + τ)/4 for τ > 1, χ = τ²(3 − τ)/4 otherwise;
//! 2. else, below O'Rourke's b_cr, reflexive separation (Ashgriz and Poo
//!    1990, J. Fluid Mech. 221) when
//!    We_d > 3 [7 (1 + Δ³)^{2/3} − 4 (1 + Δ²)] Δ (1 + Δ³)²/(Δ⁶ η₁ + η₂),
//!    η₁ = 2 (1 − ξ)² √(1 − ξ²) − 1, η₂ = 2 (Δ − ξ)² √(Δ² − ξ²) − Δ³,
//!    ξ = B (1 + Δ)/2 (the root taken as 0 where ξ > Δ), and coalescence
//!    otherwise;
//! 3. else stretching separation as O'Rourke's grazing.
//!
//! The criteria and the order are restated in Sommerfeld and Kuschel
//! (2016, Exp. Fluids 57:187, eqs. 16-20) and Magnotti (2017, thesis,
//! OSTI 1785739, eqs. 2.23-2.26). The velocities after a bounce and a
//! reflexive separation are ours (the paper was not open to us): the drops
//! keep their sizes and take the grazing velocities with the share s = B,
//! as the grazing of `coalescence = false` (a head-on one leaves both with
//! the common velocity, mass and momentum kept, energy lost).
//! No satellite drops form.
//!
//! Collision cells of their own (D350, `cell = h`): instead of the cells
//! of the flow mesh the drops meet in cubes of edge h, with V = h³. Every
//! step the cube lattice is turned by a random rotation (uniform over the
//! rotations, from a random unit quaternion, Shoemake 1992) and shifted by
//! a random vector in [0, h)³, so no lattice line stays in place and the
//! rate does not follow the flow mesh (Hou and Schmidt 2006, Int. J.
//! Multiphase Flow 32; the uniform rotated collision mesh as described by
//! Lanfaloni, MSc thesis, Politecnico di Torino 2020, p. 15). A cube cut
//! by a wall or a mirror plane is counted with its full volume (ours).
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
    /// Of the events, bounces (Post outcomes, D349).
    pub bounces: u64,
    /// Of the events, reflexive separations (Post outcomes, D349).
    pub reflexive: u64,
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
    /// The outcomes of Post and Abraham (D349) in place of O'Rourke's.
    pub post: bool,
    /// The edge of the collision cells (D350); `None`: the flow cells.
    pub cell: Option<f64>,
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

/// The Weber number of the small diameter below which two drops bounce
/// (Estrade et al. 1999), size ratio `delta` = r_s/r_l ≤ 1, impact
/// parameter `b` = b/(r_s + r_l) in [0, 1].
pub fn bounce_weber(delta: f64, b: f64) -> f64 {
    const PHI: f64 = 3.351;
    let tau = (1.0 - b) * (1.0 + delta);
    let chi = if tau > 1.0 {
        1.0 - 0.25 * (2.0 - tau).powi(2) * (1.0 + tau)
    } else {
        0.25 * tau * tau * (3.0 - tau)
    };
    let den = chi * (1.0 - b * b);
    if den <= 0.0 {
        return 0.0;
    }
    delta * (1.0 + delta * delta) * (4.0 * PHI - 12.0) / den
}

/// The Weber number of the small diameter above which two drops separate
/// reflexively (Ashgriz and Poo 1990); arguments as [`bounce_weber`].
/// Infinite where the formula has no positive value (no reflexive
/// separation).
pub fn reflexive_weber(delta: f64, b: f64) -> f64 {
    let xi = 0.5 * b * (1.0 + delta);
    let eta1 = 2.0 * (1.0 - xi).powi(2) * (1.0 - xi * xi).max(0.0).sqrt() - 1.0;
    let eta2 =
        2.0 * (delta - xi).powi(2) * (delta * delta - xi * xi).max(0.0).sqrt() - delta.powi(3);
    let d3 = 1.0 + delta.powi(3);
    let num = 3.0 * (7.0 * d3.powf(2.0 / 3.0) - 4.0 * (1.0 + delta * delta)) * delta * d3 * d3;
    let den = delta.powi(6) * eta1 + eta2;
    if den > 0.0 && num > 0.0 {
        num / den
    } else {
        f64::INFINITY
    }
}

/// The rotation (rows) and the shift (in cells, [0, 1)³) of the collision
/// lattice of the step ending at `t1` (D350): a uniform random rotation
/// from a unit quaternion of three uniform numbers (Shoemake 1992).
pub fn lattice_frame(t1: f64) -> ([[f64; 3]; 3], [f64; 3]) {
    let tw = t1.to_bits();
    let u = |e: u64| rng::uniform_keyed([stream::SPRAY_COLLISION, u64::MAX, tw, e]);
    let (u1, u2, u3) = (u(0), u(1), u(2));
    let tau = 2.0 * std::f64::consts::PI;
    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
    let (w, x, y, z) = (
        a * (tau * u2).sin(),
        a * (tau * u2).cos(),
        b * (tau * u3).sin(),
        b * (tau * u3).cos(),
    );
    let rot = [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ];
    (rot, [u(3), u(4), u(5)])
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
    // the members of each cell in the order of their identities: the flow
    // cell, or the cube of the turned and shifted collision lattice
    let lattice = inp.cell.map(|h| (h, lattice_frame(inp.t1)));
    let key_of = |i: usize| -> [i64; 3] {
        match lattice {
            None => [p.cell[i] as i64, 0, 0],
            Some((h, (rot, shift))) => {
                let x = p.x[i];
                std::array::from_fn(|a| {
                    let y = rot[a][0] * x[0] + rot[a][1] * x[1] + rot[a][2] * x[2];
                    (y / h + shift[a]).floor() as i64
                })
            }
        }
    };
    let mut members: Vec<([i64; 3], u64, usize)> = (0..p.len())
        .filter(|&i| !skip(i) && p.count[i] > 0.0 && p.d[i] > 0.0)
        .map(|i| (key_of(i), p.id[i], i))
        .collect();
    members.sort_unstable_by_key(|&(c, id, _)| (c, id));
    let mut groups: Vec<&[([i64; 3], u64, usize)]> = Vec::new();
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
            let vol = match lattice {
                None => volume[g[0].0[0] as usize],
                Some((h, _)) => h * h * h,
            };
            let t = cell(&mut drops, vol, g[0].1, inp);
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
        tally.bounces += t.bounces;
        tally.reflexive += t.reflexive;
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
        // Post and Abraham (D349): bouncing first, reflexive separation of
        // the collisions that O'Rourke would merge
        let mut bounce_or_reflex = false;
        if inp.post {
            let (delta, bb) = (r_s / r_l, xi.sqrt());
            let we_d = 2.0 * we;
            if we_d < bounce_weber(delta, bb) {
                t.bounces += 1;
                bounce_or_reflex = true;
            } else if xi < bcr2 && we_d > reflexive_weber(delta, bb) {
                t.reflexive += 1;
                bounce_or_reflex = true;
            }
        }
        if inp.coalescence && xi < bcr2 && !bounce_or_reflex {
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
            let s = if bounce_or_reflex {
                b
            } else if inp.coalescence && bc < 1.0 {
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
            post: false,
            cell: None,
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

    /// The collision lattice (D350): a proper rotation (orthonormal rows,
    /// determinant 1), a shift in [0, 1)³, a new frame every step.
    #[test]
    fn lattice_frame_is_a_rotation() {
        for s in 0..50 {
            let (r, sh) = lattice_frame(s as f64 * 1e-6);
            for a in 0..3 {
                for b in 0..3 {
                    let d: f64 = (0..3).map(|k| r[a][k] * r[b][k]).sum();
                    assert!((d - if a == b { 1.0 } else { 0.0 }).abs() < 1e-12);
                }
                assert!((0.0..1.0).contains(&sh[a]));
            }
            let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
                - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
                + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
            assert!((det - 1.0).abs() < 1e-12);
        }
        assert_ne!(lattice_frame(1e-6).0, lattice_frame(2e-6).0);
    }

    /// On the collision lattice (D350) the drop collisions of parcels
    /// spread uniformly over a box of edge L = 40 h match the pair rate of
    /// the box, N² q² σ v Δt/(2 L³), as on the flow cells (within the few
    /// per cent that the cubes cut by the box sides lose, about 3 h/(2 L));
    /// parcels in one point meet in one cube at the rate of V = h³.
    #[test]
    fn lattice_rate_matches_the_pair_rate() {
        let (r, v, q) = (10e-6, 10.0, 1.0);
        let run = |n: usize, h: f64, spread: f64, steps: usize| {
            let mut sum = 0.0;
            for s in 0..steps {
                let mut p = Parcels::default();
                for i in 0..n {
                    let w = |e: u64| rng::uniform_keyed([99, i as u64, s as u64, e]);
                    let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                    p.id.push(i as u64);
                    p.x.push([spread * w(0), spread * w(1), spread * w(2)]);
                    p.v.push([sign * v, 0.0, 0.0]);
                    p.d.push(2.0 * r);
                    p.count.push(q);
                    p.born.push(0.0);
                    p.temp.push(300.0);
                    p.eddy.push(Eddy::default());
                    p.brk.push(BreakupState::default());
                    p.cell.push(0);
                }
                let mut inp = input(s as f64 * 1e-5, false);
                inp.cell = Some(h);
                sum += collide(&mut p, &[1.0], &|_| false, inp).drop_collisions;
            }
            let half = (n / 2) as f64;
            let rate =
                half * half * q * q * std::f64::consts::PI * (2.0 * r).powi(2) * (2.0 * v) * 1e-5;
            (sum, rate * steps as f64)
        };
        let h = 1e-4;
        let l = 40.0 * h;
        let (got, rate) = run(400, h, l, 400);
        let want = rate / (l * l * l);
        let ratio = got / want;
        assert!(
            (0.92..1.03).contains(&ratio),
            "uniform: {got} vs {want}: {ratio}"
        );
        // few drops a step, so a changed velocity hardly meets again
        let h = 1e-3;
        let (got, rate) = run(20, h, 0.0, 40_000);
        let want = rate / (h * h * h);
        let ratio = got / want;
        assert!(
            (ratio - 1.0).abs() < 0.06,
            "point: {got} vs {want}: {ratio}"
        );
    }

    /// The Post criteria (D349) at known points: equal drops head-on bounce
    /// below We_d = 2 (4φ′ − 12) = 2.808 (τ = 2, χ = 1) and separate
    /// reflexively above 6 (7 · 2^{2/3} − 8) = 18.68 (η₁ = η₂ = 1), the
    /// value of Ashgriz and Poo for water; the bounce limit rises towards
    /// grazing impacts; at large impact parameters there is no reflexive
    /// separation.
    #[test]
    fn post_criteria() {
        assert!((bounce_weber(1.0, 0.0) - 2.808).abs() < 1e-9);
        let r = 6.0 * (7.0 * 2f64.powf(2.0 / 3.0) - 8.0);
        assert!(
            (reflexive_weber(1.0, 0.0) - r).abs() < 1e-9 * r,
            "{}",
            reflexive_weber(1.0, 0.0)
        );
        assert!((18.6..18.8).contains(&r));
        assert!(bounce_weber(1.0, 0.8) > bounce_weber(1.0, 0.3));
        assert!(bounce_weber(1.0, 0.3) > bounce_weber(1.0, 0.0));
        assert_eq!(reflexive_weber(1.0, 0.9), f64::INFINITY);
        assert!(reflexive_weber(0.5, 0.0).is_finite());
    }

    /// With the Post outcomes slow drops bounce instead of merging (O'Rourke
    /// merges every collision at We → 0), and mass, momentum and enthalpy
    /// still stay; fast head-on drops of the same set separate reflexively
    /// or graze, never all merge.
    #[test]
    fn post_outcomes_bounce_slow_drops() {
        for (speed, want_bounce) in [(0.05, true), (40.0, false)] {
            let list: Vec<_> = (0..30u64)
                .map(|i| {
                    let a = i as f64;
                    let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                    (
                        i,
                        10e-6,
                        50.0,
                        [sign * speed * (1.0 + 0.1 * (a * 0.7).sin()), 0.0, 0.0],
                    )
                })
                .collect();
            let mut p = parcels(&list);
            let before = totals(&p);
            let mut tally = CollisionTally::default();
            for s in 0..200 {
                let mut inp = input(s as f64, true);
                inp.post = true;
                let t = collide(&mut p, &[1e-12], &|_| false, inp);
                tally.events += t.events;
                tally.coalescences += t.coalescences;
                tally.bounces += t.bounces;
                tally.reflexive += t.reflexive;
            }
            let after = totals(&p);
            assert!(tally.events > 10, "{tally:?}");
            if want_bounce {
                assert_eq!(tally.bounces, tally.events, "{tally:?}");
                assert_eq!(tally.coalescences, 0);
            } else {
                assert!(tally.reflexive > 0, "{tally:?}");
                assert!(tally.coalescences < tally.events, "{tally:?}");
            }
            let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
            assert!(rel(after.0, before.0) < 1e-12);
            assert!(rel(after.2, before.2) < 1e-12);
            let scale = before.1.iter().map(|x| x.abs()).fold(0.0, f64::max) + before.0 * speed;
            assert!((after.1[0] - before.1[0]).abs() < 1e-12 * scale);
        }
    }
}

#[cfg(test)]
mod checks {
    use crate::config::{CollisionModel, SpraySection, Steady};

    /// The collisions need the surface tension and the transient spray.
    #[test]
    fn the_collision_table_is_checked() {
        let mut s = SpraySection::default();
        s.collision.model = CollisionModel::Ntc;
        let e = s.validate().unwrap_err().to_string();
        assert!(e.contains("surface tension"), "{e}");
        s.liquid.surface_tension = Some(0.02);
        assert!(s.validate().is_ok());
        s.steady = Some(Steady::default());
        let e = s.validate().unwrap_err().to_string();
        assert!(e.contains("transient"), "{e}");
    }
}
