//! Monte Carlo radiation (step 29, R1a): the sampling of the energy bundles
//! and their transport, one source for the CPU and the GPU (D259).
//!
//! - **Points**: Owen-scrambled Sobol points after Burley, *Practical
//!   hash-based Owen scrambling*, JCGT 9(4), 2020: the Sobol sequence in
//!   four dimensions (direction numbers by the Bratley–Fox recurrence from
//!   the primitive polynomials x + 1, x² + x + 1, x³ + x + 1, dimension 0
//!   the van der Corput sequence), its index shuffled and every coordinate
//!   scrambled by a nested uniform scramble with a hash (Laine & Karras).
//!   One emitter's bundles are the first points of one scrambled sequence,
//!   so they are stratified over position and direction.
//! - **Sampling**: a uniform point in a tetrahedron by folding the unit cube
//!   (Rocchini & Cignoni, *Generating random points in a tetrahedron*,
//!   J. Graphics Tools 5(4), 2000), in a triangle by folding the square, an
//!   isotropic direction, and a cosine-weighted direction about a normal in
//!   the branchless orthonormal basis of Duff et al., *Building an
//!   orthonormal basis, revisited*, JCGT 6(1), 2017.
//! - **Transport** ([`transport`]): a bundle walks through the cells
//!   ([`crate::walk`]) and leaves the fraction `1 − exp(−κ ℓ)` of its weight
//!   in each cell it crosses (track-length estimator, D107). Its weight is
//!   an integer number of quanta (D107 (3)): each deposit is rounded once
//!   and subtracted exactly, so emitted = absorbed + escaped holds in
//!   integers, on every call (D107 (4)). Walls are grey diffuse (absorb the
//!   fraction ε, reflect the rest with a cosine-weighted direction), mirrors
//!   (symmetry planes) or open. A weight below `w0 / 1024` plays Russian
//!   roulette; the weight it adds or removes is tallied, so the ledger still
//!   closes exactly.
//!
//! Every random number is keyed by (emitter, bundle, call, event) in the
//! radiation stream (D107 (2), D121), so a bundle's path does not depend on
//! the thread, the rank or the order of the bundles.

use crate::entry::Raw;
use crate::walk::{Walk, WalkEnd, WalkMesh, walk};
use foehncraft_core::math::{cos, exp, floor, sin, sqrt};
use foehncraft_core::rng::{keyed, stream};

/// Sobol direction numbers of one dimension (Bratley–Fox recurrence).
const fn directions(s: u32, a: u32, m: [u32; 3]) -> [u32; 32] {
    let mut v = [0u32; 32];
    let mut k = 0;
    while k < 32 {
        let j = k as u32 + 1; // 1-based
        if s == 0 {
            v[k] = 1u32 << (32 - j);
        } else if j <= s {
            v[k] = m[k] << (32 - j);
        } else {
            let mut x = v[k - s as usize] ^ (v[k - s as usize] >> s);
            let mut i = 1;
            while i < s {
                if (a >> (s - 1 - i)) & 1 == 1 {
                    x ^= v[k - i as usize];
                }
                i += 1;
            }
            v[k] = x;
        }
        k += 1;
    }
    v
}

/// Direction numbers of the four dimensions.
const SOBOL: [[u32; 32]; 4] = [
    directions(0, 0, [0, 0, 0]),
    directions(1, 0, [1, 0, 0]),
    directions(2, 1, [1, 3, 0]),
    directions(3, 1, [1, 3, 1]),
];

/// Coordinate `dim` (0 to 3) of Sobol point `index`, as 32 bits.
#[inline(always)]
pub fn sobol(index: u32, dim: usize) -> u32 {
    let v = &SOBOL[dim];
    let mut x = 0u32;
    let mut i = index;
    let mut b = 0;
    while i != 0 {
        if i & 1 == 1 {
            x ^= v[b];
        }
        i >>= 1;
        b += 1;
    }
    x
}

#[inline(always)]
fn laine_karras(mut x: u32, seed: u32) -> u32 {
    x = x.wrapping_add(seed);
    x ^= x.wrapping_mul(0x6c50_b47c);
    x ^= x.wrapping_mul(0xb82f_1e52);
    x ^= x.wrapping_mul(0xc7af_e638);
    x ^= x.wrapping_mul(0x8d22_f6e6);
    x
}

/// Nested uniform scramble of the bits of `x` (an Owen scramble by hash).
#[inline(always)]
pub fn nested_uniform_scramble(x: u32, seed: u32) -> u32 {
    laine_karras(x.reverse_bits(), seed).reverse_bits()
}

#[inline(always)]
fn hash_combine(seed: u32, v: u32) -> u32 {
    seed ^ v
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
        .wrapping_add(0x9e37_79b9)
}

/// Point `index` of the shuffled, Owen-scrambled Sobol sequence `seed` in
/// [0, 1)⁴.
#[inline(always)]
pub fn owen_sobol4(index: u32, seed: u32) -> [f64; 4] {
    let i = nested_uniform_scramble(index, seed);
    let mut x = [0.0; 4];
    let mut d = 0;
    while d < 4 {
        let s = nested_uniform_scramble(sobol(i, d), hash_combine(seed, d as u32));
        x[d] = f64::from(s) * (1.0 / 4_294_967_296.0);
        d += 1;
    }
    x
}

#[inline(always)]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline(always)]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline(always)]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// A uniform point in the tetrahedron `t` from three uniform numbers.
pub fn point_in_tet(t: &[[f64; 3]; 4], u: [f64; 3]) -> [f64; 3] {
    let [mut s, mut v, mut w] = u;
    if s + v > 1.0 {
        s = 1.0 - s;
        v = 1.0 - v;
    }
    if v + w > 1.0 {
        let tmp = w;
        w = 1.0 - s - v;
        v = 1.0 - tmp;
    } else if s + v + w > 1.0 {
        let tmp = w;
        w = s + v + w - 1.0;
        s = 1.0 - v - tmp;
    }
    let a = 1.0 - s - v - w;
    add(
        add(scale(t[0], a), scale(t[1], s)),
        add(scale(t[2], v), scale(t[3], w)),
    )
}

/// A uniform point in the triangle (a, b, c) from two uniform numbers.
pub fn point_in_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3], u: [f64; 2]) -> [f64; 3] {
    let [mut s, mut v] = u;
    if s + v > 1.0 {
        s = 1.0 - s;
        v = 1.0 - v;
    }
    add(add(scale(a, 1.0 - s - v), scale(b, s)), scale(c, v))
}

/// A direction uniform on the unit sphere from two uniform numbers.
pub fn isotropic(u: [f64; 2]) -> [f64; 3] {
    let z = 1.0 - 2.0 * u[0];
    let r = sqrt((1.0 - z * z).max(0.0));
    let phi = 2.0 * core::f64::consts::PI * u[1];
    [r * cos(phi), r * sin(phi), z]
}

/// A direction about the unit normal `n` with the density cos θ / π
/// (Lambert's law), from two uniform numbers.
pub fn cosine_about(n: [f64; 3], u: [f64; 2]) -> [f64; 3] {
    // orthonormal basis (t, b, n) after Duff et al. 2017
    let sign = if n[2] >= 0.0 { 1.0 } else { -1.0 };
    let a = -1.0 / (sign + n[2]);
    let bb = n[0] * n[1] * a;
    let t = [1.0 + sign * n[0] * n[0] * a, sign * bb, -sign * n[0]];
    let b = [bb, sign + n[1] * n[1] * a, -n[1]];
    let r = sqrt(u[0]);
    let phi = 2.0 * core::f64::consts::PI * u[1];
    let (x, y) = (r * cos(phi), r * sin(phi));
    let z = sqrt((1.0 - u[0]).max(0.0));
    add(add(scale(t, x), scale(b, y)), scale(n, z))
}

/// `1 − exp(−τ)` for τ ≥ 0, accurate also for small τ.
#[inline(always)]
pub fn absorbed_fraction(tau: f64) -> f64 {
    if tau < 1e-5 {
        tau * (1.0 - tau * (0.5 - tau * (1.0 / 6.0)))
    } else {
        1.0 - exp(-tau)
    }
}

/// The fraction `f` of `w` quanta, rounded to the nearest quantum.
#[inline(always)]
fn part(w: i64, f: f64) -> i64 {
    (floor(w as f64 * f + 0.5) as i64).clamp(0, w)
}

/// Kind of a boundary face for the radiation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum WallKind {
    /// Grey diffuse wall: absorbs the fraction ε, reflects the rest with a
    /// cosine-weighted direction (black: ε = 1).
    Grey = 0,
    /// Symmetry plane: specular reflection.
    Mirror = 1,
    /// Open boundary: the bundle leaves the domain.
    Open = 2,
}

/// What a bundle crosses: absorption coefficient per local cell and the
/// boundary faces.
#[derive(Clone, Copy)]
pub struct Medium<'a> {
    /// Absorption coefficient of each local cell (1/m).
    pub kappa: &'a [f64],
    /// Kind of each boundary face, indexed by face − first boundary face.
    pub wall: &'a [WallKind],
    /// Emissivity of each boundary face (same index).
    pub emissivity: &'a [f64],
    /// Unit normal into the domain of each boundary face (same index).
    pub normal: &'a [[f64; 3]],
    /// First boundary face.
    pub first_boundary: usize,
    /// Length of a ray segment: longer than any chord of the domain.
    pub reach: f64,
    /// Tetrahedra a walk may visit before it counts as lost.
    pub max_steps: u32,
}

/// An energy bundle in flight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bundle {
    /// Where it is and where it goes.
    pub walk: Walk,
    /// Weight in quanta.
    pub w: i64,
    /// Weight at emission, in quanta (the roulette threshold is `w0/1024`).
    pub w0: i64,
    /// Global identity of the emitter.
    pub emitter: u64,
    /// Index of the bundle among the emitter's bundles.
    pub index: u32,
    /// Events so far (wall interactions, roulette).
    pub event: u32,
    /// Batch of the statistical error estimate.
    pub batch: u32,
}

/// Where the tallies of a call go: integer quanta, so the order of the
/// additions does not matter (atomics on the GPU, atomics or per-thread
/// arrays on the CPU).
pub trait Tally {
    /// `q` quanta absorbed in local cell `cell`.
    fn cell(&self, batch: u32, cell: usize, q: i64);
    /// `incident` quanta arrive at boundary face `face` (index from the
    /// first boundary face), `absorbed` of them stay.
    fn wall(&self, batch: u32, face: usize, incident: i64, absorbed: i64);
    /// `q` quanta leave through an open boundary.
    fn escaped(&self, q: i64);
    /// Russian roulette removed (`q > 0`) or added (`q < 0`) quanta.
    fn roulette(&self, q: i64);
    /// A bundle with `q` quanta was lost (degenerate geometry).
    fn lost(&self, q: i64);
}

/// How [`transport`] ended for a bundle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fate {
    /// The bundle's weight is all tallied.
    Done,
    /// The bundle stands at a processor face, in the halo cell behind it:
    /// hand it to the rank that owns the cell.
    HandOver,
}

/// A uniform number of bundle (emitter, index) for an event and a purpose:
/// the (identity, step, event) key of D121 in the radiation stream.
#[inline(always)]
fn uniform(emitter: u64, index: u32, call: u64, event: u32, purpose: u64) -> f64 {
    let step = (call << 32) | u64::from(index);
    let event = (u64::from(event) << 8) | purpose;
    (keyed([stream::RADIATION, emitter, step, event]) >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Russian roulette on a bundle of `w < w_min` quanta: it survives with
/// probability `w / w_min` and then carries `w_min`. Returns the new
/// weight; the quanta removed (or, negative, added) go to the tally.
#[inline(always)]
fn roulette<T: Tally>(w: i64, w_min: i64, u: f64, tally: &T) -> i64 {
    if u * (w_min as f64) < w as f64 {
        tally.roulette(w - w_min);
        w_min
    } else {
        tally.roulette(w);
        0
    }
}

/// Seed of the scrambled Sobol sequence of an emitter for one purpose.
#[inline(always)]
pub fn sobol_seed(emitter: u64, call: u64, purpose: u64) -> u32 {
    (keyed([
        stream::RADIATION,
        emitter,
        (call << 32) | 0xffff_ffff,
        purpose,
    ]) >> 32) as u32
}

/// Purposes of the keyed numbers.
const ROULETTE: u64 = 1;
const REFLECT: u64 = 2;

/// Events (wall interactions and roulette) after which a bundle counts as
/// lost: a guard against endless reflections, far above any physical case
/// (a bundle loses at least the fraction ε at a grey wall).
pub const MAX_EVENTS: u32 = 1 << 20;

/// Moves bundle `b` until its weight is all tallied or it reaches a
/// processor face. `call` is the number of the radiation call.
pub fn transport<T: Tally>(
    m: &WalkMesh<'_>,
    md: &Medium<'_>,
    b: &mut Bundle,
    call: u64,
    tally: &T,
) -> Fate {
    let w_min = (b.w0 >> 10).max(1);
    let (emitter, index, batch) = (b.emitter, b.index, b.batch);
    loop {
        let d = [
            b.walk.target[0] - b.walk.origin[0],
            b.walk.target[1] - b.walk.origin[1],
            b.walk.target[2] - b.walk.origin[2],
        ];
        let len = sqrt(dot(d, d));
        let (mut w, mut event) = (b.w, b.event);
        let end = walk(m, &mut b.walk, md.max_steps, |cell, a, c| {
            let tau = md.kappa[cell as usize] * ((c - a) * len);
            let q = part(w, absorbed_fraction(tau));
            if q > 0 {
                tally.cell(batch, cell as usize, q);
                w -= q;
                if w < w_min {
                    let u = uniform(emitter, index, call, event, ROULETTE);
                    event += 1;
                    w = roulette(w, w_min, u, tally);
                }
            }
            w > 0
        });
        b.w = w;
        b.event = event;
        match end {
            WalkEnd::Stopped => return Fate::Done,
            WalkEnd::Processor { .. } => return Fate::HandOver,
            WalkEnd::Lost | WalkEnd::Reached => {
                // the segment is longer than any chord, so Reached means a
                // broken walk as well
                tally.lost(b.w);
                b.w = 0;
                return Fate::Done;
            }
            WalkEnd::Boundary { face } => {
                let f = face as usize - md.first_boundary;
                let n = md.normal[f];
                let s = 1.0 / len;
                let dir = [d[0] * s, d[1] * s, d[2] * s];
                let new_dir = match md.wall[f] {
                    WallKind::Open => {
                        // the flux through the open face, none of it stays
                        tally.wall(batch, f, b.w, 0);
                        tally.escaped(b.w);
                        b.w = 0;
                        return Fate::Done;
                    }
                    WallKind::Mirror => {
                        b.event += 1;
                        if b.event >= MAX_EVENTS {
                            tally.lost(b.w);
                            b.w = 0;
                            return Fate::Done;
                        }
                        let k = 2.0 * dot(dir, n);
                        [dir[0] - k * n[0], dir[1] - k * n[1], dir[2] - k * n[2]]
                    }
                    WallKind::Grey => {
                        let q = part(b.w, md.emissivity[f]);
                        tally.wall(batch, f, b.w, q);
                        b.w -= q;
                        if b.w == 0 {
                            return Fate::Done;
                        }
                        let u = [
                            uniform(emitter, index, call, b.event, REFLECT),
                            uniform(emitter, index, call, b.event, REFLECT | 0x80),
                        ];
                        b.event += 1;
                        if b.event >= MAX_EVENTS {
                            tally.lost(b.w);
                            b.w = 0;
                            return Fate::Done;
                        }
                        if b.w < w_min {
                            let u = uniform(emitter, index, call, b.event, ROULETTE);
                            b.event += 1;
                            b.w = roulette(b.w, w_min, u, tally);
                            if b.w == 0 {
                                return Fate::Done;
                            }
                        }
                        cosine_about(n, u)
                    }
                };
                let x = b.walk.x;
                let target = add(x, scale(new_dir, md.reach));
                b.walk = Walk::new(x, b.walk.cell, b.walk.face, b.walk.edge, target);
            }
        }
    }
}

/// Words of a bundle in the flat records of the device path
/// ([`Bundle::to_words`]).
pub const BUNDLE_WORDS: usize = 22;

impl Bundle {
    /// The bundle as 64-bit words: the bits of origin, target, λ and x, then
    /// cell, face, edge, entered, steps, crossed, the weights, the emitter,
    /// the index, the event and the batch.
    pub fn to_words(&self) -> [u64; BUNDLE_WORDS] {
        let w = &self.walk;
        [
            w.origin[0].to_bits(),
            w.origin[1].to_bits(),
            w.origin[2].to_bits(),
            w.target[0].to_bits(),
            w.target[1].to_bits(),
            w.target[2].to_bits(),
            w.lambda.to_bits(),
            w.x[0].to_bits(),
            w.x[1].to_bits(),
            w.x[2].to_bits(),
            u64::from(w.cell),
            u64::from(w.face),
            u64::from(w.edge),
            u64::from(w.entered),
            u64::from(w.steps),
            u64::from(w.crossed),
            self.w as u64,
            self.w0 as u64,
            self.emitter,
            u64::from(self.index),
            u64::from(self.event),
            u64::from(self.batch),
        ]
    }

    /// The bundle of [`Bundle::to_words`] (`v` holds at least
    /// [`BUNDLE_WORDS`] words).
    pub fn from_words(v: &[u64]) -> Self {
        let f = |k: usize| f64::from_bits(v[k]);
        Self {
            walk: Walk {
                origin: [f(0), f(1), f(2)],
                target: [f(3), f(4), f(5)],
                lambda: f(6),
                x: [f(7), f(8), f(9)],
                cell: v[10] as u32,
                face: v[11] as u32,
                edge: v[12] as u32,
                entered: v[13] as u32,
                steps: v[14] as u32,
                crossed: v[15] as u32,
            },
            w: v[16] as i64,
            w0: v[17] as i64,
            emitter: v[18],
            index: v[19] as u32,
            event: v[20] as u32,
            batch: v[21] as u32,
        }
    }
}

/// Moves bundle `i` of the flat records `words` ([`Bundle::to_words`]) to
/// its end: the body of the device path, the same on the CPU and the GPU.
pub fn transport_words<T: Tally>(
    m: &WalkMesh<'_>,
    md: &Medium<'_>,
    words: &[u64],
    i: usize,
    call: u64,
    tally: &T,
) -> Fate {
    let mut b = Bundle::from_words(&words[i * BUNDLE_WORDS..(i + 1) * BUNDLE_WORDS]);
    transport(m, md, &mut b, call, tally)
}

/// Argument block of the PTX entry `radiation_transport`
/// (`gpu/ptx-kernels/src/radiation.rs`), one GPU thread per bundle: the
/// walk mesh, the medium, the bundles and the integer tallies, as device
/// addresses ([`Raw`]). The host side is `crates/foehncraft-cuda/src/radiation.rs`.
/// The tallies are batch-major (`batch × n + i`); `ledger` holds the
/// escaped, roulette and lost quanta, all as two's-complement sums.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct TransportArgs {
    /// Node coordinates.
    pub points: Raw<[f64; 3]>,
    /// Face node offsets.
    pub face_ptr: Raw<u32>,
    /// Face nodes.
    pub face_nodes: Raw<u32>,
    /// Owner of each face.
    pub owner: Raw<u32>,
    /// Neighbour of each interior face.
    pub neighbour: Raw<u32>,
    /// Cell face offsets.
    pub cell_ptr: Raw<u32>,
    /// Packed cell faces.
    pub cell_faces: Raw<u32>,
    /// Apex of each cell.
    pub apex: Raw<[f64; 3]>,
    /// Fan point of each face.
    pub fan: Raw<[f64; 3]>,
    /// Number of local cells.
    pub n_cells: u64,
    /// First periodic face.
    pub periodic_start: u64,
    /// Absorption coefficient per cell.
    pub kappa: Raw<f64>,
    /// [`WallKind`] of each boundary face as its `u8` value.
    pub wall: Raw<u8>,
    /// Emissivity of each boundary face.
    pub emissivity: Raw<f64>,
    /// Inward unit normal of each boundary face.
    pub normal: Raw<[f64; 3]>,
    /// First boundary face.
    pub first_boundary: u64,
    /// Segment length.
    pub reach: f64,
    /// Step limit of a walk.
    pub max_steps: u64,
    /// Bundles, [`BUNDLE_WORDS`] words each.
    pub bundles: Raw<u64>,
    /// Number of bundles.
    pub n_bundles: u64,
    /// Number of the call.
    pub call: u64,
    /// Number of boundary faces.
    pub n_boundary: u64,
    /// Quanta absorbed per batch and cell.
    pub cell: Raw<u64>,
    /// Quanta arriving per batch and boundary face.
    pub incident: Raw<u64>,
    /// Quanta absorbed per batch and boundary face.
    pub absorbed: Raw<u64>,
    /// Escaped, roulette and lost quanta.
    pub ledger: Raw<u64>,
}

/// Moves bundle `i` of the argument block (the body of the PTX entry).
///
/// # Safety
/// The addresses of `a` hold the arrays of one walk mesh and medium (as
/// for [`WalkMesh`] and [`Medium`], `wall` with valid [`WallKind`] values)
/// and `n_bundles` records; `i < n_bundles`. The tallies are written only
/// through `tally`.
pub unsafe fn transport_entry<T: Tally>(i: usize, a: &TransportArgs, tally: &T) -> Fate {
    // SAFETY: guaranteed by the caller; WallKind is repr(u8) and the host
    // writes only its values into `wall`.
    unsafe {
        let m = WalkMesh {
            points: a.points.get(),
            face_ptr: a.face_ptr.get(),
            face_nodes: a.face_nodes.get(),
            owner: a.owner.get(),
            neighbour: a.neighbour.get(),
            cell_ptr: a.cell_ptr.get(),
            cell_faces: a.cell_faces.get(),
            apex: a.apex.get(),
            fan: a.fan.get(),
            n_cells: a.n_cells as usize,
            periodic_start: a.periodic_start as usize,
        };
        let wall = if a.wall.len == 0 {
            &[]
        } else {
            core::slice::from_raw_parts(a.wall.ptr as usize as *const WallKind, a.wall.len as usize)
        };
        let md = Medium {
            kappa: a.kappa.get(),
            wall,
            emissivity: a.emissivity.get(),
            normal: a.normal.get(),
            first_boundary: a.first_boundary as usize,
            reach: a.reach,
            max_steps: a.max_steps as u32,
        };
        transport_words(&m, &md, a.bundles.get(), i, a.call, tally)
    }
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;
    use std::{vec, vec::Vec};

    /// Each coordinate of the first 2^m points, unscrambled and scrambled,
    /// has exactly one point in each of the 2^m intervals of width 2^-m
    /// (a (0, 1)-sequence).
    #[test]
    fn every_coordinate_is_stratified() {
        let m = 10;
        let n = 1usize << m;
        for d in 0..4 {
            let mut hit = vec![0u32; n];
            for i in 0..n as u32 {
                hit[(sobol(i, d) >> (32 - m)) as usize] += 1;
            }
            assert!(hit.iter().all(|&h| h == 1), "dimension {d}");
        }
        for seed in [0u32, 7, 0xdead_beef] {
            for d in 0..4 {
                let mut hit = vec![0u32; n];
                for i in 0..n as u32 {
                    hit[(owen_sobol4(i, seed)[d] * n as f64) as usize] += 1;
                }
                assert!(hit.iter().all(|&h| h == 1), "seed {seed}, dimension {d}");
            }
        }
    }

    /// The first two coordinates form a (0, 2)-sequence: the first 2^m
    /// points put one point in every elementary box of area 2^-m, also
    /// after the scramble.
    #[test]
    fn first_two_coordinates_are_a_02_net() {
        let m = 8u32;
        let n = 1usize << m;
        for seed in [0u32, 3, 0x1234_5678] {
            let pts: Vec<[f64; 4]> = (0..n as u32).map(|i| owen_sobol4(i, seed)).collect();
            for a in 0..=m {
                let (nx, ny) = (1usize << a, 1usize << (m - a));
                let mut hit = vec![0u32; n];
                for p in &pts {
                    let (i, j) = ((p[0] * nx as f64) as usize, (p[1] * ny as f64) as usize);
                    hit[i * ny + j] += 1;
                }
                assert!(
                    hit.iter().all(|&h| h == 1),
                    "seed {seed}, boxes {nx} x {ny}"
                );
            }
        }
    }

    #[test]
    fn bundle_words_round_trip() {
        let b = Bundle {
            walk: Walk::new([0.1, -2.5, 3.0], 17, 99, 3, [1e300, -0.0, 5e-324]),
            w: (1 << 52) - 3,
            w0: 1 << 40,
            emitter: (1 << 50) + 7,
            index: u32::MAX,
            event: 12,
            batch: 7,
        };
        let mut c = b;
        c.walk.lambda = 0.375;
        c.walk.steps = 40_000;
        c.walk.crossed = 12;
        c.walk.entered = 2;
        for x in [b, c] {
            assert_eq!(Bundle::from_words(&x.to_words()), x);
        }
    }

    #[test]
    fn samples_have_the_right_moments() {
        let n = 1u32 << 14;
        // uniform in a tetrahedron: mean at the centroid
        let t = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.3, 0.2, 1.5],
        ];
        let mut mean = [0.0; 3];
        for i in 0..n {
            let u = owen_sobol4(i, 11);
            let x = point_in_tet(&t, [u[1], u[2], u[3]]);
            for k in 0..3 {
                mean[k] += x[k] / f64::from(n);
            }
        }
        for k in 0..3 {
            let c = (t[0][k] + t[1][k] + t[2][k] + t[3][k]) / 4.0;
            assert!((mean[k] - c).abs() < 2e-3, "{mean:?}");
        }
        // cosine-weighted about a tilted normal: all inward, mean cos θ = 2/3
        let nrm = {
            let v = [0.3f64, -0.8, -0.52];
            let l = dot(v, v).sqrt();
            [v[0] / l, v[1] / l, v[2] / l]
        };
        let mut mc = 0.0;
        let mut iso = [0.0; 3];
        for i in 0..n {
            let u = owen_sobol4(i, 5);
            let d = cosine_about(nrm, [u[0], u[1]]);
            assert!((dot(d, d) - 1.0).abs() < 1e-12);
            let c = dot(d, nrm);
            assert!(c >= 0.0);
            mc += c / f64::from(n);
            let s = isotropic([u[2], u[3]]);
            assert!((dot(s, s) - 1.0).abs() < 1e-12);
            for k in 0..3 {
                iso[k] += s[k] / f64::from(n);
            }
        }
        assert!((mc - 2.0 / 3.0).abs() < 1e-3, "{mc}");
        assert!(iso.iter().all(|v| v.abs() < 1e-3), "{iso:?}");
        // the absorbed fraction is continuous where the series takes over
        let (a, b) = (
            absorbed_fraction(1e-5 * (1.0 - 1e-12)),
            absorbed_fraction(1e-5),
        );
        assert!((a - b).abs() < 1e-15, "{a} {b}");
    }
}
