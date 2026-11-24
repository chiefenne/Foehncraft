//! Communicator abstraction.

use foehncraft_core::exact::{ExactSum, WORDS, max_total, min_total};

/// Communication between the ranks of a run.
///
/// Reductions are exact and order-independent (D31 (1), D93 (4)): every rank
/// receives the same bits, and the bits do not depend on the number of
/// ranks or on the partition. A sum is reduced as the integer image of an
/// [`ExactSum`] (one all-reduce of [`WORDS`] integers), so the reduction
/// tree of the MPI library cannot change it.
pub trait Communicator {
    /// Rank of this process, `0..size()`.
    fn rank(&self) -> usize;
    /// Number of ranks.
    fn size(&self) -> usize;
    /// Waits until every rank has arrived.
    fn barrier(&self);
    /// One value from every rank, in rank order, on every rank.
    fn all_gather_f64(&self, local: f64) -> Vec<f64>;
    /// Sends `send` to rank `dest` and receives a message of the same length
    /// from rank `source` into `recv`. The basic step of a halo exchange.
    fn send_recv_f64(&self, send: &[f64], dest: usize, recv: &mut [f64], source: usize);

    /// The values of every rank, in rank order, on every rank (lengths may
    /// differ). Built from [`Communicator::all_gather_f64`] and
    /// [`Communicator::send_recv_f64`]: the lengths first, then a ring.
    fn all_gather_vec_f64(&self, local: &[f64]) -> Vec<Vec<f64>> {
        let (me, n) = (self.rank(), self.size());
        let lens = self.all_gather_f64(local.len() as f64);
        let mut out: Vec<Vec<f64>> = lens.iter().map(|&l| vec![0.0; l as usize]).collect();
        out[me].copy_from_slice(local);
        for k in 1..n {
            let dest = (me + k) % n;
            let source = (me + n - k) % n;
            let mut recv = vec![0.0; lens[source] as usize];
            self.send_recv_f64(local, dest, &mut recv, source);
            out[source] = recv;
        }
        out
    }

    /// Element-wise integer sum of `buf` over all ranks, in place, on every
    /// rank. Integer addition is associative, so the result does not depend
    /// on the reduction tree. The sums must not overflow `i64`.
    fn all_sum_i64(&self, buf: &mut [i64]);

    /// The merge of the accumulators of all ranks, the same on every rank.
    /// For fewer than 2³⁰ ranks.
    fn merge_exact(&self, local: &ExactSum) -> ExactSum {
        let mut w = local.to_words();
        self.all_sum_i64(&mut w);
        ExactSum::from_words(&w)
    }

    /// The exact sum over all ranks of the local accumulators, rounded
    /// once: bitwise independent of the rank count and the partition.
    fn sum_exact(&self, local: &ExactSum) -> f64 {
        self.merge_exact(local).value()
    }

    /// Several exact sums in one all-reduce (for example the three
    /// components of a force), in the order of `locals`.
    fn sum_exact_many(&self, locals: &[ExactSum]) -> Vec<f64> {
        let mut w: Vec<i64> = locals.iter().flat_map(|l| l.to_words()).collect();
        self.all_sum_i64(&mut w);
        w.as_chunks::<WORDS>()
            .0
            .iter()
            .map(|a| ExactSum::from_words(a).value())
            .collect()
    }

    /// Exact sum over all ranks of one value per rank (independent of the
    /// rank order and count).
    fn sum_reproducible(&self, local: f64) -> f64 {
        let mut s = ExactSum::new();
        s.add(local);
        self.sum_exact(&s)
    }

    /// Maximum over all ranks ([`max_total`]: NaN propagates, +0 > −0).
    fn max(&self, local: f64) -> f64 {
        self.all_gather_f64(local)
            .into_iter()
            .fold(f64::NEG_INFINITY, max_total)
    }

    /// Minimum over all ranks ([`min_total`]: NaN propagates, −0 < +0).
    fn min(&self, local: f64) -> f64 {
        self.all_gather_f64(local)
            .into_iter()
            .fold(f64::INFINITY, min_total)
    }

    /// Name of the implementation, for logs.
    fn name(&self) -> &'static str;
}

/// The communicator of a run on one rank.
#[derive(Debug, Default, Clone, Copy)]
pub struct SerialComm;

impl Communicator for SerialComm {
    fn rank(&self) -> usize {
        0
    }
    fn size(&self) -> usize {
        1
    }
    fn barrier(&self) {}
    fn all_gather_f64(&self, local: f64) -> Vec<f64> {
        vec![local]
    }
    fn send_recv_f64(&self, send: &[f64], dest: usize, recv: &mut [f64], source: usize) {
        assert_eq!((dest, source), (0, 0), "a serial run has only rank 0");
        recv.copy_from_slice(send);
    }
    fn all_sum_i64(&self, _buf: &mut [i64]) {}
    fn name(&self) -> &'static str {
        "serial"
    }
}

#[cfg(feature = "mpi")]
pub mod mpi_comm {
    //! MPI implementation of [`Communicator`].

    use super::Communicator;
    use mpi::Threading;
    use mpi::collective::SystemOperation;
    use mpi::environment::Universe;
    use mpi::topology::{Communicator as _, SimpleCommunicator};
    use mpi::traits::*;

    /// Communicator over `MPI_COMM_WORLD`. MPI is finalised when the value
    /// is dropped.
    ///
    /// MPI is initialised with the threading level "funneled": the rank's
    /// worker threads (rayon) compute, and only the thread that initialised
    /// MPI calls it. The value holds no communicator handle (the world
    /// handle is taken per call), so it can be shared (`Send + Sync`) with
    /// the solver structures of that thread, like the other
    /// communicators.
    pub struct MpiComm {
        universe: Universe,
    }

    impl MpiComm {
        /// Initialises MPI. Returns `None` if MPI was already initialised.
        pub fn init() -> Option<Self> {
            let (universe, _level) = mpi::initialize_with_threading(Threading::Funneled)?;
            Some(Self { universe })
        }

        fn world(&self) -> SimpleCommunicator {
            self.universe.world()
        }

        /// Name and version of the MPI library.
        pub fn library_version() -> String {
            mpi::environment::library_version().unwrap_or_default()
        }
    }

    impl Communicator for MpiComm {
        fn rank(&self) -> usize {
            self.world().rank() as usize
        }
        fn size(&self) -> usize {
            self.world().size() as usize
        }
        fn barrier(&self) {
            self.world().barrier();
        }
        fn all_gather_f64(&self, local: f64) -> Vec<f64> {
            let mut out = vec![0.0; self.size()];
            self.world().all_gather_into(&local, &mut out[..]);
            out
        }
        fn send_recv_f64(&self, send: &[f64], dest: usize, recv: &mut [f64], source: usize) {
            let world = self.world();
            let d = world.process_at_rank(dest as i32);
            let s = world.process_at_rank(source as i32);
            mpi::request::scope(|scope| {
                let sreq = d.immediate_send(scope, send);
                s.receive_into(recv);
                sreq.wait();
            });
        }
        fn all_sum_i64(&self, buf: &mut [i64]) {
            let send = buf.to_vec();
            self.world()
                .all_reduce_into(&send[..], buf, SystemOperation::sum());
        }
        fn name(&self) -> &'static str {
            "mpi"
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serial_reductions() {
        let c = SerialComm;
        assert_eq!(c.sum_reproducible(2.5), 2.5);
        assert_eq!(c.max(-1.0), -1.0);
        let mut r = [0.0; 3];
        c.send_recv_f64(&[1.0, 2.0, 3.0], 0, &mut r, 0);
        assert_eq!(r, [1.0, 2.0, 3.0]);
        let mut acc = ExactSum::new();
        acc.add_slice(&[1e300, 1.0, -1e300]);
        assert_eq!(c.sum_exact(&acc), 1.0);
        assert_eq!(
            c.sum_exact_many(&[acc.clone(), ExactSum::new()]),
            [1.0, 0.0]
        );
        assert_eq!(c.min(3.0), 3.0);
    }
}
